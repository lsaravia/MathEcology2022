#
# ABC 
#

# Modelo
neighbors_VN = ((0,1),(0,-1),(1,0),(-1,0) )
"""
    step_bosque!(landscape::Matrix{Int}, neighborhood, parms)

Ejecuta un paso del modelo de bosques identico al proceso de contacto[1]. 
Los parametros son λ = tasa de crecimiento y δ la tasa de mortalidad. 
El punto crítico del modelo de campo medio es λ/δ=1

1. Oborny, B., Szabó, G., and Meszéna, G. (2007). “Survival of species in patchy landscapes: percolation in space and time,” 
in Scaling Biodiversity (Cambridge University Press), 409–440. Available at: http://dx.doi.org/10.1017/CBO9780511814938.022.

"""
function step_bosque!(landscape::Matrix{Int}, neighborhood, parms)
    λ, δ = parms
    
    # Total rate of events
    R = λ + δ         

    # Size of the grid 
    fil = size(landscape,1)
    col = size(landscape,2)
    
    # Number of iterations  to complete one time step
    N = fil*col*R

    for z in 1:N
        i= rand(1:fil) 
        j= rand(1:col)
        if landscape[i,j]==1
            if rand() ≤ λ/R
                x=[i,j] .+  rand(neighborhood)
                if !(x[1] > fil || x[1] <1 || x[2] > col || x[2]<1)
                    landscape[x[1],x[2]] = 1
                end
            else
                landscape[i,j] = 0
            end
        end
    end
    return
end

"""
    initialize_bosque(;fil::Int=200,col::Int=200, density::Float64=0.1)

Genera una matriz de nros enteros `fil x col` con una densidad `density` de 1

"""
function initialize_bosque(;fil::Int=200,col::Int=200, density::Float64=0.1)
	
	space = zeros(Int64,fil,col)
    ini_dens = trunc(Int, fil*col*density)
    for i in 1:fil, j in 1:col
        if rand() ≤ density        
            space[i,j] = 1
        end
    end
    return space
end


#
#  Ejecutar el modelo de bosque 
#
function run_bosque(neighborhood, parms; size=200, pasos=100, plot::Bool=false)
	pob = zeros(pasos,2) 
	m = initialize_bosque(fil=size,col=size)
	N = length(m)
	pob[1,1] = sum(m)/N 
    pob[1,2] = var(m)/pob[1,1]
	for i in 2:pasos
		step_bosque!(m,neighborhood,parms)
		if plot 
			#= 
			usar display para graficar dentro de una funcion que retorna otra cosa que no sea el plot
			=#
			display(plot_bosque(m))
		end
		pob[i,1]=sum(m)/N          # densidad en el espacio 
        pob[i,2]=var(m)/pob[i,1]   # Coeficiente de variacion = Varianza/media
	end
	return pob
end


#
# Datos
#

# * Opcion 1 generar datos a partir de una distribución
# using Distributions
# plot(rand(Laplace(50,4),100))

# * Opcion 2 generar datos a partir del mismo modelo

# * Opcion 3 tomar datos reales

using CSV
using DataFrames, DataFramesMeta


bosqData = CSV.File("Data/GlobalForestCoverThresholded.csv") |> DataFrame

# mostrar las Regiones
#
@show @chain bosqData begin 
    unique( :regsub)
    @select :regsub
end

# Filtrar la region de interes
#
bD = @chain bosqData begin
    @subset(:regsub .== "SAT1", :threshold.== 30 ) 
    @select(:regsub,:threshold,:number_of_patches, :total_patch_area, :year)
    groupby( :threshold)
    @transform :density = :total_patch_area ./ maximum(:total_patch_area)
end

# Usar statPlots para plotear los dataframes
#
using StatsPlots
@df bD plot(:year,:density, group=:threshold, legend = :outertopright)

# Elegir un estadístico / distancia

vecDens = @chain bD begin
            _.density
end

#=
function distance(data::Vector{Float64},simulation::Vector{Float64})
    x = [mean(data), var(data)/mean(data)]
    y = [mean(simulation), var(simulation)/mean(simulation)]
    return sqrt(sum((x.-y).^2.0))
end
=#

function distance(data::Vector{Float64},simulation::Vector{Float64})
    return sqrt(sum((data.-simulation).^2.0))
end

# Cutoff to reject the sample

threshold  = 10

# Parameters priors
using Distributions
d_λ = Uniform(0.2,5)
d_δ = Uniform(0.1,2)
d_size = [50,100,150,200]

# Cuanto vale el threshold???
#
ts = run_bosque(neighbors_VN, (5,.4), size=200, pasos=200)
ts1 = run_bosque(neighbors_VN, (5,.4), size=200, pasos=200)

vecDens2 = ts[1+end-length(vecDens):end,1]
vecDens1 = ts1[1+end-length(vecDens):end,1]
distance(vecDens1,vecDens2)
distance(vecDens1,vecDens)
distance(vecDens2,vecDens)

mean(vecDens2)
mean(vecDens)
#
# ABC - initialization 
# Priors
# Pre-allocate space for the retained samples
#
d_λ = Uniform(0.2,5)
d_δ = Uniform(0.1,2)
d_size = [50,100,150,200]
params= zeros(Float64, 3, 100)
successes = 0
trials = 0 
threshold  = 0.40


@time while successes < size(params, 2 )
    trials += 1 
    # Candidate sample
    t_λ = rand(d_λ)
    t_δ = rand(d_δ)
    t_size = rand(d_size)
    if t_λ/t_δ > 1.5 
        # Create timeseries
        ts = run_bosque(neighbors_VN, (t_λ,t_δ), size=t_size, pasos=200)
        last_ts = ts[1+end-length(vecDens):end,1]
        score = distance(vecDens,last_ts)
        if score < threshold
            successes += 1
            params[:,successes] .= [t_λ, t_δ, t_size]
            @info "Got $successes good samples"
        end 
    end
end

# Accepted vs Trials
# Time = 1240 sec  suc/tri = 0.01
successes/trials 

# Comparison prior/posterior
#
plot(d_λ ,lab="λ (Prior)") 
density!(params[1,:], lab= "λ (Posterior)")
xaxis!((0.2,5))

plot(d_δ ,lab="δ (Prior)") 
density!(params[2,:], lab= "δ (Posterior)")
xaxis!((0.1,2))


histogram(params[3,:])

median(params[1,:])
median(params[2,:])
median(params[1,:])/median(params[2,:])

median(params[3,:])

ts = run_bosque(neighbors_VN, (1.83,.24), size=100, pasos=200)
vecDens2 = ts[1+end-length(vecDens):end,1]
@df bD plot(:year,:density, group=:threshold, legend = :outertopright)
plot!(collect(2000:2015),vecDens2)

#
# ABC con tiempo como parametro 
# Asumimos que los datos pueden estar en el transiente
#
#
# ABC - initialization 
# Priors
# Pre-allocate space for the retained samples
#
d_λ = Uniform(0.2,5)
d_δ = Uniform(0.1,2)
d_time = collect(1:(101- length(vecDens)))
params= zeros(Float64, 3, 100)
successes = 0
trials = 0 
threshold  = 0.40


@time while successes < size(params, 2 )
    trials += 1 
    # Candidate sample
    t_λ = rand(d_λ)
    t_δ = rand(d_δ)
    t_time = rand(d_time)
    if t_λ/t_δ > 1.5 
        # Create timeseries
        ts = run_bosque(neighbors_VN, (t_λ,t_δ), size=100, pasos=200)
        last_ts = ts[t_time:t_time+15,1]
        score = distance(vecDens,last_ts)
        if score < threshold
            successes += 1
            params[:,successes] .= [t_λ, t_δ, t_time]
            @info "Got $successes good samples"
        end 
    end
end

# Accepted vs Trials
# Time = 317 sec  suc/tri = 0.15
successes/trials 

# Comparison prior/posterior
#
plot(d_λ ,lab="λ (Prior)") 
density!(params[1,:], lab= "λ (Posterior)")
xaxis!((0.2,5))

plot(d_δ ,lab="δ (Prior)") 
density!(params[2,:], lab= "δ (Posterior)")
xaxis!((0.1,2))

histogram(params[3,:],bins=10)

median(params[1,:])
median(params[2,:])
median(params[1,:])/median(params[2,:])

median(params[3,:])

ts = run_bosque(neighbors_VN, (3.74,.28), size=100, pasos=200)
vecDens2 = ts[35:35+15,1]
@df bD plot(:year,:density, group=:threshold, legend = :outertopright)
plot!(collect(2000:2015),vecDens2)
