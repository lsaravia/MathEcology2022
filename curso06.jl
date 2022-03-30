#= ABC con el modelo de bosque2D 

* Elegir un estadístico de comparación

* Elegir una distancia

* Elegir una distribución a priori para los parámetros

* Generar datos con el modelo

* Verificar que el algoritmo ajusta a los valores usados como datos

* Ajustar con datos reales"

=#

"""
    initialize_bosque(;fil::Int=200, col::Int=200, densidad::Float64=0.1)

Inicializa el modelo de bosque con `fil` filas y `col` columnas con una densidad 
evaluada de forma estocastica con el parámetro `densidad`

"""
function initialize_bosque(;fil::Int=200, col::Int=200, densidad::Float64=0.1)
	
	landscape = zeros(Int,fil,col)
	
	for i in 1:fil, j in 1:col 
		if rand() ≤ densidad
			landscape[i,j] = 1 
		end
	end
	
	return landscape
	
end

"""
    plot_bosque(m::Matrix{Int64})

Gráfico de la distribución espacial del bosque

"""
function plot_bosque(m::Matrix{Int64})
	heatmap(m,aspect_ratio=1,legend=:none,xticks=:none,
		yticks=:none,framestyle=:none,color=[ :black , :green])
end

"""
    step_bosque!(landscape::Matrix{Int64},vecindad, parms)

Realiza un paso del modelo de bosque con parametros `parms` con tasa de crecimiento λ y mortalidad δ. Este modelo es identico al llamado proceso de contacto.

"""
function step_bosque!(landscape::Matrix{Int64},vecindad, parms)
	λ, δ = parms
	
	fil = size(landscape,1)
	col = size(landscape,2)
	R = λ + δ
	N = R * fil * col
	for z in 1:N 
		i = rand(1:fil)
		j = rand(1:col)
		if landscape[i,j] == 1
			if rand() ≤ λ/R
				x= [i,j] .+ rand(vecindad)
				if !(x[1] < 1 || x[1] > fil || x[2] < 1 || x[2]>col)
					landscape[x[1],x[2]] = 1
				end
			else
				landscape[i,j] = 0
			end
		end
	end
	return nothing
end

vecindad_MP = ((0,1),(0,-1),(1,0),(-1,0) )

function run_bosque(vecindad, parms; size=100, pasos=100, plot::Bool=false)
	pob = zeros(pasos)
	m = initialize_bosque(fil=size,col=size)
	N = length(m)
	pob[1] = sum(m)/N 
	for i in 2:pasos
		step_bosque!(m,vecindad,parms)
		if plot 
			#= 
			usar display para graficar dentro de una funcion que retorna otra cosa que no sea el plot
			=#
			display(plot_bosque(m))
		end
		pob[i]=sum(m)/N 
	end
	return pob
end

#
# ABC
#

# Definir distancia

function distance(data::Vector{Float64},simulation::Vector{Float64})
    return sqrt(sum((data.-simulation).^2.0))
end


# Estimar threshold
#

ts1 = run_bosque(vecindad_MP, (2,.4), size=200, pasos=200)
ts2 = run_bosque(vecindad_MP, (2,.4), size=200, pasos=200)

vecDens = ts1[1+end-15:end]
vecDens2 = ts2[1+end-15:end]
distance(vecDens,vecDens2)

#
# ABC
#
function ABC_bosque(data,d_λ,d_δ,d_size,threshold)
	
	params = zeros(Float64, 3, 100)
	successes = 0
	trials = 0 
	
	while successes < size(params, 2 )
	    trials += 1 
	    # Candidate sample
	    t_λ = rand(d_λ)
	    t_δ = rand(d_δ)
	    t_size = rand(d_size)
	    if t_λ/t_δ > 1.5 
	        # Create timeseries
	        ts = run_bosque(vecindad_MP, (t_λ,t_δ), size=t_size, pasos=200)
	        last_ts = ts[1+end-length(data):end]
	        score = distance(data,last_ts)
	        if score < threshold
	            successes += 1
	            params[:,successes] .= [t_λ, t_δ, t_size]
	            @info "Tengo $successes exitos en $trials intentos! "
	        end 
	    end
	end
	return params
end

using Distributions
d_λ = Uniform(0.2,5)
d_δ = Uniform(0.1,2)
d_size = [50,100,150,200]
p_parms = ABC_bosque(vecDens,d_λ,d_δ,d_size, 0.1)

#
# Comparación priori/posteriori
#
using Plots, StatsPlots
plot(d_λ ,lab=" (Prior)",legend = :outertopright) 
density!(p_parms[1,:], lab= "λ (Posterior)")
xaxis!((0.2,5))

plot(d_δ ,lab=" (Prior)",legend = :outertopright) 
density!(p_parms[2,:], lab= "δ (Posterior)")
xaxis!((0.1,2))

histogram(p_parms[3,:],bins=10)

median(p_parms[1,:])
median(p_parms[2,:])
median(p_parms[1,:])/median(p_parms[2,:]) 
mode(p_parms[3,:])

ts = run_bosque(vecindad_MP, (3.36,.66), size=100, pasos=200)
vecDens3 = ts[1+end-15:end]
plot(vecDens, label="Data")
plot!(vecDens3, label="Model")


#
# Ajustar a datos reales
#
# Saravia, L. A., Doyle, S., and Bond-Lamberty, B. (2018). Power laws and critical fragmentation in global forests. 
# Scientific Reports 8, 17766. doi:10.1038/s41598-018-36120-w.
#

using CSV
using DataFrames, DataFramesMeta

#
# https://github.com/lsaravia/MathEcology2022/blob/master/Data/GlobalForestCoverThresholded.csv
# 

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

#
# Ver cual es el valor posible del umbral
#

ts2 = run_bosque(vecindad_MP, (5,.4), size=200, pasos=200)
vecDens2 = ts2[1+end-length(vecDens):end]
distance(vecDens,vecDens2)
mean(vecDens2)

# Definir apriori
#
d_λ = Uniform(1,6)
d_δ = Uniform(0.1,2)
d_size = [50,100,150,200]

# Ejecutar ABC
#
#
@time p_parms = ABC_bosque(vecDens,d_λ,d_δ,d_size, 0.4)


#
# Comparación priori/posteriori
#
plot(d_λ ,lab=" (Prior)",legend = :outertopright) 
density!(p_parms[1,:], lab= "λ (Posterior)")
xaxis!((1,6))

plot(d_δ ,lab=" (Prior)",legend = :outertopright) 
density!(p_parms[2,:], lab= "δ (Posterior)")
xaxis!((0.1,2))

histogram(p_parms[3,:],bins=10)

median(p_parms[1,:])
median(p_parms[2,:])
median(p_parms[1,:])/median(p_parms[2,:]) 
mode(p_parms[3,:])

ts = run_bosque(vecindad_MP, (4.83,.35), size=50, pasos=200)
vecDens3 = ts[1+end-length(vecDens):end]
plot(vecDens, label="Data")
plot!(vecDens3, label="Model")

plot!(ts[1:length(vecDens)], label="Transiente")


# Definir apriori
#
d_λ = Truncated(Normal(5,2), 1,6)
plot(d_λ ,lab=" (Prior)",legend = :outertopright) 

d_δ = Truncated(Normal(.5,2), .1,2)
plot(d_δ ,lab=" (Prior)",legend = :outertopright) 

d_size = [50,100,150,200]
@time p_parms = ABC_bosque(vecDens,d_λ,d_δ,d_size, 0.3)

plot(d_λ ,lab=" (Prior)",legend = :outertopright) 
density!(p_parms[1,:], lab= "λ (Posterior)")
xaxis!((1,6))

plot(d_δ ,lab=" (Prior)",legend = :outertopright) 
density!(p_parms[2,:], lab= "δ (Posterior)")
xaxis!((0.1,2))

histogram(p_parms[3,:],bins=10)

median(p_parms[1,:])
median(p_parms[2,:])
median(p_parms[1,:])/median(p_parms[2,:]) 
mode(p_parms[3,:])

ts = run_bosque(vecindad_MP, (4.67,.31), size=200, pasos=200)
vecDens3 = ts[1+end-length(vecDens):end]
plot(vecDens, label="Data")
plot!(vecDens3, label="Model")

