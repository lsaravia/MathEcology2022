#
# Modelo estocastico espacial de una población
#
using Plots
using Statistics
using Distributions

fil, col =  (100,100)
space = zeros(Int64,fil,col)

[ space[rand(1:fil) , rand(1:col)] = 1 for i in 1:10]

neighbors_VN = ((0,1),(0,-1),(1,0),(-1,0) )


for i in 1:fil
    for j in 1:col
        if space[i,j]==1
            x=(i,j) .+  rand(neighbors_VN)
            if x[1] > fil || x[1] <1 || x[2] > col || x[2]<1
                @info "Outside space $x"
            else
                space[x[1],x[2]] = 1
            end
        end
    end
end

heatmap(space,aspect_ratio=1,legend=:none,xticks=:none,
yticks=:none,framestyle=:none,color=[ :black , :green])




#
#  Construir las primeras funciones
#


function initialize_bosque(;fil::Int=200,col::Int=200, density::Float64=0.1)
	
	space = zeros(Int64,fil,col)
    ini_dens = trunc(Int, fil*col*density)
	[ space[rand(1:fil) , rand(1:col)] = 1 for i in 1:ini_dens]
    # no anda si no retorna space 
    return space
end


function step_bosque!(landscape::Matrix{Int}, neighborhood)
    @show neighborhood
    for i in 1:fil
        for j in 1:col
            if landscape[i,j]==1
                x=[i,j] .+  rand(neighborhood)
                if x[1] > fil || x[1] <1 || x[2] > col || x[2]<1
                    @info "Outside space $x"
                else
                    landscape[x[1],x[2]] = 1
                end
            end
        end
    end
    return
end

#
# 
#
bosq = initialize_bosque()
sum(bosq)
mean(bosq)
var(bosq)
typeof(bosq)
step_bosque!(bosq,neighbors_VN)
sum(bosq)
mean(bosq)
var(bosq)


#
# Add random position probabilities
#
function step_bosque!(landscape::Matrix{Int}, neighborhood)
    fil = size(landscape,1)
    col = size(landscape,2)
    N = fil*col
    for z in 1:N
        i= rand(1:fil) 
        j= rand(1:col)
        if landscape[i,j]==1
            x=[i,j] .+  rand(neighborhood)
            if x[1] > fil || x[1] <1 || x[2] > col || x[2]<1
                @info "Outside space $x"
            else
                landscape[x[1],x[2]] = 1
            end

        end
    end
    return
end

#
# Funcion para graficar el bosque
# 
plot_bosq(x) = heatmap(x,aspect_ratio=1,legend=:none,xticks=:none,
yticks=:none,framestyle=:none,color=[ :black , :green])



bosq = initialize_bosque()
plot_bosq(bosq)
step_bosque!(bosq,neighbors_VN)
plot_bosq(bosq)

#
# Add parameters 
#
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


bosq = initialize_bosque()
plot_bosq(bosq)
step_bosque!(bosq,neighbors_VN, (.16,.1))
plot_bosq(bosq)



bosq = initialize_bosque()
pbosq = zeros(1000)
vbosq = zeros(1000)
pbosq[1] = sum(bosq)/40000
vbosq[1] = var(bosq)/pbosq[1]
for i in 2:1000
    step_bosque!(bosq,neighbors_VN, (0.6,.1))
    pbosq[i] = sum(bosq)/40000
    vbosq[i] = var(bosq)/pbosq[i]
end
plot(pbosq, label="Mean density")
plot!(vbosq, label="Var / mean")


#
#  Ejecutar el modelo de bosque 
#
function run_bosque(vecindad, parms; pasos=100, plot::Bool=false)
	pob = zeros(pasos)
	m = initialize_bosque()
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
