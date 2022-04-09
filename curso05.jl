#
# Modelo estocastico espacial de una población
#
using Plots
using Statistics
using Distributions

#
# Funcion para graficar el bosque
# 
plot_bosque(x) = heatmap(x,aspect_ratio=1,legend=:none,xticks=:none,
yticks=:none,framestyle=:none,color=[ :black , :red])


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


#
# Prueba de funciones
#
bosq = initialize_bosque()
plot_bosque(bosq)
step_bosque!(bosq,neighbors_VN, (.16,.1))
plot_bosque(bosq)

#
# Version mejorada con densidad inicial mas exacta
#
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
# Ejemplo de dinámica
#
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
function run_bosque(neighborhood, parms; pasos=100, plot::Bool=false)
	pob = zeros(pasos,2) 
	m = initialize_bosque()
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

rr = run_bosque(neighbors_VN,(1.5,1))
plot(rr, label=["Densidad" "CV Espacial"])

#
# generating a gif
#
using Statistics
mm = initialize_bosque(fil=200,col=200, density=0.1 )
fint = 100
N = 200*200
pob = [sum(mm)/N]
λ = 1.7
δ = 1
@gif for i in 1:fint
    step_bosque!(mm,neighbors_VN,(λ,δ));

    plt1 = plot_bosque(mm)
    push!(pob,sum(mm)/N)         # densidad en el espacio 
    #pob[i,2]=var(mm)/pob[i,1]   # Coeficiente de variacion = Varianza/media
    plt2=plot(pob,title="λ=$(λ) δ=$(δ) μ=$(λ/δ)", labels="Density",xlim=(0,fint),ylim=(0,.2))
    plot(plt1, plt2, layout = (1, 2))
end