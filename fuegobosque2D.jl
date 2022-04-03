#
# Fuego en el bosque 2D
#

#
# Structura de datos para el modelo con el Array de posiciones y las dimensiones
#
struct fb_model
    landscape::Matrix{Int8}
    fil::Int64
    col::Int64
    function fb_model(x::Matrix{Int8})
        f = size(x,1)
        c = size(x,2)
        new(x,f,c)
    end 
end

"""
    initialize_fb2D(;fil::Int=200,col::Int=200, density::Float64=0.1)

Inicializa el modelo de Fuego en el bosque con dimensiones `fil x col` con una densidad `density` de 1

"""
function initialize_fb2D(;fil::Int=200,col::Int=200, density::Float64=0.1)
    #space = zeros(Int8,fil,col)
    #space = fill(:Ø,fil,col)
    mdl = fb_model(zeros(Int8,fil,col))
    ini_dens = trunc(Int, mdl.fil*mdl.col*density)
    for i in eachindex(mdl.landscape)
        if rand() ≤ density        
            mdl.landscape[i] = 1 # :V
        end
    end
    return mdl
end


"""
    step_fb2D!(landscape::Matrix{Symbol}, neighborhood, parms)
    
Ejecuta un paso del modelo de Fuego en el bosque identico a [1]. 
Los parametros son λv tasa de propagación de la vegetación (Bosque)
δv tasa de mortalidad de la vegetación, λf tasa de propagación del fuego
bf  tasa de ignicion espontanea del fuego, δf tasa de extinción del fuego

1. Nicoletti, G., Saravia, L., Momo, F., Maritan, A., and Suweis, S. (2021). 
The emergence of scale-free fires in Australia. arXiv:2110.10014 [cond-mat]. 
Available at: http://arxiv.org/abs/2110.10014 [Accessed October 20, 2021].


"""
function step_fb2D!(mdl::fb_model, vecindad, parms)

    # Size of the grid defined in fb_model

    # Parametros
    λv, δv,  λf,  bf,  δf = parms

    # Total rate of events
    R = λv + δv + λf + bf + δf

    # Number of iterations  to complete one time step
    evals = mdl.fil*mdl.col*R

    
    # Longitud total del Array
    l = length(mdl.landscape)

    for i in 1:evals
        i= rand(1:mdl.fil) 
        j= rand(1:mdl.col)
        # rpos = rand(1:l)   # posicion al azar dentro del Array
        y1 = rand()
        if y1 ≤ λv/R 
            #
            # toma un vecino y reproduce un V
            #
            if mdl.landscape[i,j]== 1    # :V
                x=(i,j) .+  randvecino_fb2D(neighbors_VN)
                if !fueradeborde_fb2D(mdl,x)
                    mdl.landscape[x[1],x[2]] = 1 # :V
                end
            end
        elseif y1 ≤ (λv + δv) / R
            #
            # Muere el V actual
            #
            if mdl.landscape[i,j] == 1    # :V
                mdl.landscape[i,j] = 0   #   :Ø
            end
        elseif y1 ≤ (λv + δv + λf) / R
            #
            # Se propaga F
            #
            if mdl.landscape[i,j] == 2    # :F
                x=(i,j) .+  randvecino_fb2D(neighbors_VN)
                if !fueradeborde_fb2D(mdl,x)
                    if mdl.landscape[x[1],x[2]] == 1    # :V
                        mdl.landscape[x[1],x[2]] = 2    # :F
                    end
                end
            end
        elseif y1 ≤ (λv + δv + λf + bf) / R 
            #
            # Se produce una ignición
            # 
            if mdl.landscape[i,j] == 1    # :V
                mdl.landscape[i,j] = 2    # :F
            end
        else
            #
            # Se extingue el F:
            #
            if mdl.landscape[i,j] == 2    # :F
                mdl.landscape[i,j] = 0   #   :Ø
            end
        end
    end
end

#
# Vecindad de Von Neuman = 4 vecinos mas cercanos
#
neighbors_VN = ((0,1),(0,-1),(1,0),(-1,0) )

"""
    randvecino_fb2D(vecindad::Tuple)

Una vecindad aleatoria en coordenadas 2D basado en una Tuple con coordenadas
para otros tipos de vecindades debería tomar una función
"""
function randvecino_fb2D(vecindad::Tuple)
    rand(vecindad )    
end

function fueradeborde_fb2D(mdl::fb_model,x::Tuple{Int64,Int64})
    if x[1] > mdl.fil || x[1] <1 || x[2] > mdl.col || x[2]<1
        #@info "Outside space $x"
        return true
    else
        return false
    end
end

#
#  Plot 2D Animated (Falta interactive)
#
function runanimated_fb2D(;fil=200,col=200,parms=[3,.5,500,1e-5,20],tfinal=100)
    λv, dv, λf, bf, df = parms

    fbm = initialize_fb2D(fil=fil,col=col,density=0.1)
    obs = Observable(fbm.landscape)
    t = 0 
    tobs = Observable("T=$t")
    display(Makie.heatmap(obs,  axis = (aspect = 1,title=tobs),colorrange=(0, 2), colormap=[:black, :green, :red]))
    tsF = zeros(2,tfinal)

    for i in 1:tfinal
        step_fb2D!(fbm, neighbors_VN, (λv, dv, λf, bf, df) )
        obs[] = fbm.landscape
        tsF[1,i]=count(k->k==1,fbm.landscape) 
        tsF[2,i]=count(k->k==2,fbm.landscape) 
        tobs[] = "T=$i V=$(tsF[1,i]) F=$(tsF[2,i])"
        sleep(0.001)
    end
end


#
#  Plot 2D Animated (Falta interactive)
#
function saveanimated_fb2D(;fil=200,col=200,parms=[3,.5,500,1e-5,20],tfinal=100,fname="fb2Danimation.mp4")
    λv, dv, λf, bf, df = parms

    fbm = initialize_fb2D(fil=fil,col=col,density=0.1)
    obs = Observable(fbm.landscape)
    t = 0 
    tobs = Observable("T=$t")
    fig = Makie.heatmap(obs,  axis = (aspect = 1,title=tobs),colorrange=(0, 2), colormap=[:black, :green, :red])
    tsF = zeros(2,tfinal)

    record(fig, fname, 1:tfinal, framerate=24) do frame 
        step_fb2D!(fbm, neighbors_VN, (λv, dv, λf, bf, df) )
        obs[] = fbm.landscape
        t += 1
        tsF[1,t]=count(k->k==1,fbm.landscape) / (fil*col)
        tsF[2,t]=count(k->k==2,fbm.landscape) / (fil*col)
        
        tobs[] = "T=$t V=$(round(tsF[1,t],digits=2)) F=$(round(tsF[2,t],digits=2))"
        #sleep(0.001)
    end
end

# runanimated_fb2D(parms=[3,.5,50,1e-5,2],tfinal=50)

#
# Plot time series
#
"""
    runtimeseries_fb2D(;fil=200,col=200,parms=[3,.5,50,1e-5,2],tfinal=100,plot=false)

Corre el modelo de Fuego en el bosque 2D con el tamaño de grilla `fil,col` los parametros `parms = (λv, dv, λf, bf, df)` y el tiempo total `tfinal`
Si `plot=true` grafica la serie usando GLMakie.

"""
function runtimeseries_fb2D(;fil=200,col=200,parms=[3,.5,50,1e-5,2],tfinal=100,plot=false)
    λv, dv, λf, bf, df = parms

    fbm = initialize_fb2D()
    
    tsF = zeros(2,tfinal)
    for i in 1:tfinal
        step_fb2D!(fbm, neighbors_VN, (λv, dv, λf, bf, df) )
        tsF[1,i]=count(k->k==1,fbm.landscape) 
        tsF[2,i]=count(k->k==2,fbm.landscape) 
    end
    if plot
        display(lines(tsF[1,:]/(fil*col), color=:green))
        lines!(tsF[2,:]/(fil*col), color=:red)
    end
    return tsF
end

# runtimeseries_fb2D(parms=[3,.5,50,1e-5,2],tfinal=500)
