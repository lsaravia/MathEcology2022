#= ABC con el modelo de fuegobosque2D 

* Elegir un estadístico de comparación

* Elegir una distancia

* Elegir una distribución a priori para los parámetros

* Generar datos con el modelo

* Verificar que el algoritmo ajusta a los valores usados como datos

* Ajustar con datos reales"

=#
include("fuegobosque2D.jl")

#
# ABC
#

# Definir distancia

function distance(data::Vector{Float64},simulation::Vector{Float64})
    return sqrt(sum((data.-simulation).^2.0))
end


# Estimar threshold
#
#using GLMakie
#saveanimated_fb2D(parms=[3,.87,50,3e-6,1.76],tfinal=50)

ts1 = runtimeseries_fb2D(parms=[3,.5,50,1e-5,2],tfinal=200)
ts2 = runtimeseries_fb2D(parms=[3,.5,50,1e-5,2],tfinal=200)

vecDens = ts1[1,1+end-16:end]/(40000)
vecDens2 = ts2[1,1+end-16:end]/40000
distance(vecDens,vecDens2)

lines(ts1[1,:])
lines!(ts2[1,:])
#
# ABC
#
function ABC_fb2D(data,d_λv,d_dv,d_λf, d_bf,d_df,threshold)
	
	params = zeros(Float64, 5, 100)
	successes = 0
	trials = 0 
	
	while successes < size(params, 2 )
	    trials += 1 
	    # Candidate sample
	    t_λv = rand(d_λv)
	    t_dv = rand(d_dv)
        t_λf = rand(d_λf)
        t_bf = rand(d_bf)
        t_df = rand(d_df)
	    if t_λv/t_dv > 1.5 
	        # Create timeseries
	        ts = runtimeseries_fb2D(parms=[t_λv,t_dv,t_λf, t_bf,t_df],tfinal=200)
	        last_ts = ts[1,1+end-length(data):end] / 40000
	        score = distance(data,last_ts)
	        if score < threshold
	            successes += 1
	            params[:,successes] .= [t_λv,t_dv,t_λf, t_bf,t_df]
	            @info "Tengo $successes exitos en $trials intentos! "
	        end 
	    end
	end
	return params
end

using Distributions
d_λv = Truncated(Normal(3,2),1,5)
d_dv = Truncated(Normal(.5,2),.1,2)
d_λf = Truncated(Normal(50,15),10,80)
d_bf = Truncated(Normal(1e-5,1e-4),1e-6,1e-4)
d_df = Truncated(Normal(2,1),.1,5)
plot(d_λv)

p_parms = ABC_fb2D(vecDens,d_λv,d_dv,d_λf, d_bf,d_df,0.5)
#
# Comparación priori/posteriori
#
using JLD2
@save "posteriori_ABC_fb2D" p_parms


using Plots, StatsPlots
@load "posteriori_ABC_fb2D" p_parms

plot(d_λv ,lab=" (Prior)",legend = :outertopright) 
density!(p_parms[1,:], lab= "λv (Posterior)")
xaxis!((1,5))

plot(d_dv ,lab=" (Prior)",legend = :outertopright) 
density!(p_parms[2,:], lab= "dv (Posterior)")
xaxis!((0.1,2))

plot(d_λf ,lab=" (Prior)",legend = :outertopright) 
density!(p_parms[3,:], lab= "λf (Posterior)")
xaxis!((10,80))


plot(d_bf ,lab=" (Prior)",legend = :outertopright) 
density!(p_parms[4,:], lab= "df (Posterior)")
xaxis!((1e-6,1e-4))

plot(d_df ,lab=" (Prior)",legend = :outertopright) 
density!(p_parms[5,:], lab= "λf (Posterior)")
xaxis!((.1,5))


using Statistics

median(p_parms[1,:])
median(p_parms[2,:])
median(p_parms[1,:])/median(p_parms[2,:]) 
median(p_parms[3,:])
median(p_parms[4,:])
median(p_parms[5,:])

ts = runtimeseries_fb2D(parms=[3.32,.87,50,4.5e-5,1.76],tfinal=200)

vecDens3 = ts[1+end-16:end] / 40000
plot(vecDens, label="Data")
plot!(vecDens3, label="Model")


#
# Ajustar a datos reales
#
