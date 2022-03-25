
#
# Crecimiento 
#
function sim_exponencial( N₀,fin_t,h,λ )    # N\_0 TAB
    xs = zeros(fin_t)
    xs[1]=N₀
    ti = zeros(fin_t)
    for t in 2:fin_t
        xs[t] = xs[t-1] *(1 + h*λ)
        ti[t] = ti[t-1] + h 
    end
        
    return ti , xs
end

#
# Esta funcion no simula hasta fin_t
#

using Plots
pop = sim_exponencial(10,1000,0.1,.1);

plot(pop, label="Determ.")


function sto_crecimiento(N₀,fin_t,h,λ)
    xs = zeros(fin_t)
    xs[1]=N₀
    ti = zeros(fin_t)
    for t in 2:fin_t
        if rand() < h*λ 
            xs[t] = xs[t-1] + 1
        else
            xs[t] = xs[t-1] 
        end
        ti[t] = ti[t-1] + h 
    end
        
    return ti , xs
end
pop1 = sto_crecimiento(10,1000,0.1,.1);
plot!(pop1, yscale= :log10, l=:scatter, label="Estoc.")

#
# No anda porque Hay que hacerlo xs[t-1] veces
# 

scatter(pop1)


function stochastic_exacto(N₀,fin_t,λ)
    xs = zeros(fin_t)
    xs[1]=N₀
    ti = zeros(fin_t)

    for t in 2:fin_t
        
        ti[t] = ti[t-1] - log(rand()) / (xs[t-1] * λ)
        xs[t] = xs[t-1] + 1
        
    end
        
    return ti , xs
end
pop2 = stochastic_exacto(10,1000,.1)
pop1 = sim_crecimiento(10,500,.1,.1);
scatter(pop2, label= "Stoc." )
scatter!(pop1, yaxis= :log, l= :line, label= "Determ.")

#
# Funcion que calcula la solucion de la ecuación diferencial
#
crec_exp(N₀,λ,t) = N₀ *exp(λ*t)

pop3 = [crec_exp(10,.1,t) for t in 1:50]

scatter!(pop3, l= :line)

#
# Ejercicio 1: modificar la funcion sim_crecimiento
# para que simule hasta fin_t
#
# Ejercicio resuelto 
#

function sim_crecimiento( N₀,fin_t,h,λ )
    stp = Int64(fin_t / h)
    xs = zeros(stp)
    xs[1]=N₀
    ti = zeros(stp)
    for t in 2:stp
        xs[t] = xs[t-1] *(1 + h*λ)
        ti[t] = ti[t-1] + h 
    end
        
    return ti , xs
end

pop2 = stochastic_crecimiento(10,1000,.1)
pop1 = sim_crecimiento(10,50,.1,.1);
plot(pop2 )
plot!(pop1, yaxis= :log)

#
# Ejercicio 1: Modificar la funcion sim_crecimiento
#              para que simule hasta fin_t
#
# Ejercicio 2: modificar la funcion stochastic_crecimiento
#              para que simule hasta fin_t
#
# Ejercicio 3: Modificar la funcion sto_crecimiento para que simule
#              bien el crecimiento poblacional
#