#
# Simulación crecimiento muerte lineal
#
function sim_crec_muert(N₀,fin_t,h,par)
	μ,δ = par
	Nt = [N₀]
	N  = N₀
	st = [0.0]
	t  = 0.0
	while t ≤ fin_t	
		N_next =  N + h* (μ - δ* N)
		t_next = t + h
		push!(Nt,N_next)
		push!(st,t_next)
		N = N_next
		t = t_next
	end
	return st, Nt
end


parameters = [1.8, .1]
N₀ =50
X = sim_crec_muert(N₀,100,0.1,parameters)

using Plots
plot(X)

X1 = [sim_crec_muert(i,100,0.1,parameters) for i in range(10.0,100.0,5)]

plot(X1,label=collect(range(10.0,100.0,5))')    # transpose( data ) = data'

#
# Simulación crecimiento muerte estocástica
#
function sim_crec_muert_sto(N₀,fin_t,par)
	μ , δ = par
	Ns = Float64[N₀]
	N = Float64(N₀)
	ts = [0.0]
	time = 0.0
	t = 2
	while time ≤ fin_t
		B = μ                          # Tasa de natalidad 
		D = δ * N
		R = B + D 
		y1 = rand()
		if y1 < B/R 
			N = Ns[t-1] + 1  # Calculando el próximo valor de N 
		else 
			N = Ns[t-1] - 1
		end
		push!(Ns,N)         # Agregando un elemento al final del vector Ns
		time = ts[t-1] - log(rand()) / R
		push!(ts,time)
		t += 1 
	end
	return ts, Ns
end

par = [1.8, 0.1]
X1 = [ sim_crec_muert_sto(i, 100, par) for i in range(0,100,10)]
plot(X1, label=collect(range(0,100,10))')

#
# En este caso está bien que empezando de cero vaya al atractor? 
# 
# Es cero un atractor?"
#

#
# Simulación logística determinista
#
function sim_logistico( N₀,fin_t,h, par )
    r , K = par                             # desempacar los parametros 
    stp = Int64(fin_t / h)                  # dividimos el tiempo de simulación por el intervalo de integración h
    xs = zeros(stp)
    xs[1]=N₀
    ti = zeros(stp)
    for t in 2:stp
        xs[t] = xs[t-1] + h * xs[t-1] * r *(1 - xs[t-1]/K )
        ti[t] = ti[t-1] + h 
    end
        
    return ti , xs
end

#
#  Simulación logistica estocastica exacta con fin_t pasos
#
function sim_logistico_sto1(N₀,fin_t,par)

    r , K = par           # desempacar los parametros 

    xs = [N₀]             # Serie de tiempo de la poblacion
    x = N₀
    ti = [0]              # Vector de Tiempo inicial

    for t in 2:fin_t      # El t no es el tiempo sino la evaluacion anterior
        B = r* xs[t-1]
        D = r / K * xs[t-1]^2
        R = B + R
        y1 = rand() 
        if  y1 < (B / R) 
            xs[t] = xs[t-1] + 1
        else
            xs[t] = xs[t-1] - 1
        end
        ti[t] = ti[t-1] - log(rand()) / R 
        t += 1
    end     
         
    return ti , xs
end


#
#  Simulación logistica estocastica exacta con fin_t pasos
#
function sim_logistico_sto(N₀,fin_t,par)

    r , K = par           # desempacar los parametros 

    xs = [N₀]             # Serie de tiempo de la poblacion
    x = N₀
    ti = [0]              # Vector de Tiempo inicial

    for t in 2:fin_t      # El t no es el tiempo sino la evaluacion anterior
        B = r* xs[t-1]
        D = r / K * xs[t-1]^2
        R = B + R
        y1 = rand() 
        if  y1 < (B / R) 
            x = xs[t-1] + 1
        else
            x = xs[t-1] - 1
        end
        t += 1
        push!(xs,x)
		time = ti[t-1] - log(rand()) / R 
		push!(ti,time)

    end     
         
    return ti , xs
end



#
#  Simulación logistica estocastica exacta que termina en tiempo fin_t
#
function sim_logistico_sto(N₀,fin_t,par)

    r , K = par           # desempacar los parametros 

    xs = [N₀]             # Serie de tiempo de la poblacion
    x = N₀
    ti = [0.0]              # Vector de Tiempo inicial
    time = 0.0              # suma el tiempo
	t = 2                   # indice para los vectores
    while time ≤ fin_t      # El t no es el tiempo sino la evaluacion anterior
        B = r* xs[t-1]
        D = r / K * xs[t-1]^2
        R = B + D
        y1 = rand() 
        if  y1 < (B / R) 
            x = xs[t-1] + 1
        else
            x = xs[t-1] - 1
        end
    	push!(xs,x)
		time = ti[t-1] - log(rand()) / R 

		push!(ti,time)

        t += 1
    end     
         
    return ti , xs
end
