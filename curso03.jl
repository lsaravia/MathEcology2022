#
#  Simulación SIS determinística S+I = 1
#
function sim_sis(N₀,fin_t,h,par)
	μ , α = par      
	Is = Float64[N₀]
	I = Float64(N₀)
	ts = [0.0]
	t = 0.0
	while t ≤ fin_t	
		I_next = I +  h * (-α*I + μ* I *( 1 - I ))  # Calculando el próximo valor 
		if I_next < 0                               # Si es menor que cero 
			I_next = 0                              
		end
		t_next = t + h
		push!(Is,I_next)           # Agregando un elemento al final del vector Ns
		push!(ts,t_next)
		t = t_next
		I = I_next
	end
	return ts, Is
end

α = 1; μ=1
Xd = sim_sis(.4,100,0.01,[μ,α])
plot(Xd); hline!([1-α/μ], c=:grey, ls=:dash)


#
#  Simulación SIS estocastica con 4 eventos
#
function sim_sis_sto(N₀,fin_t,par)
	μ , α, N  = par          # desempacar los parametros 

    It = [N₀]               # Serie de tiempo de la poblacion
    I  =  N₀
	S  =  N - I 
	St = [S]
    ti = [0.0]              # Vector de Tiempo inicial
    time = 0.0              # suma el tiempo
	t = 2                   # indice para los vectores
    while time ≤ fin_t      # El t no es el tiempo sino la evaluacion anterior

        
        Bi = μ * S * I  
        Di = α * I
		Bs = α * I
		Ds = μ * S * I  
        R = Bi + Di + Bs +Ds
        y1 = rand() 
        if  y1 < (Bi / R) 
            I = It[t-1] + 1
		elseif y1 < (Bi+Di)/R
            I = It[t-1] - 1
		elseif y1 < (Bi+Di+Bs)/R
            S = St[t-1] + 1
		else
			S = St[t-1] - 1
        end
    	push!(It,I)
		push!(St,S)
		time = ti[t-1] - log(rand()) / R 

		push!(ti,time)

        t += 1
    end     
         
    return ti , It, St
end


par = [0.1,1, 100];
X = sim_sis_sto(50,100,par)
using Plots
plot(X[1],X[2]);plot!(X[1],X[3])
plot(X, label=false)

#
#  Simulación SIS estocastica con 2 eventos
#
function sim_sis_sto01(N₀,fin_t,par)
	μ , α, N  = par          # desempacar los parametros 

    It = [N₀]               # Serie de tiempo de la poblacion
    I  =  N₀
	S  =  N - I 
	St = [S]
    ti = [0.0]              # Vector de Tiempo inicial
    time = 0.0              # suma el tiempo
	t = 2                   # indice para los vectores
    while time ≤ fin_t      # El t no es el tiempo sino la evaluacion anterior

        
        B = μ * S * I  
        D = α * I
        R = B +D
        y1 = rand() 
        if  y1 ≤ (B / R) 
            I = It[t-1] + 1
			S = St[t-1] - 1
		else
            I = It[t-1] - 1
            S = St[t-1] + 1
        end
    	push!(It,I)
		push!(St,S)
		time = ti[t-1] - log(rand()) / R 

		push!(ti,time)

        t += 1
    end     
         
    return ti , It, St
end

parz = [.001,.5,100]
Z = sim_sis_sto01(50,100,parz)
plot(Z[1],Z[2]);plot!(Z[1],Z[3])
hline!([parz[3]-parz[2]/parz[1]],c=:grey, ls=:dash)
