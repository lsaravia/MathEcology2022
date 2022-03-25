using Plots
using DifferentialEquations

#
#  Simulación Fuego en el Bosque determinística F+V+Ø = N
#
function sim_fb_det(N₀,fin_t,h,par)
    λv , dv, λf, bf , df, N = par                       # desempaquetamos parametros 
    F, V = N₀                                        # desempaquetamos condiciones iniciales
	Fs = Float64[F]
    Vs = Float64[V]
	ts = [0.0]
	t = 0.0
	while t ≤ fin_t	
		F_next = F +  h * (-df*F + ( λf*F + bf) * V)   # Calculando el próximo valor 
        V_next = V +  h * (-(dv+bf+λf*F) * V + λv * V * (N - F - V))
		if V_next < 0                               # Si es menor que cero 
			V_next = 0                              
		end
		t_next = t + h
		push!(Fs,F_next)                            # Agregando un elemento al final del vector
		push!(Vs,V_next)            
		push!(ts,t_next)
		t = t_next
		F = F_next
        V = V_next
	end
	return ts, Vs, Fs 
end

# λv=.05; dv=10; λf=.1; bf=0.0001; df=10
# λv=.05; dv=.5; λf=.1; bf=0.0001; df=10
# λv=.05; dv=1; λf=.13; bf=1e-04; df=10
# λv=.1; dv=.01; λf=.13; bf=1e-04; df=10
# λv=5; dv=.3; λf=.2; bf=1e-03; df=90
λv=.0002; dv=.1; λf=.1; bf=0.0001; df=10

Xd = sim_fb_det([100,500],500,0.01,[λv , dv, λf, bf , df, 1000])
plot(Xd[1],Xd[2]);plot!(Xd[1],Xd[3])
plot(Xd, label=false,  xlabel="t", ylabel="V", zlabel="F")

#
#  Simulación Fuego en el Bosque estocástico
#
function sim_fb_sto(N₀,fin_t,par)
    λv , dv, λf, bf , df, N = par                       # desempaquetamos parametros 
    F, V = N₀                                        # desempaquetamos condiciones iniciales
    Ø = N - F - V
	Fs = Float64[F]
    Vs = Float64[V]
    Øs = Float64[Ø]
    ti = [0.0]              # Vector de Tiempo inicial
    time = 0.0              # suma el tiempo
	t = 2                   # indice para los vectores
    while time ≤ fin_t      # El t no es el tiempo sino la evaluacion anterior

        Bv = λv*V*Ø
        Dv = dv * V
		Bf = λf * F * V
		Bf1 = bf * V
        Df  = df * F
        R = Bv + Dv + Bf + Bf1 +Df
        y1 = rand() 
        if  y1 ≤ (Bv / R) 
            V = Vs[t-1] + 1
            Ø = Øs[t-1] - 1
		elseif y1 ≤ (Bv+Dv)/R
            V = Vs[t-1] - 1
            Ø = Øs[t-1] + 1
		elseif y1 ≤ (Bv+Dv+Bf)/R
            F = Fs[t-1] + 1
            V = Vs[t-1] - 1
        elseif y1 ≤ (Bv+Dv+Bf+Bf1)/R
            V = Vs[t-1] - 1
            F = Fs[t-1] + 1
        else 
            F = Fs[t-1] - 1
            Ø = Øs[t-1] + 1
        end
		push!(Vs,V)
    	push!(Fs,F)
		push!(Øs,Ø)

		time = ti[t-1] - log(rand()) / R 

		push!(ti,time)

        t += 1
    end     
         
    return ti , Vs, Fs, Øs
end

λv=.0002; dv=.1; λf=.1; bf=0.0001; df=10
X = sim_fb_sto([100,500],500,[λv , dv, λf, bf , df, 1000.0])
plot(X[1],X[2]);plot!(X[1],X[3])
plot(X[1:3], label=false,  xlabel="t", ylabel="V", zlabel="F")

#
#  Simulación Fuego en el Bosque estocástico Output de 0.1*t
#
function sim_fb_sto01(N₀,fin_t,par)
    λv , dv, λf, bf , df, N = par                       # desempaquetamos parametros 
    F, V = N₀                                        # desempaquetamos condiciones iniciales
    Ø = N - F - V
	Fs = Float64[F]
    Vs = Float64[V]
    Øs = Float64[Ø]
    Fn = Vn = Øn = 0.0      # Next values
    ti = [0.0]              # Vector de Tiempo inicial
    time = 0.0 
    timep = 0.0             # suma el tiempo
	t = 1                   # indice para los vectores
    while time ≤ fin_t      # El t no es el tiempo sino la evaluacion anterior

        Bv = λv*V*Ø
        Dv = dv * V
		Bf = λf * F * V
		Bf1 = bf * V
        Df  = df * F
        R = Bv + Dv + Bf + Bf1 +Df
        y1 = rand() 
        if  y1 ≤ (Bv / R) 
            Vn = V + 1
            Øn = Ø - 1
		elseif y1 ≤ (Bv+Dv)/R
            Vn = V - 1
            Øn = Ø + 1
		elseif y1 ≤ (Bv+Dv+Bf)/R
            Fn = F + 1
            Vn = V - 1
        elseif y1 ≤ (Bv+Dv+Bf+Bf1)/R
            Vn = V - 1
            Fn = F + 1
        else 
            Fn = F - 1
            Øn = Ø + 1
        end
        

		time = timep - log(rand()) / R 
        if(time > t)
            push!(Vs,V)
            push!(Fs,F)
            push!(Øs,Ø)
    		push!(ti,time)
            t += .1
        end
        F = Fn
        V = Vn
        Ø = Øn
        timep = time
    end     
         
    return ti , Vs, Fs, Øs
end

# λv=.05; dv=10; λf=.1; bf=0.0001; df=10
# λv=2; dv=1.9; λf=1; bf=0.0001; df=.01
# λv=.05; dv=1; λf=.13; bf=1e-04; df=10
# λv=5; dv=.3; λf=.2; bf=1e-03; df=90
λv=.0002; dv=.1; λf=.1; bf=0.0001; df=20

X = sim_fb_sto01([0,500],100,[λv , dv, λf, bf , df, 1000.0])
plot(X[1],X[2]);plot!(X[1],X[3])
plotly()
plot(X[1:3], label=false,  xlabel="t", ylabel="V", zlabel="F")


#
# Definicion de ecuaciones para resolverlas con el paquete
#

function fb_ode!(du, u, p, t)  # ! since modifies du -- in-place updating (can have better performance)

    
    # unpack variables and parameters:
    λv , δv, λf, bf , δf, N = p                       # desempaquetamos parametros 
    F, V = u                                        # desempaquetamos condiciones iniciales
    
    # define differential equations:
    dV = -(δv + bf+ λf * F) * V + λv * V * (N - V - F)
    dF = -δf*F + (λf * F + bf) * V
    
    du .= (dV, dF)   # copy the values into the vector du; note the `.`
end

λv=5; dv=.3; λf=.2; bf=1e-03; df=90
parameters = [λv,dv,λf,bf,df,1000]

V₀ = 100
F₀ = 1

initial_values = [F₀, V₀]

time_span = [0.0, 200.0]  # initial and final time

# set up problem:


# λv=.05; dv=10; λf=.1; bf=0.0001; df=10
# λv=2; dv=1.9; λf=1; bf=0.0001; df=.01
# λv=.05; dv=1; λf=.13; bf=1e-04; df=10
# λv=.0005; dv=.3; λf=.2; bf=1e-03; df=90
λv=.0002; dv=.1; λf=.1; bf=0.0001; df=10

function sim_fb_ode(initial_values,time_span,parameters)
    problem = ODEProblem(fb_ode!, initial_values, time_span, parameters)
    solution = solve(problem, saveat = 0.1) 
end

Xs = sim_fb_ode(initial_values,time_span,[λv,dv,λf,bf,df,1000])
plot(Xs, label=["F" "V"])


