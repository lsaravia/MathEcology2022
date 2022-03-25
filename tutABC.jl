using StatsPlots 
using Distributions
using Random
using StatsBase

island = vec([false false false true true false true true true true false false false true true])

# 1- Colonization
# 2 - extinction
# 3 - Observation

function _transition(state::Bool; c::Float64=0.1, e::Float64=0.05)::Bool
    return state ? rand() < (1-e) : rand() < c*(1-e)
end

function islanddynamics!(timeseries::Vector{Bool}; kwargs...)
    for time in 2:length(timeseries)
        timeseries[time] = _transition(timeseries[time-1]; kwargs...)
    end
    return timeseries
end

function islanddynamics(n::Int64; kwargs...)
    timeseries = zeros(Bool, n)
    return islanddynamics!(timeseries; kwargs...)
end

function observationprocess!(timeseries::Vector{Bool}; m::Float64=0.2)
    for i in eachindex(timeseries)
        if timeseries[i]
            timeseries[i] = rand() < (1-m)
        end
    end
    return timeseries
end

function observeddynamics(n::Int64; m::Float64=0.2, kwargs...)
    timeseries = islanddynamics(n; kwargs...)
    observationprocess!(timeseries; m=m)
    return timeseries
end

[sum(observeddynamics(10)) for i in 1:1000] |>  histogram



# 
# Approximate Bayesian Computation
#

island = vec([false false false true true false true true true true false false false true true])

# Summary statistic

occupancy(timeseries::Vector{Bool}) = mean(timeseries)

transitionrate(timeseries::Vector{Bool}) = mean(timeseries[1:(end-1)] .!= timeseries[2:end])

function distance(data::Vector{Bool},simulation::Vector{Bool})
    x = [occupancy(data), transitionrate(data)]
    y = [occupancy(simulation), transitionrate(simulation)]
    return sqrt(sum((x.-y).^2.0))
end

[distance(island, observeddynamics(1000)) for i in 1:1000] |> histogram

# Cutoff to reject the sample

threshold  = 0.015

# Parameters priors

d_e = Beta(2.,5.)
plot(d_e, label="Extinction")
d_c = Beta(2., 3.) 
plot!(d_c, label="Colonization")
d_m = Truncated(InverseGaussian(3, 0.5), 0.0, 1)
support(d_m)
plot!(d_m, label="Obs Error")

#
# Pre-allocate space for the retained samples

params= zeros(Float64, 3, 100)
successes = 0
trials = 0 

while successes < size(params, 2 )
    trials += 1 
    # Candidate sample
    t_m = rand(d_m)
    t_c = rand(d_c)
    t_e = rand(d_e)
    # Create timeseries
    ts = observeddynamics(1000; m=t_m, c=t_c, e=t_e)
    score = distance(island,ts)
    if score < threshold
        successes += 1
        params[:,successes] .= [t_m, t_c, t_e]
        @info "Got $successes good samples"
    end 
end

# Comparison prior/posterior
#
plot(d_c,lab="Colonization (Prior)") 
density!(params[2,:], lab= "Colonization (Posterior)")
xaxis!((0,1))
histogram(params[2,:])

plot(d_e,lab="Extinction (Prior)") 
density!(params[3,:], lab= "Extinction (Posterior)")
xaxis!((0,1))

plot(d_m,lab="Measurement (Prior)") 
density!(params[1,:], lab= "Measurement (Posterior)")
xaxis!((0,1))

pri = plot(d_c,lab="Colonization (Prior)") 
plot!(pri,d_e,lab="Extinction (Prior)") 
plot!(pri, d_m,lab="Measurement (Prior)") 

pos = density(params[2,:], lab= "Colonization (Posterior)")
density!(pos,params[3,:], lab= "Extinction (Posterior)")
density!(pos,params[1,:], lab= "Measurement (Posterior)")

plot(pri,pos)

# Timeseries with no observation Error
#
[occupancy(islanddynamics(1000,c=params[2,i], e=params[3,i])) for i in 1:size(params,2)] |> histogram
vline!([occupancy(island)],lab="Measurement")


[transitionrate(islanddynamics(1000,c=params[2,i], e=params[3,i])) for i in 1:size(params,2)] |> histogram
vline!([transitionrate(island)],lab="Measurement")

# Parameter relationship

scatter(params[2,:], params[3,:])
xaxis!("Colonization", (0,1))
yaxis!("Extinction", (0,1))