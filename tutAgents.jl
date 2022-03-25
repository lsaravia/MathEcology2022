#
#
# 

using Agents

space= GridSpace((10,10))

mutable struct Schelling <: AbstractAgent
    id::Int
    pos::Tuple{Int,Int}
    group::Int
    happy::Bool
end

properties = Dict(:min_to_be_happy => 3)

model = AgentBasedModel(Schelling, space; properties)

scheduler = Schedulers.randomly 

model = AgentBasedModel(Schelling,space; properties,scheduler)

function initialize(; N=320,M=20, min_to_be_happy = 3)

    space= GridSpace((M,M))
    scheduler = Schedulers.randomly 
    properties = Dict(:min_to_be_happy => min_to_be_happy)
    model = AgentBasedModel(Schelling,space; properties,scheduler)

    # add Agents
    for n in 1:N
        agent = Schelling(n, (1,1), n < N/2 ? 1 : 2, false)
        add_agent_single!(agent,model)
    end
    return model
end

model = initialize()

function agent_step!(agent,model)
    agent.happy && return
    nearby_same = 0
    for neighbor in nearby_agents(agent,model)
        if agent.group == neighbor.group
            nearby_same += 1
        end
    end
    if nearby_same ≥ model.min_to_be_happy
        agent.happy = true
    else
        move_agent_single!(agent,model)
    end
    return
end

step!(model,agent_step!)

using InteractiveDynamics, GLMakie

fig, _ = abm_plot(model)

display(fig)

groupcolor(agent) = agent.group == 1 ? :blue : :orange
groupmarker(agent) = agent.group == 1 ? :circle : :rect

fig, _ = abm_plot(model, ac=groupcolor, am=groupmarker)

model = initialize()
abm_play( 
    model, agent_step!,
    ac=groupcolor, am=groupmarker, as=12

)

#
# Collect data
#
z(agent) = agent.pos[1]
adata = [:group, :happy, z]
model = initialize()
adf, _ = run!(model,agent_step!,dummystep,3; adata)

using Statistics: mean
adata = [
    (:group, sum),
    (z, mean)]

model = initialize()
adf, _ = run!(model,agent_step!,dummystep,3; adata)
    