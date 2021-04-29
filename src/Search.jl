"""
    Search

Provides the types for objects and scenario management. This module also includes plot recipes for the various object types.

# Types
- `SimObject`: the abstract type that defines a basic interface for all objects within the sim environment.
    - `Target`: the subtype of `SimObject` that specializes for the search targets.
    - `Agent`: the subtype of `SimObject` that specializes for the searching agent.
    - `Sensor`: the subtype of `SimObject` that specializes for search sensors.
- `Scenario`: the type that contains an `Agent`, `Target`, `Sensors`, and other parameters.
"""
module Search

# using Distributions
# using StatsPlots
using Interpolations
using RecipesBase
using Plots

export Target, Agent, Sensor, Scenario, run!

"""
    SimObject

Abstract SimObject, which requires the following interface:
# Fields
- `f::Function`: the policy update function.
- `x::Array`: the position history (dimension N x T where N = length(position) and T = length(t)
- `v::Array`: the velocity history (dimension N x T where N = length(position) and T = length(t)
- `t::Array`: the time history (dimension T x 1)
"""
abstract type SimObject end

"""
    Sensor(f::Function, x::Array, v::Array, t::Array, r::Number)

Subtype of SimObject. It adds the field `r`.
# Fields
- `f::Function`: the policy update function.
- `x::Array`: the position history (dimension N x T where N = length(position) and T = length(t)
- `v::Array`: the velocity history (dimension N x T where N = length(position) and T = length(t)
- `t::Array`: the time history (dimension T x 1)
- `r::Number`: the detection radius around the sensor.
"""
mutable struct Sensor <: SimObject
    f::Function
    x::Array
    v::Array
    t::Array
    r::Real
end

"""
    Agent(f::Function, x::Array, v::Array, t::Array, slims::Tuple, s::Array)

Subtype of SimObject. It adds the field `slims` and `s`.
# Fields
- `f::Function`: the policy update function.
- `x::Array`: the position history (dimension N x T where N = length(position) and T = length(t)
- `v::Array`: the velocity history (dimension N x T where N = length(position) and T = length(t)
- `t::Array`: the time history (dimension T x 1)
- `slims::Tuple`: the upper and lower limits on the Agent speed.
- `s::Array`: the array of sensors available to the Agent.
"""
mutable struct Agent <: SimObject
    f::Function
    x::Array
    v::Array
    t::Array
    slims::Tuple
    s::Array
end

"""
    Target(f::Function, x::Array, v::Array, t::Array, slims::Tuple)

Subtype of SimObject. It adds the field `slims`.
# Fields
- `f::Function`: the policy update function.
- `x::Array`: the position history (dimension N x T where N = length(position) and T = length(t).
- `v::Array`: the velocity history (dimension N x T where N = length(position) and T = length(t).
- `t::Array`: the time history (dimension T x 1).
- `slims::Tuple`: the upper and lower speed limit of the target.
"""
mutable struct Target <: SimObject
    f::Function
    x::Array
    v::Array
    t::Array
    slims::Tuple
end

"""
    Scenario(a,t,s)

Type to hold the scenario parameters and objects.
# Fields
- `agent::Agent`
- `target::Target`
- `sensors::Array`
- `t::Tuple`
- `Δt::Real`
"""
mutable struct Scenario
    agent::Agent
    target::Target
    sensors::Array
    t::Tuple
    Δt::Real
end

"`Sensor(x,r)`: Initialize a Sensor with the default policy `policyDoNothing!` at time zero."
Sensor(x,r) = Sensor(policyDoNothing!, x, r)

"`Sensor(x,t,r)`: Initialize a Sensor with the default policy `policyDoNothing!` at time `t`."
Sensor(x,t,r) = Sensor(policyDoNothing!, x, [t], r)

"`Sensor(f::Function,x,r)`: Initialize a Sensor with zero initial velocity at time zero."
Sensor(f::Function,x,r) = Sensor(f, x, x.*0, [0.0], r)

"`Sensor(f::Function,x,t,r)`: Initialize a Sensor with zero initial velocity at time `t`."
Sensor(f::Function,x,t,r) = Sensor(f, x, x.*0, [t], r)


"Initialize an Agent with zero initial velocity and default speed limits."
Agent(f::Function,x, sensors; slims = (0,120)) = Agent(f, x, x.*0, [0.0], slims, sensors)

"`Target(f::Function,x,v,slims)`: Initialize a Target at time zero."
Target(f::Function,x,v,slims) = length(x)!=length(v) ? error("Position and Speed have unequal dimension.") : Target(f,x,v,[0.0],slims)
Target(x,v,slims) = Target(policyDoNothing!, x, v, slims)
Target(x; slims = (2,6)) = Target(x[1:length(x)/2], x[length(x)/2+1:end], slims)
"""
    Target(;slims, θlims)

Randomly initalize a Target at time zero with speed and direction drawn from uniform distribution over `slims` and `θlims`
"""
function Target(;slims=(2,6), θlims = (0,2π))
    x = [0.0, 0.0]
    θ = round(rand(range(θlims[1], stop = θlims[end], length = 100)), digits = 2)
    speed = round(rand(range(slims[1],stop = slims[end], length = 100)), digits = 2)
    v = speed .* [cos(θ),sin(θ)]
    Target(x,v,slims)
end

"""
    Scenario(::Agent, ::Target, ::Array, t, Δt)

Contain the objects within a scenario.
"""
function Scenario(t,Δt)
    a = Agent
end

"""
    step!(s::SimObject, Δt)

Propogate forward the position of a `SimObject` based on the velocity and time.
```
    x = x + v⋅Δt
    v = v
    t = t + Δt    
```
"""
function step!(s::SimObject, Δt)
    s.t = vcat(s.t, s.t[end]+Δt)
    s.x = hcat(s.x, s.x[:,end] + s.v[:,end] .* Δt)
    s.v = hcat(s.v, s.v[:,end])
    return s
end

"`step!(s::Sensor,Δt)`: A specialized step! for `Sensor`, as they do not move."
function step!(s::Sensor, Δt)
    s.t = vcat(s.t, s.t[end]+Δt)
    s.x = hcat(s.x, s.x[:,end])
    s.v = hcat(s.v, s.v[:,end])
end

"""
    policyDoNothing!(s::SimObject)

The default policy to make no changes from the initial conditions.
"""
function policyDoNothing!(s::SimObject)
    return
end

function policySpiral!(a::Agent, Δt)
    a.v[:,end] = sum(a.slims)/4 * [cos(a.t[end]*2π), sin(a.t[end]*2π)]
end

"""
    policyTurnPeriodically!(t::Target)

Turn at regular intervals. Parameters are defined within the function body.
"""
function policyTurnPeriodically!(t::SimObject, Δt)
    if mod(t.t[end],.5)==0
        s = (2,6)
        Δθ = .5
        θ = atan(t.v[1,end],t.v[2,end])
        θ += 2*(rand()-.5)*Δθ
        speed = round(rand(range(s[1],stop = s[end], length = 100)), digits = 2)
        t.v = speed .* [cos(θ),sin(θ)]
    end
end

"`checkContact(t::Target,s::Sensor)::Bool`: Return `true` if the distance between `t` and `s` is less than or equal to detection range `s.r`, otherwise `false`."
checkContact(t::Target,s::Sensor)::Bool = distance(t,s) <= s.r ? true : false

"""
    distance(s1::SimObject, s2::SimObject)::Real

Return the distance between two `SimObject`'s.
# Input
- `s1::SimObject`
- `s2::SimObject`

# Output
The 2-norm distance between s1 and s2.
"""
function distance(s1::SimObject, s2::SimObject)::Real
    sqrt(sum((s1.x - s2.x).^2))
end

"""
    plot(t::Agent; interval = .25)

Plot the position history of an Agent.
"""
@recipe function plot(t::Agent; interval = .25)
    tick = t.t[1]:interval:t.t[end]
    tockx = LinearInterpolation(t.t,t.x[1,:])
    tocky = LinearInterpolation(t.t,t.x[2,:])
    @series begin
        linewidth := 0
        markershape --> :+
        tockx.(tick), tocky.(tick)
    end 
    xguide --> "x (nm)"
    yguide --> "y (nm)"
    # title --> $(t.t[end])
    aspect_ratio --> :equal
    linecolor --> cgrad(:winter,rev=true)
    # colorbar --> true
    linewidth --> 2
    line_z --> t.t
    legend := false
    background_color_inside --> :lightblue
    n = length(t.t)
    seriesalpha --> (range(0, 1, length = n)).^2
    markershape --> :square
    markeralpha --> ifelse.(t.t .≈ 0, 1, 0)
    markercolor --> :red
    @series t.x[1,:], t.x[2,:]

end

"""
    plot(t::Target; interval = .25)

Plot the position history of a Target.
"""
@recipe function plot(t::Target;interval=.25)
    tick = t.t[1]:interval:t.t[end]
    tockx = LinearInterpolation(t.t,t.x[1,:])
    tocky = LinearInterpolation(t.t,t.x[2,:])
    @series begin
        linewidth := 0
        markershape --> :+
        tockx.(tick), tocky.(tick)
    end 
    xguide --> "x (nm)"
    yguide --> "y (nm)"
    # title --> $(t.t[end])
    aspect_ratio --> :equal
    linecolor --> cgrad(:CMRmap,rev=true)
    # colorbar --> true
    linewidth --> 2
    line_z --> t.t
    legend := false
    background_color_inside --> :lightblue
    # n = length(t.t)
    # seriesalpha --> (range(0, 1, length = n))
    markershape --> :diamond
    markeralpha --> ifelse.(t.t .≈ 0, 1, 0)
    markercolor --> :yellow
    t.x[1,:], t.x[2,:]
end


"""
    run!(s::Scenario)

Loop through the Scenario time steps, apply policy and then step for each object
"""
function run!(s::Scenario)
    for t in s.t[1]:s.Δt:s.t[end]
        for obj in vcat(s.sensors, s.target, s.agent)
            policy!(obj)
        end
        for obj in vcat(s.sensors, s.target, s.agent)
            step!(obj)
        end
    end
end


# mutable struct Grid
#     x::AbstractVector
#     y::AbstractVector
#     v::AbstractMatrix
# end

# function Grid(x::AbstractVector, y::AbstractVector, d::Sampleable{Multivariate, Continuous})
#     v = [pdf(d,[i,j]) for i in x, j in y]
#     Grid(x,y,v)
# end

# mutable struct SearchScenario
#     initialDistribution::Sampleable{Multivariate, Continuous}
#     searchRadius::Number
#     searchSchedule::Array
#     AOU::Grid
# end

# function clip!(g::Grid,ϵ=1e-9)
#     g.v[abs.(g.v) ≈ 0] .= 0
# end

# function normalize!(g::Grid)
#     g.v = g.v ./ sum(g.v)
# end

# function clipNormalize!(g::Grid)
#     clip!(g)
#     normalize!(g)
# end

# import Base.+
# function +(a::Grid,b::Grid)
#     if a.x != b.x
#         error("x-axis does not match")
#     end
#     if a.y != b.x
#         error("y-axis does not match")
#     end
#     Grid(a.x, a.y, a.v .+ b.v)
# end

# import Base.*
# function *(a::Grid,b::Grid)
#     if a.x != b.x
#         error("x-axis does not match")
#     end
#     if a.y != b.x
#         error("y-axis does not match")
#     end
#     Grid(a.x, a.y, a.v .* b.v)
# end

# function initializeScenario(x::AbstractVector, y::AbstractVector, initDist::Sampleable{Multivariate, Continuous},searchRadius,searchSched::Array)::SearchScenario
#     SearchScenario(initDist,searchRadius,searchSched,Grid(x,y,initDist));
# end

# function initializeScenario(initDist::Sampleable{Multivariate, Continuous},searchRadius,searchSched::Array,model::Function;N = 1000)::SearchScenario
#     targets = Vector{Target}(undef,N)
#     for i in 1:N
#         targets[i] = Target(rand(initDist));
#     end
#     scenario = SearchScenario(initDist,searchDist,searchSched,model,targets);
# end

# function step!(scenario::SearchScenario)
#     map(scenario.targetModel!,scenario.targets)
#     dots = map(t -> t.x,targets)
#     p = scatter(dots[:,1],dots[:,2])
#     display(p)
# end

# function isFound(target::Target,search::MvNormal)
#     norm(target.x - search.μ) <= search.Σ[1]*1.177 ? true : false
# end


end # module
