struct Integrator{T<:Real} <: Simulation
    scheme::Function
    h::T
    h2::T
    tmax::T

    function Integrator(scheme::Function,h::T,tmax::T) where T <: Real
        h2::T = 0.5*h
        return new{T}(scheme,h,h2,tmax)
    end
end

function (intr::Integrator)(sim::State{T}) where T <: Real
    for i in 1:intr.tmax
        zero_out!(sim.d)
        intr.scheme(sim,intr.h,intr.h2)
    end
    return
end

function (intr::Integrator)(sim::State{T},steps::T) where T <: Real
    for i in 1:steps
        zero_out!(sim.d)
        intr.scheme(sim,intr.h,intr.h2)
    end
    return
end
#=
function (integrator::Integrator)(sim::State{T},::Type{TTV})
    # DO TTV_ELEMENTS + TTV
end =#
