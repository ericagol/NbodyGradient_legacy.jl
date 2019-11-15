mutable struct DriftArrays{T<:Real} <: PreAllocArrays
    state0::Array{T,1}
    state::Array{T,1}
    x0::Array{T,1}
    v0::Array{T,1}

    function DriftArrays(::T) where T <: Real
        state0 = zeros(T,12)
        state = zeros(T,12)
        x0 = zeros(T,3)
        v0 = zeros(T,3)
        return new{T}(state0,state,x0,v0)
    end
end

mutable struct Derivatives{T<:Real} <: PreAllocArrays
    jac_phi::Array{T,2}
    jac_kick::Array{T,2}
    jac_copy::Array{T,2}
    jac_ij::Array{T,2}
    jac_tmp1::Array{T,2}
    jac_tmp2::Array{T,2}
    dqdt_phi::Array{T,1}
    dqdt_kick::Array{T,1}
    dqdt_ij::Array{T,1}
    jac_kepler::Array{T,2}

    function Derivatives(init::InitialConditions)
        T = typeof(init.t0)
        n = 7*init.nbody
        jac_phi = Matrix{T}(I,n,n)
        jac_kick = zeros(T,n,n)
        jac_copy = zeros(T,n,n)
        jac_ij = zeros(T,14,14)
        jac_tmp1 = zeros(T,14,n)
        jac_tmp2 = zeros(T,14,n)
        dqdt_phi = zeros(T,n)
        dqdt_kick = zeros(T,n)
        dqdt_ij = zeros(T,14)
        jac_kepler = zeros(T,7,7)
        return new{T}(jac_phi,jac_kick,jac_copy,jac_ij,
        jac_tmp1,jac_tmp2,dqdt_phi,dqdt_kick,dqdt_ij,jac_kepler)
    end
end
