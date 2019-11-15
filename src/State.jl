mutable struct State{T<:Real} <: Simulation
    IC::InitialConditions
    d::Derivatives{T}
    m::Array{T,2}
    x::Array{T,2}
    v::Array{T,2}
    jac_init::Array{T,2}
    jac_step::Array{T,2}
    t::T
    pair::Array{Bool,2}
    n::Int64
    NDIM::Int64
    GNEWT::T
    DRIFT::DriftArrays{T}

    # Oribital elements constructor
    function State(init::ElementsIC{T}) where T <: Real
        pair=zeros(Bool,init.nbody,init.nbody)
        x,v,jac_init = init_nbody(init)
        m = init.m
        n = init.nbody
        t = init.t0
        NDIM = init.NDIM
        GNEWT = convert(T,(39.4845/(365.242)^2))
        jac_step = Matrix{T}(I,7*n,7*n)
        d = Derivatives(init)
        DRIFT = DriftArrays(t)
        return new{T}(init,d,m,x,v,jac_init,jac_step,t,
                      pair,n,NDIM,GNEWT,DRIFT)
    end

    # Cartesian coordinate constructor
    function State(init::CartesianIC{T}) where T <: Real
        pair=zeros(Bool,init.nbody,init.nbody)
        x = init.x
        v = init.v
        jac_init = init.jac_init
        m = init.m
        n = init.nbody
        t = init.t0
        NDIM = init.NDIM
        GNEWT = convert(T,(39.4845/(365.242)^2))
        jac_step = Matrix{T}(I,7*n,7*n)
        d = Derivatives(init)
        DRIFT = DriftArrays(init)
        return new{T}(init,d,m,x,v,jac_init,jac_step,t,
                      pair,n,NDIM,GNEWT)
    end
end

function zero_out!(a::Array{T}) where T <: Real
    @inbounds for i in 1:length(a)
        a[i] = 0.0
    end
end

function zero_out!(d::Derivatives{T}) where T <: Real
    @inbounds for i in 1:size(d.jac_phi)[1], j in 1:size(d.jac_phi)[2]
        d.jac_phi[i,j] = 0.0
        if i == j
            d.jac_phi[i,j] = 1.0
        end
    end
    @inbounds for i in 1:length(d.jac_kick); d.jac_kick[i] = 0.0; end;
    @inbounds for i in 1:length(d.jac_copy); d.jac_copy[i] = 0.0; end;
    @inbounds for i in 1:length(d.jac_ij);   d.jac_ij[i] = 0.0;   end;
    @inbounds for i in 1:length(d.jac_tmp1); d.jac_tmp1[i] = 0.0; end;
    @inbounds for i in 1:length(d.jac_tmp2); d.jac_tmp2[i] = 0.0; end;
    @inbounds for i in 1:length(d.dqdt_phi); d.dqdt_phi[i] = 0.0; end;
    @inbounds for i in 1:length(d.dqdt_kick) d.dqdt_kick[i] = 0.0;end;
    @inbounds for i in 1:length(d.dqdt_ij) d.dqdt_ij[i] = 0.0;    end;
    return
end

function comp_sum(sim::State{T},sum_error::T) where T <: Real
    sum_error += sim.h
    tmp = sim.t + sum_error
    sum_error = (sim.t - tmp) + sum_error
    sim.t = tmp
    return sum_error::T
end
