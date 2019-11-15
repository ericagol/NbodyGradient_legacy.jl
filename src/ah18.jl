function ah18!(sim::State{T},h::T,h2::T,der::Bool) where T <: Real
    drift!(sim,h2)
    kickfast!(sim,h2)
    phisalpha!(sim.x,sim.v,h,sim.m,convert(T,2),sim.n,
               sim.d.jac_phi,sim.d.dqdt_phi,sim.pair)
    @inbounds for i=sim.n-1:-1:1, j=sim.n:-1:i+1
        if ~sim.pair[i,j]
            kepler_driftij!(sim.m,sim.x,sim.v,i,j,h2,true)
        end
    end
    kickfast!(sim,h2)
    drift!(sim,h2)
    return
end

function ah18!(sim::State,h::T,h2::T) where T <: Real
    drift!(sim,h2)
    kickfast!(sim,h2,sim.d)
    jac_kick_step!(sim.jac_step,sim.d.jac_kick,sim.d.jac_copy)
    indi = 0; indj = 0
    @inbounds for i=1:sim.n-1
        indi = (i-1)*7
        for j=i+1:sim.n
            indj = (j-1)*7
            if ~sim.pair[i,j]
                kepler_driftij!(sim,sim.d,i,j,h2,true)
                jac_tmp!(sim,sim.d,indi,indj)
            end
        end
    end
    # Missing phic here [ ]
    phisalpha!(sim.x,sim.v,h,sim.m,convert(T,2.0),sim.n,
               sim.d.jac_phi,sim.d.dqdt_phi,sim.pair)
    jac_kick_step!(sim.jac_step,sim.d.jac_phi,sim.d.jac_copy)
    indi=0;indj=0
    @inbounds for i=sim.n-1:-1:1
        indi = (i-1)*7
        for j=sim.n:-1:i+1
            indj = (j-1)*7
            if ~sim.pair[i,j]
                kepler_driftij!(sim,sim.d,i,j,h2,false)
                jac_tmp!(sim,sim.d,indi,indj)
            end
        end
    end
    kickfast!(sim,h2,sim.d)
    jac_kick_step!(sim.jac_step,sim.d.jac_kick,sim.d.jac_copy)
    drift!(sim,h2)
    return
end

function jac_kick_step!(jac1::Array{T,2},jac2::Array{T,2},
                        jac_copy::Array{T,2}) where T<:Float64
    @inbounds for i in eachindex(jac1)
        jac_copy[i] = jac1[i]
    end
    BLAS.gemm!('N','N',one(T),jac2,jac_copy,zero(T),jac1)
end

function jac_kick_step!(jac1::Array{BigFloat,2},jac2::Array{BigFloat,2},
                        jac_copy::Array{BigFloat,2})
    @inbounds for i in eachindex(jac1)
        jac_copy[i] = jac1[i]
    end
    jac1 = *(jac2,jac_copy)
    return
end

function jac_tmp!(sim::State{T},d::Derivatives{T},indi::Int64,indj::Int64) where T <: Real
    @inbounds for k2=1:(7*sim.n), k1=1:7
        d.jac_tmp1[k1,k2] = sim.jac_step[indi+k1,k2]
    end
    @inbounds for k2=1:(7*sim.n), k1=1:7
        d.jac_tmp1[7+k1,k2] = sim.jac_step[indj+k1,k2]
    end
    # Carry out multiplication on the i/j components of matrix
    if T == BigFloat
        d.jac_tmp2 = *(d.jac_ij,d.jac_tmp1)
    else
        BLAS.gemm!('N','N',one(T),d.jac_ij,d.jac_tmp1,zero(T),d.jac_tmp2)
    end
    # Copy back to the Jacobian
    @inbounds for k2=1:(7*sim.n), k1=1:7
        sim.jac_step[indi+k1,k2] = d.jac_tmp2[k1,k2]
    end
    @inbounds for k2=1:(7*sim.n), k1=1:7
        sim.jac_step[indj+k1,k2] = d.jac_tmp2[7+k1,k2]
    end
    return
end

function drift!(sim::State{T},h2::T,der::Bool) where T<: Real
    @inbounds for i=1:sim.n, j=1:sim.NDIM
        sim.x[j,i] += h2*sim.v[j,i]
     end
     return
end

function drift!(sim::State{T},h2::T) where T <: Real
    indi = 0
    @inbounds for i=1:sim.n 
        for j=1:sim.NDIM
            sim.x[j,i] += h2*sim.v[j,i]
        end 
        indi = (i-1)*7
        for k=1:7*sim.n, j=1:sim.NDIM
            sim.jac_step[indi+j,k] += h2*sim.jac_step[indi+3+j,k]
        end
    end
    return
end

#=
function kickfast!(init::IC{T},sim::Sim{T},h2::T) where T <: Real
    rij = zeros(T,3)
    @inbounds for i=1:sim.n-1, j=i+1:sim.n
        if sim.pair[i,j]
            r2 = zero(T)
            for k=1:3
                rij[k] = x[k,i] - x[k,j]
                r2 += rij[k]^2
            end
            r3_inv = one(T)/(r2*sqrt(r2))
            for k=1:3
                fac = h2*sim.GNEWT*rij[k]*r3_inv
                sim.v[k,i] -= init.m[j]*fac
                sim.v[k,j] += init.m[i]*fac
            end
        end
    end
    return
end
=#
function kickfast!(sim::State{T},h2::T,d::Derivatives{T}) where T <: Real
    rij = zeros(T,3)
    d.jac_kick .= Matrix{T}(I,7*sim.n,7*sim.n)

    @inbounds for i=1:sim.n-1
        indi = (i-1)*7
        @inbounds for j=i+1:sim.n
            indj = (j-1)*7
            if sim.pair[i,j]
                @inbounds for k=1:3
                    rij[k] = sim.x[k,i] - sim.x[k,j]
                end
                r2inv = 1.0/(rij[1]*rij[1]+rij[2]*rij[2]+rij[3]*rij[3])
                r3inv = r2inv*sqrt(r2inv)
                @inbounds for k=1:3
                    fac = h2*sim.GNEWT*rij[k]*r3inv
                    sim.v[k,i] -= sim.m[j]*fac
                    sim.v[k,j] += sim.m[i]*fac
                    d.dqdt_kick[indi+3+k] -= sim.m[j]*fac/h2
                    d.dqdt_kick[indj+3+k] += sim.m[i]*fac/h2
                    d.jac_kick[indi+3+k,indj+7] -= fac
                    d.jac_kick[indj+3+k,indi+7] += fac
                    fac *= 3.0*r2inv
                    @inbounds for p=1:3
                        d.jac_kick[indi+3+k,indi+p] += fac*sim.m[j]*rij[p]
                        d.jac_kick[indi+3+k,indj+p] -= fac*sim.m[j]*rij[p]
                        d.jac_kick[indj+3+k,indj+p] += fac*sim.m[i]*rij[p]
                        d.jac_kick[indj+3+k,indi+p] -= fac*sim.m[i]*rij[p]
                    end
                    fac = h2*sim.GNEWT*r3inv
                    d.jac_kick[indi+3+k,indi+k] -= fac*sim.m[j]
                    d.jac_kick[indi+3+k,indj+k] += fac*sim.m[j]
                    d.jac_kick[indj+3+k,indj+k] -= fac*sim.m[i]
                    d.jac_kick[indj+3+k,indi+k] += fac*sim.m[i]
                end
            end
        end
    end
    return
end

function kepler_driftij!(sim::State{T},d::Derivatives{T},i::Int64,j::Int64,h::T,drift_first::Bool) where {T <: Real}

zero_out!(sim.DRIFT.state0)
zero_out!(sim.DRIFT.state)

for k=1:sim.NDIM
    sim.DRIFT.state0[1+k] = sim.x[k,i] - sim.x[k,j]
    sim.DRIFT.state0[4+k] = sim.v[k,i] - sim.v[k,j]
end

gm = sim.GNEWT*(sim.m[i]+sim.m[j])

for i in 1:size(d.jac_ij)[1], j in 1:size(d.jac_ij)[2]
    d.jac_ij[i,j] = 0.0
    if i == j
        d.jac_ij[i,j] = 1.0
    end
end

if gm != 0
    zero_out!(d.jac_kepler)

    kepler_drift_step!(gm, h, sim.DRIFT.state0, sim.DRIFT.state, d.jac_kepler,drift_first)
    mijinv =1.0/(sim.m[i] + sim.m[j])
    mi = sim.m[i]*mijinv # Normalize the masses
    mj = sim.m[j]*mijinv
    for k=1:3
        sim.x[k,i] += mj*sim.DRIFT.state[1+k]
        sim.v[k,i] += mj*sim.DRIFT.state[4+k]
        sim.x[k,j] -= mi*sim.DRIFT.state[1+k]
        sim.v[k,j] -= mi*sim.DRIFT.state[4+k]
    end
    for l=1:6, k=1:6
        d.jac_ij[  k,  l] += mj*d.jac_kepler[k,l]
        d.jac_ij[  k,7+l] -= mj*d.jac_kepler[k,l]
        d.jac_ij[7+k,  l] -= mi*d.jac_kepler[k,l]
        d.jac_ij[7+k,7+l] += mi*d.jac_kepler[k,l]
    end
    for k=1:6
        d.jac_ij[   k, 7] = -mj*sim.DRIFT.state[1+k]*mijinv + sim.GNEWT*mj*d.jac_kepler[  k,7]
        d.jac_ij[   k,14] =  mi*sim.DRIFT.state[1+k]*mijinv + sim.GNEWT*mj*d.jac_kepler[  k,7]
        d.jac_ij[ 7+k, 7] = -mj*sim.DRIFT.state[1+k]*mijinv - sim.GNEWT*mi*d.jac_kepler[  k,7]
        d.jac_ij[ 7+k,14] =  mi*sim.DRIFT.state[1+k]*mijinv - sim.GNEWT*mi*d.jac_kepler[  k,7]
    end
end
return
end
