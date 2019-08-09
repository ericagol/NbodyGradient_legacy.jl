# Initial Conditions generation
using DelimitedFiles
include("/Users/langfzac/Research/NbodyGradient.jl/src/kepler_init.jl")
include("/Users/langfzac/Research/NbodyGradient.jl/src/setup_hierarchy.jl")

abstract type InitialConditions
end

abstract type SystemState
end

struct IC <: InitialConditions
    elements::Array{<:Real,2}
    ϵ::Array{<:Real,2}
    amat::Array{<:Real,2}
    derivatives::Bool
    NDIM::Integer
    nbody::Integer
    m::Array{<:Real,2}

    ## Constructor for IC composite type
    function IC(filename::String,system::Array{Any,1};
                der::Bool = false, NDIM::Integer = 3)
        ϵ = hierarchy(system)
        elements = readdlm(filename,',',comments=true)
        nbody = length(ϵ[:,1])
        m = reshape(vcat(elements[:,1])[1:nbody],nbody,1)
        amat = amatrix(elements,ϵ,m)
        return new(elements,ϵ,amat,der,NDIM,nbody,m)
    end
end

function init_nbody(init::IC,t0::T) where T <: Real

    if !isa(t0,eltype(init.elements))
         error("t0 needs to be same type as values in elements array")
    end

    r, rdot = kepcalc(t0,init)
    Ainv = inv(init.amat)

    # Cartesian coordinates
    x = zeros(eltype(t0),init.NDIM)
    x = transpose(*(Ainv,r))
    x = convert(Array{eltype(t0),2},x)

    v = zeros(eltype(t0),init.NDIM)
    v = transpose(*(Ainv,rdot))
    v = convert(Array{eltype(t0),2},v)

    # Calculate and return mass derivatives if true
    if init.derivatives
        jac_init = d_dm(init)
        return x,v,jac_init
    else
        return x,v
    end
end

# Compute derivatives of
function d_dm(init::IC)

    N = init.nbody
    m = init.m
    jac_init = zeros(eltype(m),7*N,7*N)
    jac_kepler = zeros(eltype(m),6*N,7*N)
    dAdm = zeros(eltype(m),N,N,N)
    dxdm = zeros(eltype(m),NDIM,N)
    dvdm = zeros(eltype(m),NDIM,N)

    # Differentiate A matrix wrt the mass of each body
    for i in 1:N, i in 1:N, j in 1:N
        dAdm[i,j,k] = ((δ_(k,j)*ϵ[i,j])/Σm(m,i,j,ϵ)) -
        ((δ_(ϵ[i,j],ϵ[i,k]))*ϵ[i,j]*m[j]/(Σm(m,i,j,ϵ)^2))
    end

    # Calculate inverse of dAdm
    Ainv = inv(init.amat)
    dAinvdm = zeros(eltype(m),N,N,N)
    for k in 1:N
        dAinvdm[:,:,k] = -Ainv*dAdm[:,:,k]*Ainv
    end

    # Fill in jac_init array
    for i in 1:N
        for k in 1:N
            for j in 1:3, l in 1:7*N
                jac_init[(i-1)*7+j,l] += Ainv[i,k]*jac_kepler[(k-1)*6+j,l]
                jac_init[(i-1)*7+3+j,l] += Ainv[i,k]*jac_kepler[(k-1)*6+3+j,l]
            end
        end

        # Derivatives of cartesian coordinates wrt masses
        for k in 1:N
            dxdm = transpose(dAinvdm[:,:,k]*rkepler)
            dvdm = transpose(dAinvdm[:,:,k]*rdotkepler)
            jac_init[(i-1)*7+1:(i-1)*7+3,k*7] += dxdm[1:3,i]
            jac_init[(i-1)*7+4:(i-1)*7+6,k*7] += dvdm[1:3,i]
        end
        jac_init[i*7,i*7] = 1.0
    end
    return jac_init
end

# Computes Kepler's problem
function kepcalc(t0::T,init::IC) where T <: Real

    # Kepler position/velocity arrays (body,pos/vel)
    rkepler = zeros(eltype(t0),init.nbody,init.NDIM)
    rdotkepler = zeros(eltype(t0),init.nbody,init.NDIM)

    # Compute Kepler's problem
    for i in 1:init.nbody

        ind1 = findall(isequal(-1.0),init.ϵ[i,:])
        ind2 = findall(isequal(1.0),init.ϵ[i,:])
        m1 = sum(init.m[ind1])
        m2 = sum(init.m[ind2])

        if i < init.nbody
            r, rdot = kepler_init(t0,m1+m2,init.elements[i+1,2:7])
            for j = 1:init.NDIM
                rkepler[i,j] = r[j]
                rdotkepler[i,j] = rdot[j]
            end
        end
    end
    return rkepler, rdotkepler
end

## Creates A matrix from elements and ϵ matrix
function amatrix(elements::Array{T,2},ϵ::Array{T,2},m::Array{T,2}) where T<:Real
    A = zeros(eltype(ϵ),size(ϵ)) # Empty A matrix
    N = length(ϵ[:,1]) # Number of bodies in system

    for i in 1:N, j in 1:N
        A[i,j] = (ϵ[i,j]*m[j])/(Σm(m,i,j,ϵ))
    end
    return A
end

# Sums masses in current keplerian
function Σm(masses::Array{T,2},i::Integer,j::Integer,
            ϵ::Array{T,2}) where T <: Real
    m = 0.0
    for l in 1:size(masses)[1]
        m += masses[l]*δ_(ϵ[i,j],ϵ[i,l])
    end
    return m
end

# N x N Kronecker Delta
function δ_(i::T,j::T) where T <: Real
    if i == j
        return 1.0
    else
        return 0.0
    end
end

# Simply system test. Need to add other params (GWENT, etc.)
function test_type()
    global GNEWT = 39.4845/365.242^2
    global third = 1.0/3.0
    elements = "test/elements.txt"
    sys = [4,"1,1,1"]
    t0 = 7257.93115525
    init = IC(elements,sys)
    return init_nbody(init,t0)
end
