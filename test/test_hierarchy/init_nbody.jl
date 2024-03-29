include("kepler_init.jl")
include("setup_hierarchy.jl")

function init_nbody(elements::Array{T,2}, t0::T, IC::Array{Any,1}) where {T <: Real}

    n_body = IC[1]

    h = hierarchy(IC)

    mass = vcat(elements[:,1])

    rkepler = zeros(eltype(elements),n_body,NDIM)
    rdotkepler = zeros(eltype(elements),n_body,NDIM)

    for i in 1:n_body-1
        index1 = findall(isequal(-1.0),h[i,:])
        iter1 = length(index1)
        if iter1 == 1
            h[i,index1] .= -1.0
            mc1 = mass[index1[1]]
        elseif iter1 > 1
            mc1 = summass(mass,index1,iter1)
            h[i,index1] = -mass[index1]/mc1
        end
        index2 = findall(isequal(1.0),h[i,:])
        iter2 = length(index2)
        if iter2 == 1
            h[i,index2] .= 1.0
            mc2 = mass[index2[1]]
        elseif iter2 > 1
            mc2 = summass(mass,index2,iter2)
            h[i,index2] = mass[index2]/mc2
        end
        r,rdot = kepler_init(t0,mc1+mc2,elements[i+1,2:7])
        for j=1:3
            rkepler[i,j] = r[j]
            rdotkepler[i,j] = rdot[j]
        end
    end
    
    index = findall(isequal(-1.0),h[end,:])
    iter = length(index)
    mtot = summass(mass,index,iter)
    h[end,:] .= mass[index]/mtot
	

    ainv = inv(amat)

    x = zeros(eltype(elements),NDIM,n_body)
    v = zeros(eltype(elements),NDIM,n_body)

    x = transpose(*(ainv,rkepler)) ; x = convert(Array{Float64,2},x)
    v = transpose(*(ainv,rdotkepler)); v = convert(Array{Float64,2},x)

    return x,v

end

function init_nbody(elements::Array{Float64,2},t0::Float64,IC::Array{Any,1},jac_init::Array{Float64,2})
    n_body = IC[1]
    
    fill!(jac_init,0.0)

    h = hierarchy(IC)
    dhdm = zeros(eltype(elements),n_body,n_body,n_body)

    mass = vcat(elements[:,1])
    rkepler = zeros(eltype(elements),n_body,NDIM)
    rdotkelper = zeros(eltype(elements),n_body,NDIM)

    jac_21 = zeros(Float64,7,7)
    jac_kepler = zeros(Float64,n_body*6,n_body*7)

    
# FILL A MATRIX #
    for i in 1:n_body-1
        index1 = findall(isequal(-1.0),h[i,:])
        iter1 = length(index1)
        if iter1 == 1
            h[i,index1] .= -1.0
            mc1 = mass[index1[1]]
            dhdm[i,index1,:] .= -1.0
        elseif iter1 > 1
            mc1 = summass(mass,index1,iter1)
            h[i,index1] = -mass[index1]/mc1
	    dhdm[i,]  
        end
        index2 = findall(isequal(1.0),h[i,:])
        iter2 = length(index2)
        if iter2 == 1
            h[i,index2] .= 1.0
            mc2 = mass[index2[1]]
        elseif iter2 > 1
            mc2 = summass(mass,index2,iter2)
            h[i,index2] = mass[index2]/mc2
        end
        r,rdot = kepler_init(t0,mc1+mc2,elements[i+1,2:7])
        for j=1:NDIM
            rkepler[i,j] = r[j]
            rdotkepler[i,j] = rdot[j]
        end

# SAVE KEPLERIAN JACOBIAN TO MATRIX
        for j=1:6, k=1:6
	    jac_kepler[(i-1)*6+j,i*7+k] = jac_21[j,k]
	end
        

    end
    index = findall(isequal(-1.0),h[end,:])
    iter = length(index)
    mtot = summass(mass,index,iter)
    h[end,:] .= mass[index]/mtot

    	

end

################################################################################
# Sums the masses in a range of indices. Called in amat functions.
################################################################################
function summass(mass::Array{T,1},index::Array{Int64,1},iter::Int64) where {T<:Real}
    tmass = 0.0
    for i in 1:iter
        tmass += mass[index[i]]
    end
    return tmass
end

