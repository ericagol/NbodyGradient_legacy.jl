#Testing Dependencies
#include("setup_hierarchy.jl")
#elements = readdlm("../elements.txt",',',skipstart=1)
#IC = [8,"1,1,1,1,1,1,1"]
#h = hierarchy(IC)
#n_body = IC[1]
#jac_init = zeros(Float64,7*n_body,7*n_body)

################################################################################
# Generates the A matrix from Hamers & Portegies Zwart 2016.
#
# Generates ONLY the A matrix.
################################################################################
function amat1(h::Array{},elements::Array{},IC::Array{Any,1}) where {T <: Real}

    mass = vcat(elements[:,1])

    for i in 1:length(h[:,1])
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
    end
end

################################################################################
# Generates A matrix from Hamers & Portegies Zwart (2016). Calls kepler_init()
# to calcuate initial Cartesian positions and velocities.
################################################################################
function amat_kep(h::Array{T},elements::Array{T,2},IC::Array{Any,1},rkepler::Array{T},rdotkepler::Array{T}) where {T <: Real}

    mass = vcat(elements[:,1])

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
end

################################################################################
# INCOMPLETE
# Generates A matrix from Hamer & Portegies Zwart (2016). Calls kepler_init() to
# calculate initial Cartesian position and velocities.
#
# Need to finish adding in mass derivatives.
################################################################################
"""
function amat_kep_jac(h::Array{T},rkepler::Array{T,2},rdotkepler::Array{T,2},elements::Array{T,2},IC::Array{Any,1},jac_init::Array{Float64,2}) where {T <: Real}

    mass = vcat(elements[:,1])
    rkepler = zeros(eltype(elements),n_body,3)
    rdotkepler = zeros(eltype(elements),n_body,3)
    h = convert(eltype(elements),h)
    dh = zeros(eltype(elements),n_body,n_body,n_body)

    for i in 1:length(h[:,1])
        index1 = findall(isequal(-1.0),h[i,:])
        iter1 = length(index1)
        if iter1 == 1
            mc1 = mass[index1[1]]
            h[i,index1] .= -mass[index1]/mc1
            dh[i,index1,1:n_body] .= 0
        elseif iter1 > 1
            mc1 = summass(mass,index1,iter1)
            h[i,index1] .= -mass[index1]/mc1
            dh[i,index1,1:n_body] .= mass[index1]/mc1^2
        end
        index2 = findall(isequal(1.0),h[i,:])
        iter2 = length(index2)
        if iter2 == 1
            mc2 = mass[index2[1]]
            h[i,index2] .= -1.0
            dh[i,index1,1:n_body] .= 0
        elseif iter2 > 1
            mc2 = summass(mass,index2,iter2)
            h[i,index2] = mass[index2]/mc2
            dh[i,index2,1:n_body] .= -mass[index2]/mc2^2
        end
        r,rdot = kepler_init(t0,mc1+mc2,elements[i+1,2:7],jac_init)
        for j=1:3
          rkepler[i,j] = r[j]
          rdotkepler[i,j] = rdot[j]
        end
        for j=1:6, k=1:6
          jac_kepler[(i-1)*6+j,i*7+k] = jac_21[j,k]
        end
        for j=1:n_body
          # Check which bodies participate in current Keplerian:
          if indices[i,j] != 0
            for k=1:6
              jac_kepler[(i-1)*6+k,j*7] = jac_21[k,7]
            end
          end
        end
    end
end
"""

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
