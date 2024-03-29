using LinearAlgebra
using DelimitedFiles

include("kepler_init.jl")
include("setup_hierarchy.jl")

function init_nbody(elements::Array{T,2},t0::T,IC::Array{Any,1}) where {T <: Real}

nbody = IC[1]
rkepler = zeros(typeof(t0),nbody,NDIM)
rdotkepler = zeros(typeof(t0),nbody,NDIM)

ϵ = hierarchy(IC) # Indices matrix
m = reshape(vcat(elements[:,1])[1:nbody],nbody,1) # put masses into a vector

# Create A matrix
A = zeros(typeof(t0),size(ϵ))
for i in 1:nbody, j in 1:nbody
    A[i,j] = (ϵ[i,j]*m[j])/(Σm(m,i,j,ϵ))
end

# Computes Kepler's problem
for i in 1:nbody

    ind1 = findall(isequal(-1.0),ϵ[i,:])
    ind2 = findall(isequal(1.0),ϵ[i,:])
    m1 = sum(m[ind1])
    m2 = sum(m[ind2])

    if i < nbody
        r, rdot = kepler_init(t0,m1+m2,elements[i+1,2:7])
        for j=1:NDIM
            rkepler[i,j] = r[j]
            rdotkepler[i,j] = rdot[j]
        end
    end

end

# Inverse A matrix
Ainv = inv(A)

# Cartesian coordinates
x = zeros(eltype(t0),NDIM);
x = transpose(*(Ainv,rkepler))
x = convert(Array{eltype(t0),2},x)

v = zeros(eltype(t0),NDIM)
v = transpose(*(Ainv,rdotkepler))
v = convert(Array{eltype(t0),2},v)

return x,v

end




function init_nbody(elements::Array{Float64,2}, t0::Float64, IC::Array{Any,1}, jac_init::Array{Float64,2})

nbody = IC[1]
fill!(jac_init,0.0)
jac_21 = zeros(Float64,7,7)
jac_kepler = zeros(Float64,nbody*6,nbody*7)
rkepler = zeros(Float64,nbody,NDIM)
rdotkepler = zeros(Float64,nbody,NDIM)

ϵ = hierarchy(IC) # Indices matrix
m = reshape(vcat(elements[:,1])[1:nbody],nbody,1) # put masses into a vector

# Create A matrix
A = zeros(Float64,size(ϵ))
for i in 1:nbody, j in 1:nbody
    A[i,j] = (ϵ[i,j]*m[j])/(Σm(m,i,j,ϵ))
end
# Computes the Keplerians
for i in 1:nbody

    ind1 = findall(isequal(-1.0),ϵ[i,:])
    ind2 = findall(isequal(1.0),ϵ[i,:])
    m1 = sum(m[ind1])
    m2 = sum(m[ind2])

    if i < nbody
        # Computes Kepler's problem
		r, rdot = kepler_init(t0,m1+m2,elements[i+1,2:7],jac_21)
		for j=1:3
            rkepler[i,j] = r[j]
            rdotkepler[i,j] = rdot[j]
        end
		# Save Keplerian Jacobian to a matrix. First, postitions/velocities vs. elements
		for j=1:6, k=1:6
	    	jac_kepler[(i-1)*6+j,i*7+k] = jac_21[j,k]
		end
		# Then mass derivatives
		for j=1:nbody
	    	if ϵ[i,j] != 0
	        	for k = 1:6
		    	jac_kepler[(i-1)*6+k,j*7] = jac_21[k,7]
				end
	    	end
		end
    end
end

# Create dAdm matrix
dAdm = zeros(eltype(elements),nbody,nbody,nbody)
for k in 1:nbody, i in 1:nbody, j in 1:nbody
    dAdm[i,j,k] = ((δ_(k,j)*ϵ[i,j])/Σm(m,i,j,ϵ)) -
    ((δ_(ϵ[i,j],ϵ[i,k]))*ϵ[i,j]*m[j]/(Σm(m,i,j,ϵ)^2))
end

# Create dAinvdm
Ainv = inv(A)
dAinvdm = zeros(Float64,nbody,nbody,nbody)
for k in 1:nbody
    dAinvdm[:,:,k] = -Ainv*dAdm[:,:,k]*Ainv
end

# Cartesian coordinates
x = zeros(Float64,NDIM);
x = transpose(*(Ainv,rkepler))
x = convert(Array{Float64,2},x)

v = zeros(Float64,NDIM)
v = transpose(*(Ainv,rdotkepler))
v = convert(Array{Float64,2},v)

# Cartesian derivatives
dxdm = zeros(Float64,NDIM,nbody)
dvdm = zeros(Float64,NDIM,nbody)

# Fill in jac_init
for i=1:nbody
  for k=1:nbody
    for j=1:3, l=1:7*nbody
      jac_init[(i-1)*7+j,l] += Ainv[i,k]*jac_kepler[(k-1)*6+j,l]
      jac_init[(i-1)*7+3+j,l] += Ainv[i,k]*jac_kepler[(k-1)*6+3+j,l]
    end
  end
  for k=1:nbody
    dxdm = transpose(dAinvdm[:,:,k]*rkepler)
    dvdm = transpose(dAinvdm[:,:,k]*rdotkepler)
    jac_init[(i-1)*7+1:(i-1)*7+3,k*7] += dxdm[1:3,i]
    jac_init[(i-1)*7+4:(i-1)*7+6,k*7] += dvdm[1:3,i]
  end
  # Masses are conserved:
  jac_init[i*7,i*7] = 1.0
end

return x,v
end

# Sums masses in current keplerian
function Σm(masses::Array{T,2},i::Integer,j::Integer,ϵ::Array{Float64,2}) where T<:Real
	m = zero(eltype(masses))
    for l in 1:size(masses)[1]
        m += masses[l]*δ_(ϵ[i,j],ϵ[i,l])
    end
    return m
end

# N x N Kronecker Delta
function δ_(i::T,j::T) where {T<:Real}
    if i == j
        return 1.0
    else
        return 0.0
    end
end
