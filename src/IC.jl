"""
    InitialConditions

Abstract type for [`IC`](@ref).
"""
abstract type InitialConditions end

"""
    IC{T<:Real} <: InitialConditions

Structure for holding the initial conditions of a nested hierarchical system.

# Fields
- `elements::Array{<:Real,2}`: Array of masses and initial orbital elements of the systems.
- `ϵ::Array{<:Real,2}`: Matrix of relative position vectors for each level of the hierarchical system. (Hamers & Portegies Zwart 2016).
- `amatrix::Array{<:Real,2}`: A-matrix. Adds center of mass elements to ϵ.
- `derivatives::Bool`: Determines if derivatives are to be calculated.
- `NDIM::Integer`: Number of spatial dimensions.
- `nbody::Integer`: Number of bodies in the system
- `m::Array{<:Real,2}`: Masses of each body.
# Inner Constructor Arguements
- `filename::String`: Name of file containing initial orbital elements and masses.
- `system::Array{Int64,1}`: Array with the first index being the number of bodies in the system and the subsequent indices the number of binarys on that level of the hierarchy.
# Options
- `der::Bool`: Determines whether derivatives are to be computed. Default = `true`.
- `NDIM::Int64`: Number of spatial dimensions. Default = `3`
- `prec`: Precision for all initial conditions. Needs to be of a type <: Real. Default = `Float64`
# Example
The following uses the inner constructor to create an IC type for a simple system of 4 bodies in concentric orbits.
```julia
filename = "/my/path/to/elements.txt" # File containing orbital elements and masses
system = [4,1,1,1] # Example system array syntax. 4 bodies in concentric orbits.
init = IC(filename,system)
```
"""
mutable struct IC{T<:AbstractFloat} <: InitialConditions
    elements::Array{T,2}
    ϵ::Array{T,2}
    amat::Array{T,2}
    derivatives::Bool
    NDIM::Int64
    nbody::Int64
    m::Array{T,2}
    t0::T
    # Inner Constructor
    function IC(filename::String,system::Array{Int64,1},t0::T;
                der::Bool = true, NDIM::Int64 = 3,) where T <: AbstractFloat
        ϵ = convert(Array{T},hierarchy(system))
        elements = convert(Array{T},readdlm(filename,',',comments=true))
        nbody = system[1]
        m = reshape(vcat(elements[:,1])[1:nbody],nbody,1)
        amat = amatrix(ϵ,m)
        return new{T}(elements,ϵ,amat,der,NDIM,nbody,m,t0);
    end
end
