
#include("../src/ttv.jl")
#include("/Users/ericagol/Computer/Julia/regress.jl")

#maxabs(x) = maximum(abs.(x))

using DelimitedFiles
using LinearAlgebra
#KEPLER_TOL  = 1e-15
#TRANSIT_TOL = 1e-15

@testset "findtransit3" begin

#n = 8
n = 3
t0 = 7257.93115525
#h  = 0.12
h = 0.05
#tmax = 600.0
tmax = 100.0
#tmax = 10.0
IC = [3,"1,1"]

# Read in initial conditions:
elements = readdlm("elements.txt",',',comments=true)
# Increase masses for debugging:
elements[2,1] *= 10.0
elements[3,1] *= 10.0

# Make an array, tt,  to hold transit times:
# First, though, make sure it is large enough:
ntt = zeros(Int64,n)
for i=2:n
  global ntt[i] = ceil(Int64,tmax/elements[i,2])+3
end
tt  = zeros(n,maximum(ntt))
# Save a counter for the actual number of transit times of each planet:
count = zeros(Int64,n)
# Now, compute derivatives (with respect to initial cartesian positions/masses at
# beginning of each step):
dtdq0 = zeros(n,maximum(ntt),7,n)
dtdq0_num = zeros(BigFloat,n,maximum(ntt),7,n)
dlnq = big(1e-10)
# Make radius of star large:
rstar = 1e12
dtdelements_num = ttv_elements!(n,t0,h,tmax,elements,IC,tt,count,dtdq0,dtdq0_num,dlnq,rstar)

dtdq0_num = convert(Array{Float64,4},dtdq0_num)

mask = zeros(Bool, size(dtdq0))
for i=2:n, j=1:count[i], k=1:5, l=1:n
  global mask[i,j,k,l] = true
end
#println("Max diff log(dtdq0): ",maximum(abs.(dtdq0_num[mask]./dtdq0[mask]-1.0)))
#println("dtdq0: ", dtdq0_num[mask])
#println("Max diff asinh(dtdq0): ",maximum(abs.(asinh.(dtdq0_num[mask])-asinh.(dtdq0[mask]))))
#println("Max diff     dtdq0 : ",maximum((dtdq0_num[mask]-dtdq0[mask])))
println("Max diff asinh(dtdq0): ",maximum(abs.(asinh.(dtdq0_num[mask])-asinh.(dtdq0[mask]))))
println("Max diff     dtdq0 : ",maximum((dtdq0_num[mask]-dtdq0[mask])))
@test isapprox(asinh.(dtdq0[mask]),asinh.(dtdq0_num[mask]);norm=maxabs)
#unit = ones(dtdq0[mask])
#@test isapprox(dtdq0[mask]./convert(Array{Float64,4},dtdq0_num)[mask],unit;norm=maxabs)
end
