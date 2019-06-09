using DelimitedFiles
using LinearAlgebra

include("/Users/zllangford/Desktop/Research/NbodyGradient/test/test_hierarchy/setup_hierarchy.jl")
include("/Users/zllangford/Desktop/Research/NbodyGradient/test/test_hierarchy/amat.jl")
include("/Users/zllangford/Desktop/Research/NbodyGradient/src/kepler_init.jl")

global t0 = 7257.93115525
global NDIM = 3
const YEAR = 365.242
const GNEWT = 39.4845/YEAR^2
const TRANSIT_TOL = 10.0*sqrt(eps(1.0))
const third = 1.0/3.0
const alpha0 = 0.0

elements = readdlm("/Users/zllangford/Desktop/Research/NbodyGradient/test/elements.txt",',',skipstart=1)
IC = [8,"1,1,1,1,1,1,1"]
h = hierarchy(IC)
n_body = IC[1]
jac_init = zeros(Float64,7*n_body,7*n_body)
jac_21 = zeros(Float64,7,7)
jac_kepler = zeros(Float64,n_body*6,n_body*7)
rkepler = zeros(eltype(elements),n_body,NDIM)
rdotkepler = zeros(eltype(elements),n_body,NDIM)

amat_kep_jac(h,elements,IC,rkepler,rdotkepler,jac_init)
show(stdout, "text/plain", h)
println(" ")
