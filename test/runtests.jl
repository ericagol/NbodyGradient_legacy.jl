include("../src/ttv.jl")
maxabs(x) = maximum(abs.(x))
using Test

include("test_kepler_init.jl")

# init_nbody.jl working, test is not.
#include("test_init_nbody.jl")
include("test_elliptic_derivative.jl")
include("test_keplerij.jl")
include("test_kickfast.jl")
include("test_phisalpha.jl")
include("test_dh17.jl")

# Not working
#include("test_findtransit2.jl")
include("test_ttv_cartesian.jl")
include("test_ttv_elements.jl")
