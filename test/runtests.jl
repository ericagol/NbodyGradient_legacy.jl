include("../src/ttv.jl")
maxabs(x) = maximum(abs.(x))
using Test

include("test_kepler_init.jl")
include("test_init_nbody.jl")
include("test_elliptic_derivative.jl")
include("test_keplerij.jl")
include("test_kickfast.jl")
include("test_phisalpha.jl")
include("test_dh17.jl")
include("test_findtransit3.jl")
include("test_ttv_cartesian.jl")
include("test_ttv_elements.jl")
