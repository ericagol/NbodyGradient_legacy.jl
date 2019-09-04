using NbodyGradient
using Test

maxabs(x) = maximum(abs.(x))

const YEAR  = 365.242
const GNEWT = 39.4845/YEAR^2
const NDIM  = 3
const third = 1.0/3.0
const alpha0 = 0.0

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
