module NbodyGradient

    #################### Usings ####################

        using DelimitedFiles
        using Distributed
        using ForwardDiff
        using LinearAlgebra
        using Test
        using Statistics

    ################### Includes ###################

        include("types.jl")
        include("init_nbody.jl")
        include("g3.jl")
        include("kepler.jl")
        include("kepler_drift_solver.jl")
        include("kepler_drift_step.jl")
        include("kepler_init.jl")
        include("kepler_solver_derivative.jl")
        include("kepler_step.jl")
        include("setup_hierarchy.jl")
        include("ttv.jl")

	################### Exports ###################

		export
		InitialConditions, IC, amatrix,
        init_nbody, kepler_init, 
        ttv_elements!, ttv!,
		kep_ell_hyp!, keplerij!, dh17!,
		kickfast!, phisalpha!, ah18!

end
