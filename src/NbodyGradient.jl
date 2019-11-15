module NbodyGradient
    ################ Abstract Types ################
    
        abstract type InitialConditions end
        abstract type Simulation end
        abstract type PreAllocArrays end
   

    #################### Usings ####################
        using Reexport
        @reexport using DelimitedFiles
        using Distributed
        using ForwardDiff
        @reexport using LinearAlgebra
        using Test
        using Statistics
    
    ################### Includes ###################

        include("InitialConditions.jl")
        include("init_nbody.jl")
        include("PreAllocArrays.jl")
        include("State.jl")
        include("Integrators.jl")
        include("ah18.jl")
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
		InitialConditions, Simulation, Integrator, 
        State, ElementsIC, CartesianIC, 
        amatrix, init_nbody, kepler_init,  
        ttv_elements!, ttv!, drift!, kepler_drift_step!,
		kepler_driftij!,kep_ell_hyp!, keplerij!, 
        dh17!, kickfast!, phisalpha!, ah18!

end
