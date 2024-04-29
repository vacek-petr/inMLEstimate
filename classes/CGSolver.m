classdef CGSolver < Solver
    
    properties
    stoppingCriterion % struct 
    end
    
    methods
        function cgSolver = CGSolver(stoppingCriterion)
            cgSolver@Solver()
            cgSolver.stoppingCriterion = stoppingCriterion;
        end
        
        function  [approx,info] = apply(cgSolver,A,b,approx)
        [approx,info] = cg(A,b,approx,cgSolver.stoppingCriterion);
        end
    end
end

