classdef Solver
% SOLVER class for calling solver for Ax = b

methods
    function solver = Solver()
    end
end

methods (Abstract)
[approx,info] = apply(solver)
end

end