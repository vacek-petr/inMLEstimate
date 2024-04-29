classdef MultigridHierarchy< handle
    % MULTIGRIDHIERARCHY class for storing the hierarchy of stiffness and
    % prolongation matrices, their properties and the right-hand sides

    properties
        name    % string 
        A   % cell array - stiffness matrices, A{1} - coarsest level,  A{end} - finest level
        P   % cell array - prolongation matrices, P{1} empty, P{2} prolongation matrix from level 1 to level 2
        F   % cell array - right-hand sides
        numberOfLevels  % number
        ASmallestEigenvalues   % cell array 
        solution  % cell array - approx of A{j}^-1F{j}
    end

    methods
        function mh = MultigridHierarchy(name)
            % MULTIGRIDHIERARCHY Construct an instance of class MatrixHierarchy
            mh.name = name;
        end
    end
end