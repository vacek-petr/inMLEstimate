classdef MatrixHierarchy< handle
    % MATRIXHIERARCHY class for storing the hierarchy of stiffness and
    % prolongation matrices, their properties and the right-hand sides

    properties
        name    % string 
        A   % cell array - stiffness matrices, A{1} - coarsest level,  A{end} - finest level
        P   % cell array - prolongation matrices, P{1} empty, P{2} prolongation matrix from level 1 to level 2
        F   % cell array - right-hand sides
        numberOfLevels  % number
        ASmallestEigenvalues   % cell array 
        solution  % cell array - approx of A{j}^-1F{j} computedy using Matlab backslash 
    end

    methods
        function obj = MatrixHierarchy()
        % MATRIXHIERARCHY Construct an instance of class MatrixHierarchy
        end
        function selectLevels(obj,options)
            arguments
                obj
                options.numberOfLevels = obj.numberOfLevels;
                options.from = 1;
            end
            obj.numberOfLevels = options.numberOfLevels;
            from  = options.from;
            to = from + obj.numberOfLevels - 1;
            obj.A = obj.A(from:to);
            obj.P = obj.P(from:to);
            obj.F = obj.F(from:to);
            obj.solution = obj.solution(from:to);
            obj.ASmallestEigenvalues = obj.ASmallestEigenvalues(from:to);
        end
    end
end