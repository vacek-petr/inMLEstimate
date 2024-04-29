classdef BackSlashSolver < Solver

    methods
        function backSlashSolver = BackSlashSolver()
            backSlashSolver@Solver();
        end
        
        function [approx,info] = apply(backSlashSolver,A,b,~)
        approx = A\b;
        info = [];
        end
    end
end