classdef BackSlashSolver < Solver
    %BACKSLASHSOLVER Summary of this class goes here
    %   Detailed explanation goes here

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