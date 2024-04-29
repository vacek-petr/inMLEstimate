function [approx,coarseSolverInfo] = vcycle(A,P,j,RHS,approx,smoother,solver)
% VCYCLE Recursive implementation of the multigrid V-cycle method
%
% Inputs:
% A - cell array - stiffness matrices
% P - cell array - prolongation matrices
% j - number - indicating the level of the hierarchy
% RHS - vector - right hand side
% approx - vector - approximation of A{j}^{(-1)}RHS
% smoother - class Smoother
% solver - class Solver
%
% Outputs:
% approx - vector - approximation of A{j}^{(-1)}RHS
% coarseSolverInfo - structure containing
%   numberOfIterations
%   res2norm

if (j > 1)
    approx = smoother.apply(A{j},RHS,approx,j,'pre');
    [correction,coarseSolverInfo] = vcycle(A,P,j-1,P{j}'*(RHS-A{j}*approx),zeros(size(A{j-1},1),1),smoother,solver);
    approx = smoother.apply(A{j},RHS,approx+P{j}*correction,j,'post');
else
    [approx, coarseSolverInfo] = solver.apply(A{1},RHS,zeros(size(RHS)));
end
end