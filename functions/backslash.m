function [approx,info] = backslash(A,b)
% BACKSHLASH call Matlab backslash operator
%
% Inputs:
% A - matrix
% b - right-hand side vector
%
% Outputs:
% approx - computed approximation of A^(-1)*b
% info - structure - contains
%   res2norm - Euclidian norm of residual b-A*approx
%   numberOfIterations  = 1

approx = A\b;
info.res2norm = norm(b-A*approx);
info.numberOfIterations = 1;
end