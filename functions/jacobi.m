function approx  = jacobi(A,b,approx,numberOfIterations,dampingParameter)
% JACOBI Implementation of the (damped) Jacobi method
%
% Inputs:
% A - matrix
% b - right-hand side vector
% approx - intial approximation of A^(-1)*b
% numberOfIterations
% dampingParameter
%
% Outputs:
% approx - computed approximation to A^(-1)*b

M = (1/dampingParameter).*diag(diag(A));
N = diag(diag(A)) - A;


for i = 1:numberOfIterations
    approx   = M\(N*approx + b) + (1-dampingParameter)*approx;
end

end