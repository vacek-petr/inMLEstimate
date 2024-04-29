function approx = gs(A,b,approx,numberOfIterations)
% GS Implementation of the Gauss-Seidel method
%
% Inputs:
% A - matrix
% b - right-hand side vector
% approx - intial approximation of A^(-1)*b
% numberOfIterations - number of iterations to be performed
%
% Outputs:
% approx - computed approximation to A^(-1)*b

U = triu(A,1);
L = tril(A);

for i = 1:numberOfIterations
    approx = L \ (b - U*approx);
end
end