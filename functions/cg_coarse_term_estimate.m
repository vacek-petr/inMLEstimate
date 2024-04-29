function [estimate,info] = cg_coarse_term_estimate(A,b,approx,fineLevelTermsL2normSq,parameters)
% CG_COARSE_TERM_ESTIMATE Implementation of the CG method for computing the approximation
% of b'*A^(-1)*b, in the context of multilevel estimator described in 
% P. Vacek, J. Papez and Z. Strakos, A posteriori error estimates based on 
% multilevel decompositions with large problems on the coarsest level.
%
% The computation is stopped using the estimate on the A-norm of the error based on the Gauss-Radau
% quadrature (Meurant,Tichy 2023, second inequality in (3.5) with updating
% formula for a coefficient (3.3)]), such that 
% ||error||^2_A <= (fineLevelTermsL2norm^2 + mu^2)*ratio
% 
% Inputs:
% A - symmetric positive definite matrix
% b - right-hand side vector
% approx - initial approximation of A^(-1)*b
% parameters - structure containing
%   .smallestEigenvalueLowerBound - lower bound of the smallest eigenvalue of A
%   .ratio
% fineLevelTermsL2normSq of the MLEtimator
%   
% Outputs:
% bound  - upper bound of the coarsest-level term in MLEstimator
% info - structure containing
%   .numberOfIterations

iter = 0;
r = b - A*approx;
rrNew = r'*r;
p = r;

smallestEigenvalueLowerBound = parameters.smallestEigenvalueLowerBound;
gamma = 1/smallestEigenvalueLowerBound;
ratio = parameters.ratio;
muSq = 0;
zetaSq = gamma*rrNew;
if(zetaSq > ratio*fineLevelTermsL2normSq)
    stoppingCriterionSatisfied = 0;
else
    stoppingCriterionSatisfied = 1;
    estimate = muSq + zetaSq;
end

while(~stoppingCriterionSatisfied)
    rr = rrNew;
    iter = iter + 1;
    Ap = A*p;
    alpha = rr/(p'*Ap);
    approx = approx + alpha*p;
    r = r - alpha*Ap;
    rrNew = r'*r;
    delta = rrNew/rr;
    p = r + delta*p;

    gamma = (gamma-alpha)/(smallestEigenvalueLowerBound*(gamma-alpha) + delta);
    zetaSq = gamma*rrNew;
    muSq = muSq + alpha*rr;
    if (zetaSq > ratio*(fineLevelTermsL2normSq + muSq))
        stoppingCriterionSatisfied = 0;
    else
        stoppingCriterionSatisfied = 1;
        estimate = muSq + zetaSq;
    end
end
info.numberOfIterations = iter;
end
