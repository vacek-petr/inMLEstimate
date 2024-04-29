function [approx,info] = cg(A,b,approx,stoppingCriterion)
% CG Implementation of the CG method (Hestenes, Stiefel 1952) with
% various stopping criteria.
%
% Inputs:
% A - symmetric positive definite matrix
% b - right-hand side vector
% approx - initial approximation of A^(-1)*b
% stoppingCriterion - structure
%   implemented variants:
%   .name='errAnorm' - stopping when (approximate) absolute A norm of the error
%                      is less than given .tolerance
%   .name='res2norm' - stopping when (relative [.relative=1] or absolute [.relative=0])
%                      Euclidien norm of residual is less than given .tolerance
%   .name='GR' - stopping criterion based on Gauss-Radau quadrature
%                interpretation of CG, see (Meurant,Tichy 2023, second inequality in (3.5) 
%                with updating formula for a coefficient (3.3)]),
%                requires .tolerance and .smallestEigenvalueLowerBound
%   .name='fixNumberOfIteration' - stop after .numberOfIterations 
%
% Outputs:
% approx - approximation of A^(-1)*b
% info - structure containing
%   res2norm - Euclidean norm of the residual b-A*approx
%   numberOfIterations - number of performed CG iterations


iter = 0;
r = b - A*approx;
rrNew = r'*r;
p = r;


switch(stoppingCriterion.name)
    case 'errAnorm'
        tolerance = stoppingCriterion.tolerance;
        approxBackslash = A\b;
        err_Anorm = sqrt((approx-approxBackslash)'*A*(approx-approxBackslash));
        if (err_Anorm > tolerance)
            stoppingCriterionSatisfied = 0;
        else
            stoppingCriterionSatisfied = 1;
        end
    case 'res2norm'
        if (stoppingCriterion.relative)
            tolerance = norm(b)*stoppingCriterion.tolerance;
        else
            tolerance = stoppingCriterion.tolerance;
        end
        if(sqrt(rrNew) > tolerance)
            stoppingCriterionSatisfied = 0;
        else
            stoppingCriterionSatisfied = 1;
        end
    case 'GR'
        tolerance = stoppingCriterion.tolerance;
        smallestEigenvalueLowerBound = stoppingCriterion.smallestEigenvalueLowerBound;
        gamma = 1/smallestEigenvalueLowerBound;
        UpperBoundGRSq = gamma*rrNew;
        if (UpperBoundGRSq < 0)
            error('GR bound failed!');
        end
        if(sqrt(UpperBoundGRSq) > tolerance)
            stoppingCriterionSatisfied = 0;
        else
            stoppingCriterionSatisfied = 1;
        end 
    case 'fixNumberOfIterations'
        stoppingCriterionSatisfied = 0;
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

    switch(stoppingCriterion.name)
        case 'errAnorm'
            err_Anorm = sqrt((approx-approxBackslash)'*A*(approx-approxBackslash));
            if (err_Anorm > tolerance)
                stoppingCriterionSatisfied = 0;
            else
                stoppingCriterionSatisfied = 1;
            end
        case 'res2norm'
            if(sqrt(rrNew) > tolerance)
                stoppingCriterionSatisfied = 0;
            else
                stoppingCriterionSatisfied = 1;
            end
        case 'GR'
            gamma = (gamma-alpha)/(smallestEigenvalueLowerBound*(gamma-alpha) + delta);
            UpperBoundGRSq = gamma*rrNew;
            if (UpperBoundGRSq < 0)
                error('GR bound failed.');
            end
            if (sqrt(UpperBoundGRSq) > tolerance)
                stoppingCriterionSatisfied = 0;
            else
                stoppingCriterionSatisfied = 1;
            end
        case 'fixNumberOfIterations'
            if (iter < stoppingCriterion.numberOfIterations)
                stoppingCriterionSatisfied = 0;
            else
                stoppingCriterionSatisfied = 1;
            end
    end
end
info.res2norm = sqrt(rrNew);
info.numberOfIterations = iter;
end