function [estimate,info] = ml_estimate(A,P,F,J,approx,options)
% ML_ESTIMATE implementation of the multilevel error
% estimator on the Anorm of the error.
% See P. Vacek, J. Papez and Z. Strakos,
% A posteriori error estimates based on multilevel 
% decompositions with large problems on the coarsest level.
% 
% Inputs:
% A - cell contaning system matrices
% P - cell contaning prolongation matrices
% F - right hand side vector on the finest level
% numberOfLevels
% approx - approximation of A{J}^(-1)*F
% options
%   .coarseSolveName = 'backslash' compute the term coresponding to the
%       coarsest level using the MATLAB backslash operator
%   .coarseSolveName = 'CGAdaptive' compute the term coresponding to the coarsest
%       level using cg with a procedure for stopping described in the paper.
%       This variant requires bound on the smallest eigenvalue 
%       of the problem on the coarsest level
%       .coarseSolveParameters.smallestEigenvalueLowerBound
%       and parameter coarseSolveParameters.ratio 
%   .coarseSolveName = 'none' the coarsest level term is not included in
%       the estimate
%   .coarseSolveName = 'diag' compute the term coresponding to the
%       coarsest level by replacing the coarsest-level matrix with
%       its diagonal, the term is multiplied by .coarseSolveParameters.constant 
%   .coarseSolveName = 'CGFixNumberOfIteration' compute the term coresponding to the
%       coarsest level using .coarseSolveParameters.numberOfIterations
%   .summationStyle  = 'l1' or 'l2',  the terms corresponding to the individual level are summed
%       as l1 or l2 vector norms
% Output:
% estimate of the A-norm of the error on the finest level
% info structure containing .numberOfIterations if 'CGAdaptive' is used

arguments 
    A
    P
    F
    J
    approx
    options.summationStyle {mustBeMember(options.summationStyle,{'l1','l2'})} = 'l2';
    options.coarseSolveName {mustBeMember(options.coarseSolveName,{'backslash','diag', ...
                                                    'CGAdaptive','CGFixNumberOfIterations'})} = 'backslash';
    options.coarseSolveParameters = [];
end

info = [];

% finest level
R{J} = F - A{J}*approx;
z = diag(diag(A{J})) \ R{J};
term = R{J}'*z;

switch (options.summationStyle)
    case 'l1'
        sum = sqrt(term);
    case 'l2'
        sum = term;
end

% coarse levels
for j=J-1:-1:2
    R{j} = P{j+1}'*R{j+1};
    z = diag(diag(A{j})) \ R{j};
    term = R{j}'*z;
    switch (options.summationStyle)
        case 'l1'
            sum = sum + sqrt(term);
        case 'l2'
            sum = sum + term;
    end
end

% coarsest level
R{1} = P{2}'*R{2};

switch(options.coarseSolveName)
    case 'backslash'
        z = A{1}\R{1};
        term = R{1}'*z;
    case 'CGAdaptive'
        switch (options.summationStyle)
        case 'l1'
            error('Combination of summationStyle=''l1'' and coarseSolveName=''CGAdaptive'' is not possible.')
        case 'l2'
            [term,info] = cg_coarse_term_estimate(A{1},R{1},zeros(size(R{1})),sum,options.coarseSolveParameters);
        end
    case 'diag'
        z = diag(diag(A{1})) \ R{1};  
        term = options.coarseSolveParameters.constant*R{1}'*z;
    case 'CGFixNumberOfIterations'
        stoppingCriterion.name = 'fixNumberOfIterations';
        stoppingCriterion.numberOfIterations = options.coarseSolveParameters.numberOfIterations;
        z = cg(A{1},R{1},zeros(size(R{1})),stoppingCriterion);
        term = R{1}'*z;
    case 'none'
        term = 0;
end

switch (options.summationStyle)
        case 'l1'
            sum = sum + sqrt(term);
            estimate = sum;
        case 'l2'
            sum = sum + term;
            estimate = sqrt(sum);
end

end