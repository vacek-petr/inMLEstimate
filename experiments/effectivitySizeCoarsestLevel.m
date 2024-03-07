clc; clear all; close all;

addpath('..\functions');
addpath('..\classes');
addpath('..\data\3Dpeak');

vcycleRelativeTolerance = 10^(-11);
vcycleMaximumNumberOfIterations = 100;
smoother = Smoother('gs',3,3);
coarsestLevelSolver = Solver('cg');
coarsestLevelSolver.stoppingCriterion.name = 'res2norm';

coarsestLevelSolver.stoppingCriterion.relative = true;
coarsestLevelSolver.stoppingCriterion.tolerance = 0.1;

numberOfLevels = 2;

% The constant C is chosen as a minimal value such that the efficiency 
% index for the variant with the direct solver on the coarsest level 
% is above or equal to one for all iterates.
C = 1.100096470053534;

for j=2:7-(numberOfLevels-1)
    disp("setting "+num2str(j-1)+"/5")
    load('3Dpeak.mat','mh');
    mh.selectLevels(numberOfLevels=numberOfLevels,from=j);

    efficiencyIndexBackslash = [];
    efficiencyIndexCGFixNumberOfIterations = [];
    efficiencyIndexDiag = [];
    efficiencyIndexNewProcedure= [];
    newProcedureCGNumberOfIteration = [];

    iter = 0;
    approx = zeros(size(mh.F{numberOfLevels}));
    errAnormInit = sqrt((approx-mh.solution{numberOfLevels})'*mh.A{numberOfLevels}*(approx-mh.solution{numberOfLevels}));
    errAnorm = errAnormInit;

    while (errAnorm/errAnormInit>vcycleRelativeTolerance)&&(iter<vcycleMaximumNumberOfIterations)
        iter = iter + 1;
        disp("V-cycle iteration "+num2str(iter))
        approx = vcycle(mh.A,mh.P,numberOfLevels,mh.F{numberOfLevels},approx,smoother,coarsestLevelSolver);
        errAnorm = sqrt((approx-mh.solution{numberOfLevels})'*mh.A{numberOfLevels}*(approx-mh.solution{numberOfLevels}));

        efficiencyIndexBackslash(iter) = C*computeErrAnormMLEstimate(mh.A,mh.P,mh.F{numberOfLevels},numberOfLevels,approx,coarseSolveName='backslash')/errAnorm;

        coarseSolveParameters = [];
        coarseSolveParameters.numberOfIterations = 4;
        efficiencyIndexCGFixNumberOfIterations(iter) = C*computeErrAnormMLEstimate(mh.A,mh.P,mh.F{numberOfLevels},numberOfLevels,approx,coarseSolveName='CGFixNumberOfIterations',coarseSolveParameters=coarseSolveParameters)/errAnorm;
 
        coarseSolveParameters = [];
        coarseSolveParameters.constant = 3*(1.7321*2^(j-1))^2;
        efficiencyIndexDiag(iter) = C*computeErrAnormMLEstimate(mh.A,mh.P,mh.F{numberOfLevels},numberOfLevels,approx,coarseSolveName='diag', coarseSolveParameters=coarseSolveParameters)/errAnorm;

        coarseSolveParameters = [];
        coarseSolveParameters.smallestEigenvalueLowerBound = mh.ASmallestEigenvalues{1};
        coarseSolveParameters.ratio = 0.1;
        [errAnormMLEstimateNewProcedure,info] = computeErrAnormMLEstimate(mh.A,mh.P,mh.F{numberOfLevels},numberOfLevels,approx,coarseSolveName='CGNewProcedure',coarseSolveParameters=coarseSolveParameters);
        efficiencyIndexNewProcedure(iter) = C*errAnormMLEstimateNewProcedure/errAnorm;
        newProcedureCGNumberOfIteration(iter) = info.numberOfIterations;
    end

    figure(11)
    plot(j*ones(length(efficiencyIndexBackslash),1),efficiencyIndexBackslash,'LineStyle', 'none',Marker='o',MarkerEdgeColor='g')
    hold on
    ylabel('Efficiency index')
    title('direct solve')

    figure(12)
    plot(j*ones(length(efficiencyIndexCGFixNumberOfIterations),1),efficiencyIndexCGFixNumberOfIterations,'LineStyle', 'none',Marker='+',MarkerEdgeColor='b')
    hold on
    ylabel('Efficiency index')
    title('CG 4 iter')

    figure(13)
    plot(j*ones(length(efficiencyIndexDiag),1),efficiencyIndexDiag,'LineStyle', 'none',Marker='*',MarkerEdgeColor='b')
    hold on
    ylabel('Efficiency index')
    title('diagonal')

    figure(14)
    plot(j*ones(length(efficiencyIndexNewProcedure),1),efficiencyIndexNewProcedure,'LineStyle', 'none',Marker='x', MarkerEdgeColor='r')
    ylabel('Efficiency index')
    title('new procedure')
    hold on

    figure(15)
    plot(1:length(newProcedureCGNumberOfIteration),newProcedureCGNumberOfIteration,'LineWidth', 1,'DisplayName', 'new pocedure',Marker='x');
    xlabel('V-cycle iterations')
    ylabel('Number of CG iter')
    hold on

    efficiencyIndexBackslashTable{j} = efficiencyIndexBackslash';
    efficiencyIndexCGFixNumberOfIterationsTable{j} = efficiencyIndexCGFixNumberOfIterations';
    efficiencyIndexDiagTable{j} = efficiencyIndexDiag';
    efficiencyIndexNewProcedureTable{j} = efficiencyIndexNewProcedure';
    newProcedureCGNumberOfIterationTable{j} = newProcedureCGNumberOfIteration';
end
