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
coarsestLevelSolver.stoppingCriterion.tolerance = 1e-1;

% The constant C is chosen as a minimal value such that the efficiency 
% index is above or equal to one for all iterates.
C = 1.193071412121244;

for j = 2:6
    disp("setting "+num2str(j-1)+"/5")
    numberOfLevels = j;
    load('3Dpeak.mat','mh');
    mh.selectLevels(numberOfLevels=j,from=2);

    efficiencyIndexBackslashL1 = [];
    efficiencyIndexBackslashL2 = [];
    
    iter = 0;
    approx = zeros(size(mh.F{numberOfLevels}));
    errAnormInit = sqrt((approx-mh.solution{numberOfLevels})'*mh.A{numberOfLevels}*(approx-mh.solution{numberOfLevels}));
    errAnorm = errAnormInit;
    
    while (errAnorm/errAnormInit>vcycleRelativeTolerance)&&(iter<vcycleMaximumNumberOfIterations)
        iter = iter + 1;
        disp("V-cycle iteration "+num2str(iter))
        approx = vcycle(mh.A,mh.P,numberOfLevels,mh.F{numberOfLevels},approx,smoother,coarsestLevelSolver);        
        errAnorm = sqrt((approx-mh.solution{numberOfLevels})'*mh.A{numberOfLevels}*(approx-mh.solution{numberOfLevels}));
        efficiencyIndexBackslashL1(iter) = C*computeErrAnormMLEstimate(mh.A,mh.P,mh.F{numberOfLevels},numberOfLevels,approx,coarseSolveName='backslash',summationStyle='l1')/errAnorm;
        efficiencyIndexBackslashL2(iter) = C*computeErrAnormMLEstimate(mh.A,mh.P,mh.F{numberOfLevels},numberOfLevels,approx,coarseSolveName='backslash',summationStyle='l2')/errAnorm;    
    end
    
    figure(21)
    plot(j*ones(length(efficiencyIndexBackslashL1),1),efficiencyIndexBackslashL1,'LineStyle', 'none',Marker='o',MarkerEdgeColor='g')
    hold on
    plot(j*ones(length(efficiencyIndexBackslashL2),1),efficiencyIndexBackslashL2,'LineStyle', 'none',Marker='x',MarkerEdgeColor='b')
    ylabel('Efficiency index')
    xlabel('number of levels')
    title('l1 vs l2 summationStyle')
    
    efficiencyIndexBackslashL1Table{j} = efficiencyIndexBackslashL1';
    efficiencyIndexBackslashL2Table{j} = efficiencyIndexBackslashL2';
end
