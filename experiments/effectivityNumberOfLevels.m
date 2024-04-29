clear all; close all;
addpath('..\functions');
addpath('..\classes');
if ~exist('..\results', 'dir')
    mkdir('..\results')
end

% select problem name from "2Dpeak", "2Dshell", "3Dpeak", "3Dshell"
problemName = "3Dpeak";

disp("Solving " + problemName)
addpath('..\data\' + problemName);
% load multigrid hierarchy data stored in mh an instance of class MultilevelHierarchy
load(problemName + ".mat",'mh');

% V-cycle setting
vcycleRelativeTolerance = 1e-11;
vcycleMaximumNumberOfIterations = 1000;
smoother = GSSmoother(3,3);
% define coarsest level solver
stoppingCriterion.name = 'res2norm';
stoppingCriterion.relative = true;
switch(problemName)
    case {"2Dpeak", "2Dshell", "3Dpeak"}
       stoppingCriterion.tolerance = 0.1;
    case  "3Dshell"
        stoppingCriterion.tolerance = 0.5;
end
coarsestLevelSolver = CGSolver(stoppingCriterion);

% The constant C_num_exp is chosen as a minimal value such that the efficiency 
% index for the variant with the direct solver on the coarsest level 
% is above or equal to one for all iterates in this experiment 
% and also in the experiment effectivitySizeCoarsestLevel.
switch(problemName)
    case "2Dpeak"
        C_numexp = 1.178301355536950;
    case"2Dshell"
        C_numexp = 1.386793086895968;
    case"3Dpeak"
        C_numexp = 1.282709081580298;
    case "3Dshell"
        C_numexp = 1.84246163538;5764;
end

% loop over different setting, J=vcycleNumberOfLevels
for J = 2:mh.numberOfLevels
    disp("setting " + num2str(J-1) + "/" + num2str(mh.numberOfLevels-1))
    
    % select stiffness and prolongation matrices, the rhs vector and
    % smallest eigenvalue estimate
    A = mh.A(1:J);
    P = mh.P(1:J);
    F = mh.F{J};
    solution = mh.solution{J};
    
    iter = 0;
    approxVcycle = zeros(size(F));
    errAnormInit = sqrt((approxVcycle-solution)'*A{J}*(approxVcycle-solution));
    errAnorm = errAnormInit;
    errAnormPrevisous = errAnorm;
    disp("Vcycle iter. " + num2str(iter) + ...
        ", errorAnorm: " + num2str(errAnorm))
    
    while (errAnorm/errAnormInit>vcycleRelativeTolerance)&&(iter<vcycleMaximumNumberOfIterations)
        iter = iter + 1;
        [approxVcycle,coarsestLevelSolverInfo] = vcycle(A,P,J,F,approxVcycle,smoother,coarsestLevelSolver);
        
        errAnorm = sqrt((approxVcycle-solution)'*A{J}*(approxVcycle-solution));
        disp("Vcycle iter. " + num2str(iter) + ...
            ", coarse-level solver iter. " + num2str(coarsestLevelSolverInfo.numberOfIterations) + ...
            ", errorAnorm " + num2str(errAnorm) + ...
            ", rateAnorm: " + num2str(errAnorm/errAnormPrevisous))
        errAnormPrevisous = errAnorm;

        estimate = ml_estimate(A,P,F,J,approxVcycle, ...
            coarseSolveName='backslash', ...
            summationStyle='l1');
        efficiencyIndexL1(iter,J) = C_numexp*estimate/errAnorm;
        
        estimate = ml_estimate(A,P,F,J,approxVcycle, ...
            coarseSolveName='backslash', ...
            summationStyle='l2');
        efficiencyIndexL2(iter,J) = C_numexp*estimate/errAnorm;    
    end
        
    f = figure(1);
    plot(J*ones(iter,1),efficiencyIndexL1(1:iter,J), ...
        LineStyle='none', ...
        Marker='x', ...
        MarkerEdgeColor='#28ae80', ...
        MarkerSize=8, ...
        LineWidth=2)    
    hold on
    plot(J*ones(iter,1),efficiencyIndexL2(1:iter,J), ...
        LineStyle='none', ...
        Marker='x', ...
        MarkerEdgeColor='#7ad151', ...
        MarkerSize=8, ...
        LineWidth=2)    
    hold on
    ylabel('Efficiency index')
    xlabel('number of levels')
    title('Dependence on the number of levels')
end

fileName = "..\results\effectivityNumberOfLevels_" + problemName;
saveas( f, fileName , "png" )
saveas( f, fileName , "fig" )

results.problemName = problemName;
results.efficiencyIndexL1 = efficiencyIndexL1;
results.efficiencyIndexL2 = efficiencyIndexL2;

fileName = "..\results\effectivityNumberOfLevelsData_" + problemName;
save(fileName , "results" )
