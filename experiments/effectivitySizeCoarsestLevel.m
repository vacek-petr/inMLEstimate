clear all; close all;
addpath('..\functions');
addpath('..\classes');
if ~exist('..\results', 'dir')
    mkdir('..\results')
end

% "2Dpeak", "2Dshell", "3Dpeak", "3Dshell",
problemName = "2Dpeak"; 

disp("Solving " + problemName)
addpath('..\data\' + problemName)
% load multigrid hierarchy data stored in mh an instance of class MultilevelHierarchy
load(problemName + ".mat",'mh');

% vcycle setting
switch(problemName)
    case {"2Dpeak", "2Dshell"}
        J = 3; % vcycleNumberOfLevels
    case {"3Dpeak", "3Dshell"}
        J = 2; % vcycleNumberOfLevels
end

vcycleRelativeTolerance = 1e-11;
vcycleMaximumNumberOfIterations = 30;
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
% and also in the experiment effectivityNumberOfLevels.
switch(problemName)
    case "2Dpeak"
        C_numexp = 1.178301355536950;
        diamOmega = sqrt(2);
        meshsize_0 = 0.12856486930664487;
    case"2Dshell"
        C_numexp = 1.386793086895968;
        diamOmega = 1;
        meshsize_0 =  0.09771455584161245;
    case"3Dpeak"
        C_numexp = 1.282709081580298;
        diamOmega = sqrt(3);
        meshsize_0 = 0.2886751345948128;
    case "3Dshell"
        C_numexp = 1.842461635385764;
        diamOmega = 1;
        meshsize_0 = 0.13832808829250579;
end

% 2Dpeak: 1/0.84; 2Dshell: 1/0.72; 3Dshell: 1/0.49; %1.193071412121244;

% loop over different setting - increasing the size of the coarsest-level problem 
for j = 1:mh.numberOfLevels - (J-1)

    disp("setting " + num2str(j) + "/" + num2str(mh.numberOfLevels - (J-1)))
    % select stiffness and prolongation matrices, the rhs vector and
    % smallest eigenvalue estimate
    A = mh.A(j:j+J-1);
    P = mh.P(j:j+J-1);
    F = mh.F{j+J-1};
    smallestEigenvalueLowerBound = mh.ASmallestEigenvalues{j};
    solution = mh.solution{j+J-1};
    
    iter = 0;
    approxVcycle = zeros(size(F));
    errAnormInit = sqrt((approxVcycle-solution)'*A{J}*(approxVcycle-solution));
    errAnorm = errAnormInit;
    errAnormPrevisous = errAnorm;
    disp("Vcycle iter. " + num2str(iter) + ", errorAnorm: " + num2str(errAnorm))

    while (errAnorm/errAnormInit>vcycleRelativeTolerance)&&(iter<vcycleMaximumNumberOfIterations)
        iter = iter + 1;
        [approxVcycle,coarsestLevelSolverInfo] = vcycle(A,P,J,F,approxVcycle,smoother,coarsestLevelSolver);
        
        errAnorm = sqrt((approxVcycle-solution)'*A{J}*(approxVcycle-solution));
        disp("Vcycle iter. " + num2str(iter) + ...
            ", coarse-level solver iter. " + num2str(coarsestLevelSolverInfo.numberOfIterations) + ...
            ", errorAnorm " + num2str(errAnorm) + ...
            ", rateAnorm: " + num2str(errAnorm/errAnormPrevisous))
        errAnormPrevisous = errAnorm;

        
        estimate = ml_estimate(A,P,F,J,approxVcycle,coarseSolveName='backslash');
        efficiencyIndexBackslash(iter,j) = C_numexp*estimate/errAnorm;
        
        estimate = ml_estimate(A,P,F,J,approxVcycle, ...
            coarseSolveName='CGFixNumberOfIterations', ...
            coarseSolveParameters=struct('numberOfIterations',4));
        efficiencyIndexCGFixNumberOfIterations(iter,j) = C_numexp*estimate/errAnorm;

        meshsize = meshsize_0/(2^(j-1));
        coarseConstant = (diamOmega/meshsize)^2;
       
        estimate = ml_estimate(A,P,F,J,approxVcycle,coarseSolveName='diag', ...
            coarseSolveParameters=struct('constant',coarseConstant));
        efficiencyIndexDiag(iter,j) = C_numexp*estimate/errAnorm;

        
        [estimate,info] = ml_estimate(A,P,F,J,approxVcycle,coarseSolveName='CGAdaptive', ...
            coarseSolveParameters=struct('smallestEigenvalueLowerBound',smallestEigenvalueLowerBound, ...
                                         'ratio',0.1));
        efficiencyIndexCGAdaptive(iter,j) = C_numexp*estimate/errAnorm;
        numberOfCGIterationCGAdaptive(iter,j) = info.numberOfIterations;
    end

    coarsestLevelDoFs = size(A{1},1);

    f1 = figure(11);
    subplot(2,2,1)
    semilogx(coarsestLevelDoFs*ones(iter,1),efficiencyIndexBackslash(1:iter,j), ...
        LineStyle='none', ...
        Marker='x', ...
        MarkerEdgeColor='#7ad151', ...
        MarkerSize=8, ...
        LineWidth=2)
    hold on
    ylabel('Efficiency index')
    title('direct solver')
    xlabel('# coarsest-level DoFs')

    subplot(2,2,2)
    semilogx(coarsestLevelDoFs*ones(iter,1),efficiencyIndexCGFixNumberOfIterations(1:iter,j), ...
        LineStyle='none', ...
        Marker='x', ...
        MarkerEdgeColor='#2a788e', ...
        MarkerSize=8, ...
        LineWidth=2)
    hold on
    ylabel('Efficiency index')
    title( '4 iterations of CG')
    xlabel('# coarsest-level DoFs')

    subplot(2,2,3)
    semilogx(coarsestLevelDoFs*ones(iter,1),efficiencyIndexDiag(1:iter,j), ...
        LineStyle='none', ...
        Marker='x', ...
        MarkerEdgeColor='#414487', ...
        MarkerSize=8, ...
        LineWidth=2)
    hold on
    ylabel('Efficiency index')
    title('diagonal')
    xlabel('# coarsest-level DoFs')

    subplot(2,2,4)
    semilogx(coarsestLevelDoFs*ones(iter,1),efficiencyIndexCGAdaptive(1:iter,j), ...
        LineStyle='none', ...
        Marker='x', ...
        MarkerEdgeColor='#22a884', ...
        MarkerSize=8, ...
        LineWidth=2)
    ylabel('Efficiency index')
    xlabel('# coarsest-level DoFs')
    title('adaptive CG approximation')
    hold on

    f2 = figure(15);
    plot(1:iter,numberOfCGIterationCGAdaptive(1:iter,j), ...
        Marker='x', ...
        MarkerSize=8, ...
        LineWidth= 2)
    xlabel('V-cycle iterations')
    ylabel('Number of CG iterations')
    title('adaptive CG approximation')
    hold on
end
fileName = "..\results\effectivitySizeCoaresetLevel_" + problemName;
saveas( f1, fileName , "png" )
saveas( f1, fileName , "fig" )

fileName = "..\results\effectivitySizeCoaresetLevelCGAdaptiveIter_" + problemName;
saveas( f2, fileName , "png" )
saveas( f2, fileName , "fig" )

results.problemName = problemName;
results.efficiencyIndexBackslash = efficiencyIndexBackslash;
results.efficiencyIndexCGAdaptive = efficiencyIndexCGAdaptive;
results.numberOfCGIterationCGAdaptive = numberOfCGIterationCGAdaptive;
results.efficiencyIndexCGFixNumberOfIterations = efficiencyIndexCGFixNumberOfIterations;
results.efficiencyIndexDiag = efficiencyIndexDiag;

fileName = "..\results\effectivitySizeCoaresetLevelData_" + problemName;
save(fileName , "results" )
