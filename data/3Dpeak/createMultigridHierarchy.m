% script for saving matrices generated in Fenics 
% to a MultigridHierarchy class containing
% .name
% .A - cell array stiffness matrices
% .P - cell array prolongation matrices
% .F - cell array right-hand sides
% .numberOfLevels interger
% .solution - cell array 
% .ASmallestEigenvalues - cell array

addpath('..\..\classes\');
addpath('..\..\functions\');

% load data generated using assemble_matrices.py
load("A.mat","A");
load("F.mat","F");
load("P.mat","P");

mh = MultigridHierarchy('3Dpeak');
mh.A = A;
clear("A");
mh.numberOfLevels = size(mh.A,2);
mh.F = F;
clear("F");
mh.P = {};
mh.P(2:mh.numberOfLevels) = P;
clear("P");

for j = 1:4
    mh.solution{j} = mh.A{j}\mh.F{j};
end

% compute solution using V-cycle on level mh.numberOfLevels-1
smoother = GSSmoother(3,3);
coarsestLevelSolver = BackSlashSolver();
for J = 5:mh.numberOfLevels
    approx = zeros(size(mh.F{J}));
    for i = 1:30
        approx = vcycle(mh.A(1:J),mh.P(1:J),J,mh.F{J},approx,smoother,coarsestLevelSolver);
        disp("V-cycle iteration " + num2str(i) + ...
            ", res2norm: " + num2str(norm(mh.F{J} - mh.A{J}*approx)))
    end
    mh.solution{J} = approx;
end



% approximate smallest eigenvalue of A on levels 1 to mh.numberOfLevels-2
for j = 1:mh.numberOfLevels-2
    mh.ASmallestEigenvalues{j} = eigs(mh.A{j},1,'smallestabs');
end
% approximate smallest eigenvalue of A on levels mh.numberOfLevels-1 and mh.numberOfLevels
% we know that there is a lower bound on the smallest eigenvalue, Ch^3,
% where C is a constant so we extrapolate the eigenvalues as 
mh.ASmallestEigenvalues{mh.numberOfLevels-1} = mh.ASmallestEigenvalues{mh.numberOfLevels-2}*0.5^3;
mh.ASmallestEigenvalues{mh.numberOfLevels} = mh.ASmallestEigenvalues{mh.numberOfLevels-1}*0.5^3;
save(mh.name,"mh")