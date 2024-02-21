% script for saving matrices generated in Fenics as MatrixHierarchy class

% first load "A.mat", "F.mat", "P.mat", and "BN.mat" generated using
% assembleMatricesFenics.py

addpath('..\classes');
mh = MatrixHierarchy();
mh.name = '3Dpeak';
mh.A = A;
mh.numberOfLevels = size(mh.A,2);
mh.F = F;
mh.P = {};
mh.P(2:mh.numberOfLevels) = P;

% In FEniCS the stiffness matrix is assembled using all nodes of the mesh. The
% homogeneous Dirichlet boundary condition is then applied by setting to zero all non-
% diagonal elements in rows and columns which correspond to nodes on the boundary
% and setting to zero the corresponding elements in the right-hand side vector. We
% modify the stiffness matrices, the prolongation matrices and the right-hand side vector
% so that the Galerkin condition A{j-1} = P{j}'*A{j}*P{j} is satisfied.

for j = 1:mh.numberOfLevels
    % cut rows and collums coresponding to nodes on boundary
    mh.A{j}(BN{j},:) = [];
    mh.A{j}(:,BN{j}) = [];
    mh.F{j}(BN{j}) = [];
    mh.F{j} = mh.F{j}';
    if (j ~= 1)
        mh.P{j}(:,BN{j-1}) = [];
        mh.P{j}(BN{j},:) = [];
        mh.P{j} = mh.P{j}.*(mh.P{j} ~= 0); % drop zero entries from sparse matrix
    end
end

% compute solution using MATLAB backslash on levels 1 to 6
for j = 1:6
    mh.solution{j} = mh.A{j}\mh.F{j};
end

% compute solution using V-cycle on level 7
smoother = Smoother('gs',3,3);
coarsestLevelSolver = Solver('backslash');
approx = zeros(size(mh.F{7}));
for i = 1:30
        approx = vcycle(mh.A,mh.P,7,mh.F{7},approx,smoother,coarsestLevelSolver);
end
mh.solution{7} = approx;

% approximate smallest eigenvalue of A on levels 1 to 5
for j = 1:5
    mh.ASmallestEigenvalues{j} = eigs(A{j},1,'smallestabs');
end
% approximate smallest eigenvalue of A on levels 6 and 7
% we know that ther is a lower bound on the smallest eigenvalue, Ch^3,
% where C is a constant so we extrapolate the eigenvalues as 
mh.ASmallestEigenvalues{6} = mh.ASmallestEigenvalues{5}*0.5^3;
mh.ASmallestEigenvalues{7} = mh.ASmallestEigenvalues{6}*0.5^3;
save(mh.name,"mh")