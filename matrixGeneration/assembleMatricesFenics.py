import sys
import petsc4py
import time
petsc4py.init(sys.argv)
from petsc4py import PETSc
import numpy as np
from scipy import sparse
import scipy.io
from dolfin import *
from mshr import *

# definition of the problem
def operator(u,v):
    a = inner(grad(u), grad(v))*dx
    return a

def rightside(f,v):
    b = f*v*dx
    return b

def boundary(x, on_boundary):
    return on_boundary

n = 3 # number of cell in one direction
mesh = UnitCubeMesh(n,n,n)
exactSolution = Expression("x[0]*(x[0] - 1)*x[1]*(x[1] - 1)*x[2]*(x[2]-1)*exp(-100*(pow(x[0] - 0.5,2) + pow(x[1] - 0.5,2) + pow(x[2] - 0.5,2)))", degree=2)
gradExactSolution = Expression(["(200*pow(x[0],3) - 300*pow(x[0],2)+98*x[0]+1)*(x[1]-1)*x[1]*(x[2]-1)*x[2]*(-exp(-100*(pow(x[0] - 0.5,2) + pow(x[1] - 0.5,2) + pow(x[2] - 0.5,2))))","(200*pow(x[1],3) - 300*pow(x[1],2)+98*x[1]+1)*(x[0]-1)*x[0]*(x[2]-1)*x[2]*(-exp(-100*(pow(x[0] - 0.5,2) + pow(x[1] - 0.5,2) + pow(x[2] - 0.5,2))))","(200*pow(x[2],3) - 300*pow(x[2],2)+98*x[2]+1)*(x[0]-1)*x[0]*(x[1]-1)*x[1]*(-exp(-100*(pow(x[0] - 0.5,2) + pow(x[1] - 0.5,2) + pow(x[2] - 0.5,2))))"],degree=2)
f = Expression("-(2*(20000*pow(x[0],4)*(x[1] - 1)*x[1]*(x[2] - 1)*x[2] - 40000*pow(x[0],3)*(x[1] - 1)*x[1]*(x[2] - 1)*x[2] + pow(x[0],2)*(20000*pow(x[1],4)*(x[2] - 1)*x[2]-40000*pow(x[1],3)*(x[2] - 1)*x[2] + pow(x[1],2)*(20000*pow(x[2],4) -40000*pow(x[2],3) + 73500*pow(x[2],2) - 53500*x[2] - 99)+ x[1]*(-20000*pow(x[2],4) + 40000*pow(x[2],3) - 53500*pow(x[2],2) + 33500*x[2] + 99) - 99*(x[2] - 1)*x[2]) + x[0]*(-20000*pow(x[1],4)*(x[2] - 1)*x[2] + 40000*pow(x[1],3)*(x[2] - 1)*x[2] + pow(x[1],2)*(-20000*pow(x[2],4) + 40000*pow(x[2],3) - 53500*pow(x[2],2) + 33500*x[2] + 99) + x[1]*(20000*pow(x[2],4) - 40000*pow(x[2],3) + 33500*pow(x[2],2) - 13500*x[2] - 99) + 99*(x[2] - 1)*x[2]) - 99*(x[1] - 1)*x[1]*(x[2] - 1)*x[2])*exp(-25*(4*pow(x[0],2) - 4*x[0] + 4*pow(x[1],2) - 4*x[1] + 4*pow(x[2],2) - 4*x[2] + 3)))",degree=2)

# number of levels in multigrid hierarchy
nlevel = 7

matricesA = np.empty((nlevel,), dtype=object) # stiffness matrices
vectorsF = np.empty((nlevel,), dtype=object) # right-hand sides
matricesP = np.empty((nlevel-1,), dtype=object) # prolongation matrices
vectorsBoundaryNodes = np.empty((nlevel,), dtype=object) # vectors containing indexis of nodes on the boundary

tic = time.time()
for j in range(0,nlevel):
    P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
    VFine = FunctionSpace(mesh, P1)
    u = TrialFunction(VFine)
    v = TestFunction(VFine)
    bc = DirichletBC(VFine, 0.0, boundary)

    A, F = assemble_system(operator(u,v),rightside(f,v),bc)
    
    # solve on level j
    u = Function(VFine)
    solve(A, u.vector(),F,"cg")
    
    # save solution level j in xdmf file 
    fileName = "solution_level" + str(j) + ".xdmf"
    file = XDMFFile(fileName)
    file.write(u,0)
    
    # compute errors
    errorL2Norm = sqrt(assemble((u-exactSolution)**2*dx(mesh)))
    errorEnergyNorm = sqrt(assemble(inner(grad(u)-gradExactSolution,grad(u)-gradExactSolution)*dx(mesh)))

    print("level " + str(j) + ", number of DoFs " + str(VFine.dim()) + ", meshsize " + str(mesh.hmin()))
    print("L2 norm of the error = {0:16.8e}".format(errorL2Norm))
    print("Energy norm of the error = {0:16.8e}".format(errorEnergyNorm))
    print("-------------------------------------------")

    # put matrices and vectors into np arrays containing scipy sparse matrices
    APetsc = as_backend_type(A).mat()
    indptr, indices, data = APetsc.getValuesCSR()
    AScipy = sparse.csr_matrix((data,indices,indptr))
    matricesA[j] = AScipy
    
    FPetsc = as_backend_type(F).vec()
    vectorsF[j] = FPetsc 

    boundaryNodes = np.fromiter(bc.get_boundary_values().keys(), dtype=int)  
    boundaryNodes = boundaryNodes + np.ones(boundaryNodes.shape)
    boundaryNodes.sort()
    vectorsBoundaryNodes[j] = boundaryNodes

    if (j != 0):
        P = PETScDMCollection.create_transfer_matrix(VCoarse, VFine)
        PPetsc = P.mat()
        indptr, indices, data = PPetsc.getValuesCSR()
        PScipy = sparse.csr_matrix((data,indices,indptr))
        matricesP[j-1] = PScipy

    VCoarse = VFine
    mesh = refine(mesh)

scipy.io.savemat("A.mat", {"A": matricesA})
scipy.io.savemat("F.mat", {"F": vectorsF})
scipy.io.savemat("P.mat", {"P": matricesP})
scipy.io.savemat("BN.mat", {"BN": vectorsBoundaryNodes})

toc = time.time()
print("Total time", toc-tic)