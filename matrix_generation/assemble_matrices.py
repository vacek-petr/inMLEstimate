import petsc4py
import time
from petsc4py import PETSc
import numpy as np
from scipy import sparse
import scipy.io
from dolfin import *
from mshr import *

# Select problem name: "2Dpeak", "2Dshell", "3Dpeak", or "3Dshell" 
problemName = "2Dpeak" 
onlyFreeNodes = True
#  In FEniCS the stiffness matrix is assembled using all nodes of the mesh. The
#  homogeneous Dirichlet boundary condition is then applied by setting to zero all non-
#  diagonal elements in rows and columns which correspond to nodes on the boundary
#  and setting to zero the corresponding elements in the right-hand side vector. 
#  if (onlyFreeNodes == True)
#  we modify the stiffness matrices, the prolongation matrices and the right-hand side vector
#  so that the Galerkin condition A{j-1} = P{j}'*A{j}*P{j} holds.

if ( problemName == "2Dpeak" ):
    n = 11
    nlevel = 9
    mesh = UnitSquareMesh(n,n)
    exactSolution = Expression("x[0]*(x[0] - 1)*x[1]*(x[1] - 1)*exp(-100*(pow(x[0] - 0.5,2) + pow(x[1] - 0.5,2)))", degree=2)
    gradExactSolution = Expression(["(200*pow(x[0],3) - 300*pow(x[0],2)+98*x[0]+1)*(x[1]-1)*x[1]*(-exp(-100*(pow(x[0] - 0.5,2) + pow(x[1] - 0.5,2) )))","(200*pow(x[1],3) - 300*pow(x[1],2)+98*x[1]+1)*(x[0]-1)*x[0]*(-exp(-100*(pow(x[0] - 0.5,2) + pow(x[1] - 0.5,2))))"],degree=2)
    f = Expression("-(2*exp(-100*(pow(x[0] - 0.5,2) + pow(x[1] - 0.5,2))))*(20000*pow(x[0],4)*(x[1] - 1)*x[1] - 40000*pow(x[0],3)*(x[1] - 1)*x[1] + pow(x[0],2)*(20000*pow(x[1],4)-40000*pow(x[1],3) + 49000*pow(x[1],2)-29000*x[1] - 99) + x[0]*(-20000*pow(x[1],4) + 40000*pow(x[1],3) -29000*pow(x[1],2) + 9000*x[1] + 99)-99*(x[1]-1)*x[1])",degree=4)

if ( problemName == "2Dshell" ):
    n = 11
    nlevel = 9
    domain = Circle(Point(0,0),1) - Circle(Point(0,0),0.5) 
    mesh = generate_mesh(domain,n)
    f = Expression("0.0",degree=2)

    
if ( problemName == "3Dshell" ):
    nlevel = 5
    mesh_file_name = "3Dshell_meshfile"
    mesh = Mesh(mesh_file_name + ".xml")
    mesh = refine(mesh)
    f = Expression("0.0",degree=2)

if ( problemName == "3Dpeak" ):
    nlevel = 6
    n = 6 
    mesh = UnitCubeMesh(n,n,n)
    exactSolution = Expression("x[0]*(x[0] - 1)*x[1]*(x[1] - 1)*x[2]*(x[2]-1)*exp(-100*(pow(x[0] - 0.5,2) + pow(x[1] - 0.5,2) + pow(x[2] - 0.5,2)))", degree=2)
    gradExactSolution = Expression(["(200*pow(x[0],3) - 300*pow(x[0],2)+98*x[0]+1)*(x[1]-1)*x[1]*(x[2]-1)*x[2]*(-exp(-100*(pow(x[0] - 0.5,2) + pow(x[1] - 0.5,2) + pow(x[2] - 0.5,2))))","(200*pow(x[1],3) - 300*pow(x[1],2)+98*x[1]+1)*(x[0]-1)*x[0]*(x[2]-1)*x[2]*(-exp(-100*(pow(x[0] - 0.5,2) + pow(x[1] - 0.5,2) + pow(x[2] - 0.5,2))))","(200*pow(x[2],3) - 300*pow(x[2],2)+98*x[2]+1)*(x[0]-1)*x[0]*(x[1]-1)*x[1]*(-exp(-100*(pow(x[0] - 0.5,2) + pow(x[1] - 0.5,2) + pow(x[2] - 0.5,2))))"],degree=2)
    f = Expression("-(2*(20000*pow(x[0],4)*(x[1] - 1)*x[1]*(x[2] - 1)*x[2] - 40000*pow(x[0],3)*(x[1] - 1)*x[1]*(x[2] - 1)*x[2] + pow(x[0],2)*(20000*pow(x[1],4)*(x[2] - 1)*x[2]-40000*pow(x[1],3)*(x[2] - 1)*x[2] + pow(x[1],2)*(20000*pow(x[2],4) -40000*pow(x[2],3) + 73500*pow(x[2],2) - 53500*x[2] - 99)+ x[1]*(-20000*pow(x[2],4) + 40000*pow(x[2],3) - 53500*pow(x[2],2) + 33500*x[2] + 99) - 99*(x[2] - 1)*x[2]) + x[0]*(-20000*pow(x[1],4)*(x[2] - 1)*x[2] + 40000*pow(x[1],3)*(x[2] - 1)*x[2] + pow(x[1],2)*(-20000*pow(x[2],4) + 40000*pow(x[2],3) - 53500*pow(x[2],2) + 33500*x[2] + 99) + x[1]*(20000*pow(x[2],4) - 40000*pow(x[2],3) + 33500*pow(x[2],2) - 13500*x[2] - 99) + 99*(x[2] - 1)*x[2]) - 99*(x[1] - 1)*x[1]*(x[2] - 1)*x[2])*exp(-25*(4*pow(x[0],2) - 4*x[0] + 4*pow(x[1],2) - 4*x[1] + 4*pow(x[2],2) - 4*x[2] + 3)))",degree=2)

def operator(u,v):
    a = inner(grad(u), grad(v))*dx
    return a

def rightside(f,v):
    b = f*v*dx
    return b

def boundary(x, on_boundary):
    return on_boundary

def boundary_inside3D(x, on_boundary):
    return on_boundary and (np.sqrt(x[0]**2 + x[1]**2 + x[2]**2) < 0.9)

def boundary_outside3D(x, on_boundary):
    return on_boundary and (np.sqrt(x[0]**2 + x[1]**2+ x[2]**2) > 0.9)

def boundary_inside2D(x, on_boundary):
    return on_boundary and (np.sqrt(x[0]**2 + x[1]**2 ) < 0.9)

def boundary_outside2D(x, on_boundary):
    return on_boundary and (np.sqrt(x[0]**2 + x[1]**2) > 0.9)

matricesA = np.empty((nlevel,), dtype=object) # stiffness matrices
vectorsF = np.empty((nlevel,), dtype=object) # right-hand sides
matricesP = np.empty((nlevel-1,), dtype=object) # prolongation matrices
vectorsFreeNodes = np.empty((nlevel,), dtype=object) # vectors containing indexis of nodes on the boundary


tic = time.time()

for j in range(0,nlevel):
    P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
    VFine = FunctionSpace(mesh, P1)
    u = TrialFunction(VFine)
    v = TestFunction(VFine)

    if (problemName=="2Dpeak" or problemName=="3Dpeak"):
        bc = DirichletBC(VFine, 0.0, boundary)
    
    if (problemName=="2Dshell"):
        bc1 = DirichletBC(VFine, 100.0, boundary_inside2D)
        bc2 = DirichletBC(VFine, -10.0, boundary_outside2D)
        bc = [bc1,bc2]

    if (problemName=="3Dshell"):
        bc1 = DirichletBC(VFine, 100.0, boundary_inside3D)
        bc2 = DirichletBC(VFine, -10.0, boundary_outside3D)
        bc = [bc1,bc2]

    A, F = assemble_system(operator(u,v),rightside(f,v),bc)
    
    # solve on level j
    u = Function(VFine)
    solve(A, u.vector(),F,"cg")
    
    # save solution for first two levels in xdmf file 
    if ( j < 2 ):
        fileName = "solution_level_" + str(j) + ".xdmf"
        file = XDMFFile(fileName)
        file.write(u,0)
    print("-------------------------")
    print("generating level ", j)
    print("meshsize: ", mesh.hmin())
    print("mesh quality: ", MeshQuality.radius_ratio_min_max(mesh))  
    

    if (problemName=="2Dpeak" or problemName=="3Dpeak"):
        errorL2Norm = sqrt(assemble((u-exactSolution)**2*dx(mesh)))
        errorEnergyNorm = sqrt(assemble(inner(grad(u)-gradExactSolution,grad(u)-gradExactSolution)*dx(mesh)))
        print("error L2norm:  {0:16.8e}".format(errorL2Norm))
        print("error energy norm: {0:16.8e}".format(errorEnergyNorm))
    
    APetsc = as_backend_type(A).mat()
    indptr, indices, data = APetsc.getValuesCSR()
    AScipy = sparse.csr_matrix((data,indices,indptr))
    AScipy.eliminate_zeros()

    if (onlyFreeNodes):
        if (problemName=="2Dpeak" or problemName=="3Dpeak"):
            boundaryNodes = np.fromiter(bc.get_boundary_values().keys(), dtype=int)

        if (problemName=="2Dshell" or problemName=="3Dshell"):
            boundaryNodes1 = np.fromiter(bc1.get_boundary_values().keys(), dtype=int)  
            boundaryNodes2 = np.fromiter(bc2.get_boundary_values().keys(), dtype=int)  
            boundaryNodes = np.append(boundaryNodes1,boundaryNodes2)
        
        [nDoFs,nDoFs] = AScipy.shape 
        freeNodes = np.setdiff1d(range(0,nDoFs),boundaryNodes)
        vectorsFreeNodes[j] = freeNodes
        print("number of all DoFs: ", VFine.dim())  
        print("number of free DoFs: ", freeNodes.shape[0])  
        AScipy = AScipy[freeNodes,:]
        AScipy = AScipy[:,freeNodes]
    
    matricesA[j] = AScipy
    
    FPetsc = as_backend_type(F).vec()    
    
    FScipy = np.array(FPetsc.getValues(range(0,FPetsc.getSize())))
    FScipy = np.atleast_2d(FScipy).transpose()
    if (onlyFreeNodes):
        FScipy = FScipy[freeNodes]
    vectorsF[j] = FScipy     

    if (j != 0):
        P = PETScDMCollection.create_transfer_matrix(VCoarse, VFine)
        PPetsc = P.mat()
        indptr, indices, data = PPetsc.getValuesCSR()
        PScipy = sparse.csr_matrix((data,indices,indptr))  
        PScipy.eliminate_zeros()
        if (onlyFreeNodes):
            PScipy = PScipy[freeNodes,:]
            PScipy = PScipy[:,vectorsFreeNodes[j-1]]
        matricesP[j-1] = PScipy
        
    if (j < nlevel-1):
        VCoarse = VFine
        if (problemName =="2Dpeak"):
            mesh = UnitSquareMesh(n*(2**(j+1)),n*(2**(j+1)))
        else:
            mesh = refine(mesh)

scipy.io.savemat("A.mat", {"A": matricesA})
scipy.io.savemat("F.mat", {"F": vectorsF})
scipy.io.savemat("P.mat", {"P": matricesP})

toc = time.time()
print("Total time", toc-tic)