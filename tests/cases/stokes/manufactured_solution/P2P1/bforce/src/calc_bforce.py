# Header
import sympy as sp
import numpy as np
import vtk
import os
sp.init_printing(use_unicode = True)

class bcolors:
    HEADER  = '\033[95m'
    OKBLUE  = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL    = '\033[91m'
    ENDC    = '\033[0m'
    BOLD    = '\033[1m'
    ULINE   = '\033[4m'

# Define symbols
x, y  = sp.symbols('x y')    # coordinate axes

# Input parameters
N       = 256
mshFile = '../../mesh/N%03d/mesh-complete.mesh.vtu' % (N)
fhdr    =  '../N%03d' % (N)
os.system("mkdir -p %s/csv" % (fhdr))

#=======================================================================
# Velocity
v1   =  sp.sin(2*sp.pi*x) * sp.cos(2*sp.pi*y)
v2   = -sp.cos(2*sp.pi*x) * sp.sin(2*sp.pi*y)
v    = sp.Matrix([v1, v2])

print(bcolors.HEADER + "="*80 + bcolors.ENDC)
print("")
print(bcolors.OKBLUE + "Velocity field, v: " + bcolors.ENDC)
print("")
sp.pprint(v)
print("")

#=======================================================================
# Pressure:
p    = -4*sp.pi*sp.cos(2*sp.pi*x) * sp.cos(2*sp.pi*y)
print(bcolors.HEADER + "="*80 + bcolors.ENDC)
print("")
print(bcolors.OKBLUE + "Pressure, p: " + bcolors.ENDC)
print("")
sp.pprint(p)
print("")

#=======================================================================
# Velocity gradient
gradV = sp.Matrix([[sp.diff(v1,x), sp.diff(v1,y)],\
                   [sp.diff(v2,x), sp.diff(v2,y)]])

#=======================================================================
# Stress
sigma = -p*sp.eye(2) + gradV + gradV.T
print(bcolors.HEADER + "="*80 + bcolors.ENDC)
print("")
print(bcolors.OKBLUE + "Cauchy stress: " + bcolors.ENDC)
print("")
sp.pprint(sigma)
print("")

#=======================================================================
# Body force
fb = -(sp.diff(sigma.col(0), x) + sp.diff(sigma.col(1), y))
print(bcolors.HEADER + "="*80 + bcolors.ENDC)
print("")
print(bcolors.OKBLUE + "Body force, fb: " + bcolors.ENDC)
print("")
sp.pprint(fb)
print("")

fb_fn = sp.lambdify([x, y], fb, "numpy")

#=======================================================================
# Read VTU mesh file and load point coordinates
mshReader = vtk.vtkXMLUnstructuredGridReader()
mshReader.SetFileName(mshFile)
mshReader.Update()

msh = vtk.vtkUnstructuredGrid()
msh = mshReader.GetOutput()

msh_npts = msh.GetNumberOfPoints()
msh_pts  = np.zeros((msh_npts,3))
for ipt in range(0, msh_npts):
    msh_pts[ipt,:] = msh.GetPoint(ipt)

#=======================================================================
# Evaluate body force over the domain and write to file
bff = vtk.vtkDoubleArray()
bff.SetNumberOfComponents(3)
bff.Allocate(msh_npts)
bff.SetNumberOfTuples(msh_npts)
bff.SetName("FB")

print(bcolors.HEADER + "="*80 + bcolors.ENDC)
print("")
print(bcolors.OKBLUE + "Writing csv and vtu files" + bcolors.ENDC)
print("")

fbg = fb_fn(msh_pts[:,0], msh_pts[:,1])

fname = "%s/csv/bforce.csv" % (fhdr)
np.savetxt(fname, fbg[0][0], delimiter=',')

for j in range(0, msh_npts):
    rtmp = fbg[0][0][j]
    bff.SetTuple3(j, rtmp, 0, 0)
fname = "%s/bforce.vtu" % (fhdr)
msh.GetPointData().AddArray(bff)
mshWrite = vtk.vtkXMLUnstructuredGridWriter()
mshWrite.SetInputData(msh)
mshWrite.SetFileName(fname)
mshWrite.Write()

print(bcolors.HEADER + "="*80 + bcolors.ENDC)


#=======================================================================
# EOF

