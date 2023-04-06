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
os.system("mkdir -p %s" % (fhdr))

#=======================================================================
# Velocity
v1   =  sp.sin(2*sp.pi*x) * sp.cos(2*sp.pi*y)
v2   = -sp.cos(2*sp.pi*x) * sp.sin(2*sp.pi*y)
v    = sp.Matrix([v1, v2])
v_fn = sp.lambdify([x, y], v, "numpy")

print(bcolors.HEADER + "="*80 + bcolors.ENDC)
print("")
print(bcolors.OKBLUE + "Velocity field, v: " + bcolors.ENDC)
print("")
sp.pprint(v)
print("")

#=======================================================================
# Pressure:
p    = -4*sp.pi*sp.cos(2*sp.pi*x) * sp.cos(2*sp.pi*y)
p_fn = sp.lambdify([x, y], p, "numpy")

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
s_fn  = sp.lambdify([x, y], sigma, "numpy")

print(bcolors.HEADER + "="*80 + bcolors.ENDC)
print("")
print(bcolors.OKBLUE + "Cauchy stress: " + bcolors.ENDC)
print("")
sp.pprint(sigma)
print("")

#=======================================================================
# Read vtk mesh
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
# Evaluate data at the mesh point coordiantes
vg = v_fn(msh_pts[:,0], msh_pts[:,1])
pg = p_fn(msh_pts[:,0], msh_pts[:,1])
sg = s_fn(msh_pts[:,0], msh_pts[:,1])

#=======================================================================
# Copy data to vtk objects and write to file
vtkV = vtk.vtkDoubleArray()
vtkV.SetNumberOfComponents(3)
vtkV.Allocate(msh_npts)
vtkV.SetNumberOfTuples(msh_npts)
vtkV.SetName("Velocity")

vtkP = vtk.vtkDoubleArray()
vtkP.SetNumberOfComponents(1)
vtkP.Allocate(msh_npts)
vtkP.SetNumberOfTuples(msh_npts)
vtkP.SetName("Pressure")

vtkS = vtk.vtkDoubleArray()
vtkS.SetNumberOfComponents(3)
vtkS.Allocate(msh_npts)
vtkS.SetNumberOfTuples(msh_npts)
vtkS.SetName("Stress")

print(bcolors.HEADER + "="*80 + bcolors.ENDC)
print("")
print(bcolors.OKBLUE + "Writing vtu files" + bcolors.ENDC)
print("")

sg11 = sg[0][0]
for j in range(0, msh_npts):
    vtkP.SetTuple1(j, pg[j])
    vtkV.SetTuple3(j, vg[0,0,j], vg[1,0,j], 0.0)
    vtkS.SetTuple3(j, sg11[j], 0.0, 0.0)

fname = "%s/exact_soln.vtu" % (fhdr)
msh.GetPointData().AddArray(vtkP)
msh.GetPointData().AddArray(vtkV)
msh.GetPointData().AddArray(vtkS)
mshWrite = vtk.vtkXMLUnstructuredGridWriter()
mshWrite.SetInputData(msh)
mshWrite.SetFileName(fname)
mshWrite.Write()

print(bcolors.HEADER + "="*80 + bcolors.ENDC)

#=======================================================================
# EOF

