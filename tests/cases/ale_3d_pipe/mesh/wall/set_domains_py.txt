import os
import sys
import vtk
import time
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

# Global parameters
EPS = sys.float_info.epsilon
PI  = np.pi

# User defined inputs
MESH_FILE = "mesh-complete.mesh.vtu"
N_DOMAIN  = 5

#----------------------------------------------------------------------
def is_zero(a, b=None):

    absA = abs(a)
    if b == None:
        b = 0.0

    absA = abs(a)
    absB = abs(b)
    if absB > absA:
        absA, absB = absB, absA

    nrm  = max(absA, absB)
    flag = False
    if ((absA-absB)/nrm) < (10.0*EPS):
        flag = True

    return flag
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def norm(u, v=None):

    if v == None:
        v = u

    n = np.size(u)
    l2norm = 0.0
    for i in range(0, n):
        l2norm = l2norm + (u[i]*v[i])

    return l2norm
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def create_VTK_data_array(dType, numComp, numTupls, dName):

    if dType == "double":
        D = vtk.vtkDoubleArray()
    elif dType == "int":
        D = vtk.vtkIntArray()

    D.SetNumberOfComponents(numComp)
    D.Allocate(numTupls)
    D.SetNumberOfTuples(numTupls)
    D.SetName(dName)

    return D
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def load_mesh(file_name):

    print "   Loading mesh file   <---   %s" % (file_name)
    vtuReader = vtk.vtkXMLUnstructuredGridReader()
    vtuReader.SetFileName(file_name)
    vtuReader.Update()

    vtuMesh = vtk.vtkUnstructuredGrid()
    vtuMesh = vtuReader.GetOutput()

    return vtuMesh
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def get_z_bounds(msh):
    msh_npts = msh.GetNumberOfPoints()
    zmin = 9999999.9
    zmax = -zmin
    for ipt in range(0, msh_npts):
        x = msh.GetPoint(ipt)
        if zmin > x[2]:
            zmin = x[2]
        if zmax < x[2]:
            zmax = x[2]

    return zmin, zmax
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def set_domain_ID(msh, zmin, zmax):

    print "   Setting domain IDs at cells"
    num_cell = msh.GetNumberOfCells()
    dmn_IDs  = create_VTK_data_array("int", 1, num_cell, "DOMAIN_ID")

    z_part = np.linspace(zmin, zmax, N_DOMAIN+1)

    for icell in xrange(0, num_cell):
        curr_cell = msh.GetCell(icell)
        pts_cell = curr_cell.GetPointIds()

        p0 = msh.GetPoint(pts_cell.GetId(0))
        p1 = msh.GetPoint(pts_cell.GetId(1))
        p2 = msh.GetPoint(pts_cell.GetId(2))
        p3 = msh.GetPoint(pts_cell.GetId(3))

        ctr_z = 0.25*(p0[2] + p1[2] + p2[2] + p3[2])

        for i in range(0, N_DOMAIN):
            if (ctr_z>=z_part[i]-EPS) and (ctr_z<z_part[i+1]+EPS):
                dmn_id = i
        dmn_IDs.SetTuple1(icell, dmn_id)

    return dmn_IDs
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Main function
if __name__ == '__main__':

    t1 = time.time()
    print "========================================================"
    msh = load_mesh(MESH_FILE)

    zmin, zmax = get_z_bounds(msh)

    dmnIDs = set_domain_ID(msh, zmin, zmax)

    print "   Writing domains and fibers to VTK data structure"
    msh.GetCellData().AddArray(dmnIDs)

    fileName = "domains.vtu"
    vtuWriter = vtk.vtkXMLUnstructuredGridWriter()
    vtuWriter.SetInputData(msh)
    vtuWriter.SetFileName(fileName)
    print "   Writing to vtu file   --->   %s" % (fileName)
    vtuWriter.Write()
    print "========================================================"

#----------------------------------------------------------------------


