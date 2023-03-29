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

def genFaceBC(fileName, fdir) :
    x, y  = sp.symbols('x y')    # coordinate axes

    v1   =  sp.sin(2*sp.pi*x) * sp.cos(2*sp.pi*y)
    v2   = -sp.cos(2*sp.pi*x) * sp.sin(2*sp.pi*y)
    v    = sp.Matrix([v1, v2])
    v_fn = sp.lambdify([x, y], v, "numpy")

    vtpReader = vtk.vtkXMLPolyDataReader()
    vtpReader.SetFileName(fileName)
    vtpReader.Update()

    pdata = vtk.vtkPolyData()
    pdata = vtpReader.GetOutput()

    pdata_npts = pdata.GetNumberOfPoints()
    pdata_pts  = np.zeros((pdata_npts,3))
    pdata_nid  = np.zeros((pdata_npts,1),int)
    pdata_nid  = pdata.GetPointData().GetArray('GlobalNodeID')
    for ipt in range(0, pdata_npts):
        pdata_pts[ipt,:] = pdata.GetPoint(ipt)

    path, fname = os.path.split(fileName)
    fhdr, ext   = os.path.splitext(fname)

    fname = "%s/csv/bc_%s_nodeid.csv" % (fdir, fhdr)
    np.savetxt(fname, pdata_nid, delimiter=',', fmt='%d')

    print(bcolors.HEADER + "="*80 + bcolors.ENDC)
    print("")
    print(bcolors.OKBLUE + "Writing BC data for face %s" % \
        (fhdr) + bcolors.ENDC)
    print("")

    vg = v_fn(pdata_pts[:,0], pdata_pts[:,1])
    for j in range(0, np.size(vg,0)) :
        fname = "%s/csv/bc_%s_v%d.csv" % (fdir, fhdr, j+1)
        np.savetxt(fname, vg[j,0,:])


if __name__ == '__main__':
    N       = 256
    fhdr    =  '../N%03d' % (N)

    fname = '../../mesh/N%03d/mesh-surfaces/bottom.vtp' %(N)
    genFaceBC(fname, fhdr)

    fname = '../../mesh/N%03d/mesh-surfaces/top.vtp' %(N)
    genFaceBC(fname, fhdr)

    fname = '../../mesh/N%03d/mesh-surfaces/left.vtp' %(N)
    genFaceBC(fname, fhdr)

    fname = '../../mesh/N%03d/mesh-surfaces/right.vtp' %(N)
    genFaceBC(fname, fhdr)
