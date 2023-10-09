# Header
import numpy as np
import vtk

#----------------------------------------------------------------------
def getWallNodes(FileName):
    modelReader = vtk.vtkXMLPolyDataReader()
    modelReader.SetFileName(FileName)
    modelReader.Update()

    model = vtk.vtkPolyData()
    model = modelReader.GetOutput()

    model_npts = model.GetNumberOfPoints()
    model_ID = model.GetPointData().GetArray('GlobalNodeID')
    id_list = np.zeros((model_npts,1), dtype='int32')
    for i in range(0, model_npts):
        id_list[i] = int(model_ID.GetTuple1(i))
    return id_list
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def loadVTU(fileName):

    print ("   Loading vtu file   <---   {}".format(fileName))
    vtuReader = vtk.vtkXMLUnstructuredGridReader()
    vtuReader.SetFileName(fileName)
    vtuReader.Update()

    vtuMesh = vtk.vtkUnstructuredGrid()
    vtuMesh = vtuReader.GetOutput()

    return vtuMesh
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def getSurfaceData(msh, srf_ids, keyWord):
    msh_npts = msh.GetNumberOfPoints()
    msh_data = vtk.vtkDoubleArray()
    msh_data = msh.GetPointData().GetArray(keyWord)
    num_comp = msh_data.GetNumberOfComponents()

    srf_nno = np.size(srf_ids)
    if num_comp == 1:
        srf_data = np.zeros((srf_nno,))
    else:
        srf_data = np.zeros((srf_nno,num_comp))
    for ipt in range(0, srf_nno):
        h = msh_data.GetTuple(int(srf_ids[ipt])-1)
        if num_comp == 1:
            srf_data[ipt] = h[0]
        else:
            for j in range(0, num_comp):
                srf_data[ipt,j] = h[j]

    return srf_data
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def writeSrfTraction(h, fout, fwall):
    modelReader = vtk.vtkXMLPolyDataReader()
    modelReader.SetFileName(fwall)
    modelReader.Update()

    model = vtk.vtkPolyData()
    model = modelReader.GetOutput()
    model_npts = modelReader.GetNumberOfPoints()

    vtkH = vtk.vtkDoubleArray()
    vtkH.SetNumberOfComponents(3)
    vtkH.Allocate(model_npts)
    vtkH.SetNumberOfTuples(model_npts)
    vtkH.SetName("Traction")
    for i in range(0, model_npts):
        vtkH.SetTuple3(i, h[i,0], h[i,1], h[i,2])

    model.GetPointData().AddArray(vtkH)
    modelWrite = vtk.vtkXMLPolyDataWriter()
    modelWrite.SetInputData(model)
    modelWrite.SetFileName(fout)
    modelWrite.Write()

    return
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def writeSrfPressure(p, fout, fwall):
    modelReader = vtk.vtkXMLPolyDataReader()
    modelReader.SetFileName(fwall)
    modelReader.Update()

    model = vtk.vtkPolyData()
    model = modelReader.GetOutput()
    model_npts = modelReader.GetNumberOfPoints()

    vtkP = vtk.vtkDoubleArray()
    vtkP.SetNumberOfComponents(1)
    vtkP.Allocate(model_npts)
    vtkP.SetNumberOfTuples(model_npts)
    vtkP.SetName("Pressure")
    for i in range(0, model_npts):
        vtkP.SetTuple1(i, p[i])

    model.GetPointData().AddArray(vtkP)
    modelWrite = vtk.vtkXMLPolyDataWriter()
    modelWrite.SetInputData(model)
    modelWrite.SetFileName(fout)
    modelWrite.Write()

    return
#----------------------------------------------------------------------

if __name__ == '__main__':
    srcdir = "24-procs"
    nstart = 600
    nend   = 800
    nfreq  = 100
    dt     = 0.005
    fwall  = "../mesh/mesh-surfaces/lumen_wall.vtp"

    wall_ids = getWallNodes(fwall)
    wall_nno = np.size(wall_ids)
    mean_h = np.zeros((wall_nno,3), dtype='float64')
    mean_P = np.zeros((wall_nno,), dtype='float64')

    nframe = int((nend - nstart)/nfreq) + 1
    for i in range(0, nframe):
        ntime = nstart + i*nfreq
        time  = float(ntime)*dt
        if (ntime < 100):
            fname = "%s/result_%03d_cpp.vtu" %(srcdir, ntime)
        else:
            fname = "%s/result_%d_cpp.vtu" %(srcdir, ntime)
        print ("Reading file    <-----   {}".format(fname))
        vtuMesh = loadVTU(fname)
        mean_h = mean_h + getSurfaceData(vtuMesh, wall_ids, 'Traction')
        mean_P = mean_P + getSurfaceData(vtuMesh, wall_ids, 'Pressure')

    mean_h = mean_h / float(nframe)
    mean_P = mean_P / float(nframe)
    fout = 'rigid_wall_mean_traction.vtp'
    print ("Writing traction file    ---->   {}".format(fout))
    writeSrfTraction(mean_h, fout, fwall)
    fout = 'rigid_wall_mean_pressure.vtp'
    print ("Writing pressure file    ---->   {}".format(fout))
    writeSrfPressure(mean_P, fout, fwall)

#=======================================================================
# EOF


