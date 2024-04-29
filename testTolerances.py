import os
import vtk
import time

def compute_difference(file1, file2):
    # Read current folder .vtu
    reader1 = vtk.vtkXMLUnstructuredGridReader()
    reader1.SetFileName(file1)
    reader1.Update()
    grid1 = reader1.GetOutput()
    pressure_data1 = grid1.GetPointData().GetArray("Pressure")

    # Read reference folder .vtu
    reader2 = vtk.vtkXMLUnstructuredGridReader()
    reader2.SetFileName(file2)
    reader2.Update()
    grid2 = reader2.GetOutput()
    pressure_data2 = grid2.GetPointData().GetArray("Pressure")

    # The grids should be identical
    if grid1.GetNumberOfPoints() != grid2.GetNumberOfPoints():
        raise ValueError("Meshes have different number of points")

    # Compute the difference between the values (reference-current)
    num_points = grid1.GetNumberOfPoints()
    diff = vtk.vtkDoubleArray()
    diff.SetName("Pdiff")
    diff.SetNumberOfComponents(1)
    diff.SetNumberOfTuples(num_points)
    for i in range(num_points):
        diff.SetValue(i, pressure_data2.GetValue(i) - pressure_data1.GetValue(i))

    # Add the difference array to the first grid
    grid1.GetPointData().AddArray(diff)

    return grid1

def list_folders(directory):
    folders = [f for f in os.listdir(directory) if os.path.isdir(os.path.join(directory, f))]
    return folders

def check_folder_p(folder_name):
    return folder_name.startswith("petsc_")
def check_folder_tr(folder_name):
    return folder_name.startswith("trilinos_")
def check_folder_f(folder_name):
    return folder_name.startswith("fsils_")

def main():

    # Open log file
    file_path = "/solver/logTolerances.txt"
    with open(file_path, "w") as file:
        file.write("====================================== \n")
        file.write(" Fluid/pipe_RCR_3d test on tolerances  \n")
        file.write("====================================== \n")
        file.write("\n")

    test_cases_path = "/solver/svFSIplus/tests/cases/fluid/test_tol_ls"
    file_name = ["result_050.vtu","result_100.vtu"]
    bin_petsc_folder = "/solver/build-ubuntu20-petsc/svFSI-build/bin/svFSI"
    bin_trilinos_folder = "/solver/build-ubuntu20-trilinos/svFSI-build/bin/svFSI"

    # List all the folders in test_tol_ls 
    folders = list_folders(test_cases_path)
    print(folders)
    
    # Run Simulation
    for folder in folders:
        
        os.chdir(os.path.join(test_cases_path,folder))

        if check_folder_p(folder):
            start_time = time.time()
            os.system("mpirun -n 4 " + bin_petsc_folder + " svFSI.xml")
            end_time = time.time()
        if check_folder_tr(folder):
            start_time = time.time()
            os.system("mpirun -n 4 " + bin_trilinos_folder + " svFSI.xml")
            end_time = time.time()
        if check_folder_f(folder):
            start_time = time.time()
            os.system("mpirun -n 4 " + bin_trilinos_folder + " svFSI.xml")
            end_time = time.time()
            # Build reference folders fsils_* 12, 13, 14
            if folder.endswith("12"):
                ref_path_12=os.path.join(test_cases_path,folder)
            if folder.endswith("13"):
                ref_path_13=os.path.join(test_cases_path,folder)
            if folder.endswith("14"):
                ref_path_14=os.path.join(test_cases_path,folder)

        elapsed_time = end_time - start_time
        with open(file_path, "a") as file:
            file.write("Time for "+ folder + ": " + str(elapsed_time) +" seconds \n")
            file.write("\n")

    # PostProcess
    for folder in folders:
        # fsils_ folders are the references
        if folder.startswith("fsils_"): 
            continue

        reference_folder = os.path.join(ref_path_12,"4-procs") if folder.endswith("12") else None
        reference_folder = os.path.join(ref_path_13,"4-procs") if folder.endswith("13") else None
        reference_folder = os.path.join(ref_path_14,"4-procs") if folder.endswith("14") else None

        current_folder = os.path.join(test_cases_path,folder,"4-procs")
        output_folder = current_folder
        # Iterate over pairs of files and compute the difference
        for f in file_name:
            result_mesh = compute_difference(os.path.join(current_folder,f), os.path.join(reference_folder,f))
            # Save the result to a new .vtu file
            output_file = os.path.join(output_folder, "p.diff."+f)

            # Save the modified grid to a new .vtu file
            writer = vtk.vtkXMLUnstructuredGridWriter()
            writer.SetFileName(output_file)
            writer.SetInputData(result_mesh)
            writer.Write()


if __name__ == "__main__":
    main()
