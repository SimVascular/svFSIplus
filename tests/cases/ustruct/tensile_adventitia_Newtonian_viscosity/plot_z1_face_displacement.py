'''
Plots the displacement of the z1 face of the mesh over time
'''

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import os
import glob
import xml.etree.ElementTree as ET

def get_timestep(result_file):
    '''
    Extracts the timestep of a results_###.vtu file
    '''
    # If file path is provided, get the file name
    file_name = os.path.basename(result_file)

    # Get the ###.vtu part
    s = file_name.split('_')[-1]

    # Get the ### part
    s = s.split('.')[0]

    # Return as integer
    return int(s)

def get_start_end_step(results_folder):
    """
    Automatically determine the start timestep, end timestep, and step size based on all
    svFSI results file in results_folder.

    Args:
        results_folder: A string of absolute file path of folder containing results of 
        svFSIplus simulation. This usually ends with 16-procs/ or other 
        number.
    
    Returns:
        (start_timestep, end_timestep, step): A tuple of 3 integers, giving the first 
        timestep of results to process, the last timestep, and the step size. 

    """

    # Get list of all .vtu files in results_folder sorted by time step
    list_of_files = sorted( filter( os.path.isfile,
                            glob.glob(os.path.join(results_folder, '*.vtu')) ), key = get_timestep)

    # Get start time from the first result file (list_of_files[0])
    start_file_name = os.path.basename(list_of_files[0])
    start_timestep = int("".join([i for i in start_file_name if i.isdigit()]))
    print('Start timestep:', start_timestep)

    # Get end time from the last result file (list_of_files[-1])
    end_file_name = os.path.basename(list_of_files[-1])
    end_timestep = int("".join([i for i in end_file_name if i.isdigit()]))
    print('End timestep:', end_timestep)

    # Get step size by looking at second time step
    start_plus_one_file_name = os.path.basename(list_of_files[1])
    start_time_plus_one = int("".join([i for i in start_plus_one_file_name if i.isdigit()]))
    step = start_time_plus_one - start_timestep
    print('Step:', step)

    return (start_timestep, end_timestep, step)

def get_timestep_size(svfsiplus_xml_file):
    '''
    Reads the timestep size from an svFSIplus XML file.

    Args:
        svfsiplus_xml_file: The file path to the svFSIplus XML input file

    Returns:
        timestep_size: The size of the timestep in the input file
    '''
    # Parse the XML file
    tree = ET.parse(svfsiplus_xml_file)
    root = tree.getroot()

    # Find the element containing the Time_step_size
    timestep_size = None
    for elem in root.iter():
        if elem.tag == 'Time_step_size':
            timestep_size = float(elem.text)
            break

    return timestep_size


# Move to this script's directory
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Load the z1 face
mesh = pv.read("mesh/mesh-surfaces/Z1.vtp")

# Get the start timestep, end timestep, and step size
start_timestep, end_timestep, step = get_start_end_step("4-procs/")

# Get the timestep size
timestep_size = get_timestep_size("svFSIPlus.xml")

# Get the time values
time = np.arange(start_timestep, end_timestep + 1, step) * timestep_size

# Loop over the results files
displacements = []
for n in range(start_timestep, end_timestep + 1, step):
    # Load the results file
    result_file = f"4-procs/result_{n:03d}.vtu"
    result = pv.read(result_file)

    # Sample displacement at the z1 face
    displacement = mesh.sample(result).point_data['Displacement'][:, 2].mean()

    # Append to the list
    displacements.append(displacement)

# Save the displacements to a text file
with open("z1_face_displacement.dat", "w") as file:
    file.write(f"{len(time)}\n")
    for t, s in zip(time, displacements):
        file.write(f"{t:.3f} {s:.3f}\n")


# Plot the displacement
plt.plot(time, displacements, label="With viscosity")
plt.xlabel("Time (s)")
plt.ylabel("z Displacement (cm)")
plt.title("Z1 face displacement")

# Plot the displacement with no viscosity
if os.path.exists("z1_face_displacement_no_viscosity.dat"):
    with open("z1_face_displacement_no_viscosity.dat", "r") as file:
        lines = file.readlines()
        n = int(lines[0])
        time_no_viscosity = np.zeros(n)
        displacements_no_viscosity = np.zeros(n)
        for i, line in enumerate(lines[1:]):
            time_no_viscosity[i], displacements_no_viscosity[i] = map(float, line.split())

    plt.plot(time_no_viscosity, displacements_no_viscosity, label="Without viscosity")


plt.legend()
plt.savefig("z1_face_displacement.png")