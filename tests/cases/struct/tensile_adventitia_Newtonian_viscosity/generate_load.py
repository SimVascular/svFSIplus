import numpy as np
import os

# Go to directory of this script
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Number of timesteps and number of Fourier modes
n_timesteps = 1001
n_modes = 256

# Generate time values from 0 to 1
time = np.linspace(0, 1, n_timesteps)

# Generate sinusoidal load values with an amplitude of 100
amplitude = 100
load = amplitude * np.sin(2 * np.pi * time)

# Write the time and load values to a text file
with open("load.dat", "w") as file:
    file.write(f"{n_timesteps} {n_modes}\n")
    for t, s in zip(time, load):
        file.write(f"{t:.3f} {s:.3f}\n")

# # Plot the load values
# import matplotlib.pyplot as plt
# plt.plot(time, load)
# plt.xlabel("Time (s)")
# plt.ylabel("load (dynes/cm^2)")
# plt.title("Active load")
# plt.savefig("load.png")