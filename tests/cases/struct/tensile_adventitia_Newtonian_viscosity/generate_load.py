import numpy as np
import os

# Go to directory of this script
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Number of timesteps and number of Fourier modes
n_timesteps = 1001
n_modes = 256

# Generate time values from 0 to 1
time = np.linspace(0, 1, n_timesteps)

# Generate 0.25s ramp to max load, then hold for 0.25s, then 0.25s ramp down to 0, then hold for 0.25s
load = np.zeros(n_timesteps)
max_load= -1e2
load[time < 0.25] = max_load * time[time < 0.25] / 0.25
load[(time >= 0.25) & (time < 0.5)] = max_load
load[(time >= 0.5) & (time < 0.75)] = max_load * (1 - (time[(time >= 0.5) & (time < 0.75)] - 0.5) / 0.25)
load[time >= 0.75] = 0

# Write the time and load values to a text file
with open("load.dat", "w") as file:
    file.write(f"{n_timesteps} {n_modes}\n")
    for t, s in zip(time, load):
        file.write(f"{t:.3f} {s:.3f}\n")

# Plot the load values
import matplotlib.pyplot as plt
plt.plot(time, load)
plt.xlabel("Time (s)")
plt.ylabel("load (dynes/cm^2)")
plt.title("Active load")
plt.savefig("load.png")