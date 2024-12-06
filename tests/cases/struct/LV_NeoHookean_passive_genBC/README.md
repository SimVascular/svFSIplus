This test case simulates an idealized left ventricle (LV) with a NeoHookean material model
coupled to a lumped-parameter network (LPN), implemented in genBC. The LPN consists of a large pressure source and large resistor, which together produce an approximately constant flowrate into
the LV. This inflates the LV at an approximately constant rate of change of volume.

Before running svMultiPhysics , genBC_svMultiPhysics must be compiled. Navigate to genBC_svMultiPhysics
and run `make clean` then `make`.

The results can be post-processed by running `process_results.py`. If run for 100
timesteps, the results should be the following.

![V3D vs. V0D](V3D_vs_V0D.png)

*Plot of volume vs. time for this simulation. V_3D is computed from results.vtu
by extracting the endocardial surface, warping it according to displacement at
each timestep, capping it, then measuring the enclosed volume. V_0D is read 
from AllData (produced by genBC).*


![dVdt3D vs. dVdt0D](dVdt3D_vs_dVdt0D.png)

*Plot of flowrate (dVdt) vs. time for this simulation. dVdt_3D is computed from results.vtu
by extracting the endocardial surface, warping it according to displacement at
each timestep, then computing the velocity flux integral.
V_0D is computed from AllData (produced by genBC). As seen in the plot, we
ramp the flowrate from  0 to 100 cm^3/s over 10 timesteps.*

![Pressure vs. Volume](pv_plot.png)

*Plot of pressure vs. volume as the LV is inflated.*
