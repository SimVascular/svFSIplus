This test case simulates pulling and releasing a slab of material
described by the Guccione material model and a pseudo-potential solid viscosity model.
In this viscosity model, the so-called viscous pseudo-potential $\Psi_{vis}$ is
given by
$$
\Psi_{vis} =  \frac{\mu}{2} \text{tr}(\dot{\mathbf{E}})
$$
The viscous part of the 2nd Piola-Kirchhoff stress is then given by
$$
S_{vis} = \frac{\partial \Psi_{vis}}{\partial \dot{\mathbf{E}}} = \mu  \dot{\mathbf{E}}
$$



The load profile is a 0.5s ramp to 1e5 dynes/cm^2, then held for another 0.5s. 
The load is applied on the Z1 face in the z-direction.

![Load Profile](load.png)

The load data is defined in `load.dat`, which can be generated with
`load.py`

The resulting deformation is shown in the video below:
![Deformation](movie.gif)

The z-displacement of the Z1 face is plotted below:
![Z1 Displacement](z1_face_displacement.png)

