This test case simulates prescribed extension and compression of a slab of material
described by the Neohookean material model and a Newtonian solid viscosity model.

The load profile is a 0.5s ramp to 1e5 dynes/cm^2, then held for another 0.5s. 
The load is applied on the Z1 face in the z-direction.

![Load Profile](load.png)

The load data is defined in `load.dat`, which can be generated with
`load.py`

The resulting deformation is shown in the video below:
![Deformation](movie.gif)

The z-displacement of the Z1 face is plotted below:
![Z1 Displacement](z1_face_displacement.png)

