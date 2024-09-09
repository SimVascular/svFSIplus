This test case simulates prescribed extension and compression of a slab of material
described by the Neohookean material model and a Newtonian solid viscosity model.

The load profile is sinusoidal over a period of 1 second with an amplitude
of 100 dynes/cm^2. The load is applied on the Z1 face in the z-direction.

![Load Profile](load.png)

The load data is defined in `load.dat`, which can be generated with
`load.py`

The resulting deformation is shown in the video below:
![Deformation](movie.gif)

