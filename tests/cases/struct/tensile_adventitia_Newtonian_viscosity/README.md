This test case simulates pulling and releasing a slab of material
described by the Guccione material model and a Newtonian solid viscosity model. 
In this viscosity model, the viscous deviatoric Cauchy stress is identical to that
for a Newtonian fluid:

$$
\sigma^{dev}_{vis} = 2 \mu \mathbf{d}^{dev}
$$
where
$$
\mathbf{d}^{dev} = \frac{1}{2} (\nabla_x \mathbf{v} + (\nabla_x \mathbf{v})^T) - \frac{1}{3} (\nabla_x \cdot \mathbf{v}) \mathbf{I}
$$

The viscous part of the 2nd Piola-Kirchhoff stress is then given by a pull-back operation
$$
\mathbf{S}_{vis} = 2 \mu J \mathbf{F}^{-1} \mathbf{d}^{dev} \mathbf{F}^{-T}
$$

The load profile is a 0.25s ramp to the max value, then held for 0.25s, then a 0.25s ramp down to 0, then held for 0.25s.
The load is applied on the Z1 face in the z-direction.

![Load Profile](load.png)

The load data is defined in `load.dat`, which can be generated with
`load.py`

The resulting deformation is shown in the video below:
![Deformation](movie.gif)

The z-displacement of the Z1 face is plotted below:
![Z1 Displacement](z1_face_displacement.png)

