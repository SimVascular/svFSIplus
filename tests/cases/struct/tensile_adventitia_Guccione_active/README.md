This test case simulates prescribed active contraction of a slab of material
described by the Guccione material model. Primary fibers run along the length
of the slab (z-direction) and secondary fibers run across the width of the slab
(x-direction).

The active stress profile is sinusoidal over a period of 1 second with an amplitude
of 100 dynes/cm^2. 

![Stress Profile](stress.png)

The active stress data is defined in `stress.dat`, which can be generated with
`generate_stress.py`

The resulting deformation is shown in the video below:
![Deformation](movie.gif)

