
# **Problem Description**

Simulation of a 3D pipe flow problem using the Casson viscosity model.

Casson model is defined by the following equation

![ca_model](https://latex.codecogs.com/png.image?\inline&space;\dpi{120}\bg{white}\eta&space;=&space;\frac{1}{\dot{\gamma}}\left[&space;k_0(c)&space;&plus;&space;k_1(c)\sqrt{\dot{\gamma}}&space;\right]^2)

where $k_0(c)$ and $k_1(c)$ are functions of the hematocrit $c$.

The model parameters are specified in the `Viscosity` sub-section
```
<Viscosity model="Cassons" >
  <Asymptotic_viscosity_parameter> 0.3953 </Asymptotic_viscosity_parameter> 
  <Yield_stress_parameter> 0.22803 </Yield_stress_parameter> 
  <Low_shear_rate_threshold> 0.5 </Low_shear_rate_threshold> 
</Viscosity>
```


## Reference
Boyd, Joshua, James M. Buick, and Simon Green.  Analysis of the Casson and Carreau-Yasuda Non-Newtonian Blood Models in Steady and Oscillatory Flows Using the Lattice Boltzmann Method.  *Physics of Fluids* 19, no. 9 (September 2007): 093103. https://doi.org/10.1063/1.2772250.
