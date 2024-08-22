
# **Problem Description**

Simulation of a 3D pipe flow problem using a Carreau-Yasuda viscosity model

Carreau-Yasuda viscosity model is defined by the following equation

![cy_model](https://latex.codecogs.com/png.image?\inline&space;\dpi{120}\bg{white}\eta=\eta_\infty&plus;(\eta_0-\eta_\infty)\left[&space;1&plus;\left(&space;\lambda&space;\dot{\gamma}&space;\right)^a&space;\right]^{\frac{n-1}{a}}&space;)

where

![iamge1](https://latex.codecogs.com/png.image?\inline&space;\dpi{120}\bg{white}\eta_\infty:&space;\text{Limiting&space;high&space;shear-rate&space;viscosity})

![iamge2](https://latex.codecogs.com/png.image?\inline&space;\dpi{120}\bg{white}\eta_0:&space;\text{Limiting&space;low&space;shear-rate&space;viscosity}&space;)

![iamge3](https://latex.codecogs.com/png.image?\inline&space;\dpi{120}\bg{white}\lambda:&space;\text{Shear-rate&space;tensor&space;multiplier}&space;)

![iamge4](https://latex.codecogs.com/png.image?\inline&space;\dpi{120}\bg{white}\dot{\gamma}:&space;\text{Shear&space;rate}&space;)

![iamge5](https://latex.codecogs.com/png.image?\inline&space;\dpi{120}\bg{white}a:&space;\text{Shear-rate&space;tensor&space;exponent}&space;)

![iamge6](https://latex.codecogs.com/png.image?\inline&space;\dpi{120}\bg{white}n:&space;&space;\text{Power-law&space;index})

Casson model is defined as

![ca_model](https://latex.codecogs.com/png.image?\inline&space;\dpi{120}\bg{white}\eta&space;=&space;\frac{1}{\dot{\gamma}}\left[&space;k_0(c)&space;&plus;&space;k_1(c)\sqrt{\dot{\gamma}}&space;\right]^2)

$k_0(c)$ and $k_1(c)$ are functions of the hematocrit $c$. For more information, please refer to Ref. [1].


The model parameters are specified in the `Viscosity` sub-section

```
<Viscosity model="Carreau-Yasuda" >
  <Limiting_high_shear_rate_viscosity> 0.022 </Limiting_high_shear_rate_viscosity>
  <Limiting_low_shear_rate_viscosity> 0.22 </Limiting_low_shear_rate_viscosity> 
  <Shear_rate_tensor_multiplier> 0.11 </Shear_rate_tensor_multiplier> 
  <Shear_rate_tensor_exponent> 0.644 </Shear_rate_tensor_exponent> 
  <Power_law_index> 0.392 </Power_law_index> 
</Viscosity>
```



## Reference
Boyd, Joshua, James M. Buick, and Simon Green.  Analysis of the Casson and Carreau-Yasuda Non-Newtonian Blood Models in Steady and Oscillatory Flows Using the Lattice Boltzmann Method.  *Physics of Fluids* 19, no. 9 (September 2007): 093103. https://doi.org/10.1063/1.2772250.
