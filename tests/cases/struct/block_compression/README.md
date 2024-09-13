
# **Problem Description**

Simulate a block compression problem using a displacement-based solid equation. 

## Directional Dirichlet BC

The problem configuration requires enforcing Dirichlet boundary condition along a specific direction. 

For example the displacement of the **X0** face is constrained to be zero along the Z-axis.
```
<Add_BC name="X0" > 
  <Type> Dir </Type> 
  <Value> 0.0 </Value> 
  <Effective_direction> (0, 0, 1) </Effective_direction> 
</Add_BC>
```


## References
Liu, Ju, and Alison L. Marsden.  A Unified Continuum and Variational Multiscale Formulation for Fluids, Solids, and Fluid Structure Interaction.  *Computer Methods in Applied Mechanics and Engineering* 337 (August 2018): 549 97. https://doi.org/10.1016/j.cma.2018.03.045.
