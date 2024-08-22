
# **Problem Description**

Simulate a fluid flow with dye transport in a cylindrical tube with RCR boundary conditions at the outlet and unsteady flow at the inlet. The dye is passively transported by the flow through advection and diffusion.

## Scalar transport equation

In addition to the fluid equation, the advection-diffusion equation that governs the dye transportation is added to the input file. 

The `Coupled` keyword set to `false` under the `<Add_equation type="scalarTransport" >` specifies a one-way coupling and the dye is passively transported by the flow.
```
<Coupled> false </Coupled>
```

The **Temperature** quantities normally output for the simulation are renamed `Concentration` 

```
<Output type="Alias" >
  <Temperature> Concentration </Temperature>
</Output>
```

The inlfow boundary condition specifies a contant inflow of dye
```
<Add_BC name="lumen_inlet" >
  <Type> Dirichlet </Type>
  <Time_dependence> Steady </Time_dependence>
  <Value> 1.0 </Value>
</Add_BC>
```
