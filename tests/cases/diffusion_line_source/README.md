
# **Problem Description**

Solve heat transfer problem inside a 2D plate. A steady, uniform temperature value is specified on the left boundary, and zero temperature is prescribed on the right boundary. An option to provide unsteady temperature profile on the left boundary is also available. Adiabatic boundary conditions (zero heat flux) are employed on top and bottom boundaries. The result shows a linear distribution of temperature.

<p align="center">
   <img src="./plot-temp-flux.png" width="1000">
</p>

The input file `svFSI.inp` follows the master input file [`svFSI_master.inp`](./svFSI_master.inp) as a template. Some specific input options are discussed below:

## Steady Dirichlet BC

If a line source with constant temperature is specified on the left boundary, the following directives are used

```
   Add BC: left {
      Type: Dir
      Time dependence: Steady
      Value: 10.0
      Zero out perimeter: f
   }
```

Here, `Zero out perimeter` is only effective when used alongside Dirichlet BC. When set to `true`, this setting will enforce the temperature on top-left and bottom-left corners to be zero.

## Unsteady Dirichlet BC

In some scenarios, apply a large Dirichlet BC value at first time step with cause stability or convergence issues in the software. It is sometimes recommended to gradually ramp up the BC value to the desired quantity over a short period of time. To achieve this, the following directives are used in the input file:

```
   Add BC: left {
      Type: Dir
      Time dependence: Unsteady
      Temporal values file path: bc_left.dat
      Ramp function: t
   }
```

Once the ramp function is activated, the software will look for a user-specified file through ` Temporal values file path` setting. In this case, `bc_left.dat` has only three lines:
```
<Line 1> 2    1      # Never change for Ramp function
<Line 2> 0.0  0.0    # time_point_1     data_value_1
<Line 3> 1.0  10.0   # time_point_2     data_value_2
```

When ramp function is used, the software will linearly increase the boundary value from `data_value_1` to `data_value_2`  between `time_point_1` and `time_point_2`. After `time_point_2`, the boundary value will remain constant at `data_value_2`.
