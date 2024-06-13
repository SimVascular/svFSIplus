
# **Problem Description**

Solve fluid flow in a cylindrical tube with resistance or RCR boundary conditions at the outlet and unsteady flow at inlet.

The input file `svFSI.inp` follows the master input file [`svFSI_master.inp`](./svFSI_master.inp) as a template. Some specific input options are discussed below:

## Resistance and RCR values

If the problem is solved using resistance BC at the outlet, the following keywords are expected in `Add BC`:

`Type: Neu
 Time dependence: Resistance
 Value: << resistance_value >>`

Note that when providing resistance_value, << >> are not required.

If the problem is solved using an RCR BC at the outlet, the following keywords are expected in `Add BC`:

`Type: Neu
 Time dependence: RCR   # or Windkessel
 RCR values: << (proximal_resistance, capacitance, distal_resistance) >>
 Distal pressure: << distal_pressure_value >>`

For prescribing the RCR BC, the order of the values is fixed i.e., proximal resistance followed by capacitance and the distal resistance. These values can be separated by comma or space, and should be enclosed within parentheses (), double quotes "", or <>.

## Options for providing unsteady BCs

This example provides two ways of specifying unsteady boundary conditions:

(a) Providing time-dependent data as input using the keyword,

`Temporal values file path: lumen_inlet.flow`

(b) Providing interpolated Fourier coefficients as input using the keyword,

`Fourier coefficients file path: lumen_inlet.fcs`

### File format for time-dependent data

The format for providing time-dependent data in an ASCII formatted file is:

`<Line 1>  number_of_time_points    number_of_Fourier_modes
 <Line 2>  time_point_1     data_value_1
 <Line 3>  time_point_2     data_value_2
 .
 .
 .
 <Line n>  time_point_n     data_value_n`

### File format for Fourier coefficients

The format for providing Fourier coefficients in an ASCII formatted file is:

`<Line 1>  first_time_point
 <Line 2>  last_time_point
 <Line 3>  first_data_value
 <Line 4>  first_Fourier_mode
 <Line 5>  number_of_Fourier_modes
 <Line 6>  real_Fourier_mode_1   imag_Fourier_mode_1
 <Line 7>  real_Fourier_mode_2   imag_Fourier_mode_2
 .
 .
 .
 <Line N+5>  real_Fourier_mode_N   imag_Fourier_mode_N`

In the above format, the first_Fourier_mode on Line 4 is to be computed as,

```bash
first_Fourier_mode = (last_data_value - first_data_value) / (last_time_point - first_time_point)
```
