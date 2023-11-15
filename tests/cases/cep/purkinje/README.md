
# **Problem Description**

Simulate electric signal propagation inside the Purkinje network.

<p align="center">
   <img src="./activation_10fps.gif" width="600">
</p>

The ten-Tusscher-Panfilov model is used to describe the cell activation. For details regarding the model, please refer to the following publications:

> K. H. W. J. ten Tusscher, D. Noble, P. J. Noble, and A. V. Panfilov. A model for hu-
> man ventricular tissue. American Journal of Physiology-Heart and Circulatory Physiology,
> 286(4):H1573–H1589, apr 2004.

> K. H. W. J. ten Tusscher and A. V. Panfilov. Alternans and spiral breakup in a human
> ventricular tissue model. American Journal of Physiology-Heart and Circulatory Physiology,
> 291(3):H1088–H1100, sep 2006.

The input file `svFSI.inp` follows the master input file [`svFSI_master.inp`](./svFSI_master.inp) as a template. More on Purkinje network generation and simulation can be found here:

- SimVascular Website: https://simvascular.github.io/docsSimCardio.html#purkinje
- Youtube Tutorial: https://www.youtube.com/watch?v=TCK3SmGwBa8&ab_channel=SimVascular
