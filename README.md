**SPARTACUS** is a Unified Stochastic Particle (*USP*) solver based on SPARTA.
**SPARTACUS** stands for **SPARTA** **C**ombined with **US**P.
The *USP* method is proposed by *Fei et al.* to simulate multi-scale gas flows.

> [examples_USP](./examples_USP/) contains Several benchmark cases and tutorial examples.

> [doc/collide.html](./doc/collide.html) & [doc/collide_bgk_modify.html](./doc/collide_bgk_modify.html)  contains essential documentation about this solver.

Downloading, installing and using SPARTACUS is exactly the same as SPARTA, which is fully described within [SPARTA Users Manual](https://sparta.github.io/doc/Manual.html), unless specified in the above files and directories

*Coupled time step and grid adaptation based on local characteristic scale is available now!*

![](./doc/JPG/adaptation.gif)

**SPARTACUS Papers**

The following paper describes the *USP* method and its implementation within the SPARTACUS framework. 

> Feng, K., Tian, P., Zhang, J., Fei, F., Wen, D., 2022. [**SPARTACUS: An open-source unified stochastic particle solver for the simulation of multiscale nonequilibrium gas flows.**](https://www.sciencedirect.com/science/article/abs/pii/S0010465522003265) Computer Physics Communications 108607. https://doi.org/10.1016/j.cpc.2022.108607

The following are some references on the theory and application of the *USP* method.

> Fei, F., Zhang, J., Li, J., Liu, Z., 2020. [**A unified stochastic particle Bhatnagar-Gross-Krook method for multiscale gas flows.**](https://doi.org/10.1016/j.jcp.2019.108972) Journal of Computational Physics 400, 108972. 

> Fei, F., Jenny, P., 2021. [**A hybrid particle approach based on the unified stochastic particle Bhatnagar-Gross-Krook and DSMC methods.**](https://doi.org/10.1016/j.jcp.2020.109858) Journal of Computational Physics 424, 109858. 

> Zhang, J., John, B., Pfeiffer, M., Fei, F., Wen, D., 2019.[ **Particle-based hybrid and multiscale methods for nonequilibrium gas flows.**](https://doi.org/10.1186/s42774-019-0014-7) Adv. Aerodyn. 1, 12. 

> Fei, F., Hu, Y., Jenny, P., 2022. [**A unified stochastic particle method based on the Bhatnagar-Gross-Krook model for polyatomic gases and its combination with DSMC.**](https://doi.org/10.1016/j.jcp.2022.111640) Journal of Computational Physics 471, 111640. 

> Fei, F., Jenny, P., 2022. [**A high-order unified stochastic particle method based on the Bhatnagar-Gross-Krook model for multi-scale gas flows.**](https://doi.org/10.1016/j.cpc.2022.108303) Computer Physics Communications 274, 108303. 

> Zhang, J., Yao, S., Fei, F., Ghalambaz, M., Wen, D., 2020. [**Competition of natural convection and thermal creep in a square enclosure.**](https://doi.org/10.1063/5.0022260) Physics of Fluids 32, 102001. 


----------------------------------------------------------------------
----------------------------------------------------------------------

SPARTA on [GitHub](https://github.com/sparta/sparta): https://github.com/sparta/sparta 


***Original SPARTA README:***



This is the SPARTA software package.

SPARTA stands for Stochastic PArallel Rarefied-gas Time-accurate
Analyzer.

Copyright (2014) Sandia Corporation.  Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
certain rights in this software.  This software is distributed under
the GNU General Public License.

----------------------------------------------------------------------

SPARTA is a Direct Simulation Monte Carlo (DSMC) code designed to run
efficiently on parallel computers.  It was developed at Sandia
National Laboratories, a US Department of Energy facility, with
funding from the DOE.  It is an open-source code, distributed freely
under the terms of the GNU Public License (GPL).

The primary authors of the code are Steve Plimpton and Michael Gallis,
who can be emailed at sjplimp@sandia.gov and magalli@sandia.gov.  The
SPARTA web site at http://sparta.sandia.gov has more information about
the code and its uses.

The SPARTA distribution includes the following files and directories:

README			   this file
LICENSE			   the GNU General Public License (GPL)
bench                      benchmark problems
data                       files with species/reaction params, surface files
doc                        documentation
examples                   simple test problems
lib                        additional library files
python                     Python wrapper on SPARTA as a library
src                        source files
tools                      pre- and post-processing tools

Point your browser at any of these files to get started:

doc/Manual.html	           the SPARTA manual
doc/Section_intro.html	   hi-level introduction to SPARTA
doc/Section_start.html	   how to build and use SPARTA
doc/Developer.pdf          SPARTA developer guide
