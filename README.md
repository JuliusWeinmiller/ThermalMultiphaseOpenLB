# ThermalMultiphaseOpenLB
This is an extension for OpenLB which allows for thermal single-component phase change multiphase flow

This extension was created by Julius Weinmiller as part of his master thesis.
The final thesis will be uploaded shortly to TUDelft repository, and once avaliable will be linked here.

## Features
With this extension, it is possible to simulate phase change using OpenLB of liquids at high density ratios.

An example usecase is flow boiling in microchannHelmholtz Institute Ulmels.


### Code extensions
The extension includes:
* A new forcing algorithm for pseudopotential method
  * stable simulations high density ratios
  * thermodynamic consistent densities
* Multiphase thermal coupling DDF
  * Incorporation of latent heat via the equation of state
* Semi-hybrid thermal solver
  * Solves thermal equation with explicit euler
  * Helps with setting thermal diffusivity

## Limitations
The original extension was created for OpenLB v1.3-1.

It is currently not working in combination with OpenLB v 1.4, but it is being worked on to make sure it functions as well

## Quick Installation
1. Drop the `phaseChangeExtension` folder into the OpenLB base folder
2. In the `global.mk` file:
 * add the line `phaseChangeExtension/phaseChangeSource \` to the `SUBDIRS` and `INCLUDEDIRS` 
 * add `-lgsl -lgslcblas` to the `LIBS` after the `-lz`
3. One small change needs to be made in the file `src/dynamics/dynamics.h`
  * Navigate to the class `BounceBack` (line 651)
  * Change `private` to `protected` (line 701)
4. Compile and run the examples to make sure everything works
  * Do not forget to change `BUILDTYPE` to `generic` in `config.mk`


### Cleaner Installation
It is possible to organize the folders better, but it may break some references.

1. Place the source folder into the `src/external/` folder
  * Add the location of the source folder to `global.mk`
2. Place the python folder into the `src/external/` folder
  * The python files (Used only to generate plots) in the example program need to be given the proper reference to the `eos.py` file
3. Place the phaseChangeExamples into the OpenLB's examples folder
  * If the examples do not compile, make sure that the `ROOT` in the `definitions.mk` is correct.
  * Optional: add the examples to the `EXAMPLEDIRS` in the `global.mk`
4. Any other steps in the quick installation also need to be performed
  * Reminder to add: `-lgsl -lgslcblas`
  * Reminder to edit `BounceBack`

# Core papers used
Q. Li, K. H. Luo, and X. J. Li. 
Forcing scheme in pseudopotential lattice Boltzmann model for multiphaseflows.
Physical Review E, 86(1):016709, jul 2012. doi:10.1103/PhysRevE.86.016709.

S. Gong and P. Cheng.
A lattice Boltzmann method for simulation of liquid–vapor phase-change heat transfer.
International Journal of Heat and Mass Transfer, 55(17-18):4923–4927, aug 2012. doi:10.1016/j.ijheatmasstransfer.2012.04.037

Q. Li, P. Zhou, and H. J. Yan.
Improved thermal lattice Boltzmann model for simulation of liquid-vaporphase change.
Physical Review E, 96(6), 2017. doi:10.1103/PhysRevE.96.063303

## Relevant papers for MRT collision operator

Q. Li, K. H. Luo, and X. J. Li.
Lattice Boltzmann modeling of multiphase flows at large density ratio withan improved pseudopotential model.
Physical Review E, 87(5):053301, may 2013. doi:10.1103/PhysRevE.87.053301.

## Other interesting papers

Q. Li, Y. Yu, and Z. X. Wen.
How does boiling occur in lattice Boltzmann simulations?
093306(May),2020. doi:10.1063/5.0015491.

