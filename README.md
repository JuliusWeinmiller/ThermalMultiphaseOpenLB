# ThermalMultiphaseOpenLB
This is an extension for OpenLB which allows for thermal single-component phase change multiphase flow

This extension was created by Julius Weinmiller as part of his master thesis.
The final thesis will be uploaded shortly to TUDelft repository, and once avaliable will be linked here.

## Features
With this extension, it is possible to simulate phase change using OpenLB of liquids at high density ratios.

An example usecase is flow boiling in microchannels.


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
  * Navigate to the class `BounceBack` (line 652)
  * Change `private` to `protected`
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

