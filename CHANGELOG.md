# Changelog

	This document lists the major changes in femagtools. Please clone this project to follow bug fixes and minor enhancements.

## Release 1.2.2
	- femag: temperature dependency of magnet material improved (muer)
    - windings: check unbalanced windings
    - effmap: drive mode only

## Release 1.1.26
	- amela interface added

## Release 1.1.25
        - ts: revised ts losses

## Release 1.1.24
	- bch: added skewAngle, skewSections, parallelWdgs to machine dict

## Release 1.1.23
	- parameter ident for im: noload calc with rotate
	- efficiency and losses maps
        - sizing added: spm, ipm, eesm, im
        - new parameter in machine dict to set the number of nodes in airgap: num_agnodes
## Release 1.1.22
	- mcv loss data handling

## Release 1.1.21
	- recalc mcv for dynamic (eddy current) simulations

## Release 1.1.20
	- extra_require: vtk, dxfgrabber, networkx, meshio, pyzmq, matplotlib

## Release 1.1.19
        - multiproc engine with timestep parameter for progress logging

## Release 1.1.18
	- added templatedirs parameter for user specific mako templates
        - improved TS support
        - added progress logging

## Release 1.1.17
	- rearranged module machine: im.py, sm.py, pm.py utils.py

## Release 1.1.12
	- ts module added (vtu postprocessing, loss calculation)

## Release 1.1.10
- vtu _movie support in mult_cal_fast
- added optional filename parameter in FE simulation ['plots']

## Release 1.1.3
- eesm support improved (bch, ld_ld_fast)

## Release 1.1.2
- parstudy with induction motor and wdg material

## Release 1.1.1
- dfx fsl conversion improvements
- parstudy psid-psiq ident with example

## Release 1.1
- new module parstudy with sampling methods for parameter variation: List, Grid, Sobol, LatinHypercube
- dakota integration (sampling, uq, moga, moat)
- nc/isa7: delta_node_angle

## Release 1.0.20
- femag: run femag with model par only
- winding: custom winding def
- isa7.flux_density, demagnetization: changed parameter x,y to el, removed cosys
- conductor: material properties added

## Release 1.0.19
- added templates: srm, rotorKs2, magnetShell, com_motor_sim

## Release 1.0.18
- added losses parts stator, rotor for model specific iron subregions
- fixed cosys bug in isa/nc flux_density

## Release 1.0.17
- improved meshing with auxiliary lines in dxf fsl conversion
- added Bch reader and plot for tubular machines (R/Z Coordinate System)
- fixed asm parameter ident for delta connected windings

## Release 1.0.16
- windings with zoneplan and mmf
- AC simulation with EEC parameter identification

## Release 1.0.15

- windings (experimental): create winding from scratch, winding diagrams, current linkage
- parameter variation: added values vector as an alternative to steps/bounds in decision vars
- AC simulation (experimental): asyn_motor
- Characteristic curves: losses approximation (from psid-psiq, ld-lq identification)
- poc file handling
. demagnetization analysis enhancements: per rotor position
- TS simulation: netlist handling
- phasor plot new layout
