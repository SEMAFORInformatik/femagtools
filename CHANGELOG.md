# Changelog

This document lists the major changes in femagtools. Please clone this project to follow bug fixes and minor enhancements.

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
