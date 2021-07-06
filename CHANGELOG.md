# Changelog

This document lists the major changes in femagtools. Please clone this project to follow bug fixes and minor enhancements.


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
