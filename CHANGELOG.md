# Changelog

	This document lists the major changes in femagtools. Please clone this project to follow bug fixes and minor enhancements.
## Release 1.7.7
	- dxfsl: new machine type EESM

## Release 1.7.4
	- added optional parameter dqpar for psd_psq or ld_lq (default) identification in pm parident
	- added method transient for PMRelMachinePsidq (experimental)

## Release 1.7.0
	- displ_stator_rotor: eccentricity simulation and analysis
	- improved force density plots (contour and surface)

## Release 1.6.4
	- dxfsl: split option is deprecated

## Release 1.6.1
        - calc torque frequency spectrum in bch cogging if missing

## Release 1.6.0
        - iron loss calculation: support for Bertotti- and Modified-Steinmetz Loss-Models added.
        - New FSL templates for custom calculation via PVFE_FSL
        - new loss parameter (ce, b_beta_coeff) in MC/MCV file
        - bch: new loss section (excess losses) in the BCH file (FEMAG version >= R2024 B2).
	- parameter identification from dxf (requires i1_max, speed)
	- added fieldlines plot (from svg)
	- added fieldcalc with result airgap and field_lines
	- dxf fsl improvements

## Release 1.5.5
	- PM, SM characteristics: auto correction of max torque

## Release 1.5.2
	- VTU demag access by elements added

## Release 1.4.8
	- Core loss calc method Bertotti added

## Release 1.4.5
        - TH properties in templates statorRotor3, magnetSector, magnetIron, magnetIronV
	- netCDF version update for CVE-2023-38545

## Release 1.4.4
	- model dict key 'windings' renamed to 'winding' (by keeping backward compatibility)

## Release 1.4.2
	- custom loss calc func for pm_sym_fast, torq_calc, ld_ld, psd_psq simulations, mult_cal

## Release 1.4.1
	- convert from JMAG/JPLOT to msh
	- added parident of AFM

## Release 1.4
	- new machine type AFM (Axial Flux Machine)
        - added isa/nc methods (properties): get_areas, elements (inner,outer,center), calc_iron_loss
	- split plot module into subpackage
        - extended magnet_loss calc (IALH method)

## Release 1.3.1
	- new flags: with_mtpv, with_mtpa, with_pmconst, with_tmech for effloss and pm characteristics
        - new plot func transientsc_demag
	- new hxy file read module
	- characteristics with friction and windage losses (with_tmech)
	- effmap with multiprocessing support
	- improved dxf-fsl conversion
        - improved contour plots
	- demag_pos plot with pos at max demag
        - stator resistance: frequency dependent calculation with user defined function

## Release 1.3
	- TH module for static and dynamic thermal simulations included

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
