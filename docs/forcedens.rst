ForceDensity
************

The ForceDensity object contains the PLT file results. It has
following attributes which mostly correspond to the text sections in the file:

================  =======================================================
Attribute          Description     
================  =======================================================
project            Name of model file
filename           Name of PLT file
date               calculation date
version            FEMAG version
nodes              number of nodes
elements           number of elements
quality            meshing quality
positions          Position properties
================  =======================================================

Positions
=========

list of dictionaries for each position

=========  =======================================================
position   rotor position
X          list of angles
FN         list of normal component of force
FT         list tangential component of force
Radius     list of radius
B_N        list normal component of flux density
B_T        list tangential component of flux density
=========  =======================================================

  Example::
    
   fdens = ForceDensity()
   fdens.read('example.PLT0')
   pl.title('{}, Rotor Position {}'.format(
            fdens.title, fdens.positions[0]['position']))
   pl.plot(fdens.positions[0]['X'], fdens.positions[0]['FT'])

.. plot:: pyplots/forcedens.py
   
