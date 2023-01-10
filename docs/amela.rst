Amela
*****

AMELA can be executed as a batch process in femagtools. Amela object will search the 
AMELA batch file and the FEMAG model file in the given directories. The magnet data
for the loss calculation will be extracted from the FEMAG model file. 
  * If the FEMAG model file (nc file) is located in the same directory as AMELA.
    Only the workdir is needed. 

  Example::
   
   amela = Amela(workdir, dict(name='example'))
   r = amela()
   loss = r['pm_data_se38']['total_loss']
   print(f'loss in the superelement 38 is: {loss} W')

  * If the FEMAG model file (nc file) is located in a different directory as AMELA.
    Both the workdir (nc file) and amela_dir need to be input.

  Example::
   
   amela = Amela(workdir, dict(name='example'), \
                amela_dir)
   r = amela()
   loss = r['pm_data_se38']['total_loss']
   print(f'loss in the superelement 38 is: {loss} W')

The circumferential and axial segmentation are optional parameters and can be
passed in the dictionary. 

  Example::

   magnet_data = dict(name='example', 
                      nseglen=3)
   amela = Amela(workdir, magnet_data, \
                amela_dir)
list of all the optional parameters: 

=========  ===================================
speed      rotational speed (1/min) 
sigma      electrical conductivity (S/m)
nsegwid    number of circumferential segments
nseglen    number of axial segments
=========  ===================================

The return values of the Amela object is a dictionary that contains the loss data: 

==========  ========================================================
loss_data   time dependable loss value
total_loss  total loss in the superelement (averaged over time step)
==========  ========================================================