FEMAG
*****

FEMAG can be used in 2 ways:

* as a standalone calculation program that is invoked with the filename of
  a FSL script::

   workdir = os.path.join(os.path.expanduser('~'), 'femag')
   femag = femagtools.Femag(workdir)
   femag.run('femag.fsl')
  
* as a service running on a host that is reachable by TCP/IP::

   with open('femag.fsl') as f:
        fslcmds = f.read()
	
   femag = femagtools.femag.ZmqFemag(5555, 'localhost')
   ret = femag.send_fsl(fslcmds)
   r = [json.loads(s) for s in ret]

.. Note::
   
   The communication is based on ZMQ (see http://zeromq.org). This requires a Femag version with ZMQ support.
   Contact support@profemag.ch for further information.

The return value of the send_fsl method is a list of 2 strings which can be converted to a dict as shown
in the above example. In the first Dict has always a status and a message field. The status field
contains the string 'error' if any error occured during the fsl processing. Otherwise the value is 'ok'.
Further details can be found in the string value message field. The structure of the second dict depends on
the type of the calculation. It can also be empty.

When running in service mode the input and output files can usually not be accessed directly.
Instead these files must be transferred::

    transfer_files = ['M330-50A.MCV', 'PM_270_L8_8p.poc',
                      'PM_270_L8.AUX7', 'PM_270_L8.ISA7']

    ret = femag.upload(transfer_files)
  

If Femag has created result files during the previous FSL processing the first dict of the send_fsl return value
includes the field 'result_file' which holds a list of filenames. These files can be downloaded
with the method getfile::

    if r[0]['status'] == 'ok':
        try:
            bchfile = r[0]['result_file'][0]
            status, content = femag.getfile(bchfile)
            bch = femagtools.bch.Reader()
            bch.read(content.decode())
            print('Torque [Nm] = {}'.format(bch.machine['torque']))

        except KeyError:
            logger.warn('No result file')
    else:
        logger.error(r[0]['message'])

It is recommended to not fiddling with these files directly but to use
the ZmqFemag call method instead which handles these details internally::
  
     femag = femagtools.ZmqFemag(5555, workdir=workdir,
                            magnets=magnetmat,
                            magnetizingCurves=mcv)

     machine = dict( .. )

     simulation = dict(
         calculationMode="cogg_calc",
         magn_temp=60.0,
         num_move_steps=49,
         speed=50.0)

     r = femag(machine, simulation)

     print("Order    T/Nm      %")
     tq = r.torque_fft[-1]
     for l in zip(tq['order'], tq['torque'], tq['torque_perc']):
         print('{0:<5} {1:9.2f} {2:6.1f}'.format(*l))

A complete example can be found in the directory examples/docker-zmq.
