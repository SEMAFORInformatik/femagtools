**Winding**
-----------

============    ============================  =======
Name             Parameter                    Default
============    ============================  =======
num_phases      number of phases (m)
num_wires       number of wires per slot
coil_span       coil span
num_layers      number of layers
resistance      elect. resistance (Ohm)
cufilfact       Fill factor of copper          0.45
culength        rel length of conductor        1.4
cuconduc        conductivity (S/m)             56e6
slot_indul      insulation thickness in slot   0.0
material        name of material
============    ============================  =======

.. Note:: The material parameter references a dict item in
  the list of :ref:`conductor` material. It overrides the cuconduc value and is the preferred method.


End-Winding Leakage
^^^^^^^^^^^^^^^^^^^

Windings may contain a leakage dict: leak_dist_wind, leak_evol_wind, leak_tooth_wind (version added 0.9.9)

* leak_dist_wind

  ============    ============================  =======
  Name             Parameter                    Unit
  ============    ============================  =======
  perimrad        Radius of perimeter            m
  vbendrad        Bending radius vertical        m
  endheight       End winding height             m
  m.wiredia       Wire diameter                  m
  ============    ============================  =======

* leak_evol_wind

  ============    =============================  =======
  Name             Parameter                     Unit
  ============    =============================  =======
  evol1rad        Top radius of first evolvent   m
  evol2rad        Top radius of second evolvent  m
  botlevel        Level at bottom of evolvents   m
  toplevel        Level at top of evolvents      m
  evolbend        Bending radius                 m
  endheight       End winding height             m
  m.wiredia       Wire diameter                  m
  ============    =============================  =======

* leak_tooth_wind

  ============    ============================  =======
  Name             Parameter                    Unit
  ============    ============================  =======
  bendrad         Bending radius vertical        m
  endheight       End winding height             m
  m.wiredia       Wire diameter                  m
  ============    ============================  =======

  Example::

    winding=dict(
        num_phases=3,
        num_wires=100,
        coil_span=3.0,
        num_layers=1,
        leak_dist_wind=dict(
            perimrad=67.1e-3, # Radius of perimeter [m]
            vbendrad=5e-3,    # Bending radius vertical [m]
            endheight=20e-3,  # End winding height [m]
            wiredia=1e-3)     # Wire diameter [m]
    )
