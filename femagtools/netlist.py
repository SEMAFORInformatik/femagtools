"""Netlist creation """

def create(num_poles, model_divisor, voltage_source, rotor_cage):
    rotor_slots_gen = rotor_cage['num_slots'] // model_divisor
    num_poles_gen = num_poles//model_divisor
    NoPhases = len(voltage_source['phiq'])
    StatorNodes=2*NoPhases+1

    uq = [dict(
        Type = "Uac",
        Name = f"Uq{i+1}",
        Nodes = [i+1, 0],
        Amplitude = voltage_source['Uq_amp'],
        Frequency = voltage_source['fq'],
        Phaseangle = phiq,
        Ri = voltage_source['Rqi'])
              for i, phiq in enumerate(voltage_source['phiq'])]

    hf = [dict(
        Type = "Uac",
        Name = f"Uq{i+1}_HF",
        Nodes = [i+1,i+NoPhases+1],
        Amplitude = 0,
        Frequency = 0,
        Phaseangle = 0,
        Ri = 0) for i in range(NoPhases)]

    stator_wdg = [dict(
        Type = "Wdg",
        Name = f"PH{i+1}",
        Nodes = [i+NoPhases+1,StatorNodes],
        WdgKey = i+1) for i in range(NoPhases)]
    
    rotor_bars = [dict(
        Type = "Wdg",
        Name = f"Bar{i+1}",
        Nodes = [i+StatorNodes, i+StatorNodes+rotor_slots_gen],
        WdgKey = i+1+NoPhases ) for i in range(rotor_slots_gen)]
    
    rotor_bars[0]['Nodes'] = [0, StatorNodes+rotor_slots_gen] 

    rotor_ring_segments = [dict(
        Type = "RL",
        Name = f"Z{i+1}",
        Nodes = [i+StatorNodes, i+StatorNodes+1],
        R = rotor_cage['RR']*model_divisor,
        L = rotor_cage['LR']*model_divisor)
                               for i in range(rotor_slots_gen*2)]
        
    rotor_ring_segments[0]['Nodes'] = [0,StatorNodes+1]
    if num_poles_gen % 2 == 0:   # even number of poles in model
        i = rotor_slots_gen-1
        rotor_ring_segments[i]['Nodes'] = [i+StatorNodes, 0]
        i = 2*rotor_slots_gen-1
        rotor_ring_segments[i]['Nodes'] = [i+StatorNodes, StatorNodes + rotor_slots_gen]
    else:
        i = 2*rotor_slots_gen-1
        rotor_ring_segments[i]['Nodes'] = [i+StatorNodes, 0]

    return dict(
        Branches=uq+hf+stator_wdg+rotor_bars+rotor_ring_segments )
