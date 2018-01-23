import pytest
import numpy as np
import femagtools.isa7

def test_read_file():
    filename = 'tests/data/minimal.ISA7'
    isa = femagtools.isa7.read(filename)
    
    assert isa.NUM_PNT == 39
    assert isa.NUM_LIN == 39
    assert isa.NUM_SPEL_NDCH == 112
    assert type(isa.POINT_ISA_PT_VALID[0]) == bool
    assert type(isa.POINT_ISA_POINT_REC_PT_CO_X[0]) == float

def test_lines():
    filename = "tests/data/minimal.ISA7"
    isa = femagtools.isa7.read(filename)

    assert len(isa.lines()) == isa.NUM_LIN

def test_se_outline():
    filename = "tests/data/minimal.ISA7"
    isa = femagtools.isa7.read(filename)
    ol = isa.se_outline()

    assert np.array(ol).shape == (isa.NUM_SPEL_NDCH, 2, 2)
    