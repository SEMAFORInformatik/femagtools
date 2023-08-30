#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
import femagtools.conductor

condMat = [dict(
    name="Cu",
    spmaweight=8.96,
    elconduct=56e6,
    tempcoef=3.9e-3)
]


def test_findById():
    cond = femagtools.conductor.Conductor(condMat)
    result = cond.find('Cu')
    expected = condMat[0]
    assert result == expected
