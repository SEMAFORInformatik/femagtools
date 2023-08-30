#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
import femagtools.magnet

magMat = [dict(
    name="M395",
    remanenc=1)
]


def test_findById():
    mag = femagtools.magnet.Magnet(magMat)
    result = mag.find('M395')
    expected = magMat[0]
    assert result == expected
