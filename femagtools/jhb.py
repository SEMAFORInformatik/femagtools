#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
import sys
import os
import femagtools.mcv as mcv
import logging

logger = logging.getLogger(__name__)


class Reader(object):
    def __init__(self, filename):
        self.mc1_type = mcv.MAGCRV

        with open(filename) as f:
            self.name = os.path.splitext(
                os.path.basename(f.readline().strip()))[0]
            r = f.readline().split()
            self.angle = [float(a) for a in r[1:]]
            self.curve = [dict(hi=[], bi=[])
                          for i in range(int(r[0])-1)]
            for l in f:
                r = l.split()
                if len(r) > 1:
                    bh = [float(x.replace(',', '.')) for x in r]
                    if bh[1] > 0:
                        b = bh[0]
                        for i, h in enumerate(bh[1:]):
                            self.curve[i]['hi'].append(h)
                            self.curve[i]['bi'].append(b)
            if len(self.curve) > 1:
                self.mc1_type = mcv.ORIENT_CRV
        logger.info("JHB %s size %d", self.name, self.curve[0]['bi'])
        
    def __getitem__(self, index):
        return self.__getattribute__(index)
    
    def getValues(self):
        """return values as mcv dict"""
        return {
            'name': self.name,
            'ctype': self.mc1_type,
            'curve': self.curve}


def read(filename):
    """read JHB File and return mc dict"""
    jhb = Reader(filename)
    return jhb.getValues()


if __name__ == "__main__":
    if len(sys.argv) == 2:
        filename = sys.argv[1]
    else:
        filename = sys.stdin.readline().strip()
            
    jhb = Reader(filename)
    print(jhb.getValues())
