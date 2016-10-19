#!/usr/bin/env python
#
# convert bch to xml
#
# Author: Ronald Tanner
# 2016-10-20
#
import sys
import re
import io
from .bch import Reader
import xml.etree.ElementTree as el
import xml.dom.minidom


def usage():
    print("femagtools.bchxml <bchfilename>")
    sys.exit(1)


def list_to_xml(tag, l):
    '''
    Turn a list into XML
    '''
    elem = el.Element(re.sub(r'^(\d+)$', r'w\1', tag))
    for v in l:
        if v:
            if isinstance(v, list):
                elem.append(list_to_xml('val', v))
            elif isinstance(v, dict):
                elem.append(dict_to_xml('val', v))
            else:
                child = el.Element('val')
                child.text = str(v)
                elem.append(child)
    return elem


def dict_to_xml(tag, d):
    '''
    Turn a simple dict of key/value pairs into XML
    '''
    
    elem = el.Element(re.sub(r'^(\d+)$', r'w\1', tag))
    for key, val in d.items():
        if val:
            if isinstance(val, dict):
                elem.append(dict_to_xml(key, val))
            elif isinstance(val, list):
                elem.append(list_to_xml(key, val))
            else:
                child = el.Element(key)
                child.text = str(val)
                elem.append(child)
    return elem

if __name__ == "__main__":

    if len(sys.argv) != 2:
        usage()
        
    filename = sys.argv[1]
    bchresults = Reader()
    with io.open(filename) as f:
        bchresults.read(f.readlines())

    reparsed = xml.dom.minidom.parseString(el.tostring(
        dict_to_xml("bch", bchresults), method='xml'))

    with io.open(filename.split('.')[0]+'.xml',
                 mode='w', encoding='utf-8') as f:
        f.write(reparsed.toprettyxml(indent='  '))
