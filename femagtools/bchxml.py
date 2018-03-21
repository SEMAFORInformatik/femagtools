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

    if not isinstance(tag, str):
        try:
            tag = str(tag)
        except Exception as e:
            print("bchxml.list_to_xml():  Can not convert tag {} to str!".format(tag))
            pass

    elem = el.Element(re.sub(r'^(\d+)$', r'w\1', tag))
    for v in l:
        if v is not None:  # check against None explicitly to keep ZERO values. thomas.maier/OSWALD
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

    if not isinstance(tag, str):
        try:
            tag = str(tag)
        except Exception as e:
            print("bchxml.dict_to_xml():  Can not convert tag {} to str!".format(tag))
            pass

    elem = el.Element(re.sub(r'^(\d+)$', r'w\1', tag))
    for key, val in d.items():
        if val is not None:  # check against None explicitly to keep ZERO values. thomas.maier/OSWALD
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
    with io.open(filename, encoding='latin1', errors='ignore') as f:
        bchresults.read(f.readlines())

    reparsed = xml.dom.minidom.parseString(el.tostring(
        dict_to_xml("bch", bchresults), method='xml'))

    with io.open(filename.split('.')[0]+'.xml',
                 mode='w', encoding='utf-8') as f:
        f.write(reparsed.toprettyxml(indent='  '))
