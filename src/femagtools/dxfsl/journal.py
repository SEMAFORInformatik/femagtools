# -*- coding: utf-8 -*-
"""
  femagtools.dxfsl.concat
  ~~~~~~~~~~~~~~~~~~~~~~~

  Authors: Ronald Tanner, Beat Holm
"""
from __future__ import print_function

import numpy as np
import logging
import sys
import json
import time
from pathlib import Path

logger = logging.getLogger('femagtools.journal')
journal = None

#############################
#           concat          #
#############################


def getJournal(name=None, aktiv=False):
    global journal
    if not journal:
        if not name:
            journal = Journal(name='none', aktiv=False)
        else:
            journal = Journal(name=name, aktiv=True)
    return journal


class Journal(object):
    def __init__(self,
                 name="journal",
                 aktiv=False):
        self.name = name
        self.aktiv = aktiv
        self.journal = {}
        self.data = {}
        self.filename = Path('{}.json'.format(self.name))

    def open_journal(self):
        if self.filename.is_file():
            with open(self.filename, 'r') as f:
                try:
                    self.journal = json.loads(f.read())
                except:
                    self.journal = {}
                f.close()

    def write_journal(self):
        if not self.journal:
            return
        self.set_benchmark()

        with open(self.filename, 'w') as f:
            f.write(json.dumps(self.journal, indent=4))
            f.close()

    def get_journal(self, name):
        if not self.aktiv:
            return None

        if self.data:
            return self.data

        self.open_journal()
        self.data = self.journal.get(name, None)
        if not self.data:
            self.data = {'filename': "",
                         'elements': 0,
                         'nodes': 0,
                         'concat_lines': 0,
                         'concat_arcs': 0,
                         'appendices': 0,
                         'appendices_connected': 0,
                         'appendices_remaining': 0,
                         'areas': 0,
                         'area_errors': 0}
            self.journal[name] = self.data
        return self.data

    def set_benchmark(self):
        if self.data.get('area_errors', 0) > 0:
            self.put_warning("Problem with areas")
        if self.data.get('appendices_deleted', 0) > 0:
            self.put_warning("Problem with appendices")
          
    def put(self, name, val):
        if self.data:
            self.data[name] = val

    def put_filename(self, val):
        self.put('filename', val)

    def put_areas(self, val):
        self.put('areas', val)

    def put_area_errors(self, val):
        self.put('area_errors', val)

    def put_elements(self, val):
        self.put('elements', val)

    def put_nodes(self, val):
        self.put('nodes', val)

    def put_concat_lines(self, val):
        self.put('concat_lines', val)

    def put_concat_arcs(self, val):
        self.put('concat_arcs', val)

    def put_appendices(self, val):
        self.put('appendices', val)

    def put_appendices_connected(self, val):
        self.put('appendices_connected', val)

    def put_appendices_remaining(self, val):
        self.put('appendices_remaining', val)

    def put_appendices_deleted(self, val):
        self.put('appendices_deleted', val)
       
    def put_exception(self, msg):
        self.put('exception', msg)

    def put_warning(self, msg):
        print(msg)
        lst = self.data.get('warning', [])
        lst.append(msg)
        self.put('warning', lst)
