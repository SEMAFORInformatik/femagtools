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
        self.data = {'filename': ""}  # initialise
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

    def add(self, name, val):
        if self.data:
            self.data[name] = val + self.data.get(name, 0)

    def put_filename(self, val):
        self.put('filename', val)

    def put_areas(self, val):
        self.add('areas', val)

    def put_area_errors(self, val):
        self.add('area_errors', val)

    def put_elements(self, val):
        self.put('elements', val)

    def put_nodes(self, val):
        self.put('nodes', val)

    def put_concat_lines(self, val):
        self.add('concat_lines', val)

    def put_concat_arcs(self, val):
        self.add('concat_arcs', val)

    def put_appendices(self, val):
        self.add('appendices', val)

    def put_appendices_connected(self, val):
        self.add('appendices_connected', val)

    def put_appendices_remaining(self, val):
        self.add('appendices_remaining', val)

    def put_appendices_deleted(self, val):
        self.add('appendices_deleted', val)
       
    def put_exception(self, msg):
        self.put('exception', msg)

    def put_warning(self, msg):
        lst = self.data.get('warning', [])
        lst.append(msg)
        self.put('warning', lst)
