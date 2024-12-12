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
#          journal          #
#############################


def getJournal(name=None, aktiv=False):
    global journal
    if not journal:
        if not name:
            name = 'none'
        journal = Journal(name=name, aktiv=aktiv)
    return journal


class Journal(object):
    def __init__(self,
                 name="journal",
                 aktiv=False):
        self.name = name
        self.aktiv = aktiv
        self.journal = {}
        self.data = {}
        self.data_key = None
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
        self.sort_by_keys()
        with open(self.filename, 'w') as f:
            f.write(json.dumps(self.journal, indent=4))
            f.close()

    def sort_by_keys(self):
        myKeys = list(self.data.keys())
        myKeys.sort()
        sorted_dict = {i: self.data[i] for i in myKeys}
        self.data = sorted_dict
        self.journal[self.data_key] = sorted_dict

    def get_journal(self, name):
        if not self.aktiv:
            return None

        if self.data:
            return self.data

        self.open_journal()
        self.data = self.journal.get(name, None)
        self.data = {'filename': ""}  # initialise
        self.data_key = name
        self.journal[name] = self.data
        return self.data

    def get_total(self, name):
        val = self.data.get(name, None)
        if val is None:
            return 0
        if isinstance(val, list):
            if not val:
                return 0
            if isinstance(val[0], int) or isinstance(val[0], float):
                val = sum(val)
                return val
            return len(val)
        if isinstance(val, int) or isinstance(val, float):
            return val
        return 1  # entries

    def set_benchmark(self):
        if self.get_total('area_errors') > 0:
            self.put_warning("Problem with areas")
        if self.get_total('appendices_deleted') > 0:
            self.put_warning("Problem with appendices")

    def set(self, name, val):
        if not self.data:
            return
        self.data[name] = val

    def put(self, name, val):
        if not self.data:
            return
        data = self.data.get(name, None)
        if data is None:
            self.data[name] = val
            return
        if isinstance(data, list):
            self.data[name].append(val)
            return
        data_list = [data]
        data_list.append(val)
        self.data[name] = data_list

    def set_filename(self, val):
        self.set('filename', val)

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

    def put_nodes_connected(self, val):
        self.put('nodes_connected', val)

    def put_exception(self, msg):
        self.put('exception', msg)

    def put_warning(self, msg):
        self.put('warning', msg)
