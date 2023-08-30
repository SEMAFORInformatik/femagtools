# -*- coding: utf-8 -*-
"""
    femagtools.conductor
    ~~~~~~~~~~~~~~~~~~~~

    Creating and managing coonductor material


"""


class Conductor:
    def __init__(self, clist=[]):
        """initialize this object from a list of conductor objects"""
        self.conductors = {}
        for c in clist:
            if 'id' in c:
                self.conductors[str(c['id'])] = c
            elif 'name' in c:
                self.conductors[c['name']] = c

    def find(self, id):
        """find conductor by id or name"""
        try:
            return self.conductors[id]
        except ValueError:
            pass  # not found
        except KeyError:
            pass
        return self.find_by_name(id)

    def find_by_name(self, name):
        """find conductor by name"""
        for k in self.conductors.keys():
            if self.conductors[k]['name'] == name:
                return self.conductors[k]
            # not found
        return None

    def __str__(self):
        return str([self.conductors[c]
                    for c in self.conductors])
