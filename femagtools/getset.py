# -*- coding: utf-8 -*-
"""
    femagtools.getset
    ~~~~~~~~~~~~~~~~~

    Generic data access



"""
import re

_indxOpPattern = re.compile(r'[a-zA-Z0-9_]+(\[[-0-9]+\])')


def splitindex(name):
    "returns listname and index if name contains index"
    m = _indxOpPattern.search(name)
    if m:
        return (name.split('[')[0],
                int(''.join(list(m.group(1))[1:-1])))
    return (name, None)


class GetterSetter(object):
    """Generic data access"""
    def __init__(self, parameters):
        if isinstance(parameters, dict):
            for k in parameters.keys():
                setattr(self, k, parameters[k])

    def set_value(self, name, value, p=None):
        """set value of parameter identified by name

        Args:
            name: name of parameter (can be list or string)
            value: value to be assigned to parameter
        """
        if isinstance(name, str):
            setattr(self, name, value)
            return

        if len(name) > 1:
            k = name[0]
            if hasattr(self, k):
                self.set_value(name[1:], value, getattr(self, k))
                return
            elif p:
                self.set_value(name[1:], value, p[k])
                return
            self.set_value(name[1:], value, self)
            return
        if p:
            p[name[0]] = value
            return
        setattr(self, name[0], value)


    def get(self, name, r=None):
        """return value of key name
        name can be a list such as ['torque[1]', 'ripple']
        or a string: 'dqPar'
        """
        try:
            if isinstance(name, str):  # ignore r
                lname, indx = splitindex(name)
                if indx is None:
                    return getattr(self, name)
                return getattr(self, lname).__getitem__(indx)

            if len(name) > 1:
                lname, indx = splitindex(name[0])
                if r:
                    if indx is None:
                        return self.get(name[1:],
                                        getattr(r, lname))
                    return self.get(name[1:],
                                    getattr(r, lname).__getitem__(indx))

                if indx is None:
                    return self.get(name[1:],
                                    getattr(self, lname))
                return self.get(name[1:],
                                getattr(self, lname).__getitem__(indx))

            lname, indx = splitindex(name[0])
            if r:
                if indx is None:
                    return r.get(name[0])
                return r.get(lname).__getitem__(indx)

            if indx is None:
                return getattr(self, name[0])
            return getattr(self, lname).__getitem__(indx)
        except (KeyError, IndexError, AttributeError):
            return None
