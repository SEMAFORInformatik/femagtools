# -*- coding: utf-8 -*-
"""
    femagtools.config
    ~~~~~~~~~~~~~~~~

    Config FEMAG

    :copyright: 2016 Semafor Informatik & Energie AG, Basel
    :license: BSD, see LICENSE for more details.
"""
try:
    # Load python 2 config parser
    import ConfigParser as configparser
except ImportError:
    import configparser


class Config(dict):

    def __init__(self, defaults=None):
        dict.__init__(self, defaults or {})

    def from_ini_file(self, ini_file):
        if not ini_file:
            return
        config = configparser.ConfigParser()
        config.read(ini_file)

        # Only get config for this engine
        engine = self['ENGINE']

        try:
            config[engine.lower()]
            engine = engine.lower()
        except KeyError:
            try:
                config[engine.upper()]
                engine = engine.upper()
            except KeyError:
                engine = None

        if not engine:
            return

        for key in config[engine]:
            if ',' in config[engine][key]:
                self[key.upper()] = config[engine][key].split(',')
                continue
            self[key.upper()] = config[engine][key]

    def from_ini_file_all(self, ini_file):
        if not ini_file:
            return

        config = configparser.ConfigParser()
        config.read(ini_file)

        for group in config:
            for key in config[group]:
                self[key.upper()] = config[group][key]
