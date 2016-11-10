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

    def __init__(self, defaults=None, configfile=None):
        dict.__init__(self, defaults or {})

    def from_ini_file(self, ini_file):
        if not ini_file:
            return
        config = configparser.ConfigParser()
        config.read(ini_file)

        # Only get config for this engine
        section = self['ENGINE']

        try:
            if config.items(section.lower()):
                engine = section.lower()
            else:
                if config.items(section.upper()):
                    engine = section.upper()
                else:
                    return
        except configparser.NoSectionError:
            try:
                config[section.upper()]
                engine = section.upper()
            except:
                return

        for key, value in config.items(engine):
            if ',' in value:
                self[key.upper()] = value.split(',')
                continue
            self[key.upper()] = value

    def from_ini_file_all(self, ini_file):
        if not ini_file:
            return

        config = configparser.ConfigParser()
        config.read(ini_file)

        for group in config:
            for key in config[group]:
                self[key.upper()] = config[group][key]
