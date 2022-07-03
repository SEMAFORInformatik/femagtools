# -*- coding: utf-8 -*-
"""
    femagtools.config
    ~~~~~~~~~~~~~~~~~

    Config FEMAG



"""
import sys
import platform
import os
try:
    # Load python 2 config parser
    import ConfigParser as configparser
except ImportError:
    import configparser

executable = {
    'Linux': {
        'mag_static': 'xfemag64',
        'mag_dynamic': 'cxfemag64',
        'mag_transient': 'tsfemag64',
        'mech_static': 'mefemag64',
        'therm_static': 'thfemag64'},
    'Windows': {
        'mag_static': 'wfemagw64.exe',
        'mag_dynamic': 'wcfemagw64.exe',
        'mag_transient': 'wtsfemagw64.exe',
        'mech_static': 'wmefemagw64.exe',
        'therm_static': 'wthfemagw64.exe'}
}


def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)


def which(program):
    """returns complete pathname of program or None if not found"""
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def get_executable(stateofproblem='mag_static'):
    """returns platform specific pathname of femag executable"""
    progname = executable[platform.system()][stateofproblem]
    femag = which(progname)
    if not femag:
        raise Exception(f"{progname} not found in {os.environ['PATH']}")
        # TODO: for l in ['libstdc++-6.dll', 'libgcc_s_dw2-1.dll']:
        #    lib = which(l)
        #    if not lib:
        #        raise Exception("{} not found in {}".format(
        #            l, os.environ["PATH"]))

        #    transfer_files.append(lib)
    return femag


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
