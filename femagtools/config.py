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


def get_executable():
    """returns platform specific pathname of femag executable"""
    if platform.system() == "Windows":
        wfemag = which("wfemagw64.exe")
        if not wfemag:
            raise Exception("wfemagw64 not found in {}".format(
                os.environ["PATH"]))
        # TODO: for l in ['libstdc++-6.dll', 'libgcc_s_dw2-1.dll']:
        #    lib = which(l)
        #    if not lib:
        #        raise Exception("{} not found in {}".format(
        #            l, os.environ["PATH"]))

        #    transfer_files.append(lib)
        return wfemag

    xfemag = which("xfemag64")
    if not xfemag:
        raise Exception("xfemag64 not found in {}".format(
            os.environ["PATH"]))
    return xfemag
        

def get_femag():
    """returns femag command"""
    if sys.platform.startswith('linux'):
        if platform.machine() == 'x86_64':
            return 'xfemag64'
        return 'xfemag'

    if platform.machine() == 'AMD64':
        return 'wfemagw64'

    return 'wfemag'
    

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
