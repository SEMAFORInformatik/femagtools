try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
from os import path

here = path.abspath(path.dirname(__file__))
# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    description='Femag Tools: a Python API for FEMAG',
    long_description=long_description,
    author='Ronald Tanner',
    url='https://github.com/SEMAFORInformatik/femagtools',
    author_email='tar@semafor.ch',
    version='0.0.5',
    install_requires=['numpy', 'scipy', 'mako'],
    packages=['femagtools', 'femagtools/moo'],
    package_data={'femagtools': ['templates/*.mako']},
    license='BSD',
    name='femagtools',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.4',
        'Topic :: Scientific/Engineering']
)

