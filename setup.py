#!/usr/bin/env python
"""Python setup.py file for ratran_python
"""

from setuptools import setup, find_packages

setup(
    name='ratran_python',
    version='0.1',
    author='Magnus Persson',
    author_email='magnusp@vilhelm.nu',
    packages=find_packages(),
    license='BSD',
    description='Wrapper around RATRAN',
    #install_requires=[],
)