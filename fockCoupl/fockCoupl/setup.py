from setuptools import setup

setup(name='fockCoupl',packages=['fockCoupl'],entry_points={'console_scripts':['fockCoupl=fockCoupl.fockCoupl:main']},install_requires=['numpy','matplotlib','scipy'])