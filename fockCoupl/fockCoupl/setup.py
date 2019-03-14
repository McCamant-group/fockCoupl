from setuptools import setup

setup(name='fockCoupl_DA',packages=['fockCoupl_DA'],entry_points={'console_scripts':['fockCoupl_DA=fockCoupl_DA.fockCoupl_DA:main']},install_requires=['numpy','matplotlib','scipy'])
