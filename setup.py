# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 12:54:56 2022

@author: antoine
"""

import setuptools



setuptools.setup(
    name='steroid',                                # should match the package folder
    packages=['steroid'],                     # should match the package folder
    version='0.0.01',                                # important for updates
    license='MIT',                                  # should match your chosen license
    description='Package dedicated to asteroid photometry',
    author='Antoine Choukroun',
    author_email='antcho1@amu.edu.pl',
   
    install_requires=['numpy', 'pandas', 'photutils', 'matplotlib', 'astropy', 'opencv-python', 'glob2', 'keras', 'tqdm', 'scipy'],                  # list all packages that your package uses
    package_data={'steroid': ['best_model_alexnet_good3_99.h5']},
    include_package_data=True,
    
    keywords=["pypi", "astePhot", "tutorial"], #descriptive meta-data
    classifiers=[                                   # https://pypi.org/classifiers
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    
)