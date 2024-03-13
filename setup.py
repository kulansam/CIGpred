#!/usr/bin/env python

"""Description
Setup script for CIG-Pred: Computational Tool for Cell Identity Gene Discovery
@author:  Kulandaisamy Arulsamy
@contact: KulandaiSamy.Arulsamy@childrens.harvard.edu
"""

import sys
from setuptools import setup, find_packages


def main():
    if float(sys.version[:3]) >= 3.9:
        #sys.stderr.write("We suggested users to use python version >=3.9!\n")
        sys.exit(1)

    setup(name="CIG-Pred",
          version="1.0",
          description="CIG-Pred: Computational Tool for Cell Identity Gene Discovery",
          author='Kulandaisamy Arulsamy',
          author_email='KulandaiSamy.Arulsamy@childrens.harvard.edu',
          packages=find_packages(),
          url='https://github.com/kulansam/CIGpred',
          scripts=['src/CIG_pred',
                   ],
          include_package_data=True,
          package_data={
              '': ['data/*.pkl', 'test/*.txt'],
          },
          license='MIT',
          install_requires=[
              'numpy','qnorm','regex','sklearn','rpy2','bioinfokit','anndata','scanpy','argparse','scipy'],
          data_files=[('', ['data/*']),]
          )

if __name__ == '__main__':
    main()