import setuptools

import sys
version = sys.version_info[0:2]
if version < (3, 9):
    raise Exception('Requires 3.11 > python > 3.9, you have python {}.{} installed'.format(
        version[0], version[1]))
elif version == (3, 11):
    raise Exception('Requires 3.11 > python > 3.9, you have python {}.{} installed'.format(
        version[0], version[1]))
else:
    print('Python version 3.9 or later requirement already satisfied')

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setuptools.setup(name='topoExtend',
                 use_scm_version=True,
                 setup_requires=['setuptools_scm'],
                 description='A SimScale API wrapper with easy to set objects for External Building Aerodynamics',
                 url='https://github.com/DHLynch/topo_extend',
                 author='SimScale GmbH - AE Team',
                 author_email='dlynch@simscale.com',
                 license='MIT',
                 install_requires=requirements,
                 packages=setuptools.find_packages(
                     exclude=(
                         '.github', 'assets', 'recipe', 'plugin', 'tests', '*.egg-info'
                     )
                 ),
                 zip_safe=False)