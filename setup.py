import setuptools

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setuptools.setup(name='topo-extend',
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