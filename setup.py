from setuptools import find_packages, setup

from phenosim import __project__, __version__

setup(
    name=__project__,
    packages=find_packages(),
    version=__version__,
    description='ACMG criteria variant annotator.',
    author='Carlos Borroto',
    author_email='<cborroto@genedx.com>',
    license='',
    entry_points={
        'console_scripts': [
            '%s = acmg_annotator.__main__:main' % __project__,
        ]
    },
    install_requires=[
        'fire',
        'networkx',
        'numpy',
        'obonet',
        'pandas',
    ]
)
