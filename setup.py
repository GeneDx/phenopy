from setuptools import find_packages, setup

from phenosim import __project__, __version__

setup(
    name=__project__,
    packages=find_packages(),
    version=__version__,
    description='ACMG criteria variant annotator.',
    author='Kevin Arvai <karvai@genedx.com>, Kyle Retterer <retterer@genedx.com>, Carlos Borroto <cborroto@genedx.com>',
    author_email='<datascience@genedx.com>',
    license='',
    entry_points={
        'console_scripts': [
            f'{__project__} = {__project__}.__main__:main',
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
