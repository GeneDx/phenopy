from setuptools import find_packages, setup

from phenopy import __project__, __version__

from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name=__project__,
    packages=find_packages(),
    version=__version__,
    description='Phenotype comparison scoring by semantic similarity.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Kevin Arvai <karvai@genedx.com>, Kyle Retterer <retterer@genedx.com>, Carlos Borroto <cborroto@genedx.com>, Vlad Gainullin <vgainullin@genedx.com>, Vincent Ustach <vustach@genedx.com>, Stephen McGee <smcgee@genedx.com',
    author_email='<datascience@genedx.com>',
    license='',
    entry_points={
        'console_scripts': [
            f'{__project__} = {__project__}.__main__:main',
        ]
    },
    include_package_data=True,
    install_requires=[
        'fire==0.5.0',
        'gensim<4.0',
        'networkx==2.6.3',
        'numpy==1.21.6',
        'obonet',
        'pandas==1.3.5',
        'scipy==1.7.3',
        'requests==2.31.0',
    ]
)
