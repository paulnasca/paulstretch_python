from setuptools import setup, find_packages
from codecs import open
from os import path


here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README'), encoding='utf-8') as f:
    long_description = f.read()


setup(
    name='paulstretch_python',
    version='0.1.0',

    description='paulstretch algorithm in python',
    long_description=long_description,
    url='https://github.com/paulnasca/paulstretch_python',
    author='Nasca Octavian PAUL',
    author_email='',
    license='Public Domain',
    keywords='paulstretch sound',

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Other Audience',
        'License :: Public Domain',
        'Topic :: Multimedia :: Sound/Audio :: Sound Synthesis',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],

    packages=find_packages(exclude=['contrib', 'docs', 'tests']),
    install_requires=['numpy', 'scipy'],

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    # entry_points={
    #     'console_scripts': [
    #         'sample=sample:main',
    #     ],
    # },

    scripts=[
        'paulstretch_stereo.py',
        'paulstretch_mono.py',
        'paulstretch_newmethod.py'
    ],
)
