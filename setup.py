from setuptools import setup

setup(
    name='disruption_py',
    version='0.1.0a0',
    author='Herbert Turner',
    author_email='hmturner@mit.edu',
    packages=['disruption_py', 'test'],
    scripts=['scripts/summarize_database.py', 'scripts/train.py',
             'scripts/validate_shots.py', 'scripts/generate_datasets.py'],
    url='http://pypi.python.org/pypi/PackageName/',
    license='LICENSE',
    description='A package for plasma disruption analysis and prediction',
    long_description=open('README.md').read(),
    install_requires=['attrs',
                      'cftime',
                      'cycler',
                      'dataclasses; python_version < "3.7"',
                      'exceptiongroup',
                      'fonttools',
                      'importlib_resources; python_version < "3.7"',
                      'iniconfig',
                      'JayDeBeApi',
                      'JPype1',
                      'kiwisolver',
                      'matplotlib',
                      'netCDF4',
                      'numpy',
                      'packaging',
                      'pandas',
                      'Pillow',
                      'pluggy',
                      'PyMySQL',
                      'pyparsing',
                      'pytest',
                      'python-dateutil',
                      'pytz',
                      'scipy',
                      'scikit-learn',
                      'six',
                      'tomli',
                      ],
    include_package_data=True,
)
