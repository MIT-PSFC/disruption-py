from setuptools import setup

setup(
    name='disruption_py',
    version='0.1.0-alpha',
    author='Herbert Turner',
    author_email='hmturner@mit.edu',
    packages=['disruption_py', 'disruption_py.test'],
    scripts=['scripts/summarize_database.py'],
    url='http://pypi.python.org/pypi/PackageName/',
    license='LICENSE',
    description='A package for plasma disruption analysis and prediction',
    long_description=open('README.md').read(),
    install_requires=[
        "pytest",
    ],
    include_package_data=True,
)
