from setuptools import setup, find_packages

__version__ = "0.1"

setup(
    name='r2z',
    author='Tobias Huefner',
    author_email='thuefner@health.ucsd.edu',
    description='A python class for converting rdkit objects to z-matrices',
    version=__version__,
    license='MIT',
    platforms=['Linux'],
    zip_safe=False,
    packages=find_packages(),
    include_package_data=True
    )