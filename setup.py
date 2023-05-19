from setuptools import find_packages, setup

setup(
    name='dimensneon',
    packages=find_packages(include=['dimensneon']),
    version='0.0.1',
    description='A small tool for running t-SNE visualization on 10x scRNA-seq matrices',
    author='Aaron Sonin',
    license='MIT',
    install_requires=['scanpy', 'leidenalg'],
)
