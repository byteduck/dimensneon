from setuptools import find_packages, setup

setup(
    name='dimensneon',
    packages=find_packages(include=['dimensneon']),
    version='0.0.1',
    description='A drop-in replacement for scanpy\'s t-SNE simulation functionality',
    author='Aaron Sonin',
    license='MIT',
    install_requires=['scanpy', 'leidenalg'],
)
