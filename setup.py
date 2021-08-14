from setuptools import setup, find_packages


setup(
    name='ThreeStation',
    version='0.0.1',
    packages=find_packages(),
    python_requires='>=3.6',
    install_requires=[
        'bottleneck',
        'matplotlib',
        'numpy',
        'obspy',
        'pandas',
        'pyyaml',
        'setuptools',
        'tqdm',

        'pymodule',
    ],
    author='Shane Zhang',
    author_email='shzh3924@colorado.edu',
    description='Three-station interferometry',
    license='MIT',
    url='https://github.com/shane-d-zhang/ThreeStation',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
    )
