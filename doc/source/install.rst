Installation
============

We recommend using conda_.

First, clone this repo

.. code-block:: shell

    git clone git@github.com:NoiseCIEI/ThreeStation.git

Then create an environment:

.. code-block:: shell

    cd ThreeStation
    conda env create -f environment.yml
    conda activate threestation

And install `pymodule <https://github.com/shane-d-zhang/pymodule>`__
and other dependencies (`setup.py`) in the new environment
using either conda_ or pip_.

After installing dependencies, or if you decide to use pip_
to install dependencies, run the following will install `threestation`
and missing dependencies:

.. code-block:: shell

    pip install -e .


.. _conda: https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html#anaconda-or-miniconda
.. _pip: https://pip.pypa.io/en/stable/
