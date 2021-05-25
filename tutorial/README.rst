Tutorial on Three-station Interferometry for 2021 Lamont-doherty Seismology Student Workshop IX
-----------------------------------------------------------------------------------------------

This ``README`` is only for the tutorial. If you are interested in using the full package, please see the ``README`` in the main folder.

First, create a conda environment using the ``environment.yml`` in the same repo:

.. code-block:: shell

    cd ThreeStation
    conda env create -f environment.yml
    conda activate threestation

Then, uncompress the data:

.. code-block:: shell

    cd tutorial
    tar zxvf data.tar.gz
