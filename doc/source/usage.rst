Usage
=====

A starting point is to reproduce results in `example/`. The `param.yml`
provides a :doc:`config`, summarizing processing parameters. After
decompressing `data.tar.gz`, you will find `all.csv, receiver.csv, source.csv`,
which contain station meta data. The folder `I2/` contains raw two-station
interferograms. Folders `I3_ell/, I3_hyp/` have example output from
three-station interferometry, which should be compared with.

Detailed explanations can be seen in :doc:`api`, so only a summary of the
workflow is given here: use :doc:`config` and `config.py` for configuration,
use `interferometry.py` and `preprocess.py` for general processing functions,
and use `core.py` to run main functions. In `config.py`, you may customize as
long as the variables and functions keep the same api. After checking the
I/O paths in :doc:`config` and setting `make_doc=False` in
`config.py`, you can simply run

.. code-block:: shell

    python path/to/ThreeStation/threestation/main.py path/to/param.yml

to compare with example outputs.

Note that phase shift is not applied to source-specific interferograms in
`example/`. However, the deviation from inter-receiver distance is saved in SAC
header and can be used to de-bais in dispersion measurements.
