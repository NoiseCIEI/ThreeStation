Usage
=====

A starting point is to play with the data in the folder `example`,
and compare results to make sure the code is working properly.
The `param.yml` provides a :doc:`config`, summarizing almost all parameters.
After decompressing `data.tar.gz`, you will find
`all.txt, receiver.txt, source.txt`, which contain station meta data
in the format specified in the :doc:`config`. The folder `I2` contains
raw two-station interferograms. Folders `I3_ell, I3_hyp` are example output
from three-station interferometry, which the user should compare with.

Detailed explanations can be seen in :doc:`api`, so here only a summary of
the workflow is given. The design of this code is relatively straight-forward:
use :doc:`config` and `config.py` for configuration, use `interferometry.py`
and `preprocess.py` for general processing functions, and use `core.py`
to run main functions. In `config.py`, you may customize
as long as the variables and functions keep the same api. After checking the
input/output paths in :doc:`config` and setting `make_doc=False` in `config.py`,
you can simply run

.. code-block:: shell

    python path/to/ThreeStation/threestation/main.py path/to/param.yml

to compare with example outputs.

Note that phase shift is not applied to source-specific interferograms
in `example`. However, the deviation from inter-receiver distance is
saved in SAC header and can be used to de-bais in dispersion measurements.
