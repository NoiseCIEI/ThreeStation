#!/usr/bin/env python
import os
import sys
import logging.config
import warnings



if __name__ == '__main__':
    cwd = os.getcwd()
    if len(sys.argv) != 2:
        sys.exit(f'Usage: python {sys.argv[0]} dirname')
    else:
        os.chdir(sys.argv[1])

    from threestation.core import (
        main,
        PARAM,
    )
    import pymodule as my


    logging.basicConfig(level=PARAM['misc']['log_level'])
    if not sys.warnoptions:
        warnings.filterwarnings(action='ignore')

    my.timing.timing()

    main()

    os.chdir(cwd)
