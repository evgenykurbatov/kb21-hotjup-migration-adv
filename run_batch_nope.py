# -*- coding: utf-8 -*-
"""
Run the grid of models in parallel using multiprocessing

Inspired by https://stackoverflow.com/a/5209746/1038377
"""

import sys
import logging
from multiprocessing.pool import Pool
import numpy as np

import batch_params as bp


info = logging.getLogger(__name__).info



def callback(_):
    info("process finished")



if __name__=="__main__":

    ##
    ## Command line options

    if len(sys.argv) == 3:
        ## Max number of process running simultaneously (e.g. num of CPU cores minus one)
        processes = int(sys.argv[1])
        ## Mass loss rate by the planet [g/s]
        dotM_p_ref = float(sys.argv[2])
        print("dotM_p_ref = %.0e [g/s]" % dotM_p_ref)
    else:
        print("Use\n\t%s <processes> <dotM_p_ref [g/s]>" % sys.argv[0])
        sys.exit(1)


    ##
    ## Run

    ## Initialize INFO level logger
    logging.basicConfig(level=logging.INFO,
                        format=("%(relativeCreated)04d %(process)05d %(threadName)-10s "
                                "%(levelname)-5s %(msg)s"))

    import config_nope as config
    from main import main

    pool = Pool(processes=processes)
    for a_ini in bp.a_ini:
        for alpha in bp.alpha:
            for beta in bp.beta:
                pool.apply_async(main, args=(config, dotM_p_ref, a_ini, bp.t_ini, alpha, beta), callback=callback)
    pool.close()
    pool.join()
