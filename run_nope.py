# -*- coding: utf-8 -*-

import sys



##
## Command line options
##

if len(sys.argv) == 7:
    ## Mass loss rate by the planet [g/s]
    dotM_p_ref = float(sys.argv[1])
    print("dotM_p_ref = %.0e [g/s]" % dotM_p_ref)
    ## Initial orbit [AU]
    a_ini = float(sys.argv[2])
    print("a_ini = %g [AU]" % a_ini)
    ## Initial time [yr]
    t_ini = float(sys.argv[3])
    print("t_ini = %g [yr]" % t_ini)
    ## Alpha parameter
    alpha = float(sys.argv[4])
    print("alpha = %g" % alpha)
    ## Turbulent viscosity power index
    beta = float(sys.argv[5])
    print("beta = %g" % beta)
else:
    print("Use\n\t%s <dotM_p_ref [g/s]> <a_ini [AU]> <t_ini [yr]> <alpha> <beta>" % sys.argv[0])
    print("Example:\n\t%s 1e12 0.7 1e7 0.01 1.5" % sys.argv[0])
    sys.exit(1)


## Run
import config_nope as config
from main import main

main(config, dotM_p_ref, a_ini, t_ini, alpha, beta)
