
# The role of atmospheric outflows in the migration of hot Jupiters


## Abstract

Many of observed hot Jupiters are subject to atmospheric outflows. Numerical simulations have shown that the matter escaping from the atmosphere can accumulate outside the orbit of the planet, forming a torus. In a few 10^8 yr, the mass of the torus can become large enough to exert a significant gravitational effect on the planet. Accumulation of mass, in its own turn, is hindered by the activity of the star, which leads to the photoevaporation of the torus matter. We explore the role of these and other factors in the planet's migration in the epoch when the protoplanetary disk has already disappeared. Using HD 209458 system as an example, we show that the gravitational interaction with the torus leads to the possibility of migration of the planet to its observable position, starting from an orbit ~0.3 AU.

This code is suitable for calculating the evolution of surface density in a disk under the influence of tidal forces from the planet and viscous forces. It is based on the Pringle's model of geometrically thin viscous accretion disk. The model includes photoevaporation.


## Requirements

Python with Numpy, Scipy, H5Py, Matplotlib. Latest [Anaconda](https://www.anaconda.com/) is OK.


## Code

The core model are in the **model.py** and **pringle.py** scripts. The former contains the _Model_ class which is set-up and instantiated inside of _main_ function of the **main.py** module. The parameters of a _main_ function define the actual model. It is called by the **run\_\*.py** scripts and some parameters are given in the **config\_\*.py** files.


## Publications

E. P. Kurbatov, D. V. Bisikalo. _The role of atmospheric outflows in the migration of hot Jupiters_. Submitted to MNRAS.


## Authors

### Author of the code

**Evgeny P. Kurbatov** _Institute of Astronomy, Russian Academy of Sciences / Moscow, Russia_
- <kurbatov@inasan.ru>
- <evgeny.p.kurbatov@gmail.com>

### Co-author of the Paper

**Dmitry V. Bisikalo** _Institute of Astronomy of the Russian Academy of Sciences / Moscow, Russia_
- <bisikalo@inasan.ru>
