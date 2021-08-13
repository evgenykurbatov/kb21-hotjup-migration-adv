
# The role of atmospheric outflows in the migration of hot Jupiters


## Abstract

Many of observed hot Jupiters are subject to atmospheric outflows. Numerical simulations have shown that the matter escaping from the atmosphere can accumulate outside the orbit of the planet, forming a torus. In a few 10^8 yr, the mass of the torus can become large enough to exert a significant gravitational effect on the planet. Accumulation of mass, in its own turn, is hindered by the activity of the star, which leads to the photoevaporation of the torus matter. We explore the role of these and other factors in the planet's migration in the epoch when the protoplanetary disk has already disappeared. Using HD 209458 system as an example, we show that the gravitational interaction with the torus leads to the possibility of migration of the planet to its observable position, starting from an orbit ~0.3 AU.

This code is suitable for calculating the evolution of the surface density of an accretion disk under the action of tidal force from the planet, viscous force, and stellar ionizing radiation (by photoevaporation). It is based on the Pringle's model of geometrically thin viscous accretion disk. The planet orbit migrates due to the angular momentum exchange.


## Requirements

Python with NumPy, SciPy, H5Py, Matplotlib. Latest [Anaconda](https://www.anaconda.com/) is OK.


## Code

The core of the model are **model.py** and **pringle.py** scripts. The former contains the _Model_ class which is set-up and instantiated inside of _main_ function of the **main.py** module. The parameters of a _main_ function define the actual model. It is runned by the **run\_\*.py** scripts, while some parameters are given in the **config\_\*.py** files as well as in **aux.py**.

There are two sets of scripts for configuration and running: **\*\_nope.py** and **\*\_xray.py**. The first stands for the models with no photoevaporation, and the second is for the models with photoevaporation by X-ray stellar radiation.


## Publications

E. P. Kurbatov, D. V. Bisikalo. _The role of atmospheric outflows in the migration of hot Jupiters_, [2021, MNRAS, 506, 3128–3137](https://academic.oup.com/mnras/article-abstract/506/3/3128/6321842) | [ADS](https://ui.adsabs.harvard.edu/abs/2021MNRAS.506.3128K) | [arXiv](https://arxiv.org/abs/2101.04112).


## Authors

### Author of the code

**Evgeny P. Kurbatov** _Institute of Astronomy, Russian Academy of Sciences / Moscow, Russia_ [(ORCID iD)](https://orcid.org/0000-0002-1024-9446)
- <kurbatov@inasan.ru>
- <evgeny.p.kurbatov@gmail.com>

### Co-author of the Paper

**Dmitry V. Bisikalo** _Institute of Astronomy of the Russian Academy of Sciences / Moscow, Russia_
- <bisikalo@inasan.ru>
