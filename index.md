---
title: pomp, an R package for statistical inference using partially observed Markov processes
layout: pomp
---
**pomp** provides a very general realization of nonlinear partially-observed Markov processes (AKA nonlinear stochastic dynamical systems).
These are a generalization of linear state-space and hidden Markov models to nonlinear, non-Gaussian processes in either discrete or continuous time.
In **pomp**, one can implement a model by specifying its unobserved process and measurement components;
the package uses these functions in algorithms to simulate, analyze, and fit the model to data.
The motivation, structure, and contents of the package are described, with examples, in a *Journal of Statistical Software* paper, [an updated version of which is provided on this site](./vignettes/pompjss.pdf).

Currently, **pomp** provides support for

- basic particle filtering (AKA sequential importance sampling or sequential Monte Carlo),
- trajectory matching,
- the approximate Bayesian sequential Monte Carlo algorithm of Liu&nbsp;&amp;&nbsp;West&nbsp;(2001),
- the particle Markov chain Monte Carlo method of Andrieu et al.&nbsp;(2010),
- approximate Bayesian computation (ABC; Toni et al.&nbsp;2009)
- the iterated filtering method of Ionides,&nbsp;Breto,&nbsp;&amp;&nbsp;King&nbsp;(2006),
- the improved iterated filtering method (Ionides et al. 2015),
- probe-matching methods (e.g., Kendall et al. 1999, Wood et al. 2010)
- the nonlinear forecasting method of Ellner&nbsp;&amp;&nbsp;Kendall,
- the ensemble Kalman filter of Evensen (1994, 2009), and
- the ensemble adjustment Kalman filter of Anderson (2001), and
- power-spectrum-matching methods of Reuman et al. (2006).

Future support for a variety of other algorithms is envisioned.
In particular, **pomp** is a platform upon which general inference algorithms for partially observed Markov processes can be implemented.
We welcome contributions in the form of codes, examples, improvements to the documentation, bug reports, feature requests, and requests for help!

Please let the developers know if you find **pomp** useful, if you publish results obtained using it, [if you come up with improvements, find bugs, or have suggestions or feature requests!](http://github.com/kingaa/pomp/issues)</a>
There is also a [wiki](https://github.com/kingaa/pomp/wiki):
you are invited to contribute snippets of code, interesting results, references to papers, and so on.

The package is provided under the GNU Public License. 
**pomp** is under active development:
new features are being added and old features are being improved. 
Although the developers make efforts to preserve backward compatibility, we cannot absolutely guarantee backward compatibility.
We will be sure to include warnings of changes that break backward compatibility in the NEWS file and the [**pomp** news blog](./blog.html).
To keep abreast of new developments, subscribe to the [**pomp** RSS feed](./pomp.atom).
