---
title: panelPomp package alpha release
layout: pomp
---

This alpha release provides tools for working with panel data using partially observed Markov processes.
In particular, this package allows one to model multiple, independent units (or individuals) for each of which one has (potentially multivariate) time series data.
The basic idea driving **panelPomp** is to apply to a collection of units some of the **pomp** package facilities for implementing POMP models, simulating them, and fitting them to time series data.
Regarding fitting, as of this release, only the iterated filtering (`mif2`) algorithm has currently been extended to the **panelPomp** panel setting.
The package is authored and maintained by Carles Breto.

See [the **panelPomp** github site](http://github.com/cbreto/panelPomp) for information, source code, and downloads.
