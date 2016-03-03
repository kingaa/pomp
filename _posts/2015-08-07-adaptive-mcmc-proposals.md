---
date: 7 Aug 2015
title: Adaptive PMCMC proposals
layout: pomp
---

As of version 0.78-1, **pomp** includes facilities for adaptive proposals for particle MCMC (`pmcmc`) and approximate Bayesian computation (`abc`).
The new function `mvn.rw.adaptive` generates an multivariate normal random-walk MCMC proposal function that adapts in scale and shape.
Thanks to Sebastian Funk for contributing a patch that spurred this development.
This functionality should be regarded as experimental and subject to change.
Please tell me about your experiences with it and contribute suggestions for improvement (or pull requests!) if you can.

