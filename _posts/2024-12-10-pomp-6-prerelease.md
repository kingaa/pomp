---
layout: pomp
title: pomp version 6
---

The [current release](https://github.com/kingaa/pomp/releases/latest) contains two breaking changes over **pomp** versions 5.11.X.X.

- With version 5.8.4, the manner in which one supplies extra model elements (the so-called [`userdata`](https://kingaa.github.io/manuals/pomp/html/userdata.html), elements other than parameters, state variables, observations, and covariates that are available to the basic model components) was changed.
The old way of doing things remained available, but using it generated a warning message, with the promise that the warning would someday become an error.
That day has come.
As of version 6, to supply userdata to the basic model components, one must explicitly use [the `userdata` argument](https://kingaa.github.io/manuals/pomp/html/pomp.html#:~:text=This allows the user to pass information to the basic components outside of the usual routes of covariates) in **pomp** constructors, elementary algorithms, or inference algorithms.
Users who have already adjusted their codes to eliminate the aforementioned warning should experience no change in **pomp** behavior.

- While it has always been good practice in **pomp** function calls to pass arguments by name (as opposed to by position), with version 6.0.1, this becomes (almost) mandatory.
This is meant to prevent a class of difficult-to-trace errors in which an unnamed argument is passed to a lower-level function and potentially ignored, silently.
The chief exception to this rule is the first argument in many **pomp** functions, which is typically a *pomp* object or a data frame.
Users who are in the habit of passing named arguments will not notice any change in **pomp** behavior.
