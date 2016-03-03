---
date: 25 Jul 2015
title: Better caching
layout: pomp
---

When a `pomp` object is constructed using `Csnippet`s, the C source is compiled and dynamically loaded into the running **R** session.
Heretofore, this meant that one had to rebuild the `pomp` in each new **R** session.
Changes in [version 0.75-1](http://github.com/kingaa/pomp/releases/tag/0.75-1) now make it so that one can store and re-use `pomp` objects across **R** sessions.
<!--more-->
This is achieved by storing the source code internally to the `pomp` object.
When a `pomp` object is used, a test is first performed to see if the needed shared-object file exists.
If it does not, the source code is written to a file in the session's temporary directory, which is then compiled and loaded.
This feature allows `pomp` objects to be stored and reused across **R** sessions.
To avoid collisions, the name of the file is constructed using a hash of its contents.
