---
title: install pomp
layout: pomp
---

## Installation instructions

### Development version

The development version is usually several weeks ahead of the version on CRAN.
You can install it from the Github repository by executing the following in an **R** session:

```
install.packages("pomp",repos="http://kingaa.github.io")
```

It is also easy to install from source using the <code>devtools</code> package:

```
require(devtools)
install_github("kingaa/pomp")
```

Alternatively, you can download the source as a [tarball](http://github.com/kingaa/pomp/tarball/master)
or as a [zip file](http://github.com/kingaa/pomp/zipball/master).
Unpack the tarball or zipfile (thereby creating a directory named, say, `kingaa-pomp-xxxxxxx`) and build the package source by executing

```
R CMD build kingaa-pomp-xxxxxxx
```

Then install the package from the **R** command line:

```
install.packages("pomp_X.X.X.X.tar.gz")
```

### Stable version

Source and binaries for the [CRAN version are available here](http://cran.r-project.org/package=pomp).


### Important note for Windows users

To make use of **pomp**'s facilities for accelerated computation using compiled C code, and to compile the package from source, you will need the ability to compile C code and dynamically link it into an **R** session.
For this reason, you must install the **Rtools** suite, which can be downloaded from [cran.r-project.org](http://cran.r-project.org/bin/windows/Rtools).
**Rtools** is needed *both* to compile and install the development version from source *and* to obtain full value from any version of the package.

When installing **Rtools**, it is sufficient to choose the "Package authoring installation" option.
Also during the installation, tick the "edit system PATH" box.

Installation of **pomp** from source on Windows has been tested with versions 3.2 and 3.3 of **Rtools** under **R** version 3.2.3.

### Important note for Mac OS X users

To make use of the package facilities for accelerated computation using compiled C code, you will need the ability to compile C code and dynamically link it into an **R** session.
These facilities are provided in the <code>Xcode</code> app, which is free and can be installed via the App Store or downloaded from [developer.apple.com](https://developer.apple.com/xcode/downloads/).

Some users report problems installing **pomp** from source due to lack of an appropriate **gfortran** installation, which is not included by default in all versions of **Xcode**.
If you have this problem, see [these instructions](http://kingaa.github.io/mac-fortran.html).
