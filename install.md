---
title: install pomp
id: install
layout: pomp
---

# Installation instructions

## From Github:

The Github version is usually several weeks ahead of the version on CRAN.
You can install it from the Github repository by executing the following in an **R** session:  
```
install.packages("pomp",repos="https://kingaa.github.io/")
```  
On Windows and MacOS systems, this will cause a precompiled version of the latest release of **pomp** to be installed.

You can also [download the latest release](https://github.com/kingaa/pomp/releases/) and install it locally as you would any **R** package.

## From CRAN:

Source and binaries for the [CRAN version are available here](http://cran.r-project.org/package=pomp){:target="_blank"}.

Install **pomp** from CRAN just like any other **R** package:
```
install.packages("pomp")
```

---------------------------

## Testing your installation

To make use of **pomp**'s facilities for accelerated computation using compiled C code, and to compile the package from source, you will want the ability to compile C code and dynamically link it into an **R** session.
To test this, run the following in an **R** session:
```
source("https://kingaa.github.io/scripts/hello.R",echo=TRUE)
```
This script attempts to compile two simple programs, one in C and one in FORTRAN.
Upon success, you'll see two "Hello!" messages.

If this fails, consult the instructions below, according to your operating system.

--------------------------

## Important note for Windows users

To use **pomp**'s compilation facility, you need to have the **Rtools** suite installed.
This can be [downloaded from CRAN](http://cran.r-project.org/bin/windows/Rtools){:target="_blank"}.
**Rtools** is needed to obtain full value from **pomp**, but also to compile **R** packages from source, if you ever want to do that.

A video tutorial on installing **Rtools** is [available here](https://youtu.be/lmIhiT_QsPE){:target="_blank"}.

------------------------

## Important note for Mac OS users

To use **pomp**'s compilation facility, you need to have the <code>Xcode</code> app installed.
<code>Xcode</code> is free and can be installed via the App Store or downloaded from [developer.apple.com](https://developer.apple.com/xcode/downloads/){:target="_blank"}.

Some users report problems installing **pomp** from source due to lack of an appropriate **gfortran** installation, which is not included by default in all versions of **Xcode**.
If you have this problem, see [these instructions](https://mac.r-project.org/tools/){:target="_blank"}.

------------------------

## Important note for Linux users

You may need to install the `gfortran` compiler on your machine, if it is not already there.
For example, on a Debian-based system (e.g., Ubuntu):
```
sudo apt install gfortran
```
On an RPM-based system:
```
sudo yum install gcc-gfortran
```
