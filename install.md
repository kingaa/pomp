---
title: install pomp
id: install
layout: pomp
---

# Installation instructions

## From Github:

The github version is usually several weeks ahead of the version on CRAN.
You can install it from the Github repository in three different ways.

1. Run the following in an **R** session:  
```
install.packages("pomp",repos="https://kingaa.github.io/")
```  
On Windows and MacOSX systems, this will cause a precompiled version of the latest release of **pomp** to be installed.
1. Use **devtools**:  
```
devtools::install_github("kingaa/pomp")
```  
1. Download [the source tarball](https://github.com/kingaa/pomp/releases/) of the latest release and install it locally as you would any **R** package.

## From CRAN:

Source and binaries for the [CRAN version are available here](http://cran.r-project.org/package=pomp){:target="_blank"}.

## Changes from version 1:

**pomp** version 2 introduced a number of changes that were not compatible with earlier versions.
This means that code that ran under versions &lt;2 may break under up-to-date versions of **pomp**.
An [upgrade guide](https://kingaa.github.io/pomp/vignettes/upgrade_guide.html){:target="_blank"} is available to help you transition your codes to the new version.

## Important note for Windows users

To make use of **pomp**'s facilities for accelerated computation using compiled C code, and to compile the package from source, you will need the ability to compile C code and dynamically link it into an **R** session.
For this reason, you must install the **Rtools** suite, which can be [downloaded from CRAN](http://cran.r-project.org/bin/windows/Rtools){:target="_blank"}.
**Rtools** is needed *both* to compile and install the development version from source *and* to obtain full value from any version of the package.

Note that, after installing **Rtools**, it is necessary to edit the system path, as [described here](https://cran.r-project.org/bin/windows/Rtools/#putting-rtools-on-the-path){:target="_blank"}, so that these tools can be found by the operating system.

A video tutorial on installing **Rtools** is [available here](https://youtu.be/lmIhiT_QsPE){:target="_blank"}.


## Important note for Mac OS X users

To make use of the package facilities for accelerated computation using compiled C code, you will need the ability to compile C code and dynamically link it into an **R** session.
These facilities are provided in the <code>Xcode</code> app, which is free and can be installed via the App Store or downloaded from [developer.apple.com](https://developer.apple.com/xcode/downloads/){:target="_blank"}.

Some users report problems installing **pomp** from source due to lack of an appropriate **gfortran** installation, which is not included by default in all versions of **Xcode**.
If you have this problem, see [these instructions](https://mac.r-project.org/tools/){:target="_blank"}.
