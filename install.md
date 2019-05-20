---
title: install pomp
id: install
layout: pomp
---

## Installation instructions

### **pomp** version 2 installation instructions

Because it is not fully backward-compatible, **pomp** version 2 is being made available temporarily as the **pomp2** package.
**pomp2** is not available on CRAN but it can be installed directly from the **pomp** package repository.

**pomp2** will be phased out, replaced by version 2 of **pomp**, around the middle of 2019.
At that point, **pomp** versions <2 will no longer be supported.

***It is recommended that users begin making the transition to version 2 now***.
An [upgrade guide](https://kingaa.github.io/pomp/vignettes/upgrade_guide.html){:target="_blank"} is available to help you transition your codes to the new version.
In the interim period, you can have both **pomp** and **pomp2** installed on your system simultaneously.

There are three ways to install **pomp2**.

1. Run the following in an **R** session:  
```
install.packages("pomp2",repos="https://kingaa.github.io/")
```  
1. Use **devtools**:  
```
devtools::install_github("kingaa/pomp")
```  
Note that this installs **pomp2**.
1. Download [the source tarball](https://github.com/kingaa/pomp/releases/) and install locally.

### **pomp** versions <2, from CRAN:

Source and binaries for the [CRAN version are available here](http://cran.r-project.org/package=pomp){:target="_blank"}.


### Important note for Windows users

To make use of **pomp**'s facilities for accelerated computation using compiled C code, and to compile the package from source, you will need the ability to compile C code and dynamically link it into an **R** session.
For this reason, you must install the **Rtools** suite, which can be downloaded from [cran.r-project.org](http://cran.r-project.org/bin/windows/Rtools){:target="_blank"}.
**Rtools** is needed *both* to compile and install the development version from source *and* to obtain full value from any version of the package.

When installing **Rtools**, it is sufficient to choose the "Package authoring installation" option.
Also during the installation, tick the "edit system PATH" box.

### Important note for Mac OS X users

To make use of the package facilities for accelerated computation using compiled C code, you will need the ability to compile C code and dynamically link it into an **R** session.
These facilities are provided in the <code>Xcode</code> app, which is free and can be installed via the App Store or downloaded from [developer.apple.com](https://developer.apple.com/xcode/downloads/){:target="_blank"}.

Some users report problems installing **pomp** from source due to lack of an appropriate **gfortran** installation, which is not included by default in all versions of **Xcode**.
If you have this problem, see [these instructions](http://kingaa.github.io/mac-fortran.html){:target="_blank"}.
