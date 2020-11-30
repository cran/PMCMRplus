Installation of PMCMRplus and (external) dependencies
=====================================================

Description
-----------

In order to use the extended functions of the R package **PMCMRplus**, 
several additional R packages available from CRAN need to
be imported, i.e. **mvtnorm** (Genz and Bretz 2009, Genz et al.
2015), **multcompView** (Graves et al. 2015), **Rmpfr** (Maechler 2016)
and **gmp** (Lucas et al. 2017). This will be done automatically by R's
package management system.

However, Linux user may encounter some installation problems, as several
R packages require external libraries on the system. This is why this
README file briefly describes the installation procedure of
**PMCMRplus**.

Installation under Windows
--------------------------

As R packages for Windows are distributed in binary form, there should
not be any problem with the installation. Simply run from within R the
following function:

    install.packages("PMCMRplus")

R will automatically install all the relevant dependencies. Provided
that **PMCMRplus** is already installed on your system, simply update
the package or all installed packages with:

    # update PMCMRplus
    update.packages("PMCMRplus")

    # or update all
    update.packages()

Installation under Linux from source packages
---------------------------------------------

R packages for Unix / Linux are distributed in source form. Installation 
of R add-on packages do not require root proviliges and the installation 
directory is set in the variable `$R_LIBS_USER`. The installation
directory is in the users `$HOME` directory. 

First check, whether **PMCMRplus** can be installed or updated by running the
following function from within R:

    # update PMCMRplus
    update.packages("PMCMRplus")

    # or install
    install.packages("PMCMRplus")

Both R packages **Rmpfr** and **gmp** need compilation and are wrapper
functions for the external libraries (i.e. not shipped with R) `libmpfr`
(Fousse et al. 2007, <https://www.mpfr.org/>) and `libgmp`
(<https://gmplib.org/>). For a correct compilation, the corresponding
header files of the external libraries are required. Therefore, it is
possible that the installation process breaks up with an error message
such as:

    ...
    configure: error: GNU MP not found ...
    ...
    configure: error: Header file mpfr.h not found

However, both libraries and their header files are commonly available 
on various Linux distributions.

Ubuntu and Debian
-----------------

Check for the header files by running the following commands outside of
R from the console.

    dpkg -p libgmp-dev
    dpkg -p libmpfr-dev

If any (or both) of the above packages are missing, simply install the
missing package(s) from the repository of your Linux distribution:

    sudo apt-get install libgmp-dev
    sudo apt-get install libmpfr-dev

After successful installation of the above Linux packages, repeat with
the installation of the R package **PMCMRplus** from within R:

    install.packages("PMCMRplus")

Fedora, Redhat, CentOS, opensuse, etc.
--------------------------------------

Check for the header files by running the following commands outside of
R from the console.

    dnf info gmp-devel
    dnf info mpfr-devel

If any (or both) of the above packages are missing, simply install the
missing package(s) from the repository of your Linux distribution:

    sudo dnf install gmp-devel
    sudo dnf install mpfr-devel

After successful installation of the above Linux packages, repeat with
the installation of the R package **PMCMRplus** from within R:

    install.packages("PMCMRplus")

Installation under Linux using binary packages
----------------------------------------------

Ubuntu
------

Installation instructions for R core using an Ubuntu distribution 
can be found here:

<https://cran.r-project.org/bin/linux/ubuntu/>

Additional CRAN binary packages (>1,000) for Ubuntu are availabe 
at the CRAN2deb4ubuntu PPA that can be found here

<https://launchpad.net/~marutter/+archive/ubuntu/c2d4u>

or

<https://launchpad.net/~marutter/+archive/ubuntu/c2d4u3.5>.

Provided, that the above PPA was successfully added to the 
package source list and the user has root (or su, sudo) priviliges, 
one can try to install precompiled `r-cran*` deb packages 
outside of the R environment as

    sudo apt-get install r-cran-pmcmrplus

This will install depending dep packages for **PMCMCRplus**, too.

References
----------

L. Fousse, G. Hanrot, V. Lefevre, P. Pelissier, P. Zimmermann (2007)
MPFR: A Multiple-precision Binary Floating-point Library with Correct
Rounding. *ACM Trans. Math. Softw. 33*. 13.
<https://doi.acm.org/10.1145/1236463.1236468>.

A. Genz, F. Bretz (2009) *Computation of Multivariate Normal and t
Probabilities*. Lecture Notes in Statistics. Heidelberg: Springer.

A. Genz, F. Bretz, T. Miwa, X. Mi, F. Leisch, F. Scheipl, T. Hothorn
(2017) **mvtnorm**: Multivariate Normal and t Distributions. R package
version 1.0-6, <https://CRAN.R-project.org/package=mvtnorm>.

S. Graves, H.-P. Piepho, L. Selzer, S. Dorai-Raj (2015)
**multcompView**: Visualizations of Paired Comparisons. R package
version 0.1-7, <https://CRAN.R-project.org/package=multcompView>.

A. Lucas, I. Scholz, R. Boehme, S. Jasson, M. Maechler (2017) **gmp**:
Multiple Precision Arithmetic. R package version 0.5-13.1,
<https://CRAN.R-project.org/package=gmp>.

M. Maechler (2016) **Rmpfr**: R MPFR - Multiple Precision Floating-Point
Reliable. R package version 0.6-1.
<https://CRAN.R-project.org/package=Rmpfr>.
