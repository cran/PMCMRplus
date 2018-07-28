Installation of PMCMRplus and (external) dependencies
=====================================================

Description
-----------

In order to use the extended functions of the R package **PMCMRplus**
(&gt;= 5.0), several additional R packages available from CRAN need to
be imported, i.e. **mvtnorm** (Genz and Bretz and 2009, Genz et al.
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

Installation under Linux
------------------------

R packages for Unix / Linux are distributed in source form. First check,
whether **PMCMRplus** can be installed or updated by running the
following function from within R:

    # update PMCMRplus
    update.packages("PMCMRplus")

    # or install
    install.packages("PMCMRplus")

Both R packages **Rmpfr** and **gmp** need compilation and are wrapper
functions for the external libraries (i.e. not shipped with R) `libmpfr`
(Fousse et al. 2007, <http://www.mpfr.org/>) and `libgmp`
(<https://gmplib.org/>). For a correct compilation, the corresponding
header files of the external libraries are required. Therefore, it is
possible that the installation process breaks up with an error message
such as:

    ...
    configure: error: GNU MP not found ...
    ...
    configure: error: Header file mpfr.h not found

However, both libraries and their header files are commonly available on
various Linux distributions.

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

References
----------

L. Fousse, G. Hanrot, V. Lefevre, P. Pelissier, P. Zimmermann (2007)
MPFR: A Multiple-precision Binary Floating-point Library with Correct
Rounding. *ACM Trans. Math. Softw. 33*. 13.
<http://doi.acm.org/10.1145/1236463.1236468>.

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
