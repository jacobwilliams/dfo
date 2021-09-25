# DFO FAQs 

## What is DFO? 
DFO is a Fortran package for solving general nonlinear optimization problems that have the following characteristics:
 * they are relatively small scale (less than 100 variables),
 * their objective function is relatively expensive to compute and derivatives of such functions are not available and
 * cannot be estimated efficiently.

There also may be some noise in the function evaluation procedures. Such optimization problems arise ,for example, in engineering design, where the objective function evaluation is a simulation package treated as a black box.

## What if I have a constraint function that is expensive to evaluate and whose derivatives are not available?
A new version of DFO which handles such problems will be available soon. Meanwhile the user can try to put the "difficult" constraint into the objective function with a penalty term and handle the penalty function management themselves externally.

## What if the objective function evaluation fails for certain values of the variables?  
We refer to such situation as "virtual" constraint; i.e. constraint whose value cannot be computed at all. All the information that is available about such a constraint is whether it is violated or not. The user has an option to return a "fail=true" flag if the objective function cannot be computed at a given point. The user should not assign artificially high values to the objective function at such points because this can seriously affect the efficiency of the algorithm.

## Is it useful to use DFO if the functions are inexpensive to evaluate?
Yes, but the runs will not be very efficient if the problem size is close to 100 rather than close to 20. A version of DFO for problems with inexpensive function evaluations is under development.

## Does it matter if the underlying problem is essentially smooth or not?
Maybe. The theory behind the algorithm assumes that derivatives exist but are not available. In practise the algorithm can be applied when derivatives do not exist.

## Will DFO find me a global minimum?
Maybe, but the algorithm is designed only to find local minima. In any case the user should run it from a variety of starting points. Because of the
general approach the method has a tendency not to converge to "poor" local optima.

## Are there any publications on DFO? 
Yes.
 1. A. R. Conn, K. Scheinberg and Ph.L. Toint, Recent progress in unconstrained nonlinear optimization without derivatives. Mathematical Programming, Vol. 79 (1997), pp. 397 - 414.
 1. A. R. Conn, K. Scheinberg and Ph.L. Toint, On the convergence of derivative-free methods for unconstrained optimization. Approximation Theory and Optimization: Tributes to M. J. D. Powell , Eds. A. Iserles and M. Buhmann, (1997), pp. 83-108, Cambridge University Press.
 1. A. R. Conn, K. Scheinberg and Ph.L. Toint, A derivative free optimization algorithm in practice in Proceedings of 7th AIAA/USAF/NASA/ISSMO Symposium on Multidisciplinary Analysis and Optimization, St. Louis, MO, 1998.

## What platforms does DFO run on? 
DFO has been tested on:
 * AIX V4
 * Linux

## What do I need to install DFO? 
To run DFO the user needs to have the BLAS routines and the LAPACK routines which are not distributed with the open source. In addition DFO uses a derivative based optimization package to solve a trust region minimization subproblem on every iteration. In theory any general purpose gradient based nonlinear optimization package can be used, if an appropriate interface is provided. Currently there is are interfaces for the FORTRAN 77 package NPSOL and for the C package CFSQP. The user can independently obtain either of these packages or develop their own interface with a different package.

## How do I test if I installed DFO correctly? 
There are 7 test problems from the Hock and Schittkowski test set that are provided with DFO. See the README file for instructions on running the tests. 

## Misc

### To use with IPOPT:

Obtain IPOPT and install the IPOPT library `libipopt.a'as described in
the IPOPT instructions.  In order to use IPOPT for DFO, you do not
need the AMPL Solver Interface or CUTEr.

Once IPOPT has been compiled and installed, edit the Makefile here in
the DFO directory so that the IPOPTLIB variable points to the
libipopt.a that you installed from the IPOPT distrubtion.

In addition, you might need to specify (in the variables ADDLIBS) any
libraries that IPOPT requires, such as BLAS and LAPACK, unless you
downloaded the source code for these during the installation of IPOPT
and they have been included into the IPOPT library `libipopt.a.'

Once that is done, run command "make libdfo_ipopt.a" to make
the DFO library.

Additional files that are required are:
dasum.f, dasum.f, dgeqpf.f, dpocon.f, dpotrf.f, dsyrk.f, dtrsm.f 
from LAPACK, http://www.netlib.org/lapack/, 
and ranlux.f from
http://www.camk.edu.pl/~tomek/htmls.refs/ranlux.f.html.
The user can have LAPACK and BLAS routines in a separate library.
In that case the only file needed is  ranlux.f. Otherwise
these files can all be included as objects in 
the APPL variable in the makefile.

If using IPOPT the the BLAS routines are already included, so there
is not need to point to them again.

Once all the libraries are in place, run "make dfotest_ipopt"
to create a test example. File dfotest.f contains examples
of all user provided routines. The test set includes 8 problems
briefly described in dfotest.f. They are Hock and Schittkowski
problems with various choices of treating constraints as "easy" or
"difficult" (see the manual for meaning). To run the first example
type "dfotest_ipopt 1" and press enter, etc.

### To use with NPSOL:

Obtain NPSOL and create a library. Change the makefile to point to that 
library.


Once that is done, run command "make libdfo_npsol.a" to make
the DFO library. Run command "make libappl.a" to make a library
of auxiliary routines. Change the makefile to point to these
libraries (e.g. change line
"DFO_IPOPT  = /u/katyas/Dfo/dfo_constr/libdfo_npsol.a"
to "DFO_IPOPT  = ./libdfo_npsol.a" to point to the current directory).

The files that are needed by libappl.a are:
dasum.f, dasum.f, dgeqpf.f, dpocon.f, dpotrf.f, dsyrk.f, dtrsm.f 
from LAPACK, http://www.netlib.org/lapack/, 
and ranlux.f from
http://www.camk.edu.pl/~tomek/htmls.refs/ranlux.f.html
The user can have LAPACK routines in a separate library and
modify the libappl.a library (or completely delete it, as long
as ranlux.f is linked with thre rest of the code).

The user also needs to point to the local BLAS library.

Once all the libraries are in place, run "make dfotest_npsol"
to create a test example. File dfotest.f contains examples
of all user provided routines. The test set includes 8 problems
briefly described in dfotest.f. They are Hock and Schittkowski
problems with various choices of treating constraints as "easy" or
"difficult" (see the manual for meaning). To run the first example
type "dfotest_npsol 1" and press enter, etc.
