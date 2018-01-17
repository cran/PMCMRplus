#include <R.h>
#include <Rmath.h>
/*
    Copyright (C) 2017 Thorsten Pohlert

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

// This is to get the state of RNG
void F77_SUB(rngstart)(void) { GetRNGstate(); }

// This is to put the state of RNG
void F77_SUB(rngend)(void) { PutRNGstate(); }

// This is the random generator for standard normal variates 
double F77_SUB(normrand)(void) { return norm_rand(); }

// This is a call to R_rsort(double *x, int n)
// This call gives strange results, also it 
// does not create an access memory error.
/*void F77_SUB(rrsort)(double *x, int n)
{ 
   R_rsort(x, n);
}*/
