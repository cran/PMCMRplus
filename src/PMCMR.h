#include <R.h>
#include <Rinternals.h>
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
    along with this program.  If not, see <http://www.gnu.org/licenses/>.*/

void F77_NAME(dstat)(double *x, double *d, int *n);
		    
void F77_NAME(pd)(double *q, int *n, int*m, double *p);

void F77_NAME(pava)(double *y, double *w, int *kt, int *n);
							
/* void F77_NAME(muma)(double *x, int *nr, int *k, int *p, int *n, double *t, double *et, double *vt); */
