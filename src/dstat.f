      subroutine dstat(x, d, n)
C
C     Copyright (C) 2017 Thorsten Pohlert
C
C     This program is free software: you can redistribute it and/or modify
C     it under the terms of the GNU General Public License as published by
C     the Free Software Foundation, either version 3 of the License, or
C     (at your option) any later version.
C
C     This program is distributed in the hope that it will be useful,
C     but WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C
C     You should have received a copy of the GNU General Public License
C     along with this program.  If not, see <http://www.gnu.org/licenses/>.
C
C     DESCRIPTION
C     This subroutine calculates the statistic D of
C     the double Grubbs test.
C
      implicit none
      integer :: n
      double precision :: d
      double precision, dimension(n) :: x

C     local variables
      double precision :: qtot, q2, meantot, mean2
     
C     External functions:
      double precision :: mean, ssqr

C     in case of user interrupt
      call rchkusr()

c      xx = x
C     sort array
      call qsort3(x, 1, n)
c      xx = x
c      m = n
c      call rrsort(xx, m)

C     calculate mean and Q
      meantot = mean(x, n)
      qtot = ssqr(x, meantot, n)

C     calculate mean and Q
      mean2 = mean(x, n-2)
      q2 = ssqr(x, mean2, n-2)

      d = q2 / qtot
      
      end subroutine
