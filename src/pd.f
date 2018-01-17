      subroutine pd(stat, n, m, p)
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
C     This subroutine calculates the pvalue of the 
C     double Grubbs test with a Monte-Carlo simulation
C
      implicit none
      double precision, intent(in) :: stat
      integer, intent(in) :: n, m
      double precision, intent(out) :: p
C     local variables
      integer :: j, i
      double precision, dimension(n) :: x
      double precision, dimension(m) :: statm
      double precision :: tmp

C     External functions:
      double precision :: normrand, getpval
      external :: normrand, getpval

C     get RND state
      call rngstart()

C     Monte Carlo loop
      do j = 1, m

C     user interrupt
         call rchkusr()

C     Get normal random deviates
         do i = 1, n
            x(i) = normrand()
         end do
         
C     Get statistic D for null
         call dstat(x, tmp, n)
         statm(j) = tmp  
      end do

      call rngend()

      p = getpval(statm, stat, m)
      
      end subroutine
