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
      double precision function mean(x, n)
C     This function calculates the mean
      implicit none
      integer, intent(in) :: n
      double precision, dimension(n), intent(in):: x
      
      mean = sum(x) / real(n, kind = 8)
      return
      end function


      double precision function ssqr(x, xmean, n)
C     This function calculates Q, which is the sum of squares
      implicit none
      integer, intent(in) :: n
      double precision, dimension(n), intent(in) :: x
      double precision, intent(in) :: xmean
      
C     internal
      integer :: i
      double precision :: tmp

      tmp = 0.0d0
      do i = 1, n
         tmp = tmp + (x(i) - xmean)**2.0d0
      end do
      ssqr = tmp
      return
      end function
