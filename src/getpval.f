      double precision function getpval(x, crit, n)
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
      implicit none
      integer, intent(in):: n  
      double precision, dimension(n), intent(in):: x
      double precision, intent(in) :: crit
C
C     local
      integer :: i, cnt

      cnt = 0
      do i = 1, n
         if (x(i) > crit) then
            cnt = cnt + 1
         end if
      end do
      getpval = real(cnt, kind=8) / real(n, kind=8)
      return
      end function
