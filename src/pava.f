C Output from Public domain Ratfor, version 1.0
      subroutine pava(y,w,kt,n)
      implicit double precision(a-h,o-z)
      logical same
      dimension y(n), w(n), kt(n)
      do23000 i = 1,n 
      kt(i) = i
23000 continue
23001 continue
      if(n.eq.1)then
      return
      endif
23004 continue
      same = .true.
      do23007 i = 2,n 
      if(y(i-1) .gt. y(i))then
      k1 = kt(i)
      k2 = kt(i-1)
      do23011 j = 1,n 
      if(kt(j).eq.k1)then
      kt(j) = k2
      endif
23011 continue
23012 continue
      wnew = w(i-1) + w(i)
      ynew = (w(i-1)*y(i-1)+w(i)*y(i))/wnew
      do23015 j = 1,n 
      if(kt(j).eq.k2)then
      y(j) = ynew
      w(j) = wnew
      endif
23015 continue
23016 continue
      same = .false.
      endif
23007 continue
23008 continue
      if(same)then
      goto 23006
      endif
23005 goto 23004
23006 continue
      return
      end
