Subroutine waveorder(wave,lev,wstore,indi,indj)

  ! Read in wavelength(lev,lev) and sort.
  integer :: i,ii,j,jj,ord,icount,jcount
  real*8, intent(in) :: wave(lev,lev)
  real*8, intent(out) :: wstore((lev*(lev-1))/2)
  integer, intent(out) :: indi((lev*(lev-1))/2), &
                     &    indj((lev*(lev-1))/2)

  ord = (lev*(lev-1))/2

  do i=2,lev
    do j=1,(i-1)
     jcount = 1
     do ii=2,lev
       do jj=1,(ii-1)
         if ((i.eq.ii).and.(j.eq.jj)) then
           continue
         else
           if (wave(i,j).gt.wave(ii,jj)) then
             jcount = jcount + 1
           else
             continue
           endif
         endif
       enddo
     enddo

   ! store indexes for returning from subroutine
   wstore(jcount) = wave(i,j)
   indi(jcount) = i
   indj(jcount) = j
  enddo
enddo

End Subroutine waveorder

