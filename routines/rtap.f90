
                            program rtap

  implicit real*8(A-H,O-Z)

  ! Define NAMELIST / input variables
  NAMELIST/COLLINP/Nodens,notemp,sDens,fDens,IC,popt,grad, &
          & worder,petol,mntmp,mxtmp,wrngup,wrngdn,pops,corapp,lte
  NAMELIST/LIRA/uplev,uplevl,lolev,lolevl,uppec,lopec,ppec

  ! Initialise all input variables
  integer :: clo=0,einf=0,cup=0,worder=0,corapp=0,pops=0, &
           & petol=0,lopec=1,uppec=2,ppec=0,mntmp=1,mxtmp=-1, &
           & wrngup=0,wrngdn=0,grad=0,popt=0

  integer :: i,ii,j,jj,jjj,k,kk,kkk,order,Ne,prnt, &
           & IC,icount,dum,dumm, &
           & jcount,levels,notemp,Z,tripi,tripj,tk,Nee, &
           & Nodens,sDens,fDens,uplev,lolev,kcount,dumrr, &
           & uplevl,lolevl,lte,collprint

  real*8 :: hold,IH,dumr,start,finish,start1,finish1, &
          & start2,finish2,start3,finish3,petolr, &
          & const,boltz,factor,plank,sol,Me,zero, &
          & newDen,corappv,ltev,DE1,DE2,totAup,totAlo

  real*8, allocatable :: aij(:,:),rate(:,:,:),temp(:), &
                       & holde(:),ED(:),XJ(:),weight(:), &
                       & en(:),err(:),lratio(:,:), &
                       & COLL(:,:,:,:),LHS(:,:,:), &
                       & pop(:,:,:),wlength(:,:), &
                       & pec(:,:,:,:),wout(:)

  integer, allocatable :: ind(:),mult(:),bEll(:), &
          & ipoint(:), jpoint(:)

  character*13, allocatable :: config(:)
  character*3 :: charin,chr

  parameter(MXEN=999)         ! NUMBER OF TARGET STATES

  ! Open up files for input/output
  open(15, file='rad.inp',status='unknown', form='formatted')
  open(16, file='rad.out',status='unknown', form='formatted')
  open(404, file='const_dens',status='unknown', form='formatted')
  open(405, file='const_temp',status='unknown', form='formatted')

  call cpu_time(start)
  call cpu_time(start1)

  ! Read the input from rad.inp
  read(15,collinp)
  read(15,lira)
  if (mxtmp.eq.-1) then
    mxtmp = notemp
  endif

  if ((Nodens.lt.2).or.(Nodens.gt.999)) then
      print *, 'Warning, Nodens in radinp must be 2 &
               & < Nodens < 999'
  endif

  ! Set up electron density mesh
  allocate(ED(Nodens))
  newDen = real(sDens)

  if (Nodens.eq.2) then
    ED(1) = 10.0d0**(real(sDens))
    ED(2) = 10.0d0**(real(fDens))
  else
    if ((real(fDens) - real(sDens).lt.1.1)) then
      Me = ((real(fDens) - real(sDens)))/ &
         & (real(Nodens) - 1.0)
    else
      Me = ((real(fDens) - real(sDens)))/ &
         & (real(Nodens) - 1.0)
    endif

    ! Initialise the first value
    ED(1) = 10.0d0**(real(sDens))

    ! Continue the mesh for remaining values
    do i=2,Nodens
      newDen = newDen + Me
      ED(i) = 10.0d0**(real(newDen))
      print *, 'Density: ', i, ED(i)
    enddo
  endif

  close(15)

  ! Define constant physical quantities
  const = 2.1716d-8          ! 2pia_0 const. for rates (q_ij)
  IH = 13.605698             ! Ionization of Hydrogen
  boltz = 8.61733247d-5      ! Boltzmann constant in eV
  plank = 4.13566766225d-15  ! Plancks constant to convert to Ang
  sol = 2.99792458d8         ! Speed of light

  ! Open up adf04 file for reading atomic data
  open(10, file='adf04', status='old')
  read(10,1000) Z

  ! Calculate the number of levels in adf04 file
  levels = 0
  do i=1,MXEN
    read(10,*) dum
      if(dum.eq.(-1))then
        go to 200
      else
        continue
      endif
    levels = levels + 1
  enddo

  ! Reset input back to the start once levels are determined
  200 continue
  do i=1,(levels+1)
    backspace(10)
  enddo

  ! Start to allocate necessary arrays
  allocate(ind(levels))
  allocate(en(levels))
  allocate(mult(levels))
  allocate(bEll(levels))
  allocate(XJ(levels))
  allocate(weight(levels))
  allocate(config(levels))
  allocate(aij(levels,levels))
  allocate(rate(levels,levels,notemp))
  allocate(temp(notemp))
  allocate(holde(notemp+einf))
  allocate(err(levels))

  ! Write output to log file
  write(16,7000)
  write(16,7001) Z, notemp, Nodens, sDens, fDens, Me
  write(16,7002)
  write(16,7003) levels

  ! Calculate rates from Upsilons in ADF04
  do i=1,levels
    read(10,2001) ind(i), config(i), mult(i), bEll(i), XJ(i), &
               &  en(i)

    if (i.eq.levels) then
      if (ind(i).eq.levels) then
        print *, 'Read of levels correct, total in ADF04: ', levels
      else
        print *, 'Read of crucial input data INCORRECT.'
      endif
    endif

    ! cm-1 -> eV conversion
    en(i) = en(i)*0.000123984
    write(16,7004), ind(i), config(i), en(i)

    ! Statistical weightings => LSpi / Jpi :
    if (IC.eq.0) then
      weight(i) = (mult(i))*(2*bEll(i)+1)
    else
      weight(i) = (2*XJ(i)+1)
    endif
  enddo

  ! Read temperatures
  read(10,*) dumm
  read(10,*) dumr, (temp(i), i=1,(notemp))
  
  print *, 'Check temperatures are correct:', (temp(i), i=1,(notemp))


  ! tripi/j, hold and holde used for holding variables until complete
  if (grad.eq.1) then
    open(5, file='adf04rad', status='unknown', form='formatted')
  endif

  ! Initialize all arrays
  aij = 0.0
  rate = 0.0

  ! Number of rows to calculate
  order = ((levels**2-levels)/2)

  ! Loop over everything
  do i=1,order
    read(10,2003) tripj, tripi, hold, (holde(j), j=1,notemp+einf)
    if (grad.eq.1) then
      read(5,*) dumrr, dumm, hold
    endif
    if (hold.lt.1e-20) then
      hold = 0.0
    endif

    ! Start to fill in A values (up and down)
    aij(tripj,tripi) = hold
    aij(tripi,tripj) = hold

    ! Excitation rate depends on T only
    do k=1,notemp
      rate(tripi,tripj,k) = (const/weight(tripi))* &
           & sqrt((IH)/(boltz*temp(k)))* &
           & exp((en(tripi)-en(tripj))/ &
           & (boltz*temp(k)))*holde(k)

    ! De-excitation rates
      rate(tripj,tripi,k) = (weight(tripi)/weight(tripj))* &
           & exp(-(en(tripi)-en(tripj))/(boltz*temp(k)))* &
           & rate(tripi,tripj,k)
    enddo
  enddo

  ! Calculate energy differences
  DE1 = en(uplev)-en(uplevl)
  DE2 = en(lolev)-en(lolevl)

  ! Important formats for ADF04 file
  ! These may need changed depending
  ! on the format of the adf04 file.
  1000 format(13X,I2)
  1001 format(2X,I3,3X,A13,4X,I1,1X,I1,1X,F4.1,3X,F9.1)
  1003 format(2X,I2,2X,I2,25F8.5)
  2001 format(2X,I3,3X,A13,4X,I1,1X,I1,1X,F4.1,4X,F8.1)
  2003 format(1X,I3,1X,I3,25F8.5)
  2010 format(1X,I3,1X,I3,1X,16F8.6)

  ! Deallocated arrays that are no longer needed
  if (allocated(mult)) deallocate(mult)
  if (allocated(bEll)) deallocate(bEll)
  if (allocated(XJ)) deallocate(XJ)

  ! Read is complete -
  !   -> A-values stored in aij(i,j)
  !   -> rates stored in rate(i,j, k)

  !        ####################################
  !        # Fill the Collision matrix (COLL) #
  !        ####################################

  allocate(COLL(levels,levels,notemp,Nodens))
  COLL = 0.0d0

  print *, 'Begin filling collision matrix..'
  do Ne=1,Nodens
    edens = edens + Me
    print *, Ne, edens

    do k=mntmp,mxtmp
      do i=1,levels
        do j=1,levels
          if (i.eq.j) then
            do kk=1,j
              COLL(i,j,k,Ne) = COLL(i,j,k,Ne) - aij(j,kk)
            enddo

            do kkk=1,levels
              if (kkk.eq.j) then
                go to 99
              else
                COLL(i,j,k,Ne) = COLL(i,j,k,Ne) - &
                               & ((ED(Ne))*(rate(j,kkk,k)))
              endif
              99 continue
            enddo
          else if (i.gt.j) then
            COLL(i,j,k,Ne) = COLL(i,j,k,Ne) + ((ED(Ne))*(rate(j,i,k)))
          else if (i.lt.j) then
            COLL(i,j,k,Ne) = COLL(i,j,k,Ne) + &
                           &  aij(i,j) + ((ED(Ne))*(rate(j,i,k)))
          endif
        enddo
      enddo
    enddo
  enddo

  call cpu_time(finish1)
  print '("Basic Read = ",f6.3," seconds.")',finish1-start1
  call cpu_time(start2)


  !        ########################
  !        # Gaussian Elimination #
  !        ########################
  !

  allocate(LHS((levels),notemp,Nodens))
  LHS = 0.0d0

  ! Immediately create the Left hand side of the eq.
  do k=mntmp,mxtmp
    do kk=1,Nodens
      do kkk=1,(levels)
        LHS(kkk,k,kk) = -COLL((kkk),1,k,kk)
      enddo
    enddo
  enddo

  print *, 'Begin Guassian elimination..'
  print *, 'Density      Temperature'
  print *, '________________________'

  ! Loop over densities, then temperatures
  ! Guassian elimination, tweaking all the elements of C
  do Ne=1,Nodens
    print *,'  ', Ne
    do kk=mntmp,mxtmp
      print *, '           ', kk
      do i=2,(levels-1)
        do j=i+1,levels
          factor=COLL(j,i,kk,Ne)/COLL(i,i,kk,Ne)
          COLL(j,i,kk,Ne)=0.0d0
          LHS(j,kk,Ne) = LHS(j,kk,Ne) - factor*LHS(i,kk,Ne)

          do k=i+1,levels
            COLL(j,k,kk,Ne) = COLL(j,k,kk,Ne) - factor*COLL(i,k,kk,Ne)
          enddo
        enddo
      enddo
    enddo
  enddo

  call cpu_time(finish2)
  print '("Gaussian elimination = ",f6.3," seconds.")', &
        &  finish2-start2

  write(16,7005) levels
  jcount = 0

  do k=2,levels
    if (mod(k,8).eq.0) then
      write(16,7006) (LHS(kk+(jcount*8),mntmp,1), kk=1,8)
      jcount = jcount + 1
    else if (k.eq.levels) then
      write(16,7006) (LHS(levels-mod(levels,8)+kk,mntmp,1),kk=1,(mod(levels,8)))
      continue  ! Is this needed?
    endif
  enddo

  !        #############################
  !        # Populations & Line ratios #
  !        #############################

  ! Initialization
  allocate(pop(levels,notemp,Nodens))
  allocate(lratio(Nodens,notemp))
  pop = 0.0
  icount = 0

  do Ne=1,Nodens
    do kk=mntmp,mxtmp
      do i=levels,2,-1
        pop(i,kk,Ne) = (LHS(i,kk,Ne))/(COLL(i,i,kk,Ne))
        if (i.ne.levels) then
          do j=levels,(i+1),-1
            pop(i,kk,Ne) = pop(i,kk,Ne) - &
                        &  ((COLL(i,j,kk,Ne)*pop(j,kk,Ne))/ &
                        &  (COLL(i,i,kk,Ne)))
          enddo
        endif
      enddo
    enddo
  enddo

  ! Plot the populations if pops = 1
  if (pops.eq.1) then
    do i=1,levels
      write(charin,'(I3)') i
      charin = adjustl(charin)
      open (410 + i, file='population_'//trim(charin), &
         & status='unknown', form='formatted')
      write(410 + i,7011) (temp(j), j=mntmp,notemp)
      do Ne=1,Nodens
        write(410 + i,7012) ED(Ne), (pop(i,k,Ne), k=mntmp,mxtmp)
      enddo
      close(410+i)
    enddo
  endif


  ! Keeping Density constant...
  write(16,7007) uplev, uplevl, lolev, lolevl
  prnt = 0

  do Ne=1,Nodens
    if (prnt.eq.1) then
      write(16,7008) Ne, ED(Ne)
      do kk=mntmp,mxtmp
        lratio(Ne,kk) = (pop(uplev,kk,Ne)*aij(uplev,uplevl))/ &
                     &  (pop(lolev,kk,Ne)*aij(lolev,lolevl))
        write(16,7015) temp(kk), lratio(Ne,kk)*(DE1/DE2)
      enddo
    elseif (prnt.eq.0) then
      if (Ne.eq.1) then
        write(16,*) '-> Increasing density'
        write(16,7013) (ED(jj), jj=1,Nodens)
      else
        continue
      endif

      ! Loop over temps
      do kk=mntmp,mxtmp
        lratio(Ne,kk) = (pop(uplev,kk,Ne)*aij(uplev,uplevl))/ &
                      & (pop(lolev,kk,Ne)*aij(lolev,lolevl))
      enddo
!
      if (Ne.eq.Nodens) then
        do kk=mntmp,mxtmp
          write(16,7012) temp(kk),(lratio(ii,kk)*(DE1/DE2), ii=1,Nodens)
          write(404,*) temp(kk),(lratio(ii,kk)*(DE1/DE2), ii=1,Nodens)
        enddo
      endif
    endif
  enddo

  ! Keeping Temperature constant...
  do i=1,2
    write(16,*)
  enddo

  do kk=mntmp,mxtmp
    if (prnt.eq.1) then
      write(16,7010) kk, temp(kk)
      do Ne=1,Nodens
        lratio(Ne,kk) = (pop(uplev,kk,Ne)*aij(uplev,uplevl))/ &
                     &  (pop(lolev,kk,Ne)*aij(lolev,lolevl))
        write(16,7015) ED(Ne), lratio(Ne,kk)*(DE1/DE2)
      enddo

      if (kk.eq.mxtmp) then
        write(16,*)
        write(16,*) '----END LINE RATIOS----'
        write(16,*)
      endif
    elseif (prnt.eq.0) then
      do Ne=1,Nodens
        lratio(Ne,kk) = (pop(uplev,kk,Ne)*aij(uplev,uplevl))/ &
                     &  (pop(lolev,kk,Ne)*aij(lolev,lolevl))
      enddo
      if (kk.eq.mntmp) then
        write(16,*) '-> Increasing temperature'
        write(16,7011) (temp(jj), jj=mntmp,mxtmp)
      else
        continue
      endif

      if (kk.eq.notemp) then
        do Ne=1,Nodens
          write(16,7012) ED(Ne), (lratio(Ne,ii)*(DE1/DE2),ii=mntmp,mxtmp)
          write(405,*) ED(Ne), (lratio(Ne,ii)*(DE1/DE2),ii=mntmp,mxtmp)
        enddo

        write(16,*)
        write(16,*) '----END LINE RATIOS----'
        write(16,*)
      endif
    endif
  enddo

  !        ########################
  !        #    Coronal + LTE     #
  !        ########################

  if (lte.eq.1) then
    write(16,*) 'LTE values (High Ne): '
    do k=mntmp,mxtmp
      ltev = (weight(uplev)/weight(lolev))* &
           & exp(-(en(uplev)-en(lolev))/(boltz*temp(k)))
      write(16,7014) k, temp(k), ltev*((aij(uplev,uplevl))/ &
                  &  (aij(lolev,lolevl)))*(DE1/DE2)
    enddo
  endif

  if (corapp.eq.1) then
    totAup = 0.0
    totAlo = 0.0

    write(16,*) 'Coronal values (Low Ne):'

    do i=1,(uplev-1)
      totAup = totAup + aij(uplev,i)
    enddo

    do i=1,(lolev-1)
      totAlo = totAlo +  aij(lolev,i)
    enddo

     do k=mntmp,mxtmp
       corappv = (rate(1,uplev,k)/rate(1,lolev,k))*(totAlo/totAup)* &
              &  (aij(uplev,uplevl)/aij(lolev,lolevl))
       write(16,7014) k, temp(k),(corappv*(DE1/DE2))
     enddo
  endif

  call cpu_time(start3)

  !        ########################
  !        # Line Identifications #
  !        ########################

  allocate(wlength(levels,levels-1))
  allocate(pec(levels,levels-1,notemp,Nodens))
  petolr = real(10.0**(real(petol)))

  open(18, file='lineID', status='unknown', form='formatted')
  write(18,*) ' U(j) L(i)   Up-Config       &
            & Lo-Config   Up (eV)  Lo (eV)  wlength (Ang)'

  open(19, file='PEC_plot.out', status='unknown', form='formatted')

  do kk=mntmp,mxtmp
    if ((ppec.eq.0).and.(worder.eq.0)) then
      write(19,*) 'Temperature: ', kk
      write(19,*) '--------------------'
      write(19,*) ' Wavelength     Ind      Ne           PEC'
    endif

    do Ne=1,Nodens
      do j=2,levels
        do i=1,j-1
          if ((kk.eq.mntmp).and.(Ne.eq.1)) then
            wlength(j,i) = (sol*plank)/(abs(en(j)-en(i)))
            wlength(j,i) = wlength(j,i)*1d10
            write(18,7016) j, i, config(j), config(i), &
                         & en(j), en(i), wlength(j,i)
          endif
          pec(j,i,kk,Ne) = pop(j,kk,Ne)*aij(j,i)

          ! Check the following:
          !  1) wavelength ordering is off
          !  2) ppec is for all PEC's.
          !  3) Only PEC's above a certain tolerance
          !  4) Only wavelengths in a certain range
          if (worder.eq.0) then
            if (ppec.eq.0) then
              if ((pec(j,i,kk,Ne)).gt.petolr) then
                if ((wrngup.gt.0).or.(wrngdn.gt.0)) then
                  if ((wlength(j,i).gt.wrngdn).and. &
                   & (wlength(j,i).lt.wrngup)) then
                     write(19,7017) wlength(j,i), Ne, ED(Ne), &
                                  & pec(j,i,kk,Ne)
                  endif
                endif
              endif
            endif
          endif

        enddo
      enddo
    enddo
  enddo

  ! ppec = 1 singles out a unique j -> i transition
  if ((ppec.eq.1).and.(worder.eq.0)) then
    write(19,*) 'ppec = ', ppec
    write(19,*) 'Looking at PEC for transition', &
              &  uppec, '->', lopec

    do kk=mntmp,mxtmp
      write(19,7013) temp(kk)
      do Ne=1,Nodens
        write(19,7015) ED(Ne), pec(uppec,lopec,kk,Ne)
      enddo
    enddo

    write(19,*) 'Now loop with a constant density'
    do Ne=1,Nodens
      write(19,7011) ED(Ne)
      do kk=mntmp,mxtmp
        write(19,7015) temp(kk), pec(uppec,lopec,kk,Ne)
      enddo
    enddo
  endif

  !        ########################
  !        #  Wavelength sorting  #
  !        ########################

  if (worder.eq.1) then
     print *, 'Wavelength ordering requested, &
          &   with worder = ', worder
     write (19,*) 'The following PECs are in wavelength order'
     write (19,*) 'Tolerance is = ', petolr
     write (19,*) '------------------'

     ! Allocate arrays for wavelength sorting
     allocate(wout(order))
     allocate(ipoint(order))
     allocate(jpoint(order))

     CALL waveorder(wlength,levels,wout,ipoint,jpoint)

     do kk=mntmp,mxtmp
       write(19,*) 'Temperature: ', kk
       do Ne=1,Nodens
         write(19,*) 'Density: ', Ne
         do k=1,order
           ii = ipoint(k)
           jj = jpoint(k)
           if (pec(ii,jj,kk,Ne).gt.petolr) then
             write(19,7018) ii, jj, wout(k), pec(ii,jj,kk,Ne)
           endif
         enddo
       enddo
     enddo

   else
     print *, 'No wavelength ordering requested, &
           &   with worder = ', worder
  endif

  call cpu_time(finish3)
  print '("Rearranging Wavelength Order = ",f6.2," &
       & minutes.")', (finish3-start3)/60.0
  call cpu_time(finish)
  print '("Total Time = ",f6.2," minutes.")', &
       & (finish-start)/60.0

  ! Open up output files for matlab
  open(950, file='lratio', status='unknown', form='formatted')
  open(951, file='dens', status='unknown', form='formatted')
  open(952, file='temp', status='unknown', form='formatted')

  write(950,*) 'Electron Density | Electron Temperature | Line Ratio'
  write(950,*) '----------------------------------------------------'

  ! Arrays to be used in Matlab
  do Ne=1,Nodens
    write(951,*) ED(Ne)
  enddo
  do kk=mntmp,mxtmp
    write(952,*) temp(kk)
  enddo
  do Ne=1,Nodens
    write(950,*) (lratio(Ne,kk),kk=mntmp,mxtmp)
  enddo

  !        ########################
  !        #  Format specifiers   #
  !        ########################

  7000 format(' --- Radiative Transfer Output ---',/,/,' &
   &  - Basic Read of Input - ', /, &
   & 8X, '_____________________________________________')
  7001 format(8X, '| Atomic Number = ', 22X, I3, ' |', / ,&
   & 8X, '| Number of Temperatures = ', 13X, I3, ' |', / ,&
   & 8X, '| Number of Densities = ', 15X, I4, ' |', / ,&
   & 8X, '| Starting at 10^', I1, 'K and finishing at 10^',&
   & I2, 'K',1X, '|', /, 8X,'| with a mesh increment of ',&
   & F16.4, ' |',/,8X, &
   & '_____________________________________________', /)
  7002 format('RESERVED FOR SPECIFIC LINES NAMELIST LINE.')
  7003 format('  - Basic Level Information -  ', /, /,&
   & 'Total number of levels in ADF04 file: ', I3, /, /,&
   & '_____________________________________________', /,&
   & '| IND.|   CONFIGURATION     |  ENERGY (eV)  |', /,&
   & '_____________________________________________')
  7004 format('| ',I3,' | '3X,A13,3X,' | ',F12.4,'  |')
  7005 format('__________________________________________',&
   & /, /, 'The new left hand side values from N_i = 1 - ',&
   & I3,' for a particular temperature and density &
   & (normally 1 and 1)',/)
  7006 format(8D12.4)
  7007 format(/, 'Looking at the specific line ratio of', I3,&
   & '->',I3,'/', I3,'->',I3, 'over temperatures and densities.')
  7008 format(/, '-------------- Density ...', 4X, I3, &
   & '( ',ES11.4,')')
  7009 format(2X,E10.4,4X,F20.16)
  7010 format(/, '-------------- Temperature ...', 4X, I3,&
   & '( ',ES11.4,')')
  7011 format(3X,'Density',1X,99ES11.4)
  7012 format(999ES11.4)
  7013 format(3X,'Temper.',1X,99ES11.4)
  7014 format(10X,I3,1X,ES11.4,2X,ES11.4)
  7015 format(1X,2ES11.4)
  7016 format(1X,I3,2X,I3,3X,A13,1X,'->',1X,A13,1X,F7.4,2X,F7.4,1X,F14.2)
  7017 format(F14.2,2X,I3,3X,ES11.4,2X,ES11.4)
  7018 format(25X,I3,2X,I3,F14.2,2X,ES11.4)

end program
