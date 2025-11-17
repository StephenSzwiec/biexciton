 
        program readwfnew
	implicit none

	include "omp_lib.h"

	integer ::  nthreads, tid
! --------------------------------------------------------------------
! Reads phi(k)s from WAVECAR
! --------------------------------------------------------------------

!>>>>>>>>>>>>>>>>>>>>>>   Declare variables  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	integer     npw,nkpt, nband, ae, nLU, be, n, l, qq,nijp, jae
	integer     nbandmin, nbandmax, ndiff, nHO,ndiffU, ndiffO,ic

	real(8)     emax, kpt(3), eval
	real(8)     a1, b1, c1, A(3,3)       ! Unit-vectors of the cell

	character*5 code
        integer     POTIM, Nexcc,Ndw
!------------------------- 
! New format of WAVECAR:
!--------------------------  
	integer     irec, irecl   
	real(8)     rdum, rispin, rtag

	complex(4), allocatable :: coef(:)		
	real(8)     fweight(20000), eig(20000), gweight(20000)

	integer     i,j,k,ii,kk,iband, ikpt,nxn,ee,hh,ss
	
	real(8)     inkpt, inband, npl
!------------------------------
! For transition dipole moment	
!------------------------------
	integer  n5                     ! dimension of the operator nubla
        real(8)  gx, gy, gz
	integer, allocatable :: XX(:)
	integer, allocatable :: YY(:)
	integer, allocatable :: ZZ(:)
        
        integer, allocatable :: ccc(:,:,:)
        integer, allocatable :: Nexc(:)
        integer, allocatable :: ce(:,:)
        integer, allocatable :: ch(:,:)
       
        integer, allocatable :: cveeh(:,:,:,:)
        integer, allocatable :: cvehe(:,:,:,:)
        integer, allocatable :: cvhhe(:,:,:,:)
        integer, allocatable :: cvheh(:,:,:,:)
        
!        complex(4)  JX(9400,9400)
!        complex(4)  JY(9400,9400)
!        complex(4)  dd1,dd2,dd3
!        complex(4)  JZ(9400,9400)

        complex(8)  q,qc,curx,cury,curz,fx,fy,fz,f
        complex(8)  psia,psib,psic,Rh,Rp,Rshh,Rshp
        complex(8)  psiab,psicab,psicccb,psicb,psicc_b,psic_cb
        
!--------------------------------
! coefficients of wave functions:
!--------------------------------
        complex(8), allocatable :: ac(:,:)
        complex(8), allocatable :: dc(:,:)
        complex(4), allocatable :: psi(:,:)

        complex(4), allocatable :: VCeeh(:)
        complex(4), allocatable :: VCheh(:)
        complex(4), allocatable :: VCehe(:)
        complex(4), allocatable :: VChhe(:)

        CHARACTER(LEN=20) :: t1,t2,t3
        INTEGER :: n1,n2

!---------------------------
! DIMENSIONS OF MATRICES
!---------------------------
!        INTEGER SIZE            ! LENGTH OF coupling matrix
!        parameter(SIZE=500)
!	INTEGER sqSIZE          ! LENGTH OF transition dipole 
!        parameter(sqSIZE=90000)
	
!----------------------------	
! Normalization parameters:
!----------------------------
        INTEGER digit, du, ddo,jj
!        character*20 filename, filename1
        character*20 dumy
!       integer time,t,i,j,nxn
        double precision overlap, PV,R2
        double precision w,gamma_b,gamma_s,denom,num,denom1
        double precision NORM1,NORM2, VV, r,denom2,denomln,denomnl
        double precision repsi, impsi, deltaexc,Eaeki,Ebln,Ejiw,Ebnl
!        double precision D, DD (SIZE,SIZE)  ! Normalized coupling elements
!-----------------------
! energy:
!-----------------------
!       double precision  e(1500), OCUPn(1500), NN(1500)
        double precision, allocatable :: e(:)
        double precision, allocatable :: OCUPn(:)
        double precision, allocatable :: NN(:)
        double precision, allocatable :: Eexc(:)
!, f(50,50,50), fs(50,50,50)
        double precision  dE,dw,hiw,low,RRh,RRp,RRshh,RRshp


!>>>>>>>>>>>>>>>>>>   begining the main program  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!-------------------------------------------	
! Set bands we're interested in
!-------------------------------------------
        open(75,file='input_overlap')
	     read(75,*) nbandmin
             read(75,*) nHO
	     read(75,*) nbandmax
            ! read(75,*) POTIM
	close (75)    
        
        write(*,*) ' check active window size ! '
       
           open(425,file='active_window')
            read(425,*) ddo
            read(425,*) du
            read(425,*) Nexcc
 !           read(425,*) Ndw
 !           read(425,*) low
 !           read(425,*) hiw
            read(425,*) gamma_b
           close(425)

          write(*,*)'ddo=',ddo
          write(*,*)'du=',du

         write(*,*) ' check consistency with "energy_pop" !'
         write(*,*) 'nbandmin, has to be ge',nbandmin
!         read(*,*) nbandmin
!        write(*,*) 'nbandmin=',nbandmin

        write(*,*) 'nbandmax, has to be le',nbandmax
!	read(*,*) nbandmax

        write(*,*)'nHO', nHO

        nLU=nHO+1

        allocate(e(nbandmax-nbandmin+1))
        allocate(NN(nbandmax-nbandmin+1))
        allocate(OCUPn(nbandmax-nbandmin+1))


        write(*,*) 'nbandmin=',nbandmin
        write(*,*) 'nHO=', nHO
        write(*,*) 'nbandmax=',nbandmax
!       write(*,*) 'ndiff=', ndiff,' ndiffO=', ndiffO
!---------------------------------------------------------------------
! Read input "energy_pop" files and input parameters
!---------------------------------------------------------------------	
	write(*,*) ' check file "energy_pop" !'        
	open(11,file="energy_pop")
	do i=nbandmin, nbandmax
           read(11,*) NN(i), e(i), OCUPn(i) 
!          write(*,*) NN(i), e(i), OCUPn(i)        
        enddo

        close(11)
!---------------------------------------	
! OUTPUT FILES:	
!---------------------------------------
         open(18,file="overlap")
!-----------------------------------------
! Read in first WAVECAR: |psi(t)>
!----------------------------------------
        open(12,file="WAVECAR",status="old",form="unformatted", 
     c          access='direct',recl=1500000)

! Read in from WAVECAR:
! the number of k-points
! the number of bands
! the energy maximum
! the cell dimensions
! the type of VASP code

        read(12,rec=1) rdum, rispin, rtag
	     irecl = rdum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  note to people...
!  we have to open close 12 twice:  
!  it's direct access and you need the record length so 
!  that you can read the information, but you don't know the record length
!  until you open the file and read rdum (yes, it changes for every system
!  from vasp.).  thus you have to open it to read rdum so you can use rdum
!   
!  as a second important note- if you are having problems reading the 
!  coefficients, check rtag (use a write statement)...
!
!    rtag = 45210     the coefficients are written in double precision
!    rtag = 45200     the coefficients are written in single precision
!
!   gotta love the VASP people...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	close(12)
        write(*,*)'irecl =',irecl
        write(*,*)'rdum rispin rtag:',rdum, rispin, rtag

        open(12,file="WAVECAR",status="old",form="unformatted", 
     c          access='direct',recl=irecl)

             read(12,rec=1) rdum, rispin, rtag
             read(12,rec=2) inkpt, inband, emax, ((A(i,j),i=1,3),j=1,3)
          write(*,*)'2-nd time rdum, rispin, rtag:',rdum, rispin, rtag
          write(*,*)'inkpt, inband, emax:',inkpt, inband, emax

!----------- for CASE OF RECTANGULAR UNIT CELL a1,b1,c1 are unit-vectors: 
                 a1=A(1,1)
                 b1=A(2,2)
                 c1=A(3,3)

             write(*,*) 'a1=',a1	
             write(*,*) 'b1=',b1
             write(*,*) 'c1=',c1     
             
	     nkpt = inkpt
	     nband = inband

	     write(*,*) 'number k-points nkpt=',nkpt
	     
	!>>>>>>>>>>>>>>>>>> FOR TEST <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	!      write(19,*) 
	!      write(19,*) 'nkpt  =',nkpt
	!      write(19,*) 'nband =',nband
	!      write(19,*) 'emax  =',emax
	!      write(19,*) 'A='
	!      write(19,'(3X,3(1X,f8.3))') (A(i,1),i=1,3)
	!      write(19,'(3X,3(1X,f8.3))') (A(i,2),i=1,3)
	!      write(19,'(3X,3(1X,f8.3))') (A(i,3),i=1,3)
	!>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            irec = 3

        DO 100 ikpt=1, 1 !NOTE RESTRICTED TO GAMMA POINT ONLY

! Read in from WAVECAR:
! the number of pw (plane waves) wave fns?
! the k-point? 

         read(12,rec=irec) npl,kpt(1:3),(eig(j),gweight(j),fweight(j), 
     c                    j=1,
     c                    nband, 1)

	 irec = irec + 1	
         npw = npl


! Brad thinks vasp writes in real values only (i.e. not integer)
!        write(*,*) "length = ", irecl, "  npw = ", npw

       !>>>>>>>>>>>>>>> FOR TEST <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
       !    write(19,*) 'npw= ', npw
       !    write(19,*) 'along x kpt= ', kpt(1)
       !    write(19,*) 'along y kpt= ', kpt(2)
       !    write(19,*) 'along z kpt= ', kpt(3)
       !>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<	

	    allocate(XX(npw))
	    allocate(YY(npw))
	    allocate(ZZ(npw))

!---------------- CALCULATIONS OF an OPERATOR NUBLA(r)-----------------------
!                 FOR FINDING TRANSITION DIPOLE MOMENT

           CALL find_r_dipmom(a1,b1,c1,emax,npw,XX,YY,ZZ,gx,gy,gz)
!----------------------------------------------------------------------------
            allocate(coef(npw))
            allocate(ac(npw,nbandmax-nbandmin+1))
	    ic = 1 
!******************************************************************
            DO 89 iband = 1, nband

!----------------------- Read in from WAVECAR:-----------------------------

                read(12,rec=irec) (coef(i),i=1,npw)
		irec = irec + 1
!---------------------
! Create matrix (ac):
!---------------------

! row: coefficients for plane waves
! column: coeff of 1st plane wave for bands interested in 
! in ac(i,ic), i=1..npw, ic=1..ndiff corresponding to nbandmin..nbandmax
		
                  if(iband.ge.nbandmin.and.iband.le.nbandmax) then
	            do i = 1,npw
	    	        ac(i,ic) = coef(i)
	            enddo
	            ic = ic+1
	        end if

89	     CONTINUE         
!********************************************************************
           deallocate(coef)

100	CONTINUE

	close(unit=12)

        write (*,*)'ac-coefficients were found ic=',ic-1
!       write (*,*)'ac-coefficients'
!                do i = 410,580
!                       write(*,*) 'i=',i
!                       write(*,*) 'a(i)=',a(i,ic-5)
!	            enddo

                write(*,*) 'npw=',npw	

           do i = nbandmin,nbandmax

            overlap = 0.d0

             do k = 1,npw
	       overlap=overlap+conjg(ac(k,i))*ac(k,i)
             end do
           write(18,*) i, overlap
           end do
           close(18)
           write (*,*)'overlap file computed'

!   in ac(i,ic), i=1..npw, ic=nbandmin..nbandmax corresponding to nbandmin..nbandmax
                                             
!  exciton energies written in excE file as 
!  E_k k - exciton state energy, k=1,Nexcc
!  k=1 corresponds to the lowest energy

! exciton wavefunction written in excwf file as 
!  k - exciton state number
! (ddo+1)(du+1) X 4 lines of 
!  ...............
!  j particle level
!  i hole level
!  Re[psi^k_ji]
!  Im[psi^k_ji]
!  ...............
! TRUNCATED exciton wavefunction written in excwf_trunc 
! file where only Nk non-zero entries are recorded as 
!  k - exciton state number
! (Nk) X 4 lines of 
!  ...............
!  j particle level
!  i hole level
!  Re[psi^k_ji]
!  Im[psi^k_ji]
!  ...............
! Nk are recorded in 'excwf_length' file
! as k Nk
!    ....
            allocate(psi(Nexcc,(du+1)*(ddo+1)))
            allocate(Eexc(Nexcc))
            allocate(Nexc(Nexcc))


!            emin=nHO+1
!            hmin=nHO-ddo
            open(112,file="excE")
            open(113,file="excwf_trunc")
            open(1213,file="excwf_length")

            do k = 1,Nexcc! NOTE system specific
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<
             read(112,*) Eexc(k)
             write(*,*)Eexc(k)
 
             read(1213,*) ic,Nexc(ic)
             write(*,*) ic,Nexc(ic)

             if (ic.ne.k) then
              write(*,*)"problem",ic,k
             end if
!<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>
            end do

         ! counter of eh pairs: e,h -> integer
          allocate(ccc(Nexcc,nLU+du+1,nLU+du+1))          
          allocate(ce(Nexcc,(ddo+1)*(du+1)))
          allocate(ch(Nexcc,(ddo+1)*(du+1)))

          do k = 1,Nexcc
           do i = 1,nLU+du+1
            do j = 1,nLU+du+1
             ccc(k,i,j)=0
            end do
           end do
          end do

            do k = 1,Nexcc! NOTE system specific      
      
             overlap = 0.d0       
             read(113,*) ic

             if (ic.ne.k) then
              write(*,*)"problem",ic,k
             end if

             do i = 1,Nexc(ic) 

             read(113,*)ee,hh,repsi,impsi
             ce(ic,i)=ee
             ch(ic,i)=hh
             
             ccc(ic,ee,hh)=i
             psi(ic,i)=CMPLX(repsi,impsi) 
!            write(*,*) ee,hh,psi(ic,ee-(nHO+1)+1,hh-(nHO-ddo)+1)
             overlap=overlap+repsi*repsi+impsi*impsi

	     end do
             write(*,*)ic,' ',overlap,' ',Nexc(ic)
	    end do
            close(112)
            close(113)
            close(1213)

            write(*,*) 'Psis have been read'

           do k = 1,Nexcc! NOTE system specific
!>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            ic=0
            do i = nHO-ddo,nHO
!---------------------------------------------------------
              do j = nLU,nLU+du
!=========================================================
               if (ccc(k,j,i).ne.0) then
                ic=ic+1
               end if
!=========================================================
              end do
!---------------------------------------------------------
            end do
            
            if (ic.ne.Nexc(k)) then
             write(*,*)ic,Nexc(k)
            end if
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
           end do
           write(*,*)'ccc array checked'
!           read(*,*) ic 
!!!         deallocate(ccc)

!NOTE read only VCe1ehe2 since V(e1,h,e,e2)=V(e2,e,h,e1)^*
            open(18,file='VCe1ehe2')
! NOTE NO need to add number of ksij 4plets to the beginning of the VC files
! number of e1ehe2 
!           e1hee2 
!           matrix elements
            nijp=(du+1)*(du+1)*(ddo+1)*(du+1)
            write(*,*)'number of e1ehe2 matrix elements',nijp
            allocate(VCeeh(nijp))

            allocate(cveeh(du+1,du+1,ddo+1,du+1))

            write(*,*)'VCeeh array is  allocated'
             
             ic=0
             do qq=1,nijp
!========================================================================
!                   215         236         164         215 ( 1.329723000744372E-002,-2.224093068271295E-006)       
              read(18,*)j,l,n,k,VCeeh(qq)    
              cveeh(j-nLU+1,l-nLU+1,n-(nHO-ddo)+1,k-nLU+1)=qq
              ic=ic+1
!         write(*,*)q,ck(q),cs(q),ci(q),cj(q),VC2p(q)         
!========================================================================
             enddo

             close(18)

             write(*,*)'V(e1,e,h,e2)  matrix elements have been read'
             write(*,*)'q=',ic,' nijp=',nijp
             write(*,*)j,l,n,k,VCeeh(ic)

             open(19,file='VCh1ehh2')
! NOTE NO need to add number of ksij 4plets to the beginning of the VC files
! number of h1ehh2 
!           h1heh2 
!           matrix elements
            nijp=(ddo+1)*(ddo+1)*(ddo+1)*(du+1)
            write(*,*)'number of h1ehh2 matrix elements',nijp
            allocate(VCheh(nijp))

            allocate(cvheh(ddo+1,du+1,ddo+1,ddo+1))

            write(*,*)'VCheh array allocated'

             ic=0
             do qq=1,nijp
!========================================================================
!                   164         226         164         164 (-0.179791207988772     , 2.617794969849091E-005)      
              read(19,*)j,l,n,k,VCheh(qq)    
          cvheh(j-(nHO-ddo)+1,l-nLU+1,n-(nHO-ddo)+1,k-(nHO-ddo)+1)=qq
              ic=ic+1
!         write(*,*)q,ck(q),cs(q),ci(q),cj(q),VC2p(q)         
!========================================================================
             enddo

             close(19)

             write(*,*)'C V(h1,e,h,h2) matrix elements have been read'
             write(*,*)'q=',ic,' nijp=',nijp
             write(*,*)j,l,n,k,VCheh(ic)
           
            write(*,*)'HO-LU gap is Eg=', e(nLU)-e(nHO)
            write(*,*)'exciton gap is Eg=', Eexc(1)
            write(*,*)'2*gap is 2*Eg=', 2*(e(nLU)-e(nHO))
            write(*,*)'2*exciton gap is', 2*Eexc(1)

!           w=low*Eexc(1)!(e(nLU)-e(nHO)) !eV 
!           write(*,*)'photon energy w=',w
!
!            write(*,*)'highest photon energy in Eg units', hiw
!            dw=(hiw-low)*Eexc(1)/Ndw
!           write(*,*)'dw=',dw
                
        ! gamma_b to be read from active_window
            write(*,*)'width gamma? gamma=',gamma_b
             qq=86
             open(425,file='R_h_spec_nb_w86')
             open(423,file='M2cab_h_spec_w86')
! energy of initial exciton
           w=Eexc(qq)      
            
            write(*,*)'w=',w

            ss=300
!            write(425,*)'{'           

             RRh=0
!             RRp=0
!             RRshh=0
!             RRshp=0

		!The all Sum variable declared above here

!--------------prepare for  OMP thread start

!$OMP PARALLEL PRIVATE(nthreads, tid)


        tid = OMP_GET_THREAD_NUM()
        nthreads = OMP_GET_NUM_THREADS()

	!---not included: arrays, gamma_b, ss, w, qq

!$OMP SECTIONS PRIVATE(ae, be, dE, denom, deltaexc, ii, kk), 
!$OMP& PRIVATE( Rh,Rp,Rshh,Rshp,l,n,psib,Ebln,denomln),
!$OMP& PRIVATE(denomnl,j,i,psic,psicc_b,psicccb,jae, psia),
!$OMP& PRIVATE(Ebnl,Eaeki,denom1,Ejiw,denom2,PV,jj,f,psic_cb)

 !$OMP SECTION

  	do ae = 1, 4  !---------------------------------------------------------
              do be = 1,qq !Does not go past this? !Nexcc! exciton state sum
!---------------------------------------------------------  
! Rexc(c)-->biexc(a,b)=(|R2(a,b)|^2)*delta(Ec-(Ea+Eb))
! Rh to be calculated for each ae,be
               Rh=(0,0)
!               Rp=(0,0)
!               Rshh=(0,0)
!               Rshp=(0,0)

! w=Eexc(c)=Eexc(qq) the initial exciton energy          
! check if the two excitons ae and be satisfy energy conservation w=Ea+Eb
               dE=w-Eexc(ae)-Eexc(be)
               denom=dE*dE + gamma_b*gamma_b
               deltaexc=gamma_b/denom 
               
       if (deltaexc > 1/(gamma_b*ss)) then
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<           
               do ii = 1,Nexc(qq)
!/////////////////////////////////////////////////////////
                do kk = 1,Nexc(be)
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!==================hole spectator===================================
!    read l,n common for the Theta_{j} and Theta_{-j} parts
                l=ce(be,kk)
                n=ch(be,kk)

                psib=psi(be,kk)

!         the Theta_l Theta_{-n} part
             Ebln=e(l)-e(n)-Eexc(be)
             denomln=Ebln*Ebln+gamma_b*gamma_b

!         the Theta_{-l} Theta_{n} part
             Ebnl=e(l)-e(n)+Eexc(be)
             denomnl=Ebnl*Ebnl+gamma_b*gamma_b

!=========================================================
!=========================================================
!         the Theta_j part
!=========================================================
!  indices of exciton gamma or c
             j=ce(qq,ii)
             i=ch(qq,ii)

             psic=psi(qq,ii)

             psicc_b=conjg(psic)*psib
             psicccb=conjg(psic)*conjg(psib)
             
           do k = nLU,nLU+du
!---------------------------------------------------------
               
             jae=ccc(ae,k,i)

             if (jae.ne.0) then
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
             psia=psi(ae,jae)
            
             Eaeki=Eexc(ae)+e(k)-e(i)
             denom1=Eaeki*Eaeki+gamma_b*gamma_b
            
             Ejiw=e(j)-e(i)+w
             denom2=Ejiw*Ejiw+gamma_b*gamma_b

!         the Theta_l Theta_{-n} part            
             PV=Eaeki*Ejiw*Ebln/(denom1*denom2*denomln)

!   find the Theta_l Theta_{-n} VC matrix element
             jj=cveeh(j-nLU+1,l-nLU+1,n-(nHO-ddo)+1,k-nLU+1)
            
      f=(Eexc(ae)-e(k)+e(i))*(e(l)-e(n)-Eexc(be))*(e(j)-e(i)-Eexc(qq))

           Rh=Rh+VCeeh(jj)*PV*f*psia*psicccb!psia*conjg(psic)*conjg(psib)

! the Theta_{-l} Theta_{n} part
             PV=Eaeki*Ejiw*Ebnl/(denom1*denom2*denomnl)

!   find the Theta_{-l} Theta_{n} VC matrix element
!     NOTE use V(e2,h,e,e1)=V(e1,e,h,e2)^*
             jj=cveeh(k-nLU+1,l-nLU+1,n-(nHO-ddo)+1,j-nLU+1)
!            jj=cvehe(j-nLU+1,n-(nHO-ddo)+1,l-nLU+1,k-nLU+1)
            
!      f=(Eexc(ae)-e(k)+e(i))*(e(l)-e(n)-Eexc(be))*(e(j)-e(i)-Eexc(qq))

           Rh=Rh+conjg(VCeeh(jj))*PV*f*psia*psicc_b!psia*conjg(psic)*psib
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            end if
!---------------------------------------------------------	      
	   end do! end k loop
!=========================================================
!=========================================================
!         end of the Theta_j part
!=========================================================

!=========================================================
!         the Theta_-j part
!=========================================================
!  indices of exciton gamma or c
             i=ce(qq,ii)
             j=ch(qq,ii)

             psic=psi(qq,ii)

             psicb=psic*psib
             psic_cb=psic*conjg(psib)

           do k = nHO-ddo,nHO
!---------------------------------------------------------  
             jae=ccc(ae,i,k)

             if (jae.ne.0) then
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
             psia=psi(ae,jae)
                         
             Eaeki=Eexc(ae)-e(i)+e(k)
             denom1=Eaeki*Eaeki+gamma_b*gamma_b
            
             Ejiw=e(i)-e(j)-w
             denom2=Ejiw*Ejiw+gamma_b*gamma_b

!         the Theta_l Theta_{-n} part
            
             PV=Eaeki*Ejiw*Ebln/(denom1*denom2*denomln)

!   find the Theta_l Theta_{-n} VC matrix element
           jj=cvheh(j-(nHO-ddo)+1,l-nLU+1,n-(nHO-ddo)+1,k-(nHO-ddo)+1)
            
             f=(Eexc(ae)-e(i)+e(k))*(e(l)-e(n)-Eexc(be))*(e(i)-e(j)-w)
             
       Rh=Rh+VCheh(jj)*PV*f*conjg(psia)*psic_cb!conjg(psia)*psic*conjg(psib)

! the Theta_{-l} Theta_{n} part
             PV=Eaeki*Ejiw*Ebnl/(denom1*denom2*denomnl)

!   find the Theta_{-l} Theta_{n} VC matrix element
!   NOTE use V(h2,h,e,h1)=V(h1,e,h,h2)^*
         jj=cvheh(k-(nHO-ddo)+1,l-nLU+1,n-(nHO-ddo)+1,j-(nHO-ddo)+1)
!        jj=cvhhe(j-(nHO-ddo)+1,n-(nHO-ddo)+1,l-nLU+1,k-(nHO-ddo)+1)
            
!       f=(Eexc(ae)-e(i)+e(k))*(e(l)-e(n)-Eexc(be))*(e(i)-e(j)-w)

       Rh=Rh+conjg(VCheh(jj))*PV*f*conjg(psia)*psicb!conjg(psia)*psic*psib
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            end if
!---------------------------------------------------------	      
	   end do !end k loop
      
!=========================================================
!=========================================================
!         end of the Theta_-j part
!=========================================================

!///////////////////////////////////////////////////////////	      
	        end do !end of kk loop: b exciton index loop
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\	      
	       end do !end of ii loop: c or qq or gamma exciton index loop 
! if delta(w-Ea-Eb) is big enough there is a contribution to R2
!dimension factor to be added in .nb
         RRh=RRh+conjg(Rh)*Rh*deltaexc
!         write(423,*)'{',qq,',',ae,',',be,',',conjg(Rh)*Rh,'},'
	  !NOTE: change REALPART -> REAL
          write(423,*)qq,ae,be,REAL(conjg(Rh)*Rh)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
         end if
!--------------------------------------------------------------------------------
          end do! end of be exciton state sum 
          write(*,*)ae,deltaexc,RRh
!--------------------------------------------------------------------------------
          end do! end of ae exciton state sum
 !$OMP SECTION

  	do ae = 5, 8  !---------------------------------------------------------
              do be = 1,qq !Does not go past this? !Nexcc! exciton state sum
!---------------------------------------------------------  
! Rexc(c)-->biexc(a,b)=(|R2(a,b)|^2)*delta(Ec-(Ea+Eb))
! Rh to be calculated for each ae,be
               Rh=(0,0)
!               Rp=(0,0)
!               Rshh=(0,0)
!               Rshp=(0,0)

! w=Eexc(c)=Eexc(qq) the initial exciton energy          
! check if the two excitons ae and be satisfy energy conservation w=Ea+Eb
               dE=w-Eexc(ae)-Eexc(be)
               denom=dE*dE + gamma_b*gamma_b
               deltaexc=gamma_b/denom 
               
       if (deltaexc > 1/(gamma_b*ss)) then
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<           
               do ii = 1,Nexc(qq)
!/////////////////////////////////////////////////////////
                do kk = 1,Nexc(be)
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!==================hole spectator===================================
!    read l,n common for the Theta_{j} and Theta_{-j} parts
                l=ce(be,kk)
                n=ch(be,kk)

                psib=psi(be,kk)

!         the Theta_l Theta_{-n} part
             Ebln=e(l)-e(n)-Eexc(be)
             denomln=Ebln*Ebln+gamma_b*gamma_b

!         the Theta_{-l} Theta_{n} part
             Ebnl=e(l)-e(n)+Eexc(be)
             denomnl=Ebnl*Ebnl+gamma_b*gamma_b

!=========================================================
!=========================================================
!         the Theta_j part
!=========================================================
!  indices of exciton gamma or c
             j=ce(qq,ii)
             i=ch(qq,ii)

             psic=psi(qq,ii)

             psicc_b=conjg(psic)*psib
             psicccb=conjg(psic)*conjg(psib)
             
           do k = nLU,nLU+du
!---------------------------------------------------------
               
             jae=ccc(ae,k,i)

             if (jae.ne.0) then
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
             psia=psi(ae,jae)
            
             Eaeki=Eexc(ae)+e(k)-e(i)
             denom1=Eaeki*Eaeki+gamma_b*gamma_b
            
             Ejiw=e(j)-e(i)+w
             denom2=Ejiw*Ejiw+gamma_b*gamma_b

!         the Theta_l Theta_{-n} part            
             PV=Eaeki*Ejiw*Ebln/(denom1*denom2*denomln)

!   find the Theta_l Theta_{-n} VC matrix element
             jj=cveeh(j-nLU+1,l-nLU+1,n-(nHO-ddo)+1,k-nLU+1)
            
      f=(Eexc(ae)-e(k)+e(i))*(e(l)-e(n)-Eexc(be))*(e(j)-e(i)-Eexc(qq))

           Rh=Rh+VCeeh(jj)*PV*f*psia*psicccb!psia*conjg(psic)*conjg(psib)

! the Theta_{-l} Theta_{n} part
             PV=Eaeki*Ejiw*Ebnl/(denom1*denom2*denomnl)

!   find the Theta_{-l} Theta_{n} VC matrix element
!     NOTE use V(e2,h,e,e1)=V(e1,e,h,e2)^*
             jj=cveeh(k-nLU+1,l-nLU+1,n-(nHO-ddo)+1,j-nLU+1)
!            jj=cvehe(j-nLU+1,n-(nHO-ddo)+1,l-nLU+1,k-nLU+1)
            
!      f=(Eexc(ae)-e(k)+e(i))*(e(l)-e(n)-Eexc(be))*(e(j)-e(i)-Eexc(qq))

           Rh=Rh+conjg(VCeeh(jj))*PV*f*psia*psicc_b!psia*conjg(psic)*psib
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            end if
!---------------------------------------------------------	      
	   end do! end k loop
!=========================================================
!=========================================================
!         end of the Theta_j part
!=========================================================

!=========================================================
!         the Theta_-j part
!=========================================================
!  indices of exciton gamma or c
             i=ce(qq,ii)
             j=ch(qq,ii)

             psic=psi(qq,ii)

             psicb=psic*psib
             psic_cb=psic*conjg(psib)

           do k = nHO-ddo,nHO
!---------------------------------------------------------  
             jae=ccc(ae,i,k)

             if (jae.ne.0) then
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
             psia=psi(ae,jae)
                         
             Eaeki=Eexc(ae)-e(i)+e(k)
             denom1=Eaeki*Eaeki+gamma_b*gamma_b
            
             Ejiw=e(i)-e(j)-w
             denom2=Ejiw*Ejiw+gamma_b*gamma_b

!         the Theta_l Theta_{-n} part
            
             PV=Eaeki*Ejiw*Ebln/(denom1*denom2*denomln)

!   find the Theta_l Theta_{-n} VC matrix element
           jj=cvheh(j-(nHO-ddo)+1,l-nLU+1,n-(nHO-ddo)+1,k-(nHO-ddo)+1)
            
             f=(Eexc(ae)-e(i)+e(k))*(e(l)-e(n)-Eexc(be))*(e(i)-e(j)-w)
             
       Rh=Rh+VCheh(jj)*PV*f*conjg(psia)*psic_cb!conjg(psia)*psic*conjg(psib)

! the Theta_{-l} Theta_{n} part
             PV=Eaeki*Ejiw*Ebnl/(denom1*denom2*denomnl)

!   find the Theta_{-l} Theta_{n} VC matrix element
!   NOTE use V(h2,h,e,h1)=V(h1,e,h,h2)^*
         jj=cvheh(k-(nHO-ddo)+1,l-nLU+1,n-(nHO-ddo)+1,j-(nHO-ddo)+1)
!        jj=cvhhe(j-(nHO-ddo)+1,n-(nHO-ddo)+1,l-nLU+1,k-(nHO-ddo)+1)
            
!       f=(Eexc(ae)-e(i)+e(k))*(e(l)-e(n)-Eexc(be))*(e(i)-e(j)-w)

       Rh=Rh+conjg(VCheh(jj))*PV*f*conjg(psia)*psicb!conjg(psia)*psic*psib
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            end if
!---------------------------------------------------------	      
	   end do !end k loop
      
!=========================================================
!=========================================================
!         end of the Theta_-j part
!=========================================================

!///////////////////////////////////////////////////////////	      
	        end do !end of kk loop: b exciton index loop
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\	      
	       end do !end of ii loop: c or qq or gamma exciton index loop 
! if delta(w-Ea-Eb) is big enough there is a contribution to R2
!dimension factor to be added in .nb
         RRh=RRh+conjg(Rh)*Rh*deltaexc
!         write(423,*)'{',qq,',',ae,',',be,',',conjg(Rh)*Rh,'},'
	  !NOTE: change REALPART -> REAL
          write(423,*)qq,ae,be,REAL(conjg(Rh)*Rh)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
         end if
!--------------------------------------------------------------------------------
          end do! end of be exciton state sum 
          write(*,*)ae,deltaexc,RRh
!--------------------------------------------------------------------------------
          end do! end of ae exciton state sum
 !$OMP SECTION

  	do ae = 9, 12  !---------------------------------------------------------
              do be = 1,qq !Does not go past this? !Nexcc! exciton state sum
!---------------------------------------------------------  
! Rexc(c)-->biexc(a,b)=(|R2(a,b)|^2)*delta(Ec-(Ea+Eb))
! Rh to be calculated for each ae,be
               Rh=(0,0)
!               Rp=(0,0)
!               Rshh=(0,0)
!               Rshp=(0,0)

! w=Eexc(c)=Eexc(qq) the initial exciton energy          
! check if the two excitons ae and be satisfy energy conservation w=Ea+Eb
               dE=w-Eexc(ae)-Eexc(be)
               denom=dE*dE + gamma_b*gamma_b
               deltaexc=gamma_b/denom 
               
       if (deltaexc > 1/(gamma_b*ss)) then
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<           
               do ii = 1,Nexc(qq)
!/////////////////////////////////////////////////////////
                do kk = 1,Nexc(be)
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!==================hole spectator===================================
!    read l,n common for the Theta_{j} and Theta_{-j} parts
                l=ce(be,kk)
                n=ch(be,kk)

                psib=psi(be,kk)

!         the Theta_l Theta_{-n} part
             Ebln=e(l)-e(n)-Eexc(be)
             denomln=Ebln*Ebln+gamma_b*gamma_b

!         the Theta_{-l} Theta_{n} part
             Ebnl=e(l)-e(n)+Eexc(be)
             denomnl=Ebnl*Ebnl+gamma_b*gamma_b

!=========================================================
!=========================================================
!         the Theta_j part
!=========================================================
!  indices of exciton gamma or c
             j=ce(qq,ii)
             i=ch(qq,ii)

             psic=psi(qq,ii)

             psicc_b=conjg(psic)*psib
             psicccb=conjg(psic)*conjg(psib)
             
           do k = nLU,nLU+du
!---------------------------------------------------------
               
             jae=ccc(ae,k,i)

             if (jae.ne.0) then
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
             psia=psi(ae,jae)
            
             Eaeki=Eexc(ae)+e(k)-e(i)
             denom1=Eaeki*Eaeki+gamma_b*gamma_b
            
             Ejiw=e(j)-e(i)+w
             denom2=Ejiw*Ejiw+gamma_b*gamma_b

!         the Theta_l Theta_{-n} part            
             PV=Eaeki*Ejiw*Ebln/(denom1*denom2*denomln)

!   find the Theta_l Theta_{-n} VC matrix element
             jj=cveeh(j-nLU+1,l-nLU+1,n-(nHO-ddo)+1,k-nLU+1)
            
      f=(Eexc(ae)-e(k)+e(i))*(e(l)-e(n)-Eexc(be))*(e(j)-e(i)-Eexc(qq))

           Rh=Rh+VCeeh(jj)*PV*f*psia*psicccb!psia*conjg(psic)*conjg(psib)

! the Theta_{-l} Theta_{n} part
             PV=Eaeki*Ejiw*Ebnl/(denom1*denom2*denomnl)

!   find the Theta_{-l} Theta_{n} VC matrix element
!     NOTE use V(e2,h,e,e1)=V(e1,e,h,e2)^*
             jj=cveeh(k-nLU+1,l-nLU+1,n-(nHO-ddo)+1,j-nLU+1)
!            jj=cvehe(j-nLU+1,n-(nHO-ddo)+1,l-nLU+1,k-nLU+1)
            
!      f=(Eexc(ae)-e(k)+e(i))*(e(l)-e(n)-Eexc(be))*(e(j)-e(i)-Eexc(qq))

           Rh=Rh+conjg(VCeeh(jj))*PV*f*psia*psicc_b!psia*conjg(psic)*psib
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            end if
!---------------------------------------------------------	      
	   end do! end k loop
!=========================================================
!=========================================================
!         end of the Theta_j part
!=========================================================

!=========================================================
!         the Theta_-j part
!=========================================================
!  indices of exciton gamma or c
             i=ce(qq,ii)
             j=ch(qq,ii)

             psic=psi(qq,ii)

             psicb=psic*psib
             psic_cb=psic*conjg(psib)

           do k = nHO-ddo,nHO
!---------------------------------------------------------  
             jae=ccc(ae,i,k)

             if (jae.ne.0) then
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
             psia=psi(ae,jae)
                         
             Eaeki=Eexc(ae)-e(i)+e(k)
             denom1=Eaeki*Eaeki+gamma_b*gamma_b
            
             Ejiw=e(i)-e(j)-w
             denom2=Ejiw*Ejiw+gamma_b*gamma_b

!         the Theta_l Theta_{-n} part
            
             PV=Eaeki*Ejiw*Ebln/(denom1*denom2*denomln)

!   find the Theta_l Theta_{-n} VC matrix element
           jj=cvheh(j-(nHO-ddo)+1,l-nLU+1,n-(nHO-ddo)+1,k-(nHO-ddo)+1)
            
             f=(Eexc(ae)-e(i)+e(k))*(e(l)-e(n)-Eexc(be))*(e(i)-e(j)-w)
             
       Rh=Rh+VCheh(jj)*PV*f*conjg(psia)*psic_cb!conjg(psia)*psic*conjg(psib)

! the Theta_{-l} Theta_{n} part
             PV=Eaeki*Ejiw*Ebnl/(denom1*denom2*denomnl)

!   find the Theta_{-l} Theta_{n} VC matrix element
!   NOTE use V(h2,h,e,h1)=V(h1,e,h,h2)^*
         jj=cvheh(k-(nHO-ddo)+1,l-nLU+1,n-(nHO-ddo)+1,j-(nHO-ddo)+1)
!        jj=cvhhe(j-(nHO-ddo)+1,n-(nHO-ddo)+1,l-nLU+1,k-(nHO-ddo)+1)
            
!       f=(Eexc(ae)-e(i)+e(k))*(e(l)-e(n)-Eexc(be))*(e(i)-e(j)-w)

       Rh=Rh+conjg(VCheh(jj))*PV*f*conjg(psia)*psicb!conjg(psia)*psic*psib
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            end if
!---------------------------------------------------------	      
	   end do !end k loop
      
!=========================================================
!=========================================================
!         end of the Theta_-j part
!=========================================================

!///////////////////////////////////////////////////////////	      
	        end do !end of kk loop: b exciton index loop
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\	      
	       end do !end of ii loop: c or qq or gamma exciton index loop 
! if delta(w-Ea-Eb) is big enough there is a contribution to R2
!dimension factor to be added in .nb
         RRh=RRh+conjg(Rh)*Rh*deltaexc
!         write(423,*)'{',qq,',',ae,',',be,',',conjg(Rh)*Rh,'},'
	  !NOTE: change REALPART -> REAL
          write(423,*)qq,ae,be,REAL(conjg(Rh)*Rh)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
         end if
!--------------------------------------------------------------------------------
          end do! end of be exciton state sum 
          write(*,*)ae,deltaexc,RRh
!--------------------------------------------------------------------------------
          end do! end of ae exciton state sum
 !$OMP SECTION

  	do ae = 13, 16  !---------------------------------------------------------
              do be = 1,qq !Does not go past this? !Nexcc! exciton state sum
!---------------------------------------------------------  
! Rexc(c)-->biexc(a,b)=(|R2(a,b)|^2)*delta(Ec-(Ea+Eb))
! Rh to be calculated for each ae,be
               Rh=(0,0)
!               Rp=(0,0)
!               Rshh=(0,0)
!               Rshp=(0,0)

! w=Eexc(c)=Eexc(qq) the initial exciton energy          
! check if the two excitons ae and be satisfy energy conservation w=Ea+Eb
               dE=w-Eexc(ae)-Eexc(be)
               denom=dE*dE + gamma_b*gamma_b
               deltaexc=gamma_b/denom 
               
       if (deltaexc > 1/(gamma_b*ss)) then
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<           
               do ii = 1,Nexc(qq)
!/////////////////////////////////////////////////////////
                do kk = 1,Nexc(be)
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!==================hole spectator===================================
!    read l,n common for the Theta_{j} and Theta_{-j} parts
                l=ce(be,kk)
                n=ch(be,kk)

                psib=psi(be,kk)

!         the Theta_l Theta_{-n} part
             Ebln=e(l)-e(n)-Eexc(be)
             denomln=Ebln*Ebln+gamma_b*gamma_b

!         the Theta_{-l} Theta_{n} part
             Ebnl=e(l)-e(n)+Eexc(be)
             denomnl=Ebnl*Ebnl+gamma_b*gamma_b

!=========================================================
!=========================================================
!         the Theta_j part
!=========================================================
!  indices of exciton gamma or c
             j=ce(qq,ii)
             i=ch(qq,ii)

             psic=psi(qq,ii)

             psicc_b=conjg(psic)*psib
             psicccb=conjg(psic)*conjg(psib)
             
           do k = nLU,nLU+du
!---------------------------------------------------------
               
             jae=ccc(ae,k,i)

             if (jae.ne.0) then
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
             psia=psi(ae,jae)
            
             Eaeki=Eexc(ae)+e(k)-e(i)
             denom1=Eaeki*Eaeki+gamma_b*gamma_b
            
             Ejiw=e(j)-e(i)+w
             denom2=Ejiw*Ejiw+gamma_b*gamma_b

!         the Theta_l Theta_{-n} part            
             PV=Eaeki*Ejiw*Ebln/(denom1*denom2*denomln)

!   find the Theta_l Theta_{-n} VC matrix element
             jj=cveeh(j-nLU+1,l-nLU+1,n-(nHO-ddo)+1,k-nLU+1)
            
      f=(Eexc(ae)-e(k)+e(i))*(e(l)-e(n)-Eexc(be))*(e(j)-e(i)-Eexc(qq))

           Rh=Rh+VCeeh(jj)*PV*f*psia*psicccb!psia*conjg(psic)*conjg(psib)

! the Theta_{-l} Theta_{n} part
             PV=Eaeki*Ejiw*Ebnl/(denom1*denom2*denomnl)

!   find the Theta_{-l} Theta_{n} VC matrix element
!     NOTE use V(e2,h,e,e1)=V(e1,e,h,e2)^*
             jj=cveeh(k-nLU+1,l-nLU+1,n-(nHO-ddo)+1,j-nLU+1)
!            jj=cvehe(j-nLU+1,n-(nHO-ddo)+1,l-nLU+1,k-nLU+1)
            
!      f=(Eexc(ae)-e(k)+e(i))*(e(l)-e(n)-Eexc(be))*(e(j)-e(i)-Eexc(qq))

           Rh=Rh+conjg(VCeeh(jj))*PV*f*psia*psicc_b!psia*conjg(psic)*psib
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            end if
!---------------------------------------------------------	      
	   end do! end k loop
!=========================================================
!=========================================================
!         end of the Theta_j part
!=========================================================

!=========================================================
!         the Theta_-j part
!=========================================================
!  indices of exciton gamma or c
             i=ce(qq,ii)
             j=ch(qq,ii)

             psic=psi(qq,ii)

             psicb=psic*psib
             psic_cb=psic*conjg(psib)

           do k = nHO-ddo,nHO
!---------------------------------------------------------  
             jae=ccc(ae,i,k)

             if (jae.ne.0) then
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
             psia=psi(ae,jae)
                         
             Eaeki=Eexc(ae)-e(i)+e(k)
             denom1=Eaeki*Eaeki+gamma_b*gamma_b
            
             Ejiw=e(i)-e(j)-w
             denom2=Ejiw*Ejiw+gamma_b*gamma_b

!         the Theta_l Theta_{-n} part
            
             PV=Eaeki*Ejiw*Ebln/(denom1*denom2*denomln)

!   find the Theta_l Theta_{-n} VC matrix element
           jj=cvheh(j-(nHO-ddo)+1,l-nLU+1,n-(nHO-ddo)+1,k-(nHO-ddo)+1)
            
             f=(Eexc(ae)-e(i)+e(k))*(e(l)-e(n)-Eexc(be))*(e(i)-e(j)-w)
             
       Rh=Rh+VCheh(jj)*PV*f*conjg(psia)*psic_cb!conjg(psia)*psic*conjg(psib)

! the Theta_{-l} Theta_{n} part
             PV=Eaeki*Ejiw*Ebnl/(denom1*denom2*denomnl)

!   find the Theta_{-l} Theta_{n} VC matrix element
!   NOTE use V(h2,h,e,h1)=V(h1,e,h,h2)^*
         jj=cvheh(k-(nHO-ddo)+1,l-nLU+1,n-(nHO-ddo)+1,j-(nHO-ddo)+1)
!        jj=cvhhe(j-(nHO-ddo)+1,n-(nHO-ddo)+1,l-nLU+1,k-(nHO-ddo)+1)
            
!       f=(Eexc(ae)-e(i)+e(k))*(e(l)-e(n)-Eexc(be))*(e(i)-e(j)-w)

       Rh=Rh+conjg(VCheh(jj))*PV*f*conjg(psia)*psicb!conjg(psia)*psic*psib
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            end if
!---------------------------------------------------------	      
	   end do !end k loop
      
!=========================================================
!=========================================================
!         end of the Theta_-j part
!=========================================================

!///////////////////////////////////////////////////////////	      
	        end do !end of kk loop: b exciton index loop
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\	      
	       end do !end of ii loop: c or qq or gamma exciton index loop 
! if delta(w-Ea-Eb) is big enough there is a contribution to R2
!dimension factor to be added in .nb
         RRh=RRh+conjg(Rh)*Rh*deltaexc
!         write(423,*)'{',qq,',',ae,',',be,',',conjg(Rh)*Rh,'},'
	  !NOTE: change REALPART -> REAL
          write(423,*)qq,ae,be,REAL(conjg(Rh)*Rh)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
         end if
!--------------------------------------------------------------------------------
          end do! end of be exciton state sum 
          write(*,*)ae,deltaexc,RRh
!--------------------------------------------------------------------------------
          end do! end of ae exciton state sum
 !$OMP SECTION

  	do ae = 17, 86  !---------------------------------------------------------
              do be = 1,qq !Does not go past this? !Nexcc! exciton state sum
!---------------------------------------------------------  
! Rexc(c)-->biexc(a,b)=(|R2(a,b)|^2)*delta(Ec-(Ea+Eb))
! Rh to be calculated for each ae,be
               Rh=(0,0)
!               Rp=(0,0)
!               Rshh=(0,0)
!               Rshp=(0,0)

! w=Eexc(c)=Eexc(qq) the initial exciton energy          
! check if the two excitons ae and be satisfy energy conservation w=Ea+Eb
               dE=w-Eexc(ae)-Eexc(be)
               denom=dE*dE + gamma_b*gamma_b
               deltaexc=gamma_b/denom 
               
       if (deltaexc > 1/(gamma_b*ss)) then
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<           
               do ii = 1,Nexc(qq)
!/////////////////////////////////////////////////////////
                do kk = 1,Nexc(be)
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!==================hole spectator===================================
!    read l,n common for the Theta_{j} and Theta_{-j} parts
                l=ce(be,kk)
                n=ch(be,kk)

                psib=psi(be,kk)

!         the Theta_l Theta_{-n} part
             Ebln=e(l)-e(n)-Eexc(be)
             denomln=Ebln*Ebln+gamma_b*gamma_b

!         the Theta_{-l} Theta_{n} part
             Ebnl=e(l)-e(n)+Eexc(be)
             denomnl=Ebnl*Ebnl+gamma_b*gamma_b

!=========================================================
!=========================================================
!         the Theta_j part
!=========================================================
!  indices of exciton gamma or c
             j=ce(qq,ii)
             i=ch(qq,ii)

             psic=psi(qq,ii)

             psicc_b=conjg(psic)*psib
             psicccb=conjg(psic)*conjg(psib)
             
           do k = nLU,nLU+du
!---------------------------------------------------------
               
             jae=ccc(ae,k,i)

             if (jae.ne.0) then
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
             psia=psi(ae,jae)
            
             Eaeki=Eexc(ae)+e(k)-e(i)
             denom1=Eaeki*Eaeki+gamma_b*gamma_b
            
             Ejiw=e(j)-e(i)+w
             denom2=Ejiw*Ejiw+gamma_b*gamma_b

!         the Theta_l Theta_{-n} part            
             PV=Eaeki*Ejiw*Ebln/(denom1*denom2*denomln)

!   find the Theta_l Theta_{-n} VC matrix element
             jj=cveeh(j-nLU+1,l-nLU+1,n-(nHO-ddo)+1,k-nLU+1)
            
      f=(Eexc(ae)-e(k)+e(i))*(e(l)-e(n)-Eexc(be))*(e(j)-e(i)-Eexc(qq))

           Rh=Rh+VCeeh(jj)*PV*f*psia*psicccb!psia*conjg(psic)*conjg(psib)

! the Theta_{-l} Theta_{n} part
             PV=Eaeki*Ejiw*Ebnl/(denom1*denom2*denomnl)

!   find the Theta_{-l} Theta_{n} VC matrix element
!     NOTE use V(e2,h,e,e1)=V(e1,e,h,e2)^*
             jj=cveeh(k-nLU+1,l-nLU+1,n-(nHO-ddo)+1,j-nLU+1)
!            jj=cvehe(j-nLU+1,n-(nHO-ddo)+1,l-nLU+1,k-nLU+1)
            
!      f=(Eexc(ae)-e(k)+e(i))*(e(l)-e(n)-Eexc(be))*(e(j)-e(i)-Eexc(qq))

           Rh=Rh+conjg(VCeeh(jj))*PV*f*psia*psicc_b!psia*conjg(psic)*psib
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            end if
!---------------------------------------------------------	      
	   end do! end k loop
!=========================================================
!=========================================================
!         end of the Theta_j part
!=========================================================

!=========================================================
!         the Theta_-j part
!=========================================================
!  indices of exciton gamma or c
             i=ce(qq,ii)
             j=ch(qq,ii)

             psic=psi(qq,ii)

             psicb=psic*psib
             psic_cb=psic*conjg(psib)

           do k = nHO-ddo,nHO
!---------------------------------------------------------  
             jae=ccc(ae,i,k)

             if (jae.ne.0) then
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
             psia=psi(ae,jae)
                         
             Eaeki=Eexc(ae)-e(i)+e(k)
             denom1=Eaeki*Eaeki+gamma_b*gamma_b
            
             Ejiw=e(i)-e(j)-w
             denom2=Ejiw*Ejiw+gamma_b*gamma_b

!         the Theta_l Theta_{-n} part
            
             PV=Eaeki*Ejiw*Ebln/(denom1*denom2*denomln)

!   find the Theta_l Theta_{-n} VC matrix element
           jj=cvheh(j-(nHO-ddo)+1,l-nLU+1,n-(nHO-ddo)+1,k-(nHO-ddo)+1)
            
             f=(Eexc(ae)-e(i)+e(k))*(e(l)-e(n)-Eexc(be))*(e(i)-e(j)-w)
             
       Rh=Rh+VCheh(jj)*PV*f*conjg(psia)*psic_cb!conjg(psia)*psic*conjg(psib)

! the Theta_{-l} Theta_{n} part
             PV=Eaeki*Ejiw*Ebnl/(denom1*denom2*denomnl)

!   find the Theta_{-l} Theta_{n} VC matrix element
!   NOTE use V(h2,h,e,h1)=V(h1,e,h,h2)^*
         jj=cvheh(k-(nHO-ddo)+1,l-nLU+1,n-(nHO-ddo)+1,j-(nHO-ddo)+1)
!        jj=cvhhe(j-(nHO-ddo)+1,n-(nHO-ddo)+1,l-nLU+1,k-(nHO-ddo)+1)
            
!       f=(Eexc(ae)-e(i)+e(k))*(e(l)-e(n)-Eexc(be))*(e(i)-e(j)-w)

       Rh=Rh+conjg(VCheh(jj))*PV*f*conjg(psia)*psicb!conjg(psia)*psic*psib
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            end if
!---------------------------------------------------------	      
	   end do !end k loop
      
!=========================================================
!=========================================================
!         end of the Theta_-j part
!=========================================================

!///////////////////////////////////////////////////////////	      
	        end do !end of kk loop: b exciton index loop
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\	      
	       end do !end of ii loop: c or qq or gamma exciton index loop 
! if delta(w-Ea-Eb) is big enough there is a contribution to R2
!dimension factor to be added in .nb
         RRh=RRh+conjg(Rh)*Rh*deltaexc
!         write(423,*)'{',qq,',',ae,',',be,',',conjg(Rh)*Rh,'},'
	  !NOTE: change REALPART -> REAL
          write(423,*)qq,ae,be,REAL(conjg(Rh)*Rh)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
         end if
!--------------------------------------------------------------------------------
          end do! end of be exciton state sum 
          write(*,*)ae,deltaexc,RRh
!--------------------------------------------------------------------------------
          end do! end of ae exciton state sum
!----BOTTOM SECTION--------!


!$OMP END SECTIONS



!$OMP END PARALLEL


!        write(423,*)w,R2
         write(*,*)qq,w,RRh
         write(425,*)'{',qq,',',w,',',RRh,'},'

!         write(425,*)"}"

         close(425)
         close(423)

         end
!>>>>>>>>>>>>>>>>>>>>> END OF THE MAIN PROGRAM <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


!******************************************************************************
	SUBROUTINE find_r_dipmom(a1,b1,c1,emax,npw,XX,YY,ZZ,gx,gy,gz)
!******************************************************************************	
        implicit none

        real(8)     gx, gy, gz, gcut, twopi, gg          ! progections of k-vector in k-space
        real(8)     ggx(-300:300),ggy(-300:300),ggz(-300:300)
        real(8)     a1, b1, c1                           ! unit-vectors of the cell
        real(8)     emax

        integer     npw                                  ! max number of plane waves
        integer     ixmax, iymax, izmax, nnnn,n5 
        integer     ix(-300:300), ixy(-300:300,-300:300)
        integer     izy(-300:300,-300:300),iz(-300:300) 

        integer     i,j,k,i1,j1,k1
        integer     X,Y,Z

	integer     XX(npw),YY(npw),ZZ(npw)

!         integer     XX(100000),YY(100000),ZZ(100000)

!	integer, allocatable :: YY(npw)
!	integer, allocatable :: ZZ(npw)

!---------------------- beging the subroutine ---------------------------------------
              write(*,*) 'subroutine find_r_dipmom'
                write(*,*) 'npw=',npw

c           twopi = 2.0 * 3.1415926
	twopi=2D0*4D0*DATAN(1D0)

!       write(*,*) 'twopi=', twopi

           gx=twopi/a1
           gy=twopi/b1
           gz=twopi/c1
        !   gcut = sqrt(2.0*emax/27.2116)   !in atomic units (Hartree)
            gcut = dsqrt(emax/13.605826)/0.529177249  
!----find maximum of absolute values of vector r:
           ixmax = INT(gcut/gx)
           iymax = INT(gcut/gy)
           izmax = INT(gcut/gz)

         write(*,*) 'xmax, ymax, zmax, nnnn = ', ixmax, iymax, izmax, 
     c                                          nnnn
!----find x, y, z of the vector r (r on sphere) 
           nnnn = 0
         do i1 = 1, 2*izmax+1              ! z=0,1,...,Zmax,-Zmax,...-1
            if(i1.le.izmax+1) then
               i = i1-1         
            else
               i = i1-2*(izmax+1)
            endif
              ggz(i1)= gz * i 
              iz(i1) = 0
            do j1 = 1, 2*iymax+1           ! y=0,1,...,Ymax,-Ymax,...-1
               if(j1.le.iymax+1) then
                  j = j1-1         
               else
                  j = j1-2*(iymax+1)
               endif
                ggy(j1)= gy * j
                izy(i1,j1) = 0 
               do k1 = 1, 2*ixmax+1
                  if(k1.le.ixmax+1) then   ! x=0,1,...,Xmax,-Xmax,...-1
                     k = k1-1         
                  else
                     k = k1-2*(ixmax+1)
                  endif
                   ggx(k1)= gx * k 
                   gg = sqrt(ggx(k1)**2 + ggy(j1)**2 + ggz(i1)**2)
                  if (gg.le.gcut) then
                     if(k.gt.izy(i1,j1)) izy(i1,j1) = k
                     if(j.gt.iz(i1)) iz(i1) = j
                     nnnn = nnnn + 1
                  endif
               enddo
            enddo
          enddo

        write(*,*) 'xmax, ymax, zmax, nnnn = ', ixmax, iymax, izmax, 
     c                                          nnnn

!	    allocate(XX(npw))
!	    allocate(YY(npw))
!	    allocate(ZZ(npw))

             n5=0
!            n6=0
!            n7=0



        do i1 = 1, 2*izmax+1
            if(i1.le.izmax+1) then
               i = i1-1         
            else
               i = i1-2*(izmax+1)
            endif
            
               do j1 = 1, 2*iymax+1
               if(j1.le.iymax+1) then
                  j = j1-1         
               else
                  j = j1-2*(iymax+1)
               endif

               if(abs(j).le.iz(i1)) then

                  do k1 = 1, 2*ixmax+1
                      if(k1.le.ixmax+1) then
                          k = k1-1         
                      else
                          k = k1-2*(ixmax+1)
                      endif
                      if(abs(k).le.izy(i1,j1)) then
                            n5=n5+1
                         Z=i
                         ZZ(n5)=Z
                         Y= j
                         YY(n5)=Y
                         X= k
                         XX(n5)=X
                       if(n5.GT.npw) write(*,*) 'n5',n5,'npw',npw
                      endif
c                       if(n5.GT.npw) write(*,*) 'n5',n5,'npw',npw
                   enddo !k1
                endif  
             enddo !j1
         enddo !i1
         write(*,*) 'dimension of nubla n5 = ', n5

    !      do i = 1,300
    !       write(*,*) 'i=',i,' kx=',XX(i),' ky=',YY(i),' kz=',ZZ(i)
    !      enddo
    
   !  write(*,*) 'return from subroutine find_r_dipmom' 
        Return
       END
!>>>>>>>>>>>>>>> END OF SUBROUTINE find_r_dipmo <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
