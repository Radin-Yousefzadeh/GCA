!c    Written by MEHDI YOUSEFZADEH, June. 2018 - first edition
!c    Revised multiple time by the author since then
!c    Researchers who use this code for scientific purposes please cite the following
!c    two papers by the author  
!c    DOI: https://doi.org/10.3847/1538-4357/ac6de3 
!c    DOI: https://doi.org/10.3847/1538-4357/abd8d5
!     Feel free to send me an Email: m.yousefzadeh6@gmail.com
!
c---------------------------------------------------------------------------c
!
      program main
      use mpi
!
      IMPLICIT NONE
      integer *4 Nparticle 
      PARAMETER ( Nparticle = 2d4)
      real *8 w(4),dw(4),time,dt,maxdt,mindt,dt0,tscat,pi
      real *8 vpara,x,y,z,ek,vperp1,V(3),pitch0,pitch1,xrand,absV
      real *8 B0,V0,Q0,M0,T0,L0,E0,W0,Solar_Radius,jkl,VP ,map0,map1
      real *8 xmin_d,xmax_d,ymin_d,ymax_d,zmin_d,zmax_d,dx_d,dz_d
      integer *4 n,NTimeStep, NTimeMax, itime,itime_save,MAX_TIME
      character *8 cmy
      character *3 cpa,cpr
      integer *4 I_which_B,q_sign,t2,t1,count_rate,count_max,Interm
      integer*4 jj,io 
      integer*4 mype,npe,mpierr
      real *8 psi0,Sol_Radius0,kk,magnmom,gama,vperpO
!      ! Reading the data of the Numerical B field !
      integer i1,j1,k1
      integer, parameter :: Nx_mhd = 300 , Ny_mhd = 300 , Nz_mhd = 300
       real*8 :: Bx_mhd(Nx_mhd,Ny_mhd,Nz_mhd)
     & ,By_mhd(Nx_mhd,Ny_mhd,Nz_mhd),Bz_mhd(Nx_mhd,Ny_mhd,Nz_mhd)
      character*15 fileinx,fileiny,fileinz 
      character*55 path_file,newf !,path_mhd     
!     
      common/B_topology/psi0,Sol_Radius0,I_which_B
      common/Q_SIGNN/q_sign
      common/Mio/magnmom,gama
      common/domain/xmin_d,xmax_d,ymin_d,ymax_d,zmin_d,zmax_d,dx_d,dz_d
      common/mpis/mype,npe,mpierr
      common/color/VP
      common/DATASET/Bx_mhd,By_mhd,Bz_mhd
      common/norm/B0,V0,Q0,M0,T0,L0,E0
!c      common/MHDPARA/Nx_mhd,Ny_mhd,Nz_mhd      
!c      common/domain/xmin_d,xmax_d,ymin_d,ymax_d,zmin_d,zmax_d    
!
!c--------------------------------------------------------------c
!c                   Using MPI
      call MPI_INIT(mpierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, npe, mpierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, mype, mpierr)
!c--------------------------------------------------------------c
!c-------------------Numerical Unit System (NUS)----------------c
!c--------------------------------------------------------------c
      B0 = 1.6d-9                               !.13*1d-4 ! B = 100 Gs, unit of B, 1 Gs = 1.d-4 T
      V0 = 3.d8                                 ! speed of light, unit of speed (m/s)
      q_sign = 1
      Q0 = 1.602e-19                            ! electron charge, unit of charge
      if (q_sign.eq. -1)M0 = 9.1e-31            ! electron mass, unit of mass (kg)
      if (q_sign.eq.  1)M0 = 1.6726231e-27      ! proton mass, unit of mass (kg)
      T0 = M0/Q0/B0                             ! Gyroperiod, unit of time = (qB/m)^-1 , for proton ~ 1.d-4 (s) ; for V~0.01c , V/omga_ci ~ 313m Gyroradii  
      L0 = V0*T0                                ! unit of distance  , for proton(ion) ~ 30 (km)  ;  c/omga_ce ~ 17m , c/omga_ci ~ 31322m 
      E0 = V0*B0                                ! unit of electric field (E) , ion ~ 30 km.T 
      W0 = 1./2. *(M0 * V0**2)                  ! unit of kineitc energy ~ e-11 J ~ e8 ev ~ 100Mev (1ev = 1.6e-19 J)
      Solar_Radius = 6.96e8                     ! solar radius 6.95E8 m ~ 7.e5 km
      Solar_Radius = 100.*L0                    !(ion) ~ 3.d3 (km)
      Sol_Radius0 = Solar_Radius/L0             ! 2.32d4
      psi0 = 1.d0*B0*Solar_radius**3            ! THE CONSTANT FOR BIPOLAR MAGNETIC FLUX
      psi0 = psi0/B0/L0**3                      ! THE REMAINED CONSTANT IN THE NUS (1.d3 (km))
      print*,'Gyro Frequency (MHz)',1./T0/2/3.14/1000000 !M0*E0/(Q0*B0*B0*L0)
!      print*,'Gyro Radius (m, c*T0)',L0,Sol_Radius0 !,psi0
!
c------Simulation Box Size --------------------------c
!!!!!! 1 for Coronal loop numerical B.Field
      I_which_B = 1 
c------Initial injection type ------------c
      Interm = 0         ! 0-> delta-like , 1-> random  , 2-> continious  
!!!!!!
      if (I_which_B.eq.1)then
      xmin_d = -1.d5 - Sol_Radius0
      ymin_d = -1.d5 - Sol_Radius0
      zmin_d = -1.d5 - Sol_Radius0
      xmax_d =  1.d5 + Sol_Radius0
      ymax_d =  1.d5 + Sol_Radius0
      zmax_d =  1.d5 + Sol_Radius0
      NTimeStep = 10                            !! DATA will be saved each * time steps
      NTimeMax = 3.d4 !0000                     !! Maximum Number of Time Steps     
      endif
      MAX_TIME = 3.d4 !143 !0.0005d5               !!! in unit of T0      
!
c------Initial Time STEP and Limits of Time STEP------------c
      dt0   = 1.d-4
      maxdt = 10                                !!! MAXIMUM TIME STEP
      mindt = 1.d-8                             !!! MINIMUM TIME STEP
c----NOTE that dt will be adjusted in the BS subroutine------c
c     Division step in space      
!     dx_d = (xmax_d - xmin_d)/Nx
!      dz_d = (zmax_d - zmin_d)/Nz
c---- other parameters ----c
!
         path_file =  '/public4/home/sc51081/Mehdi/gc_run/' 
         fileinx   =  'Bxx1.dat'
         fileiny   =  'Byy1.dat'
         fileinz   =  'Bzz1.dat'          
       open(22,file=trim(path_file)//fileinx,status='old',action='read')
       open(23,file=trim(path_file)//fileiny,status='old',action='read')
       open(24,file=trim(path_file)//fileinz,status='old',action='read')
!
!!!!!!!!!! General Numerical Data Reading
!
       read(22,* ,IOSTAT=io)
     & (((Bx_mhd(i1,j1,k1),i1=1,Nx_mhd),j1=1,Ny_mhd) ,k1=1,Nz_mhd)

       read(23,*,IOSTAT=io)
     & (((By_mhd(i1,j1,k1),i1=1,Nx_mhd),j1=1,Ny_mhd) ,k1=1,Nz_mhd) 

       read(24,*,IOSTAT=io) 
     & (((Bz_mhd(i1,j1,k1),i1=1,Nx_mhd),j1=1,Ny_mhd) ,k1=1,Nz_mhd)
!
!c-------------------------Main Loop--------------------------c
!c   
      CALL SYSTEM_CLOCK(t1, count_rate, count_max)
!c      
3210   continue
!c
       write (cpa,'(I3.3)') mype
       write (cmy,'(I8.8)') MAX_TIME         
!c       
       if (Nparticle .gt. 10) then
!       open(177,file='init-'//cpa//'-'//cmy,status='unknown')
!c       
       open(961,file='final1-'//cpa//'-'//cmy 
     & ,status='new') 
       open(962,file='final2-'//cpa//'-'//cmy
     & ,status='new')
       open(963,file='final3-'//cpa//'-'//cmy
     & ,status='new')
       open(964,file='final4-'//cpa//'-'//cmy
     & ,status='new')
       endif
!
!!cc        file to save the B field that ONE particle is feeling
       if (Nparticle .eq. 1) then
       open(9009,file='bb_c-'//cpa//'-'//cmy,status='unknown')
       endif
!c
        jkl = 0
        kk  = 0
        jj  = 0
!c
       do 100 n = 1, Nparticle
       write (cpr,'(I3.3)') n
c       
c initialization of particle position & velocity
c
      call init_partcl(time,w,x,y,z,V,vperp1,vpara,pitch0,map0)
c
      !!!!!!!!!!!!! Injection Time types !!!!!!!!!!
      call injec_func(n,time,MAX_TIME,Nparticle,Interm)
!
         dt = dt0
         itime = 0
!         itime0 = 0
         itime_save = 0
!
c---------------for The purpose of diagnostics-----------------
!
        ek = vpara**2 + vperp1**2
        kk = kk + 1
!
!!!        goto 100
!
!c-------------- Loop for Present particle --------------------
10      continue 
!c         
!c    Use the BS method to calculate particle motion
!c    limit time steps if they are too large or too small
!c
        if (dt.gt.maxdt) dt = maxdt
        if (dt.lt.mindt) dt = mindt
!      
       call lorentz(time, w, dw)
!        
        itime = itime + 1
!       
       call doBS(w,dw,time,x,y,z,vperp1,vpara,ek,dt,map1,pitch1
     & ,itime)
!               
!        itime = itime + 1
!
! RANDOMLY CHANGING THE PITCH ANGLE DUE TO SCATTERING EFFECT
! CHANGING VPARA WILL BE EFFECTIVE CAUSE THAT'S THE ONLY NECESSARY
! PARAMETER FOR THE NEXT STEP
!
       tscat = 2900
       pi = 4.*atan(1.)
      if ( xrand(mype) .le. dt/tscat ) then
       absV   = dsqrt(vpara*vpara + vperp1*vperp1)
       vperpO = vperp1*vperp1
       pitch1 = xrand(mype) * pi
       vpara   = absV * dcos(pitch1)
       vperp1  = absV * dsin(pitch1)           
       magnmom = (magnmom/vperpO)*(vperp1*vperp1)
       pitch1  = pitch1 * 180/pi      !Degree
!     Put it back to equations for the next step
!     if no scattering needed then simply comment the following w(1)
!
!!!       w(1) = vpara
       endif
!
!     stop calculating if time is up or the partcle is out of the domain
        if(x.gt.Nx_mhd .or. x.lt. 1.02) goto 20
        if(y.gt.Ny_mhd .or. y.lt. 1.02) goto 20
        if(z.gt.Nz_mhd) goto 20 !.or. z.lt.1.1) goto 20
!
!        if(x.ge.xmax_d .or. x.le.xmin_d) goto 20
!        if(y.ge.zmax_d .or. y.le.zmin_d) goto 20
!        if(z.ge.zmax_d .or. z.le.zmin_d) goto 20
!         
        if( z .le. 1.02 )  then
         jj = jj+1
         goto 20
        endif
!c
       if(itime.ge.NTimeMax) goto 20
       if(time.gt.MAX_TIME)  goto 20           
!       if(time.gt.5 )  goto 20 
!c
         goto 10
!c---------------- End of loop for Present particle --------
 20      continue
!---------------for The purpose of diagnostics of each particle-----------------
!
!       itime0 = itime0 +1 
!
!c---------- Save the final characteristics of Present particle -------
!c-----------------------------------BOX1---------------------------------         
!        if (z .gt. 1.02) then 
        if ( x .ge. 149 .and. x .le. 155 .and. y .ge. 117 !Close to FP
     & .and. y .le. 129 .and. z .ge. 20 .and. z .le. 50 ) then
         write(961,'(12e15.5)')
     & time,x,y,z,vperp1,vpara,kk,pitch1,ek
!,sin(pitch1)*sin(pitch1)/B_c,ek
!         jkl =jkl + 1
!         print*, jkl
        endif
!          print*, mype !, ek
!c-----------------------------------BOX2----------------------------------        
         if ( x .ge. 150 .and. x .le. 156 .and. y .ge. 125
     & .and. y .le. 140 .and. z .ge. 50 .and. z .le. 70 ) then
         write(962,'(12e15.5)')
     & time,x,y,z,vperp1,vpara,kk,pitch1,ek
!    & time,x,y,z,vperp1,vpara,sqrt(vpara*vpara+vperp1*vperp1),kk,pitch1
!     & ,map1,ek
!sin(pitch1)*sin(pitch1)/B_c,ek
!
!!         jkl =jkl + 1
!!         print*, jkl
        endif
!
!c-----------------------------------BOX3---------------------------------------      
        if ( x .ge. 150 .and. x .le. 158 .and. y .ge. 146
     & .and. y .le. 165 .and. z .ge. 80. .and. z .le. 90 ) then
         write(963,'(12e15.5)')! File3       box3 removed replaced
     & time,x,y,z,vperp1,vpara,kk,pitch1,ek
!    & time,x,y,z,vperp1,vpara,sqrt(vpara*vpara+vperp1*vperp1),kk,pitch1
!     & ,map1,ek
!sin(pitch1)*sin(pitch1)/B_c,ek
!
!        jkl =jkl + 1
!         print*, jkl
        endif
!        endif
!c-----------------------------------BOX4---------------------------------------      
!         if (z .gt. 0.1) then 
         if ( x .ge. 150 .and. x .le. 157 .and. y .ge. 160
     & .and. y .le. 192 .and. z .ge. 90 .and. z .le. 100 ) then
!$$$ !            if (VP .lt. 0. ) then   
!           endif
         write(964,'(12e15.5)')
     & time,x,y,z,vperp1,vpara,kk,pitch1,ek
!    & time,x,y,z,vperp1,vpara,sqrt(vpara*vpara+vperp1*vperp1),kk,pitch1
!     & ,map1,ek
!
!        jkl =jkl + 1
!         print*, jkl
         endif
!        endif
!c--------------------------------------------------------------------------------
100   continue
!c-----------------------Main Loop----------------------------c      
!c end of the main loop here
       print *, "I'm process",mype,"out of",
     &      npe,"processes."
     &      ," Max_Time=",Max_time
!c    
  888  continue   
!        write (cmy,'(I4.4)') mype
!        write(222,'(10e15.5)') N_vv4 / (Nparticle*10)
!         
      write(*,*) 'Number of scaped particles', jj
      if (MAX_TIME .lt. 4300) MAX_TIME = MAX_TIME + 143
      if (MAX_TIME .gt. 4300 .and. MAX_TIME .lt. 3.d4)
     & MAX_TIME = MAX_TIME + 2850
      if ( MAX_TIME .lt. 3.d4 ) goto 3210
! If you only need one run - Only your maxtime run
!!      if ( MAX_TIME .lt. Max_TIME ) goto 3210
c
      CALL SYSTEM_CLOCK(t2, count_rate, count_max)    
      WRITE(*,*) 'Elapsed time(seconds) = ',real(t2-t1)/real(count_rate)
      call MPI_FINALIZE(mpierr)
      stop      
      end
!
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c--------------- subroutines & functions -------------------c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c------------------------------------------------------------
      subroutine init_partcl(time,w,x,y,z,V,vperp1,vpara,
     &  pitch0,map0)
      IMPLICIT NONE
      real *8 time,x,y,z,vpara,theta,phi,VP,
     & xrand ,w(4),vfactor,vperp1,absV,pitch0,cospa0,
     & pi,v0,vmax,magnmom,bb_c(3),ee(3),V(3),map0,mean,stdev,
     & term1(3),term2(3),term3(3),unite(3),dB_c(3),B_c,Ug(3),Ug2,Ug1,
     & term4(3),term5(3),term6(3),bGb(3),ueGb(3),norm_out,gama
      integer *4 I_which_B
      real *8 psi0,Sol_Radius0 
      integer*4 mype,npe,mpierr
      common/B_topology/psi0,Sol_Radius0,I_which_B
      common/Mio/magnmom,gama
      common/mpis/mype,npe,mpierr
      common/color/VP
!
         time = 0.d0
         pi = 4.*atan(1.)
c
         v0  = 0.0034d0 !0.07d0 ! 0.1c, 0.3c
         vmax= 100.d0*v0  !5.d0*v0
c
         vfactor  = xrand(mype)*(vmax/v0-1.d0) + 1.d0   ! v0 - vmax 
c         vfactor  = 1.d0
c
         theta = 2.d0*xrand(mype)-1.d0
         phi   = 2.d0*pi*xrand(mype)
c         
       if ( I_which_B.eq.1 ) then
!         
!!! POSITION ON THE UNTWISTED LOOP
! FP injection
!         x= ( 5*xrand(mype) ) + 147 !150.6
!         y= ( 5*xrand(mype) ) + 115 !118
!         z= ( 5*xrand(mype) ) + 12  !15
! LT injection
          x= ( 3*xrand(mype) ) + 151 !107 !80 !87 !156.8
          y= ( 3*xrand(mype) ) + 185 !178  !225!204 !162.9
          z= ( 3*xrand(mype) ) + 93  !37 !25 !15  !13
! FP & LT injection
!          x= ( 5*xrand(mype) ) + 151 !107 !80 !87 !156.8
!          y= ( 5*xrand(mype) ) + 155 !178  !225!204 !162.9
!          z= ( 5*xrand(mype) ) + 84 !37 !25 !15  !13
!!! POSITION ON THE TWISTED LOOP
! Bmax INJECTION
!          x= ( 5*xrand(mype) ) + 94    !236 !213 !173 !92
!          y= ( 5*xrand(mype) ) + 221   !123 !155 !165 !193
!          z= ( 5*xrand(mype) ) + 54    !09  !37  !16  !55
!
        endif
!c        
!cccccccccccccccccccccccccccccccc
!c
        call fields(x,y,z,bb_c,ee,unite,dB_c,term1,term2,term3,
     & term4,term5,term6,bGb,ueGb)        
c       
       B_c = dsqrt(bb_c(1)**2 + bb_c(2)**2 + bb_c(3)**2)
!
        mean  = 0. !0
        stdev = 0.14 ! 0.005
        call rand_normal(norm_out,mean,stdev) !( xrand(mype)*0.6 ) - 0.3 
         V(1) = norm_out
        mean  = 0.
        stdev = 0.14  ! 0.005
        call rand_normal(norm_out,mean,stdev) !( xrand(mype)*0.6 ) - 0.3 
         V(2) = norm_out 
        mean  = 0. !12
        stdev = 0.14 ! 0.005
        call rand_normal(norm_out,mean,stdev) !( xrand(mype)*0.6 ) - 0.3  
         V(3) = norm_out 
!
        absV = dsqrt(V(1)*V(1)+V(2)*V(2)+V(3)*V(3))
!
        vpara = V(1)*unite(1)+V(2)*unite(2)+V(3)*unite(3)
        vpara = 1*vpara  
!        if ( vpara .gt. 0 ) vpara = -vpara
        VP = vpara   !V(3)    !Color the different Direction
!               
        Ug(1)= V(1) - vpara*unite(3)   
        Ug(2)= V(2) - vpara*unite(2)   
        Ug(3)= V(3) - vpara*unite(1)   
!        
        Ug2 = Ug(1)*Ug(1)+Ug(2)*Ug(2)+Ug(3)*Ug(3)        
        Ug1 = dsqrt(absV*absV - vpara*vpara)
        if (vpara*vpara .eq. absV*absV .or.
     &  vpara*vpara .ge. absV*absV ) stop
!
        vpara = V(2)
        Ug1 = dsqrt(V(1)*V(1) + V(3)*V(3)) 
!
         w(1)= vpara 
         w(2)= x
         w(3)= y
         w(4)= z
c
c       !Magnetic moment is constant during the calculation
!         
!        magnmom= 0.5d0 * Ug2 / B_c   
        magnmom= 0.5d0 * Ug1*Ug1 / B_c   
!
!       Initial angel of particle's propagation
!!!!!!       vperp1 =  dsqrt(absV*absV - vpara*vpara)    !verp + dsqrt(2*magnmom*B_c)
        vperp1 = Ug1 !dsqrt(Ug2)
!
!       cospa0  = (V(1)*bb_c(1)+V(2)*bb_c(2)+V(3)*bb_c(3))/absV/B_c
        cospa0  = atan(vperp1/vpara)  !(Vt*unite(1)+Vt*unite(2)+Vt*unite(3))/Vt/B_c
!        The result will be in Degre
        pitch0 = cospa0 * 180/pi
        if (pitch0 .lt. 0) pitch0 = pitch0 + 180 
!!!
!c      IN ORDER TO CREATE LOSS CONE DISTRIBUTION USING MAXWELLIAN
!!!        if (pitch0 .lt. 30 .or. pitch0 .gt. 150 ) goto 1122
!
         map0 = sin(cospa0)*sin(cospa0)/B_c
!
        gama = 1
        if (absV .gt. 0.42) gama = 1/dsqrt(1-(absV*absV))
!
      return
      end
!------------------------------------------------------------
      subroutine injec_func(n,time,MAX_TIME,Nparticle,Interm)
      implicit none
      integer *4, INTENT(IN) :: n,MAX_TIME,Nparticle,Interm
      DOUBLE PRECISION, INTENT(OUT) :: time
      integer *4 interv
      DOUBLE PRECISION :: tintv,pintv,xrand
      integer *4 mype,npe,mpierr
      common/mpis/mype,npe,mpierr
!
      interv  = 10 
      tintv = 1500/interv   !3.d3/interv    !every 0.2sec up to the 2s(3.d3) !tintv = (MAX_TIME*0.1)/interv
      pintv = Nparticle/interv
      time = 0.
!
      IF(Interm.eq.1)THEN        !RANDOM INJECTION 
         if (n<= 1*pintv)                  time=tintv*0.10 +(0*tintv) !150  !300*xrand(mype)
         if (n> 1*pintv .and. n<= 2*pintv) time=tintv*0.24 +(1*tintv) !400  !300*xrand(mype)+300 
         if (n> 2*pintv .and. n<= 3*pintv) time=tintv*0.31 +(2*tintv) !500  !300*xrand(mype)+600
         if (n> 3*pintv .and. n<= 4*pintv) time=tintv*0.18 +(3*tintv) !800  !300*xrand(mype)+900  
         if (n> 4*pintv .and. n<= 5*pintv) time=tintv*0.34 +(4*tintv) !950  !300*xrand(mype)+1200
         if (n> 5*pintv .and. n<= 6*pintv) time=tintv*0.12 +(5*tintv) !1150 !300*xrand(mype)+1500
         if (n> 6*pintv .and. n<= 7*pintv) time=tintv*0.39 +(6*tintv) !1500 !300*xrand(mype)+1800
         if (n> 7*pintv .and. n<= 8*pintv) time=tintv*0.09 +(7*tintv) !2000 !300*xrand(mype)+2100
         if (n> 8*pintv .and. n<= 9*pintv) time=tintv*0.16 +(8*tintv) !2400 !300*xrand(mype)+2400
         if (n> 9*pintv)                   time=tintv*0.23 +(9*tintv) !2850 !300*xrand(mype)+2700 
!
       ELSE IF(Interm.eq.2)THEN  !CONSECUTIVE INJECTION
          if (n<= 1*pintv)                  time=0*tintv
          if (n> 1*pintv .and. n<= 2*pintv) time=1*tintv
          if (n> 2*pintv .and. n<= 3*pintv) time=2*tintv
          if (n> 3*pintv .and. n<= 4*pintv) time=3*tintv
          if (n> 4*pintv .and. n<= 5*pintv) time=4*tintv
          if (n> 5*pintv .and. n<= 6*pintv) time=5*tintv
          if (n> 6*pintv .and. n<= 7*pintv) time=6*tintv
          if (n> 7*pintv .and. n<= 8*pintv) time=7*tintv
          if (n> 8*pintv .and. n<= 9*pintv) time=8*tintv
          if (n> 9*pintv)                   time=9*tintv
       ELSE
          time = 0.
       ENDIF
!     
      if (time .gt. MAX_TIME) time= 0 !9*tintv !time/(MAX_TIME*0.5)
      return      
      End
!------------------------------------------------------------      
      subroutine rand_normal(norm_out,mean,stdev)
      implicit none
      DOUBLE PRECISION :: xrand,pi,theta,r,mean,stdev,norm_out !c,temp(2)
      integer *4 mype,npe,mpierr
      common/mpis/mype,npe,mpierr
!      common/sigmu/stdev,mean
      pi = 4.*atan(1.)
      IF(stdev <= 0.0d0) THEN
        WRITE(*,*) "Standard Deviation must be +ve"
      ELSE
!        CALL RANDOM_NUMBER(temp)
!        r=(-2.0d0*log(temp(1)))**0.5
!        theta = 2.0d0*pi*temp(2)
!        c= mean+stdev*r*sin(theta)
        r = (-2.0d0*log(xrand(mype)))**0.5
        theta = 2.0d0*pi*xrand(mype)
        norm_out = mean+stdev*r*sin(theta)
      END IF
      END
!----------------------------------------------------------------
!c-------------------------------------------------------------c
! This subroutine will devide the space of vpara and vperp into
! different square and give those an index in order to chategorize the
! initial velocity of particles.  "N_v" will show you the density (intensity)
! of particles within each square.
!c-------------------------------------------------------------c     
! REMOVED -- UNNESESSARY
c-------------------------------------------------------------c
! REMOVED -- UNNESESSARY
!c------------------------------------------------------------
!
      subroutine doBS(w,dw,time,x,y,z,vperp1,vpara,ek,dt,map1,pitch1
     & ,itime) 
      IMPLICIT NONE
      real *8 w(4),dw(4),w0(4),ws(4),TINY,error,time,dt,epsln,
     & dta,dte,bb_c(3),ee(3),u_x,u_y,u_z,B0,V0,Q0,M0,T0,L0,E0
      real *8 x,y,z,vdrift,vpara,ek,magnmom,vperp1,cospa1,pitch1,pi,
     & term1(3),term2(3),term3(3),unite(3),dB_c(3),B_c,map1,
     & term4(3),term5(3),term6(3),bGb(3),ueGb(3),gama
      integer*4 mype,npe,mpierr
      integer*4 ndim,i,q_sign ,itime,NTimeSte 
      external lorentz
      common/Mio/magnmom,gama
      common/Q_SIGNN/q_sign
      common/mpis/mype,npe,mpierr
      common/norm/B0,V0,Q0,M0,T0,L0,E0
       TINY  = 1.d-30 
       error = 1.d-6
       ndim  = 4     
       NTimeSte = 10
!
          pi = 4.*atan(1.)
!
         do i=1,ndim
            w0(i) = w(i)
            ws(i) = dabs(w(i)) + dabs(dt*dw(i)) + TINY
         end do 
c
         call bsstep(w,dw,ndim,time,dt,error,ws,dta,dte,lorentz)
         dt = dte !!! New Time Step
c
      vpara = w(1)
       x  = w(2)
       y  = w(3)
       z  = w(4)
!
       if ( x .gt. 300 .or. y .gt. 300 .or. z .gt. 300
     & .or. x .lt. 1.02 .or. y .lt. 1.02 .or. z .lt. 1.02 ) then 
!          print*, 'problem in dobs',x,y,z
        goto 1111
!          stop 
          endif
!                 
        call fields(x,y,z,bb_c,ee,unite,dB_c,term1,term2,term3,
     & term4,term5,term6,bGb,ueGb)
c
       B_c = dsqrt(bb_c(1)**2 + bb_c(2)**2 + bb_c(3)**2)
c
      u_x = term1(1) + q_sign*gama*((vpara**2/B_c) 
     & *term2(1)+(magnmom/B_c)*term3(1)
     & +(vpara/B_c)*(term4(1) + term5(1))
     & +(1./B_c)*term6(1))
c
      u_y = term1(2) + q_sign*gama*((vpara**2/B_c) 
     & *term2(2)+(magnmom/B_c)*term3(2)
     & +(vpara/B_c)*(term4(2) + term5(2))
     & +(1./B_c)*term6(2))      
c
      u_z = term1(3) + q_sign*gama*((vpara**2/B_c) 
     & *term2(3)+(magnmom/B_c)*term3(3)  
     & +(vpara/B_c)*(term4(3) + term5(3))
     & +(1./B_c)*term6(3))
c      
       vdrift = dsqrt(u_x**2 + u_y**2 + u_z**2)
!c
!c      !Etotal = guiding center energy + gyration energy     
!c
       ek = vpara**2 + vdrift**2 + magnmom*B_c*2.d0
!       
       vperp1 = vdrift + dsqrt(2*magnmom*B_c)
!
       cospa1  = atan(vperp1/vpara)  !(Vt*unite(1)+Vt*unite(2)+Vt*unite(3))/Vt/B_c
!        cospa1 = abs(cospa1)
!        The result will be in Degree
       pitch1 = cospa1 * 180/pi
       if (pitch1 .lt. 0) pitch1 = pitch1 + 180 
!
        map1 = sin(cospa1)*sin(cospa1)/B_c
!
!c For the purpose of Checking the GCA Code in each time step
!    rl/L << 1 -> then we can say that our approximation is valid  
!
       epsln = (vperp1/B_c)/
     & (B_c/(dsqrt(dB_c(1)**2 + dB_c(2)**2 + dB_c(3)**2)))
       if (epsln .gt. 0.02 ) then
       print*, 'Critical Approximation (epsilon):', epsln
!     &  ,vperp1*epsln*0.5
!       stop
       endif
!
1111     return
      end
c
c--------------------------------------------------      
c
      subroutine lorentz(t, w, dwdt)
      IMPLICIT NONE
      real *8 t, w(4), dwdt(4)
      real *8 x,y,z,magnmom,psi0,Sol_Radius0,  !M0,Q0,   
     & ee(3),bb_c(3),vpara,u_x,u_y,u_z, 
     & term1(3),term2(3),term3(3),unite(3),dB_c(3),B_c,
     & term4(3),term5(3),term6(3),bGb(3),ueGb(3),gama
      integer *4 I_which_B,q_sign
      common/B_topology/psi0,Sol_Radius0,I_which_B
      common/Mio/magnmom,gama
      common/Q_SIGNN/q_sign
!c      
      vpara = w(1)
       x = w(2)
       y = w(3)
       z = w(4)
c
       if ( x .gt. 300 .or. y .gt. 300 .or. z .gt. 300
     & .or. x .lt. 1.02 .or. y .lt. 1.02 .or. z .lt. 1.02 ) then 
!          print*, 'problem in lorentz',x,y,z
        goto 111
          endif
!
!       print*,'vpara',vpara
        call fields(x,y,z,bb_c,ee,unite,dB_c,term1,term2,term3,
     & term4,term5,term6,bGb,ueGb)        
c         
       B_c = dsqrt(bb_c(1)**2 + bb_c(2)**2 + bb_c(3)**2)
c
      dwdt(1) = (q_sign/gama)
     & *((ee(1)*unite(1)+ee(2)*unite(2)+ee(3)*unite(3))) 
     & - magnmom*(unite(1)*dB_c(1)+ unite(2)*dB_c(2)+unite(3)*dB_c(3))
     & + ( vpara*term1(1)*bGb(1) + vpara*term1(2)*bGb(2) 
     & + vpara*term1(3)*bGb(3) )
     & + ( term1(1)*ueGb(1) + term1(2)*ueGb(2) + term1(3)*ueGb(3) )
c 
      u_x = term1(1) + q_sign*gama*((vpara**2/B_c) 
     & *term2(1)+(magnmom/B_c)*term3(1) 
     & +(vpara/B_c)*(term4(1) + term5(1))
     & +(1./B_c)*term6(1))      
c
      u_y = term1(2) + q_sign*gama*((vpara**2/B_c) 
     & *term2(2)+(magnmom/B_c)*term3(2)
     & +(vpara/B_c)*(term4(2) + term5(2))
     & +(1./B_c)*term6(2))     
c
      u_z = term1(3) + q_sign*gama*((vpara**2/B_c) 
     & *term2(3)+(magnmom/B_c)*term3(3)
     & +(vpara/B_c)*(term4(3) + term5(3))
     & +(1./B_c)*term6(3))     
c
      dwdt(2) = u_x + vpara*unite(1)
      dwdt(3) = u_y + vpara*unite(2)
      dwdt(4) = u_z + vpara*unite(3)
c    
111    return
      end
!c 
!c----------------------------------------------------c      
!c  To give the electro-magnetic fields               c
!c----------------------------------------------------c
!c
      subroutine fields(x,y,z,bb_c,ee,unite,dB_c,term1,term2,term3,
     & term4,term5,term6,bGb,ueGb)
      IMPLICIT NONE
      real *8, intent(in) :: x,y,z
      real *8 bb_c(3),ee(3),dunitex(3),dunitey(3),dunitez(3)
      real *8 dBx_c(3), dBy_c(3), dBz_c(3),ueGb(3),B_c,gradB(9), !aax,hhy,
     & term1(3),term2(3),term3(3),unite(3),dB_c(3),bGb(3),       !rrr,rrr2
     & term4(3),term5(3),term6(3),dEx_c(3),dEy_c(3),dEz_c(3)  !,bb0,h1r,hz
      integer *4 I_which_B
      real *8 psi0,Sol_Radius0
      common/B_topology/psi0,Sol_Radius0,I_which_B
!c  
      if ( I_which_B.eq.1 ) then
!
c    Definition of E on Cartisian Coordinates
c
      ee(1) = 0.0d0   ! Ex 
      ee(2) = 0.0d0   ! Ey
      ee(3) = 0.0d0   ! Ez      
!!!!!!!!!!!!!!
      dEx_c(1) = 0   
      dEx_c(2) = 0  
      dEx_c(3) = 0
c
      dEy_c(1) = 0
      dEy_c(2) = 0
      dEy_c(3) = 0
c     
      dEz_c(1) = 0
      dEz_c(2) = 0
      dEz_c(3) = 0
!
c
      call read_mhd ( x,y,z,bb_c,gradB )
!c
      B_c = dsqrt(bb_c(1)**2 + bb_c(2)**2 + bb_c(3)**2)
!c            
      unite(1) = bb_c(1)/B_c
      unite(2) = bb_c(2)/B_c
      unite(3) = bb_c(3)/B_c
!c
      dBx_c(1) = gradB(1) !*psi0
      dBx_c(2) = gradB(2) !*psi0
      dBx_c(3) = gradB(3) !*psi0
!c
      dBy_c(1) = gradB(4) !*psi0
      dBy_c(2) = gradB(5) !*psi0
      dBy_c(3) = gradB(6) !*psi0
!c
      dBz_c(1) = gradB(7) !*psi0
      dBz_c(2) = gradB(8) !*psi0
      dBz_c(3) = gradB(9) !*psi0
!c
           endif
!c
      dB_c(1) = dBx_c(1)*unite(1) +dBy_c(1)*unite(2) +dBz_c(1)*unite(3) 
      dB_c(2) = dBx_c(2)*unite(1) +dBy_c(2)*unite(2) +dBz_c(2)*unite(3) 
      dB_c(3) = dBx_c(3)*unite(1) +dBy_c(3)*unite(2) +dBz_c(3)*unite(3) 
      dB_c(1) = dB_c(1) 
      dB_c(2) = dB_c(2) 
      dB_c(3) = dB_c(3) 
!      unite2 = dsqrt(dB(1)**2+dB(2)**2+dB(3)**2)          
!c      
!
      dunitex(1) = (1./B_c)*dBx_c(1)-(bb_c(1)/B_c*B_c)*dB_c(1) 
      dunitex(2) = (1./B_c)*dBx_c(2)-(bb_c(1)/B_c*B_c)*dB_c(2)
      dunitex(3) = (1./B_c)*dBx_c(3)-(bb_c(1)/B_c*B_c)*dB_c(3)
!c
      dunitey(1) = (1./B_c)*dBy_c(1)-(bb_c(2)/B_c*B_c)*dB_c(1)
      dunitey(2) = (1./B_c)*dBy_c(2)-(bb_c(2)/B_c*B_c)*dB_c(2)
      dunitey(3) = (1./B_c)*dBy_c(3)-(bb_c(2)/B_c*B_c)*dB_c(3)
!c
      dunitez(1) = (1./B_c)*dBz_c(1)-(bb_c(3)/B_c*B_c)*dB_c(1)
      dunitez(2) = (1./B_c)*dBz_c(2)-(bb_c(3)/B_c*B_c)*dB_c(2)
      dunitez(3) = (1./B_c)*dBz_c(3)-(bb_c(3)/B_c*B_c)*dB_c(3)
!
!c      Definition of E x B vector (First Term)
!
!      call crossab(ee,bb,UE)
        term1(1)= ( (ee(2)*unite(3))-(ee(3)*unite(2)) )/B_c
        term1(2)= ( (ee(3)*unite(1))-(ee(1)*unite(3)) )/B_c
        term1(3)= ( (ee(1)*unite(2))-(ee(2)*unite(1)) )/B_c
!c      print *, 'ue1',UE(1),'ue2',UE(2),'ue3',UE(3)
!c
!c      Definition of b x (b.Grad)b (Second Term)
!c
      term2(1)=unite(2)*unite(1)*dunitez(1)+unite(2)*unite(2)*dunitez(2)
     & + unite(2)*unite(3)*dunitez(3) -unite(3)*unite(1)*dunitey(1)
     & - unite(3)*unite(2)*dunitey(2) -unite(3)*unite(3)*dunitey(3)     
!     
      term2(2)=unite(3)*unite(1)*dunitex(1)+unite(3)*unite(2)*dunitex(2)
     & + unite(3)*unite(3)*dunitex(3) -unite(1)*unite(1)*dunitez(1)
     & - unite(1)*unite(2)*dunitez(2) -unite(1)*unite(3)*dunitez(3)
!     
      term2(3)=unite(1)*unite(1)*dunitey(1)+unite(1)*unite(2)*dunitey(2)
     & + unite(1)*unite(3)*dunitey(3) -unite(2)*unite(1)*dunitex(1)
     & - unite(2)*unite(2)*dunitex(2) -unite(2)*unite(3)*dunitex(3)
!
!c       Definition of (b.Grad)b
!
       bGb(1)= unite(1)*dunitex(1) + unite(2)*dunitex(2)
     & + unite(3)*dunitex(3)
       bGb(2)= unite(1)*dunitey(1) + unite(2)*dunitey(2)
     & + unite(3)*dunitey(3)
       bGb(3)= unite(1)*dunitez(1) + unite(2)*dunitez(2)
     & + unite(3)*dunitez(3)
!     
!c
!c      Definition of b x Gradb (Third Term)
!c
!
      term3(1) = unite(2)*dB_c(3) - unite(3)*dB_c(2)
      term3(2) = unite(3)*dB_c(1) - unite(1)*dB_c(3)
      term3(3) = unite(1)*dB_c(2) - unite(2)*dB_c(1)
!
!c
!c      Definition of b x (b.Grad)ue (Fourth Term)
!c
      term4(1) = unite(1)*( ( -term1(1)*dB_c(1) / B_c ) 
     & + ( -ee(3)*dunitey(1) + ee(2)*dunitez(1)
     & + unite(3)*dEy_c(1) - unite(2)*dEz_c(1) ) )
     
     & + unite(2)*( ( -term1(1)*dB_c(2) / B_c ) 
     & + ( -ee(3)*dunitey(2) + ee(2)*dunitez(2)
     & + unite(3)*dEy_c(2) - unite(2)*dEz_c(2) ) )
     
     & + unite(3)*( ( -term1(1)*dB_c(3) / B_c ) 
     & + ( -ee(3)*dunitey(3) + ee(2)*dunitez(3)
     & + unite(3)*dEy_c(3) - unite(2)*dEz_c(3) ) )
!!!!!--------------------------------!!
     
      term4(2) = unite(1)*( ( -term1(2)*dB_c(1) / B_c ) 
     & + ( ee(3)*dunitex(1) - ee(1)*dunitez(1)
     & - unite(3)*dEx_c(1) + unite(1)*dEz_c(1) ) )
     
     & + unite(2)*( ( -term1(2)*dB_c(2) / B_c ) 
     & + ( ee(3)*dunitex(2) - ee(1)*dunitez(2)
     & - unite(3)*dEx_c(2) + unite(1)*dEz_c(2) ) )
     
     & + unite(3)*( (-term1(2)*dB_c(3) / B_c) 
     & + ( ee(3)*dunitex(3) - ee(1)*dunitez(3)
     & - unite(3)*dEx_c(3) + unite(1)*dEz_c(3) ) )
!!!!!---------------------------------!!     

      term4(3) = unite(1)*( ( -term1(3)*dB_c(1) / B_c ) 
     & + ( -ee(2)*dunitex(1) + ee(1)*dunitey(1)
     & + unite(2)*dEx_c(1) - unite(1)*dEy_c(1) ) )
     
     & + unite(2)*( ( -term1(3)*dB_c(2) / B_c ) 
     & + ( -ee(2)*dunitex(2) + ee(1)*dunitey(2)
     & + unite(2)*dEx_c(2) - unite(1)*dEy_c(2) ) )
     
     & + unite(3)*( (-term1(3)*dB_c(3) / B_c) 
     & + ( -ee(2)*dunitex(3) + ee(1)*dunitey(3)
     & + unite(2)*dEx_c(3) - unite(1)*dEy_c(3) ) )
!
!c
!c      Definition of b x (ue.Grad)b (Fifth Term)
!c
      term5(1) = unite(2)*term1(1)*dunitez(1) 
     & + unite(2)*term1(2)*dunitez(2) + unite(2)*term1(3)*dunitez(3)
     & - unite(3)*term1(1)*dunitey(1) - unite(3)*term1(2)*dunitey(2)
     & - unite(3)*term1(3)*dunitey(3)
!
      term5(2) = unite(3)*term1(1)*dunitex(1) 
     & + unite(3)*term1(2)*dunitex(2) + unite(3)*term1(3)*dunitex(3)
     & - unite(1)*term1(1)*dunitez(1) - unite(1)*term1(2)*dunitez(2)
     & - unite(1)*term1(3)*dunitez(3) 
!
      term5(3) = unite(1)*term1(1)*dunitey(1) 
     & + unite(1)*term1(2)*dunitey(2) + unite(1)*term1(3)*dunitey(3)
     & - unite(2)*term1(1)*dunitex(1) - unite(2)*term1(2)*dunitex(2)
     & - unite(2)*term1(3)*dunitex(3)
!
!c       Definition of (ue.Grad)b
!     
      ueGb(1)= term1(1)*dunitex(1) 
     & + term1(2)*dunitex(2) + term1(3)*dunitex(3)
      ueGb(2)= term1(1)*dunitey(1) 
     & + term1(2)*dunitey(2) + term1(3)*dunitey(3)
      ueGb(3)= term1(1)*dunitez(1) 
     & + term1(2)*dunitez(2) + term1(3)*dunitez(3) 
!
!c
!c      Definition of b x (ue.Grad)ue (Sixth Term)
!c
      term6(1) = term1(1)*( ( -term1(1)*dB_c(1) / B_c ) 
     & + ( -ee(3)*dunitey(1) + ee(2)*dunitez(1)
     & + unite(3)*dEy_c(1) - unite(2)*dEz_c(1) ) )
     
     & + term1(2)*( ( -term1(1)*dB_c(2) / B_c ) 
     & + ( -ee(3)*dunitey(2) + ee(2)*dunitez(2)
     & + unite(3)*dEy_c(2) - unite(2)*dEz_c(2) ) )
     
     & + term1(3)*( ( -term1(1)*dB_c(3) / B_c ) 
     & + ( -ee(3)*dunitey(3) + ee(2)*dunitez(3)
     & + unite(3)*dEy_c(3) - unite(2)*dEz_c(3) ) )
     
!!!!!--------------------------------!!           
      term6(2) = term1(1)*( ( -term1(2)*dB_c(1) / B_c ) 
     & + ( ee(3)*dunitex(1) - ee(1)*dunitez(1)
     & - unite(3)*dEx_c(1) + unite(1)*dEz_c(1) ) )
     
     & + term1(2)*( ( -term1(2)*dB_c(2) / B_c ) 
     & + ( ee(3)*dunitex(2) - ee(1)*dunitez(2)
     & - unite(3)*dEx_c(2) + unite(1)*dEz_c(2) ) )
     
     & + term1(3)*( (-term1(2)*dB_c(3) / B_c) 
     & + ( ee(3)*dunitex(3) - ee(1)*dunitez(3)
     & - unite(3)*dEx_c(3) + unite(1)*dEz_c(3) ) )
     
!!!!!--------------------------------!!           
      term6(3) = term1(1)*( ( -term1(3)*dB_c(1) / B_c ) 
     & + ( -ee(2)*dunitex(1) + ee(1)*dunitey(1)
     & + unite(2)*dEx_c(1) - unite(1)*dEy_c(1) ) )
     
     & + term1(2)*( ( -term1(3)*dB_c(2) / B_c ) 
     & + ( -ee(2)*dunitex(2) + ee(1)*dunitey(2)
     & + unite(2)*dEx_c(2) - unite(1)*dEy_c(2) ) )
     
     & + term1(3)*( (-term1(3)*dB_c(3) / B_c) 
     & + ( -ee(2)*dunitex(3) + ee(1)*dunitey(3)
     & + unite(2)*dEx_c(3) - unite(1)*dEy_c(3) ) )
!
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!            c  read MHD data at one time frame c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_mhd ( x,y,z,bb_c,gradB )
      implicit none    
       real*8,intent(in) :: x,y,z !,time    
       real*8,intent(out) :: bb_c(3),gradB(9) !,dBx_c(3),dBy_c(3), dBz_c(3)       
       real*8 Bx_int,By_int,Bz_int
       real*8,dimension(6) :: Bx_int2,By_int2,Bz_int2!(6)
       real*8,dimension(8) :: Bxij,Byij,Bzij!(8) 
       Integer xi,yi,zi,xi_1,yi_1,zi_1
       Integer,parameter :: mhd_stp = 1      
       Integer, parameter :: Nx_mhd = 300 , Ny_mhd = 300 , Nz_mhd = 300
       real*8 :: Bx_mhd(Nx_mhd,Ny_mhd,Nz_mhd)
     & ,By_mhd(Nx_mhd,Ny_mhd,Nz_mhd),Bz_mhd(Nx_mhd,Ny_mhd,Nz_mhd)
      common/DATASET/Bx_mhd,By_mhd,Bz_mhd
!c      common/MHDPARA/Nx_mhd,Ny_mhd,Nz_mhd
!c
         if ( x .lt. 0 .and. x .ne. int(x) ) then 
             xi = int(x) -1 
         else
             xi = int(x) 
         endif
!
         if ( y .lt. 0 .and. y .ne. int(y) ) then 
             yi = int(y) -1
         else
             yi = int(y) 
         endif
!
         if ( z .lt. 0 .and. z .ne. int(z) ) then 
             zi = int(z) -1
         else
             zi = int(z) 
         endif
!                 
          xi_1 = xi +mhd_stp
          yi_1 = yi +mhd_stp
          zi_1 = zi +mhd_stp
!
          Bxij(1) = Bx_mhd(xi,yi,zi) 
          Byij(1) = By_mhd(xi,yi,zi)
          Bzij(1) = Bz_mhd(xi,yi,zi)
!
          Bxij(2) = Bx_mhd(xi_1,yi,zi) 
          Byij(2) = By_mhd(xi_1,yi,zi)
          Bzij(2) = Bz_mhd(xi_1,yi,zi)
!
          Bxij(3) = Bx_mhd(xi,yi_1,zi) 
          Byij(3) = By_mhd(xi,yi_1,zi)
          Bzij(3) = Bz_mhd(xi,yi_1,zi)
!
          Bxij(4) = Bx_mhd(xi,yi,zi_1) 
          Byij(4) = By_mhd(xi,yi,zi_1)
          Bzij(4) = Bz_mhd(xi,yi,zi_1)
!
          Bxij(5) = Bx_mhd(xi_1,yi_1,zi) 
          Byij(5) = By_mhd(xi_1,yi_1,zi)
          Bzij(5) = Bz_mhd(xi_1,yi_1,zi)
!
          Bxij(6) = Bx_mhd(xi_1,yi,zi_1) 
          Byij(6) = By_mhd(xi_1,yi,zi_1)
          Bzij(6) = Bz_mhd(xi_1,yi,zi_1)
!
          Bxij(7) = Bx_mhd(xi,yi_1,zi_1) 
          Byij(7) = By_mhd(xi,yi_1,zi_1)
          Bzij(7) = Bz_mhd(xi,yi_1,zi_1)
!
          Bxij(8) = Bx_mhd(xi_1,yi_1,zi_1) 
          Byij(8) = By_mhd(xi_1,yi_1,zi_1)
          Bzij(8) = Bz_mhd(xi_1,yi_1,zi_1)
!
!c         
      call interp3( x,y,z ,xi,yi,zi,mhd_stp ,Bxij,Bx_int )
      call interp3( x,y,z ,xi,yi,zi,mhd_stp ,Byij,By_int )
      call interp3( x,y,z ,xi,yi,zi,mhd_stp ,Bzij,Bz_int )
!
       bb_c(1) = Bx_int
       bb_c(2) = By_int
       bb_c(3) = Bz_int
!
!      x constant plane (1,3,4,7)
      call interp2( y,z, yi,zi,mhd_stp ,Bxij(1),Bxij(3)
     & ,Bxij(4),Bxij(7),Bx_int2(1) )
      call interp2( y,z, yi,zi,mhd_stp ,Byij(1),Byij(3)
     & ,Byij(4),Byij(7),By_int2(1) )
      call interp2( y,z, yi,zi,mhd_stp ,Bzij(1),Bzij(3)
     & ,Bzij(4),Bzij(7),Bz_int2(1) )
!     
!      x_1 constant plane (2,5,6,8)
      call interp2( y,z, yi,zi,mhd_stp ,Bxij(2),Bxij(5)
     & ,Bxij(6),Bxij(8),Bx_int2(2) )
      call interp2( y,z, yi,zi,mhd_stp ,Byij(2),Byij(5)
     & ,Byij(6),Byij(8),By_int2(2) )
      call interp2( y,z, yi,zi,mhd_stp ,Bzij(2),Bzij(5)
     & ,Bzij(6),Bzij(8),Bz_int2(2) )
!
!      y constant plane (1,2,4,6)
      call interp2( x,z, xi,zi,mhd_stp ,Bxij(1),Bxij(2)
     & ,Bxij(4),Bxij(6),Bx_int2(3) )
      call interp2( x,z, xi,zi,mhd_stp ,Byij(1),Byij(2)
     & ,Byij(4),Byij(6),By_int2(3) )
      call interp2( x,z, xi,zi,mhd_stp ,Bzij(1),Bzij(2)
     & ,Bzij(4),Bzij(6),Bz_int2(3) )
!     
!      y_1 constant plane (3,5,7,8)
      call interp2( x,z, xi,zi,mhd_stp ,Bxij(3),Bxij(5)
     & ,Bxij(7),Bxij(8),Bx_int2(4) )
      call interp2( x,z, xi,zi,mhd_stp ,Byij(3),Byij(5)
     & ,Byij(7),Byij(8),By_int2(4) )
      call interp2( x,z, xi,zi,mhd_stp ,Bzij(3),Bzij(5)
     & ,Bzij(7),Bzij(8),Bz_int2(4) )
!       
!      z constant plane (1,2,3,5)
      call interp2( x,y, xi,yi,mhd_stp ,Bxij(1),Bxij(2)
     & ,Bxij(3),Bxij(5),Bx_int2(5) )
      call interp2( x,y, xi,yi,mhd_stp ,Byij(1),Byij(2)
     & ,Byij(3),Byij(5),By_int2(5) )
      call interp2( x,y, xi,yi,mhd_stp ,Bzij(1),Bzij(2)
     & ,Bzij(3),Bzij(5),Bz_int2(5) )
!     
!      z_1 constant plane (4,6,7,8)
      call interp2( x,y, xi,yi,mhd_stp ,Bxij(4),Bxij(6)
     & ,Bxij(7),Bxij(8),Bx_int2(6) )
      call interp2( x,y, xi,yi,mhd_stp ,Byij(4),Byij(6)
     & ,Byij(7),Byij(8),By_int2(6) )
      call interp2( x,y, xi,yi,mhd_stp ,Bzij(4),Bzij(6)
     & ,Bzij(7),Bzij(8),Bz_int2(6) )
!
!
        call centdiff(xi,xi_1,Bx_int2(1),Bx_int2(2), gradB(1))
        call centdiff(yi,yi_1,By_int2(1),By_int2(2), gradB(2))
        call centdiff(zi,zi_1,Bz_int2(1),Bz_int2(2), gradB(3))
!
        ! For dBy/dx,y,z !
        call centdiff(xi,xi_1,Bx_int2(3),Bx_int2(4), gradB(4))
        call centdiff(yi,yi_1,By_int2(3),By_int2(4), gradB(5))
        call centdiff(zi,zi_1,Bz_int2(3),Bz_int2(4), gradB(6))
!
        ! For dBz/dx,y,z !
        call centdiff(xi,xi_1,Bx_int2(5),Bx_int2(6), gradB(7))
        call centdiff(yi,yi_1,By_int2(5),By_int2(6), gradB(8))
        call centdiff(zi,zi_1,Bz_int2(5),Bz_int2(6), gradB(9))
!
         return
        end
c
c----------------------------------------------------------
c!!!!!!!!!!!!!!! LENEAR INTERPOLATION !!!!!!!!!!!!!!!!!!
c----------------------------------------------------------        
c
      subroutine interp3( x,y,z ,xi,yi,zi,mhd_stp ,Bxij,Bx_int)
      implicit none
      Integer,Intent(IN) :: xi,yi,zi,mhd_stp          
      real*8,Intent(IN) :: x,y,z,Bxij(8) !,Byij(8),Bzij(8)
      real*8,Intent(OUT) :: Bx_int      
      real*8 Bx_old,delx,dely,delz,delt !,Bx_new
!       
       delx = ( x - xi ) / (( xi +mhd_stp ) - xi )
       dely = ( y - yi ) / (( yi +mhd_stp ) - yi )
       delz = ( z - zi ) / (( zi +mhd_stp ) - zi )
       delt = 0     !(t - ti))/(ti+1 - ti))          ! For only one snap shot of the
!       field the delt = 0 cause there is no second snapshot to be taken into account
!
           !! For BX(x,y,z)!!
         Bx_old = (1-delx)*(1-dely)*(1-delz)*Bxij(1) 
     & + delx*(1-dely)*(1-delz)*Bxij(2) +dely*(1-delx)*(1-delz)*Bxij(3)
     & + delz*(1-delx)*(1-dely)*Bxij(4) +delx*dely*(1-delz)*Bxij(5)
     & + delx*delz*(1-dely)*Bxij(6) +dely*delz*(1-delx)*Bxij(7)
     & + delx*dely*delz*Bxij(8) 
!        
!        Bx_new = (1-delx)*(1-dely)*     +delx*(1-dely)*   +dely*(1-delx)*
!     & +delx*dely*i
         Bx_int = (1-delt)*Bx_old !+ delt*B_new
!
         return
      end
!c==========================================================      
      subroutine interp2( xx,yy, xxi,yyi,mhd_stp ,B1,B2 ,B3,B4,B_int )
      implicit none
      Integer,Intent(IN) :: xxi,yyi,mhd_stp
      real*8,Intent(IN) :: xx,yy,B1,B2,B3,B4
      real*8,Intent(OUT) :: B_int
      real*8 Bx_old,delx,dely,delt !,Bx_new
!      save !B_int
!     & ,Bz_old,Bz_new,Bz_int,By_old,By_new,By_int
!
       delx = ( xx - xxi ) / (( xxi +mhd_stp ) - xxi )
       dely = ( yy - yyi ) / (( yyi +mhd_stp ) - yyi )
!       delz = ( z - zi ) / (( zi +mhd_stp ) - zi )
       delt = 0     !(t - ti))/(ti+1 - ti))          ! For only one snap shot of the
!       field the delt = 0 cause there is no second snapshot to be taken into account
!
           !! For BX(x,y,z)!!
        Bx_old = (1-delx)*(1-dely)*B1 + delx*(1-dely)*B2
     &   + dely*(1-delx)*B3 + delx*dely*B4
!c     
         B_int = (1-delt)*Bx_old !+ delt*B_new
!
         return
      end      
!
c---------------------------------------------------------
c calculate the partial difference of Bx,By,Bz
c---------------------------------------------------------      
c 
       subroutine centdiff(x_1,x_2,B_x1,B_x2,dBdx)
          real*8,Intent(IN) :: B_x1,B_x2
          integer,Intent(IN) :: x_1,x_2
          real*8,Intent(OUT) :: dBdx          
!          real*8,intent(out) :: dBdx
!           print*,x_1,x_2,B_x1,B_x2
!           pause
          dBdx = 1.d0*(B_x2 - B_x1)/(x_2 - x_1)
!
          return
        end      
c      
c------------------------------------------------------------------c
c- NOTE, the following subroutines are standard codes, DONOT REVISE-c
c------------------------------------------------------------------c
c
      SUBROUTINE BSSTEP(Y,DYDX,NV,X,HTRY,EPS,YSCAL,HDID,HNEXT,DERIVS)
      IMPLICIT NONE
      INTEGER*4 NMAX,IMAX,NUSE,NV
      REAL*8 ONE,SHRINK,GROW
      PARAMETER (NMAX=10,IMAX=11,NUSE=7,ONE=1.e0,SHRINK=.95e0,GROW=1.2e0
     *)
      REAL*8 Y(*),DYDX(*),YSCAL(*),YERR(NMAX),
     *    YSAV(NMAX),DYSAV(NMAX),YSEQ(NMAX)
      REAL*8 X,HTRY,EPS,HDID,HNEXT,H,XSAV,XEST,ERRMAX,flag
      INTEGER*4 NSEQ(IMAX),I,J
      EXTERNAL DERIVS
      DATA NSEQ /2,4,6,8,12,16,24,32,48,64,96/      ! XL
      SAVE
c      NSEQ =(/2,4,6,8,12,16,24,32,48,64,96/)

      H=HTRY
      XSAV=X
      DO 11 I=1,NV
         YSAV(I)=Y(I)
         DYSAV(I)=DYDX(I)
11    CONTINUE
 1    DO 10 I=1,IMAX
         CALL MMID(YSAV,DYSAV,NV,XSAV,H,NSEQ(I),YSEQ,DERIVS)
         XEST=(H/NSEQ(I))**2
         CALL RZEXTR(I,XEST,YSEQ,Y,YERR,NV,NUSE)
         ERRMAX=0.
        DO 12 J=1,NV
           ERRMAX=MAX(ERRMAX,DABS(YERR(J)/YSCAL(J)))
12      CONTINUE
        ERRMAX=ERRMAX/EPS
        IF(ERRMAX.LT.ONE) THEN
           X=X+H
           HDID=H
           IF(I.EQ.NUSE)THEN
              HNEXT=H*SHRINK
           ELSE IF(I.EQ.NUSE-1)THEN
              HNEXT=H*GROW
           ELSE
              HNEXT=(H*NSEQ(NUSE-1))/NSEQ(I)
           ENDIF
           RETURN
        ENDIF
10    CONTINUE
      H=0.25*H/2**((IMAX-NUSE)/2)
c      H=0.25d0*H/2.**((IMAX-NUSE)/2.)
c      IF(X+H.EQ.X) then
c	flag = 1.
c	return
c      endif
c      IF(X+H.EQ.X)PAUSE 'Step size underflow.'
      IF(X+H.EQ.X) Print*, 'Step size underflow.'
      GOTO 1
      END
      
      
      SUBROUTINE MMID(Y,DYDX,NVAR,XS,HTOT,NSTEP,YOUT,DERIVS)
      IMPLICIT NONE
      INTEGER*4 NMAX,NVAR
      PARAMETER (NMAX=10)
      REAL*8 Y(*),DYDX(*),YOUT(*),YM(NMAX),YN(NMAX)
      REAL*8 XS,HTOT,H,X,SWAP,H2
      INTEGER*4 NSTEP,I,N
      EXTERNAL DERIVS
      SAVE                   ! XL

      H=HTOT/NSTEP
      DO 11 I=1,NVAR
         YM(I)=Y(I)
         YN(I)=Y(I)+H*DYDX(I)
11    CONTINUE
      X=XS+H
      CALL DERIVS(X,YN,YOUT)
      H2=2.*H
      DO 13 N=2,NSTEP
         DO 12 I=1,NVAR
            SWAP=YM(I)+H2*YOUT(I)
            YM(I)=YN(I)
            YN(I)=SWAP
12      CONTINUE
        X=X+H
        CALL DERIVS(X,YN,YOUT)
13    CONTINUE
      DO 14 I=1,NVAR
         YOUT(I)=0.5*(YM(I)+YN(I)+H*YOUT(I))
14    CONTINUE
      RETURN
      END
     
      
      SUBROUTINE RZEXTR(IEST,XEST,YEST,YZ,DY,NV,NUSE)
      IMPLICIT NONE
      INTEGER*4 IMAX,NMAX,NCOL,NV
      PARAMETER (IMAX=11,NMAX=10,NCOL=7)
      REAL*8 X(IMAX),YEST(*),YZ(*),DY(*),D(NMAX,NCOL),FX(NCOL)
      REAL*8 XEST,YY,V,C,B,B1,DDY
      INTEGER*4 IEST,NUSE,J,M1,K
      SAVE !d,x               ! XL
      X(IEST)=XEST
      IF(IEST.EQ.1) THEN
         DO 11 J=1,NV
            YZ(J)=YEST(J)
            D(J,1)=YEST(J)
            DY(J)=YEST(J)
11      CONTINUE
      ELSE
         M1=MIN(IEST,NUSE)
         DO 12 K=1,M1-1
            FX(K+1)=X(IEST-K)/XEST
12      CONTINUE
        DO 14 J=1,NV
           YY=YEST(J)
           V=D(J,1)
           C=YY
           D(J,1)=YY
           DO 13 K=2,M1
              B1=FX(K)*V
              B=B1-C
              IF(B.NE.0.) THEN
                 B=(C-V)/B
                 DDY=C*B
                 C=B1*B
              ELSE
                 DDY=V
              ENDIF
              V=D(J,K)
              D(J,K)=DDY
              YY=YY+DDY
13        CONTINUE
          DY(J)=DDY
          YZ(J)=YY
14      CONTINUE
      ENDIF
      RETURN
      END 

      

c-------------------------------------------------      
      function xrand(np)
      
      implicit none
      real*8 xrand,ran3
      real*4 rand
      integer*4 idum,ic,np,seed(0:100)
      data ic/0/
      save idum,ic,seed

      ic=ic+1
      if( ic.eq.1 )then
        seed(np) = 100000*np +
     &  nint(     10000.0*rand(np+1000000)
     &        +1000.0*rand(np+10000)
     &        +100.0*rand(np+100)
     &        +10.0*rand(np+00)
     &        +1.0*rand(np+10))
        idum=-seed(np)
        xrand=ran3(idum)
      end if
      
      xrand=ran3(idum)

      return
      end

      
      FUNCTION ran3(idum)
      implicit none
      INTEGER*4 idum
      INTEGER*4 MBIG,MSEED,MZ
      REAL*8 ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
      INTEGER*4 i,iff,ii,inext,inextp,k
      INTEGER*4 mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/

      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif

      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC

      return
      END
