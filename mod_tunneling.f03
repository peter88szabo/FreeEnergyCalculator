     !================================================================
     ! Everything in atomic units in this module to avoid Planck const
     !================================================================
      module mod_tunneling
      private

      public :: Wigner
      public :: Bell
      public :: Skodje_Truhlar
      public :: Skodje_Truhlar_exact 
      public :: Eckart 
      private :: simpson_integrate

      contains

     !==========================================================================
      subroutine Wigner(beta,omega,V0,kappa)
     !==========================================================================
     !     omega  = imag omega
     !     V0    = barrier height
     !     beta  = 1/(kB*T)
     !==========================================================================
      implicit none
      real*8 :: beta,omega,V0,kappa,u
      real*8 :: dum0,dum1,dum2

     !in atomic units
      u=omega*beta

      kappa = 1.d0 + u*u/24.d0

      end subroutine
     !==========================================================================



     !==========================================================================
      subroutine Skodje_Truhlar(beta,omega,V0,kappa)
     !==========================================================================
     !     It assumes exoergic reaction
     !     
     !     omega  = imag omega
     !     V0    = barrier height
     !     beta  = 1/(kB*T)
     !==========================================================================
      implicit none
      real*8,parameter :: pi=4.0*atan(1.d0)
      real*8,parameter :: twopi=2.d0*pi
      real*8 :: beta,omega,V0,kappa,alpha
      real*8 :: dum0,dum1,dum2

     !in atomic units
      alpha=twopi/omega

      if(beta > alpha) then
         !when imag omega is large and low T           
         kappa = (dexp((beta-alpha)*V0)-1.d0)*beta/(beta-alpha)

      elseif(beta < alpha) then
         !when imag omega is small and high T           
         dum0 = pi*beta/alpha
         dum1 = dum0/dsin(dum0)
         dum2 = dexp((beta-alpha)*V0)*beta/(beta-alpha)
         kappa = dum1+dum2
      else
         stop "Wrong mode-argument in Skodje_Truhlar routine"
      endif

      end subroutine
     !==========================================================================

     !==========================================================================
      subroutine Skodje_Truhlar_exact(beta,omega,V0,kappa)
     !==========================================================================
     !     It assumes exoergic reaction
     !     
     !     omega  = imag omega
     !     V0    = barrier height
     !     beta  = 1/(kB*T)
     !==========================================================================
      implicit none
      real*8,parameter :: pi=4.0*atan(1.d0)
      real*8,parameter :: twopi=2.d0*pi
      integer, parameter :: nmax=100
      integer :: n
      real*8 :: beta,omega,V0,kappa,alpha
      real*8 :: res,denom,numerator,alter,dum

     !in atomic units
      alpha=twopi/omega

      res=0.d0

      do n=0,nmax

         numerator = 1.d0 - dexp((beta - (n+1)*alpha)*v0)
         denom = (n+1)*alpha - beta
         dum = 1.d0/(n*alpha + beta)

         alter = 1.d0
         if(mod(n,2) == 1) alter = -1.d0

         res = res + alter*beta*(numerator/denom + dum) 

      enddo

      kappa = res

      end subroutine
     !==========================================================================



     !==========================================================================
      subroutine Bell(beta,omega,V0,kappa)
     !==========================================================================
     !     omega  = imag omega
     !     V0    = barrier height
     !     beta  = 1/(kB*T)
     !==========================================================================
      implicit none
      real*8,parameter :: pi=4.0*atan(1.d0)
      real*8,parameter :: twopi=2.d0*pi
      real*8 :: beta,omega,V0,kappa,u
      real*8 :: dum1,dum2,dum3,a2,a3,a4,hgf

     !in atomic units (no planck constant, hbar = 1)
      u=omega*beta

      dum1 = 0.5d0*u/dsin(0.5d0*u)
      dum2 = u/(u-twopi)
      dum3 = dexp(V0*beta/dum2)

      a2 = 1.d0-u/twopi
      a3 = a2+1.d0
      a4 = -dexp(-twopi*V0*beta/u)

!      call hygfx(1.d0,a2,a3,a4,hgf)
!      call hygfx(2.d0,2.0,0.9,-2.2,hgf)

     !Together with hypergeometric functions we have the exact Bell kappa 
      kappa=abs(dum1) !-dum2*dum3!*hgf

      end subroutine
     !==========================================================================



     !==========================================================================
      subroutine Eckart(beta,omega,Vf,Vb,dE,Emax,kappa1)
     !==========================================================================
     !     omega  = imag omega
     !     Vf    = forward barrier height
     !     Vb    = backward barrier height
     !     dE    = energy grid for the numerical integration
     !     beta  = 1/(kB*T)
     !==========================================================================
      implicit none
      real*8,parameter :: pi=4.0*atan(1.d0)
      real*8,parameter :: twopi=2.d0*pi
      integer :: i, nEmax
      real*8  :: Emax
      real*8  :: beta,omega,Vf,Vb,dE
      real*8  :: kappa1,kappa2
      real*8  :: u,E,E0,csi
      real*8  :: A,B,D,alpha1,alpha2
      real*8  :: Ptun1, Ptun2,Res
      real*8,allocatable :: toInt1(:),toInt2(:)

      nEmax=nint(Emax/dE)+1
      allocate(toInt1(nEmax),toInt2(nEmax)) 

     !in atomic units (no planck constant, hbar = 1)
      u=omega*beta

      alpha1 = twopi*Vf/omega
      alpha2 = twopi*Vb/omega

      A = dsqrt(Vf) + dsqrt(Vb)
      A = A*A
      B = Vf-Vb
      D = dsqrt(((B*B - A*A)**2)/(A*A*A)/8.d0)/omega  


      E0=0.d0
      do i=1,nEmax
         E = E0 +(i-1)*dE

         csi = E/Vf 

         call Tunprop1(alpha1,alpha2,csi,Ptun1)

         !call Tunprop2(A,B,D,E,Ptun2)

         toInt1(i) = dexp(-E*beta)  * Ptun1
         !toInt2(i) = dexp(-E*beta)  * Ptun2
      enddo

      call simpson_integrate(nEmax,1,nEmax,dE,ToInt1,Res)
      !write(6,*) "Res1 after integrate Eckart: ",Res
      kappa1=Res*dexp(beta*Vf)*beta

      !call simpson_integrate(nEmax,1,nEmax,dE,ToInt2,Res)
      !write(6,*) "Res2 after integrate Eckart: ",Res
      !kappa2=Res*dexp(beta*Vf)*beta

      end subroutine
     !==========================================================================





     !==========================================================================
      subroutine Tunprop1(alpha1,alpha2,csi,Ptun)
     !==========================================================================
     !     omega  = imag omega
     !     Vf    = forward barrier height
     !     Vb    = backward barrier height
     !     beta  = 1/(kB*T)
     !==========================================================================
      implicit none
      real*8,parameter :: pi=4.0*atan(1.d0)
      real*8,parameter :: twopi=2.d0*pi
      real*8 :: alpha1,alpha2,csi,Ptun
      real*8 :: denom,twopiA,twopiB,twopiD
      real*8 :: csh_amb,csh_apb,csh_d

      denom = 1.d0/dsqrt(alpha1) + 1.d0/dsqrt(alpha2) 

      twopiA = 2.d0*dsqrt(alpha1*csi) / denom
      twopiB = 2.d0*dsqrt(alpha1*csi + dabs(alpha1-alpha2)) / denom
      twopiD = 2.d0*dsqrt(dabs(alpha1*alpha2 - pi*pi/4.d0))

      csh_amb = dcosh(twopiA-twopiB)
      csh_apb = dcosh(twopiA+twopiB)
      csh_d   = dcosh(twopiD)

!     Tunneling probability at energy (csi = E/Vf)
      Ptun = 1.d0 - (csh_amb + csh_d)/(csh_apb + csh_d)

      end subroutine
     !==========================================================================




     !==========================================================================
      subroutine Tunprop2(A,B,D,E,Ptun)
     !==========================================================================
     !     omega = imag omega
     !     Vf    = forward barrier height
     !     Vb    = backward barrier height
     !     beta  = 1/(kB*T)
     !==========================================================================
      implicit none
      real*8,parameter :: pi=4.0*atan(1.d0)
      real*8,parameter :: twopi=2.d0*pi
      real*8 :: A,B,D,E,Ptun,alpha,beta,delta
      real*8 :: csh_amb,csh_apb,csh_d


      alpha = twopi*D*dsqrt(2.d0*E)

      beta = twopi*D*dsqrt(2.d0*(E-B))

      delta = twopi*D*dsqrt(2.d0*A - 1.d0/(2.d0*D*D))

      csh_amb = dcosh(alpha-beta)
      csh_apb = dcosh(alpha+beta)
      csh_d   = dcosh(delta)

!     Tunneling probability at energy (csi = E/Vf)
      Ptun = 1.d0 - (csh_amb + csh_d)/(csh_apb + csh_d)

      end subroutine
     !==========================================================================



     !========================================================================
      subroutine simpson_integrate(ndim,nmin,nmax,step,func,res)
     !========================================================================
     !     Written by Peter Szabo, 2014
     !========================================================================
     !     Integration a numerical data series with Simpson's 1/3-rule
     !     ndim --> lenght of func array
     !     nmin --> lower limit of integration
     !     nmax --> upper limit of integration
     !     step --> equidistant step in variable
     !     func --> function for integration
     !     res  --> result, value of the integral
     !========================================================================
      implicit none
      integer, intent (in) :: nmin,nmax,ndim
      integer :: i,ndata
      double precision, intent (in) :: step
      double precision :: s0,s1,s2
      double precision, intent (out) :: res
      double precision, intent (in), dimension(ndim) :: func
      res=0.d0
      s0=0.d0
      s1=0.d0
      s2=0.d0
      ndata=nmax-nmin+1

      if(nmax > nmin) then
        do i=nmin+1,nmax-1,2
           s1 = s1+func(i-1)
           s0 = s0+func(i)
           s2 = s2+func(i+1)
        enddo
        res = step*(s1+4.d0*s0+s2)/3.d0

!       If n is even, add the last slice separately
        if(mod(ndata,2) .eq. 0) then
           res=res+                                                    &
           step*(5.d0*func(nmax)+8.d0*func(nmax-1)-func(nmax-2))/12.d0
        endif
      endif
      end subroutine
     !========================================================================

      end module
