!===================================================================
 module mod_thermo_functions
 implicit none

!------------------------------------------------------------------------------------------
!Here we define the class of "type_molecule" including its variables and the corresponding
!functions and subroutines acting on this class
!------------------------------------------------------------------------------------------
 type, public :: type_molecule
!------------------------------------------------------------------------------------------
    integer :: natom, ntrans, nrot, nvib, multiplicity,nlin

    real*8 :: sigma, chiral  !symmetry numbers
    real*8 :: freq_imag, Vforward, Vbackward !for Eckart-tunneling of TS
    real*8 :: masstot,Eelec, freq_scale
    real*8 :: U_ele, H_ele, S_ele, F_ele, G_ele
    real*8 :: U_tot, H_tot, S_tot, F_tot, G_tot, ZPE_tot
    real*8 :: U_vib_tot, H_vib_tot, S_vib_tot, F_vib_tot, G_vib_tot
    real*8 :: U_rot_tot, H_rot_tot, S_rot_tot, F_rot_tot, G_rot_tot
    real*8 :: U_tra_tot, H_tra_tot, S_tra_tot, F_tra_tot, G_tra_tot
    real*8 :: PF_ele, PF_tra, PF_rot, PF_vib, PF_tot


    real*8, allocatable :: mass(:),  wmass(:) 
    real*8, allocatable :: omega(:)
    real*8, allocatable :: Inertia(:) 
    real*8, allocatable :: qxyz(:,:)

   !Thermodynamical functions for each degrees of freedom
    real*8, allocatable :: U_vib(:), H_vib(:), S_vib(:), F_vib(:), G_vib(:), ZPE_mod(:)
    real*8, allocatable :: U_rot(:), H_rot(:), S_rot(:), F_rot(:), G_rot(:)
    real*8, allocatable :: U_tra(:), H_tra(:), S_tra(:), F_tra(:), G_tra(:)

    character(len=2), allocatable :: atom(:)
    character(len=250) :: molname
    character(len=4) :: moltype !moltyp = reac, prod, ts

   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !Here we connect(bound) several procedures (function and/or subroutines) to this class
   !acting on all of the objects which defined as "type_molecule"
   ! 
   !These procedures can be used in the main program where
   !any objects is defined as type_molecule:
   !type(type_molecule),allocatable :: Reactant(:)
   !
   !After this, the subroutines can be called on Reactant(i) as
   !call Reactant(i)%AllocateVariable()
   !or
   !call Reactant(i)%Eval_All_Thermo_Functions()
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    contains
    procedure, public  :: Eval_All_Thermo_Functions => Eval_All_Thermo_Functions_Procedure !renamed
    procedure, public  :: Read_Structure_and_Freq
    procedure, public  :: Print_thermo_funcs_long
    procedure, public  :: Print_thermo_funcs_short
    procedure, private :: Allocate_Variables_Procedure  
    procedure, private :: Electronic
    procedure, private :: Rotation_fullDim
    procedure, private :: Vibration_fullDim
   !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

 end type type_molecule



!Global parameters and variables
 real*8, parameter :: pi = 4.d0*atan(1.d0)
 real*8, parameter :: twopi = 2.d0*pi

![Anstrom]*c1=[bohr]
 real*8,parameter :: c1=1.0d0/0.5291772d0

![Hartree]*c2=[kcal/mol]
 real*8,parameter :: c2=627.51d0

![g/mol]*c3=[electron mass unit]
!real*8,parameter :: c3=1.d0/5.4857990d-4! 1836.15267343 !1838.6836605d0
 real*8,parameter :: c3=1836.15267343 !1838.6836605d0

![Hartree]*c4=[eV]
 real*8,parameter :: c4=27.2114

![Hartree]*c5=[cm-1]
 real*8,parameter :: c5=219474.d0

![Hartree]*c6=[kJ/mol]
 real*8,parameter :: c6=2625.5d0

![frequency in cm-1]/c7=[omega(bohr^(-1))]
 real*8,parameter :: c7=1.0e8*c1

![femto-sec]*c8=[time in au]
 real*8,parameter :: c8=41.341105

![Pascal]*c9=[atomic pressure unit]
 real*8,parameter :: c9=1.0d-13/2.9421912

![speed of light in atomic unit]
 real*8,parameter :: clight=137.035999074

! freq in cm-1 --> energy quantum in Hartree
 real*8, parameter :: c10=clight/c7*twopi

! freq in bohr3 --> cm3
 real*8,parameter :: c11 = 1.d-24/(c1*c1*c1)

!Gas constant
 real*8, parameter :: Rgas = 8.31446261815324/1000.0/c6 !kJ/mol/K-->atomic units

!in atomic units hbar = hplanck/2pi = 1, hence in atomic units:
 real*8, parameter :: hplanck = twopi

 real*8 :: freq_threshold
 real*8 :: Temp,Pressure
 integer :: nMolec, nTS
 character(len=100) :: which_unit, which_unit_S


!This are the classes we use
 type(type_molecule), allocatable :: Molec(:), TS(:)
!------------------------------------------------------------------------------------------


!****************************************************************************************
 private :: Rotation_Class_1D, Vibration_1D
 private :: Entropy_rot, Entropy_vib, Grimme_entropy_lowfreq_vib

 contains

   !Suborutine and functions working on the class

   !------------------------------------------------------------
    subroutine Allocate_Variables_Procedure(T,imod)
   !------------------------------------------------------------
    implicit none
    class(type_molecule) :: T !--> T = this object
    integer :: imod, nat,nv,nr,nt
    nat = T%natom        
    nv  = T%nvib 
    nr  = T%nrot 
    nt  = T%ntrans 

    if(imod == 1) then
       allocate(T%atom(nat), T%qxyz(nat,3),T%mass(nat), T%wmass(nt))
       allocate(T%omega(nv), T%Inertia(nr))
       allocate(T%U_vib(nv), T%H_vib(nv), T%S_vib(nv))
       allocate(T%F_vib(nv), T%G_vib(nv), T%ZPE_mod(nv))
       allocate(T%U_rot(nr), T%H_rot(nr), T%S_rot(nr))
       allocate(T%F_rot(nr), T%G_rot(nr))
       allocate(T%U_tra(nt), T%H_tra(nt), T%S_tra(nt))
       allocate(T%F_tra(nt), T%G_tra(nt))
    else
       deallocate(T%mass, T%wmass, T%omega)
       deallocate(T%Inertia, T%qxyz, T%atom)
       deallocate(T%U_vib, T%H_vib, T%S_vib)
       deallocate(T%F_vib, T%G_vib, T%ZPE_mod)
       deallocate(T%U_rot, T%H_rot, T%S_rot, T%F_rot, T%G_rot)
       deallocate(T%U_tra, T%H_tra, T%S_tra, T%F_tra, T%G_tra)
    endif

    end subroutine
   !------------------------------------------------------------

   !------------------------------------------------------------
    subroutine Read_Structure_and_Freq(T,ifile)
   !------------------------------------------------------------
    implicit none
    class(type_molecule) :: T !--> T = this object
    integer :: i,j,isave,iatom, ndegree, ifile
    character(len=250) :: fgeom
    character(len=7)   :: what
    character(len=100) :: inertia_unit
    real*8 :: atomic_mass, umass, req_diatom

    real*8, parameter  :: av=6.022045d0
    real*8, parameter  :: ph=6.626176d0
    real*8, parameter  :: sl=2.99792458d0
    real*8, parameter  :: tocmm1=1.0d2*(ph*av)/(8.d0*pi*pi*sl)*c3*c1*c1
    real*8, parameter  :: toghz=1.0d3*(ph*av)/(8.d0*pi*pi)*c3*c1*c1

   !Read blank line for the line sperator
    !read(ifile,*)

    read(ifile,*) T%molname
    read(ifile,*) T%moltype
    read(ifile,*) T%natom, T%nlin
    read(ifile,*) T%nvib, T%nrot
    T%ntrans = 3
      

    call Allocate_Variables_Procedure(T,1) 

    read(ifile,*) T%Eelec, T%multiplicity

    if(T%moltype .ne. "Atom") then

      read(ifile,*) T%sigma, T%chiral
      read(ifile,*) what 


      if(what == "xyzgeom") then
         read(ifile,*) fgeom
         open(11,file=fgeom)
         read(11,*) !blank read of natom
         read(11,*) !blank read of comment
         do iatom=1,T%natom
             read(11,*) T%atom(iatom), T%qxyz(iatom,1), T%qxyz(iatom,2), T%qxyz(iatom,3)
         enddo
         close(11)

         do i=1,T%natom
            T%mass(i) = atomic_mass(T%atom(i))
         enddo
         T%masstot=sum(T%mass)

         isave=9
         if(T%nrot == 3) then     
            open(isave,file="Principial_moments_of_inertia_analysis.dat")
            call Moment_of_Inertia(T%natom,T%atom,T%mass,T%qxyz,isave,T%Inertia)
            close(isave)

            !T%Inertia comes from Moment_of_Inertia() in Angstrom^2 * amu
            !We have to convert it to a.u.
            T%Inertia = T%Inertia*c3*c1*c1
         else
           write(6,*) "Wrong value of nrot=",T%nrot
           write(6,*) "When the keyword xyz-geometry is used then the number of rotational DOF must be nrot=3"
           write(6,*) "When 3 > nrot or nrot < 3 , then explicit rotational constants must be given in the input"
           call flush(6)
           stop
         endif


      elseif(what == "diatom") then

         if(T%nrot == 1) then
           read(ifile,*) T%atom(1),T%atom(2), req_diatom

           req_diatom = req_diatom*c1 !in bohr

           do i=1,T%natom
              T%mass(i) = atomic_mass(T%atom(i))
           enddo
           !write(6,*) "Diatom masses:"
           !write(6,*) T%atom(1), T%mass(1)
           !write(6,*) T%atom(2), T%mass(2)
           !write(6,*)

           T%masstot = sum(T%mass)

           umass = c3*T%mass(1)*T%mass(2)/T%masstot !reduced mass in atomic units

           T%Inertia(1) = umass*req_diatom*req_diatom
           
           !write(6,*) "diatomic molecule Intertia (au,cm-1):", T%Inertia(1),tocmm1/T%Inertia(1)!/c3*c1*c1 

         else
           write(6,*) "Wrong value of nrot=",T%nrot
           write(6,*) "When diatom is given then number of rotational DOF must be nrot=1"
           call flush(6)
           stop
         endif

         
      elseif(what == "inertia") then
         read(ifile,*) T%masstot
         read(ifile,*) inertia_unit
         do i=1,T%nrot
            read(ifile,*) T%Inertia(i)
         enddo

         if(trim(inertia_unit) == "amuAng2") then

            T%Inertia = T%Inertia*c3*c1*c1

         elseif(trim(inertia_unit) == "cm-1") then

            do j=1,T%nrot            
               T%Inertia(j) = tocmm1/T%Inertia(j)
            enddo

         elseif(trim(inertia_unit) == "GHz") then

            do j=1,T%nrot
               T%Inertia(i) = toghz/T%Inertia(j)
            enddo

         else
            stop "!!!Error!!!  Wrong Inertia unit. Only: amuAng2 or cm-1 or GHz"
         endif
      else
         stop "Wrong keyword for rotation in input. It should be either xyzgeom or inertia"
      endif



      if(T%moltype == 'TS') then
         !read(ifile,*) T%freq_imag, T%Vforward, T%Vbackward
      endif
      read(ifile,*) T%freq_scale
      do i=1,T%nvib
         read(ifile,*) T%omega(i)
      enddo

     !converting freq[cm-1] to energy quantum: hbar*omega in Hatree
      T%omega = T%omega*c10   !*clight/c7*twopi
    endif !no Atom 

   !Read blank line for the line sperator
    !read(ifile,*)

    T%masstot = T%masstot*c3 !in atomic unit
 



    write(6,"(/a)") "================================================================="
    write(6,"(a11,a100)") "Molecule:  ", T%molname
    write(6,"(a)")  "-----------------------------------------------------------------"
    write(6,"(a13, f12.3,1x,a3/)") "Temperature: ", Temp, " K "
    write(6,"(a24,i5)") "# of rotational DOF    : ", T%nrot
    write(6,"(a24,i5)") "# of vibrational DOF   : ", T%nvib
    write(6,"(a24,i5/)") "# of translational DOF : ", T%ntrans
    write(6,"(a24,i5)") "Electronic degeneracy  : ", T%multiplicity
    write(6,"(a24,f5.1/)") "Rot Symmetry number    : ", T%sigma

    call flush(ifile)
    end subroutine
   !------------------------------------------------------------



   !----------------------------------------------------------------
    subroutine Eval_All_Thermo_Functions_Procedure(T)
   !----------------------------------------------------------------
    implicit none
    class(type_molecule) :: T !--> T = this object
    real*8 :: Ezero, RT

    call Electronic(T)
    call Translation_fullDim(T)

    if(T%moltype == "Atom") then
       T%U_rot_tot = 0.d0
       T%U_vib_tot = 0.d0

       T%H_rot_tot = 0.d0
       T%H_vib_tot = 0.d0

       T%F_rot_tot = 0.d0
       T%F_vib_tot = 0.d0

       T%G_rot_tot = 0.d0
       T%G_vib_tot = 0.d0

       T%S_rot_tot = 0.d0
       T%S_vib_tot = 0.d0

       T%ZPE_tot = 0.d0
    else
      call Rotation_fullDim(T)
      call Vibration_fullDim(T) 
    endif

    Ezero = T%Eelec + T%ZPE_tot

    T%U_tot = Ezero + T%U_ele + T%U_tra_tot + T%U_rot_tot + T%U_vib_tot
    T%H_tot = Ezero + T%H_ele + T%H_tra_tot + T%H_rot_tot + T%H_vib_tot 
    T%F_tot = Ezero + T%F_ele + T%F_tra_tot + T%F_rot_tot + T%F_vib_tot
    T%G_tot = Ezero + T%G_ele + T%G_tra_tot + T%G_rot_tot + T%G_vib_tot

    T%S_tot = T%S_tra_tot + T%S_rot_tot + T%S_vib_tot
   

    RT = Rgas*Temp
    T%PF_ele = dexp(-T%F_ele/RT) 
    T%PF_tra = dexp(-T%F_tra_tot/RT) 
    T%PF_rot = dexp(-T%F_rot_tot/RT) 
    T%PF_vib = dexp(-(T%F_vib_tot+Ezero)/RT) 
    !T%PF_tot = dexp(-T%F_tot/RT) 
    T%PF_tot = dexp(-(T%F_tot - Ezero) /RT) 

    end subroutine
   !-----------------------------------------------------------------


   !----------------------------------------------------------------
    subroutine Translation_fullDim(T)
   !----------------------------------------------------------------
    implicit none
    class(type_molecule) :: T 
    integer :: i
    real*8 :: RT, Vol,lam
    real*8 :: lambda, Utra, Htra, Stra, Ftra, Gtra

    lam = sqrt(twopi*T%masstot*Rgas*Temp/(hplanck**2))
    lam = lam*lam*lam

    RT = Rgas*Temp

    Vol = RT/Pressure  

   !Each of the is defined for a 3D object
    T%U_tra_tot = 1.5d0*RT  
    T%H_tra_tot = 2.5d0*RT 
    T%S_tra_tot = Rgas*log(lam*Vol) + 2.5d0*Rgas
    T%F_tra_tot = -RT*log(lam*Vol) - RT 
    T%G_tra_tot = -RT*log(lam*Vol) 

    end subroutine
   !-----------------------------------------------------------------


   !----------------------------------------------------------------
    subroutine Rotation_fullDim(T)
   !----------------------------------------------------------------
    implicit none
    class(type_molecule) :: T 
    integer :: i
    real*8 :: Urot, Hrot, Srot, Frot, Grot
    real*8, parameter  :: av=6.022045d0
    real*8, parameter  :: ph=6.626176d0
    real*8, parameter  :: sl=2.99792458d0
    real*8, parameter  :: tocmm1=1.0d2*(ph*av)/(8.d0*pi*pi*sl)*c3*c1*c1
    real*8, parameter  :: toghz=1.0d3*(ph*av)/(8.d0*pi*pi)*c3*c1*c1
    integer :: call_first_rot = 1
    real*8 :: rot_temp

   !if(call_first_rot == 1) then
       !write(6,"(/a)") "================================================================="
       !write(6,"(a10,a)") "Molecule: ", T%molname
       !write(6,"(a13, f12.3,1x,a3/)") "Temperature: ", Temp, " K " 
       !write(6,"(a24,i5)") "# of rotational DOF    : ", T%nvib
       !write(6,"(a24,i5)") "# of vibrational DOF   : ", T%nrot
       !write(6,"(a24,i5/)") "# of translational DOF : ", T%ntrans
       !write(6,"(a24,i5)") "Electronic degeneracy  : ", T%multiplicity
       !write(6,"(a24,f5.1)") "Rot Symmetry number    : ", T%sigma
       write(6,"(/a)") "Rotational temperatures and constants:"
       write(6,"(a5,a12,3a16)") "     ", "  T[K]  ", "Brot[cm-1]", "Brot[GHz]" , "Irot[amuAng2]"
       !call_first_rot = 99
    !endif

   !In case of Atom:
    do i=1,T%nrot 


      !Carefully, it works only for odd dimensions (1D, 3D, 5D...) 
      !Check some unimolecular book (old Holbrook Pilling) of the missing
      !Gamma function factor 

       rot_temp = hplanck*hplanck/(8.d0*pi*pi*T%Inertia(i)*Rgas)

       write(6,"(i5,f12.3,3f16.5)") i, rot_temp, tocmm1/T%Inertia(i), toghz/T%Inertia(i),T%Inertia(i)/(c3*c1*c1)

       call Rotation_Class_1D(T%Inertia(i),T%sigma,Urot,Hrot,Srot,Frot,Grot)

       T%U_rot(i) = Urot
       T%H_rot(i) = Hrot
       T%S_rot(i) = Srot
       T%F_rot(i) = Frot
       T%G_rot(i) = Grot
    enddo

   !From: 
   !A Simple Method to Estimate Entropy and Free
   !Energy of Atmospheric Gases from Their Action                        
   !Ivan Kennedy, Harold Geering, Michael Rose, Angus Crossan Eq.19:    
    T%U_rot_tot = sum(T%U_rot)
    T%H_rot_tot = sum(T%H_rot)
    T%S_rot_tot = sum(T%S_rot) + Rgas*log(sqrt(pi)/T%sigma) !symn_num term added on top of the classical rotor
    T%F_rot_tot = sum(T%F_rot) - Temp*Rgas*log(sqrt(pi)/T%sigma)
    T%G_rot_tot = sum(T%G_rot) - Temp*Rgas*log(sqrt(pi)/T%sigma)

    end subroutine
   !-----------------------------------------------------------------


   !----------------------------------------------------------------
    subroutine Vibration_fullDim(T)
   !----------------------------------------------------------------
    implicit none
    class(type_molecule) :: T 
    integer :: i
    real*8 :: Uvib, Hvib, Svib, Fvib, Gvib
    integer :: call_first_vib = 1
    real*8 :: vib_temp

    !if(call_first_vib == 1) then
       write(6,"(/a)") "Vibrational temperatures and Freq[cm-1]:"
       !write(6,"(a5,a12,a16)") "     ","  T[K]  ", "Freq[cm-1]"
       call_first_vib = 99
    !endif

    do i=1,T%nvib
      vib_temp = T%omega(i)/Rgas
      write(6,"(i5, f12.3,f16.2)") i, vib_temp, T%omega(i)/c10 
      call Vibration_1D(T%omega(i), Uvib, Hvib, Svib, Fvib, Gvib)

      T%U_vib(i) = Uvib
      T%H_vib(i) = Hvib
      T%S_vib(i) = Svib
      T%F_vib(i) = Fvib
      T%G_vib(i) = Gvib
      T%ZPE_mod(i) = 0.5d0*T%omega(i)
    enddo

    T%U_vib_tot = sum(T%U_vib)
    T%H_vib_tot = sum(T%H_vib)
    T%S_vib_tot = sum(T%S_vib)
    T%F_vib_tot = sum(T%F_vib)
    T%G_vib_tot = sum(T%G_vib)

    T%ZPE_tot = sum(T%ZPE_mod)

    end subroutine
   !-----------------------------------------------------------------


   !----------------------------------------------------------------
    subroutine Electronic(T)
   !----------------------------------------------------------------
    implicit none
    class(type_molecule) :: T 
 
   !When only one electronic state is available
   !at a given Temp, then the partition function
   !approx Qele = multiplicity

    T%U_ele = 0.d0 

    T%H_ele = T%U_ele

    T%S_ele = Rgas*log(float(T%multiplicity))

    T%F_ele = T%U_ele - Temp*T%S_ele

    T%G_ele = T%F_ele 

    end subroutine
   !-----------------------------------------------------------------



   !----------------------------------------------------------------------
    subroutine Rotation_Class_1D(Inert,sigma,Urot, Hrot, Srot, Frot, Grot)
   !----------------------------------------------------------------------
   ! Here sigma, the symmetry number is not used eventually
   ! It's introduced where (right after) this routine called
   !----------------------------------------------------------------------
    implicit none
    real*8, intent(in)  :: Inert,sigma
    real*8, intent(out) :: Urot, Hrot, Srot, Frot, Grot
    real*8 :: theta
    
    theta = sqrt(8.d0*pi*pi*Inert*Rgas*Temp/hplanck/hplanck)

    Urot = 0.5d0*Rgas*Temp
   
    Hrot = Urot 
   
    Srot = Rgas*log(theta) + 0.5d0*Rgas
   
    Frot = Urot - Temp*Srot !-Rgas*Temp*log(Theta/sigma)
   
    Grot = Hrot - Temp*Srot
    
    end subroutine
   !----------------------------------------------------------------------



   !------------------------------------------------------------
    real*8 function Entropy_rot(mu) result(Sr)
   !------------------------------------------------------------
    implicit none
    real*8, intent(in) :: mu
    real*8 :: dum
   
   !1D rotational entropy: 
    dum = log(sqrt(8.d0*pi**3 * mu * Rgas*Temp))/hplanck/hplanck

    Sr = Rgas*(0.5d0 + dum)

    end function        
   !------------------------------------------------------------


   !------------------------------------------------------------
    real*8 function Entropy_vib(omega) result(Sv)
   !------------------------------------------------------------
    implicit none
    real*8, intent(in) :: omega
    real*8 :: theta, dum

    theta = omega/Rgas/Temp  !in atomic units

    dum = theta/(exp(theta)-1.d0) - log(1.d0 - exp(-theta))

    Sv = Rgas*dum

    end function
   !------------------------------------------------------------


   !------------------------------------------------------------
    real*8 function Grimme_entropy_lowfreq_vib(ome) result(Sv)
   !------------------------------------------------------------
   ! Grimme's approximation for low freq modes                                                  
   ! which are replaced with a rigid rotor mode
   ! Ref: Stefan Grimme, Vol. 18, (2012), pp. 9955-9964
   ! https://doi.org/10.1002/chem.201200497
   !------------------------------------------------------------
    implicit none
    real*8, intent(in) :: ome
    real*8 :: Wgt, mu, mu_prime,dum
    real*8, parameter :: ome_zero = 100.d0*c10 !100 cm-1 ~ 0.5kT at room temp
    real*8, parameter :: alpha = 4.d0
    real*8, parameter :: av=6.022045d0
   !real*8, parameter :: Bav = 1.d-4  *E-40 !kg*m2
   !real*8, parameter :: Bav = 0.1    *E-40 ! g*m2
   !real*8, parameter :: Bav = 1000   *E-40 ! g*cm2
   !real*8, parameter :: Bav = 1000*av/10 ! amu*Ansgtr2
    real*8, parameter :: Bav = 100**av*c3*c1*c1 ! amu*Ansgtr2


   !e1=eig(i)/av*10.d0

   !Effective inertia of rotation with same period
   !as the low-frew vibration mode
    mu = hplanck/(8.d0*pi*pi*ome)

    !mu_prime = mu*Bav / (mu + Bav)
    mu_prime = mu*Bav / (mu + Bav)

    write(88,*) ome/c10, "mu = ", mu, "Bav = ", Bav, "mu_prime = ", mu_prime

   !Weight function to smoothly switch between vib --> rot
    dum = (ome_zero/ome)**alpha
    Wgt = 1.d0 / (1.d0 + dum)

    Sv = Wgt*Entropy_vib(ome) + (1.d0-Wgt)*Entropy_rot(mu_prime)
    end function                                                             
   !------------------------------------------------------------

                                                                           
                                                                           
   !------------------------------------------------------------
    subroutine Vibration_1D(omega, Uvib, Hvib, Svib, Fvib, Gvib)                                  
   !------------------------------------------------------------
    implicit none                                                             
    real*8, intent(in) :: omega
    real*8, intent(out) :: Uvib, Hvib, Svib, Fvib, Gvib
    !real*8 :: Entropy_vib, Entropy_rot                                        
    real*8 :: theta, B, Wgt, mu                                               

    theta = omega/Rgas                                                          
                                                                           
    Uvib = Rgas*theta/(exp(theta/Temp)-1.d0) 

    Hvib = Uvib                                                             

    if(omega > freq_threshold) then                                                   
      Svib = Entropy_vib(omega)                                              
    else
     !Grimme's approximation for low freq vibration modes                                                  
      Svib = Grimme_entropy_lowfreq_vib(omega)                               
    endif
                                                                           
    Fvib = Uvib - Temp*Svib !Rgas*Temp*log(1.d0 - exp(-theta/Temp))

    Gvib = Hvib - Temp*Svib

    end subroutine
   !------------------------------------------------------------




   !------------------------------------------------------------
    subroutine Print_thermo_funcs_long(T)
   !------------------------------------------------------------
    implicit none
    class(type_molecule) :: T 
    character(len=100) :: form0,form1,form2
    real*8 :: cx


    if(which_unit == "kcal/mol") then
      cx = c2 
      which_unit_S = "cal/mol/K"
    elseif(which_unit == "kJ/mol") then
      cx = c6 
      which_unit_S = "J/mol/K"
    endif       
    which_unit = trim(which_unit)


    form0 = "(2x,a,f15.8,2x,a)"
    form1 = "(2x,a,f15.8,2x,a,f12.3,2x,a)"
    form2 = "(2x,a,f15.8,2x,a,f12.3,2x,a8,10x,a8,f12.3,2x,a12)"
   !write(6,"(/a)") "===================================================="
   !write(6,"(a,i5,5x,a)") " Molecule ", i, T%molname
   !write(6,"(a,f10.2)") "Temperature: ", Temp 
    write(6,"(a)") "================================================================"

    write(6,"(/a)")"Electronic and zero-point energy:"
    write(6,form2) "Eele    = ", T%Eelec,     "Eh  "
    write(6,form1) "ZPE     = ", T%ZPE_tot,   "Eh  ", T%ZPE_tot*cx,    which_unit 

    write(6,"(//a)") "Internal energy thermal contributions:"
    write(6,form1) "U_ele   = ", T%U_ele, "Eh  ", T%U_ele*cx,  which_unit 
    write(6,form1) "U_vib   = ", T%U_vib_tot, "Eh  ", T%U_vib_tot*cx,  which_unit 
    write(6,form1) "U_rot   = ", T%U_rot_tot, "Eh  ", T%U_rot_tot*cx,  which_unit 
    write(6,form1) "U_tra   = ", T%U_tra_tot, "Eh  ", T%U_tra_tot*cx,  which_unit 
    write(6,"(2x,a)") "-----------------------------------------------------"
    write(6,form1) "U_therm = ", T%U_tot-T%Eelec, "Eh  ",  (T%U_tot-T%Eelec)*cx, which_unit
    write(6,"(2x,a)") "-----------------------------------------------------"
    write(6,"(2x,a)") "U_therm = ZPE + U_ele + U_vib + U_rot + U_tra"
    write(6,"(2x,a)") "-----------------------------------------------------"

    write(6,"(//a)") "Enthalpy thermal contributions:"
    write(6,form1) "H_ele   = ", T%H_ele, "Eh  ", T%H_ele*cx, which_unit
    write(6,form1) "H_vib   = ", T%H_vib_tot, "Eh  ", T%H_vib_tot*cx, which_unit
    write(6,form1) "H_rot   = ", T%H_rot_tot, "Eh  ", T%H_rot_tot*cx, which_unit
    write(6,form1) "H_tra   = ", T%H_tra_tot, "Eh  ", T%H_tra_tot*cx, which_unit
    write(6,form1) "kB*T    = ", Rgas*Temp,          "Eh  ", Rgas*Temp*cx         , which_unit
    write(6,"(2x,a)") "-----------------------------------------------------"
    write(6,form1) "H_therm = ", T%H_tot-T%Eelec, "Eh  ",  (T%H_tot-T%Eelec)*cx, which_unit
    write(6,"(2x,a)") "-----------------------------------------------------"
    write(6,"(2x,a)") "H_therm = ZPE + H_ele + H_vib + H_rot + H_tra"
    write(6,"(2x,a)") "-----------------------------------------------------"


    write(6,"(//a)") "Entropy contributions in energy unit (S*T)               and     in natural unit:"
    write(6,form2) "S_ele*T = ", T%S_ele*Temp, "Eh  ",              &
       &  T%S_ele*cx*Temp, which_unit,                             &
       &  "S_ele = ", T%S_ele*cx*1000.d0, which_unit_S

    write(6,form2) "S_vib*T = ", T%S_vib_tot*Temp, "Eh  ",              &
       &  T%S_vib_tot*cx*Temp, which_unit,                             &
       &  "S_vib = ", T%S_vib_tot*cx*1000.d0, which_unit_S

    write(6,form2) "S_rot*T = ", T%S_rot_tot*Temp, "Eh  ",              &
       &  T%S_rot_tot*cx*Temp, which_unit,                             &
       &  "S_rot = ", T%S_rot_tot*cx*1000.d0, which_unit_S

    write(6,form2) "S_tra*T = ", T%S_tra_tot*Temp, "Eh  ",              &
       &  T%S_tra_tot*cx*Temp, which_unit,                             &
       &  "S_tra = ", T%S_tra_tot*cx*1000.d0, which_unit_S

    write(6,"(2x,a,10x,a)") "-----------------------------------------------------", &
       &           "-------------------------------"

    write(6,form2) "S_tot*T = ", T%S_tot*Temp, "Eh  ",              &
       &  T%S_tot*cx*Temp, which_unit,                             &
       &  "S_tot = ", T%S_tot*cx*1000.d0, which_unit_S
    write(6,"(2x,a,10x,a)") "-----------------------------------------------------", &
       &           "-------------------------------"

    write(6,"(2x,a)") "S_tot = S_ele + S_vib + S_rot + S_tra"
    write(6,"(2x,a)") "-----------------------------------------------------"



    write(6,"(//a)") "Gibbs Free Energy thermal contributions:"
    write(6,form1) "G_ele   = ", T%G_ele, "Eh  ", T%G_ele*cx, which_unit 
    write(6,form1) "G_vib   = ", T%G_vib_tot, "Eh  ", T%G_vib_tot*cx, which_unit 
    write(6,form1) "G_rot   = ", T%G_rot_tot, "Eh  ", T%G_rot_tot*cx, which_unit 
    write(6,form1) "G_tra   = ", T%G_tra_tot, "Eh  ", T%G_tra_tot*cx, which_unit 
    write(6,"(2x,a)") "-----------------------------------------------------"
    write(6,form1) "G_therm = ", T%G_tot-T%Eelec, "Eh  ",  (T%G_tot-T%Eelec)*cx, which_unit
    write(6,"(2x,a)") "-----------------------------------------------------"
    write(6,"(2x,a)") "G_therm = ZPE + G_ele + G_vib + G_rot + G_tra"
    write(6,"(2x,a)") "-----------------------------------------------------"

    write(6,"(//a)") "Helmholtz Free Energy thermal contributions:"
    write(6,form1) "F_ele   = ", T%F_ele, "Eh  ", T%F_ele*cx, which_unit
    write(6,form1) "F_vib   = ", T%F_vib_tot, "Eh  ", T%F_vib_tot*cx, which_unit
    write(6,form1) "F_rot   = ", T%F_rot_tot, "Eh  ", T%F_rot_tot*cx, which_unit
    write(6,form1) "F_tra   = ", T%F_tra_tot, "Eh  ", T%F_tra_tot*cx, which_unit
    write(6,"(2x,a)") "-----------------------------------------------------"
    write(6,form1) "F_therm = ", T%F_tot-T%Eelec, "Eh  ",  (T%F_tot-T%Eelec)*cx, which_unit
    write(6,"(2x,a)") "-----------------------------------------------------"
    write(6,"(2x,a)") "F_therm = ZPE + F_ele + F_vib + F_rot + F_tra"
    write(6,"(2x,a)") "-----------------------------------------------------"


    write(6,"(//a)")"Electronic and ZPE energy:"
    write(6,form0) "Eele         = ", T%Eelec,         "Eh"
    write(6,form0) "ZPE          = ", T%ZPE_tot,      "Eh"
    write(6,form0) "Eele+ZPE     = ", T%Eelec+T%ZPE_tot,      "Eh"
    write(6,"(/a)")"Energy functions including the electronic and ZPE energy:"
    write(6,form0) "U_tot        = ", T%U_tot,      "Eh"
    write(6,form0) "H_tot        = ", T%H_tot,      "Eh"
    write(6,form0) "F_tot        = ", T%F_tot,      "Eh"
    write(6,form0) "G_tot        = ", T%G_tot,      "Eh"
    write(6,form0) "T*S_tot      = ", T%S_tot*Temp,      "Eh"
!    write(6,form0) "G_tot'       = ", T%H_tot - T%S_tot*Temp, "Eh"
!    write(6,form0) "G_tot'-Elec  = ", (T%H_tot - T%S_tot*Temp - T%Eelec)*cx, which_unit

    write(6,*)
    write(6,"(1x,a)") "--------------------------------------------------------------"
    write(6,"(a17, 4x, a40)") "Partfunc      = ", "exp(-(F_ele + F_vib + F_rot + F_tra)/RT)"
    write(6,"(a17, 4x, a40)") "Partfunc      = ", "exp(-(F_tot - Eele - ZPE)/RT)"
    write(6,"(1x,a)") "--------------------------------------------------------------"
    write(6,"(a17, e15.5)") "Partfunc      = ", T%PF_tot
    write(6,"(1x,a/)") "--------------------------------------------------------------"
    write(6,"(a)") "==================== End of Molecule =========================="
    write(6,*)
    write(6,*)

    end subroutine


   !------------------------------------------------------------
    subroutine Print_thermo_funcs_short(T)
   !------------------------------------------------------------
    implicit none
    class(type_molecule) :: T 
    character(len=100) :: form0,form1,form2
    real*8 :: cx


    if(which_unit == "kcal/mol") then
      cx = c2 
      which_unit_S = "cal/mol/K"
    elseif(which_unit == "kJ/mol") then
      cx = c6 
      which_unit_S = "J/mol/K"
    endif       
    which_unit = trim(which_unit)


    form0 = "(2x,a,f15.8,2x,a)"
    form1 = "(2x,a,f15.8,2x,a,f12.3,2x,a)"
    form2 = "(2x,a,f15.8,2x,a,f12.3,2x,a8,10x,a8,f12.3,2x,a12)"

    write(6,"(//a)")"Electronic and ZPE energy:"
    write(6,form0) "Eele         = ", T%Eelec,         "Eh"
    write(6,form0) "ZPE          = ", T%ZPE_tot,      "Eh"
    write(6,form0) "Eele+ZPE     = ", T%Eelec+T%ZPE_tot,      "Eh"
    write(6,"(/a)")"Energy functions including the electronic and ZPE energy:"
    write(6,form0) "U_tot        = ", T%U_tot,      "Eh"
    write(6,form0) "H_tot        = ", T%H_tot,      "Eh"
    write(6,form0) "F_tot        = ", T%F_tot,      "Eh"
    write(6,form0) "G_tot        = ", T%G_tot,      "Eh"
    write(6,form0) "T*S_tot      = ", T%S_tot*Temp,      "Eh"
!    write(6,form0) "G_tot'       = ", T%H_tot - T%S_tot*Temp, "Eh"
!    write(6,form0) "G_tot'-Elec  = ", (T%H_tot - T%S_tot*Temp - T%Eelec)*cx, which_unit

    write(6,*)
    write(6,"(1x,a)") "--------------------------------------------------------------"
    write(6,"(a17, 4x, a40)") "Partfunc      = ", "exp(-(F_ele + F_vib + F_rot + F_tra)/RT)"
    write(6,"(a17, 4x, a40)") "Partfunc      = ", "exp(-(F_tot - Eele - ZPE)/RT)"
    write(6,"(1x,a)") "--------------------------------------------------------------"
    write(6,"(a17, e15.5)") "Partfunc      = ", T%PF_tot
    write(6,"(1x,a/)") "--------------------------------------------------------------"
    write(6,"(a)") "==================== End of Molecule =========================="
    write(6,*)
    write(6,*)

    end subroutine


 end module mod_thermo_functions
!=========================================================================


 program test_thermo
 use mod_thermo_functions
 implicit none
 integer :: i,iread,nEquilibria,nReaction
 character(len=100) :: form1
 integer, allocatable :: EqBetween(:,:), Reac_TS(:,:)
 real*8, allocatable :: rate(:), Keq(:)
 integer :: mola, molb, molTS
 real*8 :: RT, kappa, prefac, Vmolar, dG, G_reactant

 
                                                                       
 iread = 5

 read(iread,*) nMolec, Temp, Pressure
 read(iread,*) which_unit
 read(iread,*) freq_threshold
 freq_threshold = freq_threshold*c10 !cm-1 is converted into energy unit

 read(iread,*) nEquilibria, nReaction
 allocate(EqBetween(nEquilibria,2))
 allocate(Reac_TS(nReaction,3))

 do i=1,nEquilibria
   read(iread,*) EqBetween(i,1), EqBetween(i,2)
 enddo

 do i=1,nReaction
   read(iread,*) Reac_TS(i,1), Reac_TS(i,2),Reac_TS(i,3)
 enddo


 Pressure = Pressure*c9
 allocate(Molec(nMolec))

 do i=1,nMolec
    call Molec(i)%Read_Structure_and_Freq(iread)
    call Molec(i)%Eval_All_Thermo_Functions()
    call Molec(i)%Print_thermo_funcs_short()
 enddo


 allocate(rate(nReaction))
 allocate(Keq(nEquilibria))

 RT = Rgas*Temp
 do i=1,nEquilibria
    mola = EqBetween(i,1)
    molb = EqBetween(i,2)

    dG = Molec(mola)%G_tot - Molec(molb)%G_tot
    Keq(i) = exp(-dG/RT)
    write(6,*) "Keq = ", Keq(i), " for the equlibrium between:", Molec(mola)%molname, " and ", Molec(molb)%molname
 enddo

!kB = 1.380649
!hplanck = 6.62607015
 Vmolar = RT/Pressure*c11
 kappa = 1.d0

 prefac = Temp*1.380649/6.62607015*1.0d11*Vmolar

 do i=1,nReaction
    mola = EqBetween(i,1)
    molb = EqBetween(i,2)
    molTS = EqBetween(i,3)

    G_reactant = Molec(mola)%G_tot + Molec(molb)%G_tot
    dG = Molec(molTS)%G_tot - G_reactant


    rate(i) = kappa*prefac*exp(-dG/RT)

    write(6,*) "rate[cm3*mol-1*s-1] = ", rate(i)!, "Reaction between:", Molec(mola)%molname, " and ", Molec(molb)%molname, " through TS ", Molec(molTS)%molname
 enddo




 
  end program






