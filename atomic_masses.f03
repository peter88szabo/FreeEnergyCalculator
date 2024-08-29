! program test
! implicit none
! character(len=2) :: atom
! real*8 :: atomic_mass


! write(6,*) "give an element"
! read(5,*) atom

! write(6,"(a8,a2,a5,f12.5)") "mass of ", atom, " is: ", atomic_mass(atom)

! end program

!==========================================================================
 real*8 function atomic_mass(atom) result(mass)
!==========================================================================
! relative atomic weights of elements
!==========================================================================
 implicit none
 character(len=2), intent(in) :: atom


 select case (trim(atom))
 case("H")
     mass = 1.00784d0
 case("D")
     mass = 2.01410177d0
 case("T")
     mass = 3.10604928d0
 case("He")
     mass = 4.002600d0
 case("Li")
     mass = 6.940000d0
 case("Be")
     mass = 9.012180d0
 case("B")
     mass = 10.81000d0
 case("C") 
     mass = 12.0000d0
 case("N") 
     mass = 14.007d0
 case("O") 
     mass = 15.9994d0 
 case("F") 
     mass = 18.99840d0
 case("Ne") 
     mass = 20.17900d0 
 case("Na") 
     mass = 22.9898d0 
 case("Mg") 
     mass = 24.30500d0
 case("Al") 
     mass = 26.98154d0
 case("Si") 
     mass = 28.08550d0
 case("P") 
     mass = 30.97370d0  
 case("S") 
     mass = 32.0600d0 
 case("Cl") 
     mass = 35.45300d0
 case("Ar") 
     mass = 39.94800d0
 case("K") 
     mass = 39.09830d0
 case("Ca") 
     mass = 40.08000d0
 case("Sc") 
     mass = 44.9559d0 
 case("Ti") 
     mass = 47.90000d0 
 case("V") 
     mass = 50.94150d0
 case("Cr") 
     mass = 51.99600d0 
 case("Mn") 
     mass = 54.93800d0
 case("Fe") 
     mass = 55.8470d0 
 case("Co") 
     mass = 58.93320d0 
 case("Ni") 
     mass = 58.6934d0 
 case("Cu") 
     mass = 63.54600d0 
 case("Zn") 
     mass = 65.38000d0 
 case("Ga") 
     mass = 69.7350d0
 case("Ge") 
     mass = 72.59000d0 
 case("As") 
     mass = 74.92160d0 
 case("Se") 
     mass = 78.96000d0  
 case("Br") 
     mass = 79.90400d0
 case("Kr") 
     mass = 83.7980d0 
 case("Rb") 
     mass = 85.46780d0 
 case("Sr") 
     mass = 87.62000d0  
 case("Y") 
     mass = 88.90590d0 
 case("Zr") 
     mass = 91.22000d0 
 case("Nb") 
     mass = 92.9064d0 
 case("Mo") 
     mass = 95.94000d0 
 case("Tc") 
     mass = 98.90620d0  
 case("Ru") 
     mass = 101.0700d0 
 case("Rh") 
     mass = 102.9065d0 
 case("Pd") 
     mass = 106.400d0 
 case("Ag") 
     mass = 107.8680d0 
 case("Cd") 
     mass = 112.4100d0 
 case("In") 
     mass = 114.8200d0 
 case("Sn") 
     mass = 118.6900d0 
 case("Sb")
     mass = 121.750d0 
 case("Te")
     mass = 127.6000d0 
 case("I")
     mass = 126.9045d0 
 case("Xe")
     mass = 131.3000d0 
 case("Cs")
     mass = 132.9054d0 
 case("Ba") 
     mass = 137.327d0 
 case("Hf")
     mass = 178.4900d0 
 case("Ta")
     mass = 180.9479d0 
 case("W")
     mass = 183.8500d0 
 case("Re")
     mass = 186.207d0 
 case("Os")
     mass = 190.2300d0 
 case("Ir")
     mass = 192.2200d0 
 case("Pt")
     mass = 195.08d0 
 case("Au")
     mass = 196.9665d0 
 case("Hg")
     mass = 200.590d0 
 case("Tl")
     mass = 204.3700d0    
 case("Pb")
     mass = 207.2000d0
 case("Bi")
     mass = 208.9804d0
 case("Po")
     mass = 209.d0 
 case("At")
     mass = 210.d0
 case("Rn")
     mass = 222.d0
 case("Fr")
     mass = 223.d0
 case("Ra")
     mass = 226.d0
!F-elements
 case("La")
     mass = 138.91d0
 case("Ce")
     mass = 140.12
 case("Pr")
     mass = 140.91
 case("Nd")
     mass = 144.24 
 case("Pm")
     mass = 145.d0
 case("Sm")
     mass = 150.36d0 
 case("Eu")
     mass = 151.25d0 
 case("Gd")
     mass = 157.25d0 
 case("Tb")
     mass = 158.93d0 
 case("Dy")
     mass = 162.5d0 
 case("Ho")
     mass = 164.93d0
 case("Er")
     mass = 167.26d0 
 case("Tm")
     mass = 168.93d0 
 case("Yb")
     mass = 173.05d0 
 case("Lu")
     mass = 174.97d0 
 case("Ac")
     mass = 227.d0 
 case("Th")
     mass = 232.04d0 
 case("Pa")
     mass = 231.04 
 case("U")
     mass = 238.03 
 case("Np")
     mass = 237.d0 
 case("Pu")
     mass = 244.d0 
 case("Am")
     mass = 243.d0
 case("Cm")
     mass = 247.d0 
 case("Bk")
     mass = 247.d0 
 case("Cf")
     mass = 251.d0 
 case("Es")
     mass = 252.d0
 case("Fm")
     mass = 257.d0
 case("Md")
     mass = 258.d0
 case("No")
     mass = 259.d0
 case("Lr")
     mass = 266.d0

!New unstable elements which are not useful for chemistry
 case("Rf")
     mass = 267.d0
 case("Db")
     mass = 268.d0
 case("Sg")
     mass = 269.d0
 case("Bh")
     mass = 270.d0
 case("Hs")
     mass = 277.d0
 case("Mt")
     mass = 278.d0
 case("Ds")
     mass = 281.d0
 case("Rg")
     mass = 282.d0
 case("Cn")
     mass = 285.d0
 case("Nh")
     mass = 286.d0
 case("Fl")
     mass = 289.d0
 case("Mc")
     mass = 290.d0
 case("Lv")
     mass = 293.d0
 case("Ts")
     mass = 294.d0
 case("Og")
     mass = 294.d0

 case default
    write(6,"(/a18)") "!!!!Error!!!!"
    write(6,"(a18,2x,a3)") "Invalid atom name: ",atom
    call flush(6)
    stop 
 end select

 end function atomic_mass
