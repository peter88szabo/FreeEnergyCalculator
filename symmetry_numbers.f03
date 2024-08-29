! program test
! implicit none
! character(len=2) :: atom
! real*8 :: atomic_mass


! write(6,*) "give an element"
! read(5,*) atom

! write(6,"(a8,a2,a5,f12.5)") "mass of ", atom, " is: ", atomic_mass(atom)

! end program

!==========================================================================
 real*8 function symmetry_number(irrep) result(sym)
!==========================================================================
! relative atomic weights of elements
!==========================================================================
 implicit none
 character(len=5), intent(in) :: irrep

 select case (trim(atom))
 case("C1")
     symm = 1.d0
 case("Cs")
     symm = 1.d0
 case("Ci")
     symm = 1.d0
 case("C2")
     symm = 2.d0
 case("C3")
     symm = 3.d0
 case("C4")
     symm = 4.d0
 case("C5")
     symm = 5.d0
 case("C6")
     symm = 6.d0
 case("C7")
     symm = 7.d0
 case("C7")
     symm = 8.d0
 case("D2")
     symm = 4.d0
 case("D3")
     symm = 6.d0
 case("D4")
     symm = 8.d0
 case("D5")
     symm = 10.d0
 case("D6")
     symm = 12.d0
 case("D7")
     symm = 14.d0
 case("D8")
     symm = 16.d0
 case("C2v")
     symm = 2.d0
 case("C3v")
     symm = 3.d0
 case("C4v")
     symm = 4.d0
 case("C5v")
     symm = 5.d0
 case("C6v")
     symm = 6.d0
 case("C7v")
     symm = 7.d0
 case("C8v")
     symm = 8.d0
 case("C2h")
     symm = 2.d0
 case("C3h")
     symm = 3.d0
 case("C4h")
     symm = 4.d0
 case("C5h")
     symm = 5.d0
 case("C6h")
     symm = 6.d0
 case("C7h")
     symm = 7.d0
 case("C8h")
     symm = 8.d0
 case("D2h")
     symm = 4.d0
 case("D3h")
     symm = 6.d0
 case("D4h")
     symm = 8.d0
 case("D5h")
     symm = 10.d0
 case("D6h")
     symm = 12.d0
 case("D7h")
     symm = 14.d0
 case("D8h")
     symm = 16.d0
 case("D2d")
     symm = 4.d0
 case("D3d")
     symm = 6.d0
 case("D4d")
     symm = 8.d0
 case("D5d")
     symm = 10.d0
 case("D6d")
     symm = 12.d0
 case("D7d")
     symm = 14.d0
 case("D8d")
     symm = 16.d0
 case("S4")
     symm = 4.d0
 case("S6")
     symm = 6.d0
 case("S4")
     symm = 4.d0
 case("S6")
     symm = 6.d0
 case("S8")
     symm = 8.d0
 case("T")
     symm = 6.d0
 case("Th")
     symm = 12.d0
 case("Td")
     symm = 12.d0
 case("O")
     symm = 12.d0
 case("Oh")
     symm = 24.d0
 case("Cinfv")
     symm = 1.d0
 case("Dinfh")
     symm = 2.d0
 case("I")
     symm = 30.d0
 case("Ih")
     symm = 60.d0
 case("Kh")
     symm = 1.d0
 case default
    write(6,"(/a18)") "!!!!Error!!!!"
    write(6,"(a30,2x,a10)") "Invalid irrep name in symmetry_number(): ",irrep
    call flush(6)
    stop 
 end select

 end function symmetry_number

      
