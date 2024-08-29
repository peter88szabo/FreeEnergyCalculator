! program calcinertia
! implicit none
! integer :: i,j,natom
! real*8 :: atomic_mass
! character(len=2), allocatable :: atom(:)
! real*8, allocatable :: mass(:), coord(:,:)
! real*8, dimension(3) :: Inertia

!Read and allocate
! read(5,*) natom
! read(5,*) 
! allocate(mass(natom),atom(natom), coord(natom,3))

! do i=1,natom
!    read(5,*) atom(i),(coord(i,j),j=1,3) 
!    mass(i) = atomic_mass(atom(i)) 
! enddo

! call Moment_of_Inertia(natom,atom,mass,coord,Inertia)

! end program

!======================================================================
 subroutine Moment_of_Inertia(natom,atom,mass,coord,isave,eig)
!======================================================================
 implicit none
 real*8, parameter  :: pi=4.d0*datan(1.d0)
 real*8, parameter  :: eps=1.d-10
 real*8, parameter  :: av=6.022045d0
 real*8, parameter  :: ph=6.626176d0
 real*8, parameter  :: sl=2.99792458d0
 real*8, parameter  :: tocmm1=1.0d2*(ph*av)/(8.d0*pi*pi*sl)
 real*8, parameter  :: toghz=1.0d3*(ph*av)/(8.d0*pi*pi) 

 integer :: i,j,k,isave, natom
 real*8, dimension(3)  :: eig
 real*8, dimension(3,3):: rotm, MomIn
 character(len=2), dimension(natom) :: atom
 real*8, dimension(3)       :: cm(3)
 real*8, dimension(natom)   :: mass
 real*8, dimension(natom,3) :: coord(natom,3)
 real*8 :: e1,e2,e3


 write(isave,'(/a,f16.6,2x,a)') 'Molecular weight:',sum(mass), "g/mol"
 write(666,'(/a,f16.6,2x,a)') 'Molecular weight:',sum(mass), "g/mol"
 write(isave,'(/a)') "Coordinates before COM-shift " 
 write(isave,"(a2, 3a14,3x,a12)") "    ","  x  ",  "  y  " , "  z  ", " mass "
 do i=1,natom
    write(isave,'(a2,3f14.5,3x,f12.5)') atom(i),(coord(i,j),j=1,3),mass(i)
    write(666,'(a2,3f14.5,3x,f12.5)') atom(i),(coord(i,j),j=1,3),mass(i)
 enddo

!computation of the center of mass of the molecule     
 call cenmass(natom,mass,coord,cm)
 write(isave,'(/a)') 'Center of mass (COM):'
 write(isave,'(3a14)') ' x ',' y ',' z '
 write(isave,'(3f14.5)') cm(1),cm(2),cm(3)
      
 write(isave,'(/a)') "Coordinates after COM-shift " 
 write(isave,"(a2, 3a14,3x,a12)") "    ","  x  ",  "  y  " , "  z  ", " mass "
 do i=1,natom
    write(isave,'(a2,3f14.5,3x,f12.5)') atom(i),(coord(i,j),j=1,3),mass(i)
 enddo

!Generation of the intertial tensor of the molecule        
 call inertia(natom,mass,coord,MomIn)
 write(isave,'(/a)') 'Intertial tensor (in Angstrom and atomic mass units)'

 do i=1,3
    write(isave,"(3f16.8)") (MomIn(i,j),j=1,3)
 enddo

!Diagonalization of the inertial tensor
 call jacobi(MomIn,3,3,eig,rotm,eps,1)

 write(isave,"(/a)") 'Principial moments of inertia (eigenvectors of intertial tensor)'
 write(666,"(/a)") 'Principial moments of inertia (eigenvectors of intertial tensor)'
 write(isave,'(a4,5x,4a16)') "    ", "  [amu*A+2]  " , " [g*cm+2*E-40]", " [cm-1]  ", "   [GHz]  "  
 write(666,'(a4,5x,4a16)') "    ", "  [amu*A+2]  " , " [g*cm+2*E-40]", " [cm-1]  ", "   [GHz]  "  
 
 do i=1,3
     e1=0.d0
     e2=0.d0
     e3=0.d0
     if(eig(i) .gt. 1.d-3) then
        e1=eig(i)/av*10.d0  
        e2=tocmm1/eig(i)
        e3=toghz/eig(i)
     else
        eig(i)=0.d0
     endif

     write(isave,'(2x,i2,3x,4f16.5)') i,eig(i),e1,e2,e3
     write(666,'(2x,i2,3x,4f16.5)') i,eig(i),e1,e2,e3
 enddo
    
 write(isave,'(/a)') 'Eigenvectors of inertial tensor (Rotation Matrix)'
 do i=1,3
    write(isave,"(3f16.8)") (rotm(i,j),j=1,3)
 enddo

 call flush(isave)

 end subroutine
!===================================================================



!===================================================================
 subroutine inertia(natom,mass,coord,f)
!===================================================================
!this routine calculates the inertial tensor
!===================================================================
 implicit none 
 integer :: i,natom
 real*8, dimension(natom)   :: mass
 real*8, dimension(natom,3) :: coord
 real*8, dimension(3,3) :: f

 f=0.d0 

 do i=1,natom
    f(1,1) = f(1,1) + mass(i)*(coord(i,2)**2 + coord(i,3)**2) 
    f(2,2) = f(2,2) + mass(i)*(coord(i,3)**2 + coord(i,1)**2) 
    f(3,3) = f(3,3) + mass(i)*(coord(i,1)**2 + coord(i,2)**2) 
    f(2,1) = f(2,1) - mass(i)*coord(i,1)*coord(i,2) 
    f(3,1) = f(3,1) - mass(i)*coord(i,1)*coord(i,3) 
    f(3,2) = f(3,2) - mass(i)*coord(i,2)*coord(i,3) 
 enddo

 f(1,2)=f(2,1)
 f(1,3)=f(3,1)
 f(2,3)=f(3,2)
  
 end subroutine
!================================================================


!================================================================
 subroutine cenmass(natom,mass,coord,cm)
!================================================================
 implicit none
 integer :: i, natom
 real*8, dimension(natom)  :: mass
 real*8, dimension(natom,3):: coord
 real*8, dimension(3)      :: cm
 real*8 :: sumwx, sumwy, sumwz

 sumwx=0.d0
 sumwy=0.d0
 sumwz=0.d0

 do i=1,natom
    sumwx = sumwx + mass(i)*coord(i,1)
    sumwy = sumwy + mass(i)*coord(i,2)
    sumwz = sumwz + mass(i)*coord(i,3)
 enddo

 cm(1)=sumwx/sum(mass)
 cm(2)=sumwy/sum(mass)
 cm(3)=sumwz/sum(mass)

 do i=1,natom
    coord(i,1)=coord(i,1)-cm(1)
    coord(i,2)=coord(i,2)-cm(2)
    coord(i,3)=coord(i,3)-cm(3)
 enddo

 end subroutine
!===============================================================



!===========================================================================
      subroutine jacobi(a,nmax,n,d,x,th,iord)                          
!===========================================================================
!  threshold jacobian method for eigenvalues of symmetric matricies        !
!===========================================================================
! a(*,*)--> input matrix to be diagonalized                                !
! nmax  --> maximum order of the problem                                   !
! n     --> actual order                                                   !
! th    --> treshold                                                       !
! iord  --> sorting variable, +1: increasing order, -1: decrecasing order  !
! x(*,*)--> matrix of the eigenvectors                                     !
! d(*)  --> block of the eigenvalues                                       !
! itmax --> maximum number of the jacobian rotations                       ! 
!===========================================================================

      implicit real*8 (a-h,o-z)
      parameter (itmax=50)
      dimension a(nmax,nmax),x(nmax,nmax),d(nmax)

      do i=1,n
        do j=1,n
         x(i,j)=0.0d0
        enddo
         x(i,i)=1.d0
      enddo

      do i=1,n
        d(i)=a(i,i)
      enddo
!----------------------------------------------------
      do iter=1,itmax
       amax=0.0d0
        do i=2,n
         do 10 j=1,i-1
           aii=d(i)
           ajj=d(j)
           aij=a(i,j)
           if (dabs(aij) .gt. amax) amax=dabs(aij)
           if (dabs(aij) .le. th) goto 10 
           alpha=0.5d0*(aii-ajj)/aij
!   t=s/c     
           t=-alpha+dsqrt(1.d0+alpha*alpha)
           if (alpha .lt. 0.d0)  t=-alpha-dsqrt(1.d0+alpha*alpha)
           c=1.0d0/dsqrt(1.d0 + t*t)
           s=c*t
            
           do 20 k=1,n
              xj=c*x(k,j)-s*x(k,i)
              x(k,i)=s*x(k,j)+c*x(k,i)
              x(k,j)=xj
!             x(k,j)=s*x(k,j)+c*x(k,i)
!             x(k,i)=xj
              if (k .eq. j) goto 20
              if (k .lt. j) then
                xj=c*a(j,k)-s*a(i,k)
                a(i,k)=s*a(j,k)+c*a(i,k)
                a(j,k)=xj
                goto 20
              endif

              if (k .eq. i) goto 20
              if (k .lt. i) then
                xj=c*a(k,j)-s*a(i,k)
                a(i,k)=s*a(k,j)+c*a(i,k)
                a(k,j)=xj
                goto 20
              endif
              xj=c*a(k,j)-s*a(k,i)
              a(k,i)=s*a(k,j)+c*a(k,i)
              a(k,j)=xj
  20       continue
           d(i)=c*c*aii+s*s*ajj+2.d0*s*c*aij
           d(j)=c*c*ajj+s*s*aii-2.d0*s*c*aij
           a(i,j)=0.0d0
  10     continue
       enddo !end of i loop

      if (amax .le. th) goto 30

      enddo !end of iter loop
  30  continue
!--------------------------------------------------      
      if (iord .eq. 1) then
!     arrange eigenvalues in increasing order      
      
      do k=1,n-1
        dmn=d(k)
        kmin=k
        do j=k+1,n
          if (dmn .gt. d(j) ) then 
            kmin=j
            dmn=d(j)
          endif
        enddo
        if (k .ne. kmin) then
          do j=1,n
           call swap( x(j,kmin),x(j,k) )
          enddo
           call swap( d(kmin),d(k) )
        endif
      enddo
!------------------------------------------------
      elseif (iord .eq. -1) then
!     arrange egeinvalues in decreasing order      
      do k=1,n-1
        dmx=d(k)
        kmax=k
        do j=k+1,n
          if (dmx .lt. d(j)) then
            kmax=j
            dmx=d(j)
          endif
        enddo
         if (k .ne. kmax) then
           do j=1,n
            call swap(x(j,kmax), x(j,k))
           enddo
            call swap(d(kmax),d(k))
         endif
      enddo
      endif
!----------------------------------------------
!   restore the original a(i,j) matrix
      do i=1,n
       do j=1,i-1
         a(i,j)=a(j,i)
       enddo  
      enddo
!----------------------------------------------   
      return
      end subroutine
!====================================================================
      subroutine swap(a,b)
      implicit real*8 (a-h,o-z)
       
      temp=a
      a=b
      b=temp
      return
      end subroutine
!====================================================================



