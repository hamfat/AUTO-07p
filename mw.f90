!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

! Evaluates the algebraic equations or ODE right hand side

! Input arguments :
!      NDIM   :   Dimension of the ODE system 
!      U      :   State variables
!      ICP    :   Array indicating the free parameter(s)
!      PAR    :   Equation parameters

! Values to be returned :
!      F      :   ODE right hand side values

! Normally unused Jacobian arguments : IJAC, DFDU, DFDP (see manual)

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, IJAC, ICP(*)
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)
      
      DOUBLE PRECISION V,N,MINF,NINF,Lamda, gl, vl, gk, vk, gca, v1, v2, v3, v4, psi

       V=U(1)
       N=U(2)

       MINF=0.5*(1+tanh((V-PAR(6))/PAR(7)))
       NINF=0.5*(1+tanh((V-PAR(8))/PAR(9)))
       Lamda=cosh((V-PAR(8))/2*PAR(9))

       F(1)=-PAR(1)*(V-PAR(2))-PAR(3)*N*(V-PAR(4))-PAR(5)*MINF*(V-1)
       F(2)=PAR(10)*Lamda*(NINF-N)

      END SUBROUTINE FUNC
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----

! Input arguments :
!      NDIM   :   Dimension of the ODE system 

! Values to be returned :
!      U      :   A starting solution vector
!      PAR    :   The corresponding equation-parameter values
!      T      :	  Not used here

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

      DOUBLE PRECISION gl, vl, gk, vk, gca, v1, v2, v3, v4, psi
      !DOUBLE PRECISION, EXTERNAL :: GETP

! Initialize the equation parameters
    gl=PAR(1)
	vl=PAR(2)
	gk=PAR(3)
	vk=PAR(4)
	gca=PAR(5)
	v1=PAR(6)
	v2=PAR(7)
	v3=PAR(8)
	v4=PAR(9)
	psi=PAR(10)
	


       PAR(1)=0.25!0.1429
       PAR(2)= -0.875 ! -1.2
       PAR(3)=1.0!0.5715
       PAR(4)=-1.125
       PAR(5)=0.4997!0.2856
       PAR(6)=-0.5
       PAR(7)=0.3125
       PAR(8)=-0.1380!-0.1875!-0.1503
       PAR(9)=0.1812 
       PAR(10)=0.1665!0.0952
       

! Initialize the solution

    
      !U(1)=-0.2245
      !U(2)=0.3994
	  
  ! v1=-0.23
      U(1)=-0.18523
      U(2)=0.37255
	  
	
      END SUBROUTINE STPNT
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! The following subroutines are not used here,
! but they must be supplied as dummy routines

      SUBROUTINE BCND 
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      

    SUBROUTINE PVLS(NDIM,U,PAR)
!--------- ----

     IMPLICIT NONE
     INTEGER, INTENT(IN) :: NDIM
     DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
     DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
 ! Set PAR(3) equal to the minimum of U(2)
     !X PAR(18)=GETP('MIN',1,U)
     ! DOUBLE PRECISION, EXTERNAL :: GETP
      !set PAR(14) to be min of U(1)
      !PAR(16)=GETP('MIN', 1, U)

    END SUBROUTINE PVLS
!----------------------------------------------------------------------
!----------------------------------------------------------------------
