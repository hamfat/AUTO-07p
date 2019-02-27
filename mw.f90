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
      
      DOUBLE PRECISION V,N,MINF,NINF,Lamda, gl, vl, gk, vk, gca, va, vb, vc, vd, psi

       V=U(1)
       N=U(2)

       MINF=0.5*(1+tanh(V*PAR(6)-PAR(7)))
       NINF=0.5*(1+tanh(V*PAR(8)-PAR(9)))
       Lamda=cosh((V*PAR(8)-PAR(9))/2)

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

      DOUBLE PRECISION gl, vl, gk, vk, gca, va, vb, vc, vd, psi
      !DOUBLE PRECISION, EXTERNAL :: GETP

! Initialize the equation parameters
        gl=PAR(1)
	vl=PAR(2)
	gk=PAR(3)
	vk=PAR(4)
	gca=PAR(5)
	va=PAR(6)
	vb=PAR(7)
	vc=PAR(8)
	vd=PAR(9)
	psi=PAR(10)
	


       PAR(1)=0.25
       PAR(2)=-0.875
       PAR(3)=1.0
       PAR(4)=-1.125
       PAR(5)=0.4997
       PAR(6)=3.1999
       PAR(7)=-1.184!-0.8999 
       PAR(8)=5.5172
       PAR(9)=-0.7586!-2.5   
       PAR(10)=0.1665
       

! Initialize the solution

!vd=-2.5
       !U(1)=-0.7677
     
      
      ! U(2)=0.03014
       

!vb=-0.7088
 !U(1)=-0.7769
  !U(2)=0.002198

  ! vb=-1.184
      U(1)=-0.213
      U(2)=0.305

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
