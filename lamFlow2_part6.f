      PROGRAM CHANNEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Laminar Channel flow. Steady. Pseudo-temporal numerical method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      PARAMETER (NR=100,NT=400000)
      !!!!!!
      REAL(8) RY(NR),RDY(NR),DY2(NR),SU(NR)
      COMMON/BP/SU
      COMMON/RAD/ RY,RDY,DY2
!!!!!!
      REAL(8) Y(NR),Y1(NR),DY(NR),RP(NR),OY(NR),U(NR),
     &PDW(NR),FUN(NR),GEF(NR),SF(NR),SFF(NR),FDN(NR)
      REAL(8) T,DT,PQ1,VISM,CK,P0,H,UT,UMAX,REU
      REAL(8) RESU,RES1U
      CHARACTER C15*40
!!!!!!
      DATA T,DT,H,UB,AL/0.0D0,10.D0,1.D0,5.7D-4,1.D-8/
!     U is the flow velocity component in the x-direction, which is along the channel axis
!     UB is the bulk velocity
!     UMAX is the value of U at the channel axis
!     VISM is dynamic viscosity, here for air
!     Y is the coordinate in the direction normal to the channel wall with 0 being at the channel axis
!     DY is the step in the y-direction
!     P0 is the pressure drop in the channel defined as (-1/rho(dP/dx))
!     H is the channel half-width equal to 1 in this problem
!     REU is the characteristic Reynolds number based on the bulk flow velocity and the channel width
!     it should be less than 1350 for a channel flow be laminar
!     T and DT are time and time step
!     RESU and RES1U are residuals to monitor the solution convergence
!     AL is the parameter used to evaluate the solution converegence 
!     SU, SF, SFF, GEF, FUN, FDN, PDW, PQ1 are functions/parameters used in the numerical method     
!     other functions and parameters are auxilary 
      VISM=1.568D-5
      REU=1000.D0
      UB=(0.5D0*REU*VISM)/H  ! WTF GOING ON HERE =>  RE*VISCOSITY/2H, BP
      P0=3.D0*VISM*UB/(H**2.) ! LINEARLY DEPENDENT ON VISCOSITY, CHANGE IN P BECAUSE WALLS RESISITING FLOW, BP
!!!!!!THE BEGINNING OF CALCULATION 	 
	  N=NR   ! WHY REDEFINE? NR DEFINED AS PARAMETER (LINE 5) MESSY CODE, BP
	  NM=N-1 ! WHY? BP
!!!!!! GRID GENERATION
	  CK=1.03D0       ! COEFFICIENT IN GRID GENERATION
	  !CK=1.00D0 
	  y(1)=dble(0.0)  ! SPACE DIRECTION
        y(2)=dble(1.0)
      DO j=1,n-2
      y(j+2)=dble(y(j+1)+CK*(y(j+1)-y(j)))  ! MAKES NODES MORE DISTANTLY SPACED, BP
      enddo
      DO j=1,n
      y(j)=dble(y(j)/y(n)) ! NORMALIZES SO LARGEST Y-VALUE IS EQUAL TO 1, BP
      ENDDO 
      DO j=1,n
          y(j)=1.-y(j) !FLIPS IT,(WEIRD, WHY?) IN FACT, NOW, THE SMALLER NODE VALUES ARE MORE DISTANTLY SPACED, BP
      ENDDO        ! DO THIS CORRECTLY TO BEGIN WITH, FIXING WHAT THEY BROKE, BP
	  do j=1,n
      y1(j)=y(N+1-j)  ! CHANGES THE INDEXING TO PUT Y'S IN ASCENDING ORDER, AGAIN ONLY NEEDED BECAUSE
	  enddo         ! NOT DONE CORRECTLY IN THE BEGINNING. Y=0 START POINT, NOT SO IMPORTANT, BP
	  y=y1*H  ! OUR H IS JUST 1, BP
!!!!!!! THE END OF GRID GENERATION
!!!!!!!!!!!!!! initial profiles
      U=UB ! approximation for the flow velocity profile required to start simulations
	  U(N)=0.D0  !non-slip condition on the channel wall
!!!!!!!!!!!!!! auxilary functions
      CALL RADI(NR,NM,Y,DY,RP,OY) ! ONE OF 5 SUBROUTINES AT BOTTOM OF FILE, NO COMMENT, BP
                                  ! ALL 5 ONLY CALLED ONCE, WHY MAKE A SUBROUTINE, BP
      CALL MIZES(NR,NM,PDW,Y)     ! ONE OF 5 SUBROUTINES AT BOTTOM OF FILE, NO COMMENT, ONLY CALLED ONCE "", BP
     	PQ1=(Y(3)/Y(2))**2. 
      OPEN(4,FILE='RES.DAT')
      KM=1   ! WHAT IS THIS? BP
!!!!!! START THE MAIN CYCLE
     	DO 7 M=1,NT    ! WHY DON'T PEOPLE INDENT, STARTS LOOP, BP
      UMAX=MAXVAL(U)
	  WRITE(*,*) M,UMAX
	  SU=0.D0
      RESU=0.D0
      DO I=1,NM
      RESU=RESU+U(I)
      ENDDO
      DO I=1,N
      SU(I)=P0
      ENDDO   
!     EQUIVALENTLY, YOU CAN JUST WRITE SU=P0 IF YOUR PLATFORM CAN TAKE IT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!THE BEGINNING OF CONTROL VOLUME BLOCK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL COEGEF(NR,N,GEF,VISM,RP)
      CALL COESF(NR,NM,SF)
      SFF=0.D0
      DO I=1,N
      FUN(I)=U(I)
      ENDDO
      CALL SOLVE(NR,NM,GEF,SF,SFF,FUN,DT,PDW,PQ1,FDN)
      DO I=1,N
      U(I)=FDN(I)
      ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!THE END OF CONTROL VOLUME BLOCK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      T=T+DT  !!!! TRANSITION TO THE NEXT STEP IN TIME
!!!!!! CALCULATION OF RESIDUALS      
      RES1U=0.D0
      DO I=1,NM
      RES1U=RES1U+U(I)
      ENDDO
      RESU=DABS((RES1U-RESU)/RESU)
      IF(M.EQ.KM)THEN
      WRITE(4,99) M,RESU
      KM=KM+100
      ENDIF
      IF(RESU.LE.AL)GOTO 10  !!! CACULATIONS ARE TERMINATED IF RESIDUALS HAVE REACHED THE LIMIT
!!!!!! THE END OF CALCULATION OF RESIDUALS       
!!!!!!
  7   CONTINUE
!!!!!! THE END OF THE MAIN CYCLE
  10  CONTINUE     
      CLOSE(4)
      C15='OUTPUT.DAT'
      UMAX=MAXVAL(U)
      OPEN(3,FILE=C15)
	  WRITE(3,48) (Y1(I),U(I)/UMAX,I=1,N)
      CLOSE(3)
      C15='OUTPUT2.DAT'
      OPEN(2,FILE=C15)
	  WRITE(2,48) (Y1(I), 1.0 - Y1(I)**2.0 ,I=1,N)   !!! Y1 = y/H (NORMALIZED), BP
      CLOSE(2)
  48  FORMAT(2F16.8)
  99  FORMAT(I8,6F25.20)  
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FOR GRID CALCS, CREATES A WHOLE SPRAY OF VARIABLES, NOTHING INSIGHTFUL, BP
      SUBROUTINE RADI(NR,NM,Y,DY,RP,OY) ! SO MESSY, WHY MAKE SOME OUTPUTS COMMON AND PLACE OTHERS IN SUBROUTINE DEF,BP
      REAL(8) RY(100),RDY(100),DY2(100)
	  COMMON/RAD/ RY,RDY,DY2
      REAL(8) DY(NR),RP(NR),Y(NR),OY(NR)
      N=NM+1   ! LINE 135 IS ONLY TIME THIS USED, RIDICULOUS TO DEFINE, BP
	  OY=1.-Y
      DO 1 I=1,NM
      IP=I+1
      RY(IP)=1./Y(IP)
      DY(I)=Y(IP)-Y(I)
  1   RDY(I)=1./DY(I)
      DO 2 I=2,N
  2   RP(I)=RDY(I-1)
      DO 3 I=2,NM
  3   DY2(I)=1./(DY(I-1)+DY(I))
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! FOR GRID CALCS, CREATES PDW, MAYBE USEFUL TO NUM ROUTINE, BP 
      SUBROUTINE MIZES(NR,NM,PDW,Y)
      REAL(8) PDW(NR),Y(NR),W(NR)
      DO I=1,NM
	  W(I)=(Y(I+1)+Y(I))/2.0D0
	  ENDDO
      DO I=2,NM 
      PDW(I)=W(I)-W(I-1)
	  ENDDO
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      SUBROUTINE COEGEF(NR,N,GEF,VISM,RP)
      REAL(8) GEF(NR),RP(NR),VISM
      DO 1 I=1,N
      GEF(I)=VISM*RP(I)
  1   CONTINUE
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE COESF(NR,NM,SF)
      REAL(8) SU(100)     
      COMMON/BP/SU
	  REAL(8) SF(NR)
      DO 1 I=2,NM
      SF(I)=SU(I)
   1  CONTINUE
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE SOLVE(NR,NM,GEF,SF,SFF,FU,DT,PDW,PQ1,FD)
      REAL(8) GEF(NR),SF(NR),SFF(NR),FU(NR),PDW(NR),FD(NR)
      REAL(8) A(NR),B(NR),C(NR),D(NR),P(NR),Q(NR)
     	REAL(8) DT,PQ1,R,RDT,R1
      N=NM+1
      DO 1 I=2,NM
      A(I)=GEF(I+1)
      B(I)=GEF(I)
  1   CONTINUE
      RDT=1.D0/DT
      DO 2 I=2,NM
      C(I)=PDW(I)*(FU(I)*RDT+SF(I))
  2   D(I)=A(I)+B(I)+PDW(I)*(RDT-SFF(I))
      P(1)=0.D0    
      Q(1)=0.D0   
      R1=1.D0/((PQ1-1.D0)*A(2)-B(2))   
      P(1)=(PQ1*A(2)-D(2))*R1
      Q(1)=C(2)*R1
      !!!!!!DIRECT CALCULATION
      DO 3 I=2,NM
      IM=I-1
      R=1.D0/(D(I)-B(I)*P(IM))
      P(I)=A(I)*R
  3   Q(I)=(C(I)+B(I)*Q(IM))*R
!!!!!!REVERSE CALCULATION
      FD(N)=FU(N)
      DO 4 IR=1,NM
      I=NM-IR+1
  4   FD(I)=P(I)*FD(I+1)+Q(I)
      RETURN
      END
