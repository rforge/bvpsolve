C KARLINE: change name SOL -> SOL_GAM
C   dec -> DEC_GAM, interp -> INTERP_GAM
C 

       SUBROUTINE   gam(R,FCN,T0,Y0,TEND,H,
     &                  RTOL,ATOL,ITOL,
     &                  JAC ,IJAC,MLJAC,MUJAC,
     &                  SOLOUT,IOUT,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)




C
C     VERSION: NOVEMBER 25, 1999
C
C     PURPOSE: THE GAM CODE NUMERICALLY SOLVES A (POSSIBLY STIFF)
C              SYSTEM OF FIRST 0RDER ORDINARY DIFFERENTIAL
C              EQUATIONS IN THE FORM  Y'=F(X,Y), WITH A GIVEN
C              INITIAL CONDITION.
C
C     AUTHORS: F. IAVERNARO AND F. MAZZIA
C              UNIVERSITA' DEGLI STUDI DI BARI,
C             DIPARTIMENTO DI MATEMATICA
C             VIA ORABONA 4, 70125 BARI, ITALY
C               E-MAIL:  LABOR@ALPHAMATH.DM.UNIBA.IT
C
C     METHODS: THE METHODS USED ARE IN THE CLASS OF BOUNDARY VALUE
C              METHODS (BVMs), NAMELY THE GENERALIZED ADAMS METHODS
C              (GAMs) OF ORDER 3-5-7-9 WITH STEP SIZE CONTROL.
C
C  REFERENCES: L.BRUGNANO, D.TRIGIANTE,  Solving Differential Problems
C              by Multistep Initial and Boundary Value Methods,
C              Gordon & Breach.
C
C              F.IAVERNARO, F.MAZZIA,  Block-Boundary Value Methods
C              for the solution of Ordinary Differential Equations,
C              (submitted).
C
C              F.IAVERNARO, F.MAZZIA,  Solving Ordinary Differential
C              Equations by Generalized Adams Methods: properties and
C              implementation techniques,
C              proceedings of NUMDIFF8, Appl. Numer. Math. (to appear).
C
C DESCRIPTION: THE GAM CODE CONSISTS OF THREE FILES:
C            - gam.f CONTAINS THE MAIN SUBROUTINES THAT IMPLEMENT THE
C              INTEGRATION PROCEDURE;
C            - subgam.f  CONTAINS THE ADDITIONAL LINEAR ALGEBRA ROUTINES
C              REQUIRED BY GAM.F PLUS SOME OTHER SUBROUTINES PROPER OF
C              THE USED METHODS;
C            - param.dat CONTAINS THE DEFINITION OF ALL THE COSTANTS.
C              THIS FILE MUST BE INSERTED IN THE SAME DIRECTORY OF THE MAIN 
C              PROGRAM; IT IS INCLUDED IN ALL THE SUBROUTINE WITH THE
C              FORTRAN 77 INSTRUCTION
C              INCLUDE "param.dat" .
C
C    COMMENTS: THE PHILOSOFY AND THE STYLE USED IN WRITING THE CODE ARE VERY
C              SIMILAR TO THOSE CHARACTERIZING THE CODE  RADAU5.
C              INDEED THE AUTHORS IMPORTED FROM RADAU5 SOME SUBROUTINES,
C              COMMENTS AND IMPLEMENTATION THECNIQES LEAVING UNCHANGED
C              THE NAME AND THE MEANING OF  A NUMBER OF VARIABLES.
C              THE AUTHORS ARE VERY GRATEFUL TO ANYONE USES THE CODE AND
C              WOULD APPRECIATE ANY CRITICISM AND REMARKS ON HOW IT PERFORMS.
C
C
C
C -------------------------------------------------------------------------
C     INPUT PARAMETERS
C -------------------------------------------------------------------------
C     R           DIMENSION OF THE SYSTEM
C
C     FCN         NAME (EXTERNAL) OF THE SUBROUTINE COMPUTING THE
C                 VALUE OF F(T,Y):
C                    SUBROUTINE FCN(R,T,Y,F,RPAR,IPAR)
C                    DOUBLE PRECISION X,Y(R),F(R)
C                    F(1)=...   ETC.
C                 (RPAR, IPAR    SEE BELOW)
C
C     T0          INITIAL T-VALUE
C
C     Y0          INITIAL VALUES FOR Y
C
C     TEND        FINAL T-VALUE (TEND-T MUST BE POSITIVE)
C
C     H           INITIAL STEP SIZE GUESS;
C                 FOR STIFF EQUATIONS WITH INITIAL TRANSIENT,
C                 H=1.D0/(NORM OF F'), USUALLY 1.D-3 OR 1.D-5, IS GOOD.
C                 THIS CHOICE IS NOT VERY IMPORTANT, THE STEP SIZE IS
C                 QUICKLY ADAPTED. (IF H=0.D0, THE CODE PUTS H=1.D-6).
C
C     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. THEY
C                 CAN BE BOTH SCALARS OR ELSE BOTH VECTORS OF LENGTH R.
C
C     ITOL        SWITCH FOR RTOL AND ATOL:
C                   ITOL=0: BOTH RTOL AND ATOL ARE SCALARS.
C                     THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF
C                     Y(I) BELOW RTOL*ABS(Y(I))+ATOL
C                   ITOL=1: BOTH RTOL AND ATOL ARE VECTORS.
C                     THE CODE KEEPS THE LOCAL ERROR OF Y(I) BELOW
C                     RTOL(I)*ABS(Y(I))+ATOL(I).
C
C     JAC         NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES
C                 THE PARTIAL DERIVATIVES OF F(T,Y) WITH RESPECT TO Y
C                 (THIS ROUTINE IS ONLY CALLED IF IJAC=1; SUPPLY
C                 A DUMMY SUBROUTINE IN THE CASE IJAC=0).
C                 FOR IJAC=1, THIS SUBROUTINE MUST HAVE THE FORM
C                    SUBROUTINE JAC(R,T,Y,DFY,LDFY,RPAR,IPAR)
C                    DOUBLE PRECISION T,Y(R),DFY(LDFY,R)
C                    DFY(1,1)= ...
C                 LDFY, THE COLUMN-LENGTH OF THE ARRAY, IS
C                 FURNISHED BY THE CALLING PROGRAM.
C                 IF (MLJAC.EQ.R) THE JACOBIAN IS SUPPOSED TO
C                    BE FULL AND THE PARTIAL DERIVATIVES ARE
C                    STORED IN DFY AS
C                       DFY(I,J) = PARTIAL F(I) / PARTIAL Y(J)
C                 ELSE, THE JACOBIAN IS TAKEN AS BANDED AND
C                    THE PARTIAL DERIVATIVES ARE STORED
C                    DIAGONAL-WISE AS
C                       DFY(I-J+MUJAC+1,J) = PARTIAL F(I) / PARTIAL Y(J).
C
C     IJAC        SWITCH FOR THE COMPUTATION OF THE JACOBIAN:
C                    IJAC=0: JACOBIAN IS COMPUTED INTERNALLY BY FINITE
C                       DIFFERENCES, SUBROUTINE "JAC" IS NEVER CALLED.
C                    IJAC=1: JACOBIAN IS SUPPLIED BY SUBROUTINE JAC.
C
C     MLJAC       SWITCH FOR THE BANDED STRUCTURE OF THE JACOBIAN:
C                    MLJAC=R: JACOBIAN IS A FULL MATRIX. THE LINEAR
C                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION.
C                       0<=MLJAC<R: MLJAC IS THE LOWER BANDWITH OF JACOBIAN
C                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW
C                       THE MAIN DIAGONAL).
C
C     MUJAC       UPPER BANDWITH OF JACOBIAN  MATRIX (>= NUMBER OF NON-
C                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL).
C                 NEED NOT BE DEFINED IF MLJAC=R.
C
C
C     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE
C                 NUMERICAL SOLUTION DURING INTEGRATION.
C                 IF IOUT=1, IT IS CALLED AFTER EVERY SUCCESSFUL STEP.
C                 SUPPLY A DUMMY SUBROUTINE IF IOUT=0.
C                 IT MUST HAVE THE FORM
C                    SUBROUTINE SOLOUT(R,TP,YP,FF,NT,DBLK,ORD,RPAR,IPAR,IRTRN)
C                    INTEGER R, DBLK, ORD, IPAR(*), IRTRN, NT
C                    DOUBLE PRECISION TP(*), YP(R,*), RPAR(*), FF(R,*)
C                    ....
C                 SOLOUT FURNISHES THE SOLUTION "YP" AT THE
C                    GRID-POINTS "TP(*)".
C                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN
C                    IS SET <0, GAM  RETURNS TO THE CALLING PROGRAM.
C
C                 CONTINUOUS OUTPUT:
C                 DURING CALLS TO "SOLOUT", A CONTINUOUS SOLUTION
C                 FOR THE INTERVAL [TP(1),TP(DBLK+1)] IS AVAILABLE THROUGH
C                 THE FUNCTION
C                        >>>   CONTR(I,R,T,TP,FF,DBLK,NT)   <<<
C                 WHICH PROVIDES AN APPROXIMATION TO THE I-TH
C                 COMPONENT OF THE SOLUTION AT THE POINT T. THE VALUE
C                 T SHOULD LIE IN THE INTERVAL [TP(1),TP(DBLK+1)] ON
C                 WHICH THE SOLUTION IS COMPUTED AT CURRENT STEP.
C                 DO NOT CHANGE THE ENTRIES OF FF(R,*) and NT, IF THE
C                 DENSE OUTPUT FUNCTION IS USED.
C
C     IOUT        SWITCH FOR CALLING THE SUBROUTINE SOLOUT:
C                    IOUT=0: SUBROUTINE IS NEVER CALLED
C                    IOUT=1: SUBROUTINE IS AVAILABLE FOR OUTPUT.
C
C     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK".
C                 WORK(1), WORK(2),.., WORK(19) SERVE AS PARAMETERS
C                 FOR THE CODE. FOR STANDARD USE OF THE CODE
C                 WORK(1),..,WORK(18) MUST BE SET TO ZERO BEFORE
C                 CALLING. SEE BELOW FOR A MORE SOPHISTICATED USE.
C                 WORK(19),..,WORK(LWORK) SERVE AS WORKING SPACE
C                 FOR ALL VECTORS AND MATRICES.
C                 "LWORK" MUST BE AT LEAST:
C
C                 IN THE USUAL CASE WHERE THE JACOBIAN IS FULL
C                 STORAGE REQUIREMENT IS
C                             LWORK = 2*R*R+42*R+18.
C
C                 IN THE CASE WHERE THE JACOBIAN IS SPARSE
C                 STORAGE REQUIREMENT IS
C                             LWORK = (3*MLJAC+2*MUJAC+44)*R+18
C
C
C     LWORK       DECLARED LENGTH OF ARRAY "WORK".
C
C     IWORK       INTEGER WORKING SPACE OF LENGTH "LIWORK".
C                 IWORK(1),IWORK(2),...,IWORK(24) SERVE AS PARAMETERS
C                 FOR THE CODE. FOR STANDARD USE, SET IWORK(1),..,
C                 IWORK(8) TO ZERO BEFORE CALLING.
C                 IWORK(10),...,IWORK(LIWORK) SERVE AS WORKING AREA.
C                 "LIWORK" MUST BE AT LEAST
C
C                          LIWORK = 24 + R
C
C     LIWORK      DECLARED LENGTH OF ARRAY "IWORK".
C
C     RPAR, IPAR  REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS) WHICH
C                 CAN BE USED FOR COMMUNICATION BETWEEN YOUR CALLING
C                 PROGRAM AND THE FCN, JAC SUBROUTINES.
C
C -------------------------------------------------------------------------
C     SOPHISTICATED SETTING OF PARAMETERS
C -------------------------------------------------------------------------
C              SEVERAL PARAMETERS OF THE CODE ARE TUNED TO MAKE IT WORK
C              WELL. THEY MAY BE DEFINED BY SETTING WORK(1),...
C              AS WELL AS IWORK(1),... DIFFERENT FROM ZERO.
C              FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES:
C
C    IWORK(1)  NOT USED
C
C    IWORK(2)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS.
C              THE DEFAULT VALUE (FOR IWORK(2)=0) IS 100000.
C
C    IWORK(3)  ORDMIN, 3 <= ORDMIN <= 9,
C
C    IWORK(4)  ORDMAX, ORDMIN <= ORDMAX <= 9
C
C    IWORK(5)  THE MAXIMUM NUMBER OF SPLITTING-NEWTON ITERATIONS FOR THE
C              SOLUTION OF THE IMPLICIT SYSTEM IN EACH STEP FOR ORDER 3.
C              THE DEFAULT VALUE (FOR IWORK(5)=0) IS 10.
C
C    IWORK(6)  THE MAXIMUM NUMBER OF SPLITTING-NEWTON ITERATION FOR
C              ORDER 5, THE DEFAULT VALUE (FOR IWORK(6)=0) IS 18.
C
C    IWORK(7)  THE MAXIMUM NUMBER OF SPLITTING-NEWTON ITERATION FOR
C              ORDER 7, THE DEFAULT VALUE (FOR IWORK(7)=0) IS 26.
C
C    IWORK(8)  THE MAXIMUM NUMBER OF SPLITTING-NEWTON ITERATION FOR
C              ORDER 9, THE DEFAULT VALUE (FOR IWORK(5)=0) IS 36.
C
C
C    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 1.D-16.
C
C    WORK(2)   HMAX  MAXIMAL STEP SIZE, DEFAULT TEND-T0.
C
C    WORK(3)   THET DECIDE WHETHER THE JACOBIAN SHOULD BE RECOMPUTED
C
C    WORK(4)   FACNEWT:  stopping criterion for splitting-Newton method
C                 for small values of min(abs(y_i)) and min(abs(f_j)).
C
C    WORK(5)   TETAK0 stopping criterium for the splitting-Newton method
C                the iterates must be decreasing by a factor tetak0
C
C    WORK(6)   CS(2): EMPIRICAL COMPUTATIONAL COST FOR ORDER  5 METHOD
C              USED IN THE ORDER VARIATION STRATEGY
C              (DEFAULT WORK(6) = 2.4D0)
C
C    WORK(7)   CS(3): EMPIRICAL COMPUTATIONAL COST FOR ORDER  7 METHOD
C              USED IN THE ORDER VARIATION STRATEGY
C              (DEFAULT WORK(6) = 4.0D0)
C
C    WORK(8)   CS(4): EMPIRICAL COMPUTATIONAL COST FOR ORDER  9 METHOD
C              USED IN THE ORDER VARIATION STRATEGY
C              (DEFAULT WORK(6) = 7.2D0)
C
C    WORK(9)-WORK(10)   FACL-FACR: PARAMETERS FOR STEP SIZE SELECTION
C               THE NEW STEPSIZE IS CHOSEN SUBJECT TO THE RESTRICTION
C               FACL  <=  HNEW/HOLD <= FACR
C               (DEFAULT WORK(9) = 0.12, WORK(10) = 10 )
C
C    WORK(11)  SFDOWN:SAFETY FACTOR IN STEP SIZE PREDICTION
C                  USED FOR THE LOWER ORDER METHOD
C                  (DEFAULT WORK(11) = 20D0)
C
C    WORK(12)  SFUP:SAFETY FACTOR IN STEP SIZE PREDICTION
C                  USED FOR THE UPPER ORDER METHOD
C                  (DEFAULT WORK(12) = 40D0)
C
C    WORK(13)  SFSAME: SAFETY FACTOR IN STEP SIZE PREDICTION
C                  USED FOR THE CURRENT ORDER METHOD
C                  (DEFAULT WORK(13) = 18D0)
C
C    WORK(14)  SF: SAFETY FACTOR IN STEP SIZE PREDICTION
C                  USED FOR THE CURRENT ORDER METHOD WHEN IS
C                  FAILED THE ERROR CONTROL TEST (DEFAULT WORK(14) = 15D0)
C
C
C
C    WORK(15)  FACNEWT stopping criterion for splitting-Newton method ORDER 3
C                  (DEFAULT WORK(15) = 1.0D-3)
C
C    WORK(16)  FACNEWT stopping criterion for splitting-Newton method ORDER 5
C                  (DEFAULT WORK(16) = 9.0D-2)
C
C    WORK(17)  FACNEWT stopping criterion for splitting-Newton method ORDER 7
C                  (DEFAULT WORK(17) = 9.0D-1)
C
C    WORK(18)  FACNEWT stopping criterion for splitting-Newton method ORDER 9
C                  (DEFAULT WORK(18) = 9.9D-1)
C
C -------------------------------------------------------------------------
C     OUTPUT PARAMETERS
C -------------------------------------------------------------------------
C     T0          T-VALUE FOR WHICH THE SOLUTION HAS BEEN COMPUTED
C                 (AFTER SUCCESSFUL RETURN T0=TEND).
C
C     Y(N)        NUMERICAL SOLUTION AT T0
C
C     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP
C
C     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN:
C                   IDID= 1  COMPUTATION SUCCESSFUL,
C                   IDID=-1  INPUT IS NOT CONSISTENT,
C                   IDID=-2  LARGER NMAX IS NEEDED,
C                   IDID=-3  STEP SIZE BECOMES TOO SMALL,
C                   IDID=-4  MATRIX IS REPEATEDLY SINGULAR.
C
C   IWORK(10)  NFCN    NUMBER OF FUNCTION EVALUATIONS (THOSE FOR NUMERICAL
C                      EVALUATION OF THE JACOBIAN ARE NOT COUNTED) 
C KS CHANGED THAT:NFCN = total nr of evaluations
C   IWORK(11)  NJAC    NUMBER OF JACOBIAN EVALUATIONS (EITHER ANALYTICALLY
C                      OR NUMERICALLY)
C   IWORK(12)  NSTEP(1)  NUMBER OF COMPUTED STEPS   ORD 3
C   IWORK(13)  NSTEP(2)  NUMBER OF COMPUTED STEPS   ORD 5
C   IWORK(14)  NSTEP(3)  NUMBER OF COMPUTED STEPS   ORD 7
C   IWORK(15)  NSTEP(4)  NUMBER OF COMPUTED STEPS   ORD 9
C   IWORK(16)  NNEWT(1)  NUMBER OF REJECTED STEPS (DUE TO NEWTON CONVERGENCE) 3
C   IWORK(17)  NNEWT(2)  NUMBER OF REJECTED STEPS (DUE TO NEWTON CONVERGENCE) 5
C   IWORK(18)  NNEWT(3)  NUMBER OF REJECTED STEPS (DUE TO NEWTON CONVERGENCE) 7
C   IWORK(19)  NNEWT(4)  NUMBER OF REJECTED STEPS (DUE TO NEWTON CONVERGENCE) 9
C   IWORK(20)  NERR(1)   NUMBER OF REJECTED STEPS (DUE TO ERROR TEST) 3
C   IWORK(21)  NERR(2)   NUMBER OF REJECTED STEPS (DUE TO ERROR TEST) 5
C   IWORK(22)  NERR(3)   NUMBER OF REJECTED STEPS (DUE TO ERROR TEST) 7
C   IWORK(23)  NERR(4)   NUMBER OF REJECTED STEPS (DUE TO ERROR TEST) 9
C   IWORK(24)  NDEC      NUMBER OF LU-DECOMPOSITIONS
C-----------------------------------------------------------------------
C     DECLARATIONS
C -------------------------------------------------------------------------
       IMPLICIT NONE
C
C   INPUT VARIABLES
C------------------------------------
       INTEGER R, ORDMIN, ORDMAX, IDID, ITOL, IJAC, MLJAC, MUJAC, IOUT,
     &         IPAR(*), ITINT(4), ITMAX, IJOB, NMAX, LDJAC, LDLU
   
       DOUBLE PRECISION TEND, ATOL(*), RTOL(*), RPAR(*), FACNORD(4),
     &                  HMAX, THET, FACNEWT, TETAK0, CS(4), FACL, FACR,
     &                  SFDOWN, SFUP, SFSAME, SF, UROUND

C
C   OUTPUT VARIABLES
C------------------------------------
       INTEGER NDEC, NFCN, NJAC, NSTEP(4), NNEWT(4), NERR(4)
C
C   INPUT/OUTPUT VARIABLES
C------------------------------------
       INTEGER  LIWORK, LWORK, IWORK(LIWORK)
       DOUBLE PRECISION T0, Y0(R), H, WORK(LWORK)

C
C   LOCAL VARIABLES
C------------------------------------
       INTEGER  IEYP, IEFP, IEDN,IEF,IEF1, IEJF0,IELU,ISTORE,
     &          I, IEIPIV, IESC
       LOGICAL  ARRET, JBAND
      CHARACTER*80 MSG

C
C   EXTERNAL FUNCTIONS
C------------------------------------
       EXTERNAL FCN,JAC, SOLOUT
C -------------------------------------------------------------------------
C     SETTING THE PARAMETERS
C -------------------------------------------------------------------------
      NFCN    =0
      NJAC    =0
      NSTEP(1)=0
      NSTEP(2)=0
      NSTEP(3)=0
      NSTEP(4)=0
      NNEWT(1)=0
      NNEWT(2)=0
      NNEWT(3)=0
      NNEWT(4)=0
      NERR(1) =0
      NERR(2) =0
      NERR(3) =0
      NERR(4) =0
      NDEC    =0
      ARRET   = .FALSE.
C -------- NMAX := THE MAXIMAL NUMBER OF STEPS -----
      IF (IWORK(2).EQ.0) THEN
         NMAX=100000
      ELSE
         NMAX=IWORK(2)
         IF (NMAX.LE.0) THEN
            WRITE(msg,*)' WRONG INPUT IWORK(2)=',IWORK(2)
            CALL rwarn(msg)
            ARRET=.TRUE.
         END IF
      END IF
C--------- ORDMIN  :=  MINIMAL ORDER
      IF (IWORK(3).EQ.0) THEN
         ORDMIN = 1
      ELSE
         ORDMIN=(IWORK(3)-1)/2
         IF ((ORDMIN.LE.0).OR.(ORDMIN.GT.4)) THEN
            WRITE(msg,*)' CURIOUS INPUT IWORK(3)=',IWORK(3)
            CALL rwarn(msg)
            ARRET=.TRUE.
         END IF
      END IF
C--------- ORDMAX :=  MAXIMAL ORDER
      IF (IWORK(4).EQ.0) THEN
         ORDMAX = 4
      ELSE
         ORDMAX=(IWORK(4)-1)/2
         IF ((ORDMAX.LE.0).OR.(ORDMAX.GT.4).OR.(ORDMAX.LT.ORDMIN)) THEN
            WRITE(msg,*)' CURIOUS INPUT IWORK(4)=',IWORK(4)
            CALL rwarn(msg)
            ARRET=.TRUE.
         END IF
      END IF
C -------- ITINT(1) :=  NUMBER OF SPLITTING-NEWTON ITERATIONS ORD 3
      IF (IWORK(5).EQ.0) THEN
         ITINT(1)=12
      ELSE
         ITINT(1)=IWORK(5)
         IF (ITINT(1).LE.0) THEN
            WRITE(msg,*)' CURIOUS INPUT IWORK(5)=',IWORK(5)
            CALL rwarn(msg)
            ARRET=.TRUE.
         END IF
      END IF
      ITMAX = ITINT(1)
C -------- ITINT(2) :=  NUMBER OF SPLITTING-NEWTON ITERATIONS ORD 5
      IF (IWORK(6).EQ.0) THEN
         ITINT(2)=18
      ELSE
         ITINT(2)=IWORK(6)
         IF (ITINT(2).LT.0) THEN
            WRITE(msg,*)' CURIOUS INPUT IWORK(6)=',IWORK(6)
            CALL rwarn(msg)
            ARRET=.TRUE.
         END IF
      END IF
C -------- ITINT(3) :=  NUMBER OF SPLITTING-NEWTON ITERATIONS ORD 7
      IF (IWORK(7).EQ.0) THEN
         ITINT(3)= 26
      ELSE
         ITINT(3)=IWORK(7)
         IF (ITINT(3).LT.0) THEN
            WRITE(msg,*)' CURIOUS INPUT IWORK(8)=',IWORK(7)
            CALL rwarn(msg)
            ARRET=.TRUE.
         END IF
      END IF
C -------- ITINT(4) :=  NUMBER OF SPLITTING-NEWTON ITERATIONS ORD 9
      IF (IWORK(8).EQ.0) THEN
         ITINT(4)= 36
      ELSE
         ITINT(4)=IWORK(8)
         IF (ITINT(4).LT.0) THEN
            WRITE(msg,*)' CURIOUS INPUT IWORK(8)=',IWORK(8)
            CALL rwarn(msg)
            ARRET=.TRUE.
         END IF
      END IF

C -------- UROUND :=  SMALLEST NUMBER SATISFYING 1.0D0+UROUND>1.0D0
      IF (WORK(1).EQ.0.0D0) THEN
         UROUND=1.0D-16
      ELSE
         UROUND=WORK(1)
         IF (UROUND.LE.1.0D-19.OR.UROUND.GE.1.0D0) THEN
            WRITE(msg,*)'COEFFICIENTS HAVE 20 DIGITS, UROUND=',WORK(1)
            CALL rwarn(msg)
            ARRET=.TRUE.
         END IF
      END IF
C -------- HMAX := MAXIMAL STEP SIZE
      IF (WORK(2).EQ.0.D0) THEN
         HMAX=TEND-T0
      ELSE
         HMAX=WORK(2)
         IF (HMAX.GT.TEND-T0) THEN
            HMAX=TEND-T0
         END IF
      END IF
C -------- THET  DECIDE WHETHER THE JACOBIAN SHOULD BE RECOMPUTED
      IF (WORK(3).EQ.0.D0) THEN
         THET = 0.005
      ELSE
         THET=WORK(3)
         IF (THET .GT. 1d0) THEN
            WRITE(msg,*)' CURIOUS INPUT WORK(3)=',WORK(3)
            CALL rwarn(msg)
            ARRET=.TRUE.
         END IF
      END IF
C--------- FACNEWT: STOPPING CRITERION FOR SPLITTING-NEWTON METHOD
C--------           FOR SMALL VALUES OF min(abs(y_i)) and min(abs(f_j))
      IF (WORK(4).EQ.0.D0) THEN
          FACNEWT=5d-3
          FACNEWT=DMAX1(FACNEWT,UROUND/RTOL(1) )
      ELSE
         FACNEWT=WORK(4)
         FACNEWT=DMAX1(FACNEWT,UROUND/RTOL(1) )
         IF (FACNEWT.GE.1.0D0) THEN
            WRITE(msg,*)'WRONG INPUT FOR WORK(4) ',WORK(4)
            CALL rwarn(msg)
            ARRET=.TRUE.
         END IF
      END IF
C--------- FACNORD(1): STOPPING CRITERION FOR SPLITTING-NEWTON METHOD
C---------             ORDER 3
      IF (WORK(15).EQ.0.D0) THEN
          FACNORD(1) = 1d-3
          FACNORD(1) = DMAX1(FACNORD(1) ,UROUND/RTOL(1) )
      ELSE
         FACNORD(1) = WORK(15)
         FACNORD(1) = DMAX1(FACNORD(1) ,UROUND/RTOL(1) )
         IF (FACNEWT.GE.1.0D0) THEN
            WRITE(msg,*)'WRONG INPUT FOR WORK(15) ',WORK(15)
            CALL rwarn(msg)
            ARRET=.TRUE.
         END IF
      END IF
C--------- FACNORD(2): STOPPING CRITERION FOR SPLITTING-NEWTON METHOD
C--------              ORDER 5
      IF (WORK(16).EQ.0.D0) THEN
          FACNORD(2) = 9d-2
          FACNORD(2) = DMAX1(FACNORD(2) ,UROUND/RTOL(1) )
      ELSE
         FACNORD(2) = WORK(16)
         FACNORD(2) = DMAX1(FACNORD(2) ,UROUND/RTOL(1) )
         IF (FACNORD(2).GE.1.0D0) THEN
            WRITE(msg,*)'WRONG INPUT FOR WORK(16) ',WORK(16)
            CALL rwarn(msg)
            ARRET=.TRUE.
         END IF
      END IF
C--------- FACNORD(3): STOPPING CRITERION FOR SPLITTING-NEWTON METHOD
C---------             ORDER 7
      IF (WORK(17).EQ.0.D0) THEN
          FACNORD(3) = 9d-1
          FACNORD(3) = DMAX1(FACNORD(3), UROUND/RTOL(1) )
      ELSE
         FACNORD(3) = WORK(17)
         FACNORD(3) = DMAX1(FACNORD(3), UROUND/RTOL(1) )
         IF (FACNORD(3).GE.1.0D0) THEN
            WRITE(msg,*)'WRONG INPUT FOR WORK(17) ',WORK(17)
            CALL rwarn(msg)
            ARRET=.TRUE.
         END IF
      END IF
C--------- FACNORD(4): STOPPING CRITERION FOR SPLITTING-NEWTON METHOD
C---------             ORDER 9
      IF (WORK(18).EQ.0.D0) THEN
          FACNORD(4) = 9.9d-1
          FACNORD(4) = DMAX1(FACNORD(4),UROUND/RTOL(1) )
      ELSE
         FACNORD(4) = WORK(18)
         FACNORD(4) = DMAX1(FACNORD(4),UROUND/RTOL(1) )
         IF (FACNORD(4).GE.1.0D0) THEN
            WRITE(msg,*)'WRONG INPUT FOR WORK(18) ',WORK(18)
            CALL rwarn(msg)
            ARRET=.TRUE.
         END IF
      END IF

C--------- TETAK0: STOPPING CRITERIUM FOR THE SPLITTING-NEWTON METHOD
C---------         THE ERROR IN THE ITERATES MUST BE DECREASING
C---------         BY A FACTOR TETAK0
      IF (WORK(5).EQ.0.D0) THEN
         TETAK0 = 0.9D0
      ELSE
         TETAK0 = WORK(5)
         IF (TETAK0.LE.0.0D0) THEN
            WRITE(msg,*)'WRONG INPUT FOR WORK(5) ',WORK(5)
            CALL rwarn(msg)
            ARRET=.TRUE.
         END IF
      END IF
      CS(1) = 1.0D0
C--------- CS(2): EMPIRICAL COMPUTATIONAL COST FOR ORDER 5
C---------        USED IN THE ORDER VARIATION STRATEGY.
      IF (WORK(6).EQ.0.D0) THEN
         CS(2) = 2.4D0
      ELSE
         CS(2) = WORK(6)
         IF (CS(2).LE.0.0D0) THEN
            WRITE(msg,*)'WRONG INPUT FOR WORK(6) ',WORK(6)
            CALL rwarn(msg)
            ARRET=.TRUE.
         END IF
      END IF
C---------  CS(3): EMPIRICAL COMPUTATIONAL COST FOR ORDER 7
C---------         USED IN THE ORDER VARIATION STRATEGY.
      IF (WORK(7).EQ.0.D0) THEN
         CS(3) = 4.0D0
      ELSE
         CS(3) = WORK(7)
         IF (CS(3).LE.0.0D0) THEN
            WRITE(msg,*)'WRONG INPUT FOR WORK(7) ',WORK(7)
            CALL rwarn(msg)
            ARRET=.TRUE.
         END IF
      END IF
C--------- CS(4): EMPIRICAL COMPUTATIONAL COST FOR ORDER 9
C---------        USED IN THE ORDER VARIATION STRATEGY.
      IF (WORK(8).EQ.0.D0) THEN
         CS(4) =7.2D0
      ELSE
         CS(4) = WORK(8)
         IF (CS(4).LE.0.0D0) THEN
            WRITE(msg,*)'WRONG INPUT FOR WORK(8) ',WORK(8)
            CALL rwarn(msg)
            ARRET=.TRUE.
         END IF
      END IF
C--------- FACL: PARAMETER FOR STEP SIZE SELECTION
C---------       THE NEW STEPSIZE IS CHOSEN SUBJECT TO THE RESTRICTION
C---------       FACL <= HNEW/HOLD
      IF (WORK(9).EQ.0.D0) THEN
         FACL = 0.12D0
      ELSE
         FACL = WORK(9)
         IF (FACL.LE.0.0D0) THEN
            WRITE(msg,*)'WRONG INPUT FOR WORK(9) ',WORK(9)
            CALL rwarn(msg)
            ARRET=.TRUE.
         END IF
      END IF
C--------- FACR: PARAMETER FOR STEP SIZE SELECTION
C---------       THE NEW STEPSIZE IS CHOSEN SUBJECT TO THE RESTRICTION
C---------       HNEW/HOLD <= FACR
      IF (WORK(10).EQ.0.D0) THEN
         FACR = 10D0
      ELSE
         FACR = WORK(10)
         IF (FACR.LE.0.0D0) THEN
            WRITE(msg,*)'WRONG INPUT FOR WORK(10) ',WORK(10)
            CALL rwarn(msg)
            ARRET=.TRUE.
         END IF
      END IF
C--------- SFDOWN: SAFETY FACTOR IN STEP SIZE PREDICTION
C---------         USED FOR THE LOWER ORDER METHOD
      IF (WORK(11).EQ.0.D0) THEN
         SFDOWN = 20.0D0
      ELSE
         SFDOWN = WORK(11)
         IF (SFDOWN.LE.0.0D0) THEN
            WRITE(msg,*)'WRONG INPUT FOR WORK(11) ',WORK(11)
            CALL rwarn(msg)
            ARRET=.TRUE.
         END IF
      END IF
C--------- SFUP:  SAFETY FACTOR IN STEP SIZE PREDICTION
C---------        USED FOR THE UPPER ORDER METHOD
      IF (WORK(12).EQ.0.D0) THEN
         SFUP = 40.0D0
      ELSE
         SFUP = WORK(12)
         IF (SFUP.LE.0.0D0) THEN
            WRITE(msg,*)'WRONG INPUT FOR WORK(12) ',WORK(12)
            CALL rwarn(msg)
            ARRET=.TRUE.
         END IF
      END IF
C--------- SFSAME: SAFETY FACTOR IN STEP SIZE PREDICTION
C---------         USED FOR THE CURRENT ORDER METHOD
      IF (WORK(13).EQ.0.D0) THEN
         SFSAME = 18.0D0
      ELSE
         SFSAME = WORK(13)
         IF (SFSAME.LE.0.0D0) THEN
            WRITE(msg,*)'WRONG INPUT FOR WORK(13) ',WORK(13)
            CALL rwarn(msg)
            ARRET=.TRUE.
         END IF
      END IF
C--------- SF: SAFETY FACTOR IN STEP SIZE PREDICTION
C---------     USED FOR THE CURRENT ORDER METHOD WHEN
C---------     THE ERROR CONTROL TEST fails
      IF (WORK(14).EQ.0.D0) THEN
         SF = 15.0D0
      ELSE
         SF = WORK(14)
         IF (SF.LE.0.0D0) THEN
            WRITE(msg,*)'WRONG INPUT FOR WORK(14) ',WORK(14)
            CALL rwarn(msg)
            ARRET=.TRUE.
         END IF
      END IF
C -------- CHECK  THE TOLERANCES
      IF (ITOL.EQ.0) THEN
          IF (ATOL(1).LE.0.D0.OR.RTOL(1).LE. UROUND) THEN
              WRITE (msg,*) ' TOLERANCES ARE TOO SMALL'
              CALL rwarn(msg)
              ARRET=.TRUE.
          END IF
      ELSE
          DO I=1,R
          IF (ATOL(I).LE.0.D0.OR.RTOL(I).LE. UROUND) THEN
              WRITE (msg,*) ' TOLERANCES(',I,') ARE TOO SMALL'
              CALL rwarn(msg)
              ARRET=.TRUE.
          END IF
          END DO
      END IF
C -------------------------------------------------------------------------
C     COMPUTATION OF ARRAY ENTRIES
C -------------------------------------------------------------------------
C -------- BANDED OR NOT 
      JBAND=MLJAC.LT.R
C -------- COMPUTATION OF THE ROW-DIMENSIONS OF THE 2-ARRAYS
      IF (JBAND) THEN
        LDJAC = MLJAC+MUJAC+1
        LDLU  = MLJAC+LDJAC
      ELSE
        LDJAC = R
        LDLU  = R
      END IF
      IF (JBAND) THEN
         IJOB=2
      ELSE
         IJOB=1
      END IF
C -------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK
      IEYP  = 19
      IEFP  = IEYP  + 10*R
      IEDN  = IEFP  + 10*R
      IEF   = IEDN  + R
      IEF1  = IEF   + 10*R
      IESC  = IEF1  + 10*R
      IEJF0 = IESC  + R
      IELU  = IEJF0 + R*LDJAC
C--------- TOTAL STORAGE REQUIREMENT
      ISTORE = IELU + R*LDLU - 1
      IF(ISTORE.GT.LWORK)THEN
        WRITE(msg,*)' INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=',ISTORE
        CALL rwarn(msg)
        ARRET=.TRUE.
      END IF
C -------- ENTRY POINTS FOR INTEGER WORKSPACE
      IEIPIV=25
C -------- TOTAL REQUIREMENT
      ISTORE=IEIPIV+R-1
      IF (ISTORE.GT.LIWORK) THEN
         WRITE(msg,*)' INSUFF. STORAGE FOR IWORK, MIN. LIWORK=',ISTORE
         CALL rwarn(msg)
         ARRET=.TRUE.
      END IF

C -------- WHEN A FAIL HAS OCCURED, GAM RETURNs WITH IDID=-1
      IF (ARRET) THEN
         IDID=-1
         RETURN
      END IF

      DO I=1,4
        NSTEP(I) = 0
        NNEWT(I) = 0
        NERR(I)  = 0
      END DO
      NDEC = 0
      NFCN = 0
      NJAC = 0


C -------------------------------------------------------------------------
C     CALL TO CORE INTEGRATOR
C -------------------------------------------------------------------------
      CALL ETRO(R, FCN, T0, Y0, TEND, HMAX, H, RTOL, ATOL, ITOL,
     &   JAC, IJAC, MLJAC, MUJAC, SOLOUT, IOUT, IDID, NMAX,
     &   UROUND, THET, FACNEWT, FACNORD, TETAK0, CS, FACL, FACR, SFDOWN,
     &   SFUP, SFSAME, SF, ORDMIN, ORDMAX, ITINT, ITMAX,
     &   JBAND, IJOB, LDJAC, LDLU, WORK(IEYP), WORK(IEFP),
     &   WORK(IEDN), WORK(IEF), WORK(IEF1), WORK(IESC),
     &   WORK(IEJF0), WORK(IELU), IWORK(IEIPIV),
     &   NFCN, NJAC, NSTEP, NNEWT, NERR, NDEC, RPAR, IPAR)

      IWORK(10)= NFCN
      IWORK(11)= NJAC
      IWORK(12)= NSTEP(1)
      IWORK(13)= NSTEP(2)
      IWORK(14)= NSTEP(3)
      IWORK(15)= NSTEP(4)
      IWORK(16)= NNEWT(1)
      IWORK(17)= NNEWT(2)
      IWORK(18)= NNEWT(3)
      IWORK(19)= NNEWT(4)
      IWORK(20)= NERR(1)
      IWORK(21)= NERR(2)
      IWORK(22)= NERR(3)
      IWORK(23)= NERR(4)
      IWORK(24)= NDEC

      RETURN
      END
C
C--------- END OF SUBROUTINE GAM
C
C -------------------------------------------------------------------------
C     SUBROUTINE  ETRO (Extended trapezoidal Rules of Odd order,
C                       that is GAMs)
C -------------------------------------------------------------------------
      SUBROUTINE  ETRO(R,FCN,T0,Y0,TEND,HMAX,H,RTOL,ATOL,ITOL,
     &   JAC,IJAC,MLJAC,MUJAC,SOLOUT,IOUT,IDID,
     &   NMAX,UROUND,THET,FACNEWT,FACNORD,TETAK0,CS,FACL,FACR,SFDOWN,
     &   SFUP,SFSAME,SF, ORDMIN,ORDMAX,ITINT,ITMAX,
     &   JBAND,IJOB,LDJAC,LDLU,YP,FP,
     &   DN,F,F1,SCAL, JF0, LU, IPIV,
     &   NFCN,NJAC,NSTEP,NNEWT,NERR,NDEC,RPAR,IPAR)
C -------------------------------------------------------------------------
C     CORE INTEGRATOR FOR GAM
C     PARAMETERS SAME AS IN GAM WITH ADDED WORKSPACE
C -------------------------------------------------------------------------
C     DECLARATIONS
C ----------------------------------------------------------
      IMPLICIT NONE
C
C   COMMON
C------------------------------------
      COMMON/LINAL/MLLU,MULU,MDIAG
C
C   INPUT VARIABLES
C------------------------------------
       INTEGER R, ORDMIN, ORDMAX, IDID, ITOL, IJAC, MLJAC, MUJAC, IOUT,
     &         IPAR(*), ITINT(4), ITMAX, IJOB, NMAX, LDJAC, LDLU
   
       DOUBLE PRECISION TEND, ATOL(*), RTOL(*), RPAR(*), FACNORD(4),
     &                  HMAX, THET, FACNEWT, TETAK0, CS(4), FACL, FACR,
     &                  SFDOWN, SFUP, SFSAME, SF, UROUND

C
C   OUTPUT VARIABLES
C------------------------------------
       INTEGER NDEC, NFCN, NJAC, NSTEP(4), NNEWT(4), NERR(4), IER
C
C   INPUT/OUTPUT VARIABLES
C---- 
       DOUBLE PRECISION T0, Y0(R), H, SCAL(R), YP(R,10), FP(R,10), 
C     &                  F(R,9),DN(R), F1(R,9), JF0(LDJAC,R), LU(LDLU,R)  CHANGED...
     &                F(R,10),DN(R), F1(R,10), JF0(LDJAC,R), LU(LDLU,R)
     
C
C   LOCAL VARIABLES
C------------------------------------
       INTEGER  I, J, IPIV(R), NSING,
     &          FAILNI, FAILEI, ORDOLD, ORD, ORD2, ORDN, 
     &          IT, DBL(4), DBLK, DBLKOLD,
     &          NSTEPS, IRTRN, NT1, MLLU, MULU, MDIAG

       DOUBLE PRECISION ERRV(10), TP(10), T1(10), YSAFE, DELT, 
     &                  THETA, TETAK, TETAKOLD,
     &                  HOLD, ERRNEWT, ERRNEWT0,  ERRNEWT1, ESP,
     &                  ERRUP, ERRSAME, ERRDOWN, RR, RRN, TH, THN, FACN
      LOGICAL  JBAND, CALJAC, NEWJAC, JVAI, TER, EXTRAP
C
C   EXTERNAL FUNCTIONS
C------------------------------------
       EXTERNAL FCN,JAC, SOLOUT
       CHARACTER*80 MSG

C
C   INCLUDE
C------------------------------------
       DOUBLE PRECISION 
     & B3511, B3512, B3513, B3514, B3515, B3521, B3522, B3523, B3524, 
     & B3525, B3531, B3532, B3533, B3534, B3535, B3541, B3542, B3543, 
     & B3544, B3545, 
     & L321, L331, L332, L341, L342, L343, B311, B312, B313, 
     & B511, B512, B513, B514, B515, B521, B522, B523, B524, B525, 
     & L521, L531, L541, L551, L561, L532, L542, L552, L562, L543, 
     & L553, L563, L554, L564, L565, 
     & B5711, B5712, B5713, B5714, B5715, B5716, B5717, B5721, B5722, 
     & B5723, B5724, B5725, B5726, B5727, B5731, B5732, B5733, B5734, 
     & B5735, B5736, B5737, B5741, B5742, B5743, B5744, B5745, B5746, 
     & B5747, 
     & B711, B712, B713, B714, B715, B716, B717, B721, B722, B723, 
     & B724, B725, B726, B727, B731, B732, B733, B734, B735, B736, B737, 
     & B7911, B7912, B7913, B7914, B7915, B7916, B7917, B7918, B7919, 
     & B7921, B7922, B7923, B7924, B7925, B7926, B7927, B7928, B7929, 
     & B7931, B7932, B7933, B7934, B7935, B7936, B7937, B7938, B7939, 
     & B7941, B7942, B7943, B7944, B7945, B7946, B7947, B7948, B7949, 
     & B7951, B7952, B7953, B7954, B7955, B7956, B7957, B7958, B7959, 
     & L721, L731, L741, L751, L761, L771, L781, L732, L742, L752, 
     & L762, L772, L782, L743, L753, L763, L773, L783, L754, L764, 
     & L774, L784, L765, L775, L785, L776, L786, L787, 
     & CP31, CP51, CP52, CP71, CP72, CP73, 
     & CP91, CP92, CP93, CP94, 
     & B911, B912, B913, B914, B915, B916, B917, B918, B919, B921, 
     & B922, B923, B924, B925, B926, B927, B928, B929, B931, B932, 
     & B933, B934, B935, B936, B937, B938, B939, B941, B942, B943, 
     & B944, B945, B946, B947, B948, B949, 
     & L921, L932, L943, L954, L965, L976, L987, L998, 
     & B91011, B91012, B91013, B91014, B91015, B91021, B91022, B91023, 
     & B91024, B91025, B91031, B91032, B91033, B91034, B91035, B91041, 
     & B91042, B91043, B91044, B91045 
       PARAMETER( B311 = 5d0/12d0, 
     &            B312 = 8d0/12d0, 
     &            B313 = -1d0/12d0) 
 
       PARAMETER( Cp31 = 1d0/24d0) 
ccccc L3 dimensione 4 
       PARAMETER( L321 =  4.012395124208693d-01, 
     &            L331 =  8.819910099032529d-04, 
     &            L341 =  1.728116022258560d-04, 
     &            L332 =  3.680857287181573d-01, 
     &            L342 =  1.635381132422046d-03, 
     &            L343 =  3.688541178419062d-01) 
c--------- B3-B5 OF DIMENSION  4 
       PARAMETER( B3511 = 49d0/720d0, 
     &            B3512 = -83d0/360d0, 
     &            B3513 = 17d0/60d0, 
     &            B3514 = -53d0/360d0, 
     &            B3515 = 19d0/720d0, 
     &            B3521 = 19d0/720d0, 
     &            B3522 = -23d0/360d0, 
     &            B3523 = 1d0/30d0, 
     & B3524 = 7d0/360d0, 
     & B3525 = -11d0/720d0, 
     & B3531 = -11d0/720d0, 
     & B3532 = 37d0/360d0, 
     & B3533 = -13d0/60d0, 
     & B3534 = 67d0/360d0, 
     & B3535 = -41d0/720d0, 
     & B3541 = 19d0/720d0, 
     & B3542 = -53d0/360d0, 
     & B3543 =  17d0/60d0, 
     & B3544 = -83d0/360d0, 
     & B3545 =  49d0/720d0) 
 
c[ 49/720, -83/360,  17/60, -53/360,  19/720] 
c[ 19/720, -23/360,   1/30,   7/360, -11/720] 
c[-11/720,  37/360, -13/60,  67/360, -41/720] 
c[ 19/720, -53/360,  17/60, -83/360,  49/720] 
 
C--------- A5, B5, B53, B56 Cp5 := MATRICES DEFINING  GAM5 
 
       PARAMETER( B511 = 251d0/720d0, 
     & B512 = 323d0/360d0, 
     & B513 = - 11d0/30d0, 
     & B514 = 53d0/360d0, 
     & B515 = -19d0/720d0, 
     & B521 = -19d0/720d0, 
     & B522 =  173d0/360d0, 
     & B523 = 19d0/30d0, 
     & B524 = -37d0/360d0, 
     & B525 = 11d0/720d0) 
 
       PARAMETER( B5711 = 1997d0/60480d0, 
     & B5712 = -113d0/630d0, 
     & B5713 = 1619d0/4032d0, 
     & B5714 = -715d0/1512d0, 
     & B5715 = 1241d0/4032d0, 
     & B5716 = -263d0/2520d0, 
     & B5717 = 863d0/60480d0, 
     & B5721 = -733d0/60480d0, 
     & B5722 = 41d0/630d0, 
     & B5723 = -193d0/1344d0, 
     & B5724 = 251d0/1512d0, 
     & B5725 = -425d0/4032d0, 
     & B5726 = 29d0/840d0, 
     & B5727 = -271d0/60480d0, 
     & B5731 = -271d0/60480d0, 
     & B5732 = 97d0/5040d0, 
     & B5733 = -13d0/448d0, 
     & B5734 = 5d0/378d0, 
     & B5735 = 37d0/4032d0, 
     & B5736 = -19d0/1680d0, 
     & B5737 = 191d0/60480d0, 
     & B5741 = 191d0/60480d0, 
     & B5742 = -67d0/2520d0, 
     & B5743 = 115d0/1344d0, 
     & B5744 = -211d0/1512d0, 
     & B5745 = 499d0/4032d0, 
     & B5746 = -2d0/35d0, 
     & B5747 = 653d0/60480d0) 
 
c--------- B5 - B7 
c 
c[1997/60480,  -113/630, 1619/4032, -715/1512, 1241/4032, -263/2520,  863/60480] 
c[-733/60480,    41/630, -193/1344,  251/1512, -425/4032,    29/840, -271/60480] 
c[-271/60480,   97/5040,   -13/448,     5/378,   37/4032,  -19/1680,  191/60480] 
c[ 191/60480,  -67/2520,  115/1344, -211/1512,  499/4032,     -2/35,  653/60480] 
c[-271/60480,    29/840, -425/4032,  251/1512, -193/1344,    41/630, -733/60480] 
c[ 863/60480, -263/2520, 1241/4032, -715/1512, 1619/4032,  -113/630, 1997/60480] 
c 
c 
       PARAMETER( L521 =  3.668340831928216D-01, 
     & L531 =  2.477905683677308D-03, 
     & L541 = -1.919925047010838D-03, 
     & L551 =  2.218385581234200D-03, 
     & L561 = -5.442189351609260D-03, 
     & L532 =  3.216639533696728D-01, 
     & L542 =  1.231925763308414D-03, 
     & L552 =  7.841944627374794D-03, 
     & L562 =  1.002485104590053D-03, 
     & L543 =  3.375100828961925D-01, 
     & L553 = -2.614300734741796D-04, 
     & L563 =  1.066631182323580D-03, 
     & L554 =  3.523137378783708D-01, 
     & L564 = -3.596681121610224D-04, 
     & L565 =  3.617716171655064D-01, 
     & CP51 =   3D0/160D0, 
     & CP52 = -11D0/1440D0) 
 
C--------- A7, B7, B57, B58 Cp7 := MATRICES DEFINING GAM7 
 
        PARAMETER( B711 = 19087d0/60480d0, 
     & B712 = 2713d0/2520d0, 
     & B713 = -15487d0/20160d0, 
     & B714 = 586d0/945d0, 
     & B715 = -6737d0/20160d0, 
     & B716 = 263d0/2520d0, 
     & B717 = -863d0/60480d0, 
     & B721 = -863d0/60480d0, 
     & B722 = 349d0/840d0, 
     & B723 = 5221d0/6720d0, 
     & B724 = -254d0/945d0, 
     & B725 = 811d0/6720d0, 
     & B726 = -29d0/840d0, 
     & B727 = 271d0/60480d0, 
     & B731 = 271d0/60480d0, 
     & B732 = -23d0/504d0, 
     & B733 = 10273d0/20160d0, 
     & B734 = 586d0/945d0, 
     & B735 = -2257d0/20160d0, 
     & B736 = 67d0/2520d0, 
     & B737 = -191d0/60480d0) 
 
C 
C--------- THE LAST THREE ROWS ARE THE REVERSE OF THE FIRST THREE 
C 
c B79 = 
c 
c[ 75203/3628800, -280187/1814400, 129781/259200, -238937/259200, 27289/25920, -197687/259200,  88531/259200, -156437/1814400,  33953/3628800] 
c[-17827/3628800,   66043/1814400, -30389/259200,   55513/259200, -6281/25920,   44983/259200, -19859/259200,   34453/1814400,  -7297/3628800] 
c[  8963/3628800,  -32987/1814400,  15061/259200,  -27257/259200,  3049/25920,  -21527/259200,   9331/259200,  -15797/1814400,   3233/3628800] 
c[  3233/3628800,  -10067/1814400,   3601/259200,   -4337/259200,     23/3240,    1393/259200,  -2129/259200,    7123/1814400,  -2497/3628800] 
c[ -2497/3628800,   12853/1814400,  -7859/259200,   18583/259200, -2681/25920,   24313/259200, -13589/259200,   30043/1814400,  -8227/3628800] 
c[  3233/3628800,  -15797/1814400,   9331/259200,  -21527/259200,  3049/25920,  -27257/259200,  15061/259200,  -32987/1814400,   8963/3628800] 
c[ -7297/3628800,   34453/1814400, -19859/259200,   44983/259200, -6281/25920,   55513/259200, -30389/259200,   66043/1814400, -17827/3628800] 
c[ 33953/3628800, -156437/1814400,  88531/259200, -197687/259200, 27289/25920, -238937/259200, 129781/259200, -280187/1814400,  75203/3628800] 
c 
       PARAMETER( B7911 =   75203d0/3628800d0, 
     & B7912 = -280187d0/1814400d0, 
     & B7913 =  129781d0/259200d0, 
     & B7914 = -238937d0/259200d0, 
     & B7915 =   27289d0/25920d0, 
     & B7916 = -197687d0/259200d0, 
     & B7917 =   88531d0/259200d0, 
     & B7918 = -156437d0/1814400d0, 
     & B7919 =   33953d0/3628800d0) 
 
       PARAMETER(B7921 = -17827d0/3628800d0, 
     & B7922 =  66043d0/1814400d0, 
     & B7923 = -30389d0/259200d0, 
     & B7924 =  55513d0/259200d0, 
     & B7925 =  -6281d0/25920d0, 
     & B7926 =  44983d0/259200d0, 
     & B7927 = -19859d0/259200d0, 
     & B7928 =  34453d0/1814400d0, 
     & B7929 =  -7297d0/3628800d0) 
 
       PARAMETER(B7931 =   8963d0/3628800d0, 
     & B7932 = -32987d0/1814400d0, 
     & B7933 =  15061d0/259200d0, 
     & B7934 = -27257d0/259200d0, 
     & B7935 =   3049d0/25920d0, 
     & B7936 = -21527d0/259200d0, 
     & B7937 =   9331d0/259200d0, 
     & B7938 = -15797d0/1814400d0, 
     & B7939 =   3233d0/3628800d0) 
 
       PARAMETER(B7941 =   3233d0/3628800d0, 
     & B7942 = -10067d0/1814400d0, 
     & B7943 =   3601d0/259200d0, 
     & B7944 =  -4337d0/259200d0, 
     & B7945 =     23d0/3240d0, 
     & B7946 =   1393d0/259200d0, 
     & B7947 =  -2129d0/259200d0, 
     & B7948 =   7123d0/1814400d0, 
     & B7949 =  -2497d0/3628800d0) 
 
 
       PARAMETER(B7951 =  -2497d0/3628800d0, 
     & B7952 =  12853d0/1814400d0, 
     & B7953 =  -7859d0/259200d0, 
     & B7954 =  18583d0/259200d0, 
     & B7955 =  -2681d0/25920d0, 
     & B7956 =  24313d0/259200d0, 
     & B7957 = -13589d0/259200d0, 
     & B7958 =  30043d0/1814400d0, 
     & B7959 =  -8227d0/3628800d0) 
 
       PARAMETER(L721 = 3.023839891568610D-01, 
     & L731 = 3.201698610574002D-05, 
     & L741 = 4.193101163680004D-04, 
     & L751 = 1.686924996069667D-04, 
     & L761 = 4.806043527549464D-05, 
     & L771 = 3.598347048026785D-06, 
     & L781 = 7.892534649789167D-04, 
     & L732 = 2.559868364091398D-01, 
     & L742 = 1.336896192287030D-04, 
     & L752 = 3.080994719931695D-03, 
     & L762 = 1.457177183563680D-04, 
     & L772 = 9.259360509484074D-04, 
     & L782 = 2.397658879381223D-04, 
     & L743 = 2.639734712170458D-01, 
     & L753 = 1.734338929611258D-04, 
     & L763 = 6.704398263264620D-03, 
     & L773 = 4.559927214651730D-05, 
     & L783 = 6.396418554053151D-05, 
     & L754 = 2.817729090368562D-01, 
     & L764 = 2.877761776030408D-04, 
     & L774 = 1.810919475521773D-04, 
     & L784 = 1.009049833235848D-03, 
     & L765 = 2.993040718034231D-01, 
     & L775 = 2.009850887505898D-03, 
     & L785 = 1.748065618845750D-03, 
     & L776 = 3.150349043479135D-01, 
     & L786 = 3.243816792609449D-05, 
     & L787 = 3.271307059448932D-01) 
 
       PARAMETER(CP71 = 103D0/9061D0, 
     & CP72 = -13D0/4480D0, 
     & CP73 =  67D0/42431D0) 
 
C--------- A8, B8, B86, B810 Cp8 := MATRICES DEFINING GAM9 
 
       PARAMETER(B911 = 1070017D0/3628800D0, 
     & B912 = 2233547D0/1814400D0, 
     & B913 = -2302297D0/1814400D0, 
     & B914 = 2797679D0/1814400D0, 
     & B915 = -31457D0/22680D0, 
     & B916 = 1573169D0/1814400D0, 
     & B917 = -645607D0/1814400D0, 
     & B918 = 156437D0/1814400D0, 
     & B919 = -33953D0/3628800D0, 
     & B921 = -33953D0/3628800D0, 
     & B922 = 687797D0/1814400D0, 
     & B923 =  1622393D0/1814400D0, 
     & B924 = -876271D0/1814400D0, 
     & B925 =   8233D0/22680D0, 
     & B926 =    -377521D0/1814400D0, 
     & B927 =   147143D0/1814400D0, 
     & B928 =  -34453D0/1814400D0, 
     & B929 =   7297D0/3628800D0, 
     & B931 = 7297D0/3628800D0, 
     & B932 =  -49813D0/1814400D0, 
     & B933 =  819143D0/1814400D0, 
     & B934 =  1315919D0/1814400D0, 
     & B935 = -5207D0/22680D0, 
     & B936 =  198929D0/1814400D0, 
     & B937 =  -71047D0/1814400D0, 
     & B938 =  15797D0/1814400D0, 
     & B939 = -3233D0/3628800D0, 
     & B941 = -3233D0/3628800D0, 
     & B942 = 18197D0/1814400D0, 
     & B943 =  -108007D0/1814400D0, 
     & B944 =  954929D0/1814400D0, 
     & B945 = 13903D0/22680D0, 
     & B946 = -212881D0/1814400D0, 
     & B947 =  63143D0/1814400D0, 
     & B948 =  -12853D0/1814400D0, 
     & B949 =  2497D0/3628800D0) 
 
 
c   B910=B9-B10; 
c               prime 4 righe : la 5 e' uguale alla 4; 
c     le ultime 4 uguali alle prime 4 negate (simmetriche); 
c     le prime  5 colonne  : le altre sono simmetriche e negate; 
c 
c[ 8183/1036800, -8183/115200,   8183/28800, -57281/86400,  57281/57600, -57281/57600,  57281/86400,  -8183/28800,  8183/115200, -8183/1036800] 
c[  -425/290304,    425/32256,    -425/8064,     425/3456,    -425/2304,     425/2304,    -425/3456,     425/8064,   -425/32256,    425/290304] 
c[      7/12800,    -63/12800,      63/3200,    -147/3200,     441/6400,    -441/6400,     147/3200,     -63/3200,     63/12800,      -7/12800] 
c[-2497/7257600,  2497/806400, -2497/201600,   2497/86400,  -2497/57600,   2497/57600,  -2497/86400,  2497/201600, -2497/806400,  2497/7257600] 
c[-2497/7257600,  2497/806400, -2497/201600,   2497/86400,  -2497/57600,   2497/57600,  -2497/86400,  2497/201600, -2497/806400,  2497/7257600] 
c[ 2497/7257600, -2497/806400,  2497/201600,  -2497/86400,   2497/57600,  -2497/57600,   2497/86400, -2497/201600,  2497/806400, -2497/7257600] 
c[     -7/12800,     63/12800,     -63/3200,     147/3200,    -441/6400,     441/6400,    -147/3200,      63/3200,    -63/12800,       7/12800] 
c[   425/290304,   -425/32256,     425/8064,    -425/3456,     425/2304,    -425/2304,     425/3456,    -425/8064,    425/32256,   -425/290304] 
c[-8183/1036800,  8183/115200,  -8183/28800,  57281/86400, -57281/57600,  57281/57600, -57281/86400,   8183/28800, -8183/115200,  8183/1036800] 
 
       PARAMETER(B91011 =   8183d0/1036800d0, 
     & B91012 =  -8183d0/115200d0, 
     & B91013 =   8183d0/28800d0, 
     & B91014 = -57281d0/86400d0, 
     & B91015 =  57281d0/57600d0) 
 
       PARAMETER(B91021 =   -425d0/290304d0, 
     & B91022 =    425d0/32256d0, 
     & B91023 =   -425d0/8064d0, 
     & B91024 =    425d0/3456d0, 
     & B91025 =   -425d0/2304d0) 
 
       PARAMETER(B91031 =     7d0/12800d0, 
     & B91032 =   -63d0/12800d0, 
     & B91033 =    63d0/3200d0, 
     & B91034 =  -147d0/3200d0, 
     & B91035 =   441d0/6400d0) 
 
       PARAMETER(B91041 = -2497d0/7257600d0, 
     & B91042 =  2497d0/806400d0, 
     & B91043 = -2497d0/201600d0, 
     & B91044 =  2497d0/86400d0, 
     & B91045 = -2497d0/57600d0) 
 
       PARAMETER(Cp91 =  7.892554012345216d-03, 
     & Cp92 = -1.463982583774219d-03, 
     & Cp93 =  5.468749999999983d-04, 
     & Cp94 = -3.440531305114634d-04) 
 
C--------- THE OTHERS ARE THE SAME WITH CHANGED SIGN 
 
       PARAMETER( L921 = 2.590721934790442d-01, 
     & L932   =   2.077575545359853d-01, 
     & L943   =   2.032874698558627d-01, 
     & L954   =   2.036384888660128d-01, 
     & L965   =   2.039599505779785d-01, 
     & L976   =   2.034044409161703d-01, 
     & L987   =   2.017245408702437d-01, 
     & L998   =   1.986549276295617d-01) 
 
c -------- CONSTANTS
      MLLU=MLJAC
      MULU=MUJAC
      MDIAG=MLLU + MULU +1
C--------- DBL(1:4) := SIZE OF THE COEFFICIENT MATRICES DEFINING THE GAMs

      DBL(1) = 4
      DBL(2) = 6
      DBL(3) = 8
      DBL(4) = 9
      ORD  = ORDMIN
      NSTEPS = 0
C--------- STARTING VALUES FOR NEWTON ITERATION
      DBLK = DBL(ORD)
      H    = MIN( H, ABS(TEND-T0)/DBLK )
      DBLKOLD = DBLK
      CALL FCN(R,T0,Y0,FP(1,1), RPAR,IPAR)
      NFCN = NFCN + 1
C -------- NUMBER OF FAILURES IN THE SPLITTING-NEWTON SCHEME
      FAILNI = 0
C -------- NUMBER OF FAILURES DUE TO THE ERROR TEST
      FAILEI = 0
      NSING  = 0
      ORDOLD = 2*ORD
      HOLD   = 2*H
      CALJAC = .TRUE.
      EXTRAP = .FALSE.
      ITMAX  = ITINT(ORD)
C--------- MAIN LOOP (ADVANCING IN TIME)
100   CONTINUE

C--------- (EVENTUALLY) COMPUTE THE JACOBIAN MATRIX NUMERICALLY
       NEWJAC = .FALSE.
       IF (CALJAC) THEN
          IF (IJAC.EQ.0) THEN
            DO I=1,R
               YSAFE=Y0(I)
               DELT=DSQRT(UROUND*DMAX1(1.D-5,DABS(YSAFE)))
               Y0(I)=YSAFE+DELT
               CALL FCN(R,T0,Y0,F,RPAR,IPAR)
               IF (JBAND) THEN
                 DO J=MAX(1,I-MUJAC),MIN(R,I+MLJAC)
                   JF0(J-I+MUJAC+1,I) = (F(J,1)-FP(J,1))/DELT
                 END DO
               ELSE
                 DO J=1,R
                   JF0(J,I)=(F(J,1)-FP(J,1))/DELT
                 END DO
               END IF
               Y0(I)=YSAFE
            END DO
C KS: added that
          NFCN = NFCN + R
         ELSE
C -------- COMPUTE JACOBIAN MATRIX ANALYTICALLY
           CALL JAC(R,T0,Y0,JF0,LDJAC,RPAR,IPAR)
        END IF
         NJAC = NJAC + 1
         NEWJAC = .TRUE.
       END IF

C--------- FACTORIZE THE ITERATION MATRIX
      IF ((ORDOLD.NE.ORD).OR.(HOLD.NE.H).OR.(NEWJAC) ) THEN
        HOLD = H
        ORDOLD = ORD
        IER = 1
        DO WHILE ( IER .NE. 0)
          CALL  DECLU(R,JF0,H,LDJAC,LU,LDLU,IPIV,ORD,IER,IJOB)
          IF (IER.NE.0) THEN
            NSING = NSING + 1
            IF (NSING.GT.5) THEN
              WRITE(msg,*) 'MATRIX IS REPEATEDLY SINGULAR, IER= ',IER
              CALL rwarn(msg)
              WRITE(msg,900) T0
              CALL rwarn(msg)
              IDID=-4
              GOTO 800
            ELSE
              H = H/2D0
            END IF
          END IF
          NDEC = NDEC + 1
        END DO
      END IF

C--------- DEFINE TP AND YP
      IF (EXTRAP) THEN
         T1(1) = T0+H
         DO I=2,DBLK+1
           T1(I) = T1(I-1)+H
         END DO
         CALL INTERP_GAM(R,TP,YP,T1,F1,NT1,DBLKOLD,DBLK,T0,Y0,ORD)
      ELSE
         TP(1) = T0
         DO J=1,R
           YP(J,1) = Y0(J)
         END DO
         DO I=2,DBLK+1
           DO J=1,R
             YP(J,I) = Y0(J)
           END DO
           TP(I) = TP(I-1)+H
         END DO
      END IF

C--------- DEFINE SCAL AND FACN

      THN = 1d0
      J    = 0
      IF (ITOL.EQ.0) THEN
          DO I=1,R
            RRN = DABS(Y0(I))
            SCAL(I)=ATOL(1)+RTOL(1)*RRN
            IF (RRN .LT. THN) THEN
               J = I
               THN = RRN
            ENDIF
          END DO
      ELSE
          DO I=1,R
            RRN = DABS( Y0(I) )
            SCAL(I)=ATOL(I)+RTOL(I)*RRN
            IF (RRN .LT. THN) THEN
               J = I
               THN = RRN
            ENDIF
          END DO
      END IF



       FACN = FACNORD(ORD)
       IF (THN .LT. 1d-1) THEN
          IF (DABS(FP(J,1)) .LT. 1d-5) THEN
              FACN = MIN(FACNEWT, FACNORD(ORD) )
          END IF
       END IF


C---------- COMPUTE THE NUMERICAL SOLUTION AT TIMES T1(1)...T1(DBLK)
C---------- DEFINE VARIABLES NEEDED IN THE ITERATION
        ERRNEWT  = FACN+1d0
        ERRNEWT0 = FACN+1d0
C----------
        TETAK    = 1.0D0
        THETA    = 1.0D0
        TETAKOLD = 1.0D0
        ITMAX    = ITINT(ORD)
        IT  = 0
        DO J = 2, DBLK+1
           CALL FCN(R,TP(J),YP(1,J),FP(1,J), RPAR,IPAR)
        END DO
        NFCN = NFCN + DBLK
C
C--------- SPLITTING NEWTON LOOP
C
 300    CONTINUE
          ERRNEWT1 = ERRNEWT0
          ERRNEWT0 = ERRNEWT
          ERRNEWT  = 0D0
C--------- COMPUTE ONE ITERATION FOR THE SELECTED ORDER
          GOTO (101,201,301,401) ORD
101     CALL      TERMNOT3(R,FCN,H,IT,DN, F,FP,YP,TP,NFCN,
     &  ERRNEWT,ERRNEWT0,TETAK0,LU, LDLU,IPIV, SCAL,IJOB,TER,
     &  RPAR,IPAR)
        GOTO 501
201     CALL      TERMNOT5(R,FCN,H,IT,DN, F,FP,YP,TP,NFCN,
     &  ERRNEWT,ERRNEWT0,TETAK0,LU, LDLU,IPIV, SCAL,IJOB,TER,
     &  RPAR,IPAR)
        GOTO 501
301     CALL      TERMNOT7(R,FCN,H,IT,DN, F,FP,YP,TP,NFCN,
     &  ERRNEWT,ERRNEWT0,TETAK0,LU, LDLU,IPIV, SCAL,IJOB,TER,
     &  RPAR,IPAR)
        GOTO 501
401     CALL      TERMNOT9(R,FCN,H,IT,DN, F,FP,YP,TP,NFCN,
     &  ERRNEWT,ERRNEWT0,TETAK0,LU, LDLU,IPIV, SCAL,IJOB,TER,
     &  RPAR,IPAR)
501     CONTINUE
         IF (TER)  THEN
           ERRNEWT = FACN + 1
           GOTO 999
         END IF
C--------- COMPUTE TETAK, ETAK

           TETAKOLD=TETAK
           TETAK   = ERRNEWT/SQRT(ERRNEWT0*ERRNEWT1)
           IF (IT.LE.2) THEN
              THETA=THET/2d0
           ELSE IF (IT .GT. 2) THEN
              THETA = SQRT(TETAK*TETAKOLD)
           END IF

           IT = IT+1

           JVAI = (IT .LE. ITMAX).AND.(ERRNEWT.GT.FACN).AND.
     &((THETA.LT.TETAK0).OR.(IT.LE.2)).AND. (ERRNEWT.GT.0d0)
           IF (JVAI) GO TO 300
C
C--------- END OF NEWTON LOOP
C
999   CONTINUE
      IF ((ERRNEWT.GT.FACN).OR.(.not.(ERRNEWT.GT.0d0))) THEN
C--------- THE ITERATION DOES NOT CONVERGE
         FAILNI = FAILNI + 1
         NNEWT(ORD) = NNEWT(ORD)+1
C--------- CHOICE OF THE NEW STEPSIZE
         H=H/2d0
         DBLKOLD = DBLK
         EXTRAP = .FALSE.
         IF (FAILNI .EQ. 1) THEN
            CALJAC = .NOT. NEWJAC
         ELSE
            CALJAC = .FALSE.
         END IF
C--------- RETURN TO THE MAIN LOOP
       ELSE
C--------- THE ITERATION CONVERGES
C--------- ERROR ESTIMATION
        CALL  ESTERR(ERRV, ERRSAME, ERRUP, ERRDOWN, FP,
     &     R, H, ORD, DBLK, LU, LDLU,
     &     IPIV, F, F1, SCAL, ORDMAX,ORDMIN,IJOB)
         IF ( ERRSAME .GT. 1D0 ) THEN
           FAILEI = FAILEI + 1
           NERR(ORD) = NERR(ORD) + 1
           CALJAC = (THETA .GT. THET)
C--------- NEW STEPSIZE SELECTION
           ORD2 = 2*ORD
           ESP = 1D0/(ORD2+1D0)
           RRN=DMAX1(FACL,DMIN1(FACR,(SF*ERRSAME)**ESP))
           H = H/RRN
           DBLKOLD = DBLK
           DO I=1,DBLKOLD+1
            DO J=1,R
              F1(J,I) = YP(J,I)
            END DO
           END DO
           call DIFFDIV(TP,F1,R,DBLK,NT1)
           EXTRAP = .TRUE.

C--------- RETURN TO THE MAIN LOOP
         ELSE
C--------- THE STEPSIZE IS ACCEPTED
           NSTEP(ORD) = NSTEP(ORD)+1
c          write(55,*) H
           T0 = TP(DBLK+1)
           DO I=1, R
             Y0(I) = YP(I,DBLK+1)
             FP(I,1) = FP(I,DBLK+1)
           END DO
C--------- NEW STEPSIZE SELECTION
            ORD2 = 2*ORD
            ESP = 1D0/(ORD2+1d0)
            RRN=DMAX1(FACL,DMIN1(FACR,(SFSAME*ERRSAME)**ESP))
            THN=DBL(ORD)/(CS(ORD)*RRN)
            ORDN = ORD
            IF  (ORD.LT.ORDMAX) THEN
              ESP = 1D0/(ORD2+3D0)
              RR=DMAX1(FACL,DMIN1(FACR,(SFUP*ERRUP)**ESP))
              TH=DBL(ORD+1)/(CS(ORD+1)*RR )
              IF (TH .GT. THN ) THEN
                ORDN = ORD + 1
                RRN  = RR
                THN  = TH
              END IF
            END IF

            IF ( ORD.GT.ORDMIN)  THEN
             ESP = 1D0/(ORD2-1d0)
             RR=DMAX1(FACL,DMIN1(FACR,(SFDOWN*ERRDOWN)**ESP))
             TH=DBL(ORD-1)/(CS(ORD-1)*RR )
             IF (TH .GT. THN ) THEN
                ORDN = ORD - 1
                RRN  = RR
             END IF
            END IF
           HOLD = H
c
c
           IF (ORDN.GT.ORD) THEN
               H = MIN(H/RRN, HOLD)
           ELSE
               H = H/RRN
           END IF
           ORDOLD = ORD
           ORD = ORDN
           DBLKOLD = DBLK
           DBLK = DBL(ORD)

           CALJAC = (THETA .GT. THET)
           EXTRAP = .TRUE.

           IF ((FAILNI.NE.0).OR.(FAILEI.NE.0)) THEN
             H = DMIN1( H, HOLD)
           END IF
           IF  (.NOT. CALJAC) THEN
             IF ((H/HOLD.LE.1.1D0 ).AND.(H/HOLD.GE.0.9D0)) THEN
                H = HOLD
             END IF
           END IF
           H = DMIN1( H, DMIN1(HMAX, (TEND-T0)/DBLK) )


           DO I=1,DBLKOLD+1
            DO J=1,R
              F1(J,I) = YP(J,I)
            END DO
           END DO
           CALL DIFFDIV(TP,F1,R,DBLKOLD,NT1)
           EXTRAP = .TRUE.
c
        IF (IOUT.NE.0) THEN
C--------- CALL SOLOUT
          CALL SOLOUT(R,TP,YP,F1,NT1,DBLKOLD,ORDOLD,RPAR,IPAR,IRTRN)
          IF (IRTRN.LT.0) GOTO 800
        END IF
          IF (NSTEPS .EQ. 0) THEN
            FAILNI = 0
            FAILEI = 0
          ELSE
            FAILNI = MAX(FAILNI-1,0)
            FAILEI = MAX(FAILEI-1,0)
          END IF
          NSING  = 0
        END IF
c--------- END IF ERRSAME > 1
      END IF
c--------- END IF ERRNEWT > 1
           NSTEPS = NSTEPS + 1
           IF (0.1d0*DABS(T0-TEND)/DBLK .GT. dabs(T0)*UROUND ) THEN
            IF (0.1d0*DABS(H).LE.DABS(T0)*UROUND) THEN
              WRITE(msg,*) ' STEPSIZE TOO SMALL, H=',H
              CALL rwarn(msg)
              WRITE(msg,900) T0
              CALL rwarn(msg)
              IDID=-3
              GOTO 800
            END IF
            IF (NSTEPS.GT.NMAX) THEN
              WRITE(msg,*) ' MORE THAN NMAX =',NMAX,'STEPS ARE NEEDED'
              CALL rwarn(msg)
              WRITE(msg,900) T0
              CALL rwarn(msg)
              IDID=-2
              GOTO 800
            END IF

           GOTO 100
C
C---------- END WHILE T0 < T
           ELSE
             H    = HOLD
             IDID = 1
           END IF
 900  FORMAT(' EXIT OF GAM AT T=',E18.4)
 800  RETURN
      END




c-----------------------------------------------------------------------
c     ADDITIONAL LINEAR ALGEBRA ROUTINES REQUIRED BY GAM
c-----------------------------------------------------------------------
C     VERSION OF AUGUST 20, 1997
c-----------------------------------------------------------------------
C
      SUBROUTINE DECLU(R,JF0,H,LDJAC,LU,LDLU,IPIV,ORD,IER,IJOB)
      IMPLICIT NONE
C
C   COMMON
C------------------------------------
      COMMON/LINAL/MLLU,MULU,MDIAG
C
C   INPUT VARIABLES
C------------------------------------
      INTEGER R, LDJAC, LDLU, ORD, MLLU, MULU, MDIAG, IJOB
      DOUBLE PRECISION  JF0(LDJAC,*), H
C
C   OUTPUT VARIABLES
C------------------------------------
       INTEGER IER, IPIV(R)
       DOUBLE PRECISION LU(LDLU,*)
C
C   LOCAL VARIABLES
C------------------------------------
      INTEGER I,J
      DOUBLE PRECISION  FAC, L31, L51, L71, L91
      PARAMETER(L31  =  6.411501944628007d-01,
     &          L51  =  6.743555662880509D-01,
     &          L71  =  7.109158294404152D-01,
     &          L91  =  7.440547954061898d-01)
C
C   EXECUTABLE STATEMENTS
C---------------------------------
C
      GOTO (10,20,30,40), ORD
  10   FAC = -(L31*H)
      GOTO 50
  20   FAC = -(L51*H)
      GOTO 50
  30   FAC = -(L71*H)
      GOTO 50
  40   FAC = -(L91*H)
  50  CONTINUE

      GO TO (1,2) IJOB

  1   CONTINUE

C -------- JACOBIAN A FULL MATRIX

      DO J=1,R
         DO  I=1,R
            LU(I,J)= FAC*JF0(I,J)
         END DO
         LU(J,J)=LU(J,J)+1d0
      END DO
      CALL DEC_GAM (R,LDLU,LU,IPIV,IER)
      RETURN

  2   CONTINUE

C -------- JACOBIAN A BAND MATRIX

      DO J=1,R
         DO I=1,MDIAG
            LU(I+MLLU,J)= FAC*JF0(I,J)
         END DO
         LU(MDIAG,J)=LU(MDIAG,J)+1d0
      END DO
      CALL DECB_gam (R,LDLU,LU,MLLU,MULU,IPIV,IER)
      RETURN

      END
C
C  SUBROUTINE SOLLU
C
      SUBROUTINE SOLLU(R,LU,LDLU,F,IPIV,IJOB)
      IMPLICIT NONE
C
C   COMMON
C------------------------------------
      COMMON/LINAL/MLLU,MULU,MDIAG
C
C   INPUT VARIABLES
C------------------------------------
      INTEGER R, LDLU, IPIV(R), MLLU, MULU, MDIAG, IJOB
      DOUBLE PRECISION  LU(LDLU,*)
C
C   INPUT/OUTPUT VARIABLES
C------------------------------------
       DOUBLE PRECISION F(R)
C
C   EXECUTABLE STATEMENTS
C---------------------------------
C

      GO TO (1,2) IJOB

  1   CONTINUE

C -------- JACOBIAN A FULL MATRIX

        CALL SOL_GAM (R,LDLU,LU,F,IPIV)
      RETURN

  2   CONTINUE


C -------- JACOBIAN A BAND MATRIX

      CALL SOLB_gam (R,LDLU,LU,MLLU,MULU,F,IPIV)
      RETURN

      END

C
C  SUBROUTINE NEWTGS
C
      SUBROUTINE NEWTGS(R,DBLK,LU,LDLU,IPIV,F,DN,IJOB)
      IMPLICIT NONE
C
C   INPUT VARIABLES
C------------------------------------
      INTEGER R, LDLU, IPIV(R), IJOB, DBLK
      DOUBLE PRECISION  LU(LDLU,*), F(R,DBLK)
C
C   OUTPUT VARIABLES
C------------------------------------
       DOUBLE PRECISION DN(R,DBLK)
C
C   LOCAL VARIABLES
C------------------------------------
      INTEGER I, J
C
C   EXECUTABLE STATEMENTS
C---------------------------------
C
         DO I=1, R
           DN(I,1) = -F(I,1)
         END DO
        CALL  SOLLU(R,LU,LDLU,DN(1,1),IPIV,IJOB)
        DO J=2,DBLK
           DO I=1, R
             DN(I,J) =  -F(I,J)+DN(I,J-1)
           END DO
           CALL  SOLLU(R,LU,LDLU,DN(1,J),IPIV,IJOB)
        END DO
        RETURN
        END
C
C    SUBROUTINE INTERP_GAM
C
      SUBROUTINE INTERP_GAM(R,TP,YP,T1,F1,NT1,DBLKOLD,DBLK,T0,Y0,ORD)
      IMPLICIT NONE
C
C   INPUT VARIABLES
C------------------------------------
      INTEGER R,  DBLK, DBLKOLD, ORD, NT1
      DOUBLE PRECISION  T0, Y0(R), T1(*), F1(R,*)
C
C   OUTPUT VARIABLES
C------------------------------------
       DOUBLE PRECISION TP(*), YP(R,*)
C
C   LOCAL VARIABLES
C------------------------------------
      INTEGER I, J, N, IT1, NT2
C
C   EXECUTABLE STATEMENTS
C---------------------------------
C
      NT2 = NT1+1
      N = DBLKOLD+1


      DO IT1=2,DBLK+1
        DO I=1,R
          YP(I,IT1) = F1(I,NT1)
          DO J=NT2,N
            YP(I,IT1) = YP(I,IT1)*(T1(IT1-1)-TP(J)) + F1(I,J)
          END DO
        END DO
      END DO

      DO J=1,R
         YP(J,1) = Y0(J)
      END DO

      GOTO (10,20,30,40) ORD
10    CONTINUE
        TP(1) = T0
        TP(2) = T1(1)
        TP(3) = T1(2)
        TP(4) = T1(3)
        TP(5) = T1(4)
      RETURN
20    CONTINUE
        TP(1) = T0
        TP(2) = T1(1)
        TP(3) = T1(2)
        TP(4) = T1(3)
        TP(5) = T1(4)
        TP(6) = T1(5)
        TP(7) = T1(6)
      RETURN
30    CONTINUE
        TP(1) = T0
        TP(2) = T1(1)
        TP(3) = T1(2)
        TP(4) = T1(3)
        TP(5) = T1(4)
        TP(6) = T1(5)
        TP(7) = T1(6)
        TP(8) = T1(7)
        TP(9) = T1(8)
      RETURN
40    CONTINUE
        TP(1) = T0
        TP(2) = T1(1)
        TP(3) = T1(2)
        TP(4) = T1(3)
        TP(5) = T1(4)
        TP(6) = T1(5)
        TP(7) = T1(6)
        TP(8) = T1(7)
        TP(9) = T1(8)
        TP(10) = T1(9)
      RETURN
      END
C
C    SUBROUTINE DIFFDIV
C

      SUBROUTINE DIFFDIV(TP,YP,R,DBLK,NT1)
      IMPLICIT NONE
C
C   INPUT VARIABLES
C------------------------------------
      INTEGER R,  DBLK
C
C   INPUT/OUTPUT VARIABLES
C------------------------------------
       INTEGER NT1
       DOUBLE PRECISION TP(1), YP(R,1)
C
C   LOCAL VARIABLES
C------------------------------------
      INTEGER I, J, K, N
C
C   EXECUTABLE STATEMENTS
C---------------------------------
C

      N = DBLK+1
      NT1 = 3


      DO J=N-1,NT1,-1
        DO K=1,J
         DO I=1,R
           YP(I,K)= ( YP(I,K)- YP(I,K+1) )/( TP(K)-TP(K+N-J))
         END DO
        END DO
      END DO
      RETURN
      END
C
C Karline: made subroutine from CONTR called CONTOUT
C

      SUBROUTINE contout(R,T,TP,FF,DBLK,NT1, CONTR)
C ----------------------------------------------------------
C     THIS FUNCTION CAN BE USED FOR CONTINUOUS OUTPUT. IT PROVIDES AN
C     APPROXIMATION TO THE I-TH COMPONENT OF THE SOLUTION AT T.
C     IT GIVES THE VALUE OF THE INTERPOLATION POLYNOMIAL, DEFINED FOR
C     THE LAST SUCCESSFULLY COMPUTED STEP.
C ----------------------------------------------------------
      IMPLICIT NONE
C
C   INPUT VARIABLES
C------------------------------------
      INTEGER R, I, DBLK, NT1
      DOUBLE PRECISION T, TP(*), FF(R,*), CONTR(*)
C
C   LOCAL VARIABLES
C------------------------------------
      INTEGER J, N, NT2
      DOUBLE PRECISION YP
C
C   EXECUTABLE STATEMENTS
C---------------------------------
C
      N = DBLK+1
      NT2=NT1+1
      DO I = 1, R
	      YP = FF(I,NT1)
        DO J=NT2,N
           YP = YP*(T-TP(J)) + FF(I,J)
        END DO
        CONTR(I) = YP
      ENDDO

      RETURN
      END

C
C     FUNCTION  CONTR
C
      DOUBLE PRECISION FUNCTION CONTR(I,R,T,TP,FF,DBLK,NT1)
C ----------------------------------------------------------
C     THIS FUNCTION CAN BE USED FOR CONTINUOUS OUTPUT. IT PROVIDES AN
C     APPROXIMATION TO THE I-TH COMPONENT OF THE SOLUTION AT T.
C     IT GIVES THE VALUE OF THE INTERPOLATION POLYNOMIAL, DEFINED FOR
C     THE LAST SUCCESSFULLY COMPUTED STEP.
C ----------------------------------------------------------
      IMPLICIT NONE
C
C   INPUT VARIABLES
C------------------------------------
      INTEGER R, I, DBLK, NT1
      DOUBLE PRECISION T, TP(*), FF(R,*)
C
C   LOCAL VARIABLES
C------------------------------------
      INTEGER J, N, NT2
      DOUBLE PRECISION YP
C
C   EXECUTABLE STATEMENTS
C---------------------------------
C
      N = DBLK+1
      NT2=NT1+1
        YP = FF(I,NT1)
        DO J=NT2,N
           YP = YP*(T-TP(J)) + FF(I,J)
        END DO
        CONTR = YP
      RETURN
      END

C
C  SUBROUTINE TERMNOT3  (ORDER 3)
C
      SUBROUTINE  TERMNOT3(R,FCN,H,IT,DN, F1,FP,YP,TP,NFCN,
     &  ERRNEWT,ERRNEWT0,TETAK0,LU, LDLU,IPIV, SCAL,IJOB,TER,
     &  RPAR,IPAR)
       IMPLICIT NONE
C
C   INCLUDE
C------------------------------------
       DOUBLE PRECISION 
     & B3511, B3512, B3513, B3514, B3515, B3521, B3522, B3523, B3524, 
     & B3525, B3531, B3532, B3533, B3534, B3535, B3541, B3542, B3543, 
     & B3544, B3545, 
     & L321, L331, L332, L341, L342, L343, B311, B312, B313, 
     & B511, B512, B513, B514, B515, B521, B522, B523, B524, B525, 
     & L521, L531, L541, L551, L561, L532, L542, L552, L562, L543, 
     & L553, L563, L554, L564, L565, 
     & B5711, B5712, B5713, B5714, B5715, B5716, B5717, B5721, B5722, 
     & B5723, B5724, B5725, B5726, B5727, B5731, B5732, B5733, B5734, 
     & B5735, B5736, B5737, B5741, B5742, B5743, B5744, B5745, B5746, 
     & B5747, 
     & B711, B712, B713, B714, B715, B716, B717, B721, B722, B723, 
     & B724, B725, B726, B727, B731, B732, B733, B734, B735, B736, B737, 
     & B7911, B7912, B7913, B7914, B7915, B7916, B7917, B7918, B7919, 
     & B7921, B7922, B7923, B7924, B7925, B7926, B7927, B7928, B7929, 
     & B7931, B7932, B7933, B7934, B7935, B7936, B7937, B7938, B7939, 
     & B7941, B7942, B7943, B7944, B7945, B7946, B7947, B7948, B7949, 
     & B7951, B7952, B7953, B7954, B7955, B7956, B7957, B7958, B7959, 
     & L721, L731, L741, L751, L761, L771, L781, L732, L742, L752, 
     & L762, L772, L782, L743, L753, L763, L773, L783, L754, L764, 
     & L774, L784, L765, L775, L785, L776, L786, L787, 
     & CP31, CP51, CP52, CP71, CP72, CP73, 
     & CP91, CP92, CP93, CP94, 
     & B911, B912, B913, B914, B915, B916, B917, B918, B919, B921, 
     & B922, B923, B924, B925, B926, B927, B928, B929, B931, B932, 
     & B933, B934, B935, B936, B937, B938, B939, B941, B942, B943, 
     & B944, B945, B946, B947, B948, B949, 
     & L921, L932, L943, L954, L965, L976, L987, L998, 
     & B91011, B91012, B91013, B91014, B91015, B91021, B91022, B91023, 
     & B91024, B91025, B91031, B91032, B91033, B91034, B91035, B91041, 
     & B91042, B91043, B91044, B91045 
       PARAMETER( B311 = 5d0/12d0, 
     &            B312 = 8d0/12d0, 
     &            B313 = -1d0/12d0) 
 
       PARAMETER( Cp31 = 1d0/24d0) 
ccccc L3 dimensione 4 
       PARAMETER( L321 =  4.012395124208693d-01, 
     &            L331 =  8.819910099032529d-04, 
     &            L341 =  1.728116022258560d-04, 
     &            L332 =  3.680857287181573d-01, 
     &            L342 =  1.635381132422046d-03, 
     &            L343 =  3.688541178419062d-01) 
c--------- B3-B5 OF DIMENSION  4 
       PARAMETER( B3511 = 49d0/720d0, 
     &            B3512 = -83d0/360d0, 
     &            B3513 = 17d0/60d0, 
     &            B3514 = -53d0/360d0, 
     &            B3515 = 19d0/720d0, 
     &            B3521 = 19d0/720d0, 
     &            B3522 = -23d0/360d0, 
     &            B3523 = 1d0/30d0, 
     & B3524 = 7d0/360d0, 
     & B3525 = -11d0/720d0, 
     & B3531 = -11d0/720d0, 
     & B3532 = 37d0/360d0, 
     & B3533 = -13d0/60d0, 
     & B3534 = 67d0/360d0, 
     & B3535 = -41d0/720d0, 
     & B3541 = 19d0/720d0, 
     & B3542 = -53d0/360d0, 
     & B3543 =  17d0/60d0, 
     & B3544 = -83d0/360d0, 
     & B3545 =  49d0/720d0) 
 
c[ 49/720, -83/360,  17/60, -53/360,  19/720] 
c[ 19/720, -23/360,   1/30,   7/360, -11/720] 
c[-11/720,  37/360, -13/60,  67/360, -41/720] 
c[ 19/720, -53/360,  17/60, -83/360,  49/720] 
 
C--------- A5, B5, B53, B56 Cp5 := MATRICES DEFINING  GAM5 
 
       PARAMETER( B511 = 251d0/720d0, 
     & B512 = 323d0/360d0, 
     & B513 = - 11d0/30d0, 
     & B514 = 53d0/360d0, 
     & B515 = -19d0/720d0, 
     & B521 = -19d0/720d0, 
     & B522 =  173d0/360d0, 
     & B523 = 19d0/30d0, 
     & B524 = -37d0/360d0, 
     & B525 = 11d0/720d0) 
 
       PARAMETER( B5711 = 1997d0/60480d0, 
     & B5712 = -113d0/630d0, 
     & B5713 = 1619d0/4032d0, 
     & B5714 = -715d0/1512d0, 
     & B5715 = 1241d0/4032d0, 
     & B5716 = -263d0/2520d0, 
     & B5717 = 863d0/60480d0, 
     & B5721 = -733d0/60480d0, 
     & B5722 = 41d0/630d0, 
     & B5723 = -193d0/1344d0, 
     & B5724 = 251d0/1512d0, 
     & B5725 = -425d0/4032d0, 
     & B5726 = 29d0/840d0, 
     & B5727 = -271d0/60480d0, 
     & B5731 = -271d0/60480d0, 
     & B5732 = 97d0/5040d0, 
     & B5733 = -13d0/448d0, 
     & B5734 = 5d0/378d0, 
     & B5735 = 37d0/4032d0, 
     & B5736 = -19d0/1680d0, 
     & B5737 = 191d0/60480d0, 
     & B5741 = 191d0/60480d0, 
     & B5742 = -67d0/2520d0, 
     & B5743 = 115d0/1344d0, 
     & B5744 = -211d0/1512d0, 
     & B5745 = 499d0/4032d0, 
     & B5746 = -2d0/35d0, 
     & B5747 = 653d0/60480d0) 
 
c--------- B5 - B7 
c 
c[1997/60480,  -113/630, 1619/4032, -715/1512, 1241/4032, -263/2520,  863/60480] 
c[-733/60480,    41/630, -193/1344,  251/1512, -425/4032,    29/840, -271/60480] 
c[-271/60480,   97/5040,   -13/448,     5/378,   37/4032,  -19/1680,  191/60480] 
c[ 191/60480,  -67/2520,  115/1344, -211/1512,  499/4032,     -2/35,  653/60480] 
c[-271/60480,    29/840, -425/4032,  251/1512, -193/1344,    41/630, -733/60480] 
c[ 863/60480, -263/2520, 1241/4032, -715/1512, 1619/4032,  -113/630, 1997/60480] 
c 
c 
       PARAMETER( L521 =  3.668340831928216D-01, 
     & L531 =  2.477905683677308D-03, 
     & L541 = -1.919925047010838D-03, 
     & L551 =  2.218385581234200D-03, 
     & L561 = -5.442189351609260D-03, 
     & L532 =  3.216639533696728D-01, 
     & L542 =  1.231925763308414D-03, 
     & L552 =  7.841944627374794D-03, 
     & L562 =  1.002485104590053D-03, 
     & L543 =  3.375100828961925D-01, 
     & L553 = -2.614300734741796D-04, 
     & L563 =  1.066631182323580D-03, 
     & L554 =  3.523137378783708D-01, 
     & L564 = -3.596681121610224D-04, 
     & L565 =  3.617716171655064D-01, 
     & CP51 =   3D0/160D0, 
     & CP52 = -11D0/1440D0) 
 
C--------- A7, B7, B57, B58 Cp7 := MATRICES DEFINING GAM7 
 
        PARAMETER( B711 = 19087d0/60480d0, 
     & B712 = 2713d0/2520d0, 
     & B713 = -15487d0/20160d0, 
     & B714 = 586d0/945d0, 
     & B715 = -6737d0/20160d0, 
     & B716 = 263d0/2520d0, 
     & B717 = -863d0/60480d0, 
     & B721 = -863d0/60480d0, 
     & B722 = 349d0/840d0, 
     & B723 = 5221d0/6720d0, 
     & B724 = -254d0/945d0, 
     & B725 = 811d0/6720d0, 
     & B726 = -29d0/840d0, 
     & B727 = 271d0/60480d0, 
     & B731 = 271d0/60480d0, 
     & B732 = -23d0/504d0, 
     & B733 = 10273d0/20160d0, 
     & B734 = 586d0/945d0, 
     & B735 = -2257d0/20160d0, 
     & B736 = 67d0/2520d0, 
     & B737 = -191d0/60480d0) 
 
C 
C--------- THE LAST THREE ROWS ARE THE REVERSE OF THE FIRST THREE 
C 
c B79 = 
c 
c[ 75203/3628800, -280187/1814400, 129781/259200, -238937/259200, 27289/25920, -197687/259200,  88531/259200, -156437/1814400,  33953/3628800] 
c[-17827/3628800,   66043/1814400, -30389/259200,   55513/259200, -6281/25920,   44983/259200, -19859/259200,   34453/1814400,  -7297/3628800] 
c[  8963/3628800,  -32987/1814400,  15061/259200,  -27257/259200,  3049/25920,  -21527/259200,   9331/259200,  -15797/1814400,   3233/3628800] 
c[  3233/3628800,  -10067/1814400,   3601/259200,   -4337/259200,     23/3240,    1393/259200,  -2129/259200,    7123/1814400,  -2497/3628800] 
c[ -2497/3628800,   12853/1814400,  -7859/259200,   18583/259200, -2681/25920,   24313/259200, -13589/259200,   30043/1814400,  -8227/3628800] 
c[  3233/3628800,  -15797/1814400,   9331/259200,  -21527/259200,  3049/25920,  -27257/259200,  15061/259200,  -32987/1814400,   8963/3628800] 
c[ -7297/3628800,   34453/1814400, -19859/259200,   44983/259200, -6281/25920,   55513/259200, -30389/259200,   66043/1814400, -17827/3628800] 
c[ 33953/3628800, -156437/1814400,  88531/259200, -197687/259200, 27289/25920, -238937/259200, 129781/259200, -280187/1814400,  75203/3628800] 
c 
       PARAMETER( B7911 =   75203d0/3628800d0, 
     & B7912 = -280187d0/1814400d0, 
     & B7913 =  129781d0/259200d0, 
     & B7914 = -238937d0/259200d0, 
     & B7915 =   27289d0/25920d0, 
     & B7916 = -197687d0/259200d0, 
     & B7917 =   88531d0/259200d0, 
     & B7918 = -156437d0/1814400d0, 
     & B7919 =   33953d0/3628800d0) 
 
       PARAMETER(B7921 = -17827d0/3628800d0, 
     & B7922 =  66043d0/1814400d0, 
     & B7923 = -30389d0/259200d0, 
     & B7924 =  55513d0/259200d0, 
     & B7925 =  -6281d0/25920d0, 
     & B7926 =  44983d0/259200d0, 
     & B7927 = -19859d0/259200d0, 
     & B7928 =  34453d0/1814400d0, 
     & B7929 =  -7297d0/3628800d0) 
 
       PARAMETER(B7931 =   8963d0/3628800d0, 
     & B7932 = -32987d0/1814400d0, 
     & B7933 =  15061d0/259200d0, 
     & B7934 = -27257d0/259200d0, 
     & B7935 =   3049d0/25920d0, 
     & B7936 = -21527d0/259200d0, 
     & B7937 =   9331d0/259200d0, 
     & B7938 = -15797d0/1814400d0, 
     & B7939 =   3233d0/3628800d0) 
 
       PARAMETER(B7941 =   3233d0/3628800d0, 
     & B7942 = -10067d0/1814400d0, 
     & B7943 =   3601d0/259200d0, 
     & B7944 =  -4337d0/259200d0, 
     & B7945 =     23d0/3240d0, 
     & B7946 =   1393d0/259200d0, 
     & B7947 =  -2129d0/259200d0, 
     & B7948 =   7123d0/1814400d0, 
     & B7949 =  -2497d0/3628800d0) 
 
 
       PARAMETER(B7951 =  -2497d0/3628800d0, 
     & B7952 =  12853d0/1814400d0, 
     & B7953 =  -7859d0/259200d0, 
     & B7954 =  18583d0/259200d0, 
     & B7955 =  -2681d0/25920d0, 
     & B7956 =  24313d0/259200d0, 
     & B7957 = -13589d0/259200d0, 
     & B7958 =  30043d0/1814400d0, 
     & B7959 =  -8227d0/3628800d0) 
 
       PARAMETER(L721 = 3.023839891568610D-01, 
     & L731 = 3.201698610574002D-05, 
     & L741 = 4.193101163680004D-04, 
     & L751 = 1.686924996069667D-04, 
     & L761 = 4.806043527549464D-05, 
     & L771 = 3.598347048026785D-06, 
     & L781 = 7.892534649789167D-04, 
     & L732 = 2.559868364091398D-01, 
     & L742 = 1.336896192287030D-04, 
     & L752 = 3.080994719931695D-03, 
     & L762 = 1.457177183563680D-04, 
     & L772 = 9.259360509484074D-04, 
     & L782 = 2.397658879381223D-04, 
     & L743 = 2.639734712170458D-01, 
     & L753 = 1.734338929611258D-04, 
     & L763 = 6.704398263264620D-03, 
     & L773 = 4.559927214651730D-05, 
     & L783 = 6.396418554053151D-05, 
     & L754 = 2.817729090368562D-01, 
     & L764 = 2.877761776030408D-04, 
     & L774 = 1.810919475521773D-04, 
     & L784 = 1.009049833235848D-03, 
     & L765 = 2.993040718034231D-01, 
     & L775 = 2.009850887505898D-03, 
     & L785 = 1.748065618845750D-03, 
     & L776 = 3.150349043479135D-01, 
     & L786 = 3.243816792609449D-05, 
     & L787 = 3.271307059448932D-01) 
 
       PARAMETER(CP71 = 103D0/9061D0, 
     & CP72 = -13D0/4480D0, 
     & CP73 =  67D0/42431D0) 
 
C--------- A8, B8, B86, B810 Cp8 := MATRICES DEFINING GAM9 
 
       PARAMETER(B911 = 1070017D0/3628800D0, 
     & B912 = 2233547D0/1814400D0, 
     & B913 = -2302297D0/1814400D0, 
     & B914 = 2797679D0/1814400D0, 
     & B915 = -31457D0/22680D0, 
     & B916 = 1573169D0/1814400D0, 
     & B917 = -645607D0/1814400D0, 
     & B918 = 156437D0/1814400D0, 
     & B919 = -33953D0/3628800D0, 
     & B921 = -33953D0/3628800D0, 
     & B922 = 687797D0/1814400D0, 
     & B923 =  1622393D0/1814400D0, 
     & B924 = -876271D0/1814400D0, 
     & B925 =   8233D0/22680D0, 
     & B926 =    -377521D0/1814400D0, 
     & B927 =   147143D0/1814400D0, 
     & B928 =  -34453D0/1814400D0, 
     & B929 =   7297D0/3628800D0, 
     & B931 = 7297D0/3628800D0, 
     & B932 =  -49813D0/1814400D0, 
     & B933 =  819143D0/1814400D0, 
     & B934 =  1315919D0/1814400D0, 
     & B935 = -5207D0/22680D0, 
     & B936 =  198929D0/1814400D0, 
     & B937 =  -71047D0/1814400D0, 
     & B938 =  15797D0/1814400D0, 
     & B939 = -3233D0/3628800D0, 
     & B941 = -3233D0/3628800D0, 
     & B942 = 18197D0/1814400D0, 
     & B943 =  -108007D0/1814400D0, 
     & B944 =  954929D0/1814400D0, 
     & B945 = 13903D0/22680D0, 
     & B946 = -212881D0/1814400D0, 
     & B947 =  63143D0/1814400D0, 
     & B948 =  -12853D0/1814400D0, 
     & B949 =  2497D0/3628800D0) 
 
 
c   B910=B9-B10; 
c               prime 4 righe : la 5 e' uguale alla 4; 
c     le ultime 4 uguali alle prime 4 negate (simmetriche); 
c     le prime  5 colonne  : le altre sono simmetriche e negate; 
c 
c[ 8183/1036800, -8183/115200,   8183/28800, -57281/86400,  57281/57600, -57281/57600,  57281/86400,  -8183/28800,  8183/115200, -8183/1036800] 
c[  -425/290304,    425/32256,    -425/8064,     425/3456,    -425/2304,     425/2304,    -425/3456,     425/8064,   -425/32256,    425/290304] 
c[      7/12800,    -63/12800,      63/3200,    -147/3200,     441/6400,    -441/6400,     147/3200,     -63/3200,     63/12800,      -7/12800] 
c[-2497/7257600,  2497/806400, -2497/201600,   2497/86400,  -2497/57600,   2497/57600,  -2497/86400,  2497/201600, -2497/806400,  2497/7257600] 
c[-2497/7257600,  2497/806400, -2497/201600,   2497/86400,  -2497/57600,   2497/57600,  -2497/86400,  2497/201600, -2497/806400,  2497/7257600] 
c[ 2497/7257600, -2497/806400,  2497/201600,  -2497/86400,   2497/57600,  -2497/57600,   2497/86400, -2497/201600,  2497/806400, -2497/7257600] 
c[     -7/12800,     63/12800,     -63/3200,     147/3200,    -441/6400,     441/6400,    -147/3200,      63/3200,    -63/12800,       7/12800] 
c[   425/290304,   -425/32256,     425/8064,    -425/3456,     425/2304,    -425/2304,     425/3456,    -425/8064,    425/32256,   -425/290304] 
c[-8183/1036800,  8183/115200,  -8183/28800,  57281/86400, -57281/57600,  57281/57600, -57281/86400,   8183/28800, -8183/115200,  8183/1036800] 
 
       PARAMETER(B91011 =   8183d0/1036800d0, 
     & B91012 =  -8183d0/115200d0, 
     & B91013 =   8183d0/28800d0, 
     & B91014 = -57281d0/86400d0, 
     & B91015 =  57281d0/57600d0) 
 
       PARAMETER(B91021 =   -425d0/290304d0, 
     & B91022 =    425d0/32256d0, 
     & B91023 =   -425d0/8064d0, 
     & B91024 =    425d0/3456d0, 
     & B91025 =   -425d0/2304d0) 
 
       PARAMETER(B91031 =     7d0/12800d0, 
     & B91032 =   -63d0/12800d0, 
     & B91033 =    63d0/3200d0, 
     & B91034 =  -147d0/3200d0, 
     & B91035 =   441d0/6400d0) 
 
       PARAMETER(B91041 = -2497d0/7257600d0, 
     & B91042 =  2497d0/806400d0, 
     & B91043 = -2497d0/201600d0, 
     & B91044 =  2497d0/86400d0, 
     & B91045 = -2497d0/57600d0) 
 
       PARAMETER(Cp91 =  7.892554012345216d-03, 
     & Cp92 = -1.463982583774219d-03, 
     & Cp93 =  5.468749999999983d-04, 
     & Cp94 = -3.440531305114634d-04) 
 
C--------- THE OTHERS ARE THE SAME WITH CHANGED SIGN 
 
       PARAMETER( L921 = 2.590721934790442d-01, 
     & L932   =   2.077575545359853d-01, 
     & L943   =   2.032874698558627d-01, 
     & L954   =   2.036384888660128d-01, 
     & L965   =   2.039599505779785d-01, 
     & L976   =   2.034044409161703d-01, 
     & L987   =   2.017245408702437d-01, 
     & L998   =   1.986549276295617d-01) 
 
C
C   INPUT VARIABLES
C------------------------------------
      INTEGER R, IT, IJOB, IPIV(R), LDLU, IPAR(*)
      DOUBLE PRECISION  H, SCAL(R), TP(*), ERRNEWT0,    
     &                   LU(LDLU,R),RPAR(*)
C
C   INPUT/OUTPUT VARIABLES
C------------------------------------
       INTEGER NFCN
       DOUBLE PRECISION  ERRNEWT, TETAK0, YP(R,*), FP(R,*), F1(R,*),
     &                   DN(R) 
       LOGICAL TER
C
C   LOCAL VARIABLES
C------------------------------------
      INTEGER  J 
      DOUBLE PRECISION  ERRVJ, SUM
C
C   EXTERNAL FUNCTIONS
C------------------------------------

      EXTERNAL FCN

C
C   EXECUTABLE STATEMENTS
C---------------------------------
C          
C--------- ONE STEP OF THE ITERATION PROCEDURE
      TER = .FALSE.
C
C  
      DO J=1,R
         SUM = B311*FP(J,1)+B312*FP(J,2)+B313*FP(J,3)
         DN(J)=YP(J,2)-YP(J,1)-H*SUM
      END DO
      CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
      ERRVJ = 0D0
      DO J=1,R
         YP(J,2)=YP(J,2)-DN(J)
         SUM = (DN(J)/SCAL(J))
         ERRVJ =  ERRVJ + SUM*SUM
      END DO
      ERRVJ = dsqrt(ERRVJ/R)
      ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
      IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
         TER = .TRUE.
         RETURN
      END IF
      CALL FCN(R,TP(2),YP(1,2),F1(1,1), RPAR,IPAR)
      NFCN = NFCN + 1


      DO J=1,R
         SUM = L321*(F1(J,1)-FP(J,2))+B311*FP(J,2)
     &                  +B312*FP(J,3)+B313*FP(J,4)
         DN(J)=YP(J,3)-YP(J,2)-H*SUM
      END DO
      CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
      ERRVJ = 0D0
      DO J=1,R
        YP(J,3)=YP(J,3)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
      END DO
      ERRVJ = DSQRT(ERRVJ/R)
      ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
      IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
         TER = .TRUE.
         RETURN
      END IF
      CALL FCN(R,TP(3),YP(1,3),F1(1,2), RPAR,IPAR)
      NFCN = NFCN + 1


      DO J=1,R
        SUM = L331*(F1(J,1)-FP(J,2))+L332*(F1(J,2)-FP(J,3))
     &              +B311*FP(J,3)+B312*FP(J,4)+B313*FP(J,5)
        DN(J)=YP(J,4)-YP(J,3)-H*SUM
      END DO
      CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
      ERRVJ = 0D0
      DO J=1,R
         YP(J,4)=YP(J,4)-DN(J)
         SUM = (DN(J)/SCAL(J))
         ERRVJ =  ERRVJ + SUM*SUM
      END DO
      ERRVJ = DSQRT(ERRVJ/R)
      ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
      IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
         TER = .TRUE.
         RETURN
       END IF
       CALL FCN(R,TP(4),YP(1,4),F1(1,3), RPAR,IPAR)
       NFCN = NFCN + 1

       DO J=1,R
         SUM = L341*(F1(J,1)-FP(J,2))+L342*(F1(J,2)-FP(J,3))
     &                               +L343*(F1(J,3)-FP(J,4))
     &               +B313*FP(J,3)+B312*FP(J,4)+B311*FP(J,5)
         DN(J)=YP(J,5)-YP(J,4)-H*SUM
       END DO
       CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
       ERRVJ = 0D0
       DO J=1,R
          YP(J,5)=YP(J,5)-DN(J)
          SUM = (DN(J)/SCAL(J))
          ERRVJ =  ERRVJ + SUM*SUM
       END DO
       ERRVJ = DSQRT(ERRVJ/R)
       ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
       IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
          TER = .TRUE.
          RETURN
       END IF
       CALL FCN(R,TP(5),YP(1,5),F1(1,4), RPAR,IPAR)
       NFCN = NFCN + 1

       DO J=1,R
         FP(J,2) = F1(J,1)
         FP(J,3) = F1(J,2)
         FP(J,4) = F1(J,3)
         FP(J,5) = F1(J,4)
       END DO


       RETURN
       END
C
C  SUBROUTINE TERMNOT5  (ORDER 5)
C
      SUBROUTINE  TERMNOT5(R,FCN,H,IT,DN, F1,FP,YP,TP,NFCN,
     &  ERRNEWT,ERRNEWT0,TETAK0,LU, LDLU,IPIV, SCAL,IJOB,TER,
     &  RPAR,IPAR)
       IMPLICIT NONE
C
C   INCLUDE
C------------------------------------
       DOUBLE PRECISION 
     & B3511, B3512, B3513, B3514, B3515, B3521, B3522, B3523, B3524, 
     & B3525, B3531, B3532, B3533, B3534, B3535, B3541, B3542, B3543, 
     & B3544, B3545, 
     & L321, L331, L332, L341, L342, L343, B311, B312, B313, 
     & B511, B512, B513, B514, B515, B521, B522, B523, B524, B525, 
     & L521, L531, L541, L551, L561, L532, L542, L552, L562, L543, 
     & L553, L563, L554, L564, L565, 
     & B5711, B5712, B5713, B5714, B5715, B5716, B5717, B5721, B5722, 
     & B5723, B5724, B5725, B5726, B5727, B5731, B5732, B5733, B5734, 
     & B5735, B5736, B5737, B5741, B5742, B5743, B5744, B5745, B5746, 
     & B5747, 
     & B711, B712, B713, B714, B715, B716, B717, B721, B722, B723, 
     & B724, B725, B726, B727, B731, B732, B733, B734, B735, B736, B737, 
     & B7911, B7912, B7913, B7914, B7915, B7916, B7917, B7918, B7919, 
     & B7921, B7922, B7923, B7924, B7925, B7926, B7927, B7928, B7929, 
     & B7931, B7932, B7933, B7934, B7935, B7936, B7937, B7938, B7939, 
     & B7941, B7942, B7943, B7944, B7945, B7946, B7947, B7948, B7949, 
     & B7951, B7952, B7953, B7954, B7955, B7956, B7957, B7958, B7959, 
     & L721, L731, L741, L751, L761, L771, L781, L732, L742, L752, 
     & L762, L772, L782, L743, L753, L763, L773, L783, L754, L764, 
     & L774, L784, L765, L775, L785, L776, L786, L787, 
     & CP31, CP51, CP52, CP71, CP72, CP73, 
     & CP91, CP92, CP93, CP94, 
     & B911, B912, B913, B914, B915, B916, B917, B918, B919, B921, 
     & B922, B923, B924, B925, B926, B927, B928, B929, B931, B932, 
     & B933, B934, B935, B936, B937, B938, B939, B941, B942, B943, 
     & B944, B945, B946, B947, B948, B949, 
     & L921, L932, L943, L954, L965, L976, L987, L998, 
     & B91011, B91012, B91013, B91014, B91015, B91021, B91022, B91023, 
     & B91024, B91025, B91031, B91032, B91033, B91034, B91035, B91041, 
     & B91042, B91043, B91044, B91045 
       PARAMETER( B311 = 5d0/12d0, 
     &            B312 = 8d0/12d0, 
     &            B313 = -1d0/12d0) 
 
       PARAMETER( Cp31 = 1d0/24d0) 
ccccc L3 dimensione 4 
       PARAMETER( L321 =  4.012395124208693d-01, 
     &            L331 =  8.819910099032529d-04, 
     &            L341 =  1.728116022258560d-04, 
     &            L332 =  3.680857287181573d-01, 
     &            L342 =  1.635381132422046d-03, 
     &            L343 =  3.688541178419062d-01) 
c--------- B3-B5 OF DIMENSION  4 
       PARAMETER( B3511 = 49d0/720d0, 
     &            B3512 = -83d0/360d0, 
     &            B3513 = 17d0/60d0, 
     &            B3514 = -53d0/360d0, 
     &            B3515 = 19d0/720d0, 
     &            B3521 = 19d0/720d0, 
     &            B3522 = -23d0/360d0, 
     &            B3523 = 1d0/30d0, 
     & B3524 = 7d0/360d0, 
     & B3525 = -11d0/720d0, 
     & B3531 = -11d0/720d0, 
     & B3532 = 37d0/360d0, 
     & B3533 = -13d0/60d0, 
     & B3534 = 67d0/360d0, 
     & B3535 = -41d0/720d0, 
     & B3541 = 19d0/720d0, 
     & B3542 = -53d0/360d0, 
     & B3543 =  17d0/60d0, 
     & B3544 = -83d0/360d0, 
     & B3545 =  49d0/720d0) 
 
c[ 49/720, -83/360,  17/60, -53/360,  19/720] 
c[ 19/720, -23/360,   1/30,   7/360, -11/720] 
c[-11/720,  37/360, -13/60,  67/360, -41/720] 
c[ 19/720, -53/360,  17/60, -83/360,  49/720] 
 
C--------- A5, B5, B53, B56 Cp5 := MATRICES DEFINING  GAM5 
 
       PARAMETER( B511 = 251d0/720d0, 
     & B512 = 323d0/360d0, 
     & B513 = - 11d0/30d0, 
     & B514 = 53d0/360d0, 
     & B515 = -19d0/720d0, 
     & B521 = -19d0/720d0, 
     & B522 =  173d0/360d0, 
     & B523 = 19d0/30d0, 
     & B524 = -37d0/360d0, 
     & B525 = 11d0/720d0) 
 
       PARAMETER( B5711 = 1997d0/60480d0, 
     & B5712 = -113d0/630d0, 
     & B5713 = 1619d0/4032d0, 
     & B5714 = -715d0/1512d0, 
     & B5715 = 1241d0/4032d0, 
     & B5716 = -263d0/2520d0, 
     & B5717 = 863d0/60480d0, 
     & B5721 = -733d0/60480d0, 
     & B5722 = 41d0/630d0, 
     & B5723 = -193d0/1344d0, 
     & B5724 = 251d0/1512d0, 
     & B5725 = -425d0/4032d0, 
     & B5726 = 29d0/840d0, 
     & B5727 = -271d0/60480d0, 
     & B5731 = -271d0/60480d0, 
     & B5732 = 97d0/5040d0, 
     & B5733 = -13d0/448d0, 
     & B5734 = 5d0/378d0, 
     & B5735 = 37d0/4032d0, 
     & B5736 = -19d0/1680d0, 
     & B5737 = 191d0/60480d0, 
     & B5741 = 191d0/60480d0, 
     & B5742 = -67d0/2520d0, 
     & B5743 = 115d0/1344d0, 
     & B5744 = -211d0/1512d0, 
     & B5745 = 499d0/4032d0, 
     & B5746 = -2d0/35d0, 
     & B5747 = 653d0/60480d0) 
 
c--------- B5 - B7 
c 
c[1997/60480,  -113/630, 1619/4032, -715/1512, 1241/4032, -263/2520,  863/60480] 
c[-733/60480,    41/630, -193/1344,  251/1512, -425/4032,    29/840, -271/60480] 
c[-271/60480,   97/5040,   -13/448,     5/378,   37/4032,  -19/1680,  191/60480] 
c[ 191/60480,  -67/2520,  115/1344, -211/1512,  499/4032,     -2/35,  653/60480] 
c[-271/60480,    29/840, -425/4032,  251/1512, -193/1344,    41/630, -733/60480] 
c[ 863/60480, -263/2520, 1241/4032, -715/1512, 1619/4032,  -113/630, 1997/60480] 
c 
c 
       PARAMETER( L521 =  3.668340831928216D-01, 
     & L531 =  2.477905683677308D-03, 
     & L541 = -1.919925047010838D-03, 
     & L551 =  2.218385581234200D-03, 
     & L561 = -5.442189351609260D-03, 
     & L532 =  3.216639533696728D-01, 
     & L542 =  1.231925763308414D-03, 
     & L552 =  7.841944627374794D-03, 
     & L562 =  1.002485104590053D-03, 
     & L543 =  3.375100828961925D-01, 
     & L553 = -2.614300734741796D-04, 
     & L563 =  1.066631182323580D-03, 
     & L554 =  3.523137378783708D-01, 
     & L564 = -3.596681121610224D-04, 
     & L565 =  3.617716171655064D-01, 
     & CP51 =   3D0/160D0, 
     & CP52 = -11D0/1440D0) 
 
C--------- A7, B7, B57, B58 Cp7 := MATRICES DEFINING GAM7 
 
        PARAMETER( B711 = 19087d0/60480d0, 
     & B712 = 2713d0/2520d0, 
     & B713 = -15487d0/20160d0, 
     & B714 = 586d0/945d0, 
     & B715 = -6737d0/20160d0, 
     & B716 = 263d0/2520d0, 
     & B717 = -863d0/60480d0, 
     & B721 = -863d0/60480d0, 
     & B722 = 349d0/840d0, 
     & B723 = 5221d0/6720d0, 
     & B724 = -254d0/945d0, 
     & B725 = 811d0/6720d0, 
     & B726 = -29d0/840d0, 
     & B727 = 271d0/60480d0, 
     & B731 = 271d0/60480d0, 
     & B732 = -23d0/504d0, 
     & B733 = 10273d0/20160d0, 
     & B734 = 586d0/945d0, 
     & B735 = -2257d0/20160d0, 
     & B736 = 67d0/2520d0, 
     & B737 = -191d0/60480d0) 
 
C 
C--------- THE LAST THREE ROWS ARE THE REVERSE OF THE FIRST THREE 
C 
c B79 = 
c 
c[ 75203/3628800, -280187/1814400, 129781/259200, -238937/259200, 27289/25920, -197687/259200,  88531/259200, -156437/1814400,  33953/3628800] 
c[-17827/3628800,   66043/1814400, -30389/259200,   55513/259200, -6281/25920,   44983/259200, -19859/259200,   34453/1814400,  -7297/3628800] 
c[  8963/3628800,  -32987/1814400,  15061/259200,  -27257/259200,  3049/25920,  -21527/259200,   9331/259200,  -15797/1814400,   3233/3628800] 
c[  3233/3628800,  -10067/1814400,   3601/259200,   -4337/259200,     23/3240,    1393/259200,  -2129/259200,    7123/1814400,  -2497/3628800] 
c[ -2497/3628800,   12853/1814400,  -7859/259200,   18583/259200, -2681/25920,   24313/259200, -13589/259200,   30043/1814400,  -8227/3628800] 
c[  3233/3628800,  -15797/1814400,   9331/259200,  -21527/259200,  3049/25920,  -27257/259200,  15061/259200,  -32987/1814400,   8963/3628800] 
c[ -7297/3628800,   34453/1814400, -19859/259200,   44983/259200, -6281/25920,   55513/259200, -30389/259200,   66043/1814400, -17827/3628800] 
c[ 33953/3628800, -156437/1814400,  88531/259200, -197687/259200, 27289/25920, -238937/259200, 129781/259200, -280187/1814400,  75203/3628800] 
c 
       PARAMETER( B7911 =   75203d0/3628800d0, 
     & B7912 = -280187d0/1814400d0, 
     & B7913 =  129781d0/259200d0, 
     & B7914 = -238937d0/259200d0, 
     & B7915 =   27289d0/25920d0, 
     & B7916 = -197687d0/259200d0, 
     & B7917 =   88531d0/259200d0, 
     & B7918 = -156437d0/1814400d0, 
     & B7919 =   33953d0/3628800d0) 
 
       PARAMETER(B7921 = -17827d0/3628800d0, 
     & B7922 =  66043d0/1814400d0, 
     & B7923 = -30389d0/259200d0, 
     & B7924 =  55513d0/259200d0, 
     & B7925 =  -6281d0/25920d0, 
     & B7926 =  44983d0/259200d0, 
     & B7927 = -19859d0/259200d0, 
     & B7928 =  34453d0/1814400d0, 
     & B7929 =  -7297d0/3628800d0) 
 
       PARAMETER(B7931 =   8963d0/3628800d0, 
     & B7932 = -32987d0/1814400d0, 
     & B7933 =  15061d0/259200d0, 
     & B7934 = -27257d0/259200d0, 
     & B7935 =   3049d0/25920d0, 
     & B7936 = -21527d0/259200d0, 
     & B7937 =   9331d0/259200d0, 
     & B7938 = -15797d0/1814400d0, 
     & B7939 =   3233d0/3628800d0) 
 
       PARAMETER(B7941 =   3233d0/3628800d0, 
     & B7942 = -10067d0/1814400d0, 
     & B7943 =   3601d0/259200d0, 
     & B7944 =  -4337d0/259200d0, 
     & B7945 =     23d0/3240d0, 
     & B7946 =   1393d0/259200d0, 
     & B7947 =  -2129d0/259200d0, 
     & B7948 =   7123d0/1814400d0, 
     & B7949 =  -2497d0/3628800d0) 
 
 
       PARAMETER(B7951 =  -2497d0/3628800d0, 
     & B7952 =  12853d0/1814400d0, 
     & B7953 =  -7859d0/259200d0, 
     & B7954 =  18583d0/259200d0, 
     & B7955 =  -2681d0/25920d0, 
     & B7956 =  24313d0/259200d0, 
     & B7957 = -13589d0/259200d0, 
     & B7958 =  30043d0/1814400d0, 
     & B7959 =  -8227d0/3628800d0) 
 
       PARAMETER(L721 = 3.023839891568610D-01, 
     & L731 = 3.201698610574002D-05, 
     & L741 = 4.193101163680004D-04, 
     & L751 = 1.686924996069667D-04, 
     & L761 = 4.806043527549464D-05, 
     & L771 = 3.598347048026785D-06, 
     & L781 = 7.892534649789167D-04, 
     & L732 = 2.559868364091398D-01, 
     & L742 = 1.336896192287030D-04, 
     & L752 = 3.080994719931695D-03, 
     & L762 = 1.457177183563680D-04, 
     & L772 = 9.259360509484074D-04, 
     & L782 = 2.397658879381223D-04, 
     & L743 = 2.639734712170458D-01, 
     & L753 = 1.734338929611258D-04, 
     & L763 = 6.704398263264620D-03, 
     & L773 = 4.559927214651730D-05, 
     & L783 = 6.396418554053151D-05, 
     & L754 = 2.817729090368562D-01, 
     & L764 = 2.877761776030408D-04, 
     & L774 = 1.810919475521773D-04, 
     & L784 = 1.009049833235848D-03, 
     & L765 = 2.993040718034231D-01, 
     & L775 = 2.009850887505898D-03, 
     & L785 = 1.748065618845750D-03, 
     & L776 = 3.150349043479135D-01, 
     & L786 = 3.243816792609449D-05, 
     & L787 = 3.271307059448932D-01) 
 
       PARAMETER(CP71 = 103D0/9061D0, 
     & CP72 = -13D0/4480D0, 
     & CP73 =  67D0/42431D0) 
 
C--------- A8, B8, B86, B810 Cp8 := MATRICES DEFINING GAM9 
 
       PARAMETER(B911 = 1070017D0/3628800D0, 
     & B912 = 2233547D0/1814400D0, 
     & B913 = -2302297D0/1814400D0, 
     & B914 = 2797679D0/1814400D0, 
     & B915 = -31457D0/22680D0, 
     & B916 = 1573169D0/1814400D0, 
     & B917 = -645607D0/1814400D0, 
     & B918 = 156437D0/1814400D0, 
     & B919 = -33953D0/3628800D0, 
     & B921 = -33953D0/3628800D0, 
     & B922 = 687797D0/1814400D0, 
     & B923 =  1622393D0/1814400D0, 
     & B924 = -876271D0/1814400D0, 
     & B925 =   8233D0/22680D0, 
     & B926 =    -377521D0/1814400D0, 
     & B927 =   147143D0/1814400D0, 
     & B928 =  -34453D0/1814400D0, 
     & B929 =   7297D0/3628800D0, 
     & B931 = 7297D0/3628800D0, 
     & B932 =  -49813D0/1814400D0, 
     & B933 =  819143D0/1814400D0, 
     & B934 =  1315919D0/1814400D0, 
     & B935 = -5207D0/22680D0, 
     & B936 =  198929D0/1814400D0, 
     & B937 =  -71047D0/1814400D0, 
     & B938 =  15797D0/1814400D0, 
     & B939 = -3233D0/3628800D0, 
     & B941 = -3233D0/3628800D0, 
     & B942 = 18197D0/1814400D0, 
     & B943 =  -108007D0/1814400D0, 
     & B944 =  954929D0/1814400D0, 
     & B945 = 13903D0/22680D0, 
     & B946 = -212881D0/1814400D0, 
     & B947 =  63143D0/1814400D0, 
     & B948 =  -12853D0/1814400D0, 
     & B949 =  2497D0/3628800D0) 
 
 
c   B910=B9-B10; 
c               prime 4 righe : la 5 e' uguale alla 4; 
c     le ultime 4 uguali alle prime 4 negate (simmetriche); 
c     le prime  5 colonne  : le altre sono simmetriche e negate; 
c 
c[ 8183/1036800, -8183/115200,   8183/28800, -57281/86400,  57281/57600, -57281/57600,  57281/86400,  -8183/28800,  8183/115200, -8183/1036800] 
c[  -425/290304,    425/32256,    -425/8064,     425/3456,    -425/2304,     425/2304,    -425/3456,     425/8064,   -425/32256,    425/290304] 
c[      7/12800,    -63/12800,      63/3200,    -147/3200,     441/6400,    -441/6400,     147/3200,     -63/3200,     63/12800,      -7/12800] 
c[-2497/7257600,  2497/806400, -2497/201600,   2497/86400,  -2497/57600,   2497/57600,  -2497/86400,  2497/201600, -2497/806400,  2497/7257600] 
c[-2497/7257600,  2497/806400, -2497/201600,   2497/86400,  -2497/57600,   2497/57600,  -2497/86400,  2497/201600, -2497/806400,  2497/7257600] 
c[ 2497/7257600, -2497/806400,  2497/201600,  -2497/86400,   2497/57600,  -2497/57600,   2497/86400, -2497/201600,  2497/806400, -2497/7257600] 
c[     -7/12800,     63/12800,     -63/3200,     147/3200,    -441/6400,     441/6400,    -147/3200,      63/3200,    -63/12800,       7/12800] 
c[   425/290304,   -425/32256,     425/8064,    -425/3456,     425/2304,    -425/2304,     425/3456,    -425/8064,    425/32256,   -425/290304] 
c[-8183/1036800,  8183/115200,  -8183/28800,  57281/86400, -57281/57600,  57281/57600, -57281/86400,   8183/28800, -8183/115200,  8183/1036800] 
 
       PARAMETER(B91011 =   8183d0/1036800d0, 
     & B91012 =  -8183d0/115200d0, 
     & B91013 =   8183d0/28800d0, 
     & B91014 = -57281d0/86400d0, 
     & B91015 =  57281d0/57600d0) 
 
       PARAMETER(B91021 =   -425d0/290304d0, 
     & B91022 =    425d0/32256d0, 
     & B91023 =   -425d0/8064d0, 
     & B91024 =    425d0/3456d0, 
     & B91025 =   -425d0/2304d0) 
 
       PARAMETER(B91031 =     7d0/12800d0, 
     & B91032 =   -63d0/12800d0, 
     & B91033 =    63d0/3200d0, 
     & B91034 =  -147d0/3200d0, 
     & B91035 =   441d0/6400d0) 
 
       PARAMETER(B91041 = -2497d0/7257600d0, 
     & B91042 =  2497d0/806400d0, 
     & B91043 = -2497d0/201600d0, 
     & B91044 =  2497d0/86400d0, 
     & B91045 = -2497d0/57600d0) 
 
       PARAMETER(Cp91 =  7.892554012345216d-03, 
     & Cp92 = -1.463982583774219d-03, 
     & Cp93 =  5.468749999999983d-04, 
     & Cp94 = -3.440531305114634d-04) 
 
C--------- THE OTHERS ARE THE SAME WITH CHANGED SIGN 
 
       PARAMETER( L921 = 2.590721934790442d-01, 
     & L932   =   2.077575545359853d-01, 
     & L943   =   2.032874698558627d-01, 
     & L954   =   2.036384888660128d-01, 
     & L965   =   2.039599505779785d-01, 
     & L976   =   2.034044409161703d-01, 
     & L987   =   2.017245408702437d-01, 
     & L998   =   1.986549276295617d-01) 
 
C
C   INPUT VARIABLES
C------------------------------------
      INTEGER R, IT, IJOB, IPIV(R), LDLU, IPAR(*)
      DOUBLE PRECISION  H, SCAL(R), TP(*), ERRNEWT0,    
     &                   LU(LDLU,R),RPAR(*)
C
C   INPUT/OUTPUT VARIABLES
C------------------------------------
       INTEGER NFCN
       DOUBLE PRECISION  ERRNEWT, TETAK0, YP(R,*), FP(R,*), F1(R,*),
     &                   DN(R) 
       LOGICAL TER
C
C   LOCAL VARIABLES
C------------------------------------
      INTEGER  J 
      DOUBLE PRECISION  ERRVJ, SUM
C
C   EXTERNAL FUNCTIONS
C------------------------------------

      EXTERNAL FCN

C
C   EXECUTABLE STATEMENTS
C---------------------------------
C          
C--------- ONE STEP OF THE ITERATION PROCEDURE

      TER = .FALSE.
      DO J=1,R
         SUM = B511*FP(J,1)+B512*FP(J,2)+B513*FP(J,3)+B514*FP(J,4)
     &                                               +B515*FP(J,5)
         DN(J)=YP(J,2)-YP(J,1)-H*SUM
      END DO
              CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
              ERRVJ = 0D0
              DO J=1,R
                YP(J,2)=YP(J,2)-DN(J)
                SUM = (DN(J)/SCAL(J))
                ERRVJ =  ERRVJ + SUM*SUM
              END DO
              ERRVJ = DSQRT(ERRVJ/R)
              ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
              IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
                 TER = .TRUE.
                 RETURN
              END IF
              CALL FCN(R,TP(2),YP(1,2),F1(1,1), RPAR,IPAR)
              NFCN = NFCN + 1

      DO J=1,R
         SUM = L521*(F1(J,1)-FP(J,2))+B521*FP(J,1)+B522*FP(J,2)
     &                  +B523*FP(J,3)+B524*FP(J,4)+B525*FP(J,5)
         DN(J)=YP(J,3)-YP(J,2)-H*SUM
      END DO
              CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
              ERRVJ = 0D0
              DO J=1,R
                YP(J,3)=YP(J,3)-DN(J)
                SUM = (DN(J)/SCAL(J))
                ERRVJ =  ERRVJ + SUM*SUM
              END DO
              ERRVJ = DSQRT(ERRVJ/R)
              ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
              IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
                 TER = .TRUE.
                 RETURN
              END IF

              CALL FCN(R,TP(3),YP(1,3),F1(1,2), RPAR,IPAR)
              NFCN = NFCN + 1

      DO J=1,R
          SUM = L531*(F1(J,1)-FP(J,2))+L532*(F1(J,2)-FP(J,3))
     &                +B521*FP(J,2)+B522*FP(J,3)+B523*FP(J,4)
     &                             +B524*FP(J,5)+B525*FP(J,6)
          DN(J)=YP(J,4)-YP(J,3)-H*SUM
      END DO
              CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
              ERRVJ = 0D0
              DO J=1,R
                YP(J,4)=YP(J,4)-DN(J)
                SUM = (DN(J)/SCAL(J))
                ERRVJ =  ERRVJ + SUM*SUM
              END DO
              ERRVJ = DSQRT(ERRVJ/R)
              ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
              IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
                 TER = .TRUE.
                 RETURN
              END IF
              CALL FCN(R,TP(4),YP(1,4),F1(1,3), RPAR,IPAR)
              NFCN = NFCN + 1

       DO J=1,R
          SUM = L541*(F1(J,1)-FP(J,2))+L542*(F1(J,2)-FP(J,3))
     &      +L543*(F1(J,3)-FP(J,4))+B521*FP(J,3)+B522*FP(J,4)
     &                +B523*FP(J,5)+B524*FP(J,6)+B525*FP(J,7)
          DN(J)=YP(J,5)-YP(J,4)-H*SUM
       END DO
              CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
              ERRVJ = 0D0
              DO J=1,R
                YP(J,5)=YP(J,5)-DN(J)
               SUM = (DN(J)/SCAL(J))
               ERRVJ =  ERRVJ + SUM*SUM
              END DO
              ERRVJ = DSQRT(ERRVJ/R)
              ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
              IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT.TETAK0)) THEN
                 TER = .TRUE.
                 RETURN
              END IF

              CALL FCN(R,TP(5),YP(1,5),F1(1,4), RPAR,IPAR)
              NFCN = NFCN + 1

      DO J=1,R
         SUM =  L551*(F1(J,1)-FP(J,2))+L552*(F1(J,2)-FP(J,3))
     &         +L553*(F1(J,3)-FP(J,4))+L554*(F1(J,4)-FP(J,5))
     &+B525*FP(J,3)+B524*FP(J,4)+B523*FP(J,5)+B522*FP(J,6)+B521*FP(J,7)
         DN(J)=YP(J,6)-YP(J,5)-H*SUM
      END DO
              CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
              ERRVJ = 0D0
              DO J=1,R
                YP(J,6)=YP(J,6)-DN(J)
                SUM = (DN(J)/SCAL(J))
                ERRVJ =  ERRVJ + SUM*SUM
              END DO
              ERRVJ = DSQRT(ERRVJ/R)
              ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
              IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0)) THEN
                 TER = .TRUE.
                 RETURN
              END IF
              CALL FCN(R,TP(6),YP(1,6),F1(1,5), RPAR,IPAR)
              NFCN = NFCN + 1


      DO J=1,R
         SUM =  L561*(F1(J,1)-FP(J,2))+L562*(F1(J,2)-FP(J,3))
     &         +L563*(F1(J,3)-FP(J,4))+L564*(F1(J,4)-FP(J,5))
     &                   +L565*(F1(J,5)-FP(J,6))+B515*FP(J,3)
     &   +B514*FP(J,4)+B513*FP(J,5)+B512*FP(J,6)+B511*FP(J,7)
         DN(J)=YP(J,7)-YP(J,6)-H*SUM
      END DO
              CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
              ERRVJ = 0D0
              DO J=1,R
                YP(J,7)=YP(J,7)-DN(J)
                SUM = (DN(J)/SCAL(J))
                ERRVJ =  ERRVJ + SUM*SUM
              END DO
              ERRVJ = dsqrt(ERRVJ/R)
              ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
              if ((it.ge.1).and.(errnewt/errnewt0 .gt.tetak0)) then
                 TER = .TRUE.
                 RETURN
              end if
              CALL FCN(R,TP(7),YP(1,7),F1(1,6), RPAR,IPAR)
              NFCN = NFCN + 1



               DO J=1,R
                 FP(J,2) = F1(J,1)
                 FP(J,3) = F1(J,2)
                 FP(J,4) = F1(J,3)
                 FP(J,5) = F1(J,4)
                 FP(J,6) = F1(J,5)
                 FP(J,7) = F1(J,6)
               END DO
          RETURN
         END

C
C  SUBROUTINE TERMNOT7  (ORDER 7)
C
      SUBROUTINE  TERMNOT7(R,FCN,H,IT,DN, F1,FP,YP,TP,NFCN,
     &  ERRNEWT,ERRNEWT0,TETAK0,LU, LDLU,IPIV, SCAL,IJOB,TER,
     &  RPAR,IPAR)
       IMPLICIT NONE
C
C   INCLUDE
C------------------------------------
       DOUBLE PRECISION 
     & B3511, B3512, B3513, B3514, B3515, B3521, B3522, B3523, B3524, 
     & B3525, B3531, B3532, B3533, B3534, B3535, B3541, B3542, B3543, 
     & B3544, B3545, 
     & L321, L331, L332, L341, L342, L343, B311, B312, B313, 
     & B511, B512, B513, B514, B515, B521, B522, B523, B524, B525, 
     & L521, L531, L541, L551, L561, L532, L542, L552, L562, L543, 
     & L553, L563, L554, L564, L565, 
     & B5711, B5712, B5713, B5714, B5715, B5716, B5717, B5721, B5722, 
     & B5723, B5724, B5725, B5726, B5727, B5731, B5732, B5733, B5734, 
     & B5735, B5736, B5737, B5741, B5742, B5743, B5744, B5745, B5746, 
     & B5747, 
     & B711, B712, B713, B714, B715, B716, B717, B721, B722, B723, 
     & B724, B725, B726, B727, B731, B732, B733, B734, B735, B736, B737, 
     & B7911, B7912, B7913, B7914, B7915, B7916, B7917, B7918, B7919, 
     & B7921, B7922, B7923, B7924, B7925, B7926, B7927, B7928, B7929, 
     & B7931, B7932, B7933, B7934, B7935, B7936, B7937, B7938, B7939, 
     & B7941, B7942, B7943, B7944, B7945, B7946, B7947, B7948, B7949, 
     & B7951, B7952, B7953, B7954, B7955, B7956, B7957, B7958, B7959, 
     & L721, L731, L741, L751, L761, L771, L781, L732, L742, L752, 
     & L762, L772, L782, L743, L753, L763, L773, L783, L754, L764, 
     & L774, L784, L765, L775, L785, L776, L786, L787, 
     & CP31, CP51, CP52, CP71, CP72, CP73, 
     & CP91, CP92, CP93, CP94, 
     & B911, B912, B913, B914, B915, B916, B917, B918, B919, B921, 
     & B922, B923, B924, B925, B926, B927, B928, B929, B931, B932, 
     & B933, B934, B935, B936, B937, B938, B939, B941, B942, B943, 
     & B944, B945, B946, B947, B948, B949, 
     & L921, L932, L943, L954, L965, L976, L987, L998, 
     & B91011, B91012, B91013, B91014, B91015, B91021, B91022, B91023, 
     & B91024, B91025, B91031, B91032, B91033, B91034, B91035, B91041, 
     & B91042, B91043, B91044, B91045 
       PARAMETER( B311 = 5d0/12d0, 
     &            B312 = 8d0/12d0, 
     &            B313 = -1d0/12d0) 
 
       PARAMETER( Cp31 = 1d0/24d0) 
ccccc L3 dimensione 4 
       PARAMETER( L321 =  4.012395124208693d-01, 
     &            L331 =  8.819910099032529d-04, 
     &            L341 =  1.728116022258560d-04, 
     &            L332 =  3.680857287181573d-01, 
     &            L342 =  1.635381132422046d-03, 
     &            L343 =  3.688541178419062d-01) 
c--------- B3-B5 OF DIMENSION  4 
       PARAMETER( B3511 = 49d0/720d0, 
     &            B3512 = -83d0/360d0, 
     &            B3513 = 17d0/60d0, 
     &            B3514 = -53d0/360d0, 
     &            B3515 = 19d0/720d0, 
     &            B3521 = 19d0/720d0, 
     &            B3522 = -23d0/360d0, 
     &            B3523 = 1d0/30d0, 
     & B3524 = 7d0/360d0, 
     & B3525 = -11d0/720d0, 
     & B3531 = -11d0/720d0, 
     & B3532 = 37d0/360d0, 
     & B3533 = -13d0/60d0, 
     & B3534 = 67d0/360d0, 
     & B3535 = -41d0/720d0, 
     & B3541 = 19d0/720d0, 
     & B3542 = -53d0/360d0, 
     & B3543 =  17d0/60d0, 
     & B3544 = -83d0/360d0, 
     & B3545 =  49d0/720d0) 
 
c[ 49/720, -83/360,  17/60, -53/360,  19/720] 
c[ 19/720, -23/360,   1/30,   7/360, -11/720] 
c[-11/720,  37/360, -13/60,  67/360, -41/720] 
c[ 19/720, -53/360,  17/60, -83/360,  49/720] 
 
C--------- A5, B5, B53, B56 Cp5 := MATRICES DEFINING  GAM5 
 
       PARAMETER( B511 = 251d0/720d0, 
     & B512 = 323d0/360d0, 
     & B513 = - 11d0/30d0, 
     & B514 = 53d0/360d0, 
     & B515 = -19d0/720d0, 
     & B521 = -19d0/720d0, 
     & B522 =  173d0/360d0, 
     & B523 = 19d0/30d0, 
     & B524 = -37d0/360d0, 
     & B525 = 11d0/720d0) 
 
       PARAMETER( B5711 = 1997d0/60480d0, 
     & B5712 = -113d0/630d0, 
     & B5713 = 1619d0/4032d0, 
     & B5714 = -715d0/1512d0, 
     & B5715 = 1241d0/4032d0, 
     & B5716 = -263d0/2520d0, 
     & B5717 = 863d0/60480d0, 
     & B5721 = -733d0/60480d0, 
     & B5722 = 41d0/630d0, 
     & B5723 = -193d0/1344d0, 
     & B5724 = 251d0/1512d0, 
     & B5725 = -425d0/4032d0, 
     & B5726 = 29d0/840d0, 
     & B5727 = -271d0/60480d0, 
     & B5731 = -271d0/60480d0, 
     & B5732 = 97d0/5040d0, 
     & B5733 = -13d0/448d0, 
     & B5734 = 5d0/378d0, 
     & B5735 = 37d0/4032d0, 
     & B5736 = -19d0/1680d0, 
     & B5737 = 191d0/60480d0, 
     & B5741 = 191d0/60480d0, 
     & B5742 = -67d0/2520d0, 
     & B5743 = 115d0/1344d0, 
     & B5744 = -211d0/1512d0, 
     & B5745 = 499d0/4032d0, 
     & B5746 = -2d0/35d0, 
     & B5747 = 653d0/60480d0) 
 
c--------- B5 - B7 
c 
c[1997/60480,  -113/630, 1619/4032, -715/1512, 1241/4032, -263/2520,  863/60480] 
c[-733/60480,    41/630, -193/1344,  251/1512, -425/4032,    29/840, -271/60480] 
c[-271/60480,   97/5040,   -13/448,     5/378,   37/4032,  -19/1680,  191/60480] 
c[ 191/60480,  -67/2520,  115/1344, -211/1512,  499/4032,     -2/35,  653/60480] 
c[-271/60480,    29/840, -425/4032,  251/1512, -193/1344,    41/630, -733/60480] 
c[ 863/60480, -263/2520, 1241/4032, -715/1512, 1619/4032,  -113/630, 1997/60480] 
c 
c 
       PARAMETER( L521 =  3.668340831928216D-01, 
     & L531 =  2.477905683677308D-03, 
     & L541 = -1.919925047010838D-03, 
     & L551 =  2.218385581234200D-03, 
     & L561 = -5.442189351609260D-03, 
     & L532 =  3.216639533696728D-01, 
     & L542 =  1.231925763308414D-03, 
     & L552 =  7.841944627374794D-03, 
     & L562 =  1.002485104590053D-03, 
     & L543 =  3.375100828961925D-01, 
     & L553 = -2.614300734741796D-04, 
     & L563 =  1.066631182323580D-03, 
     & L554 =  3.523137378783708D-01, 
     & L564 = -3.596681121610224D-04, 
     & L565 =  3.617716171655064D-01, 
     & CP51 =   3D0/160D0, 
     & CP52 = -11D0/1440D0) 
 
C--------- A7, B7, B57, B58 Cp7 := MATRICES DEFINING GAM7 
 
        PARAMETER( B711 = 19087d0/60480d0, 
     & B712 = 2713d0/2520d0, 
     & B713 = -15487d0/20160d0, 
     & B714 = 586d0/945d0, 
     & B715 = -6737d0/20160d0, 
     & B716 = 263d0/2520d0, 
     & B717 = -863d0/60480d0, 
     & B721 = -863d0/60480d0, 
     & B722 = 349d0/840d0, 
     & B723 = 5221d0/6720d0, 
     & B724 = -254d0/945d0, 
     & B725 = 811d0/6720d0, 
     & B726 = -29d0/840d0, 
     & B727 = 271d0/60480d0, 
     & B731 = 271d0/60480d0, 
     & B732 = -23d0/504d0, 
     & B733 = 10273d0/20160d0, 
     & B734 = 586d0/945d0, 
     & B735 = -2257d0/20160d0, 
     & B736 = 67d0/2520d0, 
     & B737 = -191d0/60480d0) 
 
C 
C--------- THE LAST THREE ROWS ARE THE REVERSE OF THE FIRST THREE 
C 
c B79 = 
c 
c[ 75203/3628800, -280187/1814400, 129781/259200, -238937/259200, 27289/25920, -197687/259200,  88531/259200, -156437/1814400,  33953/3628800] 
c[-17827/3628800,   66043/1814400, -30389/259200,   55513/259200, -6281/25920,   44983/259200, -19859/259200,   34453/1814400,  -7297/3628800] 
c[  8963/3628800,  -32987/1814400,  15061/259200,  -27257/259200,  3049/25920,  -21527/259200,   9331/259200,  -15797/1814400,   3233/3628800] 
c[  3233/3628800,  -10067/1814400,   3601/259200,   -4337/259200,     23/3240,    1393/259200,  -2129/259200,    7123/1814400,  -2497/3628800] 
c[ -2497/3628800,   12853/1814400,  -7859/259200,   18583/259200, -2681/25920,   24313/259200, -13589/259200,   30043/1814400,  -8227/3628800] 
c[  3233/3628800,  -15797/1814400,   9331/259200,  -21527/259200,  3049/25920,  -27257/259200,  15061/259200,  -32987/1814400,   8963/3628800] 
c[ -7297/3628800,   34453/1814400, -19859/259200,   44983/259200, -6281/25920,   55513/259200, -30389/259200,   66043/1814400, -17827/3628800] 
c[ 33953/3628800, -156437/1814400,  88531/259200, -197687/259200, 27289/25920, -238937/259200, 129781/259200, -280187/1814400,  75203/3628800] 
c 
       PARAMETER( B7911 =   75203d0/3628800d0, 
     & B7912 = -280187d0/1814400d0, 
     & B7913 =  129781d0/259200d0, 
     & B7914 = -238937d0/259200d0, 
     & B7915 =   27289d0/25920d0, 
     & B7916 = -197687d0/259200d0, 
     & B7917 =   88531d0/259200d0, 
     & B7918 = -156437d0/1814400d0, 
     & B7919 =   33953d0/3628800d0) 
 
       PARAMETER(B7921 = -17827d0/3628800d0, 
     & B7922 =  66043d0/1814400d0, 
     & B7923 = -30389d0/259200d0, 
     & B7924 =  55513d0/259200d0, 
     & B7925 =  -6281d0/25920d0, 
     & B7926 =  44983d0/259200d0, 
     & B7927 = -19859d0/259200d0, 
     & B7928 =  34453d0/1814400d0, 
     & B7929 =  -7297d0/3628800d0) 
 
       PARAMETER(B7931 =   8963d0/3628800d0, 
     & B7932 = -32987d0/1814400d0, 
     & B7933 =  15061d0/259200d0, 
     & B7934 = -27257d0/259200d0, 
     & B7935 =   3049d0/25920d0, 
     & B7936 = -21527d0/259200d0, 
     & B7937 =   9331d0/259200d0, 
     & B7938 = -15797d0/1814400d0, 
     & B7939 =   3233d0/3628800d0) 
 
       PARAMETER(B7941 =   3233d0/3628800d0, 
     & B7942 = -10067d0/1814400d0, 
     & B7943 =   3601d0/259200d0, 
     & B7944 =  -4337d0/259200d0, 
     & B7945 =     23d0/3240d0, 
     & B7946 =   1393d0/259200d0, 
     & B7947 =  -2129d0/259200d0, 
     & B7948 =   7123d0/1814400d0, 
     & B7949 =  -2497d0/3628800d0) 
 
 
       PARAMETER(B7951 =  -2497d0/3628800d0, 
     & B7952 =  12853d0/1814400d0, 
     & B7953 =  -7859d0/259200d0, 
     & B7954 =  18583d0/259200d0, 
     & B7955 =  -2681d0/25920d0, 
     & B7956 =  24313d0/259200d0, 
     & B7957 = -13589d0/259200d0, 
     & B7958 =  30043d0/1814400d0, 
     & B7959 =  -8227d0/3628800d0) 
 
       PARAMETER(L721 = 3.023839891568610D-01, 
     & L731 = 3.201698610574002D-05, 
     & L741 = 4.193101163680004D-04, 
     & L751 = 1.686924996069667D-04, 
     & L761 = 4.806043527549464D-05, 
     & L771 = 3.598347048026785D-06, 
     & L781 = 7.892534649789167D-04, 
     & L732 = 2.559868364091398D-01, 
     & L742 = 1.336896192287030D-04, 
     & L752 = 3.080994719931695D-03, 
     & L762 = 1.457177183563680D-04, 
     & L772 = 9.259360509484074D-04, 
     & L782 = 2.397658879381223D-04, 
     & L743 = 2.639734712170458D-01, 
     & L753 = 1.734338929611258D-04, 
     & L763 = 6.704398263264620D-03, 
     & L773 = 4.559927214651730D-05, 
     & L783 = 6.396418554053151D-05, 
     & L754 = 2.817729090368562D-01, 
     & L764 = 2.877761776030408D-04, 
     & L774 = 1.810919475521773D-04, 
     & L784 = 1.009049833235848D-03, 
     & L765 = 2.993040718034231D-01, 
     & L775 = 2.009850887505898D-03, 
     & L785 = 1.748065618845750D-03, 
     & L776 = 3.150349043479135D-01, 
     & L786 = 3.243816792609449D-05, 
     & L787 = 3.271307059448932D-01) 
 
       PARAMETER(CP71 = 103D0/9061D0, 
     & CP72 = -13D0/4480D0, 
     & CP73 =  67D0/42431D0) 
 
C--------- A8, B8, B86, B810 Cp8 := MATRICES DEFINING GAM9 
 
       PARAMETER(B911 = 1070017D0/3628800D0, 
     & B912 = 2233547D0/1814400D0, 
     & B913 = -2302297D0/1814400D0, 
     & B914 = 2797679D0/1814400D0, 
     & B915 = -31457D0/22680D0, 
     & B916 = 1573169D0/1814400D0, 
     & B917 = -645607D0/1814400D0, 
     & B918 = 156437D0/1814400D0, 
     & B919 = -33953D0/3628800D0, 
     & B921 = -33953D0/3628800D0, 
     & B922 = 687797D0/1814400D0, 
     & B923 =  1622393D0/1814400D0, 
     & B924 = -876271D0/1814400D0, 
     & B925 =   8233D0/22680D0, 
     & B926 =    -377521D0/1814400D0, 
     & B927 =   147143D0/1814400D0, 
     & B928 =  -34453D0/1814400D0, 
     & B929 =   7297D0/3628800D0, 
     & B931 = 7297D0/3628800D0, 
     & B932 =  -49813D0/1814400D0, 
     & B933 =  819143D0/1814400D0, 
     & B934 =  1315919D0/1814400D0, 
     & B935 = -5207D0/22680D0, 
     & B936 =  198929D0/1814400D0, 
     & B937 =  -71047D0/1814400D0, 
     & B938 =  15797D0/1814400D0, 
     & B939 = -3233D0/3628800D0, 
     & B941 = -3233D0/3628800D0, 
     & B942 = 18197D0/1814400D0, 
     & B943 =  -108007D0/1814400D0, 
     & B944 =  954929D0/1814400D0, 
     & B945 = 13903D0/22680D0, 
     & B946 = -212881D0/1814400D0, 
     & B947 =  63143D0/1814400D0, 
     & B948 =  -12853D0/1814400D0, 
     & B949 =  2497D0/3628800D0) 
 
 
c   B910=B9-B10; 
c               prime 4 righe : la 5 e' uguale alla 4; 
c     le ultime 4 uguali alle prime 4 negate (simmetriche); 
c     le prime  5 colonne  : le altre sono simmetriche e negate; 
c 
c[ 8183/1036800, -8183/115200,   8183/28800, -57281/86400,  57281/57600, -57281/57600,  57281/86400,  -8183/28800,  8183/115200, -8183/1036800] 
c[  -425/290304,    425/32256,    -425/8064,     425/3456,    -425/2304,     425/2304,    -425/3456,     425/8064,   -425/32256,    425/290304] 
c[      7/12800,    -63/12800,      63/3200,    -147/3200,     441/6400,    -441/6400,     147/3200,     -63/3200,     63/12800,      -7/12800] 
c[-2497/7257600,  2497/806400, -2497/201600,   2497/86400,  -2497/57600,   2497/57600,  -2497/86400,  2497/201600, -2497/806400,  2497/7257600] 
c[-2497/7257600,  2497/806400, -2497/201600,   2497/86400,  -2497/57600,   2497/57600,  -2497/86400,  2497/201600, -2497/806400,  2497/7257600] 
c[ 2497/7257600, -2497/806400,  2497/201600,  -2497/86400,   2497/57600,  -2497/57600,   2497/86400, -2497/201600,  2497/806400, -2497/7257600] 
c[     -7/12800,     63/12800,     -63/3200,     147/3200,    -441/6400,     441/6400,    -147/3200,      63/3200,    -63/12800,       7/12800] 
c[   425/290304,   -425/32256,     425/8064,    -425/3456,     425/2304,    -425/2304,     425/3456,    -425/8064,    425/32256,   -425/290304] 
c[-8183/1036800,  8183/115200,  -8183/28800,  57281/86400, -57281/57600,  57281/57600, -57281/86400,   8183/28800, -8183/115200,  8183/1036800] 
 
       PARAMETER(B91011 =   8183d0/1036800d0, 
     & B91012 =  -8183d0/115200d0, 
     & B91013 =   8183d0/28800d0, 
     & B91014 = -57281d0/86400d0, 
     & B91015 =  57281d0/57600d0) 
 
       PARAMETER(B91021 =   -425d0/290304d0, 
     & B91022 =    425d0/32256d0, 
     & B91023 =   -425d0/8064d0, 
     & B91024 =    425d0/3456d0, 
     & B91025 =   -425d0/2304d0) 
 
       PARAMETER(B91031 =     7d0/12800d0, 
     & B91032 =   -63d0/12800d0, 
     & B91033 =    63d0/3200d0, 
     & B91034 =  -147d0/3200d0, 
     & B91035 =   441d0/6400d0) 
 
       PARAMETER(B91041 = -2497d0/7257600d0, 
     & B91042 =  2497d0/806400d0, 
     & B91043 = -2497d0/201600d0, 
     & B91044 =  2497d0/86400d0, 
     & B91045 = -2497d0/57600d0) 
 
       PARAMETER(Cp91 =  7.892554012345216d-03, 
     & Cp92 = -1.463982583774219d-03, 
     & Cp93 =  5.468749999999983d-04, 
     & Cp94 = -3.440531305114634d-04) 
 
C--------- THE OTHERS ARE THE SAME WITH CHANGED SIGN 
 
       PARAMETER( L921 = 2.590721934790442d-01, 
     & L932   =   2.077575545359853d-01, 
     & L943   =   2.032874698558627d-01, 
     & L954   =   2.036384888660128d-01, 
     & L965   =   2.039599505779785d-01, 
     & L976   =   2.034044409161703d-01, 
     & L987   =   2.017245408702437d-01, 
     & L998   =   1.986549276295617d-01) 
 
C
C   INPUT VARIABLES
C------------------------------------
      INTEGER R, IT, IJOB, IPIV(R), LDLU, IPAR(*)
      DOUBLE PRECISION  H, SCAL(R), TP(*), ERRNEWT0,    
     &                   LU(LDLU,R),RPAR(*)
C
C   INPUT/OUTPUT VARIABLES
C------------------------------------
       INTEGER NFCN
       DOUBLE PRECISION  ERRNEWT, TETAK0, YP(R,*), FP(R,*), F1(R,*),
     &                   DN(R) 
       LOGICAL TER
C
C   LOCAL VARIABLES
C------------------------------------
      INTEGER  J 
      DOUBLE PRECISION  ERRVJ, SUM
C
C   EXTERNAL FUNCTIONS
C------------------------------------

      EXTERNAL FCN

C
C   EXECUTABLE STATEMENTS
C---------------------------------
C          
C--------- ONE STEP OF THE ITERATION PROCEDURE
      TER = .FALSE.
      DO J=1,R
         SUM= B711*FP(J,1)+B712*FP(J,2)+B713*FP(J,3)
     &       +B714*FP(J,4)+B715*FP(J,5)+B716*FP(J,6)+B717*FP(J,7)
         DN(J) =YP(J,2)-YP(J,1)-H*SUM
      END DO
          CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
          ERRVJ = 0D0
          DO J=1,R
            YP(J,2)=YP(J,2) - DN(J)
            SUM = (DN(J)/SCAL(J))
            ERRVJ =  ERRVJ + SUM*SUM
          END DO
          ERRVJ = DSQRT(ERRVJ/R)
          ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
          IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
              TER = .TRUE.
              RETURN
          END IF
          CALL FCN(R,TP(2),YP(1,2),F1(1,1), RPAR,IPAR)
          NFCN = NFCN + 1

      DO J=1,R
         SUM = L721*(F1(J,1)-FP(J,2))+B721*FP(J,1)+B722*FP(J,2)+
     &B723*FP(J,3)+B724*FP(J,4)+B725*FP(J,5)+B726*FP(J,6)+B727*FP(J,7)
         DN(J) = YP(J,3)-YP(J,2)-H*SUM
      END DO
              CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
              ERRVJ = 0D0
              DO J=1,R
                YP(J,3)=YP(J,3)-DN(J)
               SUM = (DN(J)/SCAL(J))
               ERRVJ =  ERRVJ + SUM*SUM
              END DO
              ERRVJ = DSQRT(ERRVJ/R)
              ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
              IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
                 TER = .TRUE.
                 RETURN
              END IF
              CALL FCN(R,TP(3),YP(1,3),F1(1,2), RPAR,IPAR)
              NFCN = NFCN + 1

       DO J=1,R
          SUM = +L731*(F1(J,1)-FP(J,2))+L732*(F1(J,2)-FP(J,3))
     &                 +B731*FP(J,1)+B732*FP(J,2)+B733*FP(J,3)
     &    +B734*FP(J,4)+B735*FP(J,5)+B736*FP(J,6)+B737*FP(J,7)
          DN(J) =YP(J,4)-YP(J,3)-H*SUM
       END DO
              CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
              ERRVJ = 0D0
              DO J=1,R
                YP(J,4)=YP(J,4)-DN(J)
               SUM = (DN(J)/SCAL(J))
               ERRVJ =  ERRVJ + SUM*SUM
              END DO
              ERRVJ = DSQRT(ERRVJ/R)
              ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
              IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
                 TER = .TRUE.
                 RETURN
              END IF
              CALL FCN(R,TP(4),YP(1,4),F1(1,3), RPAR,IPAR)
              NFCN = NFCN + 1

      DO J=1,R
         SUM =L741*(F1(J,1)-FP(J,2))+L742*(F1(J,2)-FP(J,3))
     &                              +L743*(F1(J,3)-FP(J,4))
     &              +B731*FP(J,2)+B732*FP(J,3)+B733*FP(J,4)
     & +B734*FP(J,5)+B735*FP(J,6)+B736*FP(J,7)+B737*FP(J,8)
         DN(J) =YP(J,5)-YP(J,4)-H*SUM
      END DO
              CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
              ERRVJ = 0D0
              DO J=1,R
                YP(J,5)=YP(J,5)-DN(J)
               SUM = (DN(J)/SCAL(J))
               ERRVJ =  ERRVJ + SUM*SUM
              END DO
              ERRVJ = DSQRT(ERRVJ/R)
              ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
              IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
                 TER = .TRUE.
                 RETURN
              END IF

              CALL FCN(R,TP(5),YP(1,5),F1(1,4), RPAR,IPAR)
              NFCN = NFCN + 1

       DO J=1,R
          SUM = L751*(F1(J,1)-FP(J,2))+L752*(F1(J,2)-FP(J,3))
     &         +L753*(F1(J,3)-FP(J,4))+L754*(F1(J,4)-FP(J,5))
     &                +B731*FP(J,3)+B732*FP(J,4)+B733*FP(J,5)
     &   +B734*FP(J,6)+B735*FP(J,7)+B736*FP(J,8)+B737*FP(J,9)
          DN(J)=YP(J,6)-YP(J,5)-H*SUM
       END DO
       CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
              ERRVJ = 0D0
              DO J=1,R
                YP(J,6)=YP(J,6)-DN(J)
               SUM = (DN(J)/SCAL(J))
               ERRVJ =  ERRVJ + SUM*SUM
              END DO
              ERRVJ = DSQRT(ERRVJ/R)
              ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
              IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
                 TER = .TRUE.
                 RETURN
              end if
              CALL FCN(R,TP(6),YP(1,6),F1(1,5), RPAR,IPAR)
              NFCN = NFCN + 1

      DO J=1,R
         SUM = L761*(F1(J,1)-FP(J,2))+L762*(F1(J,2)-FP(J,3))
     &        +L763*(F1(J,3)-FP(J,4))+L764*(F1(J,4)-FP(J,5))
     &                               +L765*(F1(J,5)-FP(J,6))
     &               +B737*FP(J,3)+B736*FP(J,4)+B735*FP(J,5)
     &  +B734*FP(J,6)+B733*FP(J,7)+B732*FP(J,8)+B731*FP(J,9)
         DN(J) = YP(J,7)-YP(J,6)-H*SUM
      END DO
              CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
              ERRVJ = 0D0
              DO J=1,R
                YP(J,7)=YP(J,7)-DN(J)
               SUM = (DN(J)/SCAL(J))
               ERRVJ =  ERRVJ + SUM*SUM
              END DO
              ERRVJ = DSQRT(ERRVJ/R)
              ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
              IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
                 TER = .TRUE.
                 RETURN
              END IF

              CALL FCN(R,TP(7),YP(1,7),F1(1,6), RPAR,IPAR)
              NFCN = NFCN + 1

       DO J=1,R
          SUM = L771*(F1(J,1)-FP(J,2))+L772*(F1(J,2)-FP(J,3))
     &         +L773*(F1(J,3)-FP(J,4))+L774*(F1(J,4)-FP(J,5))
     &         +L775*(F1(J,5)-FP(J,6))+L776*(F1(J,6)-FP(J,7))
     &                +B727*FP(J,3)+B726*FP(J,4)+B725*FP(J,5)
     &   +B724*FP(J,6)+B723*FP(J,7)+B722*FP(J,8)+B721*FP(J,9)
          DN(J) = YP(J,8)-YP(J,7)-H*SUM
       END DO
              CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
              ERRVJ = 0D0
              DO J=1,R
                YP(J,8)=YP(J,8)-DN(J)
               SUM = (DN(J)/SCAL(J))
               ERRVJ =  ERRVJ + SUM*SUM
              END DO
              ERRVJ = DSQRT(ERRVJ/R)
              ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
              IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
                 TER = .TRUE.
                 RETURN
              END IF
              CALL FCN(R,TP(8),YP(1,8),F1(1,7), RPAR,IPAR)
              NFCN = NFCN + 1


      DO J=1,R
         SUM = L781*(F1(J,1)-FP(J,2))+L782*(F1(J,2)-FP(J,3))
     &        +L783*(F1(J,3)-FP(J,4))+L784*(F1(J,4)-FP(J,5))
     &        +L785*(F1(J,5)-FP(J,6))+L786*(F1(J,6)-FP(J,7))
     &                               +L787*(F1(J,7)-FP(J,8))
     &               +B717*FP(J,3)+B716*FP(J,4)+B715*FP(J,5)
     &  +B714*FP(J,6)+B713*FP(J,7)+B712*FP(J,8)+B711*FP(J,9)
         DN(J) = YP(J,9)-YP(J,8)-H*SUM
      END DO
              CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
              ERRVJ = 0D0
              DO J=1,R
                YP(J,9)=YP(J,9)-DN(J)
               SUM = (DN(J)/SCAL(J))
               ERRVJ =  ERRVJ + SUM*SUM
              END DO
              ERRVJ = DSQRT(ERRVJ/R)
              ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
              IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
                 TER = .TRUE.
                 RETURN
              END IF
              CALL FCN(R,TP(9),YP(1,9),F1(1,8), RPAR,IPAR)
              NFCN = NFCN + 1

               DO J=1,R
                 FP(J,2) = F1(J,1)
                 FP(J,3) = F1(J,2)
                 FP(J,4) = F1(J,3)
                 FP(J,5) = F1(J,4)
                 FP(J,6) = F1(J,5)
                 FP(J,7) = F1(J,6)
                 FP(J,8) = F1(J,7)
                 FP(J,9) = F1(J,8)
               END DO

          RETURN
          END
C
C  SUBROUTINE TERMNOT9  (ORDER 9)
C
      SUBROUTINE  TERMNOT9(R,FCN,H,IT,DN, F1,FP,YP,TP,NFCN,
     &  ERRNEWT,ERRNEWT0,TETAK0,LU, LDLU,IPIV, SCAL,IJOB,TER,
     &  RPAR,IPAR)
       IMPLICIT NONE
C
C   INCLUDE
C------------------------------------
       DOUBLE PRECISION 
     & B3511, B3512, B3513, B3514, B3515, B3521, B3522, B3523, B3524, 
     & B3525, B3531, B3532, B3533, B3534, B3535, B3541, B3542, B3543, 
     & B3544, B3545, 
     & L321, L331, L332, L341, L342, L343, B311, B312, B313, 
     & B511, B512, B513, B514, B515, B521, B522, B523, B524, B525, 
     & L521, L531, L541, L551, L561, L532, L542, L552, L562, L543, 
     & L553, L563, L554, L564, L565, 
     & B5711, B5712, B5713, B5714, B5715, B5716, B5717, B5721, B5722, 
     & B5723, B5724, B5725, B5726, B5727, B5731, B5732, B5733, B5734, 
     & B5735, B5736, B5737, B5741, B5742, B5743, B5744, B5745, B5746, 
     & B5747, 
     & B711, B712, B713, B714, B715, B716, B717, B721, B722, B723, 
     & B724, B725, B726, B727, B731, B732, B733, B734, B735, B736, B737, 
     & B7911, B7912, B7913, B7914, B7915, B7916, B7917, B7918, B7919, 
     & B7921, B7922, B7923, B7924, B7925, B7926, B7927, B7928, B7929, 
     & B7931, B7932, B7933, B7934, B7935, B7936, B7937, B7938, B7939, 
     & B7941, B7942, B7943, B7944, B7945, B7946, B7947, B7948, B7949, 
     & B7951, B7952, B7953, B7954, B7955, B7956, B7957, B7958, B7959, 
     & L721, L731, L741, L751, L761, L771, L781, L732, L742, L752, 
     & L762, L772, L782, L743, L753, L763, L773, L783, L754, L764, 
     & L774, L784, L765, L775, L785, L776, L786, L787, 
     & CP31, CP51, CP52, CP71, CP72, CP73, 
     & CP91, CP92, CP93, CP94, 
     & B911, B912, B913, B914, B915, B916, B917, B918, B919, B921, 
     & B922, B923, B924, B925, B926, B927, B928, B929, B931, B932, 
     & B933, B934, B935, B936, B937, B938, B939, B941, B942, B943, 
     & B944, B945, B946, B947, B948, B949, 
     & L921, L932, L943, L954, L965, L976, L987, L998, 
     & B91011, B91012, B91013, B91014, B91015, B91021, B91022, B91023, 
     & B91024, B91025, B91031, B91032, B91033, B91034, B91035, B91041, 
     & B91042, B91043, B91044, B91045 
       PARAMETER( B311 = 5d0/12d0, 
     &            B312 = 8d0/12d0, 
     &            B313 = -1d0/12d0) 
 
       PARAMETER( Cp31 = 1d0/24d0) 
ccccc L3 dimensione 4 
       PARAMETER( L321 =  4.012395124208693d-01, 
     &            L331 =  8.819910099032529d-04, 
     &            L341 =  1.728116022258560d-04, 
     &            L332 =  3.680857287181573d-01, 
     &            L342 =  1.635381132422046d-03, 
     &            L343 =  3.688541178419062d-01) 
c--------- B3-B5 OF DIMENSION  4 
       PARAMETER( B3511 = 49d0/720d0, 
     &            B3512 = -83d0/360d0, 
     &            B3513 = 17d0/60d0, 
     &            B3514 = -53d0/360d0, 
     &            B3515 = 19d0/720d0, 
     &            B3521 = 19d0/720d0, 
     &            B3522 = -23d0/360d0, 
     &            B3523 = 1d0/30d0, 
     & B3524 = 7d0/360d0, 
     & B3525 = -11d0/720d0, 
     & B3531 = -11d0/720d0, 
     & B3532 = 37d0/360d0, 
     & B3533 = -13d0/60d0, 
     & B3534 = 67d0/360d0, 
     & B3535 = -41d0/720d0, 
     & B3541 = 19d0/720d0, 
     & B3542 = -53d0/360d0, 
     & B3543 =  17d0/60d0, 
     & B3544 = -83d0/360d0, 
     & B3545 =  49d0/720d0) 
 
c[ 49/720, -83/360,  17/60, -53/360,  19/720] 
c[ 19/720, -23/360,   1/30,   7/360, -11/720] 
c[-11/720,  37/360, -13/60,  67/360, -41/720] 
c[ 19/720, -53/360,  17/60, -83/360,  49/720] 
 
C--------- A5, B5, B53, B56 Cp5 := MATRICES DEFINING  GAM5 
 
       PARAMETER( B511 = 251d0/720d0, 
     & B512 = 323d0/360d0, 
     & B513 = - 11d0/30d0, 
     & B514 = 53d0/360d0, 
     & B515 = -19d0/720d0, 
     & B521 = -19d0/720d0, 
     & B522 =  173d0/360d0, 
     & B523 = 19d0/30d0, 
     & B524 = -37d0/360d0, 
     & B525 = 11d0/720d0) 
 
       PARAMETER( B5711 = 1997d0/60480d0, 
     & B5712 = -113d0/630d0, 
     & B5713 = 1619d0/4032d0, 
     & B5714 = -715d0/1512d0, 
     & B5715 = 1241d0/4032d0, 
     & B5716 = -263d0/2520d0, 
     & B5717 = 863d0/60480d0, 
     & B5721 = -733d0/60480d0, 
     & B5722 = 41d0/630d0, 
     & B5723 = -193d0/1344d0, 
     & B5724 = 251d0/1512d0, 
     & B5725 = -425d0/4032d0, 
     & B5726 = 29d0/840d0, 
     & B5727 = -271d0/60480d0, 
     & B5731 = -271d0/60480d0, 
     & B5732 = 97d0/5040d0, 
     & B5733 = -13d0/448d0, 
     & B5734 = 5d0/378d0, 
     & B5735 = 37d0/4032d0, 
     & B5736 = -19d0/1680d0, 
     & B5737 = 191d0/60480d0, 
     & B5741 = 191d0/60480d0, 
     & B5742 = -67d0/2520d0, 
     & B5743 = 115d0/1344d0, 
     & B5744 = -211d0/1512d0, 
     & B5745 = 499d0/4032d0, 
     & B5746 = -2d0/35d0, 
     & B5747 = 653d0/60480d0) 
 
c--------- B5 - B7 
c 
c[1997/60480,  -113/630, 1619/4032, -715/1512, 1241/4032, -263/2520,  863/60480] 
c[-733/60480,    41/630, -193/1344,  251/1512, -425/4032,    29/840, -271/60480] 
c[-271/60480,   97/5040,   -13/448,     5/378,   37/4032,  -19/1680,  191/60480] 
c[ 191/60480,  -67/2520,  115/1344, -211/1512,  499/4032,     -2/35,  653/60480] 
c[-271/60480,    29/840, -425/4032,  251/1512, -193/1344,    41/630, -733/60480] 
c[ 863/60480, -263/2520, 1241/4032, -715/1512, 1619/4032,  -113/630, 1997/60480] 
c 
c 
       PARAMETER( L521 =  3.668340831928216D-01, 
     & L531 =  2.477905683677308D-03, 
     & L541 = -1.919925047010838D-03, 
     & L551 =  2.218385581234200D-03, 
     & L561 = -5.442189351609260D-03, 
     & L532 =  3.216639533696728D-01, 
     & L542 =  1.231925763308414D-03, 
     & L552 =  7.841944627374794D-03, 
     & L562 =  1.002485104590053D-03, 
     & L543 =  3.375100828961925D-01, 
     & L553 = -2.614300734741796D-04, 
     & L563 =  1.066631182323580D-03, 
     & L554 =  3.523137378783708D-01, 
     & L564 = -3.596681121610224D-04, 
     & L565 =  3.617716171655064D-01, 
     & CP51 =   3D0/160D0, 
     & CP52 = -11D0/1440D0) 
 
C--------- A7, B7, B57, B58 Cp7 := MATRICES DEFINING GAM7 
 
        PARAMETER( B711 = 19087d0/60480d0, 
     & B712 = 2713d0/2520d0, 
     & B713 = -15487d0/20160d0, 
     & B714 = 586d0/945d0, 
     & B715 = -6737d0/20160d0, 
     & B716 = 263d0/2520d0, 
     & B717 = -863d0/60480d0, 
     & B721 = -863d0/60480d0, 
     & B722 = 349d0/840d0, 
     & B723 = 5221d0/6720d0, 
     & B724 = -254d0/945d0, 
     & B725 = 811d0/6720d0, 
     & B726 = -29d0/840d0, 
     & B727 = 271d0/60480d0, 
     & B731 = 271d0/60480d0, 
     & B732 = -23d0/504d0, 
     & B733 = 10273d0/20160d0, 
     & B734 = 586d0/945d0, 
     & B735 = -2257d0/20160d0, 
     & B736 = 67d0/2520d0, 
     & B737 = -191d0/60480d0) 
 
C 
C--------- THE LAST THREE ROWS ARE THE REVERSE OF THE FIRST THREE 
C 
c B79 = 
c 
c[ 75203/3628800, -280187/1814400, 129781/259200, -238937/259200, 27289/25920, -197687/259200,  88531/259200, -156437/1814400,  33953/3628800] 
c[-17827/3628800,   66043/1814400, -30389/259200,   55513/259200, -6281/25920,   44983/259200, -19859/259200,   34453/1814400,  -7297/3628800] 
c[  8963/3628800,  -32987/1814400,  15061/259200,  -27257/259200,  3049/25920,  -21527/259200,   9331/259200,  -15797/1814400,   3233/3628800] 
c[  3233/3628800,  -10067/1814400,   3601/259200,   -4337/259200,     23/3240,    1393/259200,  -2129/259200,    7123/1814400,  -2497/3628800] 
c[ -2497/3628800,   12853/1814400,  -7859/259200,   18583/259200, -2681/25920,   24313/259200, -13589/259200,   30043/1814400,  -8227/3628800] 
c[  3233/3628800,  -15797/1814400,   9331/259200,  -21527/259200,  3049/25920,  -27257/259200,  15061/259200,  -32987/1814400,   8963/3628800] 
c[ -7297/3628800,   34453/1814400, -19859/259200,   44983/259200, -6281/25920,   55513/259200, -30389/259200,   66043/1814400, -17827/3628800] 
c[ 33953/3628800, -156437/1814400,  88531/259200, -197687/259200, 27289/25920, -238937/259200, 129781/259200, -280187/1814400,  75203/3628800] 
c 
       PARAMETER( B7911 =   75203d0/3628800d0, 
     & B7912 = -280187d0/1814400d0, 
     & B7913 =  129781d0/259200d0, 
     & B7914 = -238937d0/259200d0, 
     & B7915 =   27289d0/25920d0, 
     & B7916 = -197687d0/259200d0, 
     & B7917 =   88531d0/259200d0, 
     & B7918 = -156437d0/1814400d0, 
     & B7919 =   33953d0/3628800d0) 
 
       PARAMETER(B7921 = -17827d0/3628800d0, 
     & B7922 =  66043d0/1814400d0, 
     & B7923 = -30389d0/259200d0, 
     & B7924 =  55513d0/259200d0, 
     & B7925 =  -6281d0/25920d0, 
     & B7926 =  44983d0/259200d0, 
     & B7927 = -19859d0/259200d0, 
     & B7928 =  34453d0/1814400d0, 
     & B7929 =  -7297d0/3628800d0) 
 
       PARAMETER(B7931 =   8963d0/3628800d0, 
     & B7932 = -32987d0/1814400d0, 
     & B7933 =  15061d0/259200d0, 
     & B7934 = -27257d0/259200d0, 
     & B7935 =   3049d0/25920d0, 
     & B7936 = -21527d0/259200d0, 
     & B7937 =   9331d0/259200d0, 
     & B7938 = -15797d0/1814400d0, 
     & B7939 =   3233d0/3628800d0) 
 
       PARAMETER(B7941 =   3233d0/3628800d0, 
     & B7942 = -10067d0/1814400d0, 
     & B7943 =   3601d0/259200d0, 
     & B7944 =  -4337d0/259200d0, 
     & B7945 =     23d0/3240d0, 
     & B7946 =   1393d0/259200d0, 
     & B7947 =  -2129d0/259200d0, 
     & B7948 =   7123d0/1814400d0, 
     & B7949 =  -2497d0/3628800d0) 
 
 
       PARAMETER(B7951 =  -2497d0/3628800d0, 
     & B7952 =  12853d0/1814400d0, 
     & B7953 =  -7859d0/259200d0, 
     & B7954 =  18583d0/259200d0, 
     & B7955 =  -2681d0/25920d0, 
     & B7956 =  24313d0/259200d0, 
     & B7957 = -13589d0/259200d0, 
     & B7958 =  30043d0/1814400d0, 
     & B7959 =  -8227d0/3628800d0) 
 
       PARAMETER(L721 = 3.023839891568610D-01, 
     & L731 = 3.201698610574002D-05, 
     & L741 = 4.193101163680004D-04, 
     & L751 = 1.686924996069667D-04, 
     & L761 = 4.806043527549464D-05, 
     & L771 = 3.598347048026785D-06, 
     & L781 = 7.892534649789167D-04, 
     & L732 = 2.559868364091398D-01, 
     & L742 = 1.336896192287030D-04, 
     & L752 = 3.080994719931695D-03, 
     & L762 = 1.457177183563680D-04, 
     & L772 = 9.259360509484074D-04, 
     & L782 = 2.397658879381223D-04, 
     & L743 = 2.639734712170458D-01, 
     & L753 = 1.734338929611258D-04, 
     & L763 = 6.704398263264620D-03, 
     & L773 = 4.559927214651730D-05, 
     & L783 = 6.396418554053151D-05, 
     & L754 = 2.817729090368562D-01, 
     & L764 = 2.877761776030408D-04, 
     & L774 = 1.810919475521773D-04, 
     & L784 = 1.009049833235848D-03, 
     & L765 = 2.993040718034231D-01, 
     & L775 = 2.009850887505898D-03, 
     & L785 = 1.748065618845750D-03, 
     & L776 = 3.150349043479135D-01, 
     & L786 = 3.243816792609449D-05, 
     & L787 = 3.271307059448932D-01) 
 
       PARAMETER(CP71 = 103D0/9061D0, 
     & CP72 = -13D0/4480D0, 
     & CP73 =  67D0/42431D0) 
 
C--------- A8, B8, B86, B810 Cp8 := MATRICES DEFINING GAM9 
 
       PARAMETER(B911 = 1070017D0/3628800D0, 
     & B912 = 2233547D0/1814400D0, 
     & B913 = -2302297D0/1814400D0, 
     & B914 = 2797679D0/1814400D0, 
     & B915 = -31457D0/22680D0, 
     & B916 = 1573169D0/1814400D0, 
     & B917 = -645607D0/1814400D0, 
     & B918 = 156437D0/1814400D0, 
     & B919 = -33953D0/3628800D0, 
     & B921 = -33953D0/3628800D0, 
     & B922 = 687797D0/1814400D0, 
     & B923 =  1622393D0/1814400D0, 
     & B924 = -876271D0/1814400D0, 
     & B925 =   8233D0/22680D0, 
     & B926 =    -377521D0/1814400D0, 
     & B927 =   147143D0/1814400D0, 
     & B928 =  -34453D0/1814400D0, 
     & B929 =   7297D0/3628800D0, 
     & B931 = 7297D0/3628800D0, 
     & B932 =  -49813D0/1814400D0, 
     & B933 =  819143D0/1814400D0, 
     & B934 =  1315919D0/1814400D0, 
     & B935 = -5207D0/22680D0, 
     & B936 =  198929D0/1814400D0, 
     & B937 =  -71047D0/1814400D0, 
     & B938 =  15797D0/1814400D0, 
     & B939 = -3233D0/3628800D0, 
     & B941 = -3233D0/3628800D0, 
     & B942 = 18197D0/1814400D0, 
     & B943 =  -108007D0/1814400D0, 
     & B944 =  954929D0/1814400D0, 
     & B945 = 13903D0/22680D0, 
     & B946 = -212881D0/1814400D0, 
     & B947 =  63143D0/1814400D0, 
     & B948 =  -12853D0/1814400D0, 
     & B949 =  2497D0/3628800D0) 
 
 
c   B910=B9-B10; 
c               prime 4 righe : la 5 e' uguale alla 4; 
c     le ultime 4 uguali alle prime 4 negate (simmetriche); 
c     le prime  5 colonne  : le altre sono simmetriche e negate; 
c 
c[ 8183/1036800, -8183/115200,   8183/28800, -57281/86400,  57281/57600, -57281/57600,  57281/86400,  -8183/28800,  8183/115200, -8183/1036800] 
c[  -425/290304,    425/32256,    -425/8064,     425/3456,    -425/2304,     425/2304,    -425/3456,     425/8064,   -425/32256,    425/290304] 
c[      7/12800,    -63/12800,      63/3200,    -147/3200,     441/6400,    -441/6400,     147/3200,     -63/3200,     63/12800,      -7/12800] 
c[-2497/7257600,  2497/806400, -2497/201600,   2497/86400,  -2497/57600,   2497/57600,  -2497/86400,  2497/201600, -2497/806400,  2497/7257600] 
c[-2497/7257600,  2497/806400, -2497/201600,   2497/86400,  -2497/57600,   2497/57600,  -2497/86400,  2497/201600, -2497/806400,  2497/7257600] 
c[ 2497/7257600, -2497/806400,  2497/201600,  -2497/86400,   2497/57600,  -2497/57600,   2497/86400, -2497/201600,  2497/806400, -2497/7257600] 
c[     -7/12800,     63/12800,     -63/3200,     147/3200,    -441/6400,     441/6400,    -147/3200,      63/3200,    -63/12800,       7/12800] 
c[   425/290304,   -425/32256,     425/8064,    -425/3456,     425/2304,    -425/2304,     425/3456,    -425/8064,    425/32256,   -425/290304] 
c[-8183/1036800,  8183/115200,  -8183/28800,  57281/86400, -57281/57600,  57281/57600, -57281/86400,   8183/28800, -8183/115200,  8183/1036800] 
 
       PARAMETER(B91011 =   8183d0/1036800d0, 
     & B91012 =  -8183d0/115200d0, 
     & B91013 =   8183d0/28800d0, 
     & B91014 = -57281d0/86400d0, 
     & B91015 =  57281d0/57600d0) 
 
       PARAMETER(B91021 =   -425d0/290304d0, 
     & B91022 =    425d0/32256d0, 
     & B91023 =   -425d0/8064d0, 
     & B91024 =    425d0/3456d0, 
     & B91025 =   -425d0/2304d0) 
 
       PARAMETER(B91031 =     7d0/12800d0, 
     & B91032 =   -63d0/12800d0, 
     & B91033 =    63d0/3200d0, 
     & B91034 =  -147d0/3200d0, 
     & B91035 =   441d0/6400d0) 
 
       PARAMETER(B91041 = -2497d0/7257600d0, 
     & B91042 =  2497d0/806400d0, 
     & B91043 = -2497d0/201600d0, 
     & B91044 =  2497d0/86400d0, 
     & B91045 = -2497d0/57600d0) 
 
       PARAMETER(Cp91 =  7.892554012345216d-03, 
     & Cp92 = -1.463982583774219d-03, 
     & Cp93 =  5.468749999999983d-04, 
     & Cp94 = -3.440531305114634d-04) 
 
C--------- THE OTHERS ARE THE SAME WITH CHANGED SIGN 
 
       PARAMETER( L921 = 2.590721934790442d-01, 
     & L932   =   2.077575545359853d-01, 
     & L943   =   2.032874698558627d-01, 
     & L954   =   2.036384888660128d-01, 
     & L965   =   2.039599505779785d-01, 
     & L976   =   2.034044409161703d-01, 
     & L987   =   2.017245408702437d-01, 
     & L998   =   1.986549276295617d-01) 
 
C
C   INPUT VARIABLES
C------------------------------------
      INTEGER R, IT, IJOB, IPIV(R), LDLU, IPAR(*)
      DOUBLE PRECISION  H, SCAL(R), TP(*), ERRNEWT0,    
     &                   LU(LDLU,R),RPAR(*)
C
C   INPUT/OUTPUT VARIABLES ksks: ,1->,*
C------------------------------------
       INTEGER NFCN
       DOUBLE PRECISION  ERRNEWT, TETAK0, YP(R,*), FP(R,*), F1(R,*),
     &                   DN(R) 
       LOGICAL TER
C
C   LOCAL VARIABLES
C------------------------------------
      INTEGER  J 
      DOUBLE PRECISION  ERRVJ, SUM
C
C   EXTERNAL FUNCTIONS
C------------------------------------

      EXTERNAL FCN

C
C   EXECUTABLE STATEMENTS
C---------------------------------
C          
C--------- ONE STEP OF THE ITERATION PROCEDURE
      TER = .FALSE.

      DO J=1,R
         SUM = B911*FP(J,1)+B912*FP(J,2)+B913*FP(J,3)+B914*FP(J,4)
     &+B915*FP(J,5)+B916*FP(J,6)+B917*FP(J,7)+B918*FP(J,8)+B919*FP(J,9)
         DN(J) = YP(J,2)-YP(J,1)-H*SUM
      END DO
              CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
              ERRVJ = 0D0
              DO J=1,R
                YP(J,2)=YP(J,2)-DN(J)
                SUM = (DN(J)/SCAL(J))
                ERRVJ =  ERRVJ + SUM*SUM
              END DO
              ERRVJ = DSQRT(ERRVJ/R)
              ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
              IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
                 TER = .TRUE.
                 RETURN
              END IF
              CALL FCN(R,TP(2),YP(1,2),F1(1,1), RPAR,IPAR)
              NFCN = NFCN + 1


      DO J=1,R
         SUM = L921*(F1(J,1)-FP(J,2))+B921*FP(J,1)+B922*FP(J,2)
     &                  +B923*FP(J,3)+B924*FP(J,4)+B925*FP(J,5)
     &     +B926*FP(J,6)+B927*FP(J,7)+B928*FP(J,8)+B929*FP(J,9)
         DN(J) = YP(J,3)-YP(J,2)-H*SUM
      END DO
              CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
              ERRVJ = 0D0
              DO J=1,R
                YP(J,3)=YP(J,3)-DN(J)
                SUM = (DN(J)/SCAL(J))
                ERRVJ =  ERRVJ + SUM*SUM
              END DO
              ERRVJ = DSQRT(ERRVJ/R)
              ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
              IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
                 TER = .TRUE.
                 RETURN
              END IF

              CALL FCN(R,TP(3),YP(1,3),F1(1,2), RPAR,IPAR)
              NFCN = NFCN + 1

      DO J=1,R
         SUM = L932*(F1(J,2)-FP(J,3))+B931*FP(J,1)
     &     +B932*FP(J,2)+B933*FP(J,3)+B934*FP(J,4)+B935*FP(J,5)
     &     +B936*FP(J,6)+B937*FP(J,7)+B938*FP(J,8)+B939*FP(J,9)
         DN(J) = YP(J,4)-YP(J,3)-H*SUM
      END DO
              CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
              ERRVJ = 0D0
              DO J=1,R
                YP(J,4)=YP(J,4)-DN(J)
                SUM = (DN(J)/SCAL(J))
                ERRVJ =  ERRVJ + SUM*SUM
              END DO
              ERRVJ = DSQRT(ERRVJ/R)
              ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
              IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
                 TER = .TRUE.
                 RETURN
              END IF

              CALL FCN(R,TP(4),YP(1,4),F1(1,3), RPAR,IPAR)
              NFCN = NFCN + 1


       DO J=1,R
          SUM = L943*(F1(J,3)-FP(J,4))+B941*FP(J,1)+B942*FP(J,2)
     &      +B943*FP(J,3)+B944*FP(J,4)+B945*FP(J,5)+B946*FP(J,6)
     &                   +B947*FP(J,7)+B948*FP(J,8)+B949*FP(J,9)
          DN(J) =YP(J,5)-YP(J,4)-H*SUM
       END DO
              CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
              ERRVJ = 0D0
              DO J=1,R
                YP(J,5)=YP(J,5)-DN(J)
                SUM = (DN(J)/SCAL(J))
                ERRVJ =  ERRVJ + SUM*SUM
              END DO
              ERRVJ = DSQRT(ERRVJ/R)
              ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
              IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
                 TER = .TRUE.
                 RETURN
              END IF

              CALL FCN(R,TP(5),YP(1,5),F1(1,4), RPAR,IPAR)
              NFCN = NFCN + 1


       DO J=1,R
          SUM = L954*(F1(J,4)-FP(J,5))+B941*FP(J,2)+B942*FP(J,3)
     &      +B943*FP(J,4)+B944*FP(J,5)+B945*FP(J,6)+B946*FP(J,7)
     &                  +B947*FP(J,8)+B948*FP(J,9)+B949*FP(J,10)
          DN(J) =YP(J,6)-YP(J,5)-H*SUM
       END DO
              CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
              ERRVJ = 0D0
              DO J=1,R
                YP(J,6)=YP(J,6)-DN(J)
                SUM = (DN(J)/SCAL(J))
                ERRVJ =  ERRVJ + SUM*SUM
              END DO
              ERRVJ = DSQRT(ERRVJ/R)
              ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
              IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
                 TER = .TRUE.
                 RETURN
              END IF

              CALL FCN(R,TP(6),YP(1,6),F1(1,5), RPAR,IPAR)
              NFCN = NFCN + 1

       DO J=1,R
          SUM = L965*(F1(J,5)-FP(J,6))+B949*FP(J,2)+B948*FP(J,3)
     &      +B947*FP(J,4)+B946*FP(J,5)+B945*FP(J,6)+B944*FP(J,7)
     &                  +B943*FP(J,8)+B942*FP(J,9)+B941*FP(J,10)
          DN(J) =YP(J,7)-YP(J,6)-H*SUM
       END DO
              CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
              ERRVJ = 0D0
              DO J=1,R
                YP(J,7)=YP(J,7)-DN(J)
                SUM = (DN(J)/SCAL(J))
                ERRVJ =  ERRVJ + SUM*SUM
              END DO
              ERRVJ = DSQRT(ERRVJ/R)
              ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
              IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
                 TER = .TRUE.
                 RETURN
              END IF

              CALL FCN(R,TP(7),YP(1,7),F1(1,6), RPAR,IPAR)
              NFCN = NFCN + 1

       DO J=1,R
          SUM = L976*(F1(J,6)-FP(J,7))+B939*FP(J,2)+B938*FP(J,3)
     &      +B937*FP(J,4)+B936*FP(J,5)+B935*FP(J,6)+B934*FP(J,7)
     &                  +B933*FP(J,8)+B932*FP(J,9)+B931*FP(J,10)
          DN(J) =YP(J,8)-YP(J,7)-H*SUM
       END DO
              CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
              ERRVJ = 0D0
              DO J=1,R
                YP(J,8)=YP(J,8)-DN(J)
                SUM = (DN(J)/SCAL(J))
                ERRVJ =  ERRVJ + SUM*SUM
              END DO
              ERRVJ = DSQRT(ERRVJ/R)
              ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
              IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
                 TER = .TRUE.
                 RETURN
              END IF
              CALL FCN(R,TP(8),YP(1,8),F1(1,7), RPAR,IPAR)
              NFCN = NFCN + 1

       DO J=1,R
          SUM = L987*(F1(J,7)-FP(J,8))+B929*FP(J,2)+B928*FP(J,3)
     &      +B927*FP(J,4)+B926*FP(J,5)+B925*FP(J,6)+B924*FP(J,7)
     &                  +B923*FP(J,8)+B922*FP(J,9)+B921*FP(J,10)
          DN(J) =YP(J,9)-YP(J,8)-H*SUM
       END DO
              CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
              ERRVJ = 0D0
              DO J=1,R
                YP(J,9)=YP(J,9)-DN(J)
                SUM = (DN(J)/SCAL(J))
                ERRVJ =  ERRVJ + SUM*SUM
              END DO
              ERRVJ = DSQRT(ERRVJ/R)
              ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
              IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
                 TER = .TRUE.
                 RETURN
              END IF
              CALL FCN(R,TP(9),YP(1,9),F1(1,8), RPAR,IPAR)
              NFCN = NFCN + 1

       DO J=1,R
          SUM = L998*(F1(J,8)-FP(J,9))+B919*FP(J,2)+B918*FP(J,3)
     &      +B917*FP(J,4)+B916*FP(J,5)+B915*FP(J,6)+B914*FP(J,7)
     &                  +B913*FP(J,8)+B912*FP(J,9)+B911*FP(J,10)
          DN(J) =YP(J,10)-YP(J,9)-H*SUM
       END DO

              CALL  SOLLU(R,LU,LDLU,DN,IPIV,IJOB)
              ERRVJ = 0D0
              DO J=1,R
                YP(J,10)=YP(J,10)-DN(J)
                SUM = (DN(J)/SCAL(J))
                ERRVJ =  ERRVJ + SUM*SUM
              END DO
              ERRVJ = DSQRT(ERRVJ/R)
              ERRNEWT = DMAX1( ERRNEWT, ERRVJ )
              IF ((IT.GE.1).AND.(ERRNEWT/ERRNEWT0 .GT. TETAK0 )) THEN
                 TER = .TRUE.
                 RETURN
              END IF
              CALL FCN(R,TP(10),YP(1,10),F1(1,9), RPAR,IPAR)
              NFCN = NFCN + 1

               DO J=1,R
                 FP(J,2) = F1(J,1)
                 FP(J,3) = F1(J,2)
                 FP(J,4) = F1(J,3)
                 FP(J,5) = F1(J,4)
                 FP(J,6) = F1(J,5)
                 FP(J,7) = F1(J,6)
                 FP(J,8) = F1(J,7)
                 FP(J,9) = F1(J,8)
                 FP(J,10) = F1(J,9)
               END DO
         RETURN
         END


C
C  SUBROUTINE ESTERR
C  ERRORS ESTIMATION.  ERRSAME: THE CURRENT ORDER
C                         ERRUP: GREATER ORDER (THAN ERRSAME)
C                       ERRDOWN: LOWER ORDER
C
      SUBROUTINE  ESTERR(ERRV, ERRSAME, ERRUP, ERRDOWN, FP,
     &     R, H, ORD, DBLK, LU, LDLU,
     &     IPIV, F, DN,SCAL,ORDMAX,ORDMIN,IJOB)
      IMPLICIT NONE
C
C   INCLUDE
C------------------------------------
       DOUBLE PRECISION 
     & B3511, B3512, B3513, B3514, B3515, B3521, B3522, B3523, B3524, 
     & B3525, B3531, B3532, B3533, B3534, B3535, B3541, B3542, B3543, 
     & B3544, B3545, 
     & L321, L331, L332, L341, L342, L343, B311, B312, B313, 
     & B511, B512, B513, B514, B515, B521, B522, B523, B524, B525, 
     & L521, L531, L541, L551, L561, L532, L542, L552, L562, L543, 
     & L553, L563, L554, L564, L565, 
     & B5711, B5712, B5713, B5714, B5715, B5716, B5717, B5721, B5722, 
     & B5723, B5724, B5725, B5726, B5727, B5731, B5732, B5733, B5734, 
     & B5735, B5736, B5737, B5741, B5742, B5743, B5744, B5745, B5746, 
     & B5747, 
     & B711, B712, B713, B714, B715, B716, B717, B721, B722, B723, 
     & B724, B725, B726, B727, B731, B732, B733, B734, B735, B736, B737, 
     & B7911, B7912, B7913, B7914, B7915, B7916, B7917, B7918, B7919, 
     & B7921, B7922, B7923, B7924, B7925, B7926, B7927, B7928, B7929, 
     & B7931, B7932, B7933, B7934, B7935, B7936, B7937, B7938, B7939, 
     & B7941, B7942, B7943, B7944, B7945, B7946, B7947, B7948, B7949, 
     & B7951, B7952, B7953, B7954, B7955, B7956, B7957, B7958, B7959, 
     & L721, L731, L741, L751, L761, L771, L781, L732, L742, L752, 
     & L762, L772, L782, L743, L753, L763, L773, L783, L754, L764, 
     & L774, L784, L765, L775, L785, L776, L786, L787, 
     & CP31, CP51, CP52, CP71, CP72, CP73, 
     & CP91, CP92, CP93, CP94, 
     & B911, B912, B913, B914, B915, B916, B917, B918, B919, B921, 
     & B922, B923, B924, B925, B926, B927, B928, B929, B931, B932, 
     & B933, B934, B935, B936, B937, B938, B939, B941, B942, B943, 
     & B944, B945, B946, B947, B948, B949, 
     & L921, L932, L943, L954, L965, L976, L987, L998, 
     & B91011, B91012, B91013, B91014, B91015, B91021, B91022, B91023, 
     & B91024, B91025, B91031, B91032, B91033, B91034, B91035, B91041, 
     & B91042, B91043, B91044, B91045 
       PARAMETER( B311 = 5d0/12d0, 
     &            B312 = 8d0/12d0, 
     &            B313 = -1d0/12d0) 
 
       PARAMETER( Cp31 = 1d0/24d0) 
ccccc L3 dimensione 4 
       PARAMETER( L321 =  4.012395124208693d-01, 
     &            L331 =  8.819910099032529d-04, 
     &            L341 =  1.728116022258560d-04, 
     &            L332 =  3.680857287181573d-01, 
     &            L342 =  1.635381132422046d-03, 
     &            L343 =  3.688541178419062d-01) 
c--------- B3-B5 OF DIMENSION  4 
       PARAMETER( B3511 = 49d0/720d0, 
     &            B3512 = -83d0/360d0, 
     &            B3513 = 17d0/60d0, 
     &            B3514 = -53d0/360d0, 
     &            B3515 = 19d0/720d0, 
     &            B3521 = 19d0/720d0, 
     &            B3522 = -23d0/360d0, 
     &            B3523 = 1d0/30d0, 
     & B3524 = 7d0/360d0, 
     & B3525 = -11d0/720d0, 
     & B3531 = -11d0/720d0, 
     & B3532 = 37d0/360d0, 
     & B3533 = -13d0/60d0, 
     & B3534 = 67d0/360d0, 
     & B3535 = -41d0/720d0, 
     & B3541 = 19d0/720d0, 
     & B3542 = -53d0/360d0, 
     & B3543 =  17d0/60d0, 
     & B3544 = -83d0/360d0, 
     & B3545 =  49d0/720d0) 
 
c[ 49/720, -83/360,  17/60, -53/360,  19/720] 
c[ 19/720, -23/360,   1/30,   7/360, -11/720] 
c[-11/720,  37/360, -13/60,  67/360, -41/720] 
c[ 19/720, -53/360,  17/60, -83/360,  49/720] 
 
C--------- A5, B5, B53, B56 Cp5 := MATRICES DEFINING  GAM5 
 
       PARAMETER( B511 = 251d0/720d0, 
     & B512 = 323d0/360d0, 
     & B513 = - 11d0/30d0, 
     & B514 = 53d0/360d0, 
     & B515 = -19d0/720d0, 
     & B521 = -19d0/720d0, 
     & B522 =  173d0/360d0, 
     & B523 = 19d0/30d0, 
     & B524 = -37d0/360d0, 
     & B525 = 11d0/720d0) 
 
       PARAMETER( B5711 = 1997d0/60480d0, 
     & B5712 = -113d0/630d0, 
     & B5713 = 1619d0/4032d0, 
     & B5714 = -715d0/1512d0, 
     & B5715 = 1241d0/4032d0, 
     & B5716 = -263d0/2520d0, 
     & B5717 = 863d0/60480d0, 
     & B5721 = -733d0/60480d0, 
     & B5722 = 41d0/630d0, 
     & B5723 = -193d0/1344d0, 
     & B5724 = 251d0/1512d0, 
     & B5725 = -425d0/4032d0, 
     & B5726 = 29d0/840d0, 
     & B5727 = -271d0/60480d0, 
     & B5731 = -271d0/60480d0, 
     & B5732 = 97d0/5040d0, 
     & B5733 = -13d0/448d0, 
     & B5734 = 5d0/378d0, 
     & B5735 = 37d0/4032d0, 
     & B5736 = -19d0/1680d0, 
     & B5737 = 191d0/60480d0, 
     & B5741 = 191d0/60480d0, 
     & B5742 = -67d0/2520d0, 
     & B5743 = 115d0/1344d0, 
     & B5744 = -211d0/1512d0, 
     & B5745 = 499d0/4032d0, 
     & B5746 = -2d0/35d0, 
     & B5747 = 653d0/60480d0) 
 
c--------- B5 - B7 
c 
c[1997/60480,  -113/630, 1619/4032, -715/1512, 1241/4032, -263/2520,  863/60480] 
c[-733/60480,    41/630, -193/1344,  251/1512, -425/4032,    29/840, -271/60480] 
c[-271/60480,   97/5040,   -13/448,     5/378,   37/4032,  -19/1680,  191/60480] 
c[ 191/60480,  -67/2520,  115/1344, -211/1512,  499/4032,     -2/35,  653/60480] 
c[-271/60480,    29/840, -425/4032,  251/1512, -193/1344,    41/630, -733/60480] 
c[ 863/60480, -263/2520, 1241/4032, -715/1512, 1619/4032,  -113/630, 1997/60480] 
c 
c 
       PARAMETER( L521 =  3.668340831928216D-01, 
     & L531 =  2.477905683677308D-03, 
     & L541 = -1.919925047010838D-03, 
     & L551 =  2.218385581234200D-03, 
     & L561 = -5.442189351609260D-03, 
     & L532 =  3.216639533696728D-01, 
     & L542 =  1.231925763308414D-03, 
     & L552 =  7.841944627374794D-03, 
     & L562 =  1.002485104590053D-03, 
     & L543 =  3.375100828961925D-01, 
     & L553 = -2.614300734741796D-04, 
     & L563 =  1.066631182323580D-03, 
     & L554 =  3.523137378783708D-01, 
     & L564 = -3.596681121610224D-04, 
     & L565 =  3.617716171655064D-01, 
     & CP51 =   3D0/160D0, 
     & CP52 = -11D0/1440D0) 
 
C--------- A7, B7, B57, B58 Cp7 := MATRICES DEFINING GAM7 
 
        PARAMETER( B711 = 19087d0/60480d0, 
     & B712 = 2713d0/2520d0, 
     & B713 = -15487d0/20160d0, 
     & B714 = 586d0/945d0, 
     & B715 = -6737d0/20160d0, 
     & B716 = 263d0/2520d0, 
     & B717 = -863d0/60480d0, 
     & B721 = -863d0/60480d0, 
     & B722 = 349d0/840d0, 
     & B723 = 5221d0/6720d0, 
     & B724 = -254d0/945d0, 
     & B725 = 811d0/6720d0, 
     & B726 = -29d0/840d0, 
     & B727 = 271d0/60480d0, 
     & B731 = 271d0/60480d0, 
     & B732 = -23d0/504d0, 
     & B733 = 10273d0/20160d0, 
     & B734 = 586d0/945d0, 
     & B735 = -2257d0/20160d0, 
     & B736 = 67d0/2520d0, 
     & B737 = -191d0/60480d0) 
 
C 
C--------- THE LAST THREE ROWS ARE THE REVERSE OF THE FIRST THREE 
C 
c B79 = 
c 
c[ 75203/3628800, -280187/1814400, 129781/259200, -238937/259200, 27289/25920, -197687/259200,  88531/259200, -156437/1814400,  33953/3628800] 
c[-17827/3628800,   66043/1814400, -30389/259200,   55513/259200, -6281/25920,   44983/259200, -19859/259200,   34453/1814400,  -7297/3628800] 
c[  8963/3628800,  -32987/1814400,  15061/259200,  -27257/259200,  3049/25920,  -21527/259200,   9331/259200,  -15797/1814400,   3233/3628800] 
c[  3233/3628800,  -10067/1814400,   3601/259200,   -4337/259200,     23/3240,    1393/259200,  -2129/259200,    7123/1814400,  -2497/3628800] 
c[ -2497/3628800,   12853/1814400,  -7859/259200,   18583/259200, -2681/25920,   24313/259200, -13589/259200,   30043/1814400,  -8227/3628800] 
c[  3233/3628800,  -15797/1814400,   9331/259200,  -21527/259200,  3049/25920,  -27257/259200,  15061/259200,  -32987/1814400,   8963/3628800] 
c[ -7297/3628800,   34453/1814400, -19859/259200,   44983/259200, -6281/25920,   55513/259200, -30389/259200,   66043/1814400, -17827/3628800] 
c[ 33953/3628800, -156437/1814400,  88531/259200, -197687/259200, 27289/25920, -238937/259200, 129781/259200, -280187/1814400,  75203/3628800] 
c 
       PARAMETER( B7911 =   75203d0/3628800d0, 
     & B7912 = -280187d0/1814400d0, 
     & B7913 =  129781d0/259200d0, 
     & B7914 = -238937d0/259200d0, 
     & B7915 =   27289d0/25920d0, 
     & B7916 = -197687d0/259200d0, 
     & B7917 =   88531d0/259200d0, 
     & B7918 = -156437d0/1814400d0, 
     & B7919 =   33953d0/3628800d0) 
 
       PARAMETER(B7921 = -17827d0/3628800d0, 
     & B7922 =  66043d0/1814400d0, 
     & B7923 = -30389d0/259200d0, 
     & B7924 =  55513d0/259200d0, 
     & B7925 =  -6281d0/25920d0, 
     & B7926 =  44983d0/259200d0, 
     & B7927 = -19859d0/259200d0, 
     & B7928 =  34453d0/1814400d0, 
     & B7929 =  -7297d0/3628800d0) 
 
       PARAMETER(B7931 =   8963d0/3628800d0, 
     & B7932 = -32987d0/1814400d0, 
     & B7933 =  15061d0/259200d0, 
     & B7934 = -27257d0/259200d0, 
     & B7935 =   3049d0/25920d0, 
     & B7936 = -21527d0/259200d0, 
     & B7937 =   9331d0/259200d0, 
     & B7938 = -15797d0/1814400d0, 
     & B7939 =   3233d0/3628800d0) 
 
       PARAMETER(B7941 =   3233d0/3628800d0, 
     & B7942 = -10067d0/1814400d0, 
     & B7943 =   3601d0/259200d0, 
     & B7944 =  -4337d0/259200d0, 
     & B7945 =     23d0/3240d0, 
     & B7946 =   1393d0/259200d0, 
     & B7947 =  -2129d0/259200d0, 
     & B7948 =   7123d0/1814400d0, 
     & B7949 =  -2497d0/3628800d0) 
 
 
       PARAMETER(B7951 =  -2497d0/3628800d0, 
     & B7952 =  12853d0/1814400d0, 
     & B7953 =  -7859d0/259200d0, 
     & B7954 =  18583d0/259200d0, 
     & B7955 =  -2681d0/25920d0, 
     & B7956 =  24313d0/259200d0, 
     & B7957 = -13589d0/259200d0, 
     & B7958 =  30043d0/1814400d0, 
     & B7959 =  -8227d0/3628800d0) 
 
       PARAMETER(L721 = 3.023839891568610D-01, 
     & L731 = 3.201698610574002D-05, 
     & L741 = 4.193101163680004D-04, 
     & L751 = 1.686924996069667D-04, 
     & L761 = 4.806043527549464D-05, 
     & L771 = 3.598347048026785D-06, 
     & L781 = 7.892534649789167D-04, 
     & L732 = 2.559868364091398D-01, 
     & L742 = 1.336896192287030D-04, 
     & L752 = 3.080994719931695D-03, 
     & L762 = 1.457177183563680D-04, 
     & L772 = 9.259360509484074D-04, 
     & L782 = 2.397658879381223D-04, 
     & L743 = 2.639734712170458D-01, 
     & L753 = 1.734338929611258D-04, 
     & L763 = 6.704398263264620D-03, 
     & L773 = 4.559927214651730D-05, 
     & L783 = 6.396418554053151D-05, 
     & L754 = 2.817729090368562D-01, 
     & L764 = 2.877761776030408D-04, 
     & L774 = 1.810919475521773D-04, 
     & L784 = 1.009049833235848D-03, 
     & L765 = 2.993040718034231D-01, 
     & L775 = 2.009850887505898D-03, 
     & L785 = 1.748065618845750D-03, 
     & L776 = 3.150349043479135D-01, 
     & L786 = 3.243816792609449D-05, 
     & L787 = 3.271307059448932D-01) 
 
       PARAMETER(CP71 = 103D0/9061D0, 
     & CP72 = -13D0/4480D0, 
     & CP73 =  67D0/42431D0) 
 
C--------- A8, B8, B86, B810 Cp8 := MATRICES DEFINING GAM9 
 
       PARAMETER(B911 = 1070017D0/3628800D0, 
     & B912 = 2233547D0/1814400D0, 
     & B913 = -2302297D0/1814400D0, 
     & B914 = 2797679D0/1814400D0, 
     & B915 = -31457D0/22680D0, 
     & B916 = 1573169D0/1814400D0, 
     & B917 = -645607D0/1814400D0, 
     & B918 = 156437D0/1814400D0, 
     & B919 = -33953D0/3628800D0, 
     & B921 = -33953D0/3628800D0, 
     & B922 = 687797D0/1814400D0, 
     & B923 =  1622393D0/1814400D0, 
     & B924 = -876271D0/1814400D0, 
     & B925 =   8233D0/22680D0, 
     & B926 =    -377521D0/1814400D0, 
     & B927 =   147143D0/1814400D0, 
     & B928 =  -34453D0/1814400D0, 
     & B929 =   7297D0/3628800D0, 
     & B931 = 7297D0/3628800D0, 
     & B932 =  -49813D0/1814400D0, 
     & B933 =  819143D0/1814400D0, 
     & B934 =  1315919D0/1814400D0, 
     & B935 = -5207D0/22680D0, 
     & B936 =  198929D0/1814400D0, 
     & B937 =  -71047D0/1814400D0, 
     & B938 =  15797D0/1814400D0, 
     & B939 = -3233D0/3628800D0, 
     & B941 = -3233D0/3628800D0, 
     & B942 = 18197D0/1814400D0, 
     & B943 =  -108007D0/1814400D0, 
     & B944 =  954929D0/1814400D0, 
     & B945 = 13903D0/22680D0, 
     & B946 = -212881D0/1814400D0, 
     & B947 =  63143D0/1814400D0, 
     & B948 =  -12853D0/1814400D0, 
     & B949 =  2497D0/3628800D0) 
 
 
c   B910=B9-B10; 
c               prime 4 righe : la 5 e' uguale alla 4; 
c     le ultime 4 uguali alle prime 4 negate (simmetriche); 
c     le prime  5 colonne  : le altre sono simmetriche e negate; 
c 
c[ 8183/1036800, -8183/115200,   8183/28800, -57281/86400,  57281/57600, -57281/57600,  57281/86400,  -8183/28800,  8183/115200, -8183/1036800] 
c[  -425/290304,    425/32256,    -425/8064,     425/3456,    -425/2304,     425/2304,    -425/3456,     425/8064,   -425/32256,    425/290304] 
c[      7/12800,    -63/12800,      63/3200,    -147/3200,     441/6400,    -441/6400,     147/3200,     -63/3200,     63/12800,      -7/12800] 
c[-2497/7257600,  2497/806400, -2497/201600,   2497/86400,  -2497/57600,   2497/57600,  -2497/86400,  2497/201600, -2497/806400,  2497/7257600] 
c[-2497/7257600,  2497/806400, -2497/201600,   2497/86400,  -2497/57600,   2497/57600,  -2497/86400,  2497/201600, -2497/806400,  2497/7257600] 
c[ 2497/7257600, -2497/806400,  2497/201600,  -2497/86400,   2497/57600,  -2497/57600,   2497/86400, -2497/201600,  2497/806400, -2497/7257600] 
c[     -7/12800,     63/12800,     -63/3200,     147/3200,    -441/6400,     441/6400,    -147/3200,      63/3200,    -63/12800,       7/12800] 
c[   425/290304,   -425/32256,     425/8064,    -425/3456,     425/2304,    -425/2304,     425/3456,    -425/8064,    425/32256,   -425/290304] 
c[-8183/1036800,  8183/115200,  -8183/28800,  57281/86400, -57281/57600,  57281/57600, -57281/86400,   8183/28800, -8183/115200,  8183/1036800] 
 
       PARAMETER(B91011 =   8183d0/1036800d0, 
     & B91012 =  -8183d0/115200d0, 
     & B91013 =   8183d0/28800d0, 
     & B91014 = -57281d0/86400d0, 
     & B91015 =  57281d0/57600d0) 
 
       PARAMETER(B91021 =   -425d0/290304d0, 
     & B91022 =    425d0/32256d0, 
     & B91023 =   -425d0/8064d0, 
     & B91024 =    425d0/3456d0, 
     & B91025 =   -425d0/2304d0) 
 
       PARAMETER(B91031 =     7d0/12800d0, 
     & B91032 =   -63d0/12800d0, 
     & B91033 =    63d0/3200d0, 
     & B91034 =  -147d0/3200d0, 
     & B91035 =   441d0/6400d0) 
 
       PARAMETER(B91041 = -2497d0/7257600d0, 
     & B91042 =  2497d0/806400d0, 
     & B91043 = -2497d0/201600d0, 
     & B91044 =  2497d0/86400d0, 
     & B91045 = -2497d0/57600d0) 
 
       PARAMETER(Cp91 =  7.892554012345216d-03, 
     & Cp92 = -1.463982583774219d-03, 
     & Cp93 =  5.468749999999983d-04, 
     & Cp94 = -3.440531305114634d-04) 
 
C--------- THE OTHERS ARE THE SAME WITH CHANGED SIGN 
 
       PARAMETER( L921 = 2.590721934790442d-01, 
     & L932   =   2.077575545359853d-01, 
     & L943   =   2.032874698558627d-01, 
     & L954   =   2.036384888660128d-01, 
     & L965   =   2.039599505779785d-01, 
     & L976   =   2.034044409161703d-01, 
     & L987   =   2.017245408702437d-01, 
     & L998   =   1.986549276295617d-01) 
 
C
C   INPUT VARIABLES
C------------------------------------
      INTEGER R, ORD,ORDMIN, ORDMAX, DBLK, IJOB, IPIV(R), LDLU
      DOUBLE PRECISION  H, SCAL(R),  LU(LDLU,R)  
 
C
C   INPUT/OUTPUT VARIABLES  ksks:,1->,*
C------------------------------------
      DOUBLE PRECISION  ERRV(DBLK), ERRSAME, ERRUP, ERRDOWN,
     &                  FP(R,*), F(R,*), DN(R,*)
C
C   LOCAL VARIABLES
C------------------------------------
      INTEGER  I, J, ORDDOWN
      DOUBLE PRECISION ERRVJ,          
     &                 FP1,FP2,FP3,FP4,FP5,FP6,FP7,FP8,FP9,FP10
C
C   EXECUTABLE STATEMENTS
C---------------------------------
C          
C--------- ERRSAME ESTIMATION
       GOTO (10,20,30,40),  ORD
 10   CONTINUE
       DO I=1,R
        FP1= FP(I,1)
        FP2= FP(I,2)
        FP3= FP(I,3)
        FP4= FP(I,4)
        FP5= FP(I,5)
        F(I,1) = H*(B3511*FP1+B3512*FP2+B3513*FP3+B3514*FP4+B3515*FP5)
        F(I,2) = H*(B3521*FP1+B3522*FP2+B3523*FP3+B3524*FP4+B3525*FP5)
        F(I,3) = H*(B3531*FP1+B3532*FP2+B3533*FP3+B3534*FP4+B3535*FP5)
        F(I,4) = H*(B3541*FP1+B3542*FP2+B3543*FP3+B3544*FP4+B3545*FP5)
       END DO

      GOTO 50
 20   CONTINUE
        DO J=1,R
        FP1= FP(J,1)
        FP2= FP(J,2)
        FP3= FP(J,3)
        FP4= FP(J,4)
        FP5= FP(J,5)
        FP6= FP(J,6)
        FP7= FP(J,7)
       F(J,1) = H*(B5711*FP1+B5712*FP2+B5713*FP3
     &            +B5714*FP4+B5715*FP5+B5716*FP6+B5717*FP7)
       F(J,2) = H*(B5721*FP1+B5722*FP2+B5723*FP3
     &            +B5724*FP4+B5725*FP5+B5726*FP6+B5727*FP7)
       F(J,3) = H*(B5731*FP1+B5732*FP2+B5733*FP3
     &            +B5734*FP4+B5735*FP5+B5736*FP6+B5737*FP7)
       F(J,4) = H*(B5741*FP1+B5742*FP2+B5743*FP3
     &            +B5744*FP4+B5745*FP5+B5746*FP6+B5747*FP7)
       F(J,5) = H*(B5727*FP1+B5726*FP2+B5725*FP3
     &            +B5724*FP4+B5723*FP5+B5722*FP6+B5721*FP7)
       F(J,6) = H*(B5717*FP1+B5716*FP2+B5715*FP3
     &            +B5714*FP4+B5713*FP5+B5712*FP6+B5711*FP7)
        END DO

      GOTO 50
  30  CONTINUE
        DO J=1,R
          FP1= FP(J,1)
          FP2= FP(J,2)
          FP3= FP(J,3)
          FP4= FP(J,4)
          FP5= FP(J,5)
          FP6= FP(J,6)
          FP7= FP(J,7)
          FP8= FP(J,8)
          FP9= FP(J,9)
       F(J,1) = H*(B7911*FP1+B7912*FP2+B7913*FP3+B7914*FP4
     & +B7915*FP5+B7916*FP6+B7917*FP7+B7918*FP8+B7919*FP9)
       F(J,2) = H*(B7921*FP1+B7922*FP2+B7923*FP3+B7924*FP4
     & +B7925*FP5+B7926*FP6+B7927*FP7+B7928*FP8+B7929*FP9)
       F(J,3) = H*(B7931*FP1+B7932*FP2+B7933*FP3+B7934*FP4
     & +B7935*FP5+B7936*FP6+B7937*FP7+B7938*FP8+B7939*FP9)
       F(J,4) = H*(B7941*FP1+B7942*FP2+B7943*FP3+B7944*FP4
     & +B7945*FP5+B7946*FP6+B7947*FP7+B7948*FP8+B7949*FP9)
       F(J,5) = H*(B7951*FP1+B7952*FP2+B7953*FP3+B7954*FP4
     & +B7955*FP5+B7956*FP6+B7957*FP7+B7958*FP8+B7959*FP9)
       F(J,6) = H*(B7939*FP1+B7938*FP2+B7937*FP3+B7936*FP4
     & +B7935*FP5+B7934*FP6+B7933*FP7+B7932*FP8+B7931*FP9)
       F(J,7) = H*(B7929*FP1+B7928*FP2+B7927*FP3+B7926*FP4
     & +B7925*FP5+B7924*FP6+B7923*FP7+B7922*FP8+B7921*FP9)
       F(J,8) = H*(B7919*FP1+B7918*FP2+B7917*FP3+B7916*FP4
     & +B7915*FP5+B7914*FP6+B7913*FP7+B7912*FP8+B7911*FP9)
       END DO

      GOTO 50
  40  CONTINUE
        DO J=1,R
          FP1= FP(J,1)
          FP2= FP(J,2)
          FP3= FP(J,3)
          FP4= FP(J,4)
          FP5= FP(J,5)
          FP6= FP(J,6)
          FP7= FP(J,7)
          FP8= FP(J,8)
          FP9= FP(J,9)
          FP10= FP(J,10)
       F(J,1) = H*(B91011*(FP1-FP10)+B91012*(FP2-FP9)+B91013*(FP3-FP8)
     & +B91014*(FP4-FP7)+B91015*(FP5-FP6) )
       F(J,2) = H*(B91021*(FP1-FP10)+B91022*(FP2-FP9)+B91023*(FP3-FP8)
     & +B91024*(FP4-FP7)+B91025*(FP5-FP6) )
       F(J,3) = H*(B91031*(FP1-FP10)+B91032*(FP2-FP9)+B91033*(FP3-FP8)
     & +B91034*(FP4-FP7)+B91035*(FP5-FP6) )
       F(J,4) = H*(B91041*(FP1-FP10)+B91042*(FP2-FP9)+B91043*(FP3-FP8)
     & +B91044*(FP4-FP7)+B91045*(FP5-FP6) )
       F(J,5) =  F(J,4)
       F(J,6) = -F(J,4)
       F(J,7) = -F(J,3)
       F(J,8) = -F(J,2)
       F(J,9) = -F(J,1)
       END DO

  50  CONTINUE

C--------- A SINGLE SPLITTING-NEWTON ITERATION FOR ERRSAME
           CALL NEWTGS(R,DBLK,LU,LDLU,IPIV,F,DN,IJOB)


C--------- COMPUTE  ERRSAME AND ERRV (VECTOR ERROR)
              ERRSAME = 0D0
              DO J=1,DBLK
                ERRV(J) = 0D0
                DO I=1,R
                  FP1 = (DN(I,J)/SCAL(I) )
                  ERRV(J) =  ERRV(J)+ FP1*FP1
                END DO
                ERRV(J) = DSQRT(ERRV(J)/R)
                ERRSAME = DMAX1( ERRSAME, ERRV(J) )
              END DO
              ERRSAME = DMAX1(ERRSAME, 1d-15)
       ERRDOWN = 0D0
       ERRUP   = 0D0
       IF ( ERRSAME .LE. 1d0) THEN

       IF (ORD .LT. ORDMAX) THEN
C--------- ERRUP ESTIMATION
         GOTO (11,21,31), ORD
 11    CONTINUE
          DO I=1,R
             FP1 =  F(I,1)/CP31
             FP2 =  F(I,2)/CP31
             FP3 =  F(I,3)/CP31
             FP4 = -F(I,4)/CP31
             F(I,1) =  (-FP1 + 2d0*FP2 - FP3)*CP51
             F(I,2) =  (-FP1 + 2d0*FP2 - FP3)*CP52
             F(I,3) = -(-FP2 + 2d0*FP3 - FP4)*CP52
             F(I,4) = -(-FP2 + 2d0*FP3 - FP4)*CP51
          END DO

          GOTO 41
 21    CONTINUE
           DO I = 1, R
             FP1 =  F(I,1)/CP51
             FP2 =  F(I,2)/CP52
             FP3 =  F(I,3)/CP52
             FP4 =  F(I,4)/CP52
             FP5 = -F(I,5)/CP52
             FP6 = -F(I,6)/CP51
             F(I,1) =  (-FP1 + 2d0*FP2 - FP3)*CP71
             F(I,2) =  (-FP1 + 2d0*FP2 - FP3)*CP72
             F(I,3) =  (-FP2 + 2d0*FP3 - FP4)*CP73
             F(I,4) = -(-FP3 + 2d0*FP4 - FP5)*CP73
             F(I,5) = -(-FP4 + 2d0*FP5 - FP6)*CP72
             F(I,6) = -(-FP4 + 2d0*FP5 - FP6)*CP71

          END DO
          GOTO 41
  31   CONTINUE
           DO I = 1, R
             FP1 =  F(I,1)/CP71
             FP2 =  F(I,2)/CP72
             FP3 =  F(I,3)/CP73
             FP4 =  F(I,4)/CP73
             FP5 =  F(I,5)/CP73
             FP6 = -F(I,6)/CP73
             FP7 = -F(I,7)/CP72
             FP8 = -F(I,8)/CP71
             F(I,1) =  (-FP1 + 2d0*FP2 - FP3)*CP91
             F(I,2) =  (-FP1 + 2d0*FP2 - FP3)*CP92
             F(I,3) =  (-FP2 + 2d0*FP3 - FP4)*CP93
             F(I,4) =  (-FP3 + 2d0*FP4 - FP5)*CP94
             F(I,5) = -(-FP4 + 2d0*FP5 - FP6)*CP94
             F(I,6) = -(-FP5 + 2d0*FP6 - FP7)*CP93
             F(I,7) = -(-FP6 + 2d0*FP7 - FP8)*CP92
             F(I,8) = -(-FP6 + 2d0*FP7 - FP8)*CP91
           END DO

  41   CONTINUE

C--------- A SINGLE SPLITTING-NEWTON ITERATION FOR ERRUP

          CALL NEWTGS(R,DBLK,LU,LDLU,IPIV,F,DN,IJOB)


C-------- COMPUTE ERRUP
            ERRUP = 0D0
              DO J=1,DBLK
                ERRVJ = 0D0
                DO I=1,R
                  FP1 = (DN(I,J)/SCAL(I) )
                  ERRVJ =  ERRVJ + FP1*FP1
                END DO
                 ERRVJ = DSQRT(ERRVJ/R)
                 ERRUP = DMAX1( ERRUP, ERRVJ )
              END DO
              ERRUP = dmax1(ERRUP, 1d-15)
        END IF
        IF (ORD .GT. ORDMIN) THEN
C--------- ERRDOWN ESTIMATION
          ORDDOWN = ORD-1
          GOTO (13, 23, 33), ORDDOWN
 13   CONTINUE
       DO J=1,R
         FP1= FP(J,1)
         FP2= FP(J,2)
         FP3= FP(J,3)
         FP4= FP(J,4)
         FP5= FP(J,5)
         FP6= FP(J,6)
         FP7= FP(J,7)
       F(J,1) = -H*(B3511*FP1+B3512*FP2+B3513*FP3+B3514*FP4+B3515*FP5)
       F(J,2) = -H*(B3521*FP1+B3522*FP2+B3523*FP3+B3524*FP4+B3525*FP5)
       F(J,3) = -H*(B3521*FP2+B3522*FP3+B3523*FP4+B3524*FP5+B3525*FP6)
       F(J,4) = -H*(B3521*FP3+B3522*FP4+B3523*FP5+B3524*FP6+B3525*FP7)
       F(J,5) = -H*(B3531*FP3+B3532*FP4+B3533*FP5+B3534*FP6+B3535*FP7)
       F(J,6) = -H*(B3541*FP3+B3542*FP4+B3543*FP5+B3544*FP6+B3545*FP7)
       END DO

      GOTO 43
 23   CONTINUE
       DO J=1,R
          FP1= FP(J,1)
          FP2= FP(J,2)
          FP3= FP(J,3)
          FP4= FP(J,4)
          FP5= FP(J,5)
          FP6= FP(J,6)
          FP7= FP(J,7)
          FP8= FP(J,8)
          FP9= FP(J,9)
       F(J,1) = -H*(B5711*FP1+B5712*FP2+B5713*FP3
     &            +B5714*FP4+B5715*FP5+B5716*FP6+B5717*FP7)
       F(J,2) = -H*(B5721*FP1+B5722*FP2+B5723*FP3
     &            +B5724*FP4+B5725*FP5+B5726*FP6+B5727*FP7)
       F(J,3) = -H*(B5731*FP1+B5732*FP2+B5733*FP3
     &            +B5734*FP4+B5735*FP5+B5736*FP6+B5737*FP7)
       F(J,4) = -H*(B5731*FP2+B5732*FP3+B5733*FP4
     &            +B5734*FP5+B5735*FP6+B5736*FP7+B5737*FP8)
       F(J,5) = -H*(B5731*FP3+B5732*FP4+B5733*FP5
     &            +B5734*FP6+B5735*FP7+B5736*FP8+B5737*FP9)
       F(J,6) = -H*(B5741*FP3+B5742*FP4+B5743*FP5
     &            +B5744*FP6+B5745*FP7+B5746*FP8+B5747*FP9)
       F(J,7) = -H*(B5727*FP3+B5726*FP4+B5725*FP5
     &            +B5724*FP6+B5723*FP7+B5722*FP8+B5721*FP9)
       F(J,8) = -H*(B5717*FP3+B5716*FP4+B5715*FP5
     &            +B5714*FP6+B5713*FP7+B5712*FP8+B5711*FP9)

        END DO

      GOTO 43
  33  CONTINUE
       DO J=1,R
          FP1= FP(J,1)
          FP2= FP(J,2)
          FP3= FP(J,3)
          FP4= FP(J,4)
          FP5= FP(J,5)
          FP6= FP(J,6)
          FP7= FP(J,7)
          FP8= FP(J,8)
          FP9= FP(J,9)
          FP10= FP(J,10)
       F(J,1) = H*(B7911*FP1+B7912*FP2+B7913*FP3+B7914*FP4
     & +B7915*FP5+B7916*FP6+B7917*FP7+B7918*FP8+B7919*FP9)

       F(J,2) = H*(B7921*FP1+B7922*FP2+B7923*FP3+B7924*FP4
     & +B7925*FP5+B7926*FP6+B7927*FP7+B7928*FP8+B7929*FP9)

       F(J,3) = H*(B7931*FP1+B7932*FP2+B7933*FP3+B7934*FP4
     & +B7935*FP5+B7936*FP6+B7937*FP7+B7938*FP8+B7939*FP9)

       F(J,4) = H*(B7941*FP1+B7942*FP2+B7943*FP3+B7944*FP4
     & +B7945*FP5+B7946*FP6+B7947*FP7+B7948*FP8+B7949*FP9)

       F(J,5) = -H*(B7941*FP2+B7942*FP3+B7943*FP4+B7944*FP5
     & +B7945*FP6+B7946*FP7+B7947*FP8+B7948*FP9+B7949*FP10)


       F(J,6) = -H*(B7951*FP2+B7952*FP3+B7953*FP4+B7954*FP5
     & +B7955*FP6+B7956*FP7+B7957*FP8+B7958*FP9+B7959*FP10)

       F(J,7) = -H*(B7939*FP2+B7938*FP3+B7937*FP4+B7936*FP5
     & +B7935*FP6+B7934*FP7+B7933*FP8+B7932*FP9+B7931*FP10)

       F(J,8) = -H*(B7929*FP2+B7928*FP3+B7927*FP4+B7926*FP5
     & +B7925*FP6+B7924*FP7+B7923*FP8+B7922*FP9+B7921*FP10)

       F(J,9) = -H*(B7919*FP2+B7918*FP3+B7917*FP4+B7916*FP5
     & +B7915*FP6+B7914*FP7+B7913*FP8+B7912*FP9+B7911*FP10)

            END DO

  43  CONTINUE

C--------- A SINGLE SPLITTING-NEWTON ITERATION FOR ERRDOWN
           CALL NEWTGS(R,DBLK,LU,LDLU,IPIV,F,DN,IJOB)


C--------- COMPUTE ERRDOWN
            ERRDOWN = 0D0
              DO J=1,DBLK
                ERRVJ = 0D0
                DO I=1,R
                  FP1 = (DN(I,J)/SCAL(I) )
                  ERRVJ =  ERRVJ + FP1*FP1
                END DO
                 ERRVJ = DSQRT(ERRVJ/R)
                 ERRDOWN = DMAX1( ERRDOWN, ERRVJ )
              END DO
              ERRDOWN = dmax1(ERRDOWN, 1d-15)
        END IF
      END IF
      RETURN
      END
C
C     SUBROUTINE DEC_GAM
C
      SUBROUTINE DEC_GAM (N, NDIM, A, IP, IER)
C VERSION REAL DOUBLE PRECISION
      INTEGER N,NDIM,IP,IER,NM1,K,KP1,M,I,J
      DOUBLE PRECISION A,T
      DIMENSION A(NDIM,N), IP(N)
C-----------------------------------------------------------------------
C  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION.
C  INPUT..
C     N = ORDER OF MATRIX.
C     NDIM = DECLARED DIMENSION OF ARRAY  A .
C     A = MATRIX TO BE TRIANGULARIZED.
C  OUTPUT..
C     A(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U .
C     A(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
C     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW.
C     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O .
C     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE
C           SINGULAR AT STAGE K.
C  USE  SOL_GAM  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
C  DETERM(A) = IP(N)*A(1,1)*A(2,2)*...*A(N,N).
C  IF IP(N)=O, A IS SINGULAR, SOL_GAM WILL DIVIDE BY ZERO.
C
C  REFERENCE..
C     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
C     C.A.C.M. 15 (1972), P. 274.
C-----------------------------------------------------------------------
      IER = 0
      IP(N) = 1
      IF (N .EQ. 1) GO TO 70
      NM1 = N - 1
      DO 60 K = 1,NM1
        KP1 = K + 1
        M = K
        DO 10 I = KP1,N
          IF (DABS(A(I,K)) .GT. DABS(A(M,K))) M = I
 10     CONTINUE
        IP(K) = M
        T = A(M,K)
        IF (M .EQ. K) GO TO 20
        IP(N) = -IP(N)
        A(M,K) = A(K,K)
        A(K,K) = T
 20     CONTINUE
        IF (T .EQ. 0.D0) GO TO 80
        T = 1.D0/T
        DO 30 I = KP1,N
 30       A(I,K) = -A(I,K)*T
        DO 50 J = KP1,N
          T = A(M,J)
          A(M,J) = A(K,J)
          A(K,J) = T
          IF (T .EQ. 0.D0) GO TO 45
          DO 40 I = KP1,N
 40         A(I,J) = A(I,J) + A(I,K)*T
 45       CONTINUE
 50       CONTINUE
 60     CONTINUE
 70   K = N
      IF (A(N,N) .EQ. 0.D0) GO TO 80
      RETURN
 80   IER = K
      IP(N) = 0
      RETURN
C----------------------- END OF SUBROUTINE DEC_GAM -------------------------
      END
C
C     SUBROUTINE SOL_GAM
C
      SUBROUTINE SOL_GAM (N, NDIM, A, B, IP)
C VERSION REAL DOUBLE PRECISION
      INTEGER N,NDIM,IP,NM1,K,KP1,M,I,KB,KM1
      DOUBLE PRECISION A,B,T
      DIMENSION A(NDIM,N), B(N), IP(N)
C-----------------------------------------------------------------------
C  SOLUTION OF LINEAR SYSTEM, A*X = B .
C  INPUT..
C    N = ORDER OF MATRIX.
C    NDIM = DECLARED DIMENSION OF ARRAY  A .
C    A = TRIANGULARIZED MATRIX OBTAINED FROM DEC_GAM.
C    B = RIGHT HAND SIDE VECTOR.
C    IP = PIVOT VECTOR OBTAINED FROM DEC_GAM.
C  DO NOT USE IF DEC_GAM HAS SET IER .NE. 0.
C  OUTPUT..
C    B = SOLUTION VECTOR, X .
C-----------------------------------------------------------------------
      IF (N .EQ. 1) GO TO 50
      NM1 = N - 1
      DO 20 K = 1,NM1
        KP1 = K + 1
        M = IP(K)
        T = B(M)
        B(M) = B(K)
        B(K) = T
        DO 10 I = KP1,N
 10       B(I) = B(I) + A(I,K)*T
 20     CONTINUE
      DO 40 KB = 1,NM1
        KM1 = N - KB
        K = KM1 + 1
        B(K) = B(K)/A(K,K)
        T = -B(K)
        DO 30 I = 1,KM1
 30       B(I) = B(I) + A(I,K)*T
 40     CONTINUE
 50   B(1) = B(1)/A(1,1)
      RETURN
C----------------------- END OF SUBROUTINE SOL_GAM -------------------------
      END
C
C     SUBROUTINE DECB_gam
C
      SUBROUTINE DECB_gam (N, NDIM, A, ML, MU, IP, IER)
      REAL*8 A,T
      DIMENSION A(NDIM,N), IP(N)
C-----------------------------------------------------------------------
C  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A BANDED
C  MATRIX WITH LOWER BANDWIDTH ML AND UPPER BANDWIDTH MU
C  INPUT..
C     N       ORDER OF THE ORIGINAL MATRIX A.
C     NDIM    DECLARED DIMENSION OF ARRAY  A.
C     A       CONTAINS THE MATRIX IN BAND STORAGE.   THE COLUMNS
C                OF THE MATRIX ARE STORED IN THE COLUMNS OF  A  AND
C                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS
C                ML+1 THROUGH 2*ML+MU+1 OF  A.
C     ML      LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C     MU      UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C  OUTPUT..
C     A       AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND
C                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
C     IP      INDEX VECTOR OF PIVOT INDICES.
C     IP(N)   (-1)**(NUMBER OF INTERCHANGES) OR O .
C     IER     = 0 IF MATRIX A IS NONSINGULAR, OR  = K IF FOUND TO BE
C                SINGULAR AT STAGE K.
C  USE  SOLB_gam  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
C  DETERM(A) = IP(N)*A(MD,1)*A(MD,2)*...*A(MD,N)  WITH MD=ML+MU+1.
C  IF IP(N)=O, A IS SINGULAR, SOLB_gam WILL DIVIDE BY ZERO.
C
C  REFERENCE..
C     THIS IS A MODIFICATION OF
C     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
C     C.A.C.M. 15 (1972), P. 274.
C-----------------------------------------------------------------------
      IER = 0
      IP(N) = 1
      MD = ML + MU + 1
      MD1 = MD + 1
      JU = 0
      IF (ML .EQ. 0) GO TO 70
      IF (N .EQ. 1) GO TO 70
      IF (N .LT. MU+2) GO TO 7
      DO 5 J = MU+2,N
      DO 5 I = 1,ML
  5   A(I,J) = 0.D0
  7   NM1 = N - 1
      DO 60 K = 1,NM1
        KP1 = K + 1
        M = MD
        MDL = MIN(ML,N-K) + MD
        DO 10 I = MD1,MDL
          IF (DABS(A(I,K)) .GT. DABS(A(M,K))) M = I
 10     CONTINUE
        IP(K) = M + K - MD
        T = A(M,K)
        IF (M .EQ. MD) GO TO 20
        IP(N) = -IP(N)
        A(M,K) = A(MD,K)
        A(MD,K) = T
 20     CONTINUE
        IF (T .EQ. 0.D0) GO TO 80
        T = 1.D0/T
        DO 30 I = MD1,MDL
 30       A(I,K) = -A(I,K)*T
        JU = MIN0(MAX0(JU,MU+IP(K)),N)
        MM = MD
        IF (JU .LT. KP1) GO TO 55
        DO 50 J = KP1,JU
          M = M - 1
          MM = MM - 1
          T = A(M,J)
          IF (M .EQ. MM) GO TO 35
          A(M,J) = A(MM,J)
          A(MM,J) = T
 35       CONTINUE
          IF (T .EQ. 0.D0) GO TO 45
          JK = J - K
          DO 40 I = MD1,MDL
            IJK = I - JK
 40         A(IJK,J) = A(IJK,J) + A(I,K)*T
 45       CONTINUE
 50       CONTINUE
 55     CONTINUE
 60     CONTINUE
 70   K = N
      IF (A(MD,N) .EQ. 0.D0) GO TO 80
      RETURN
 80   IER = K
      IP(N) = 0
      RETURN
C----------------------- END OF SUBROUTINE DECB_gam ------------------------
      END
C
C     SUBROUTINE SOLB_gam
C
      SUBROUTINE SOLB_gam (N, NDIM, A, ML, MU, B, IP)
      REAL*8 A,B,T
      DIMENSION A(NDIM,N), B(N), IP(N)
C-----------------------------------------------------------------------
C  SOLUTION OF LINEAR SYSTEM, A*X = B .
C  INPUT..
C    N      ORDER OF MATRIX A.
C    NDIM   DECLARED DIMENSION OF ARRAY  A .
C    A      TRIANGULARIZED MATRIX OBTAINED FROM DECB_gam.
C    ML     LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C    MU     UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C    B      RIGHT HAND SIDE VECTOR.
C    IP     PIVOT VECTOR OBTAINED FROM DECB_gam.
C  DO NOT USE IF DECB_gam HAS SET IER .NE. 0.
C  OUTPUT..
C    B      SOLUTION VECTOR, X .
C-----------------------------------------------------------------------
      MD = ML + MU + 1
      MD1 = MD + 1
      MDM = MD - 1
      NM1 = N - 1
      IF (ML .EQ. 0) GO TO 25
      IF (N .EQ. 1) GO TO 50
      DO 20 K = 1,NM1
        M = IP(K)
        T = B(M)
        B(M) = B(K)
        B(K) = T
        MDL = MIN(ML,N-K) + MD
        DO 10 I = MD1,MDL
          IMD = I + K - MD
 10       B(IMD) = B(IMD) + A(I,K)*T
 20     CONTINUE
 25   CONTINUE
      DO 40 KB = 1,NM1
        K = N + 1 - KB
        B(K) = B(K)/A(MD,K)
        T = -B(K)
        KMD = MD - K
        LM = MAX0(1,KMD+1)
        DO 30 I = LM,MDM
          IMD = I - KMD
 30       B(IMD) = B(IMD) + A(I,K)*T
 40     CONTINUE
 50   B(1) = B(1)/A(MD,1)
      RETURN
C----------------------- END OF SUBROUTINE SOLB_gam ------------------------
      END

