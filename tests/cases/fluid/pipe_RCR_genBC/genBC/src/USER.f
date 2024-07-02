c     MUSC 2 Immediate postop
c     Ethan Kung  keo@ucsd.edu

c     Created by Mahdi Esmaily Moghadam 12-01-2010
c     Please report any problem to mesmaily@ucsd.edu, memt63@yahoo.com

c This subroutine initializes the parameters, you need to read the
c comments and specify them carefully
c--------------------------------------------------------------------
c This is an example for RCR boundary condition with parameters
c Rd, R, and C which are distal, and proximal resistance and capacitor.

      SUBROUTINE INITIALIZE(nTimeStep)
      USE COM
      IMPLICIT NONE
      INTENT(OUT) nTimeStep

      LOGICAl ierr
      INTEGER i, nTimeStep
      REAL(KIND=8), ALLOCATABLE :: tZeroX(:)
c
c********************************************************************
c For instance if pressure in 3D solver is in cgs and here mmHg
c pConv=1334=1.334D3, also if flowrate in 3D solver is mm^3/s and
c here is mL/s qConv=1000=1D3. In the case both solver are using same
c unites you can set these two conversion coefficients to 1D0
      pConv = 1D0
      qConv = 1D0

c Only when all the surfaces of you model are coupled with NeumannSrfs
c you may set this to .TRUE.
      pCorr = .FALSE.
      qCorr = .FALSE.

c********************************************************************
c Block for your inputs

c These two value should match to that of defined in solver.inp
      nDirichletSrfs = 0
      nNeumannSrfs   = 1
c Number of unknowns that you need inside your lumped parameter network
      nUnknowns      = 1
c Number of time step between N and N+alpha, for higher values this code
c would be more accurate and more costly
      nTimeStep = 1000

c Number of parameters to be printed in AllData file (the first
c nUnknowns columns are the "X" values)
      nXprint = 3

c--------------------------------------------------------------------
c You don't need to change this part                                 !
      ALLOCATE (tZeroX(nUnknowns), srfToXdPtr(nDirichletSrfs))       !
      ALLOCATE (srfToXPtr(nNeumannSrfs))                             !
      tZeroX = 0D0                                                   !
c--------------------------------------------------------------------

c Value of your unknown at time equal to zero (This is going to be used 
c ONLY when you are initiating simulation)
       
C INITIALIZE INLET FLOW AND RCR CAPACITOR PRESSURE
      tZeroX(1) = 0D0

c--------------------------------------------------------------------
c You don't need to change this part                                 !
      INQUIRE (FILE='InitialData', EXIST=ierr)                       !
      IF (.NOT.ierr) THEN                                            !
         PRINT *, 'Initializing unknowns in LPM'                     !
         OPEN (1, FILE='InitialData',STATUS='NEW',FORM='UNFORMATTED')!
         WRITE (1) 0D0                                               !
         DO i=1, nUnknowns                                           !
            WRITE (1) tZeroX(i)                                      !
         END DO                                                      !
         CLOSE(1)                                                    !
      END IF                                                         !
c--------------------------------------------------------------------

c Surface to X pointer: this defines which Unknown corresponds to which
c suface in "List of Neumann Surfaces" inside solver.inp
c For example, if you have "List of Neumann Surfaces= 2 8 4 ...."
c and you set "srfToXPtr = (/5,3,9,.../)"
C this means X(5) corresponds to surface 2, X(3) <-> surface 8,
c and X(9) <-> surface 4
c      srfToXPtr = (/7/)
c This is exactly the same for srfToXdPtr only srfToXdPtr is for Dirichlet 
c surfaces
c      srfToXdPtr = (/1/)
      srfToXPtr = (/1/)

      END SUBROUTINE INITIALIZE

c####################################################################
c Here you should find the f_i=dx_i/dt, based on the following parameters:
c  Current x_i:                   x(i)
c  Current time:                  t
c  Flowrates of Neumann faces:    Q(i)
c  Pressure of Dirichlet faces: P(i)

      SUBROUTINE FINDF(t, x, f, Q, P)
      USE COM
      IMPLICIT NONE
      INTENT(IN) t, Q
      INTENT(OUT) f

      REAL(KIND=8) t, x(nUnknowns), f(nUnknowns), Q(nNeumannSrfs),
     2   P(nDirichletSrfs), pi

c RCR parameters
      REAL(8) Rp, C, Rd

      pi = ATAN(1D0)*4D0
       
      Rp = 121D0
      C  = 1.5D-4
      Rd = 1212D0

      f(1) = (1D0/C) * (Q(1) - x(1)/Rd)
      offset(1) = Q(1)*Rp
       
      Xprint(1)=t
      Xprint(2)= Q(1)
      Xprint(3) = offset(1)

      RETURN
      END SUBROUTINE FINDF


