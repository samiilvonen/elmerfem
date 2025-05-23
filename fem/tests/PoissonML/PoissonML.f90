SUBROUTINE PoissonSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve the Poisson equation without constructing a global matrix!
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh, materials, BCs, etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear & nonlinear equation solver options
!
!  REAL(KIND=dp) :: dt,
!     INPUT: Timestep size for time dependent simulations
!
!  LOGICAL :: TransientSimulation
!     INPUT: Steady state or transient simulation
!
!******************************************************************************
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  LOGICAL :: AllocationsDone = .FALSE., Found
  TYPE(Element_t),POINTER :: Element

  REAL(KIND=dp) :: Norm, d_val, diag, rt, ct, eps
  INTEGER :: sz, i, j, k, l, m, n, nelem, nb, nd, t, istat, active, maxi
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: BodyForce, BC
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), LOAD(:), FORCE(:)

  REAL(KIND=dp), POINTER :: b(:), x(:)

  INTEGER, ALLOCATABLE :: inds(:)

  ! structure listing degrees of freedom
  TYPE np_t
    REAL(KIND=dp) :: diag = 0._dp
    INTEGER :: eCount = 0
    INTEGER, ALLOCATABLE :: Elems(:)
  END TYPE np_t

  ! structure listing elements
  TYPE ed_t 
    INTEGER :: nd
    INTEGER, ALLOCATABLE :: dofIndeces(:)
    REAL(KIND=dp), ALLOCATABLE :: stiff(:,:), force(:)
  END TYPE ed_t

  TYPE(np_t), ALLOCATABLE :: np(:)
  INTEGER, ALLOCATABLE :: tp(:)
  TYPE(ed_t), ALLOCATABLE, TARGET :: ed(:)

  SAVE STIFF, LOAD, FORCE, AllocationsDone
!------------------------------------------------------------------------------

  !Allocate some permanent storage, this is done first time only:
  !--------------------------------------------------------------
  Mesh => GetMesh()

  IF ( .NOT. AllocationsDone ) THEN
     N = Solver % Mesh % MaxElementDOFs  ! just big enough for elemental arrays
     ALLOCATE( FORCE(N), LOAD(N), STIFF(N,N), STAT=istat )
     IF ( istat /= 0 ) THEN
        CALL Fatal( 'PoissonSolve', 'Memory allocation error.' )
     END IF
     AllocationsDone = .TRUE.
  END IF

   Active = GetNOFActive()

   n = SIZE(Solver % Variable % Values)
   ALLOCATE(ed(Active), np(n))
   DO i=1,n
     ALLOCATE(np(i) % Elems(4))
   END DO

   !System assembly:
   !----------------
   CALL DefaultInitialize()
   DO t=1,Active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      nb = GetElementNOFBDOFs()

      LOAD = 0.0d0
      BodyForce => GetBodyForce()
      IF ( ASSOCIATED(BodyForce) ) THEN
         Load(1:n) = GetReal( BodyForce, 'Source', Found )
      END IF

      !Get element local matrix and rhs vector:
      !----------------------------------------
      CALL LocalMatrix(  STIFF, FORCE, LOAD, Element, n, nd+nb )
      CALL LCondensate( nd, nb, STIFF, FORCE )

      ed(t) % Stiff = STIFF(1:nd,1:nd)
      ed(t) % Force = FORCE(1:nd)
      ed(t) % dofIndeces = [(i,i=1,nd)]
      nd =  GetElementDOFs(ed(t) % dofIndeces, Element)
      ed(t) % dofIndeces = solver % variable % perm(ed(t) % dofIndeces)
      DO i=1,nd
        j = ed(t) % dofIndeces(i)
        sz = np(j) % eCount
        IF ( sz>=SIZE(np(j) % Elems) ) THEN
          tp = np(j) % Elems(1:sz)
          DEALLOCATE(np(j) % Elems)
          ALLOCATE(np(j) % Elems(sz*2))
          np(j) % Elems(1:sz) = tp
        END IF
        np(j) % Elems(sz+1)=t
        np(j) % eCount = sz+1
        np(j) % diag = np(j) % diag + Stiff(i,i)
      END DO
   END DO

   ! Dirichlet-conditions
!$omp parallel do private(t,n,nelem,nd,nb,i,j,k,l,m,inds,found,element,BC,d_val,diag)
   DO t=1,GetNOFBoundaryElements()
     Element => GetBoundaryElement(t)

     BC => GetBC()
     IF (.NOT.ASSOCIATED(BC)) CYCLE

     d_val = GetCReal( BC, Solver % Variable % Name, Found )
     IF (.NOT. Found ) CYCLE

     n  = GetElementNOFNodes(Element)
     nd = GetElementNofDOFs(Element)

     inds = [(i,i=1,nd)]
     nd = GetElementDOFs(inds, Element)
     inds = solver % variable % perm(inds)

     DO j=1,nd
       nelem = np(inds(j)) % eCount
       diag = np(inds(j)) % diag
       DO k=1,nelem
         m = np(inds(j)) % elems(k)
         DO l=1,SIZE(ed(m) % dofIndeces)
           IF(ed(m) % dofIndeces(l)==inds(j)) EXIT
         END DO
         IF(l>SIZE(ed(m) % dofIndeces)) stop 'l'

         IF (j<=n) THEN
           ed(m) % force = ed(m) % force - ed(m) % stiff(:,l)*d_val
         END IF
         ed(m) % stiff(:,l) = 0._dp
         ed(m) % stiff(l,:) = 0._dp
         ed(m) % stiff(l,l) = diag/nelem
         IF(j<=n) THEN
           ed(m) % force(l) = d_val * diag/nelem
         ELSE
           ed(m) % force(l) = 0
         ENDIF
       END DO
     END DO
   END DO
!$omp end parallel do

   b => Solver % Matrix % RHS
   x => Solver % Variable % Values
   n = SIZE(Solver % Variable % Values)

   b = 0._dp
   DO i=1,SIZE(ed)
     nd = SIZE(ed(i) % dofIndeces)
     DO j=1,nd
       k = ed(i) % dofIndeces(j)
       b(k) = b(k) + ed(i) % force(j)
     END DO
   END DO

   ct = Cputime(); rt=RealTime();

   eps  = GetCReal( GetSolverParams(), 'Linear System Convergence Tolerance', Found )
   maxi = GetInteger( GetSolverParams(), 'Linear System Max Iterations', Found )
   CALL PCG( n, x, b, eps, maxi )

   print*,'Linear System Timing: ', cputime()-ct,realtime()-rt

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, LOAD, Element, n, nd )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:), LOAD(:)
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    LOGICAL :: Stat
    INTEGER :: i,t
    TYPE(GaussIntegrationPoints_t) :: IP

    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------
    CALL GetElementNodes( Nodes )
    STIFF = 0.0d0
    FORCE = 0.0d0

    !Numerical integration:
    !----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
       IP % W(t), detJ, Basis, dBasisdx )

      ! The source term at the integration point:
      !------------------------------------------
      LoadAtIP = SUM( Basis(1:n) * LOAD(1:n) )

      ! Finally, the elemental matrix & vector:
      !----------------------------------------
      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + IP % s(t) * DetJ * &
             MATMUL( dBasisdx, TRANSPOSE( dBasisdx ) )
      FORCE(1:nd) = FORCE(1:nd) + IP % s(t) * DetJ * LoadAtIP * Basis(1:nd)
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

! Tailored local CG algo for speed testing (somewhat faster than any of the 
! library routines but not so much...)
!-------------------------------------------------------------------------
  SUBROUTINE  pcg( n, x, b, eps, maxiter )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: x(n),b(n), eps
    INTEGER :: n, maxiter
    TYPE(Matrix_t), POINTER :: A

    REAL(KIND=dp):: alpha, beta, rho, oldrho
    REAL(KIND=dp) :: r(n), p(n), q(n), z(n), s
    INTEGER :: iter, i, j, k
    REAL(KIND=dp) :: residual, eps2,st

    eps2 = eps*eps
    call mv(n,x,r)
    r = b - r

    residual = SUM(r*r)
    IF(residual<eps2) RETURN

    call prec(n,r,z)
    p = z
    rho = SUM(r*z)

    DO iter=1,maxiter
      oldrho = rho
      call mv(n,p,q)
      alpha = rho/SUM(p*q)

      x = x + alpha * p
      r = r - alpha * q

      residual = SUM(r*r)
      IF(MOD(iter,10)==0) WRITE (*, '(I8, E11.4)') iter, SQRT(residual)
      IF (residual < eps2) EXIT

      call prec(n,r,z)
      rho = SUM(r*z)
      beta = rho / oldrho
      p = z + beta * p
    END DO

    call mv(n,x,r)
    r = b - r
    residual = SQRT(SUM(r*r))
    WRITE (*, '(I8, E11.4)') iter, residual
  END SUBROUTINE pcg


  SUBROUTINE mv(n,u,v)
    REAL(KIND=dp) :: u(n), v(n)
    INTEGER :: i,j,k,l,m,n,nd
    INTEGER, POINTER :: inds(:)

     v(1:n) = 0._dp
!$omp parallel do private(i,inds) shared(ed,u,v)
     DO i=1,SIZE(ed)
       inds => ed(i)  % dofIndeces
       v(inds) = v(inds) + MATMUL(ed(i) % stiff,u(inds))
     END DO
!$omp end parallel do
  END SUBROUTINE mv


  SUBROUTINE prec(n,u,v)
    INTEGER :: i,j,k,l,n
    REAL(KIND=dp) :: u(n), v(n), dval

!!omp parallel do private(i,j,inds) shared(ed,u,v)
!      DO i=1,SIZE(ed)
!        inds = ed(i) % dofIndeces
!        DO j=1,size(inds)
!          v(inds(j)) = u(inds(j)) / ed(i) % stiff(j,j)
!        END DO
!      END DO
!!omp end parallel do

!$omp parallel do private(i) shared(np)
    DO i=1,n
      v(i) = u(i) / np(i) % diag
    END DO
!$omp end parallel do
  END SUBROUTINE prec

!------------------------------------------------------------------------------
  SUBROUTINE LCondensate( N, Nb, K, F )
!------------------------------------------------------------------------------
    USE LinearAlgebra
    INTEGER :: N, Nb
    REAL(KIND=dp) :: K(:,:),F(:),Kbb(Nb,Nb), &
         Kbl(Nb,N), Klb(N,Nb), Fb(Nb)

    INTEGER :: m, i, j, l, p, Ldofs(N), Bdofs(Nb)

    IF ( Nb <= 0 ) RETURN

    Ldofs = (/ (i, i=1,n) /)
    Bdofs = (/ (i, i=n+1,n+nb) /)

    Kbb = K(Bdofs,Bdofs)
    Kbl = K(Bdofs,Ldofs)
    Klb = K(Ldofs,Bdofs)
    Fb  = F(Bdofs)

    CALL InvertMatrix( Kbb,nb )

    F(1:n) = F(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    K(1:n,1:n) = &
         K(1:n,1:n) - MATMUL( Klb, MATMUL( Kbb, Kbl ) )
!------------------------------------------------------------------------------
  END SUBROUTINE LCondensate
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------

           

!------------------------------------------------------------------------------
END SUBROUTINE PoissonSolver
!------------------------------------------------------------------------------
