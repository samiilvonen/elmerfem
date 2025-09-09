!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Module for solving the nodal preconditioning equation associated to AV solver.
! *  It is assumed that the matrix is asssembled by the AV solver and here it is
! *  ready to be used.
! *
! *  Authors: Peter Råback, Mika Malinen
! *  Email:   Peter.Raback@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 3.9.2025
! *
! *****************************************************************************/

!> \ingroup Solvers
!> \{

!------------------------------------------------------------------------------
SUBROUTINE APrecSolver_Init( Model,Solver,dt,Transient ) ! {{{
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver       !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model         !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt            !< Timestep size for time dependent simulations
  LOGICAL :: Transient           !< Steady state or transient simulation
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Found, ElectroDynamics
  LOGICAL :: Monolithic
  CHARACTER(*), PARAMETER :: Caller = 'APrecSolver_Init'
  
  Params => GetSolverParams()
  CALL ListAddNewLogical( Params,'Apply Mortar BCs',.TRUE.)
  CALL ListAddNewLogical( Params,'Use Global Mass Matrix',.TRUE.)
  CALL ListAddNewString( Params,'Exec Solver','never')
  CALL ListAddNewLogical( Params,'Skip Compute Nonlinear Change',.TRUE.)

  Monolithic = ListGetLogical( Params,'Monolithic Solver',Found )  
  IF( Monolithic ) THEN
    CALL ListAddNewString( Params,'Variable','-dofs 3 Nodal A' )
  ELSE
    CALL ListAddNewString( Params,'Variable','-nooutput nodal A tmp')
    CALL ListAddString( Params,&
        NextFreeKeyword('Exported Variable', Params),'-dofs 3 nodal A')
    CALL ListAddString( Params,&
        NextFreeKeyword('Exported Variable', Params),'-dofs 3 nodal A cum')
    CALL ListAddString( Params,&
        NextFreeKeyword('Exported Variable', Params),'-dofs 3 nodal A rhs')
  END IF
    
!------------------------------------------------------------------------------
END SUBROUTINE APrecSolver_Init ! }}}
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Solve the magnetic vector potential expressed in terms of a single component.
!> The solver may take into account rotating boundary conditions.
!> Also optionally compute moments and inertia. 
!------------------------------------------------------------------------------
SUBROUTINE APrecSolver( Model,Solver,dt,Transient ) ! {{{
!------------------------------------------------------------------------------
  USE DefUtils
  USE CircuitUtils
  USE ZirkaUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver       !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model         !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt            !< Timestep size for time dependent simulations
  LOGICAL :: Transient           !< Steady state or transient simulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  LOGICAL :: Found
  TYPE(Element_t), POINTER :: Element
  REAL(KIND=dp) :: Norm
  INTEGER :: i,j,k,n, nb, nd, t
  TYPE(ValueList_t), POINTER :: BC
  TYPE(Mesh_t),   POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Variable_t), POINTER :: Avar, Svar, NodeResVar, EdgeSolVar, EdgeResVar, pVar
  TYPE(Matrix_t), POINTER, SAVE :: Proj => NULL()
  LOGICAL :: Monolithic, SecondFamily, SecondOrder, PiolaVersion
  TYPE(ValueList_t), POINTER :: EdgeSolverParams
  CHARACTER(LEN=MAX_NAME_LEN) :: sname
  LOGICAL, SAVE :: Visited = .FALSE.
  REAL(KIND=dp), POINTER :: allrhs(:) => NULL(), BulkValues(:) => NULL()
  INTEGER :: comps, compi
  INTEGER :: NoVisited = 0
  LOGICAL, POINTER, SAVE :: NodeSkip(:)
  
  CHARACTER(*), PARAMETER :: Caller = 'APrecSolver'


  SAVE :: allrhs, BulkValues, NoVisited 
  
!------------------------------------------------------------------------------

  CALL Info( Caller,'-------------------------------------------------------', Level=10 )
  CALL Info( Caller,'Solving preconditioning equation for vector potential', Level=6 )
  CALL Info( Caller,'-------------------------------------------------------', Level=10 )

  Mesh => Solver % Mesh 
  SolverParams => Solver % Values
  SVar => Solver % Variable
  
  Monolithic = ListGetLogical( SolverParams,'Monolithic Solver', Found )
  NoVisited = NoVisited + 1

  IF(Monolithic) THEN
    ! Solve all 3 components at the same time!
    ! Note that the current assembly in AVSolver is not compatible with this!
    comps = 1
    allrhs => Solver % Matrix % rhs
  ELSE   
    ! Solve one component at a time -> faster!
    comps = 3
    NodeResVar => VariableGet( Mesh % Variables,'nodal a rhs')
    allrhs => NodeResVar % Values
    IF(.NOT. ASSOCIATED(Solver % Matrix % BulkValues)) THEN
      ALLOCATE(Solver % Matrix % BulkValues(SIZE(Solver % Matrix % Values)))
    END IF
    Solver % Matrix % BulkValues = Solver % Matrix % Values
  END IF
  BulkValues => Solver % Matrix % BulkValues

  
  IF(Monolithic) THEN
    AVar => SVar
  ELSE
    Avar => VariableGet( Mesh % Variables,'Nodal A')        
    IF(.NOT. ASSOCIATED(Avar)) THEN
      CALL Fatal(Caller,'Could not find variable "Nodal A"')
    END IF
    IF(SVar % dofs /= 1) THEN
      CALL Fatal(Caller,'Componentwise solver size should be 1!')
    END IF
  END IF
  IF(AVar % dofs /= 3) THEN
    CALL Fatal(Caller,'Full solution size should be 3!')
  END IF

  EdgeSolVar => NULL()
  sname = ListGetString( SolverParams, 'Edge Update Name', Found)
  IF (Found) THEN
    EdgeSolVar => VariableGet(Mesh % Variables, sname)
    IF (.NOT. ASSOCIATED(EdgeSolVar)) CALL Fatal(Caller, 'Could not found field: '//TRIM(sname))
  ELSE
    CALL Fatal(Caller, 'Give "Edge Update Name" to enable the use as a preconditioner')
  END IF

  EdgeResVar => NULL()  
  sname = ListGetString( SolverParams,'Edge Residual Name',Found)
  IF(.NOT. Found) THEN
    CALL Fatal(Caller, 'Give Edge Residual Name to enable the use as a preconditioner')
  END IF

  EdgeResVar => VariableGet( Mesh % Variables, sname )
  IF(.NOT. ASSOCIATED( EdgeResVar ) ) CALL Fatal(Caller,'Could not find field: '//TRIM(sname))

  EdgeSolverParams => GetSolverParams(EdgeResVar % Solver)
  CALL EdgeElementStyle(EdgeSolverParams, PiolaVersion, SecondFamily, SecondOrder, Check = .TRUE.)

  IF (.NOT. ASSOCIATED(Proj)) THEN
    CALL Info(Caller,'Creating projection matrix to map a nodal solution into vector element space', Level=10)
    CALL NodalToNedelecInterpolation_GlobalMatrix(Mesh, Avar, EdgeSolVar, Proj, cdim=3)
  END IF

 
  ! Now EdgeResVar represents the residual with respect
  ! to the basis for H(curl). We need to apply a transformation so that
  ! we may solve the residual correction equation by using the nodal basis.
  !-----------------------------------------------------------------------------
  CALL Info(Caller,'Using Transposed Projection Matrix: H(curl) -> H1', Level=10)
  CALL CRS_TransposeMatrixVectorMultiply(Proj, EdgeResVar % Values, allrhs )           

  
  IF(comps == 3 ) THEN  
    
    IF(NoVisited == 1 ) THEN
      n = SIZE(Solver % Matrix % rhs)
      ALLOCATE(NodeSkip(n))    
      CALL CreateNodeSkipMask(NodeSkip,Solver % Variable)
      n = COUNT(NodeSkip)
      IF(n==0) DEALLOCATE(NodeSkip)
    END IF
      
    DO compi = 1, comps      
      Solver % Matrix % Values = BulkValues
      Solver % Matrix % rhs = allrhs(compi::comps)

      ! By construction do not apply any residual to the mortar boundary. 
      IF(ASSOCIATED(NodeSkip)) THEN
        WHERE(NodeSkip)
          Solver % Matrix % rhs = 0.0_dp
        END WHERE
      END IF
        
      ! Different components will have different Dirichlet BC's!
      ! Note that Solver % Variable now points to correct component of Avar so there is no
      ! need to copy values to it!
      Solver % Variable => VariableGet( Model % Variables,TRIM(Avar % name)//' '//I2S(compi))
      IF(.NOT. ASSOCIATED(Solver % Variable)) THEN
        CALL Fatal(Caller,'Could not find variable for component :'//I2S(compi))
      END IF
      
      IF(ALLOCATED(Solver % Matrix % ConstrainedDOF ) ) &
          Solver % Matrix % ConstrainedDOF = .FALSE.
      CALL DefaultDirichletBCs()
      Norm = DefaultSolve()
    END DO
    Solver % Variable => SVar

    pVar => VariableGet( Mesh % Variables,'nodal a cum')
    IF(ASSOCIATED(pVar)) THEN
      pVar % Values = pVar % Values + AVar % Values
    END IF
  ELSE          
    CALL DefaultDirichletBCs()
    Norm = DefaultSolve()
  END IF

  CALL Info(Caller,'Projecting nodal solution to vector element space', Level=20)
  CALL CRS_MatrixVectorMultiply(Proj, Avar % Values, EdgeSolVar % Values ) 

  CALL Info(Caller,'Auxiliary space nodal solution finished!',Level=10)
    
END SUBROUTINE APrecSolver
!------------------------------------------------------------------------------

!> \}

