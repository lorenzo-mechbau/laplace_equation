PROGRAM LAPLACE_EQUATION

  USE OpenCMISS
  USE OpenCMISS_Iron

#ifndef NOMPIMOD
  USE MPI
#endif

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  !-----------------------------------------------------------------------------------------------------------
  ! PROGRAM VARIABLES AND TYPES
  !-----------------------------------------------------------------------------------------------------------
  
  !Test program parameters
  REAL(CMISSRP), PARAMETER :: HEIGHT=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: WIDTH=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: LENGTH=1.0_CMISSRP

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: LinearBasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: HermiteBasisUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=12
 
  !Program types
  
  !Program variables
  INTEGER(CMISSIntg) :: RowNo, ColumnNo, GlobalElementNo, GlobalNodeNo, I
  INTEGER(CMISSIntg) :: NumberOfElements, NumberOfNodes

  !CMISS variables
  TYPE(cmfe_BasisType) :: LinearBasis
  TYPE(cmfe_BasisType) :: HermiteBasis
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: Equations
  TYPE(cmfe_EquationsSetType) :: EquationsSet
  TYPE(cmfe_FieldType) :: GeometricField,EquationsSetField,DependentField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_NodesType) :: Nodes
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_SolverType) :: Solver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations
  TYPE(cmfe_MeshElementsType) :: MeshElements

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: FirstNodeNumber,LastNodeNumber
  INTEGER(CMISSIntg) :: FirstNodeDomain,LastNodeDomain
  INTEGER(CMISSIntg) :: Err
  LOGICAL :: SetDecompositionDistributed
  
  ! my variables
  INTEGER(CMISSIntg) :: count_components = 1  
  
#ifdef WIN32
  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif

  !Intialise OpenCMISS
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)

  !CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_OUTPUT_ERROR,Err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)

  CALL cmfe_RandomSeedsSet(9999,Err)
  
  !CALL cmfe_DiagnosticsSetOn(CMFE_IN_DIAG_TYPE,[1,2],"",&
  !  & ["DOMAIN_MAPPINGS_NODES_CALCULATE"],Err)

  !CALL cmfe_DiagnosticsSetOn(CMFE_FROM_DIAG_TYPE,[3],"",&
  !  & ["DECOMPOSITION_CREATE_FINISH"],Err)

  !CALL cmfe_DiagnosticsSetOn(CMFE_IN_DIAG_TYPE,[1,2],"",&
  !  & ["DECOMPOSITION_ADJACENT_DOMAINS_CALCULATE"],Err)

  CALL cmfe_DiagnosticsSetOn(CMFE_IN_DIAG_TYPE,[2],"",&
    & ["FIELD_MAPPINGS_CALCULATE"],Err)

  
  
  !CALL cmfe_OutputSetOn("diagnostics.txt",Err)

  !Get the computational nodes information
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  !-----------------------------------------------------------------------------------------------------------
  !COORDINATE SYSTEM
  !-----------------------------------------------------------------------------------------------------------  
 
  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,2,Err)
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !-----------------------------------------------------------------------------------------------------------
  !REGION
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL cmfe_Region_LabelSet(Region,"laplace_equation",Err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(Region,Err)

  !-----------------------------------------------------------------------------------------------------------
  !BASIS
  !-----------------------------------------------------------------------------------------------------------

  !Create linear basis
  CALL cmfe_Basis_Initialise(LinearBasis,Err)
  CALL cmfe_Basis_CreateStart(LinearBasisUserNumber,LinearBasis,Err)
  CALL cmfe_Basis_TypeSet(LinearBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(LinearBasis,2,Err)
  CALL cmfe_Basis_InterpolationXiSet(LinearBasis,&
    & [CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION,CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(LinearBasis,[2,2],Err)
  CALL cmfe_Basis_CreateFinish(LinearBasis,Err)
  
  !Create hermite basis
  CALL cmfe_Basis_Initialise(HermiteBasis,Err)
  CALL cmfe_Basis_CreateStart(HermiteBasisUserNumber,HermiteBasis,Err)
  CALL cmfe_Basis_TypeSet(HermiteBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(HermiteBasis,2,Err)
  CALL cmfe_Basis_InterpolationXiSet(HermiteBasis,&
    & [CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION,CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION],Err)
  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(HermiteBasis,[3,3],Err)
  CALL cmfe_Basis_CreateFinish(HermiteBasis,Err)
   
  !-----------------------------------------------------------------------------------------------------------
  !MESH
  !-----------------------------------------------------------------------------------------------------------

 ! !Start the creation of a generated mesh in the region
 ! CALL cmfe_GeneratedMesh_Initialise(GeneratedMesh,Err)
 ! CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
 ! !Set up a regular x*y*z mesh
 ! CALL cmfe_GeneratedMesh_TypeSet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
 ! !Set the default basis
 ! CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,Basis,Err)   
 ! !Define the mesh on the region
 ! IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
 !   CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT],Err)
 !   CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS],Err)
 ! ELSE
 !   CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT,LENGTH],Err)
 !   CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
 !     & NUMBER_GLOBAL_Z_ELEMENTS],Err)
 ! ENDIF    
 ! !Finish the creation of a generated mesh in the region
 ! CALL cmfe_Mesh_Initialise(Mesh,Err)
 ! CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

  ! Create mesh as specified in "PartitionedMesh_4x4.pdf"
  
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_Mesh_CreateStart(MeshUserNumber,Region,2,Mesh,Err)
  CALL cmfe_Mesh_NumberOfComponentsSet(Mesh,1,Err)
  NumberOfElements = 16
  IF (NumberOfComputationalNodes == 3) NumberOfElements = 3
  CALL cmfe_Mesh_NumberOfElementsSet(Mesh,NumberOfElements,Err)
  
  ! Create elements 
  CALL cmfe_MeshElements_Initialise(MeshElements,Err)
  CALL cmfe_MeshElements_CreateStart(Mesh,1,LinearBasis,MeshElements,Err)
  
  ! Define nodes for the mesh
  CALL cmfe_Nodes_Initialise(Nodes,Err)
  NumberOfNodes = 25
! This line does not make much sense:  
!  IF (NumberOfComputationalNodes == 3) NumberOfElements = 4
  IF (NumberOfComputationalNodes == 3) NumberOfNodes = 8
  CALL cmfe_Nodes_CreateStart(Region,NumberOfNodes,Nodes,Err)
  CALL cmfe_Nodes_CreateFinish(Nodes,Err)
  
  ! Set adjacent nodes for each element
  IF (NumberOfComputationalNodes == 3) THEN
    CALL cmfe_MeshElements_NodesSet(MeshElements, 1, [1,2,5,6], Err)
    CALL cmfe_MeshElements_NodesSet(MeshElements, 2, [2,3,6,7], Err)
    CALL cmfe_MeshElements_NodesSet(MeshElements, 3, [3,4,7,8], Err)
    CALL cmfe_MeshElements_BasisSet(MeshElements,1,LinearBasis,Err)
    CALL cmfe_MeshElements_BasisSet(MeshElements,2,HermiteBasis,Err)
    CALL cmfe_MeshElements_BasisSet(MeshElements,3,LinearBasis,Err)
  ELSE
    DO RowNo = 1,4
      DO ColumnNo = 1,4
        GlobalElementNo = (RowNo-1)*4 + ColumnNo
        GlobalNodeNo = (RowNo-1)*5 + ColumnNo
        ! OK
        CALL cmfe_MeshElements_NodesSet(MeshElements, GlobalElementNo, &
          & [GlobalNodeNo,GlobalNodeNo+1,GlobalNodeNo+5,GlobalNodeNo+6], Err)
         ! PRINT *, "Element"
         ! PRINT *, GlobalElementNo
         ! PRINT *, "Nodes"
         ! PRINT *, [GlobalNodeNo,GlobalNodeNo+1,GlobalNodeNo+5,GlobalNodeNo+6]          
      ENDDO
    ENDDO

    ! Set the basis
    DO GlobalElementNo=1,4
      CALL cmfe_MeshElements_BasisSet(MeshElements,GlobalElementNo,HermiteBasis,Err)
    ENDDO
    DO GlobalElementNo=5,12
      CALL cmfe_MeshElements_BasisSet(MeshElements,GlobalElementNo,LinearBasis,Err)
    ENDDO
    DO GlobalElementNo=13,16
      CALL cmfe_MeshElements_BasisSet(MeshElements,GlobalElementNo,HermiteBasis,Err)
    ENDDO
  ENDIF
  
  CALL cmfe_MeshElements_CreateFinish(MeshElements,Err)
  
  CALL cmfe_Mesh_CreateFinish(Mesh,Err)
  
  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
! Set the decomposition to be a general decomposition with the specified number of domains
  
  !CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_USER_DEFINED_TYPE,Err)
  
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  
  ! global element numbers
  !13 14 15 16
  ! 9 10 11 12
  ! 5  6  7  8
  ! 1  2  3  4
  
  PRINT *, "Computational node number:"
  PRINT *, ComputationalNodeNumber
  
  SetDecompositionDistributed = .FALSE.   ! true only works with new implementation, false works with both
  IF (NumberOfComputationalNodes == 1) THEN
    DO I=1,16
      CALL cmfe_Decomposition_ElementDomainSet(Decomposition,I,0,Err)
    ENDDO
  ELSEIF (NumberOfComputationalNodes == 2) THEN
    CALL cmfe_Decomposition_ElementDomainSet(Decomposition,1,0,Err)
    CALL cmfe_Decomposition_ElementDomainSet(Decomposition,2,0,Err)
    CALL cmfe_Decomposition_ElementDomainSet(Decomposition,5,0,Err)
    CALL cmfe_Decomposition_ElementDomainSet(Decomposition,6,0,Err)
    CALL cmfe_Decomposition_ElementDomainSet(Decomposition,9,0,Err)
    CALL cmfe_Decomposition_ElementDomainSet(Decomposition,10,0,Err)
    CALL cmfe_Decomposition_ElementDomainSet(Decomposition,13,0,Err)
    CALL cmfe_Decomposition_ElementDomainSet(Decomposition,14,0,Err)
    CALL cmfe_Decomposition_ElementDomainSet(Decomposition,3,1,Err)
    CALL cmfe_Decomposition_ElementDomainSet(Decomposition,4,1,Err)
    CALL cmfe_Decomposition_ElementDomainSet(Decomposition,7,1,Err)
    CALL cmfe_Decomposition_ElementDomainSet(Decomposition,8,1,Err)
    CALL cmfe_Decomposition_ElementDomainSet(Decomposition,11,1,Err)
    CALL cmfe_Decomposition_ElementDomainSet(Decomposition,12,1,Err)
    CALL cmfe_Decomposition_ElementDomainSet(Decomposition,15,1,Err)
    CALL cmfe_Decomposition_ElementDomainSet(Decomposition,16,1,Err)
  ELSEIF (NumberOfComputationalNodes == 3) THEN
    CALL cmfe_Decomposition_ElementDomainSet(Decomposition,1,0,Err)
    CALL cmfe_Decomposition_ElementDomainSet(Decomposition,2,1,Err)
    CALL cmfe_Decomposition_ElementDomainSet(Decomposition,3,2,Err)
  ELSEIF (NumberOfComputationalNodes >= 4) THEN
    IF (ComputationalNodeNumber == 0 .OR..NOT.SetDecompositionDistributed) THEN
      !                                             global el., domain
      CALL cmfe_Decomposition_ElementDomainSet(Decomposition,9,2,Err)
      CALL cmfe_Decomposition_ElementDomainSet(Decomposition,5,0,Err)
      CALL cmfe_Decomposition_ElementDomainSet(Decomposition,14,2,Err)
      CALL cmfe_Decomposition_ElementDomainSet(Decomposition,11,3,Err)
      CALL cmfe_Decomposition_ElementDomainSet(Decomposition,6,0,Err)
      CALL cmfe_Decomposition_ElementDomainSet(Decomposition,16,3,Err)
    ENDIF
    
    IF (ComputationalNodeNumber == 1 .OR..NOT.SetDecompositionDistributed) THEN
      CALL cmfe_Decomposition_ElementDomainSet(Decomposition,7,0,Err)
      CALL cmfe_Decomposition_ElementDomainSet(Decomposition,15,3,Err)
      CALL cmfe_Decomposition_ElementDomainSet(Decomposition,1,0,Err)
    ENDIF
    
    IF (ComputationalNodeNumber == 2 .OR..NOT.SetDecompositionDistributed) THEN
      CALL cmfe_Decomposition_ElementDomainSet(Decomposition,10,2,Err)
      CALL cmfe_Decomposition_ElementDomainSet(Decomposition,8,1,Err)
    ENDIF
    
    IF (ComputationalNodeNumber == 3 .OR..NOT.SetDecompositionDistributed) THEN
      CALL cmfe_Decomposition_ElementDomainSet(Decomposition,2,0,Err)
      CALL cmfe_Decomposition_ElementDomainSet(Decomposition,3,1,Err)
      CALL cmfe_Decomposition_ElementDomainSet(Decomposition,12,3,Err)
      CALL cmfe_Decomposition_ElementDomainSet(Decomposition,4,1,Err)
      CALL cmfe_Decomposition_ElementDomainSet(Decomposition,13,2,Err)
      
    ENDIF
  ENDIF
  
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)
 
  !Destory the mesh now that we have decomposed it
! CALL cmfe_Mesh_Destroy(Mesh,Err)

  !Variable, MaxDepth, MaxArrayLength
  !CALL cmfe_PrintDecomposition(Decomposition,3,100,Err)
  
  
  
  !PRINT *, "Abort program in laplace_equation.f90:329"
  !STOP
 
  !-----------------------------------------------------------------------------------------------------------
  !GEOMETRIC FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Start to create a default (geometric) field on the region
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the domain to be used by the field components.
 
  !!CALL cmfe_Field_NumberOfComponentsSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,7,Err)
 
  !1 vbl, 2 components: 2x Mesh (node-dofs)
  !CALL cmfe_Field_NumberOfComponentsSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,Err)
 
  !1 vbl, 3 components: 2x Mesh (node-dofs), 1x Interpolation (element-dofs)
  CALL cmfe_Field_NumberOfComponentsSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,Err)

  !1 vbl, 4 components: 2x Mesh (node-dofs), 2x Interpolation (1 element-dofs, 1 constant = 1 dof)
  !CALL cmfe_Field_NumberOfComponentsSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,4,Err)

  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,count_components,1,Err)
  count_components = count_components +1
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,count_components,1,Err)
  count_components = count_components +1
  CALL cmfe_Field_ComponentInterpolationSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,count_components, &
  & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
!  count_components = count_components +1
!  CALL cmfe_Field_ComponentInterpolationSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,count_components,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)

  !!CALL cmfe_Field_ComponentInterpolationSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,5,CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  !!CALL cmfe_Field_ComponentInterpolationSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,6,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  !!CALL cmfe_Field_ComponentInterpolationSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,7,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(GeometricField,Err)

  PRINT *, ""
  !PRINT *, "Print field GeometricField"
  !CALL cmfe_PrintField(GeometricField,3,5,Err)
  
  PRINT *, "Abort program in laplace_equation.f90:360"
  STOP
 
  
  !Update the geometric field parameters
  !CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)
  
  !-----------------------------------------------------------------------------------------------------------
  !EQUATIONS SETS
  !-----------------------------------------------------------------------------------------------------------

  !Create the Standard Laplace Equations set
  CALL cmfe_EquationsSet_Initialise(EquationsSet,Err)
  CALL cmfe_Field_Initialise(EquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumber,Region,GeometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_LAPLACE_EQUATION_TYPE,CMFE_EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE],EquationsSetFieldUserNumber, &
    & EquationsSetField,EquationsSet,Err)
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(EquationsSet,Err)

  !-----------------------------------------------------------------------------------------------------------
  ! DEPENDENT FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set dependent field variables
  CALL cmfe_Field_Initialise(DependentField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  !Set the DOFs to be contiguous across components
  CALL cmfe_Field_DOFOrderTypeSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_SEPARATED_COMPONENT_DOF_ORDER,Err)
  CALL cmfe_Field_DOFOrderTypeSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,CMFE_FIELD_SEPARATED_COMPONENT_DOF_ORDER,Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Initialise the field with an initial guess
  CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,0.5_CMISSRP, &
    & Err)

  !-----------------------------------------------------------------------------------------------------------
  ! EQUATIONS
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(Equations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
! CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_FULL_MATRICES,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_NO_OUTPUT,Err)
! CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_TIMING_OUTPUT,Err)
! CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_MATRIX_OUTPUT,Err)
! CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet,Err)

  !-----------------------------------------------------------------------------------------------------------
  !PROBLEM
  !-----------------------------------------------------------------------------------------------------------  

  !Start the creation of a problem.
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_CLASSICAL_FIELD_CLASS,CMFE_PROBLEM_LAPLACE_EQUATION_TYPE, &
    & CMFE_PROBLEM_STANDARD_LAPLACE_SUBTYPE],Problem,Err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(Problem,Err)

  !Start the creation of the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)
 
  !-----------------------------------------------------------------------------------------------------------
  !SOLVER
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the problem solvers
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_NO_OUTPUT,Err)
! CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
! CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_TIMING_OUTPUT,Err)
! CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_SOLVER_OUTPUT,Err)
! CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_MATRIX_OUTPUT,Err)
  
! CALL cmfe_Solver_LinearTypeSet(Solver,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
! CALL cmfe_Solver_LinearIterativeAbsoluteToleranceSet(Solver,1.0E-12_CMISSRP,Err)
! CALL cmfe_Solver_LinearIterativeRelativeToleranceSet(Solver,1.0E-12_CMISSRP,Err)

  CALL cmfe_Solver_LinearTypeSet(Solver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  
! CALL cmfe_Solver_LinearTypeSet(Solver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
! CALL cmfe_Solver_LibraryTypeSet(Solver,CMFE_SOLVER_MUMPS_LIBRARY,Err)
! CALL cmfe_Solver_LibraryTypeSet(Solver,CMFE_SOLVER_LAPACK_LIBRARY,Err)
! CALL cmfe_Solver_LibraryTypeSet(Solver,CMFE_SOLVER_SUPERLU_LIBRARY,Err)
! CALL cmfe_Solver_LibraryTypeSet(Solver,CMFE_SOLVER_PASTIX_LIBRARY,Err)
  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

  !-----------------------------------------------------------------------------------------------------------
  !SOLVER EQUATIONS
  !-----------------------------------------------------------------------------------------------------------  

  !Start the creation of the problem solver equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquations,Err)
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
  !Get the solve equations
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
  !Set the solver equations sparsity
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_SPARSE_MATRICES,Err)
  !CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_FULL_MATRICES,Err)  
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  !-----------------------------------------------------------------------------------------------------------
  !BOUNDARY CONDITIONS
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the equations set boundary conditions
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  !Set the first node to 0.0 and the last node to 1.0
  FirstNodeNumber=1
  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Region_NodesGet(Region,Nodes,Err)
  CALL cmfe_Nodes_NumberOfNodesGet(Nodes,LastNodeNumber,Err)
  CALL cmfe_Decomposition_NodeDomainGet(Decomposition,FirstNodeNumber,1,FirstNodeDomain,Err)
  CALL cmfe_Decomposition_NodeDomainGet(Decomposition,LastNodeNumber,1,LastNodeDomain,Err)
  IF(FirstNodeDomain==ComputationalNodeNumber) THEN
    CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,FirstNodeNumber,1, &
      & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
  ENDIF
  IF(LastNodeDomain==ComputationalNodeNumber) THEN
    CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,LastNodeNumber,1, &
      & CMFE_BOUNDARY_CONDITION_FIXED,1.0_CMISSRP,Err)
  ENDIF
  !Finish the creation of the equations set boundary conditions
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !-----------------------------------------------------------------------------------------------------------
  !SOLVE
  !-----------------------------------------------------------------------------------------------------------

  !Solve the problem
  CALL cmfe_Problem_Solve(Problem,Err)

  !-----------------------------------------------------------------------------------------------------------
  !OUTPUT
  !-----------------------------------------------------------------------------------------------------------

  !Export results
  CALL cmfe_Fields_Initialise(Fields,Err)
  CALL cmfe_Fields_Create(Region,Fields,Err)
  CALL cmfe_Fields_NodesExport(Fields,"laplace_equation","FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(Fields,"laplace_equation","FORTRAN",Err)
  CALL cmfe_Fields_Finalise(Fields,Err)
  
  !Finialise CMISS
  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
CONTAINS

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)

    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP

  END SUBROUTINE HANDLE_ERROR
    
END PROGRAM LAPLACE_EQUATION
