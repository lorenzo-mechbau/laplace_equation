PROGRAM LaplaceEquation

  USE OpenCMISS
  USE OpenCMISS_Iron
  
  IMPLICIT NONE

  !-----------------------------------------------------------------------------------------------------------
  ! PROGRAM VARIABLES AND TYPES
  !-----------------------------------------------------------------------------------------------------------

  !Test program parameters
  REAL(CMISSRP), PARAMETER :: HEIGHT=1.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: WIDTH=2.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: LENGTH=3.0_CMISSRP
 
  INTEGER(CMISSIntg), PARAMETER :: COORDINATE_SYSTEM_USER_NUMBER=1
  INTEGER(CMISSIntg), PARAMETER :: REGION_USER_NUMBER=2
  INTEGER(CMISSIntg), PARAMETER :: BASIS_USER_NUMBER=3
  INTEGER(CMISSIntg), PARAMETER :: GENERATED_MESH_USER_NUMBER=4
  INTEGER(CMISSIntg), PARAMETER :: MESH_USER_NUMBER=5
  INTEGER(CMISSIntg), PARAMETER :: DECOMPOSITION_USER_NUMBER=6
  INTEGER(CMISSIntg), PARAMETER :: GEOMETRIC_FIELD_USER_NUMBER=7
  INTEGER(CMISSIntg), PARAMETER :: EQUATIONS_SET_FIELD_USER_NUMBER=8
  INTEGER(CMISSIntg), PARAMETER :: DEPENDENT_FIELD_USER_NUMBER=9
  INTEGER(CMISSIntg), PARAMETER :: EQUATIONS_SET_USER_NUMBER=10
  INTEGER(CMISSIntg), PARAMETER :: PROBLEM_USER_NUMBER=11

  !Program types

  !Program variables
  INTEGER(CMISSIntg) :: numberOfArguments,argumentLength,status
  INTEGER(CMISSIntg) :: numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements, &
    & interpolationType,numberOfGaussXi
  CHARACTER(LEN=255) :: commandArgument,filename

  !CMISS variables

  !Generic CMISS variables
  INTEGER(CMISSIntg) :: contextUserNumber,worldRegionUserNumber
  INTEGER(CMISSIntg) :: numberOfComputationalNodes,computationalNodeNumber
  INTEGER(CMISSIntg) :: equationsSetIndex
  INTEGER(CMISSIntg) :: firstNodeNumber,lastNodeNumber
  INTEGER(CMISSIntg) :: firstNodeDomain,lastNodeDomain
  INTEGER(CMISSIntg) :: err

  !-----------------------------------------------------------------------------------------------------------
  ! PROBLEM CONTROL PANEL
  !-----------------------------------------------------------------------------------------------------------

  numberOfArguments = COMMAND_ARGUMENT_COUNT()
  IF(numberOfArguments >= 4) THEN
    !If we have enough arguments then use the first four for setting up the problem. The subsequent arguments may be used to
    !pass flags to, say, PETSc.
    CALL GET_COMMAND_ARGUMENT(1,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 1.")
    READ(commandArgument(1:argumentLength),*) numberOfGlobalXElements
    IF(numberOfGlobalXElements<=0) CALL HandleError("Invalid number of X elements.")
    CALL GET_COMMAND_ARGUMENT(2,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 2.")
    READ(commandArgument(1:argumentLength),*) numberOfGlobalYElements
    IF(numberOfGlobalYElements<=0) CALL HandleError("Invalid number of Y elements.")
    CALL GET_COMMAND_ARGUMENT(3,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 3.")
    READ(commandArgument(1:argumentLength),*) numberOfGlobalZElements
    IF(numberOfGlobalZElements<0) CALL HandleError("Invalid number of Z elements.")
    CALL GET_COMMAND_ARGUMENT(4,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 4.")
    READ(commandArgument(1:argumentLength),*) interpolationType
    IF(interpolationType<=0) CALL HandleError("Invalid Interpolation specification.")
  ELSE
    !If there are not enough arguments default the problem specification
    numberOfGlobalXElements=1
    numberOfGlobalYElements=3
    numberOfGlobalZElements=1
    interpolationType=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
!    interpolationType=CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
!    interpolationType=CMFE_BASIS_CUBIC_LAGRANGE_INTERPOLATION
  ENDIF

  !Intialise OpenCMISS
  CALL cmfe_Initialise(contextUserNumber,err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,err)
  CALL cmfe_Context_WorldRegionGet(contextUserNumber,worldRegionUserNumber,err)
  CALL cmfe_Context_RandomSeedsSet(contextUserNumber,9999,err)
  !CALL cmfe_DiagnosticsSetOn(CMFE_IN_DIAG_TYPE,[1,2,3,4,5],"Diagnostics",["Laplace_FiniteElementCalculate"],err)

  WRITE(filename,'(A,"_",I0,"x",I0,"x",I0,"_",I0)') "Laplace",numberOfGlobalXElements,numberOfGlobalYElements, &
    & numberOfGlobalZElements,interpolationType

  CALL cmfe_OutputSetOn(filename,err)

  !Get the computational nodes information
  CALL cmfe_ComputationEnvironment_NumberOfWorldNodesGet(contextUserNumber,numberOfComputationalNodes,err)
  CALL cmfe_ComputationEnvironment_WorldNodeNumberGet(contextUserNumber,computationalNodeNumber,err)

  !-----------------------------------------------------------------------------------------------------------
  ! COORDINATE SYSTEM
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,contextUserNumber,err)
  IF(numberOfGlobalZElements==0) THEN
    !Set the coordinate system to be 2D
    CALL cmfe_CoordinateSystem_DimensionSet(contextUserNumber,COORDINATE_SYSTEM_USER_NUMBER,2,err)
  ELSE
    !Set the coordinate system to be 3D
    CALL cmfe_CoordinateSystem_DimensionSet(contextUserNumber,COORDINATE_SYSTEM_USER_NUMBER,3,err)
  ENDIF
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(contextUserNumber,COORDINATE_SYSTEM_USER_NUMBER,err)

  !-----------------------------------------------------------------------------------------------------------
  ! REGION
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the region
  CALL cmfe_Region_CreateStart(REGION_USER_NUMBER,contextUserNumber,worldRegionUserNumber,err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL cmfe_Region_CoordinateSystemSet(contextUserNumber,REGION_USER_NUMBER,COORDINATE_SYSTEM_USER_NUMBER,err)
  !Set the region label
  CALL cmfe_Region_LabelSet(contextUserNumber,REGION_USER_NUMBER,"LaplaceEquation",err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(contextUserNumber,REGION_USER_NUMBER,err)

  !-----------------------------------------------------------------------------------------------------------
  ! BASIS
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a basis (default is trilinear lagrange)
  CALL cmfe_Basis_CreateStart(BASIS_USER_NUMBER,contextUserNumber,err)
  SELECT CASE(interpolationType)
  CASE(1,2,3,4)
    CALL cmfe_Basis_TypeSet(contextUserNumber,BASIS_USER_NUMBER,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
  CASE(7,8,9)
    CALL cmfe_Basis_TypeSet(contextUserNumber,BASIS_USER_NUMBER,CMFE_BASIS_SIMPLEX_TYPE,err)
  CASE DEFAULT
    CALL HandleError("Invalid interpolation type.")
  END SELECT
  SELECT CASE(interpolationType)
  CASE(1)
    numberOfGaussXi=2
  CASE(2)
    numberOfGaussXi=3
  CASE(3,4)
    numberOfGaussXi=4
  CASE DEFAULT
    numberOfGaussXi=0 !Don't set number of Gauss points for tri/tet
  END SELECT
  IF(numberOfGlobalZElements==0) THEN
    !Set the basis to be a bi-interpolation basis
    CALL cmfe_Basis_NumberOfXiSet(contextUserNumber,BASIS_USER_NUMBER,2,err)
    CALL cmfe_Basis_InterpolationXiSet(contextUserNumber,BASIS_USER_NUMBER,[interpolationType,interpolationType],err)
    IF(numberOfGaussXi>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(contextUserNumber,BASIS_USER_NUMBER, &
        & [numberOfGaussXi,numberOfGaussXi],err)
    ENDIF
  ELSE
    !Set the basis to be a tri-interpolation basis
    CALL cmfe_Basis_NumberOfXiSet(contextUserNumber,BASIS_USER_NUMBER,3,err)
    CALL cmfe_Basis_InterpolationXiSet(contextUserNumber,BASIS_USER_NUMBER, &
      & [interpolationType,interpolationType,interpolationType],err)
    IF(numberOfGaussXi>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(contextUserNumber,BASIS_USER_NUMBER, &
        & [numberOfGaussXi,numberOfGaussXi,numberOfGaussXi],err)
    ENDIF
  ENDIF
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(contextUserNumber,BASIS_USER_NUMBER,err)

  !-----------------------------------------------------------------------------------------------------------
  ! MESH
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a generated mesh in the region
  CALL cmfe_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,contextUserNumber,REGION_USER_NUMBER,err)
  !Set up a regular x*y*z mesh
  CALL cmfe_GeneratedMesh_TypeSet(contextUserNumber,REGION_USER_NUMBER,GENERATED_MESH_USER_NUMBER, &
    & CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,err)
  !Set the default basis
  CALL cmfe_GeneratedMesh_BasisSet(contextUserNumber,REGION_USER_NUMBER,GENERATED_MESH_USER_NUMBER,BASIS_USER_NUMBER,err)
  !Define the mesh on the region
  IF(numberOfGlobalZElements==0) THEN
    CALL cmfe_GeneratedMesh_ExtentSet(contextUserNumber,REGION_USER_NUMBER,GENERATED_MESH_USER_NUMBER,[WIDTH,HEIGHT],err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(contextUserNumber,REGION_USER_NUMBER,GENERATED_MESH_USER_NUMBER, &
      & [numberOfGlobalXElements,numberOfGlobalYElements],err)
  ELSE
    CALL cmfe_GeneratedMesh_ExtentSet(contextUserNumber,REGION_USER_NUMBER,GENERATED_MESH_USER_NUMBER,[WIDTH,HEIGHT,LENGTH],err)
    CALL cmfe_GeneratedMesh_NumberOfElementsSet(contextUserNumber,REGION_USER_NUMBER,GENERATED_MESH_USER_NUMBER, &
      & [numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements],err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL cmfe_GeneratedMesh_CreateFinish(contextUserNumber,REGION_USER_NUMBER,GENERATED_MESH_USER_NUMBER,MESH_USER_NUMBER,err)

  !Create a decomposition
  CALL cmfe_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,contextUserNumber,REGION_USER_NUMBER, &
    & MESH_USER_NUMBER,err)
! Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(contextUserNumber,REGION_USER_NUMBER,MESH_USER_NUMBER,DECOMPOSITION_USER_NUMBER, &
    & CMFE_DECOMPOSITION_CALCULATED_TYPE,err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(contextUserNumber,REGION_USER_NUMBER,MESH_USER_NUMBER,DECOMPOSITION_USER_NUMBER, &
    & numberOfComputationalNodes,err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(contextUserNumber,REGION_USER_NUMBER,MESH_USER_NUMBER,DECOMPOSITION_USER_NUMBER,err)

  !Destory the mesh now that we have decomposed it
! CALL cmfe_Mesh_Destroy(mesh,err)

  !-----------------------------------------------------------------------------------------------------------
  ! GEOMETRIC FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Start to create a default (geometric) field on the region
  CALL cmfe_Field_CreateStart(GEOMETRIC_FIELD_USER_NUMBER,contextUserNumber,REGION_USER_NUMBER,err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(contextUserNumber,REGION_USER_NUMBER,GEOMETRIC_FIELD_USER_NUMBER, &
    & MESH_USER_NUMBER,DECOMPOSITION_USER_NUMBER,err)
  !Set the domain to be used by the field components.
  CALL cmfe_Field_ComponentMeshComponentSet(contextUserNumber,REGION_USER_NUMBER,GEOMETRIC_FIELD_USER_NUMBER, &
    & CMFE_FIELD_U_VARIABLE_TYPE,1,1,err)
  CALL cmfe_Field_ComponentMeshComponentSet(contextUserNumber,REGION_USER_NUMBER,GEOMETRIC_FIELD_USER_NUMBER, &
    & CMFE_FIELD_U_VARIABLE_TYPE,2,1,err)
  IF(numberOfGlobalZElements/=0) THEN
    CALL cmfe_Field_ComponentMeshComponentSet(contextUserNumber,REGION_USER_NUMBER,GEOMETRIC_FIELD_USER_NUMBER, &
      & CMFE_FIELD_U_VARIABLE_TYPE,3,1,err)
  ENDIF
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(contextUserNumber,REGION_USER_NUMBER,GEOMETRIC_FIELD_USER_NUMBER,err)

  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(contextUserNumber,REGION_USER_NUMBER,GENERATED_MESH_USER_NUMBER, &
    & GEOMETRIC_FIELD_USER_NUMBER,err)

  !-----------------------------------------------------------------------------------------------------------
  ! EQUATIONS SETS
  !-----------------------------------------------------------------------------------------------------------

  !Create the Standard Laplace equations set
  CALL cmfe_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,contextUserNumber,REGION_USER_NUMBER,GEOMETRIC_FIELD_USER_NUMBER, &
    & [CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS,CMFE_EQUATIONS_SET_LAPLACE_EQUATION_TYPE, &
    & CMFE_EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE],EQUATIONS_SET_FIELD_USER_NUMBER,err)
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(contextUserNumber,REGION_USER_NUMBER,EQUATIONS_SET_USER_NUMBER,err)

  !-----------------------------------------------------------------------------------------------------------
  ! DEPENDENT FIELD
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateStart(contextUserNumber,REGION_USER_NUMBER,EQUATIONS_SET_USER_NUMBER, &
    & DEPENDENT_FIELD_USER_NUMBER,err)
  !Set the DOFs to be contiguous across components
  CALL cmfe_Field_DOFOrderTypeSet(contextUserNumber,REGION_USER_NUMBER,DEPENDENT_FIELD_USER_NUMBER, &
    & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_SEPARATED_COMPONENT_DOF_ORDER,err)
  CALL cmfe_Field_DOFOrderTypeSet(contextUserNumber,REGION_USER_NUMBER,DEPENDENT_FIELD_USER_NUMBER, &
    & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,CMFE_FIELD_SEPARATED_COMPONENT_DOF_ORDER,err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(contextUserNumber,REGION_USER_NUMBER,EQUATIONS_SET_USER_NUMBER,err)

  !Initialise the field with an initial guess
  CALL cmfe_Field_ComponentValuesInitialise(contextUserNumber,REGION_USER_NUMBER,DEPENDENT_FIELD_USER_NUMBER, &
    & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,0.5_CMISSRP,err)

  !-----------------------------------------------------------------------------------------------------------
  ! EQUATIONS
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateStart(contextUserNumber,REGION_USER_NUMBER,EQUATIONS_SET_USER_NUMBER,err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(contextUserNumber,REGION_USER_NUMBER,EQUATIONS_SET_USER_NUMBER, &
    & CMFE_EQUATIONS_SPARSE_MATRICES,err)
  !CALL cmfe_Equations_SparsityTypeSet(contextUserNumber,REGION_USER_NUMBER,EQUATIONS_SET_USER_NUMBER, &
  !  & CMFE_EQUATIONS_FULL_MATRICES,err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(contextUserNumber,REGION_USER_NUMBER,EQUATIONS_SET_USER_NUMBER, &
    & CMFE_EQUATIONS_NO_OUTPUT,err)
  !CALL cmfe_Equations_OutputTypeSet(contextUserNumber,REGION_USER_NUMBER,EQUATIONS_SET_USER_NUMBER, &
  !  & CMFE_EQUATIONS_TIMING_OUTPUT,err)
  !CALL cmfe_Equations_OutputTypeSet(contextUserNumber,REGION_USER_NUMBER,EQUATIONS_SET_USER_NUMBER, &
  !  & CMFE_EQUATIONS_MATRIX_OUTPUT,err)
  !CALL cmfe_Equations_OutputTypeSet(contextUserNumber,REGION_USER_NUMBER,EQUATIONS_SET_USER_NUMBER, &
  !  & CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(contextUserNumber,REGION_USER_NUMBER,EQUATIONS_SET_USER_NUMBER,err)

  !-----------------------------------------------------------------------------------------------------------
  ! PROBLEM
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a problem.
  CALL cmfe_Problem_CreateStart(PROBLEM_USER_NUMBER,contextUserNumber,[CMFE_PROBLEM_CLASSICAL_FIELD_CLASS, &
    & CMFE_PROBLEM_LAPLACE_EQUATION_TYPE,CMFE_PROBLEM_STANDARD_LAPLACE_SUBTYPE],err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(contextUserNumber,PROBLEM_USER_NUMBER,err)

  !-----------------------------------------------------------------------------------------------------------
  ! CONTROL LOOPS
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(contextUserNumber,PROBLEM_USER_NUMBER,err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(contextUserNumber,PROBLEM_USER_NUMBER,err)

  !-----------------------------------------------------------------------------------------------------------
  ! SOLVER
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the problem solvers
  CALL cmfe_Problem_SolversCreateStart(contextUserNumber,PROBLEM_USER_NUMBER,err)
  !CALL cmfe_Solver_OutputTypeSet(contextUserNumber,PROBLEM_USER_NUMBER,CMFE_CONTROL_LOOP_NODE,1, &
  !  & CMFE_SOLVER_NO_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(contextUserNumber,PROBLEM_USER_NUMBER,CMFE_CONTROL_LOOP_NODE,1, &
  !  & CMFE_SOLVER_PROGRESS_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(contextUserNumber,PROBLEM_USER_NUMBER,CMFE_CONTROL_LOOP_NODE,1, &
  !  & CMFE_SOLVER_TIMING_OUTPUT,err)
  !CALL cmfe_Solver_OutputTypeSet(contextUserNumber,PROBLEM_USER_NUMBER,CMFE_CONTROL_LOOP_NODE,1, &
  !  & CMFE_SOLVER_SOLVER_OUTPUT,err)
  CALL cmfe_Solver_OutputTypeSet(contextUserNumber,PROBLEM_USER_NUMBER,CMFE_CONTROL_LOOP_NODE,1, &
    & CMFE_SOLVER_MATRIX_OUTPUT,err)
  
  CALL cmfe_Solver_LinearTypeSet(contextUserNumber,PROBLEM_USER_NUMBER,CMFE_CONTROL_LOOP_NODE,1, &
    & CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,err)
  CALL cmfe_Solver_LinearIterativeAbsoluteToleranceSet(contextUserNumber,PROBLEM_USER_NUMBER,CMFE_CONTROL_LOOP_NODE,1, &
    & 1.0E-12_CMISSRP,err)
  CALL cmfe_Solver_LinearIterativeRelativeToleranceSet(contextUserNumber,PROBLEM_USER_NUMBER,CMFE_CONTROL_LOOP_NODE,1, &
    & 1.0E-12_CMISSRP,err)
  
  !CALL cmfe_Solver_LinearTypeSet(contextUserNumber,PROBLEM_USER_NUMBER,CMFE_CONTROL_LOOP_NODE,1, &
  !  & CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,err)
  !CALL cmfe_Solver_LibraryTypeSet(contextUserNumber,PROBLEM_USER_NUMBER,CMFE_CONTROL_LOOP_NODE,1, &
  !  & CMFE_SOLVER_MUMPS_LIBRARY,err)
  !CALL cmfe_Solver_LibraryTypeSet(contextUserNumber,PROBLEM_USER_NUMBER,CMFE_CONTROL_LOOP_NODE,1, &
  !  & CMFE_SOLVER_LAPACK_LIBRARY,err)
  !CALL cmfe_Solver_LibraryTypeSet(contextUserNumber,PROBLEM_USER_NUMBER,CMFE_CONTROL_LOOP_NODE,1, &
  !  & CMFE_SOLVER_SUPERLU_LIBRARY,err)
  !CALL cmfe_Solver_LibraryTypeSet(contextUserNumber,PROBLEM_USER_NUMBER,CMFE_CONTROL_LOOP_NODE,1, &
  !  & CMFE_SOLVER_PASTIX_LIBRARY,err)
  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(contextUserNumber,PROBLEM_USER_NUMBER,err)

  !-----------------------------------------------------------------------------------------------------------
  ! SOLVER EQUATIONS
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateStart(contextUserNumber,PROBLEM_USER_NUMBER,err)
  !Set the solver equations sparsity
  CALL cmfe_SolverEquations_SparsityTypeSet(contextUserNumber,PROBLEM_USER_NUMBER,CMFE_CONTROL_LOOP_NODE,1, &
    & CMFE_SOLVER_SPARSE_MATRICES,err)
  !CALL cmfe_SolverEquations_SparsityTypeSet(contextUserNumber,PROBLEM_USER_NUMBER,CMFE_CONTROL_LOOP_NODE,1, &
  !  & CMFE_SOLVER_FULL_MATRICES,err)
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(contextUserNumber,PROBLEM_USER_NUMBER,CMFE_CONTROL_LOOP_NODE,1, &
    & REGION_USER_NUMBER,EQUATIONS_SET_USER_NUMBER,equationsSetIndex,err)
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(contextUserNumber,PROBLEM_USER_NUMBER,err)

  !-----------------------------------------------------------------------------------------------------------
  ! BOUNDARY CONDITIONS
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the equations set boundary conditions
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(contextUserNumber,PROBLEM_USER_NUMBER,CMFE_CONTROL_LOOP_NODE,1,err)
  !Set the first node to 0.0 and the last node to 1.0
  firstNodeNumber=1
  CALL cmfe_Nodes_NumberOfNodesGet(contextUserNumber,REGION_USER_NUMBER,lastNodeNumber,err)
  CALL cmfe_Decomposition_NodeDomainGet(contextUserNumber,REGION_USER_NUMBER,MESH_USER_NUMBER,DECOMPOSITION_USER_NUMBER, &
    & firstNodeNumber,1,firstNodeDomain,err)
  CALL cmfe_Decomposition_NodeDomainGet(contextUserNumber,REGION_USER_NUMBER,MESH_USER_NUMBER,DECOMPOSITION_USER_NUMBER, &
    & lastNodeNumber,1,lastNodeDomain,err)
  IF(firstNodeDomain==computationalNodeNumber) THEN
    CALL cmfe_BoundaryConditions_SetNode(contextUserNumber,PROBLEM_USER_NUMBER,CMFE_CONTROL_LOOP_NODE,1, &
      & REGION_USER_NUMBER,DEPENDENT_FIELD_USER_NUMBER,CMFE_FIELD_U_VARIABLE_TYPE,1,1,firstNodeNumber,1, &
      & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,err)
  ENDIF
  IF(lastNodeDomain==computationalNodeNumber) THEN
    CALL cmfe_BoundaryConditions_SetNode(contextUserNumber,PROBLEM_USER_NUMBER,CMFE_CONTROL_LOOP_NODE,1, &
      & REGION_USER_NUMBER,DEPENDENT_FIELD_USER_NUMBER,CMFE_FIELD_U_VARIABLE_TYPE,1,1,lastNodeNumber,1, &
      & CMFE_BOUNDARY_CONDITION_FIXED,1.0_CMISSRP,err)
  ENDIF
  !Finish the creation of the equations set boundary conditions
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(contextUserNumber,PROBLEM_USER_NUMBER,CMFE_CONTROL_LOOP_NODE,1,err)

  !-----------------------------------------------------------------------------------------------------------
  ! SOLVE
  !-----------------------------------------------------------------------------------------------------------

  !Solve the problem
  CALL cmfe_Problem_Solve(contextUserNumber,PROBLEM_USER_NUMBER,err)

  !-----------------------------------------------------------------------------------------------------------
  ! OUTPUT
  !-----------------------------------------------------------------------------------------------------------

  !Export results
  CALL cmfe_Fields_NodesExport(contextUserNumber,REGION_USER_NUMBER,"LaplaceEquation","FORTRAN",err)
  CALL cmfe_Fields_ElementsExport(contextUserNumber,REGION_USER_NUMBER,"LaplaceEquation","FORTRAN",err)

  !Finialise OpenCMISS
  CALL cmfe_Finalise(contextUserNumber,err)
  WRITE(*,'(A)') "Program successfully completed."
  STOP

CONTAINS

  SUBROUTINE HandleError(errorString)
    CHARACTER(LEN=*), INTENT(IN) :: errorString
    WRITE(*,'(">>ERROR: ",A)') errorString(1:LEN_TRIM(errorString))
    STOP
  END SUBROUTINE HandleError

END PROGRAM LaplaceEquation
