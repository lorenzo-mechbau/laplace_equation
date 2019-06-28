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
  REAL(CMISSRP), PARAMETER :: WIDTH=2.0_CMISSRP
  REAL(CMISSRP), PARAMETER :: LENGTH=3.0_CMISSRP

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: LagrangeBasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: HermiteBasisUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=13
 
  !Program types
  
  !Program variables
  INTEGER(CMISSIntg) :: RowNo, ColumnNo, GlobalElementNo, GlobalNodeNo, I, J
  INTEGER(CMISSIntg) :: NumberOfElements, NumberOfNodes
  INTEGER(CMISSIntg) :: equationSparsity, solverSparsity
  INTEGER(CMISSIntg) :: numberOfArguments,argumentLength,status
  INTEGER(CMISSIntg) :: numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements, &
    & lagrangeInterpolationType,numberOfGaussXi
  CHARACTER(LEN=255) :: commandArgument,filename

  LOGICAL                               :: directory_exists = .FALSE.

  !CMISS variables
  TYPE(cmfe_BasisType) :: LagrangeBasis
  TYPE(cmfe_BasisType) :: HermiteBasis
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: Equations
  TYPE(cmfe_EquationsSetType) :: EquationsSet
  TYPE(cmfe_FieldType) :: GeometricField,EquationsSetField,DependentField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_GeneratedMeshType) :: generatedMesh
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
  
  LOGICAL              :: useGeneratedMesh ! FALSE: User-defined mesh
  INTEGER(CMISSIntg)   :: whichBasis 
                                   ! For MANUAL mesh, decide if all linear (2), all Hermite (3)
                                   ! or Hermite/linear (1) (cf. Benjamin's example)
                                   ! Generated mesh: only 2 and 3 (NO mix)!!!  
  INTEGER(CMISSIntg) :: decompositionType, numberOfComponents    ! Geometric Field components:
  ! my variables 
  INTEGER(CMISSIntg) :: count_components = 1
  INTEGER(CMISSIntg) :: domain, localNodesNumber, derivativeNumber
  INTEGER(CMISSIntg), ALLOCATABLE :: elementUserNodes(:), hermiteNodes(:)
  LOGICAL, ALLOCATABLE :: trackNodes(:)
  LOGICAL :: updateDerivativeForNode
  CHARACTER(LEN=60) :: diagFilename
  REAL(CMISSRP) :: offSetVer,offSetHor

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

  !-----------------------------------------------------------------------------------------------------------
  ! PROBLEM CONTROL PANEL
  !-----------------------------------------------------------------------------------------------------------

  ! No input arguments for now, but keep consistent with the example @develop. 
  numberOfArguments = COMMAND_ARGUMENT_COUNT()
  IF(numberOfArguments >= 8) THEN
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
    READ(commandArgument(1:argumentLength),*) lagrangeInterpolationType
    IF(lagrangeInterpolationType<=0) CALL HandleError("Invalid Interpolation specification.")
    CALL GET_COMMAND_ARGUMENT(5,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 5.")
    READ(commandArgument(1:argumentLength),*) whichBasis
    IF(whichBasis<=0) CALL HandleError("Invalid choice of basis.")
    CALL GET_COMMAND_ARGUMENT(6,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 6.")
    READ(commandArgument(1:argumentLength),*) equationSparsity
    IF(equationSparsity<=0 .OR. equationSparsity>=3) CALL HandleError("Invalid choice of eq. sparsity.")
    CALL GET_COMMAND_ARGUMENT(7,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 6.")
    READ(commandArgument(1:argumentLength),*) solverSparsity
    IF(solverSparsity<=0 .OR. solverSparsity>=3) CALL HandleError("Invalid choice of solver sparsity.")
    CALL GET_COMMAND_ARGUMENT(8,commandArgument,argumentLength,status)
    IF(status>0) CALL HandleError("Error for command argument 8.")
    READ(commandArgument(1:argumentLength),*) useGeneratedMesh
  ELSE

    !If there are not enough arguments default the problem specification
    numberOfGlobalXElements=4!4!2!2!2!
    numberOfGlobalYElements=4!4!1!2!1!
    numberOfGlobalZElements=0!0!1!2!1!
    ! interpolation type for the Lagrange basis
    lagrangeInterpolationType = 1!CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION ! i.e. 1
    whichBasis = 2                 ! For MANUAL mesh, decide if
                                   ! - all Lagrange linear/quad (2),
                                   ! - all Hermite (3)
                                   ! - mixed Hermite/Lagrange LINEAR (1)
                                   ! Generated mesh: only 2 and 3 (NO mix)!!!
    !If Hermite elements 3D we need LINES computation. Default 3D is false.
    !Set below to TRUE manually, will return an error msg if -np>1.
    !3D lines computation still missing!

    equationSparsity = CMFE_EQUATIONS_FULL_MATRICES ! Sparsity of the EQUATIONS matrix

    solverSparsity   = CMFE_SOLVER_SPARSE_MATRICES !CMFE_SOLVER_FULL_MATRICES
                                                   !Sparsity of the SOLVER matrix
    numberOfComponents = 2
    ! -2 Standard 2 node-based components (automatically increased to 3 if 3D).
    ! -3 add 1 element-based
    ! -4 add 1 constant-based
    ! Decomposition OK also with >2. cmfe_Field_CreateFinish must be tested!
    ! Solution can be computed ONLY with node-based components though!!!

    useGeneratedMesh = .TRUE. ! Manual mesh (.FALSE.) only for 2D case, see picture in docs/.

    ! User-defined or automatic decomposition
    decompositionType = CMFE_DECOMPOSITION_USER_DEFINED_TYPE!CMFE_DECOMPOSITION_CALCULATED_TYPE!!

    SetDecompositionDistributed = .TRUE. !.FALSE.! Set to TRUE ONLY for new implementation (FV_merged).
                                         ! True only works with new implementation, false works with both.
                                         ! Affects user-defined decomposition for -np 4 nodes (see below).
  END IF


  !Intialise OpenCMISS
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)

  !CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_OUTPUT_ERROR,Err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)

  CALL cmfe_RandomSeedsSet(9999,Err)
  
  diagFilename = "Diagnostics" !Develop"

  CALL cmfe_DiagnosticsSetOn(CMFE_IN_DIAG_TYPE,[1,2,3],TRIM(diagFileName),&
    & ["DOMAIN_MAPPINGS_NODES_CALCULATE               ", &
    &  "DOMAIN_MAPPINGS_INITIALISE                    ", &
    &  "MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH          ", &
    &  "DECOMPOSITION_TOPOLOGY_LINES_CALCULATE        ", &
    &  "FIELD_MAPPINGS_CALCULATE                      ", &
    &  "SOLVER_MATRIX_STRUCTURE_CALCULATE             ", &
    &  "SOLVER_MAPPING_CALCULATE                      ", &
    &  "DistributedMatrix_StorageLocationsSet         ", &
    &  "DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE          ", &
    &  "FIELD_SCALINGS_CALCULATE                      ", &
    &  "Field_GeometricParametersLineLengthsCalculate ", &
    &  "DomainMappings_2DLinesCalculate               ", &
    &  "DomainMappings_3DLinesCalculate               ", &
    &  "MeshTopology_ElementsAdjacentElementsCalculate", &
    &  "DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE   ", &
    &  "DECOMPOSITION_ELEMENTS_MAPPING_MAP            ", &
    &  "DECOMPOSITION_ADJACENT_DOMAINS_CALCULATE      ", &
    &  "DOMAIN_MAPPINGS_ELEMENTS_CALCULATE            ", &
    &  "DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE        ", &
    &  "DECOMPOSITION_ELEMENTS_MAPPING_MAP_CALCULATE  ", &
    &  "DECOMPOSITION_ELEMENTS_MAPPING_CALCULATE      "], Err)

  !CALL cmfe_OutputSetOn("diagnostics.txt",Err)

  !Get the computational nodes information
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  IF (NumberOfComputationalNodes > 4) THEN
    CALL HandleError("Max. 4 computational nodes!")
  END IF

  !-----------------------------------------------------------------------------------------------------------
  !COORDINATE SYSTEM
  !-----------------------------------------------------------------------------------------------------------  
 
  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  IF(numberOfGlobalZElements==0) THEN
    !Set the coordinate system to be 2D
    CALL cmfe_CoordinateSystem_DimensionSet(coordinateSystem,2,err)
  ELSE
    !Set the coordinate system to be 3D
    CALL cmfe_CoordinateSystem_DimensionSet(coordinateSystem,3,err)
  END IF

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
  CALL cmfe_Basis_Initialise(LagrangeBasis,Err)
  CALL cmfe_Basis_CreateStart(LagrangeBasisUserNumber,LagrangeBasis,Err)

  SELECT CASE(lagrangeInterpolationType)
  CASE(1,2,3,4)
    CALL cmfe_Basis_TypeSet(LagrangeBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
  CASE(7,8,9)
    CALL cmfe_Basis_TypeSet(LagrangeBasis,CMFE_BASIS_SIMPLEX_TYPE,err)
    CALL HandleError("Simplex basis should be checked for this example first!")
  CASE DEFAULT
    CALL HandleError("Invalid interpolation type.")
  END SELECT
  SELECT CASE(lagrangeInterpolationType)
  CASE(1) ! linear Lagrange
    numberOfGaussXi=2
  CASE(2) ! quadratic Lagrange
    numberOfGaussXi=3
  CASE(3,4) ! cubic Hermite/Lagrange
    numberOfGaussXi=4
  CASE DEFAULT
    numberOfGaussXi=0 !Don't set number of Gauss points for tri/tet
  END SELECT

  IF(numberOfGlobalZElements==0) THEN
    !Set the basis to be a bi-interpolation basis
    CALL cmfe_Basis_NumberOfXiSet(LagrangeBasis,2,err)
    CALL cmfe_Basis_InterpolationXiSet(LagrangeBasis,[lagrangeInterpolationType,lagrangeInterpolationType],err)
    IF(numberOfGaussXi>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(LagrangeBasis,[numberOfGaussXi,numberOfGaussXi],err)
    END IF
  ELSE
    !Set the basis to be a tri-interpolation basis
    CALL cmfe_Basis_NumberOfXiSet(LagrangeBasis,3,err)
    CALL cmfe_Basis_InterpolationXiSet(LagrangeBasis, &
      & [lagrangeInterpolationType,lagrangeInterpolationType,lagrangeInterpolationType],err)
    IF(numberOfGaussXi>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(LagrangeBasis,[numberOfGaussXi,numberOfGaussXi,numberOfGaussXi],err)
    END IF
  END IF

  !CALL cmfe_Basis_NumberOfXiSet(LagrangeBasis,2,Err)
  !CALL cmfe_Basis_InterpolationXiSet(LagrangeBasis,&
  !  & [CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION,CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err) !i.e. 1
  !CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(LagrangeBasis,[2,2],Err)
  
  CALL cmfe_Basis_CreateFinish(LagrangeBasis,Err)
  
  !Create hermite basis
  CALL cmfe_Basis_Initialise(HermiteBasis,Err)
  CALL cmfe_Basis_CreateStart(HermiteBasisUserNumber,HermiteBasis,Err)
  CALL cmfe_Basis_TypeSet(HermiteBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)

  IF(numberOfGlobalZElements==0) THEN
    !Set the basis to be a bi-interpolation basis
    CALL cmfe_Basis_NumberOfXiSet(HermiteBasis,2,err)
    CALL cmfe_Basis_InterpolationXiSet(HermiteBasis, & 
      & [CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION,CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION],err)
    IF(numberOfGaussXi>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(HermiteBasis,[4,4],err)
    END IF
  ELSE
    !CALL HandleError("No Hermite basis defined for 3D case.")
    !Set the basis to be a tri-interpolation basis
    CALL cmfe_Basis_NumberOfXiSet(HermiteBasis,3,err)
    CALL cmfe_Basis_InterpolationXiSet(HermiteBasis,[CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION, &
      & CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION,CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION],err)
    IF(numberOfGaussXi>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(HermiteBasis,[4,4,4],err)
    END IF
  END IF
 
  CALL cmfe_Basis_CreateFinish(HermiteBasis,Err)

  !-----------------------------------------------------------------------------------------------------------
  !MESH
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a generated mesh in the region
  IF (useGeneratedMesh) THEN
    WRITE(*,*) "Using a generated mesh"
    CALL cmfe_GeneratedMesh_Initialise(generatedMesh,err)
    CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,generatedMesh,err)
    !Set up a regular x*y*z mesh
    CALL cmfe_GeneratedMesh_TypeSet(generatedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,err)
    !Set the default basis
    !Using different basis types is not supported for generated meshes! ==> ALL LINEAR/HERMITE!!!
    SELECT CASE (whichBasis)
    CASE (2)
      CALL cmfe_GeneratedMesh_BasisSet(generatedMesh,LagrangeBasis,err) 
    CASE (3)
      CALL cmfe_GeneratedMesh_BasisSet(generatedMesh,HermiteBasis,err) 
    CASE DEFAULT
      CALL HandleError("Mixed basis not possible for generated mesh!")
    END SELECT

    !Define the mesh on the region
    IF(numberOfGlobalZElements==0) THEN
      CALL cmfe_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,HEIGHT],err)
      CALL cmfe_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements,numberOfGlobalYElements],err)
      NumberOfElements = numberOfGlobalXElements*numberOfGlobalYElements
    ELSE
      CALL cmfe_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,HEIGHT,LENGTH],err)
      CALL cmfe_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements,numberOfGlobalYElements, &
        & numberOfGlobalZElements],err)
      NumberOfElements = numberOfGlobalXElements*numberOfGlobalYElements*numberOfGlobalZElements
    END IF
    !Finish the creation of a generated mesh in the region
    CALL cmfe_Mesh_Initialise(Mesh,err)
    CALL cmfe_GeneratedMesh_CreateFinish(generatedMesh,MeshUserNumber,Mesh,err)

    ! Get the number of nodes
    CALL cmfe_Nodes_Initialise(Nodes,Err)
    CALL cmfe_Region_NodesGet(Region,Nodes,Err)
    CALL cmfe_Nodes_NumberOfNodesGet(Nodes,NumberOfNodes,Err)

  ELSE
  ! Create MANUAL mesh as specified in "PartitionedMesh_4x4.pdf"
  
    CALL cmfe_Mesh_Initialise(Mesh,Err)
    IF (numberOfGlobalZElements == 0) THEN
      CALL cmfe_Mesh_CreateStart(MeshUserNumber,Region,2,Mesh,Err) ! 2 is number of dim
    ELSE
      CALL cmfe_Mesh_CreateStart(MeshUserNumber,Region,3,Mesh,Err) ! 3 is number of dim
      CALL HandleError("User-defined mesh only 2D!")
    END IF
    CALL cmfe_Mesh_NumberOfComponentsSet(Mesh,1,Err)

    NumberOfElements = 16
    ! Special case 3 nodes
    IF (numberOfGlobalXElements == 3) NumberOfElements = 3

    CALL cmfe_Mesh_NumberOfElementsSet(Mesh,NumberOfElements,Err)
  
    ! Create elements 
    CALL cmfe_MeshElements_Initialise(MeshElements,Err)
    CALL cmfe_MeshElements_CreateStart(Mesh,1,LagrangeBasis,MeshElements,Err)

    PRINT *, "Basis set start!"
    ! Set the basis
    SELECT CASE (whichBasis)
    CASE(1) ! mix linear-Hermite  
      SELECT CASE (numberOfElements) 
      CASE (16)
        DO GlobalElementNo=1,4
          CALL cmfe_MeshElements_BasisSet(MeshElements,GlobalElementNo,HermiteBasis,Err)
        END DO
        DO GlobalElementNo=5,12
          CALL cmfe_MeshElements_BasisSet(MeshElements,GlobalElementNo,LagrangeBasis,Err)
        END DO
        DO GlobalElementNo=13,16
          CALL cmfe_MeshElements_BasisSet(MeshElements,GlobalElementNo,HermiteBasis,Err)
        END DO
      CASE (3)
        CALL cmfe_MeshElements_BasisSet(MeshElements,1,LagrangeBasis,Err)
        CALL cmfe_MeshElements_BasisSet(MeshElements,2,HermiteBasis,Err)
        CALL cmfe_MeshElements_BasisSet(MeshElements,3,LagrangeBasis,Err)
      CASE DEFAULT
        CALL HandleError("Number of elements not supported for user-defined mesh!!!") 
      END SELECT        
    CASE(2) ! all linear
      DO GlobalElementNo=1,NumberOfElements
        CALL cmfe_MeshElements_BasisSet(MeshElements,GlobalElementNo,LagrangeBasis,Err)
      END DO
    CASE(3) ! all Hermite
      DO GlobalElementNo=1,NumberOfElements
        CALL cmfe_MeshElements_BasisSet(MeshElements,GlobalElementNo,HermiteBasis,Err)
      END DO
    CASE DEFAULT
      CALL HandleError("Basis choice not supported!!!")
    END SELECT

    PRINT *, "Basis set finish!"

    ! Define nodes for the mesh
    CALL cmfe_Nodes_Initialise(Nodes,Err)

    IF (lagrangeInterpolationType==1 .OR. whichBasis==3) THEN ! linear Lagrange or cubic Hermite
      NumberOfNodes = (numberOfGlobalXElements+1)*(numberOfGlobalXElements+1)
    ELSEIF (lagrangeInterpolationType==2) THEN
      NumberOfNodes = (numberOfGlobalXElements*2+1)*(numberOfGlobalXElements*2+1)
    END IF
    ! Special case 3 nodes
    IF (numberOfGlobalXElements == 3) NumberOfNodes = 8

    CALL cmfe_Nodes_CreateStart(Region,NumberOfNodes,Nodes,Err)
    CALL cmfe_Nodes_CreateFinish(Nodes,Err)

  ! Special case 3 nodes
  ! Set adjacent nodes for each element
    IF (numberOfGlobalXElements == 3) THEN
      CALL cmfe_MeshElements_NodesSet(MeshElements, 1, [1,2,5,6], Err)
      CALL cmfe_MeshElements_NodesSet(MeshElements, 2, [2,3,6,7], Err)
      CALL cmfe_MeshElements_NodesSet(MeshElements, 3, [3,4,7,8], Err)
    ELSE

    ! Distribute nodes among elements for a SQUARE(cf. picture)
      DO RowNo = 1,numberOfGlobalXElements
        DO ColumnNo = 1,numberOfGlobalXElements
          GlobalElementNo = (RowNo-1)*numberOfGlobalXElements + ColumnNo
          IF (lagrangeInterpolationType==1 .OR. whichBasis==3) THEN     !linear Lagrange or cubic Hermite
            GlobalNodeNo = (RowNo-1)*5 + ColumnNo
            CALL cmfe_MeshElements_NodesSet(MeshElements, GlobalElementNo, &
              & [GlobalNodeNo,  GlobalNodeNo+1, &
              &  GlobalNodeNo+5,GlobalNodeNo+6], Err)
          ELSE IF (lagrangeInterpolationType==2) THEN
            GlobalNodeNo = (RowNo-1)*18+(ColumnNo*2-1)
            CALL cmfe_MeshElements_NodesSet(MeshElements, GlobalElementNo, &
              & [GlobalNodeNo,   GlobalNodeNo+1, GlobalNodeNo+2,  &
              &  GlobalNodeNo+9, GlobalNodeNo+10,GlobalNodeNo+11, &
              &  GlobalNodeNo+18,GlobalNodeNo+19,GlobalNodeNo+20], Err)
          END IF
        END DO
      END DO
    END IF ! number of comp. nodes (3 vs. other)
  
    CALL cmfe_MeshElements_CreateFinish(MeshElements,Err)
    PRINT *, "MeshElements create finish!"
    CALL cmfe_Mesh_CreateFinish(Mesh,Err)
    PRINT *, "Mesh create finish!"

    ! Collect Hermite nodes to update geometric field in case of a mixed basis
    IF (whichBasis==1) THEN
      ALLOCATE(hermiteNodes(0), stat=err)
      IF(err/=0) CALL HandleError("Could not allocate array of Hermite elements.")
      ALLOCATE(elementUserNodes(4),STAT=Err)
      IF(ERR/=0) CALL HandleError("Could not allocate nodes on element array.")
      SELECT CASE (numberOfElements) 
      CASE (16)
        DO GlobalElementNo=1,4
          CALL cmfe_MeshElements_NodesGet(MeshElements,GlobalElementNo,elementUserNodes,Err)
          hermiteNodes = [hermiteNodes,elementUserNodes]
        END DO
        DO GlobalElementNo=13,16
          CALL cmfe_MeshElements_NodesGet(MeshElements,GlobalElementNo,elementUserNodes,Err)
          hermiteNodes = [hermiteNodes,elementUserNodes]
        END DO
      CASE (3)
        CALL cmfe_MeshElements_NodesGet(MeshElements,2,elementUserNodes,Err)
        hermiteNodes = [hermiteNodes,elementUserNodes]
      CASE DEFAULT
        CALL HandleError("Number of elements not supported for user-defined mesh!!!") 
      END SELECT        
    END IF ! mixed basis

  END IF ! generated mesh or manual mesh

  ! Check the nodes distribution among elements:
  ! Use diagnostic MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH!
  ! global element numbers
  !13 14 15 16
  ! 9 10 11 12
  ! 5  6  7  8
  ! 1  2  3  4
  
  ! OR (-np 3)
  ! 1 2 3
  
  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  ! Set the decomposition to be a general decomposition with the specified number of domains

  !If Hermite elements 3D we need LINES computation. Default 3D is false.
  !Fixed for now only for -np 1!
  !IF(whichBasis/=2)?
  !IF (numberOfGlobalZElements >=1 .AND. NumberOfComputationalNodes==1) &
  !Compute lines also in 3D (default false) - there will be an error msg if -np >1:
  CALL cmfe_Decomposition_CalculateLinesSet(Decomposition, .TRUE., Err)

  CALL cmfe_Decomposition_TypeSet(decomposition,decompositionType,err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
   
  IF(decompositionType==CMFE_DECOMPOSITION_USER_DEFINED_TYPE) THEN
    IF (numberOfGlobalZElements >=1) THEN
     ! 3D, max -np 2 for 2 elements
      SELECT CASE (NumberOfComputationalNodes)
      CASE (1)
        DO I=1,NumberOfElements
        ! All nodes to the same rank
          CALL cmfe_Decomposition_ElementDomainSet(Decomposition,I,0,Err)
        END DO
      CASE (2)
        IF (NumberOfElements ==2) THEN
          CALL cmfe_Decomposition_ElementDomainSet(Decomposition,1,0,Err)
          CALL cmfe_Decomposition_ElementDomainSet(Decomposition,2,1,Err)
        ELSE
          CALL HandleError("3D case 2 ranks requires 2 elements!!!")
        END IF
      CASE(4)
!        IF (NumberOfElements ==8) THEN
!          CALL cmfe_Decomposition_ElementDomainSet(Decomposition,1,0,Err)
!          CALL cmfe_Decomposition_ElementDomainSet(Decomposition,2,0,Err)
!          CALL cmfe_Decomposition_ElementDomainSet(Decomposition,3,1,Err)
!          CALL cmfe_Decomposition_ElementDomainSet(Decomposition,4,1,Err)
!          CALL cmfe_Decomposition_ElementDomainSet(Decomposition,5,2,Err)
!          CALL cmfe_Decomposition_ElementDomainSet(Decomposition,6,2,Err)
!          CALL cmfe_Decomposition_ElementDomainSet(Decomposition,7,3,Err)
!          CALL cmfe_Decomposition_ElementDomainSet(Decomposition,8,3,Err)
!        ELSE
          CALL HandleError("3D case 4 ranks requires calculated decomposition!!!")
!        END IF
      CASE DEFAULT
        CALL HandleError("3D case supported for 2 or 4 ranks!!!")
      END SELECT
    ELSE
    ! 2D case
      SELECT CASE (NumberOfComputationalNodes)
      CASE (1)
        DO I=1,NumberOfElements
        ! All nodes to the same rank
          CALL cmfe_Decomposition_ElementDomainSet(Decomposition,I,0,Err)
        END DO
      CASE (2)
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
      CASE (3)
        CALL cmfe_Decomposition_ElementDomainSet(Decomposition,1,0,Err)
        CALL cmfe_Decomposition_ElementDomainSet(Decomposition,2,1,Err)
        CALL cmfe_Decomposition_ElementDomainSet(Decomposition,3,2,Err)
      CASE (4)
        IF (ComputationalNodeNumber == 0 .OR..NOT.SetDecompositionDistributed) THEN
      !                                             global el., domain
          CALL cmfe_Decomposition_ElementDomainSet(Decomposition,9,2,Err)
          CALL cmfe_Decomposition_ElementDomainSet(Decomposition,5,0,Err)
          CALL cmfe_Decomposition_ElementDomainSet(Decomposition,14,2,Err)
          CALL cmfe_Decomposition_ElementDomainSet(Decomposition,11,3,Err)
          CALL cmfe_Decomposition_ElementDomainSet(Decomposition,6,0,Err)
          CALL cmfe_Decomposition_ElementDomainSet(Decomposition,16,3,Err)
        END IF
    
        IF (ComputationalNodeNumber == 1 .OR..NOT.SetDecompositionDistributed) THEN
          CALL cmfe_Decomposition_ElementDomainSet(Decomposition,7,0,Err)
          CALL cmfe_Decomposition_ElementDomainSet(Decomposition,15,3,Err)
          CALL cmfe_Decomposition_ElementDomainSet(Decomposition,1,0,Err)
        END IF
    
        IF (ComputationalNodeNumber == 2 .OR..NOT.SetDecompositionDistributed) THEN
          CALL cmfe_Decomposition_ElementDomainSet(Decomposition,10,2,Err)
          CALL cmfe_Decomposition_ElementDomainSet(Decomposition,8,1,Err)
        END IF

        IF (ComputationalNodeNumber == 3 .OR..NOT.SetDecompositionDistributed) THEN
          CALL cmfe_Decomposition_ElementDomainSet(Decomposition,2,0,Err)
          CALL cmfe_Decomposition_ElementDomainSet(Decomposition,3,1,Err)
          CALL cmfe_Decomposition_ElementDomainSet(Decomposition,12,3,Err)
          CALL cmfe_Decomposition_ElementDomainSet(Decomposition,4,1,Err)
          CALL cmfe_Decomposition_ElementDomainSet(Decomposition,13,2,Err)
        END IF
      CASE DEFAULT
        CALL HandleError("Number of nodes NOT supported (max 4)!!!")
      END SELECT
    END IF ! number of dimensions
  END IF ! decomposition type
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)
 
  !Destroy the mesh now that we have decomposed it
  !CALL cmfe_Mesh_Destroy(Mesh,Err)
  !Variable, MaxDepth, MaxArrayLength
  !CALL cmfe_PrintDecomposition(Decomposition,3,100,Err)

  !Possibility to test cmfe_Decomposition_ElementDomainGet for new implementation.
  !Returns -1 if element is not local.
  IF (.FALSE.) THEN
    DO GlobalElementNo = 1,NumberOfElements
      CALL cmfe_Decomposition_ElementDomainGet(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,&
        & GlobalElementNo, domain, err)
        ! Only LOCAL elements can be updated
      IF (ComputationalNodeNumber == domain) THEN
        PRINT *, "Rank that owns element:"
        PRINT *, domain, GlobalElementNo
      END IF
    END DO
  END IF
 
  !-----------------------------------------------------------------------------------------------------------
  !GEOMETRIC FIELD
  !-----------------------------------------------------------------------------------------------------------


  !Start to create a default (geometric) field on the region
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the domain to be used by the field components.
 
  !Components of the geometric field:
  !ndim components nodal-based
  !numberOfComponents = 2
  !Option: add 1 element-based +1 constant but NO solve if these components added! -> Moved above.
  IF(numberOfGlobalZElements/=0) numberOfComponents = numberOfComponents+1
  ! Just for decomposition test purposes!!

  !1 vbl, components: 2(3)x Mesh (node-dofs), 1 x Element-based, 1 x constant 
  CALL cmfe_Field_NumberOfComponentsSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,numberOfComponents,Err)
 
  DO count_components=1,2
    CALL cmfe_Field_ComponentMeshComponentSet &
      & (GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,count_components,1,Err)
  END DO
  count_components=count_components-1

  IF(numberOfGlobalZElements/=0) THEN
    count_components = count_components +1
    CALL cmfe_Field_ComponentMeshComponentSet &
      & (geometricField,CMFE_FIELD_U_VARIABLE_TYPE,count_components,1,err)
  END IF

  ! other interpolations
  IF (count_components < numberOfComponents) THEN
    count_components = count_components +1
    CALL cmfe_Field_ComponentInterpolationSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,count_components, &
      & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
  END IF
  IF (count_components < numberOfComponents) THEN
    count_components = count_components +1
    CALL cmfe_Field_ComponentInterpolationSet &  
      & (GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,count_components,CMFE_FIELD_CONSTANT_INTERPOLATION,Err)
  END IF

  !Finish creating the field
  CALL cmfe_Field_CreateFinish(GeometricField,Err)

  WRITE (*,*) "Created geometric field components"
  WRITE (*,*) count_components

  PRINT *, "Ready to update geometric field parameters!"
  !PRINT *, "Print field GeometricField"
  !CALL cmfe_PrintField(GeometricField,3,5,Err)
  
  !PRINT *, "Abort program in laplace_equation.f90:360"
  !STOP
  
  ! Update the geometric field parameters for a parallel simulation
  IF (useGeneratedMesh) THEN
    CALL cmfe_GeneratedMesh_GeometricParametersCalculate(generatedMesh,GeometricField,Err)
    WRITE(*,*) "Update completed."
  ELSE

    ! Update the geometric field parameters manually if mesh is NOT a generated mesh
    ! cf. procedure in GeneratedMesh_RegularGeometricParametersCalculate:
    ! update nodes + update start + update finish

    ! ELEMENTS
    ! Makes no sense for a geometric field in FE!
    ! WHAT is the value to be updated? Crashes in any case at updatefinish below!
    ! Same for a constant-interpolation.
    IF ((numberOfComponents >= 3 .AND. numberOfGlobalZElements==0) .OR. &
      & (numberOfComponents >= 4 .AND. numberOfGlobalZElements/=0)) THEN
    
      CALL HandleError("No meaningful update for element-based interpolation!")

      DO GlobalElementNo = 1,NumberOfElements 
        CALL cmfe_Decomposition_ElementDomainGet(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,&
          & GlobalElementNo, domain, err)
        ! Only LOCAL elements can be updated
        IF (ComputationalNodeNumber == domain) THEN
         !SUBROUTINE cmfe_Field_ParameterSetUpdateElementIntgObj (field,variableType,fieldSetType,userElementNumber,componentNumber, &
         !  & value,err)
        PRINT *, "Updating element number in rank"
        PRINT *, domain, GlobalElementNo
        CALL cmfe_Field_ParameterSetUpdateElement(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
            & CMFE_FIELD_VALUES_SET_TYPE,GlobalElementNo,3,0.0_CMISSRP,Err)
        END IF
      END DO
    END IF 

    ! NODES
    ! For Hermite basis: Tangent vector must have unit magnitude, update as follows:
    ! Component 1: dx/ds1 =1 dx/ds2 =0 d2x/ds1ds2 =0 
    ! Component 2: dy/ds1 =0 dy/ds2 =1 d2y/ds1ds2 =0
    ! For linear basis and no-derivative, update position of nodes. 
    ! CMFE_NO_GLOBAL_DERIV = 1

    ! Case quadratic Lagrange (ALL elements!!!)
    IF (lagrangeInterpolationType==2 .AND. whichBasis==2) THEN !Basis is all quadratic Lagrange
      derivativeNumber=1
      DO RowNo = 1,numberOfGlobalXElements
        offSetVer = HEIGHT/4.0_CMISSRP*(REAL(RowNo,CMISSRP)-1)
        DO ColumnNo = 1,numberOfGlobalXElements
          offSetHor = WIDTH/4.0_CMISSRP*(REAL(ColumnNo,CMISSRP)-1)
          GlobalNodeNo = (RowNo-1)*18+(ColumnNo*2-1)
     ! To which rank this node belongs
          CALL cmfe_Decomposition_NodeDomainGet(RegionUserNumber,MeshUserNumber,DecompositionUserNumber, &
            & GlobalNodeNo,1,domain,Err)
     ! Only LOCAL nodes can be updated
          IF (ComputationalNodeNumber == domain) THEN
            CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
              & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo,2, &
              & offSetVer,Err)
            CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
              & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo,1, &
              & offSetHor,Err)
          END IF
          CALL cmfe_Decomposition_NodeDomainGet(RegionUserNumber,MeshUserNumber,DecompositionUserNumber, &
            & GlobalNodeNo+9,1,domain,Err)
     ! Only LOCAL nodes can be updated
          IF (ComputationalNodeNumber == domain) THEN
            CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
              & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo+9,2, &
              & offSetVer+HEIGHT/8.0_CMISSRP,Err)
            CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
              & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo+9,1, &
              & offSetHor,Err)
          END IF
          CALL cmfe_Decomposition_NodeDomainGet(RegionUserNumber,MeshUserNumber,DecompositionUserNumber, &
            & GlobalNodeNo+18,1,domain,Err)
     ! Only LOCAL nodes can be updated
          IF (ComputationalNodeNumber == domain) THEN
            CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
              & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo+18,2, &
              & offSetVer+2.0_CMISSRP*HEIGHT/8.0_CMISSRP,Err)
            CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
              & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo+18,1, &
              & offSetHor,Err)
          END IF
          CALL cmfe_Decomposition_NodeDomainGet(RegionUserNumber,MeshUserNumber,DecompositionUserNumber, &
            & GlobalNodeNo+1,1,domain,Err)
     ! Only LOCAL nodes can be updated
          IF (ComputationalNodeNumber == domain) THEN
            CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
              & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo+1,2, &
              & offSetVer,Err)
            CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
              & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo+1,1, &
              & offSetHor+WIDTH/8.0_CMISSRP,Err)
          END IF
          CALL cmfe_Decomposition_NodeDomainGet(RegionUserNumber,MeshUserNumber,DecompositionUserNumber, &
            & GlobalNodeNo+10,1,domain,Err)
     ! Only LOCAL nodes can be updated
          IF (ComputationalNodeNumber == domain) THEN
            CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
              & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo+10,2, &
              & offSetVer+HEIGHT/8.0_CMISSRP,Err)
            CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
              & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo+10,1, &
              & offSetHor+WIDTH/8.0_CMISSRP,Err)
          END IF
          CALL cmfe_Decomposition_NodeDomainGet(RegionUserNumber,MeshUserNumber,DecompositionUserNumber, &
            & GlobalNodeNo+19,1,domain,Err)
     ! Only LOCAL nodes can be updated
          IF (ComputationalNodeNumber == domain) THEN
            CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
              & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo+19,2, &
              & offSetVer+2.0_CMISSRP*HEIGHT/8.0_CMISSRP,Err)
            CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
              & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo+19,1, &
              & offSetHor+WIDTH/8.0_CMISSRP,Err)
          END IF
          CALL cmfe_Decomposition_NodeDomainGet(RegionUserNumber,MeshUserNumber,DecompositionUserNumber, &
            & GlobalNodeNo+2,1,domain,Err)
     ! Only LOCAL nodes can be updated
          IF (ComputationalNodeNumber == domain) THEN
            CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
              & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo+2,2, &
              & offSetVer,Err)
            CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
              & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo+2,1, &
              & offSetHor+2.0_CMISSRP*WIDTH/8.0_CMISSRP,Err)
          END IF
          CALL cmfe_Decomposition_NodeDomainGet(RegionUserNumber,MeshUserNumber,DecompositionUserNumber, &
            & GlobalNodeNo+11,1,domain,Err)
     ! Only LOCAL nodes can be updated
          IF (ComputationalNodeNumber == domain) THEN
            CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
              & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo+11,2, &
              & offSetVer+HEIGHT/8.0_CMISSRP,Err)
            CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
              & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo+11,1, &
              & offSetHor+2.0_CMISSRP*WIDTH/8.0_CMISSRP,Err)
          END IF
          CALL cmfe_Decomposition_NodeDomainGet(RegionUserNumber,MeshUserNumber,DecompositionUserNumber, &
            & GlobalNodeNo+20,1,domain,Err)
     ! Only LOCAL nodes can be updated
          IF (ComputationalNodeNumber == domain) THEN
            CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
              & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo+20,2, &
              & offSetVer+2.0_CMISSRP*HEIGHT/8.0_CMISSRP,Err)
            CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
              & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo+20,1, &
              & offSetHor+2.0_CMISSRP*WIDTH/8.0_CMISSRP,Err)
          END IF
            !CALL cmfe_MeshElements_NodesSet(MeshElements, GlobalElementNo, &
            !  & [GlobalNodeNo,   GlobalNodeNo+1, GlobalNodeNo+2,  &
            !  &  GlobalNodeNo+9, GlobalNodeNo+10,GlobalNodeNo+11, &
            !  &  GlobalNodeNo+18,GlobalNodeNo+19,GlobalNodeNo+20], Err)
        END DO
      END DO
    ! IF linear lagrange or Hermite or mixed
    ELSE

      DO GlobalNodeNo = 1,NumberOfNodes

     ! To which rank this node belongs
        CALL cmfe_Decomposition_NodeDomainGet(RegionUserNumber,MeshUserNumber,DecompositionUserNumber, &
          & GlobalNodeNo,1,domain,Err)

     ! Only LOCAL nodes can be updated
        IF (ComputationalNodeNumber == domain) THEN

          derivativeNumber=1
          SELECT CASE (numberOfElements)
          CASE (16)

            SELECT CASE (GlobalNodeNo)
            CASE (1,6,11,16,21)
    ! node 1
    !SUBROUTINE FIELD_PARAMETER_SET_UPDATE_NODE_INTG(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,VERSION_NUMBER,DERIVATIVE_NUMBER, &
    !  & USER_NODE_NUMBER,COMPONENT_NUMBER,VALUE,ERR,ERROR,*)
              CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
                & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo,1,0.0_CMISSRP,Err)
            CASE (2,7,12,17,22)
              CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
                & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo,1,WIDTH/4.0_CMISSRP,Err)
            CASE (3,8,13,18,23)
              CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
                & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo,1,WIDTH/2.0_CMISSRP,Err)
            CASE (4,9,14,19,24)
              CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
                & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo,1,3.0_CMISSRP*WIDTH/4.0_CMISSRP,Err)
            CASE (5,10,15,20,25)
              CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
                & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo,1,WIDTH,Err)
            CASE DEFAULT
              CALL HandleError("Node out of reach!")
            END SELECT

            SELECT CASE (GlobalNodeNo)
            CASE (1,2,3,4,5)
              CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
                & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo,2,0.0_CMISSRP,Err)
            CASE (6,7,8,9,10)
              CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
                & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo,2,HEIGHT/4.0_CMISSRP,Err)
            CASE (11,12,13,14,15)
               CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
                 & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo,2,HEIGHT/2.0_CMISSRP,Err)
            CASE (16,17,18,19,20)
              CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
                & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo,2,3.0_CMISSRP*HEIGHT/4.0_CMISSRP,Err)
            CASE (21,22,23,24,25)
              CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
                & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo,2,HEIGHT,Err)
            CASE DEFAULT
              CALL HandleError("Node out of reach!")
            END SELECT

          CASE(3) ! special case 3 nodes

            SELECT CASE (GlobalNodeNo)
            CASE (1,5)
              CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
                & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo,1,0.0_CMISSRP,Err)
            CASE (2,6)
              CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
                & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo,1,WIDTH/3.0_CMISSRP,Err)
            CASE (3,7)
              CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
                & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo,1,2.0_CMISSRP*WIDTH/3.0_CMISSRP,Err)
            CASE (4,8)
              CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
                & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo,1,WIDTH,Err)
            CASE DEFAULT
              CALL HandleError("Node out of reach!")
            END SELECT

            SELECT CASE (GlobalNodeNo)
            CASE (1,2,3,4)
              CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
                & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo,2,0.0_CMISSRP,Err)
            CASE (5,6,7,8)
              CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
                & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo,2,HEIGHT,Err)
            CASE DEFAULT
              CALL HandleError("Node out of reach!")
            END SELECT

          CASE DEFAULT
            CALL HandleError("Number of elements not supported!")
          END SELECT ! number of elements
        !END DO ! derivatives

          updateDerivativeForNode = .FALSE.

          SELECT CASE (whichBasis)
          CASE(1) ! MIX: look up if node is a Hermite node
            DO I = 1,size(hermiteNodes)
              IF (GlobalNodeNo == hermiteNodes(I)) updateDerivativeForNode = .TRUE.
            END DO
          CASE(2) ! all linear: do not set any derivative
            updateDerivativeForNode = .FALSE.
          CASE(3) ! all Hermite: update all derivatives for all nodes
            updateDerivativeForNode = .TRUE.
          CASE DEFAULT
            CALL HandleError("Choice of basis incorrect!")
          END SELECT

          IF (updateDerivativeForNode) THEN
            derivativeNumber = 2 ! x derivative
            CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
              & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo,1,1.0_CMISSRP,Err)
            CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
              & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo,2,0.0_CMISSRP,Err)
            derivativeNumber = 3 ! y derivative
            CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
              & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo,1,0.0_CMISSRP,Err)
            CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
              & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo,2,1.0_CMISSRP,Err)
            derivativeNumber = 4 ! xy derivative
            CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
              & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo,1,0.0_CMISSRP,Err)
            CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
              & CMFE_FIELD_VALUES_SET_TYPE,1,derivativeNumber,GlobalNodeNo,2,0.0_CMISSRP,Err)
          END IF

          PRINT *, "Updated node number in rank"
          PRINT *, domain, GlobalNodeNo

        END IF ! Comp. node number == domain of the node

      END DO ! Global nodes

    END IF ! Lagrange Int Type

    IF(numberOfGlobalZElements/=0) THEN
      CALL HandleError("No manual update implemented for 3D!")
    END IF

    CALL cmfe_Field_ParameterSetUpdateStart(geometricField, &
      & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err) 
    ! Posts irecv and isend from and to adj domains for local ghosts
    ! Operations on internal dofs could happen in between
    CALL cmfe_Field_ParameterSetUpdateFinish(geometricField, &
      & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
    ! Waits for end of ghost communication

    WRITE(*,*) "Update completed."
    
  END IF ! Generated or manual mesh
  

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

  !Set the DOFs to be contiguous(?????) across components
  !FIELD_SEPARATED_COMPONENT_DOF_ORDER vs. FIELD_CONTIGUOUS_COMPONENT_DOF_ORDER
  !Note: for contiguous component DOF ordering all the components of the field variable 
  ! must have the same interpolation type.
  CALL cmfe_Field_DOFOrderTypeSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_SEPARATED_COMPONENT_DOF_ORDER,Err)
  CALL cmfe_Field_DOFOrderTypeSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,CMFE_FIELD_SEPARATED_COMPONENT_DOF_ORDER,Err)
  !CALL cmfe_Field_DOFOrderTypeSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_CONTIGUOUS_COMPONENT_DOF_ORDER,Err)
  !CALL cmfe_Field_DOFOrderTypeSet(DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,CMFE_FIELD_CONTIGUOUS_COMPONENT_DOF_ORDER,Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Initialise the field with an initial guess
  CALL cmfe_Field_ComponentValuesInitialise(DependentField, &
    & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,0.5_CMISSRP, Err)

  !-----------------------------------------------------------------------------------------------------------
  ! EQUATIONS
  !-----------------------------------------------------------------------------------------------------------

  !Create the equations set equations
  CALL cmfe_Equations_Initialise(Equations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  ! Affects backsubstitution and deludeln output!
  CALL cmfe_Equations_SparsityTypeSet(Equations,equationSparsity,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_NO_OUTPUT,Err)
! CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_TIMING_OUTPUT,Err)
  IF(ComputationalNodeNumber==1) &
    & CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_MATRIX_OUTPUT,Err)
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
! CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_NO_OUTPUT,Err)
! CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
! CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_TIMING_OUTPUT,Err)
! CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_SOLVER_OUTPUT,Err)
 CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_MATRIX_OUTPUT,Err)

!  CASE FULL
!  CALL cmfe_Solver_LinearTypeSet(Solver,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)
!  CALL cmfe_Solver_LinearIterativeAbsoluteToleranceSet(Solver,1.0E-12_CMISSRP,Err)
!  CALL cmfe_Solver_LinearIterativeRelativeToleranceSet(Solver,1.0E-12_CMISSRP,Err)
!  CASE SPARSE
  CALL cmfe_Solver_LinearTypeSet(Solver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)

! Do not select any for now
! CASE SPARSE  
! CALL cmfe_Solver_LibraryTypeSet(Solver,CMFE_SOLVER_MUMPS_LIBRARY,Err)

! CASE FULL 
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
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,solverSparsity,Err)
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
  !Set the first node to 0.0 and the last node to 1.0 (X-direction)
  FirstNodeNumber=1
  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Region_NodesGet(Region,Nodes,Err)
  CALL cmfe_Nodes_NumberOfNodesGet(Nodes,LastNodeNumber,Err)
  CALL cmfe_Decomposition_NodeDomainGet(Decomposition,FirstNodeNumber,1,FirstNodeDomain,Err)
  CALL cmfe_Decomposition_NodeDomainGet(Decomposition,LastNodeNumber,1,LastNodeDomain,Err)
  IF(FirstNodeDomain==ComputationalNodeNumber) THEN
    CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,FirstNodeNumber,1, &
      & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    !WRITE(*,*) "BC set 0"
    !CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,FirstNodeNumber,1, &
    !  & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
  END IF
  IF(LastNodeDomain==ComputationalNodeNumber) THEN
    CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,LastNodeNumber,1, &
      & CMFE_BOUNDARY_CONDITION_FIXED,1.0_CMISSRP,Err)
    !WRITE(*,*) "BC set 1"
    !CALL cmfe_BoundaryConditions_AddNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,LastNodeNumber,1, &
    !  & CMFE_BOUNDARY_CONDITION_FIXED,1.0_CMISSRP,Err)
  END IF
  !Finish the creation of the equations set boundary conditions
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !-----------------------------------------------------------------------------------------------------------
  !SOLVE
  !-----------------------------------------------------------------------------------------------------------

  !WRITE(*,*) "Stop program in laplace_equation.F90 before cmfe_Problem_Solve"
  !STOP
  !Solve the problem
  CALL cmfe_Problem_Solve(Problem,Err)
  !-----------------------------------------------------------------------------------------------------------
  !OUTPUT
  !-----------------------------------------------------------------------------------------------------------
  !Export results
  ! Set export file name
  WRITE(filename, "(A29,I1,A1,I1,A1,I1,A2,I1,A1,I1,A1,I1,A2,I1,A3,I1,A7,I1,A7,I1)") & !,A2,I1,A3,I1,A3,I1,A3,I1)") &
    & "results_current/current_run/l", &
    & INT(WIDTH),"x",INT(HEIGHT),"x",INT(LENGTH), &
    & "_n", &
    & numberOfGlobalXElements,"x",numberOfGlobalYElements,"x",numberOfGlobalZElements, &
    & "_i",lagrangeInterpolationType, "_np", NumberOfComputationalNodes, &
    & "_eqSpar", equationSparsity, &
    & "_soSpar", solverSparsity !,"_s",SolverIsDirect, &
    !& "_fd",JACOBIAN_FD,"_gm",useGeneratedMesh,"_bc",bcDirichlet
  ! make sure directories exist
  INQUIRE(file="./results_current/", exist=directory_exists)
  IF (.NOT.directory_exists) THEN
    CALL execute_command_line ("mkdir ./results_current/")
  END IF
  INQUIRE(file="./results_current/current_run/", exist=directory_exists)
  IF (.NOT.directory_exists) THEN
    CALL execute_command_line ("mkdir ./results_current/current_run/")
  END IF
  INQUIRE(file=trim(filename), exist=directory_exists)
  IF (.NOT.directory_exists) THEN
    CALL execute_command_line ("mkdir ./"//trim(filename))
  END IF

  ! Export solution
  CALL cmfe_Fields_Initialise(Fields,Err)
  CALL cmfe_Fields_Create(Region,Fields,Err)
  filename=trim(filename)//trim("/laplace_equation")
  CALL cmfe_Fields_NodesExport(Fields,filename,"FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(Fields,filename,"FORTRAN",Err)
  CALL cmfe_Fields_Finalise(Fields,Err)
  
  !Finialise CMISS
  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  ! CALL Get_MemoryConsumption()
  ! Output must be multiplied X number of cores!!!

  STOP
  
CONTAINS

  SUBROUTINE HandleError(errorString)
    CHARACTER(LEN=*), INTENT(IN) :: errorString
    WRITE(*,'(">>ERROR: ",A)') errorString(1:LEN_TRIM(errorString))
    STOP
  END SUBROUTINE HandleError

  SUBROUTINE Get_MemoryConsumption()
    CHARACTER(LEN=100) :: GetMemoryConsumption
    CHARACTER(LEN=100) ::  &
      & pid, comm, state, ppid, pgrp, session, tty_nr, &
      & tpgid, flags, minflt, cminflt, majflt, cmajflt, &
      & utime, stime, cutime, cstime, priority, nice, &
      & O, itrealvalue, starttime, Description, Limit
    CHARACTER(LEN=10000) :: Debug
    INTEGER(CMISSIntg) :: I,Stat
    INTEGER(CMISSLintg) :: MemoryConsumption, &
      & VmSize, VmRSS, Shared, Text, Lib, Data, Dt, RssLimBytes, RssAnon, Pagesize, VSizeBytes

    ! Critical Section
    DO I=0,NumberOfComputationalNodes
    !CALL MPI_BARRIER(MPI_COMM_WORLD, Err)
      IF (I == ComputationalNodeNumber) THEN

      ! read memory page size
        Stat = 0
        IF (I==0) CALL SYSTEM("getconf PAGESIZE > pagesize", Stat)
        IF (Stat == 0) THEN
          OPEN(UNIT=10, FILE="pagesize", ACTION="read", IOSTAT=Stat)
          IF (Stat == 0) THEN
            READ(10,*, IOSTAT=Stat) Pagesize
            CLOSE(UNIT=10)
          ELSE
            PRINT*, "Error opening pagesize."
          ENDIF
        ELSE
          PRINT*, "Error calling 'getconf PAGESIZE'"
          Pagesize = 4096
        ENDIF

      ! read from /proc/self/stat
      ! see http://man7.org/linux/man-pages/man5/proc.5.html for reference
        OPEN(UNIT=10, FILE="/proc/self/stat", ACTION="read", IOSTAT=stat)
        IF (STAT /= 0) THEN
          PRINT*, "Could not read memory consumption from /proc/self/stat."
        ELSE
        !READ(10,"(A)",IOSTAT=stat,advance='no') Debug
        !PRINT*, "proc/self/stat: "//TRIM(Debug)
          READ(10,*, IOSTAT=stat) pid, comm, state, ppid, pgrp, session, tty_nr, &
            & tpgid, flags, minflt, cminflt, majflt, cmajflt, &
            & utime, stime, cutime, cstime, priority, nice, &
            & O, itrealvalue, starttime, VSizeBytes, VmRSS, RssLimBytes
          CLOSE(UNIT=10)

          MemoryConsumption = VSizeBytes

          IF (ComputationalNodeNumber == 0) THEN
            WRITE(*, "(A,F7.3,A,I11,A,F7.3,A)") "     MemoryConsumption:", (MemoryConsumption/1e9), " GB (", &
              & MemoryConsumption, " B), Resident: ", (VmRSS*PageSize/1e9), " GB"
          ENDIF

          WRITE(GetMemoryConsumption, *) MemoryConsumption
        ENDIF

      ! read from /proc/self/limits
        IF (.TRUE.) THEN
          OPEN(UNIT=10, FILE="/proc/self/limits", ACTION="read", IOSTAT=stat)
          IF (STAT /= 0) THEN
            PRINT*, "Could not read limits from /proc/self/limits."
          ELSE
            DO
              READ(10, "(A26,A21)", IOSTAT=Stat) Description, Limit
              !PRINT*, "Description:["//TRIM(Description)//"], Limit:["//TRIM(Limit)//"]"
              IF (Stat /= 0) EXIT
              IF (INDEX(Description, "Max resident set") /= 0) THEN
                IF (TRIM(Limit) == "unlimited") THEN
                  IF (ComputationalNodeNumber == 0) PRINT*, "    (Resident has no soft limit)"
                ELSE
                  IF (ComputationalNodeNumber == 0) PRINT*, "    (Resident is limited to ", Limit,")"
                ENDIF
              ENDIF
            ENDDO
            CLOSE(UNIT=10)
          ENDIF
        ENDIF

      ! read from /proc/self/statm
        OPEN(UNIT=10, FILE="/proc/self/statm", ACTION="read", IOSTAT=stat)
        IF (STAT /= 0) THEN
          PRINT*, "Could not read memory consumption from /proc/self/statm."
        ELSE
          READ(10,*, IOSTAT=stat) VmSize, VmRSS, Shared, Text, Lib, Data, Dt
          CLOSE(UNIT=10)

        ! VmRSS = RssAnon + Shared (all as number of pages)
        ! RssAnon is the resident set size of anonymous memory (real memory in RAM, not laid out in files)
          RssAnon = VmRSS - Shared

          IF (ComputationalNodeNumber == 0) THEN
          !PRINT*, "VmSize: ", VmSize, ", VmRSS: ", VmRSS, ", Shared: ", Shared, ", Text: ", Text, ", Lib:", Lib, ", Data:", Data
          !PRINT*, "RssAnon: ", RSSAnon, ", RssLimBytes: ", RssLimBytes, ", Pagesize: ", Pagesize

          ! Output Percentage
          !WRITE(*, "(3(A,F7.3),A,F5.1,A)") "     VmSize:", (VmSize*Pagesize/1.e9), " GB, RssAnon:", (RssAnon*Pagesize/1.e9), &
          !  & " GB, RssLimBytes: ", (RssLimBytes/1.e9), " GB (", (REAL(RssAnon*Pagesize) / RssLimBytes * 100.0), "%)"
          ENDIF
        ENDIF

      !PRINT*, TRIM(cmfe_CustomProfilingGetInfo(Err))
      ENDIF
    ENDDO
  END SUBROUTINE Get_MemoryConsumption

    
END PROGRAM LAPLACE_EQUATION
