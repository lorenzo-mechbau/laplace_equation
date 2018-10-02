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
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=13
 
  !Program types
  
  !Program variables
  INTEGER(CMISSIntg) :: RowNo, ColumnNo, GlobalElementNo, GlobalNodeNo, I
  INTEGER(CMISSIntg) :: NumberOfElements, NumberOfNodes

  INTEGER(CMISSIntg) :: numberOfArguments,argumentLength,status
  INTEGER(CMISSIntg) :: numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements, &
    & interpolationType,numberOfGaussXi
  CHARACTER(LEN=255) :: commandArgument,filename

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
  INTEGER(CMISSIntg) :: count_components = 1  
  CHARACTER(LEN=60) :: diagFilename
  INTEGER(CMISSIntg), ALLOCATABLE ::elementUserNodes(:)
  LOGICAL, PARAMETER             :: useGeneratedMesh = .FALSE.
  LOGICAL, PARAMETER             :: allLinear = .TRUE. ! For MANUAL mesh, decide if all linear 
                                                       ! or Hermite/linear (cf. Benjamin's example)  
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
    numberOfGlobalXElements=4
    numberOfGlobalYElements=4
    numberOfGlobalZElements=0
    interpolationType=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
  ENDIF


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

  ! Print out of nodes/els mapping (old implementation) + dofs (new implementation)
  !CALL cmfe_DiagnosticsSetOn(CMFE_IN_DIAG_TYPE,[2],"Diagnostics",&
  !  & ["FIELD_MAPPINGS_CALCULATE"],Err)

  ! Uncomment this if on one of the merge branches! 
  ! Determine name of diagnostics file according to branch
  !SELECT CASE(CMFE_BENJAMIN)
  !CASE (CMFE_B_ORIGINAL)
  !  diagFilename = "Diag_original"
  !CASE (CMFE_B_MERGE)
  !  diagFilename = "Diag_merge"
  !CASE (CMFE_B_FACES)
  !  diagFilename = "Diag_faces"
  !CASE DEFAULT
  !  diagFilename = "Diag_other"
  !END SELECT
  diagFilename = "Diagnostics"

  CALL cmfe_DiagnosticsSetOn(CMFE_IN_DIAG_TYPE,[1,2,3],TRIM(diagFileName),&
    & ["DOMAIN_MAPPINGS_NODES_CALCULATE", &
    &  "DOMAIN_MAPPINGS_INITIALISE     "], Err)
!    &  "FIELD_MAPPINGS_CALCULATE       "],Err)
  
  
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
  IF(numberOfGlobalZElements==0) THEN
    !Set the coordinate system to be 2D
    CALL cmfe_CoordinateSystem_DimensionSet(coordinateSystem,2,err)
  ELSE
    !Set the coordinate system to be 3D
    CALL cmfe_CoordinateSystem_DimensionSet(coordinateSystem,3,err)
  ENDIF

  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !-----------------------------------------------------------------------------------------------------------
  !REGION
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of the region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  CALL cmfe_Region_LabelSet(Region,"LaplaceEquation",Err)
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

  SELECT CASE(interpolationType)
  CASE(1,2,3,4)
    CALL cmfe_Basis_TypeSet(LinearBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,err)
  CASE(7,8,9)
    CALL cmfe_Basis_TypeSet(LinearBasis,CMFE_BASIS_SIMPLEX_TYPE,err)
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
    CALL cmfe_Basis_NumberOfXiSet(LinearBasis,2,err)
    CALL cmfe_Basis_InterpolationXiSet(LinearBasis,[interpolationType,interpolationType],err)
    IF(numberOfGaussXi>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(LinearBasis,[numberOfGaussXi,numberOfGaussXi],err)
    ENDIF
  ELSE
    !Set the basis to be a tri-interpolation basis
    CALL cmfe_Basis_NumberOfXiSet(LinearBasis,3,err)
    CALL cmfe_Basis_InterpolationXiSet(LinearBasis,[interpolationType,interpolationType,interpolationType],err)
    IF(numberOfGaussXi>0) THEN
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(LinearBasis,[numberOfGaussXi,numberOfGaussXi,numberOfGaussXi],err)
    ENDIF
  ENDIF

  !CALL cmfe_Basis_NumberOfXiSet(LinearBasis,2,Err)
  !CALL cmfe_Basis_InterpolationXiSet(LinearBasis,&
  !  & [CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION,CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err) !i.e. 1
  !CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(LinearBasis,[2,2],Err)
  
  CALL cmfe_Basis_CreateFinish(LinearBasis,Err)
  
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
      CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(HermiteBasis,[3,3],err)
    ENDIF
  ELSE
    CALL HandleError("No Hermite basis defined for 3D case.")
    !Set the basis to be a tri-interpolation basis
    !CALL cmfe_Basis_NumberOfXiSet(HermiteBasis,3,err)
    !CALL cmfe_Basis_InterpolationXiSet(HermiteBasis,[interpolationType,interpolationType,interpolationType],err)
    !IF(numberOfGaussXi>0) THEN
    !  CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(HermiteBasis,[numberOfGaussXi,numberOfGaussXi,numberOfGaussXi],err)
    !ENDIF
  ENDIF


  !CALL cmfe_Basis_NumberOfXiSet(HermiteBasis,2,Err)
  !CALL cmfe_Basis_InterpolationXiSet(HermiteBasis,&
  !  & [CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION,CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION],Err)
  !CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(HermiteBasis,[3,3],Err)
  
  CALL cmfe_Basis_CreateFinish(HermiteBasis,Err)
   
  !-----------------------------------------------------------------------------------------------------------
  !MESH
  !-----------------------------------------------------------------------------------------------------------

  !Start the creation of a generated mesh in the region
  IF (useGeneratedMesh) THEN

    CALL cmfe_GeneratedMesh_Initialise(generatedMesh,err)
    CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,generatedMesh,err)
    !Set up a regular x*y*z mesh
    CALL cmfe_GeneratedMesh_TypeSet(generatedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,err)
    !Set the default basis
    !Using different basis types is not supported for generated meshes! ==> ALL LINEAR!!!
    CALL cmfe_GeneratedMesh_BasisSet(generatedMesh,LinearBasis,err) 
    !Define the mesh on the region
    IF(numberOfGlobalZElements==0) THEN
      CALL cmfe_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,HEIGHT],err)
      CALL cmfe_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements,numberOfGlobalYElements],err)
    ELSE
      CALL cmfe_GeneratedMesh_ExtentSet(generatedMesh,[WIDTH,HEIGHT,LENGTH],err)
      CALL cmfe_GeneratedMesh_NumberOfElementsSet(generatedMesh,[numberOfGlobalXElements,numberOfGlobalYElements, &
        & numberOfGlobalZElements],err)
    ENDIF
    !Finish the creation of a generated mesh in the region
    CALL cmfe_Mesh_Initialise(Mesh,err)
    CALL cmfe_GeneratedMesh_CreateFinish(generatedMesh,MeshUserNumber,Mesh,err)

  ELSE
  ! Create MANUAL mesh as specified in "PartitionedMesh_4x4.pdf"
  
    CALL cmfe_Mesh_Initialise(Mesh,Err)
    CALL cmfe_Mesh_CreateStart(MeshUserNumber,Region,2,Mesh,Err) ! 2 is number of dim
    CALL cmfe_Mesh_NumberOfComponentsSet(Mesh,1,Err)

    NumberOfElements = 16
    ! Special case 3 nodes
    IF (NumberOfComputationalNodes == 3) NumberOfElements = 3

    CALL cmfe_Mesh_NumberOfElementsSet(Mesh,NumberOfElements,Err)
  
    ! Create elements 
    CALL cmfe_MeshElements_Initialise(MeshElements,Err)
    CALL cmfe_MeshElements_CreateStart(Mesh,1,LinearBasis,MeshElements,Err)
  
    ! Define nodes for the mesh
    CALL cmfe_Nodes_Initialise(Nodes,Err)

    NumberOfNodes = 25
    ! Special case 3 nodes
    IF (NumberOfComputationalNodes == 3) NumberOfNodes = 8

    CALL cmfe_Nodes_CreateStart(Region,NumberOfNodes,Nodes,Err)
    CALL cmfe_Nodes_CreateFinish(Nodes,Err)

  ! Special case 3 nodes
  ! Set adjacent nodes for each element
    IF (NumberOfComputationalNodes == 3) THEN
      CALL cmfe_MeshElements_NodesSet(MeshElements, 1, [1,2,5,6], Err)
      CALL cmfe_MeshElements_NodesSet(MeshElements, 2, [2,3,6,7], Err)
      CALL cmfe_MeshElements_NodesSet(MeshElements, 3, [3,4,7,8], Err)
      CALL cmfe_MeshElements_BasisSet(MeshElements,1,LinearBasis,Err)
      CALL cmfe_MeshElements_BasisSet(MeshElements,2,HermiteBasis,Err)
      CALL cmfe_MeshElements_BasisSet(MeshElements,3,LinearBasis,Err)
    ELSE

    ! Distribute nodes among elements (cf. picture)
    DO RowNo = 1,4
      DO ColumnNo = 1,4
        GlobalElementNo = (RowNo-1)*4 + ColumnNo
        GlobalNodeNo = (RowNo-1)*5 + ColumnNo
        ! OK
        CALL cmfe_MeshElements_NodesSet(MeshElements, GlobalElementNo, &
          & [GlobalNodeNo,GlobalNodeNo+1,GlobalNodeNo+5,GlobalNodeNo+6], Err)

        PRINT *, "Element"
        PRINT *, GlobalElementNo
        PRINT *, "Nodes"
        PRINT *, [GlobalNodeNo,GlobalNodeNo+1,GlobalNodeNo+5,GlobalNodeNo+6]          
      END DO
    END DO

      ! Set the basis
      !IF (.FALSE.) THEN
    IF (.NOT.allLinear) THEN
      DO GlobalElementNo=1,4
        CALL cmfe_MeshElements_BasisSet(MeshElements,GlobalElementNo,HermiteBasis,Err)
      ENDDO
      DO GlobalElementNo=5,12
        CALL cmfe_MeshElements_BasisSet(MeshElements,GlobalElementNo,LinearBasis,Err)
      ENDDO
      DO GlobalElementNo=13,16
        CALL cmfe_MeshElements_BasisSet(MeshElements,GlobalElementNo,HermiteBasis,Err)
      ENDDO
    ELSE    ! all linear
      DO GlobalElementNo=1,16
        CALL cmfe_MeshElements_BasisSet(MeshElements,GlobalElementNo,LinearBasis,Err)
      ENDDO
    ENDIF

  END IF ! number of comp. nodes (3 vs. other)
  
  CALL cmfe_MeshElements_CreateFinish(MeshElements,Err)
  
  CALL cmfe_Mesh_CreateFinish(Mesh,Err)

  END IF

  ! Check the nodes distribution
  IF (.NOT. useGeneratedMesh) THEN
    ALLOCATE(elementUserNodes(4),STAT=Err)
    DO RowNo = 1,4
      DO ColumnNo = 1,4
        GlobalElementNo = (RowNo-1)*4 + ColumnNo
        ! Check with nodes get
        CALL cmfe_MeshElements_NodesGet(MeshElements,GlobalElementNo,elementUserNodes,Err)
        PRINT *, "Element"
        PRINT *, GlobalElementNo
        PRINT *, "Nodes"
        PRINT *, elementUserNodes          
      ENDDO
    ENDDO
  END IF 
  
  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  ! Set the decomposition to be a general decomposition with the specified number of domains

  ! Standard (old) decomposition
  ! No need
  !CALL cmfe_Decomposition_TypeSet(decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,err)
  !CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)

  ! User-defined decomposition
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_USER_DEFINED_TYPE,Err)
  
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  
  ! global element numbers
  !13 14 15 16
  ! 9 10 11 12
  ! 5  6  7  8
  ! 1  2  3  4
  
  PRINT *, "Computational node number:"
  PRINT *, ComputationalNodeNumber
  
  SetDecompositionDistributed = .FALSE. ! True only works with new implementation, false works with both.
                                        ! Affects decomposition for 4 nodes (below).
                                        ! NOT SURE WHAT IT MEANS THOUGH!!!!!!!!!!
  SELECT CASE (NumberOfComputationalNodes)
  CASE (1)
    DO I=1,16
      CALL cmfe_Decomposition_ElementDomainSet(Decomposition,I,0,Err)
    ENDDO
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
  CASE DEFAULT
    CALL HandleError("Number of nodes NOT supported (max 4)!!!")
  END SELECT
  
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)
 
  !Destroy the mesh now that we have decomposed it
  !CALL cmfe_Mesh_Destroy(Mesh,Err)

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

  !1 vbl, 1 component: 1x Mesh (node-dofs)
  ! CALL cmfe_Field_NumberOfComponentsSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,Err)
 
  !1 vbl, 2 components: 2x Mesh (node-dofs)
  CALL cmfe_Field_NumberOfComponentsSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,Err)
  IF(numberOfGlobalZElements/=0) THEN
    CALL cmfe_Field_NumberOfComponentsSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,Err)
  ENDIF
 
  !1 vbl, 3 components: 2x Mesh (node-dofs), 1x Interpolation (element-dofs)
  !CALL cmfe_Field_NumberOfComponentsSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,Err)

  !1 vbl, 4 components: 2x Mesh (node-dofs), 2x Interpolation (1 element-dofs, 1 constant = 1 dof)
  !CALL cmfe_Field_NumberOfComponentsSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,4,Err)

  !SUBROUTINE FIELD_COMPONENT_MESH_COMPONENT_SET(FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,MESH_COMPONENT_NUMBER,ERR,ERROR,*)

  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,count_components,1,Err)
  count_components = count_components +1
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,count_components,1,Err)
  IF(numberOfGlobalZElements/=0) THEN
    count_components = count_components +1
    CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,count_components,1,err)
  ENDIF
!  count_components = count_components +1
!  CALL cmfe_Field_ComponentInterpolationSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,count_components, &
!  & CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
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
  
  !PRINT *, "Abort program in laplace_equation.f90:360"
  !STOP
 
  
  ! Update the geometric field parameters
  IF (useGeneratedMesh) THEN
    CALL cmfe_GeneratedMesh_GeometricParametersCalculate(generatedMesh,GeometricField,Err)
  ELSE
    ! Update the geometric field parameters manually if mesh is NOT a generated mesh
    CALL cmfe_Field_ParameterSetUpdateStart(geometricField, &
      & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
    ! node 1
    !SUBROUTINE FIELD_PARAMETER_SET_UPDATE_NODE_INTG(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,VERSION_NUMBER,DERIVATIVE_NUMBER, &
    !  & USER_NODE_NUMBER,COMPONENT_NUMBER,VALUE,ERR,ERROR,*)
    DO I=1,NumberOfNodes
      SELECT CASE (I)
      CASE (1,6,11,16,21)
        CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_FIELD_VALUES_SET_TYPE,1,1,I,1,0.0_CMISSRP,Err)
      CASE (2,7,12,17,22)
        CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_FIELD_VALUES_SET_TYPE,1,1,I,1,WIDTH/4.0_CMISSRP,Err)
      CASE (3,8,13,18,23)
        CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_FIELD_VALUES_SET_TYPE,1,1,I,1,WIDTH/2.0_CMISSRP,Err)
      CASE (4,9,14,19,24)
        CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_FIELD_VALUES_SET_TYPE,1,1,I,1,3.0_CMISSRP*WIDTH/4.0_CMISSRP,Err)
      CASE (5,10,15,20,25)
        CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_FIELD_VALUES_SET_TYPE,1,1,I,1,WIDTH,Err)
      CASE DEFAULT
        CALL HandleError("Node out of reach!")
      END SELECT

      SELECT CASE (I)
      CASE (1,2,3,4,5)
        CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_FIELD_VALUES_SET_TYPE,1,1,I,2,0.0_CMISSRP,Err)
      CASE (6,7,8,9,10)
        CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_FIELD_VALUES_SET_TYPE,1,1,I,2,WIDTH/4.0_CMISSRP,Err)
      CASE (11,12,13,14,15)
        CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_FIELD_VALUES_SET_TYPE,1,1,I,2,WIDTH/2.0_CMISSRP,Err)      
      CASE (16,17,18,19,20)
        CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_FIELD_VALUES_SET_TYPE,1,1,I,2,3.0_CMISSRP*WIDTH/4.0_CMISSRP,Err)
      CASE (21,22,23,24,25)
        CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_FIELD_VALUES_SET_TYPE,1,1,I,2,WIDTH,Err)
      CASE DEFAULT
        CALL HandleError("Node out of reach!")
      END SELECT

    END DO

    IF(numberOfGlobalZElements/=0) THEN
      CALL HandleError("Node out of reach!")
    END IF
 
    CALL cmfe_Field_ParameterSetUpdateFinish(geometricField, &
      & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
    
    !CALL HandleError("No geometric field parameters update for manual mesh! Cannot solve!")

  END IF
  

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
  CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_MATRIX_OUTPUT,Err)
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
  CALL cmfe_Fields_NodesExport(Fields,"LaplaceEquation","FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(Fields,"LaplaceEquation","FORTRAN",Err)
  CALL cmfe_Fields_Finalise(Fields,Err)
  
  !Finialise CMISS
  CALL cmfe_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
CONTAINS

  SUBROUTINE HandleError(errorString)
    CHARACTER(LEN=*), INTENT(IN) :: errorString
    WRITE(*,'(">>ERROR: ",A)') errorString(1:LEN_TRIM(errorString))
    STOP
  END SUBROUTINE HandleError
    
END PROGRAM LAPLACE_EQUATION
