/*
 * 
 * This is an example program to solve Laplace's equation using OpenCMISS calls from C.
 * by Chris Bradley
 *
 */
#include <stdlib.h>
#include <stdio.h>

#include "opencmiss/iron.h"

#define STRING_SIZE 255

#define HEIGHT 1.0
#define WIDTH 2.0
#define LENGTH 3.0

#define COORDINATE_SYSTEM_USER_NUMBER 1
#define REGION_USER_NUMBER 2
#define BASIS_USER_NUMBER 3
#define GENERATED_MESH_USER_NUMBER 4
#define MESH_USER_NUMBER 5
#define DECOMPOSITION_USER_NUMBER 6
#define GEOMETRIC_FIELD_USER_NUMBER 7
#define EQUATIONS_SET_USER_NUMBER 8
#define EQUATIONS_SET_FIELD_USER_NUMBER 9
#define DEPENDENT_FIELD_USER_NUMBER 10
#define PROBLEM_USER_NUMBER 11

#define MAX_COORDINATES 3

#define CHECK_ERROR(S) \
  if(err != CMFE_NO_ERROR) { \
    if(err == CMFE_ERROR_CONVERTING_POINTER) { \
      fprintf(stderr,"Error: %s: Error converting pointer.\n",(S)); \
    } \
    else if(err == CMFE_POINTER_IS_NULL) { \
      fprintf(stderr,"Error: %s: Pointer is null.\n",(S)); \
    } \
    else if(err == CMFE_POINTER_NOT_NULL) { \
      fprintf(stderr,"Error: %s: Pointer is not null.\n",(S)); \
    } \
    else if(err == CMFE_COULD_NOT_ALLOCATE_POINTER) { \
      fprintf(stderr,"Error: %s: Could not allocate pointer.\n",(S)); \
    } \
    exit(err); \
  }

int main(int argc, char *argv[])
{
  cmfe_BasisType basis = (cmfe_BasisType)NULL;
  cmfe_BoundaryConditionsType boundaryConditions=(cmfe_BoundaryConditionsType)NULL;
  cmfe_ComputationEnvironmentType computationEnvironment=(cmfe_ComputationEnvironmentType)NULL;
  cmfe_CoordinateSystemType coordinateSystem=(cmfe_CoordinateSystemType)NULL,worldCoordinateSystem=(cmfe_CoordinateSystemType)NULL;
  cmfe_DecompositionType decomposition=(cmfe_DecompositionType)NULL;
  cmfe_EquationsType equations=(cmfe_EquationsType)NULL;
  cmfe_EquationsSetType equationsSet=(cmfe_EquationsSetType)NULL;
  cmfe_FieldsType fields=(cmfe_FieldsType)NULL;
  cmfe_FieldType geometricField=(cmfe_FieldType)NULL,dependentField=(cmfe_FieldType)NULL,equationsSetField=(cmfe_FieldType)NULL;
  cmfe_GeneratedMeshType generatedMesh=(cmfe_GeneratedMeshType)NULL;
  cmfe_MeshType mesh=(cmfe_MeshType)NULL;
  cmfe_ProblemType problem=(cmfe_ProblemType)NULL;
  cmfe_RegionType region=(cmfe_RegionType)NULL,worldRegion=(cmfe_RegionType)NULL;
  cmfe_SolverType solver=(cmfe_SolverType)NULL;
  cmfe_SolverEquationsType solverEquations=(cmfe_SolverEquationsType)NULL;

  int numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements,interpolationType;
  
  int numberOfComputationalNodes,computationalNodeNumber;
  int equationsSetIndex;
  int firstNodeNumber,lastNodeNumber;
  int firstNodeDomain,lastNodeDomain;

  int basisInterpolation[MAX_COORDINATES];
  int numberOfDimensions;
  int numberOfGauss[MAX_COORDINATES];
  int numberOfGaussXi;
  int numberXiElements[MAX_COORDINATES];
  int controlLoopIdentifier[1];
  double meshExtent[MAX_COORDINATES];

  int equationsSetSpecification[3];
  int problemSpecification[3];

  char filename[STRING_SIZE],format[STRING_SIZE],regionName[STRING_SIZE];

  int err;

  controlLoopIdentifier[0]=CMFE_CONTROL_LOOP_NODE;

  if(argc >= 4)
    {
      numberOfGlobalXElements = atoi(argv[1]);
      numberOfGlobalYElements = atoi(argv[2]);
      numberOfGlobalZElements = atoi(argv[3]);
      interpolationType = atoi(argv[4]);
    }
  else
    {
      numberOfGlobalXElements = 1;
      numberOfGlobalYElements = 3;
      numberOfGlobalZElements = 1;
      interpolationType = CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION;
   }

  err = cmfe_CoordinateSystem_Initialise(&worldCoordinateSystem);
  CHECK_ERROR("Initialising world coordinate system");
  err = cmfe_Region_Initialise(&worldRegion);
  CHECK_ERROR("Initialising world region");
  err = cmfe_Initialise(worldCoordinateSystem,worldRegion);
  CHECK_ERROR("Initialising OpenCMISS-Iron");
  err = cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR);

  sprintf(filename,"%s_%dx%dx%d_%d","Laplace",numberOfGlobalXElements,numberOfGlobalYElements,numberOfGlobalZElements,interpolationType);

  err = cmfe_OutputSetOn(STRING_SIZE,filename);

  /* Get the computational nodes information */
  err = cmfe_ComputationEnvironment_Initialise(&computationEnvironment);
  err = cmfe_ComputationEnvironment_NumberOfWorldNodesGet(computationEnvironment,&numberOfComputationalNodes);
  err = cmfe_ComputationEnvironment_WorldNodeNumberGet(computationEnvironment,&computationalNodeNumber);

  /* Start the creation of a new RC coordinate system */
  err = cmfe_CoordinateSystem_Initialise(&coordinateSystem);
  err = cmfe_CoordinateSystem_CreateStart(COORDINATE_SYSTEM_USER_NUMBER,coordinateSystem);
  if(numberOfGlobalZElements == 0)
    {
      /* Set the coordinate system to be 2D */
      numberOfDimensions = 2;
    }
  else
    {
      /* Set the coordinate system to be 3D */
      numberOfDimensions = 3;
    }
  err = cmfe_CoordinateSystem_DimensionSet(coordinateSystem,numberOfDimensions);
  /* Finish the creation of the coordinate system */
  err = cmfe_CoordinateSystem_CreateFinish(coordinateSystem);

  /* Start the creation of the region */
  sprintf(regionName,"%s","LaplaceEquation");
  err = cmfe_Region_Initialise(&region);
  err = cmfe_Region_CreateStart(REGION_USER_NUMBER,worldRegion,region);
  /* Set the regions coordinate system to the 2D RC coordinate system that we have created */
  err = cmfe_Region_CoordinateSystemSet(region,coordinateSystem);
  /* Set the regions name */
  err = cmfe_Region_LabelSet(region,STRING_SIZE,regionName);
  /* Finish the creation of the region */
  err = cmfe_Region_CreateFinish(region);

  /* Start the creation of a basis (default is trilinear lagrange) */
  err = cmfe_Basis_Initialise(&basis);
  err = cmfe_Basis_CreateStart(BASIS_USER_NUMBER,basis);
  switch(interpolationType)
    {
    case 1:
    case 2:
    case 3:
    case 4:
      err = cmfe_Basis_TypeSet(basis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE);
      break;
    case 7:
    case 8:
    case 9:
      err = cmfe_Basis_TypeSet(basis,CMFE_BASIS_SIMPLEX_TYPE);
      break;
    default:
      CHECK_ERROR("Invalid interpolation type");
    }
  switch(interpolationType)
    {
    case 1:
      numberOfGaussXi = 2;
      break;
    case 2:
      numberOfGaussXi = 3;
      break;
    case 3:
    case 4:
      numberOfGaussXi = 4;
      break;
    default:
      numberOfGaussXi = 0;
    }
  basisInterpolation[0] = interpolationType;
  basisInterpolation[1] = interpolationType;
  numberOfGauss[0] = numberOfGaussXi;
  numberOfGauss[1] = numberOfGaussXi;
  if(numberOfGlobalZElements != 0)
    {
      basisInterpolation[2] = interpolationType;
      numberOfGauss[2] = numberOfGaussXi;
    }
  err = cmfe_Basis_NumberOfXiSet(basis,numberOfDimensions);
  err = cmfe_Basis_InterpolationXiSet(basis,numberOfDimensions,basisInterpolation);
  err = cmfe_Basis_QuadratureNumberOfGaussXiSet(basis,numberOfDimensions,numberOfGauss);
  /* Finish the creation of the basis */
  err = cmfe_Basis_CreateFinish(basis);

  /* Start the creation of a generated mesh in the region */
  err = cmfe_GeneratedMesh_Initialise(&generatedMesh);
  err = cmfe_GeneratedMesh_CreateStart(GENERATED_MESH_USER_NUMBER,region,generatedMesh);
  /* Set up a regular x*y*z mesh */
  err = cmfe_GeneratedMesh_TypeSet(generatedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE);
  /* Set the default basis */
  err = cmfe_GeneratedMesh_BasisSet(generatedMesh,1,&basis);
  /* Define the mesh on the region */
  meshExtent[0] = WIDTH;
  meshExtent[1] = HEIGHT;
  numberXiElements[0] = numberOfGlobalXElements;
  numberXiElements[1] = numberOfGlobalYElements;
  if(numberOfGlobalZElements != 0)
    {
      meshExtent[2] = LENGTH;
      numberXiElements[2] = numberOfGlobalZElements;
    }
  err = cmfe_GeneratedMesh_ExtentSet(generatedMesh,MAX_COORDINATES,meshExtent);
  err = cmfe_GeneratedMesh_NumberOfElementsSet(generatedMesh,MAX_COORDINATES,numberXiElements);
  /* Finish the creation of a generated mesh in the region */
  err = cmfe_Mesh_Initialise(&mesh);
  /* Finish the creation of a generated mesh in the region */
  err = cmfe_GeneratedMesh_CreateFinish(generatedMesh,MESH_USER_NUMBER,mesh);

  /* Create a decomposition */
  err = cmfe_Decomposition_Initialise(&decomposition);
  err = cmfe_Decomposition_CreateStart(DECOMPOSITION_USER_NUMBER,mesh,decomposition);
  /* Set the decomposition to be a general decomposition with the specified number of domains */
  err = cmfe_Decomposition_TypeSet(decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE);
  err = cmfe_Decomposition_NumberOfDomainsSet(decomposition,numberOfComputationalNodes);
  /* Finish the decomposition */
  err = cmfe_Decomposition_CreateFinish(decomposition);

  /* Start to create a default (geometric) field on the region */
  err = cmfe_Field_Initialise(&geometricField);
  err = cmfe_Field_CreateStart(GEOMETRIC_FIELD_USER_NUMBER,region,geometricField);
  /* Set the decomposition to use */
  err = cmfe_Field_MeshDecompositionSet(geometricField,decomposition);
  /* Set the domain to be used by the field components. */
  err = cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,1);
  err = cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,1);
  if(numberOfGlobalZElements != 0)
    {
      err = cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,1);
    }
  /* Finish creating the field */
  err = cmfe_Field_CreateFinish(geometricField);

  /* Update the geometric field parameters */
  err = cmfe_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField);

  /* Create the equations_set */
  err = cmfe_EquationsSet_Initialise(&equationsSet);
  err = cmfe_Field_Initialise(&equationsSetField);
  equationsSetSpecification[0] = CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS;
  equationsSetSpecification[1] = CMFE_EQUATIONS_SET_LAPLACE_EQUATION_TYPE;
  equationsSetSpecification[2] = CMFE_EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE;
  err = cmfe_EquationsSet_CreateStart(EQUATIONS_SET_USER_NUMBER,region,geometricField, \
    3,equationsSetSpecification,EQUATIONS_SET_FIELD_USER_NUMBER, \
    equationsSetField,equationsSet);
  /* Finish creating the equations set */
  err = cmfe_EquationsSet_CreateFinish(equationsSet);

  /* Create the equations set dependent field variables */
  err = cmfe_Field_Initialise(&dependentField);
  err = cmfe_EquationsSet_DependentCreateStart(equationsSet,DEPENDENT_FIELD_USER_NUMBER,dependentField);
  /* Finish the equations set dependent field variables */
  err = cmfe_EquationsSet_DependentCreateFinish(equationsSet);

  /* Create the equations set equations */
  err = cmfe_Equations_Initialise(&equations);
  err = cmfe_EquationsSet_EquationsCreateStart(equationsSet,equations);
  /* Set the equations matrices sparsity type */
  err = cmfe_Equations_SparsityTypeSet(equations,CMFE_EQUATIONS_SPARSE_MATRICES);
  /* Set the equations set output */
  /* err = cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_NO_OUTPUT); */
  err = cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_TIMING_OUTPUT);
  /* err = cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_MATRIX_OUTPUT); */
  /* err = cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT); */
  /* Finish the equations set equations */
  err = cmfe_EquationsSet_EquationsCreateFinish(equationsSet);

  /* Start the creation of a problem, setting the problem to be a standard Laplace problem. */
  err = cmfe_Problem_Initialise(&problem);
  problemSpecification[0] = CMFE_PROBLEM_CLASSICAL_FIELD_CLASS;
  problemSpecification[1] = CMFE_PROBLEM_LAPLACE_EQUATION_TYPE;
  problemSpecification[2] = CMFE_PROBLEM_STANDARD_LAPLACE_SUBTYPE;
  err = cmfe_Problem_CreateStart(PROBLEM_USER_NUMBER,3,problemSpecification,problem);
  /* Finish the creation of a problem. */
  err = cmfe_Problem_CreateFinish(problem);

  /* Start the creation of the problem control loop */
  err = cmfe_Problem_ControlLoopCreateStart(problem);
  /* Finish creating the problem control loop */
  err = cmfe_Problem_ControlLoopCreateFinish(problem);

  /* Start the creation of the problem solvers */
  err = cmfe_Solver_Initialise(&solver);
  err = cmfe_Problem_SolversCreateStart(problem);
  err = cmfe_Problem_SolverGet(problem,1,controlLoopIdentifier,1,solver);
  /* err = cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_NO_OUTPUT); */
  /* err = cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_PROGRESS_OUTPUT); */
  /* err = cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_TIMING_OUTPUT); */
  /* err = cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_SOLVER_OUTPUT); */
  err = cmfe_Solver_OutputTypeSet(solver,CMFE_SOLVER_MATRIX_OUTPUT);
  err = cmfe_Solver_LinearTypeSet(solver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE);
  err = cmfe_Solver_LibraryTypeSet(solver,CMFE_SOLVER_MUMPS_LIBRARY);
  /* Finish the creation of the problem solver */
  err = cmfe_Problem_SolversCreateFinish(problem);

  /* Start the creation of the problem solver equations */
  solver=(cmfe_SolverType)NULL;
  err = cmfe_Solver_Initialise(&solver);
  err = cmfe_SolverEquations_Initialise(&solverEquations);
  err = cmfe_Problem_SolverEquationsCreateStart(problem);
  /* Get the solve equations */
  err = cmfe_Problem_SolverGet(problem,1,controlLoopIdentifier,1,solver);
  err = cmfe_Solver_SolverEquationsGet(solver,solverEquations);
  /* Set the solver equations sparsity */
  err = cmfe_SolverEquations_SparsityTypeSet(solverEquations,CMFE_SOLVER_SPARSE_MATRICES);
  /* err = cmfe_SolverEquations_SparsityTypeSet(solverEquations,CMFE_SOLVER_FULL_MATRICES);  */
  /* Add in the equations set */
  err = cmfe_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,&equationsSetIndex);
  /* Finish the creation of the problem solver equations */
  err = cmfe_Problem_SolverEquationsCreateFinish(problem);

  /* Start the creation of the equations set boundary conditions */
  err = cmfe_BoundaryConditions_Initialise(&boundaryConditions);
  err = cmfe_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions);
  /* Set the first node to 0.0 and the last node to 1.0 */
  firstNodeNumber = 1;
  if(numberOfGlobalZElements == 0)
    {
      lastNodeNumber = (numberOfGlobalXElements+1)*(numberOfGlobalYElements+1);
    }
  else
    {
      lastNodeNumber = (numberOfGlobalXElements+1)*(numberOfGlobalYElements+1)*(numberOfGlobalZElements+1);
    }
  err = cmfe_Decomposition_NodeDomainGet(decomposition,firstNodeNumber,1,&firstNodeDomain);
  err = cmfe_Decomposition_NodeDomainGet(decomposition,lastNodeNumber,1,&lastNodeDomain);
  if(firstNodeDomain == computationalNodeNumber)
    {
      err = cmfe_BoundaryConditions_SetNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,firstNodeNumber,1, \
        CMFE_BOUNDARY_CONDITION_FIXED,0.0);
    }
  if(lastNodeDomain == computationalNodeNumber)
    {
      err = cmfe_BoundaryConditions_SetNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,lastNodeNumber,1, \
        CMFE_BOUNDARY_CONDITION_FIXED,1.0);
    }
  /* Finish the creation of the equations set boundary conditions */
  err = cmfe_SolverEquations_BoundaryConditionsCreateFinish(solverEquations);

  /* Solve the problem */
  err = cmfe_Problem_Solve(problem);

  /* Output results */
  sprintf(filename,"%s","LaplaceEquation");
  sprintf(format,"%s","FORTRAN");
  err = cmfe_Fields_Initialise(&fields);
  err = cmfe_Fields_CreateRegion(region,fields);
  err = cmfe_Fields_NodesExport(fields,STRING_SIZE,filename,STRING_SIZE,format);
  err = cmfe_Fields_ElementsExport(fields,STRING_SIZE,filename,STRING_SIZE,format);
  err = cmfe_Fields_Finalise(&fields);

  /* Finalise OpenCMISS */
  err = cmfe_Finalise();

  return err;
}
