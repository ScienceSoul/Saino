//
//  FEMKernel.h
//  Saino
//
//  Created by Seddik hakime on 15/04/11.
//  Copyright 2011 ScienceSoul. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <OpenCL/OpenCL.h>

#import "FEMModel.h"
#import "FEMSolution.h"
#import "FEMMesh.h"

@interface FEMKernel : NSObject {
    
    int _normalTangentialNumberOfNodes;
    int _size1NtZeroingDone;
    int _size2NtZeroingDone;
    int _sizeBoundaryReorder;
    int _size1NtElement;
    int _size2NtElement;
    int _size1boundaryNormals;
    int _size2boundaryNormals;
    int _size1boundaryTangent1;
    int _size2boundaryTangent1;
    int _size1boundaryTangent2;
    int _size2boundaryTangent2;
    int _sizeIndexStore;
    int *_indexStore;
    bool **_ntZeroingDone;
    int *_boundaryReorder;
    int **_ntElement;
    double **_boundaryNormals;
    double **_boundaryTangent1;
    double **_boundaryTangent2;
    NSMutableString *_normalTangentialName;
    NSMutableArray *_outputLevelMask;
    BOOL _outputPrefix;
    BOOL _outputCaller;
    int _maxOutputLevel;
    int _minOutputLevel;
    int _outputPE;
}

@property(nonatomic, assign) int normalTangentialNumberOfNodes;
@property(nonatomic, assign) int size1NtZeroingDone;
@property(nonatomic, assign) int size2NtZeroingDone;
@property(nonatomic, assign) int sizeBoundaryReorder;
@property(nonatomic, assign) int size1boundaryNormals;
@property(nonatomic, assign) int size1NtElement;
@property(nonatomic, assign) int size2NtElement;
@property(nonatomic, assign) int size2boundaryNormals;
@property(nonatomic, assign) int size1boundaryTangent1;
@property(nonatomic, assign) int size2boundaryTangent1;
@property(nonatomic, assign) int size1boundaryTangent2;
@property(nonatomic, assign) int size2boundaryTangent2;
@property(nonatomic, assign) int sizeIndexStore;
@property(nonatomic, assign) int *indexStore;
@property(nonatomic, assign) bool **ntZeroingDone;
@property(nonatomic, assign) int *boundaryReorder;
@property(nonatomic, assign) int **ntElement;
@property(nonatomic, assign) double **boundaryNormals;
@property(nonatomic, assign) double **boundaryTangent1;
@property(nonatomic, assign) double **boundaryTangent2;
@property(nonatomic, strong) NSMutableString *normalTangentialName;
@property(nonatomic, strong) NSMutableArray *outputLevelMask;
@property(nonatomic, assign, getter = isOutputPrefix) BOOL outputPrefix;
@property(nonatomic, assign, getter = isOutputCaller) BOOL outputCaller;
@property(nonatomic, assign) int maxOutputLevel;
@property(nonatomic, assign) int minOutputLevel;
@property(nonatomic, assign) int outputPE;

// This class method retuns an instance (singleton) of FEMKernel
+(id)sharedKernel;

-(void)deallocation;

-(void)initializeTimeStepInSolution:(FEMSolution *)solution model:(FEMModel *)model;
-(void)initializeToZeroMatrix:(FEMMatrix *)matrix forceVector:(double *)forceVector sizeForceVector:(int)sizeForceVector model:(FEMModel *)model solution:(FEMSolution *)solution;
-(void)defaultInitializeSolution:(FEMSolution *)solution model:(FEMModel *)model;
-(BOOL)getReal:(FEMModel *)model forElement:(Element_t *)element inArray:(NSArray *)array variableName:(NSString *)name buffer:(listBuffer *)result;
-(int)isPElement:(Element_t *)element;
-(void)getNodes:(FEMSolution *)solution model:(FEMModel *)model inElement:(Element_t *)element resultNodes:(Nodes_t *)nodes numberOfNodes:(int *)nd;
-(int)getElementFamily:(Element_t *)element;
-(int)getElementDofsSolution:(FEMSolution *)uSolution model:(FEMModel *)model forElement:(Element_t *)element atIndexes:(int *)indexes;
-(int)sgetElementDofsSolution:(FEMSolution *)uSolution model:(FEMModel *)model forElement:(Element_t *)element atIndexes:(int *)indexes;
-(Element_t *)getActiveElement:(int) t solution:(FEMSolution *)solution model:(FEMModel *)model;
-(int)getNumberOfNodesForElement:(Element_t *)element;
-(int)getBoundaryConditionID:(FEMModel *)model forElement:(Element_t *)element;
-(NSArray *)getBoundaryCondition:(FEMModel *)model forElement:(Element_t *)element;
-(Element_t *)getBoundaryElement:(FEMSolution *)solution atIndex:(int)index;
-(void)getScalarLocalField:(double *)field sizeField:(int)sizeField name:(NSString *)name element:(Element_t *)element solution:(FEMSolution *)solution model:(FEMModel *)model timeStep:(int *)tStep;
-(void)getVectorLocalField:(double **)field size1Field:(int)size1Field size2Field:(int)size2Field name:(NSString *)name element:(Element_t *)element solution:(FEMSolution *)solution model:(FEMModel *)model timeStep:(int *)tStep;
-(int **)getEdgeMap:(int)elementFamily;
-(void)getBoundaryIndexes:(FEMMesh *)mesh forBoundaryElement:(Element_t *)element withParentElement:(Element_t *)parent resultVector:(int *)indexes sizeVector:(int)size indexSize:(int *)indexSize;
-(int)getBodyForceIDForElement:(Element_t *)element model:(FEMModel *)model;
-(int)getMaterialIDForElement:(Element_t *)element model:(FEMModel *)model;
-(int)getEquationIDForElement:(Element_t *)element model:(FEMModel *)model;
-(BOOL)getParentMaterialProperty:(NSString *)name forElement:(Element_t *)element parentElement:(Element_t *)parentElement model:(FEMModel *)model buffer:(listBuffer *)result;
-(BOOL)isActiveBoundaryElement:(Element_t *)element inSolution:(FEMSolution *)solution model:(FEMModel *)model;
-(void)checkNormalTangential:(FEMModel *)model inSolution:(FEMSolution *)solution forElementNumber:(int)elno numberofNodes:(int)n atIndexes:(int *)indexes atBoundary:(int)bc variableName:(NSMutableString *)name orderOfDofs:(int)dof activeCondition:(BOOL)conditional conditionName:(NSString *)condName permutationOffset:(int)offset;
-(void)zeroTheNumberOfRows:(int)n inSolutionMatrix:(FEMSolution *)solution;
-(void)setMatrixElementForSolution:(FEMSolution *)solution atIndex:(int)i andIndex:(int)j value:(double)value;
-(void)addToMatrixElementForSolution:(FEMSolution *)solution atIndex:(int)i andIndex:(int)j value:(double)value;
-(void)localBoundaryIntegral:(FEMModel *)model inSolution:(FEMSolution *)solution atBoundary:(NSArray *)bc forElement:(Element_t *)element withNumberOfNodes:(int)nd andParent:(Element_t *)parent withNumberOfNodes:(int)np boundaryName:(NSString *)name functionIntegral:(double *)integral;
-(void)localBoundaryBDOFs:(FEMModel *)model inSolution:(FEMSolution *)solution atBoundary:(NSArray *)bc forElement:(Element_t *)element withNumberOfNodes:(int)nd boundaryName:(NSMutableString *)name resultMatrix:(double **)stiff resultVector:(double *)force;
-(void)solveWithLapackMatrix:(double *)a andVector:(double *)x size:(int)n leadingDimension:(int)lda;
-(void)solveLinearSystemWithMatrix:(double **)a andVector:(double *)x size:(int)n leadingDimension:(int)lda;
-(void)setNodalLoads:(FEMModel *)model inSolution:(FEMSolution *)solution variableName:(NSString *)name orderOfDofs:(int)dof;
-(void)setDirichletBoundaries:(FEMModel *)model inSolution:(FEMSolution *)solution variableName:(NSMutableString *)name orderOfDofs:(int)dof permutationOffset:(int *)offset offDiaginalMatrix:(BOOL *)offDiaginalMatrix;
-(void)scaleLinearSystem:(FEMSolution *)solution matrix:(FEMMatrix *)matrix rhs:(double *)b result:(double *)x diagScaling:(double *)diagScaling applyScaling:(BOOL *)applyScaling rhsScaling:(BOOL *)rhsScaling;
-(void)backScaleLinearSystem:(FEMSolution *)solution matrix:(FEMMatrix *)matrix rhs:(double *)b result:(double *)x diagScaling:(double *)diagScaling sizeOFDiagScaling:(int *)sizeOfDiagScaling;
-(void)matrixVectorMultplyInSolution:(FEMSolution *)solution multiplyVector:(double *)u resultVector:(double *)v;
-(void)matrixVectorMultplyInMatrix:(FEMMatrix *)matrix multiplyVector:(double *)u resultVector:(double *)v;
-(void)invalidateVariableInTopMesh:(NSArray *)topMesh primaryMesh:(FEMMesh *)primaryMesh name:(NSString *)name model:(FEMModel *)model;
-(void)getPassiveBoundaryAtIndex:(int)bcID model:(FEMModel *)model mesh:(FEMMesh *)mesh solution:(FEMSolution *)solution;
-(void)computeNodalWeightsInSolution:(FEMSolution *)solution model:(FEMModel *)model weightBoundary:(BOOL)weightBoundary perm:(int *)perm sizePerm:(int *)sizePerm variableName:(NSString *)variableName;
-(void)condensateStiff:(double **)stiff force:(double *)force numberOfNodes:(int)n force1:(double *)force1;

// Manipulate matrix coefficients for time dependent simulations
-(void)addFirstOrderTimeModel:(FEMModel *)model solution:(FEMSolution *)solution element:(Element_t *)element massMatrix:(double **)massMatrix stiffMatrix:(double **)stiffMatrix force:(double *)force dt:(double)dt size:(int)n dofs:(int)dofs nodeIndexes:(int *)nodeIndexes rows:(int *)rows cols:(int *)cols;

-(void)updateMassMatrixModel:(FEMModel *)model inSolution:(FEMSolution *)solution localMassMatrix:(double **)localMassMatrix element:(Element_t *)element numberOfNodes:(int)n dofs:(int)dofs nodeIndexes:(int *)nodeIndexes;
-(void)defaultUpdateMass:(double **)mass element:(Element_t *)element solution:(FEMSolution *)solution model:(FEMModel *)model;
-(void)defaultUpdateComplexMass:(double complex **)cmass element:(Element_t *)element solution:(FEMSolution *)solution model:(FEMModel *)model;
-(void)defaultUpdateDamp:(double **)damp element:(Element_t *)element solution:(FEMSolution *)solution model:(FEMModel *)model;
-(void)defaultUpdateComplexDamp:(double complex **)cdamp element:(Element_t *)element solution:(FEMSolution *)solution model:(FEMModel *)model;
-(void)defaultFirstOrderTime:(FEMModel *)model inSolution:(FEMSolution *)solution forElement:(Element_t *)element realMass:(double **)mass realStiff:(double **)stiff realForce:(double *)force stiffRows:(int *)rows stiffCols:(int *)cols;
-(void)defaultFirstOrderTime:(FEMModel *)model inSolution:(FEMSolution *)solution forElement:(Element_t *)element complexMass:(double complex **)cmass complexStiff:(double complex **)cstiff complexForce:(double complex *)cforce stiffRows:(int *)rows stiffCols:(int *)cols;

-(void)defaultFirstOrderTimeGlobalModel:(FEMModel *)model inSolution:(FEMSolution *)solution;
-(void)defaultSecondOrderTimeGlobalModel:(FEMModel *)model inSolution:(FEMSolution *)solution;

// Update global equations
-(void)updateGlobalEquationsModel:(FEMModel *)model inSolution:(FEMSolution *)solution element:(Element_t *)element localStiffMatrix:(double **)localStiffMatrix forceVector:(double *)forceVector localForce:(double *)localForce size:(int)n dofs:(int)dofs nodeIndexes:(int *)nodeIndexes rows:(int *)rows cols:(int *)cols rotateNT:(BOOL *)rotateNT;
-(void)defaultUpdateEquations:(FEMModel *)model inSolution:(FEMSolution *)solution forElement:(Element_t *)element realStiff:(double **)stiff realForce:(double *)force stiffRows:(int *)rows stiffCols:(int *)cols requestBulkUpdate:(BOOL *)bulkUpdate;
-(void)defaultUpdateEquations:(FEMModel *)model inSolution:(FEMSolution *)solution forElement:(Element_t *)element complexStiff:(double **)cstiff complexForce:(double *)cforce stiffRows:(int *)rows stiffCols:(int *)cols requestBulkUpdate:(BOOL *)bulkUpdate;
-(void)defaultFinishAssemblySolution:(FEMSolution *)solution model:(FEMModel *)model;
-(void)dirichletBoundaryConditions:(FEMModel *)model inSolution:(FEMSolution *)solution usingOffset:(int *)offset offDiaginalMatrix:(BOOL *)offDiaginalMatrix;

-(void)computeChange:(FEMSolution *)solution model:(FEMModel *)aModel isSteadyState:(BOOL)steadyState nsize:(int*)nsize values:(double *)values values0:(double *)values0 sizeValues0:(int *)sizeValues0;
-(void)iterativeSolveMatrix:(FEMMatrix *)matrix result:(double *)x rhs:(double *)b dimensions:(int *)ndim solution:(FEMSolution *)solution;
-(void)solveSystemMatrix:(FEMMatrix *)matrix rhs:(double *)b result:(double *)x norm:(double *)norm dofs:(int)dofs solution:(FEMSolution *)solution  model:(FEMModel *)model;
-(double)findSolution:(FEMSolution *)solution model:(FEMModel *)model backRorateNT:(BOOL *)backRorateNT;

// Utility method in matrix-vector multiplication
-(double)stopc:(FEMMatrix *)matrix multiplyVector:(double *)x righHandSide:(double *)b ipar:(int *)ipar;

// Activating solution and starting the calculation
-(void)activateSolution:(FEMSolution *)solution model:(FEMModel *)model timeStep:(double)dt transientSimulation:(BOOL)transient;
-(void)solveEquationsModel:(FEMModel *)model timeStep:(double *)dt transientSimulation:(BOOL)transient coupledMinIteration:(int)coupledMinIter coupleMaxIteration:(int)coupleMaxIter steadyStateReached:(BOOL *)steadyStateReached realTimeStep:(int *)realTimeStep;

// OpenCL stuff
-(dispatch_queue_t)getDispatchQueueAndInfoForDeviceType:(NSString *)deviceType;

@end
