//
//  FEMCore.h
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
#import "FEMMatrixCRS.h"
#import "FEMMatrixBand.h"
#import "FEMListUtilities.h"
#import "FEMUtilities.h"
#import "FEMTimeIntegration.h"

@interface FEMCore : NSObject {
    
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
@property(nonatomic, assign, nonnull) int *indexStore;
@property(nonatomic, assign) bool * __nullable * __nullable ntZeroingDone;
@property(nonatomic, assign, nullable) int *boundaryReorder;
@property(nonatomic, assign) int * __nullable * __nullable ntElement;
@property(nonatomic, assign) double * __nullable * __nullable boundaryNormals;
@property(nonatomic, assign) double * __nullable * __nullable boundaryTangent1;
@property(nonatomic, assign) double * __nullable * __nullable boundaryTangent2;
@property(nonatomic, strong, nonnull) NSMutableString *normalTangentialName;
@property(nonatomic, strong, nonnull) NSMutableArray <NSNumber *> *outputLevelMask;
@property(nonatomic, assign, getter = isOutputPrefix) BOOL outputPrefix;
@property(nonatomic, assign, getter = isOutputCaller) BOOL outputCaller;
@property(nonatomic, assign) int maxOutputLevel;
@property(nonatomic, assign) int minOutputLevel;
@property(nonatomic, assign) int outputPE;

// This class method retuns an instance (singleton) of FEMCore
+(id __nonnull)sharedCore;
+(void)selfDestruct;

-(void)deallocation;

-(void)initializeTimeStepInSolution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model;
-(void)initializeToZeroMatrix:(FEMMatrix * __nonnull)matrix forceVector:(double * __nonnull)forceVector sizeForceVector:(int)sizeForceVector model:(FEMModel * __nonnull)model solution:(FEMSolution * __nonnull)solution ;
-(void)defaultInitializeSolution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model;
-(BOOL)getReal:(FEMModel * __nonnull)model forElement:(Element_t * __nullable)element inArray:(NSArray * __nonnull)array variableName:(NSString * __nonnull)name buffer:(listBuffer * __nonnull)result listUtilities:(FEMListUtilities * __nonnull)listUtil;
-(int)isPElement:(Element_t * __nonnull)element;
-(void)getNodes:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model inElement:(Element_t * __nonnull)element resultNodes:(Nodes_t * __nonnull)nodes numberOfNodes:(int * __nullable)nd mesh:(FEMMesh * __nullable)mesh;
-(BOOL)isFluxElement:(Element_t * __nonnull)element mesh:(FEMMesh * __nonnull)mesh ;
-(int)getElementFamily:(Element_t * __nonnull)element;
-(int)getElementDofsSolution:(FEMSolution * __nullable)uSolution model:(FEMModel * __nonnull)model forElement:(Element_t * __nonnull)element atIndexes:(int * __nonnull)indexes disableDiscontinuousGalerkin:(BOOL * __nullable)disableDiscontinuousGalerkin;
-(int)sgetElementDofsSolution:(FEMSolution * __nullable)uSolution model:(FEMModel * __nonnull)model forElement:(Element_t * __nonnull)element atIndexes:(int * __nonnull)indexes;
-(int)getNumberOfBubbleDofsElement:(Element_t * __nonnull)element solution:(FEMSolution * __nonnull)solution;
-(Element_t * __nonnull)getActiveElement:(int)t solution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model;
-(int)getNumberOfNodesForElement:(Element_t * __nonnull)element;
-(int)getBoundaryConditionID:(FEMModel * __nonnull)model forElement:(Element_t * __nonnull)element;
-(NSArray * __nullable)getBoundaryCondition:(FEMModel * __nonnull)model forElement:(Element_t * __nonnull)element ;
-(Element_t * __nonnull)getBoundaryElement:(FEMSolution * __nonnull)solution atIndex:(int)index ;
-(void)getScalarLocalField:(double * __nonnull)field sizeField:(int)sizeField name:(NSString * __nullable)name element:(Element_t * __nonnull)element solution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model timeStep:(int * __nullable)tStep;
-(void)getVectorLocalField:(double * __nonnull * __nonnull)field size1Field:(int)size1Field size2Field:(int)size2Field name:(NSString * __nullable)name element:(Element_t * __nonnull)element solution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model timeStep:(int * __nullable)tStep;
-(int * __nonnull * __nonnull)getEdgeMap:(int)elementFamily mapSize:(int * __nullable)mapSize;
-(void)getBoundaryIndexes:(FEMMesh * __nonnull)mesh forBoundaryElement:(Element_t * __nonnull)element withParentElement:(Element_t * __nonnull)parent resultVector:(int * __nonnull)indexes sizeVector:(int)size indexSize:(int * __nonnull)indexSize;
-(int)getBodyForceIDForElement:(Element_t * __nonnull)element model:(FEMModel * __nonnull)model;
-(int)getMaterialIDForElement:(Element_t * __nonnull)element model:(FEMModel * __nonnull)model;
-(int)getEquationIDForElement:(Element_t * __nonnull)element model:(FEMModel * __nonnull)model;
-(BOOL)getParentMaterialProperty:(NSString * __nonnull)name forElement:(Element_t * __nonnull)element parentElement:(Element_t * __nullable)parentElement model:(FEMModel * __nonnull)model listUtilities:(FEMListUtilities * __nonnull)listUtilities buffer:(listBuffer * __nonnull)result;
-(BOOL)isActiveBoundaryElement:(Element_t * __nonnull)element inSolution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model;
-(void)checkNormalTangential:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution forElementNumber:(int)elno numberofNodes:(int)n atIndexes:(int * __nonnull)indexes atBoundary:(int)bc variableName:(NSMutableString * __nonnull)name orderOfDofs:(int)dof activeCondition:(BOOL)conditional conditionName:(NSString * __nonnull)condName permutationOffset:(int)offset;
-(void)zeroTheNumberOfRows:(int)n inSolutionMatrix:(FEMSolution * __nonnull)solution;
-(void)setMatrixElementForSolution:(FEMSolution * __nonnull)solution atIndex:(int)i andIndex:(int)j value:(double)value;
-(void)addToMatrixElementForSolution:(FEMSolution * __nonnull)solution atIndex:(int)i andIndex:(int)j value:(double)value;
-(void)localBoundaryIntegral:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution atBoundary:(NSArray * __nonnull)bc forElement:(Element_t * __nonnull)element withNumberOfNodes:(int)nd andParent:(Element_t * __nonnull)parent withNumberOfNodes:(int)np boundaryName:(NSString * __nonnull)name functionIntegral:(double * __nonnull)integral;
-(void)localBoundaryBDOFs:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution atBoundary:(NSArray * __nonnull)bc forElement:(Element_t * __nonnull)element withNumberOfNodes:(int)nd boundaryName:(NSMutableString * __nonnull)name resultMatrix:(double * __nonnull * __nonnull)stiff resultVector:(double * __nonnull)force;
-(void)solveWithLapackMatrix:(double * __nonnull)a andVector:(double * __nonnull)x size:(int)n leadingDimension:(int)lda;
-(void)solveLinearSystemWithMatrix:(double * __nonnull * __nonnull)a andVector:(double * __nonnull)x size:(int)n leadingDimension:(int)lda;
-(void)setDirichletBoundaries:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution variableName:(NSMutableString * __nonnull)name orderOfDofs:(int)dof permutationOffset:(int * __nullable)offset offDiaginalMatrix:(BOOL * __nullable)offDiaginalMatrix listUtilites:(FEMListUtilities * __nonnull)listUtilites crsMatrix:(FEMMatrixCRS * __nonnull)crsMatrix bandMatrix:(FEMMatrixBand * __nonnull)bandMatrix;
-(void)scaleLinearSystem:(FEMSolution * __nonnull)solution matrix:(FEMMatrix * __nonnull)matrix rhs:(double * __nullable)b result:(double * __nullable)x diagScaling:(double * __nullable)diagScaling applyScaling:(BOOL * __nullable)applyScaling rhsScaling:(BOOL * __nullable)rhsScaling;
-(void)backScaleLinearSystem:(FEMSolution * __nonnull)solution matrix:(FEMMatrix * __nonnull)matrix rhs:(double * __nullable)b result:(double * __nullable)x diagScaling:(double * __nullable)diagScaling sizeOFDiagScaling:(int * __nullable)sizeOfDiagScaling;
-(void)matrixVectorMultplyInSolution:(FEMSolution * __nonnull)solution multiplyVector:(double * __nonnull)u resultVector:(double * __nonnull)v;
-(void)matrixVectorMultplyInMatrix:(FEMMatrix * __nonnull)matrix multiplyVector:(double * __nonnull)u resultVector:(double * __nonnull)v;
-(void)invalidateVariableInTopMesh:(NSArray * __nonnull)topMesh primaryMesh:(FEMMesh * __nonnull)primaryMesh name:(NSString * __nonnull)name model:(FEMModel * __nonnull)model;
-(void)getPassiveBoundaryAtIndex:(int)bcID model:(FEMModel * __nonnull)model mesh:(FEMMesh * __nonnull)mesh solution:(FEMSolution * __nonnull)solution;
-(void)computeNodalWeightsInSolution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model weightBoundary:(BOOL)weightBoundary perm:(int * __nullable)perm sizePerm:(int * __nullable)sizePerm variableName:(NSString * __nullable)variableName;
-(void)nsCondensateStiff:(double * __nonnull * __nonnull)stiff force:(double * __nonnull)force numberOfNodes:(int)n numberOfBubbles:(int)nb dimension:(int)dim force1:(double * __nonnull)force1;
-(void)condensateStiff:(double * __nonnull * __nonnull)stiff force:(double * __nonnull)force numberOfNodes:(int)n force1:(double * __nullable)force1;

// Manipulate matrix coefficients for time dependent simulations
-(void)addFirstOrderTimeModel:(FEMModel * __nonnull)model solution:(FEMSolution * __nonnull)solution element:(Element_t * __nonnull)element massMatrix:(double * __nonnull * __nonnull)massMatrix stiffMatrix:(double * __nonnull * __nonnull)stiffMatrix force:(double * __nonnull)force dt:(double)dt size:(int)n dofs:(int)dofs nodeIndexes:(int * __nonnull)nodeIndexes rows:(int * __nonnull)rows cols:(int * __nonnull)cols timeIntegration:(FEMTimeIntegration * __nonnull)timeIntegration utilities:(FEMUtilities * __nonnull)utilities;
-(void)addSecondOrderTimeModel:(FEMModel * __nonnull)model solution:(FEMSolution * __nonnull)solution element:(Element_t * __nonnull)element massMatrix:(double * __nonnull * __nonnull)massMatrix dampMatrix:(double * __nonnull * __nonnull)dampMatrix stiffMatrix:(double * __nonnull * __nonnull)stiffMatrix force:(double * __nonnull)force dt:(double)dt size:(int)n dofs:(int)dofs nodeIndexes:(int * __nonnull)nodeIndexes rows:(int * __nonnull)rows cols:(int * __nonnull)cols timeIntegration:(FEMTimeIntegration * __nonnull)timeIntegration;
-(void)updateMassMatrixModel:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution localMassMatrix:(double * __nonnull * __nonnull)localMassMatrix element:(Element_t * __nonnull)element numberOfNodes:(int)n dofs:(int)dofs nodeIndexes:(int * __nonnull)nodeIndexes;
-(void)defaultUpdateMass:(double * __nonnull * __nonnull)mass element:(Element_t * __nonnull)element solution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model;
-(void)defaultUpdateComplexMass:(double complex * __nonnull * __nonnull)cmass element:(Element_t * __nonnull)element solution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model;
-(void)defaultUpdateDamp:(double * __nonnull * __nonnull)damp element:(Element_t * __nonnull)element solution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model;
-(void)defaultUpdateComplexDamp:(double complex * __nonnull * __nonnull)cdamp element:(Element_t * __nonnull)element solution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model;
-(void)defaultFirstOrderTime:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution forElement:(Element_t * __nonnull)element realMass:(double * __nonnull * __nonnull)mass realStiff:(double * __nonnull * __nonnull)stiff realForce:(double * __nonnull)force stiffRows:(int * __nonnull)rows stiffCols:(int * __nonnull)cols timeIntegration:(FEMTimeIntegration * __nonnull)timeIntegration utilities:(FEMUtilities * __nonnull)utilities;
-(void)defaultFirstOrderTime:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution forElement:(Element_t * __nonnull)element complexMass:(double complex * __nonnull * __nonnull)cmass complexStiff:(double complex * __nonnull * __nonnull)cstiff complexForce:(double complex * __nonnull)cforce stiffRows:(int * __nonnull)rows stiffCols:(int * __nonnull)cols timeIntegration:(FEMTimeIntegration * __nonnull)timeIntegration utilities:(FEMUtilities * __nonnull)utilities;
-(void)defaultSecondOrderTime:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution forElement:(Element_t * __nonnull)element realMass:(double * __nonnull * __nonnull)mass realDamp:(double * __nonnull * __nonnull)damp realStiff:(double * __nonnull * __nonnull)stiff realForce:(double * __nonnull)force stiffRows:(int * __nonnull)rows stiffCols:(int * __nonnull)cols timeIntegration:(FEMTimeIntegration * __nonnull)timeIntegration;
-(void)defaultSecondOrderTime:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution forElement:(Element_t * __nonnull)element complexMass:(double complex * __nonnull * __nonnull)cmass complexDamp:(double complex * __nonnull * __nonnull)cdamp complexStiff:(double complex * __nonnull * __nonnull)cstiff complexForce:(double complex * __nonnull)cforce stiffRows:(int * __nonnull)rows stiffCols:(int * __nonnull)cols timeIntegration:(FEMTimeIntegration * __nonnull)timeIntegration;


-(void)defaultFirstOrderTimeGlobalModel:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution timeIntegration:(FEMTimeIntegration * __nonnull)timeIntegration utilities:(FEMUtilities * __nonnull)utilities;
-(void)defaultSecondOrderTimeGlobalModel:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution timeIntegration:(FEMTimeIntegration * __nonnull)timeIntegration;

// Update global equations
-(void)updateGlobalEquationsModel:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution element:(Element_t * __nonnull)element localStiffMatrix:(double * __nonnull * __nonnull)localStiffMatrix forceVector:(double * __nonnull)forceVector localForce:(double * __nonnull)localForce size:(int)n dofs:(int)dofs nodeIndexes:(int * __nonnull)nodeIndexes rows:(int * __nonnull)rows cols:(int * __nonnull)cols rotateNT:(BOOL * __nullable)rotateNT crsMatrix:(FEMMatrixCRS * __nonnull)crsMatrix bandMatrix:(FEMMatrixBand * __nonnull)bandMatrix;
-(void)defaultUpdateEquations:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution forElement:(Element_t * __nonnull)element realStiff:(double * __nonnull * __nonnull)stiff realForce:(double * __nonnull)force stiffRows:(int * __nonnull)rows stiffCols:(int * __nonnull)cols crsMatrix:(FEMMatrixCRS * __nonnull)crsMatrix bandMatrix:(FEMMatrixBand * __nonnull)bandMatrix;
-(void)defaultUpdateEquations:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution forElement:(Element_t * __nonnull)element complexStiff:(double complex * __nonnull * __nonnull)cstiff complexForce:(double complex * __nonnull )cforce stiffRows:(int * __nonnull)rows stiffCols:(int * __nonnull)cols crsMatrix:(FEMMatrixCRS * __nonnull)crsMatrix bandMatrix:(FEMMatrixBand * __nonnull)bandMatrix;
-(void)defaultFinishBulkAssemblySolution:(FEMSolution * __nonnull)solution bulkUpdate:(BOOL * __nullable)bulkUpdate;
-(void)defaultFinishAssemblySolution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model timeIntegration:(FEMTimeIntegration * __nonnull)timeIntegration utilities:(FEMUtilities * __nonnull)utilities;
-(void)dirichletBoundaryConditions:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution usingOffset:(int * __nullable)offset offDiaginalMatrix:(BOOL * __nullable)offDiaginalMatrix;

-(void)computeChange:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)aModel isSteadyState:(BOOL)steadyState nsize:(int * __nullable)nsize values:(double * __nullable)values values0:(double * __nullable)values0 sizeValues0:(int * __nullable)sizeValues0;
-(void)iterativeSolveMatrix:(FEMMatrix * __nonnull)matrix result:(double * __nonnull)x rhs:(double * __nonnull)b dimensions:(int * __nullable)ndim solution:(FEMSolution * __nonnull)solution;
-(void)solveSystemMatrix:(FEMMatrix * __nonnull)matrix rhs:(double * __nonnull)b result:(double * __nonnull)x norm:(double * __nonnull)norm dofs:(int)dofs solution:(FEMSolution * __nonnull)solution  model:(FEMModel * __nonnull)model;
-(double)findSolution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model backRorateNT:(BOOL * __nullable)backRorateNT;

// Utility method in matrix-vector multiplication
-(double)stopc:(FEMMatrix * __nonnull)matrix multiplyVector:(double * __nonnull)x righHandSide:(double * __nonnull)b ipar:(int * __nonnull)ipar;

// Activating solution and starting the calculation
-(void)activateSolution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model timeStep:(double)dt transientSimulation:(BOOL)transient;
-(void)solveEquationsModel:(FEMModel * __nonnull)model timeStep:(double * __nonnull)dt transientSimulation:(BOOL)transient coupledMinIteration:(int)coupledMinIter coupleMaxIteration:(int)coupleMaxIter steadyStateReached:(BOOL * __nonnull)steadyStateReached realTimeStep:(int * __nonnull)realTimeStep;

// OpenCL stuff
-(dispatch_queue_t __nonnull)getDispatchQueueAndInfoForDeviceType:(NSString * __nonnull)deviceType;

@end
