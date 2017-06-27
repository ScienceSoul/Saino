//===----------------------------------------------------------------------===//
//  FEMCore.h
//  Saino
//
//  Created by Seddik hakime on 15/04/11.
//  Copyright (c) 2011 ScienceSoul. All rights reserved.
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//===----------------------------------------------------------------------===//

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
@property(nonatomic, assign) bool * _Nullable * _Nullable ntZeroingDone;
@property(nonatomic, assign, nullable) int *boundaryReorder;
@property(nonatomic, assign) int * _Nullable * _Nullable ntElement;
@property(nonatomic, assign) double * _Nullable * _Nullable boundaryNormals;
@property(nonatomic, assign) double * _Nullable * _Nullable boundaryTangent1;
@property(nonatomic, assign) double * _Nullable * _Nullable boundaryTangent2;
@property(nonatomic, strong, nonnull) NSMutableString *normalTangentialName;
@property(nonatomic, strong, nonnull) NSMutableArray <NSNumber *> *outputLevelMask;
@property(nonatomic, assign, getter = isOutputPrefix) BOOL outputPrefix;
@property(nonatomic, assign, getter = isOutputCaller) BOOL outputCaller;
@property(nonatomic, assign) int maxOutputLevel;
@property(nonatomic, assign) int minOutputLevel;
@property(nonatomic, assign) int outputPE;

// This class method retuns an instance (singleton) of FEMCore
+(id _Nonnull)sharedCore;
+(void)selfDestruct;

-(void)deallocation;

-(void)initializeTimeStepInSolution:(FEMSolution * _Nonnull)solution model:(FEMModel * _Nonnull)model;
-(void)initializeToZeroMatrix:(FEMMatrix * _Nonnull)matrix forceVector:(double * _Nonnull)forceVector sizeForceVector:(int)sizeForceVector model:(FEMModel * _Nonnull)model solution:(FEMSolution * _Nonnull)solution ;
-(void)defaultInitializeSolution:(FEMSolution * _Nonnull)solution model:(FEMModel * _Nonnull)model;
-(BOOL)getReal:(FEMModel * _Nonnull)model forElement:(Element_t * _Nullable)element inArray:(NSArray * _Nonnull)array variableName:(NSString * _Nonnull)name buffer:(listBuffer * _Nonnull)result listUtilities:(FEMListUtilities * _Nonnull)listUtil;
-(int)isPElement:(Element_t * _Nonnull)element;
-(void)getNodes:(FEMSolution * _Nonnull)solution model:(FEMModel * _Nonnull)model inElement:(Element_t * _Nonnull)element resultNodes:(Nodes_t * _Nonnull)nodes numberOfNodes:(int * _Nullable)nd mesh:(FEMMesh * _Nullable)mesh;
-(BOOL)isFluxElement:(Element_t * _Nonnull)element mesh:(FEMMesh * _Nonnull)mesh ;
-(int)getElementFamily:(Element_t * _Nonnull)element;
-(int)getElementDofsSolution:(FEMSolution * _Nullable)uSolution model:(FEMModel * _Nonnull)model forElement:(Element_t * _Nonnull)element atIndexes:(int * _Nonnull)indexes disableDiscontinuousGalerkin:(BOOL * _Nullable)disableDiscontinuousGalerkin;
-(int)sgetElementDofsSolution:(FEMSolution * _Nullable)uSolution model:(FEMModel * _Nonnull)model forElement:(Element_t * _Nonnull)element atIndexes:(int * _Nonnull)indexes;
-(int)getNumberOfBubbleDofsElement:(Element_t * _Nonnull)element solution:(FEMSolution * _Nonnull)solution;
-(Element_t * _Nonnull)getActiveElement:(int)t solution:(FEMSolution * _Nonnull)solution model:(FEMModel * _Nonnull)model;
-(int)getNumberOfNodesForElement:(Element_t * _Nonnull)element;
-(int)getBoundaryConditionID:(FEMModel * _Nonnull)model forElement:(Element_t * _Nonnull)element;
-(NSArray * _Nullable)getBoundaryCondition:(FEMModel * _Nonnull)model forElement:(Element_t * _Nonnull)element ;
-(Element_t * _Nonnull)getBoundaryElement:(FEMSolution * _Nonnull)solution atIndex:(int)index ;
-(void)getScalarLocalField:(double * _Nonnull)field sizeField:(int)sizeField name:(NSString * _Nullable)name element:(Element_t * _Nonnull)element solution:(FEMSolution * _Nonnull)solution model:(FEMModel * _Nonnull)model timeStep:(int * _Nullable)tStep;
-(void)getVectorLocalField:(double * _Nonnull * _Nonnull)field size1Field:(int)size1Field size2Field:(int)size2Field name:(NSString * _Nullable)name element:(Element_t * _Nonnull)element solution:(FEMSolution * _Nonnull)solution model:(FEMModel * _Nonnull)model timeStep:(int * _Nullable)tStep;
-(int * _Nonnull * _Nonnull)getEdgeMap:(int)elementFamily mapSize:(int * _Nullable)mapSize;
-(void)getBoundaryIndexes:(FEMMesh * _Nonnull)mesh forBoundaryElement:(Element_t * _Nonnull)element withParentElement:(Element_t * _Nonnull)parent resultVector:(int * _Nonnull)indexes sizeVector:(int)size indexSize:(int * _Nonnull)indexSize;
-(int)getBodyForceIDForElement:(Element_t * _Nonnull)element model:(FEMModel * _Nonnull)model;
-(int)getMaterialIDForElement:(Element_t * _Nonnull)element model:(FEMModel * _Nonnull)model;
-(int)getEquationIDForElement:(Element_t * _Nonnull)element model:(FEMModel * _Nonnull)model;
-(BOOL)getParentMaterialProperty:(NSString * _Nonnull)name forElement:(Element_t * _Nonnull)element parentElement:(Element_t * _Nullable)parentElement model:(FEMModel * _Nonnull)model listUtilities:(FEMListUtilities * _Nonnull)listUtilities buffer:(listBuffer * _Nonnull)result;
-(BOOL)isActiveBoundaryElement:(Element_t * _Nonnull)element inSolution:(FEMSolution * _Nonnull)solution model:(FEMModel * _Nonnull)model;
-(void)checkNormalTangential:(FEMModel * _Nonnull)model inSolution:(FEMSolution * _Nonnull)solution forElementNumber:(int)elno numberofNodes:(int)n atIndexes:(int * _Nonnull)indexes atBoundary:(int)bc variableName:(NSMutableString * _Nonnull)name orderOfDofs:(int)dof activeCondition:(BOOL)conditional conditionName:(NSString * _Nonnull)condName permutationOffset:(int)offset;
-(void)zeroTheNumberOfRows:(int)n inSolutionMatrix:(FEMSolution * _Nonnull)solution;
-(void)setMatrixElementForSolution:(FEMSolution * _Nonnull)solution atIndex:(int)i andIndex:(int)j value:(double)value;
-(void)addToMatrixElementForSolution:(FEMSolution * _Nonnull)solution atIndex:(int)i andIndex:(int)j value:(double)value;
-(void)localBoundaryIntegral:(FEMModel * _Nonnull)model inSolution:(FEMSolution * _Nonnull)solution atBoundary:(NSArray * _Nonnull)bc forElement:(Element_t * _Nonnull)element withNumberOfNodes:(int)nd andParent:(Element_t * _Nonnull)parent withNumberOfNodes:(int)np boundaryName:(NSString * _Nonnull)name functionIntegral:(double * _Nonnull)integral;
-(void)localBoundaryBDOFs:(FEMModel * _Nonnull)model inSolution:(FEMSolution * _Nonnull)solution atBoundary:(NSArray * _Nonnull)bc forElement:(Element_t * _Nonnull)element withNumberOfNodes:(int)nd boundaryName:(NSMutableString * _Nonnull)name resultMatrix:(double * _Nonnull * _Nonnull)stiff resultVector:(double * _Nonnull)force;
-(void)solveWithLapackMatrix:(double * _Nonnull)a andVector:(double * _Nonnull)x size:(int)n leadingDimension:(int)lda;
-(void)solveLinearSystemWithMatrix:(double * _Nonnull * _Nonnull)a andVector:(double * _Nonnull)x size:(int)n leadingDimension:(int)lda;
-(void)setDirichletBoundaries:(FEMModel * _Nonnull)model inSolution:(FEMSolution * _Nonnull)solution variableName:(NSMutableString * _Nonnull)name orderOfDofs:(int)dof permutationOffset:(int * _Nullable)offset offDiaginalMatrix:(BOOL * _Nullable)offDiaginalMatrix listUtilites:(FEMListUtilities * _Nonnull)listUtilites crsMatrix:(FEMMatrixCRS * _Nonnull)crsMatrix bandMatrix:(FEMMatrixBand * _Nonnull)bandMatrix;
-(void)scaleLinearSystem:(FEMSolution * _Nonnull)solution matrix:(FEMMatrix * _Nonnull)matrix rhs:(double * _Nullable)b result:(double * _Nullable)x diagScaling:(double * _Nullable)diagScaling applyScaling:(BOOL * _Nullable)applyScaling rhsScaling:(BOOL * _Nullable)rhsScaling;
-(void)backScaleLinearSystem:(FEMSolution * _Nonnull)solution matrix:(FEMMatrix * _Nonnull)matrix rhs:(double * _Nullable)b result:(double * _Nullable)x diagScaling:(double * _Nullable)diagScaling sizeOFDiagScaling:(int * _Nullable)sizeOfDiagScaling;
-(void)matrixVectorMultplyInSolution:(FEMSolution * _Nonnull)solution multiplyVector:(double * _Nonnull)u resultVector:(double * _Nonnull)v;
-(void)matrixVectorMultplyInMatrix:(FEMMatrix * _Nonnull)matrix multiplyVector:(double * _Nonnull)u resultVector:(double * _Nonnull)v;
-(void)invalidateVariableInTopMesh:(NSArray * _Nonnull)topMesh primaryMesh:(FEMMesh * _Nonnull)primaryMesh name:(NSString * _Nonnull)name model:(FEMModel * _Nonnull)model;
-(void)getPassiveBoundaryAtIndex:(int)bcID model:(FEMModel * _Nonnull)model mesh:(FEMMesh * _Nonnull)mesh solution:(FEMSolution * _Nonnull)solution;
-(void)computeNodalWeightsInSolution:(FEMSolution * _Nonnull)solution model:(FEMModel * _Nonnull)model weightBoundary:(BOOL)weightBoundary perm:(int * _Nullable)perm sizePerm:(int * _Nullable)sizePerm variableName:(NSString * _Nullable)variableName;
-(void)nsCondensateStiff:(double * _Nonnull * _Nonnull)stiff force:(double * _Nonnull)force numberOfNodes:(int)n numberOfBubbles:(int)nb dimension:(int)dim force1:(double * _Nonnull)force1;
-(void)condensateStiff:(double * _Nonnull * _Nonnull)stiff force:(double * _Nonnull)force numberOfNodes:(int)n force1:(double * _Nullable)force1;

// Manipulate matrix coefficients for time dependent simulations
-(void)addFirstOrderTimeModel:(FEMModel * _Nonnull)model solution:(FEMSolution * _Nonnull)solution element:(Element_t * _Nonnull)element massMatrix:(double * _Nonnull * _Nonnull)massMatrix stiffMatrix:(double * _Nonnull * _Nonnull)stiffMatrix force:(double * _Nonnull)force dt:(double)dt size:(int)n dofs:(int)dofs nodeIndexes:(int * _Nonnull)nodeIndexes rows:(int * _Nonnull)rows cols:(int * _Nonnull)cols timeIntegration:(FEMTimeIntegration * _Nonnull)timeIntegration utilities:(FEMUtilities * _Nonnull)utilities;
-(void)addSecondOrderTimeModel:(FEMModel * _Nonnull)model solution:(FEMSolution * _Nonnull)solution element:(Element_t * _Nonnull)element massMatrix:(double * _Nonnull * _Nonnull)massMatrix dampMatrix:(double * _Nonnull * _Nonnull)dampMatrix stiffMatrix:(double * _Nonnull * _Nonnull)stiffMatrix force:(double * _Nonnull)force dt:(double)dt size:(int)n dofs:(int)dofs nodeIndexes:(int * _Nonnull)nodeIndexes rows:(int * _Nonnull)rows cols:(int * _Nonnull)cols timeIntegration:(FEMTimeIntegration * _Nonnull)timeIntegration;
-(void)updateMassMatrixModel:(FEMModel * _Nonnull)model inSolution:(FEMSolution * _Nonnull)solution localMassMatrix:(double * _Nonnull * _Nonnull)localMassMatrix element:(Element_t * _Nonnull)element numberOfNodes:(int)n dofs:(int)dofs nodeIndexes:(int * _Nonnull)nodeIndexes;
-(void)defaultUpdateMass:(double * _Nonnull * _Nonnull)mass element:(Element_t * _Nonnull)element solution:(FEMSolution * _Nonnull)solution model:(FEMModel * _Nonnull)model;
-(void)defaultUpdateComplexMass:(double complex * _Nonnull * _Nonnull)cmass element:(Element_t * _Nonnull)element solution:(FEMSolution * _Nonnull)solution model:(FEMModel * _Nonnull)model;
-(void)defaultUpdateDamp:(double * _Nonnull * _Nonnull)damp element:(Element_t * _Nonnull)element solution:(FEMSolution * _Nonnull)solution model:(FEMModel * _Nonnull)model;
-(void)defaultUpdateComplexDamp:(double complex * _Nonnull * _Nonnull)cdamp element:(Element_t * _Nonnull)element solution:(FEMSolution * _Nonnull)solution model:(FEMModel * _Nonnull)model;
-(void)defaultFirstOrderTime:(FEMModel * _Nonnull)model inSolution:(FEMSolution * _Nonnull)solution forElement:(Element_t * _Nonnull)element realMass:(double * _Nonnull * _Nonnull)mass realStiff:(double * _Nonnull * _Nonnull)stiff realForce:(double * _Nonnull)force stiffRows:(int * _Nonnull)rows stiffCols:(int * _Nonnull)cols timeIntegration:(FEMTimeIntegration * _Nonnull)timeIntegration utilities:(FEMUtilities * _Nonnull)utilities;
-(void)defaultFirstOrderTime:(FEMModel * _Nonnull)model inSolution:(FEMSolution * _Nonnull)solution forElement:(Element_t * _Nonnull)element complexMass:(double complex * _Nonnull * _Nonnull)cmass complexStiff:(double complex * _Nonnull * _Nonnull)cstiff complexForce:(double complex * _Nonnull)cforce stiffRows:(int * _Nonnull)rows stiffCols:(int * _Nonnull)cols timeIntegration:(FEMTimeIntegration * _Nonnull)timeIntegration utilities:(FEMUtilities * _Nonnull)utilities;
-(void)defaultSecondOrderTime:(FEMModel * _Nonnull)model inSolution:(FEMSolution * _Nonnull)solution forElement:(Element_t * _Nonnull)element realMass:(double * _Nonnull * _Nonnull)mass realDamp:(double * _Nonnull * _Nonnull)damp realStiff:(double * _Nonnull * _Nonnull)stiff realForce:(double * _Nonnull)force stiffRows:(int * _Nonnull)rows stiffCols:(int * _Nonnull)cols timeIntegration:(FEMTimeIntegration * _Nonnull)timeIntegration;
-(void)defaultSecondOrderTime:(FEMModel * _Nonnull)model inSolution:(FEMSolution * _Nonnull)solution forElement:(Element_t * _Nonnull)element complexMass:(double complex * _Nonnull * _Nonnull)cmass complexDamp:(double complex * _Nonnull * _Nonnull)cdamp complexStiff:(double complex * _Nonnull * _Nonnull)cstiff complexForce:(double complex * _Nonnull)cforce stiffRows:(int * _Nonnull)rows stiffCols:(int * _Nonnull)cols timeIntegration:(FEMTimeIntegration * _Nonnull)timeIntegration;


-(void)defaultFirstOrderTimeGlobalModel:(FEMModel * _Nonnull)model inSolution:(FEMSolution * _Nonnull)solution timeIntegration:(FEMTimeIntegration * _Nonnull)timeIntegration utilities:(FEMUtilities * _Nonnull)utilities;
-(void)defaultSecondOrderTimeGlobalModel:(FEMModel * _Nonnull)model inSolution:(FEMSolution * _Nonnull)solution timeIntegration:(FEMTimeIntegration * _Nonnull)timeIntegration;

// Update global equations
-(void)updateGlobalEquationsModel:(FEMModel * _Nonnull)model inSolution:(FEMSolution * _Nonnull)solution element:(Element_t * _Nonnull)element localStiffMatrix:(double * _Nonnull * _Nonnull)localStiffMatrix forceVector:(double * _Nonnull)forceVector localForce:(double * _Nonnull)localForce size:(int)n dofs:(int)dofs nodeIndexes:(int * _Nonnull)nodeIndexes rows:(int * _Nonnull)rows cols:(int * _Nonnull)cols rotateNT:(BOOL * _Nullable)rotateNT crsMatrix:(FEMMatrixCRS * _Nonnull)crsMatrix bandMatrix:(FEMMatrixBand * _Nonnull)bandMatrix;
-(void)defaultUpdateEquations:(FEMModel * _Nonnull)model inSolution:(FEMSolution * _Nonnull)solution forElement:(Element_t * _Nonnull)element realStiff:(double * _Nonnull * _Nonnull)stiff realForce:(double * _Nonnull)force stiffRows:(int * _Nonnull)rows stiffCols:(int * _Nonnull)cols crsMatrix:(FEMMatrixCRS * _Nonnull)crsMatrix bandMatrix:(FEMMatrixBand * _Nonnull)bandMatrix;
-(void)defaultUpdateEquations:(FEMModel * _Nonnull)model inSolution:(FEMSolution * _Nonnull)solution forElement:(Element_t * _Nonnull)element complexStiff:(double complex * _Nonnull * _Nonnull)cstiff complexForce:(double complex * _Nonnull )cforce stiffRows:(int * _Nonnull)rows stiffCols:(int * _Nonnull)cols crsMatrix:(FEMMatrixCRS * _Nonnull)crsMatrix bandMatrix:(FEMMatrixBand * _Nonnull)bandMatrix;
-(void)defaultFinishBulkAssemblySolution:(FEMSolution * _Nonnull)solution bulkUpdate:(BOOL * _Nullable)bulkUpdate;
-(void)defaultFinishAssemblySolution:(FEMSolution * _Nonnull)solution model:(FEMModel * _Nonnull)model timeIntegration:(FEMTimeIntegration * _Nonnull)timeIntegration utilities:(FEMUtilities * _Nonnull)utilities;
-(void)dirichletBoundaryConditions:(FEMModel * _Nonnull)model inSolution:(FEMSolution * _Nonnull)solution usingOffset:(int * _Nullable)offset offDiaginalMatrix:(BOOL * _Nullable)offDiaginalMatrix;

-(void)computeChange:(FEMSolution * _Nonnull)solution model:(FEMModel * _Nonnull)aModel isSteadyState:(BOOL)steadyState nsize:(int * _Nullable)nsize values:(double * _Nullable)values values0:(double * _Nullable)values0 sizeValues0:(int * _Nullable)sizeValues0;
-(void)iterativeSolveMatrix:(FEMMatrix * _Nonnull)matrix result:(double * _Nonnull)x rhs:(double * _Nonnull)b dimensions:(int * _Nullable)ndim solution:(FEMSolution * _Nonnull)solution;
-(void)solveSystemMatrix:(FEMMatrix * _Nonnull)matrix rhs:(double * _Nonnull)b result:(double * _Nonnull)x norm:(double * _Nonnull)norm dofs:(int)dofs solution:(FEMSolution * _Nonnull)solution  model:(FEMModel * _Nonnull)model;
-(double)findSolution:(FEMSolution * _Nonnull)solution model:(FEMModel * _Nonnull)model backRorateNT:(BOOL * _Nullable)backRorateNT;

// Utility method in matrix-vector multiplication
-(double)stopc:(FEMMatrix * _Nonnull)matrix multiplyVector:(double * _Nonnull)x righHandSide:(double * _Nonnull)b ipar:(int * _Nonnull)ipar;

// Activating solution and starting the calculation
-(void)activateSolution:(FEMSolution * _Nonnull)solution model:(FEMModel * _Nonnull)model timeStep:(double)dt transientSimulation:(BOOL)transient;
-(void)solveEquationsModel:(FEMModel * _Nonnull)model timeStep:(double * _Nonnull)dt transientSimulation:(BOOL)transient coupledMinIteration:(int)coupledMinIter coupleMaxIteration:(int)coupleMaxIter steadyStateReached:(BOOL * _Nonnull)steadyStateReached realTimeStep:(int * _Nonnull)realTimeStep;

// OpenCL stuff
-(dispatch_queue_t _Nonnull)getDispatchQueueAndInfoForDeviceType:(NSString * _Nonnull)deviceType;

@end
