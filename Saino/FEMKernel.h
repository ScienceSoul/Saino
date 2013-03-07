//
//  FEMKernel.h
//  Saino
//
//  Created by Seddik hakime on 15/04/11.
//  Copyright 2011 ScienceSoul. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMModel.h"
#import "FEMSolution.h"
#import "FEMMesh.h"
#import "FEMHUTIter.h"
#import "FEMPrecondition.h"
#import "FEMParallelMPI.h"
#import "FEMTimeIntegration.h"
#import "FEMMatrixCRS.h"
#import "FEMMatrixBand.h"
#import "FEMElementDescription.h"
#import "FEMNumericIntegration.h"
#import "FEMUtilities.h"

#import "Constructors.h"

@interface FEMKernel : NSObject {
    
    int _coordinateSystemDimension;
    int _normalTangentialNumberOfNodes;
    int _sizeBoundaryReorder;
    int _size1boundaryNormals;
    int _size2boundaryNormals;
    int _size1boundaryTangent1;
    int _size2boundaryTangent1;
    int _size1boundaryTangent2;
    int _size2boundaryTangent2;
    int *_boundaryReorder;
    double **_boundaryNormals;
    double **_boundaryTangent1;
    double **_boundaryTangent2;
}

@property(nonatomic, assign) int coordinateSystemDimension;
@property(nonatomic, assign) int normalTangentialNumberOfNodes;
@property(nonatomic, assign) int sizeBoundaryReorder;
@property(nonatomic, assign) int size1boundaryNormals;
@property(nonatomic, assign) int size2boundaryNormals;
@property(nonatomic, assign) int size1boundaryTangent1;
@property(nonatomic, assign) int size2boundaryTangent1;
@property(nonatomic, assign) int size1boundaryTangent2;
@property(nonatomic, assign) int size2boundaryTangent2;
@property(nonatomic, assign) int *boundaryReorder;
@property(nonatomic, assign) double **boundaryNormals;
@property(nonatomic, assign) double **boundaryTangent1;
@property(nonatomic, assign) double **boundaryTangent2;

// This class method retuns an instance (singleton) of FEMKernel
+(id)sharedKernel;

-(void)deallocation;
-(BOOL)getReal:(FEMModel *)model forElement:(Element_t *)element inList:(NSArray *)list variableName:(NSString *)name buffer:(listBuffer *)result;
-(int)isPElement:(Element_t *)element;
-(void)getNodes:(FEMSolution *)solution inElement:(Element_t *)element resultNodes:(Nodes_t *)nodes numberOfNodes:(int)nd;
-(int)getElementFamily:(Element_t *)element;
-(int)getElementDofs:(FEMSolution *)solution forElement:(Element_t *)element atIndexes:(int *)indexes;
-(int)sgetElementDofs:(FEMSolution *)solution forElement:(Element_t *)element atIndexes:(int *)indexes;
-(int)getNumberOfNodesForElement:(Element_t *)element;
-(int)getBoundaryConditionID:(FEMModel *)model forElement:(Element_t *)element;
-(NSArray *)getBoundaryCondition:(FEMModel *)model forElement:(Element_t *)element;
-(Element_t *)getBoundaryElement:(FEMSolution *)solution atIndex:(int)index;
-(BOOL)isActiveElement:(Element_t *)element inSolution:(FEMSolution *)solution;
-(void)checkNormalTangential:(FEMModel *)model inSolution:(FEMSolution *)solution forElementNumber:(int)elno numberofNodes:(int)n atIndexes:(int *)indexes atBoundary:(int)bc variableName:(NSMutableString *)name orderOfDofs:(int)dof activeCondition:(BOOL)conditional conditionName:(NSString *)condName permutationOffset:(int)offset;
-(void)getBoundaryIndexes:(FEMMesh *)mesh forBoundaryElement:(Element_t *)element withParentElement:(Element_t *)parent resultVector:(int *)indexes resultSize:(int *)indSize;
-(int **)getEdgeMap:(int)elementFamily;
-(void)zeroTheNumberOfRows:(int)n inSolutionMatrix:(FEMSolution *)solution;
-(void)setMatrixElementForSolution:(FEMSolution *)solution atIndex:(int)i andIndex:(int)j value:(double)value;
-(void)addToMatrixElementForSolution:(FEMSolution *)solution atIndex:(int)i andIndex:(int)j value:(double)value;
-(void)localBoundaryIntegral:(FEMModel *)model inSolution:(FEMSolution *)solution atBoundary:(NSArray *)bc forElement:(Element_t *)element withNumberOfNodes:(int)nd andParent:(Element_t *)parent withNumberOfNodes:(int)np boundaryName:(NSMutableString *)name functionIntegral:(double)integral;
-(void)localBoundaryBDOFs:(FEMModel *)model inSolution:(FEMSolution *)solution atBoundary:(NSArray *)bc forElement:(Element_t *)element withNumberOfNodes:(int)nd boundaryName:(NSMutableString *)name resultMatrix:(double **)stiff resultVector:(double *)force;
-(void)solveWithLapackMatrix:(double *)a andVector:(double *)x size:(int)n leadingDimension:(int)lda;
-(void)solveLinearSystemWithMatrix:(double **)a andVector:(double *)x size:(int)n leadingDimension:(int)lda;
-(void)setNodalLoads:(FEMModel *)model inSolution:(FEMSolution *)solution variableName:(NSMutableString *)name orderOfDofs:(int)dof;
-(void)setDirichletBoundaries:(FEMModel *)model inSolution:(FEMSolution *)solution variableName:(NSMutableString *)name orderOfDofs:(int)dof permutationOffset:(int *)offset;
-(void)scaleLinearSystem:(FEMSolution *)solution scaleSolutionMatrixRHS:(BOOL)scaleSolutionMatrixRHS scaleSolutionVariable:(BOOL)scaleSolutionVariable diagScaling:(double *)diagScaling applyScaling:(BOOL *)applyScaling rhsScaling:(BOOL *)rhsScaling;
-(void)backScaleLinearSystem:(FEMSolution *)solution scaleSolutionMatrixRHS:(BOOL)scaleSolutionMatrixRHS scaleSolutionVariable:(BOOL)scaleSolutionVariable diagScaling:(double *)diagScaling sizeOFDiagScaling:(int *)sizeOfDiagScaling;
-(void)matrixVectorMultplyInSolution:(FEMSolution *)solution multiplyVector:(double *)u resultVector:(double *)v;
-(void)matrixVectorMultplyInMatrix:(FEMMatrix *)matrix multiplyVector:(double *)u resultVector:(double *)v;
-(void)invalidateVariableInTopMesh:(NSArray *)topMesh primaryMesh:(FEMMesh *)primaryMesh name:(NSString *)name model:(FEMModel *)model;

-(void)defaultFirstOrderTime:(FEMModel *)model inSolution:(FEMSolution *)solution forElement:(Element_t *)element realMass:(double **)mass realStiff:(double **)stiff realForce:(double *)force stiffRows:(int *)rows stiffCols:(int *)cols;
-(void)defaultFirstOrderTime:(FEMModel *)model inSolution:(FEMSolution *)solution forElement:(Element_t *)element complexMass:(double complex **)cmass complexStiff:(double complex **)cstiff complexForce:(double complex *)cforce stiffRows:(int *)rows stiffCols:(int *)cols;

-(void)defaultUpdateEquations:(FEMModel *)model inSolution:(FEMSolution *)solution forElement:(Element_t *)element realStiff:(double **)stiff realForce:(double *)force stiffRows:(int *)rows stiffCols:(int *)cols requestBulkUpdate:(BOOL *)bulkUpdate;
-(void)defaultUpdateEquations:(FEMModel *)model inSolution:(FEMSolution *)solution forElement:(Element_t *)element complexStiff:(double **)cstiff complexForce:(double *)cforce stiffRows:(int *)rows stiffCols:(int *)cols requestBulkUpdate:(BOOL *)bulkUpdate;

-(void)dirichletBoundaryConditions:(FEMModel *)model inSolution:(FEMSolution *)solution usingOffset:(int *)offset;

-(void)iterativeSolve:(FEMSolution *)solution;
-(double)findSolution:(FEMSolution *)solution model:(FEMModel *)aModel;
-(double)stopc:(FEMSolution *)solution multiplyVector:(double *)x righHandSide:(double *)b ipar:(int *)ipar;

@end
