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
@private
    int *indexStore, sizeIndexStore;
    double **kernStiff, *kernWork;     // kernStiff(maxElementDofs)(maxElementDofs), kernWork(maxElementDofs)
    int *g_Ind, *l_Ind;                // g_Ind(maxElementDofs), l_Ind(maxElementDofs)
    int size1kernStiff, size2kernStiff, sizekernWork, size_g_Ind, size_l_Ind;
    
    int **lineEM;
    int **triangleEM;
    int **quadEM;
    int **tetraEM;
    int **prismEM;
    int **wedgeEM;
    int **brickEM;
    
    BOOL initialized[8];
  
    
}

-(void)deallocation;
-(BOOL)getReal:(FEMModel *)model forElement:(Element_t *)element inList:(NSArray *)list variableName:(NSString *)name resultArray:(double *)result;
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
-(void)checkNormalTangentiality:(FEMModel *)model inSolution:(FEMSolution *)solution forElementNumber:(int)elno numberofNodes:(int)n atIndexes:(int *)indexes atBoundary:(int)bc variableName:(NSMutableString *)name orderOfDofs:(int)dof activeCondition:(BOOL)conditional conditionName:(NSString *)condName permutationOffset:(int)offset;
-(void)getBoundaryIndexes:(FEMMesh *)mesh forBoundaryElement:(Element_t *)element withParentElement:(Element_t *)parent resultVector:(int *)indexes resultSize:(int)indSize;
-(int **)getEdgeMap:(int)elementFamily;
-(void)zeroTheNumberOfRows:(int)n inSolutionMatrix:(FEMSolution *)solution;
-(void)setMatrixElement:(FEMSolution *)solution: (int)i: (int)j: (double)value;
-(void)addToMatrixElement:(FEMSolution *)solution: (int)i: (int)j: (double)value;
-(void)localBoundaryIntegral:(FEMModel *)model inSolution:(FEMSolution *)solution atBoundary:(NSArray *)bc forElement:(Element_t *)element withNumberOfNodes:(int)nd andParent:(Element_t *)parent withNumberOfNodes:(int)np boundaryName:(NSMutableString *)name functionIntegral:(double)integral;
-(void)localBoundaryBDOFs:(FEMModel *)model inSolution:(FEMSolution *)solution atBoundary:(NSArray *)bc forElement:(Element_t *)element withNumberOfNodes:(int)nd boundaryName:(NSMutableString *)name resultMatrix:(double **)stiff resultVector:(double *)force;
-(void)solveWithLapackMatrix:(double *)a andVector:(double *)x size:(int)n;
-(void)solveLinearSystemWithMatrix:(double **)a andVector:(double *)x size:(int)n;
-(void)setNodalLoads:(FEMModel *)model inSolution:(FEMSolution *)solution variableName:(NSMutableString *)name orderOfDofs:(int)dof;
-(void)setDirichletBoundaries:(FEMModel *)model inSolution:(FEMSolution *)solution variableName:(NSMutableString *)name orderOfDofs:(int)dof permutationOffset:(int *)offset;

-(void)defaultFirstOrderTime:(FEMModel *)model inSolution:(FEMSolution *)solution forElement:(Element_t *)element realMass:(double **)mass realStiff:(double **)stiff realForce:(double *)force stiffRows:(int *)rows stiffCols:(int *)cols;
-(void)defaultFirstOrderTime:(FEMModel *)model inSolution:(FEMSolution *)solution forElement:(Element_t *)element complexMass:(double complex **)cmass complexStiff:(double complex **)cstiff complexForce:(double complex *)cforce stiffRows:(int *)rows stiffCols:(int *)cols;

-(void)defaultUpdateEquations:(FEMModel *)model inSolution:(FEMSolution *)solution forElement:(Element_t *)element realStiff:(double **)stiff realForce:(double *)force stiffRows:(int *)rows stiffCols:(int *)cols requestBulkUpdate:(BOOL *)bulkUpdate;
-(void)defaultUpdateEquations:(FEMModel *)model inSolution:(FEMSolution *)solution forElement:(Element_t *)element complexStiff:(double **)cstiff complexForce:(double *)cforce stiffRows:(int *)rows stiffCols:(int *)cols requestBulkUpdate:(BOOL *)bulkUpdate;

-(void)dirichletBoundaryConditions:(FEMModel *)model inSolution:(FEMSolution *)solution usingOffset:(int *)offset;

-(double)findSolution: (FEMSolution *)solution;
-(double)stopc:(FEMSolution *)solution: (double *)x: (double *)b: (double *)r: (int *)ipar;












// One-dimensional tables memory allocation methods
-(BOOL)allocateIntegersVector:(int *)Vector from:(long) nl to:(long) nh;
-(BOOL)allocateFloatsVector:(float *)Vector from:(long) nl to:(long) nh;
-(BOOL)allocateDoublesVector:(double *)Vector from:(long) nl to:(long) nh;

// Two-dimensional tables memory allocation methods
-(BOOL)allocateIntegersMatrix:(int **)Matrix fromRow:(long) nrl toRow:(long) nrh fromColumn:(long) ncl toColumn:(long) nch;
-(BOOL)allocateFloatsMatrix:(float **)Matrix fromRow:(long) nrl toRow:(long) nrh fromColumn:(long) ncl toColumn:(long) nch;
-(BOOL)allocateDoublesMatrix:(double **)Matrix fromRow:(long) nrl toRow:(long) nrh fromColumn:(long) ncl toColumn:(long) nch;

// Three-dimensional tables memory allocation methods
-(BOOL)allocateIntegersTensor:(int ***)Tensor fromRow:(long) nrl toRow:(long) nrh fromColumn:(long) ncl toColumn:(long) nch fromDepth:(long) ndl toDepth:(long) ndh;
-(BOOL)allocateFloatsTensor:(float ***)Tensor fromRow:(long) nrl toRow:(long) nrh fromColumn:(long) ncl toColumn:(long) nch fromDepth:(long) ndl toDepth:(long) ndh;
-(BOOL)allocateDoublesTensor:(double ***)Tensor fromRow:(long) nrl toRow:(long) nrh fromColumn:(long) ncl toColumn:(long) nch fromDepth:(long) ndl toDepth:(long) ndh;

// Memory deallocation methods for one-dimensional, two-dimensional and three-dimensional tables
-(void)freeIntegersVector:(int *)Vector from:(long) nl to:(long) nh;
-(void)freeFloatsVector:(float *)Vector from:(long) nl to:(long) nh;
-(void)freeDoublesVector:(double *)Vector from:(long) nl to:(long) nh;

-(void)freeIntegersMatrix:(int **)Matrix fromRow:(long) nrl toRow:(long) nrh fromColumn:(long) ncl toColumn:(long) nch;
-(void)freeFloatsMatrix:(float **)Matrix fromRow:(long) nrl toRow:(long) nrh fromColumn:(long) ncl toColumn:(long) nch;
-(void)freeDoublesMatrix:(double **)Matrix fromRow:(long) nrl toRow:(long) nrh fromColumn:(long) ncl toColumn:(long) nch;

-(void)freeIntergersTensor:(int ***)Tensor fromRow:(long) nrl toRow:(long) nrh fromColumn:(long) ncl toColumn:(long) nch fromDepth:(long) ndl toDepth:(long) ndh;
-(void)freeFloatsTensor:(float ***)Tensor fromRow:(long) nrl toRow:(long) nrh fromColumn:(long) ncl toColumn:(long) nch fromDepth:(long) ndl toDepth:(long) ndh;
-(void)freeDoublesTensor:(double ***)Tensor fromRow:(long) nrl toRow:(long) nrh fromColumn:(long) ncl toColumn:(long) nch fromDepth:(long) ndl toDepth:(long) ndh;
 

@end
