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

#import "Constructors.h"


@interface FEMKernel : NSObject {
@private
    int *indexStore;
  
    
}

-(int)isPElement:(Element_t *)element;
-(int)getElementFamily:(Element_t *)element;
-(int)getElementDofs:(FEMSolution *)solution forElement:(Element_t *)element atIndexes:(int *)indexes;
-(int)getNumberOfNodesForElement:(Element_t *)element;
-(int)getBoundaryConditionID:(FEMModel *)model forElement:(Element_t *)element;
-(NSArray *)getBoundaryCondition:(FEMModel *)model forElement:(Element_t *)element;

-(void)defaultFirstOrderTime:(FEMModel *)model inSolution:(FEMSolution *)solution forElement:(Element_t *)element realMass:(double **)mass realStiff:(double **)stiff realForce:(double *)force stiffRows:(int *)rows stiffCols:(int *)cols;
-(void)defaultFirstOrderTime:(FEMModel *)model inSolution:(FEMSolution *)solution forElement:(Element_t *)element complexMass:(double complex **)cmass complexStiff:(double complex **)cstiff complexForce:(double complex *)cforce stiffRows:(int *)rows stiffCols:(int *)cols;

-(void)defaultUpdateEquations:(FEMModel *)model inSolution:(FEMSolution *)solution forElement:(Element_t *)element realStiff:(double **)stiff realForce:(double *)force stiffRows:(int *)rows stiffCols:(int *)cols requestBulkUpdate:(BOOL *)bulkUpdate;
-(void)defaultUpdateEquations:(FEMModel *)model inSolution:(FEMSolution *)solution forElement:(Element_t *)element complexStiff:(double **)cstiff complexForce:(double *)cforce stiffRows:(int *)rows stiffCols:(int *)cols requestBulkUpdate:(BOOL *)bulkUpdate;

-(void)boundaryConditionsDefaultDirichlet:(FEMModel *)model toSolution:(FEMSolution *)solution usingOffset:(int)offset;

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
