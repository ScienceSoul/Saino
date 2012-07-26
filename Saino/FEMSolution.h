//
//  FEMSolution.h
//  Saino
//
//  Created by Hakime Seddik on 18/01/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMMesh.h"
#import "FEMMatrix.h"

#import "Constructors.h"
#import "memory.h"
#import "Utils.h"

@interface FEMSolution : NSObject {
    
    NSMutableDictionary *solutionInfo;
    NSString *normalTangentialName;
    
    int simulationID;
    
    int normalTangentialNOFNodes;
    int coordinateSystemDimension;
    
    int **ntElement;
    BOOL **ntZeroingDone;
    int *boundaryReorder;
    double **boundaryNormals, **boundaryTangent1, **boundaryTangent2;
    
    int sizeBoundaryReorder, size1boundaryNormals, size2boundaryNormals;
    
    int timeOrder, doneTime, order, nofEigenValues;
    int lineBeforeProc, linAfterProc;
    double alpha, beta, dt;
    
    int solutionComputeWhen;
    
    int numberOfActiveElements;
    int *acticeElements, *defDofs;
    int sizeDefDofs;
    
    ValueList_t values;
    FEMMatrix *_matrix;
    Variable_t variable;
    Variable_t *exportedVariable;
    FEMMesh *_mesh;
    
}

@property(nonatomic, strong) FEMMatrix *matrix;
@property(nonatomic, strong) FEMMesh *mesh;

#pragma mark Setters and getters for FEMSolution variables

-(id)solutionInfoForKey:(NSString *)key;
-(NSString *)normalTangentialName;
-(int)simulationID;
-(int)normalTangentialNOFNodes;
-(int)coordinateSystemDimension;
-(int)ntElement:(int)i: (int)j;
-(BOOL)ntZeroingDone:(int)i: (int)j;
-(int)boundaryReorder:(int)i;
-(double)boundaryNormals:(int)i: (int)j;
-(double)boundaryTangent1:(int)i: (int)j;
-(double)boundaryTangent2:(int)i: (int)j;
-(int)size1boundaryNormals;
-(int)size2boundaryNormals;
-(int)timeOrder;
-(int)doneTime;
-(int)order;
-(int)nofEigenValues;
-(double)alpha;
-(double)beta;
-(double)dt;
-(int)defDofs:(int)i;
-(int)sizeOfDefDofs;

-(void)setNormalTangentialName:(NSString *)string;
-(void)setSimulationID:(int)i;
-(void)setNormalTangentialNOFNodes:(int)i;
-(void)setCoordinateSystemDimension:(int)i;
-(void)setNtElement:(int)i: (int)j: (int)n;
-(void)setNtZeroingDone:(int)i: (int)j: (BOOL)n;
-(void)setBoundaryReorder:(int)i: (int)n;
-(void)setBoundaryNormals:(int)i: (int)j: (double)n;
-(void)setBoundaryTangent1:(int)i: (int)j: (double)n;
-(void)setBoundaryTangent2:(int)i: (int)j: (double)n;
-(void)setSize1boundaryNormals:(int)n;
-(void)setSize2boundaryNormals:(int)n;
-(void)setTimeOrder:(int)n;
-(void)setDoneTime:(int)n;
-(void)setOrder:(int)n;
-(void)setAlpha:(double)n;
-(void)setBeta:(double)n;
-(void)setDt:(double)n;
-(void)setDefDofs:(int)i: (int)n;
-(void)setSizeOfDefDofs:(int)n;

#pragma mark Setters and getters for Variable structure

-(char *)variableName;
-(int)variableDofs;
-(int)variablePerm:(int)i;
-(double)variableNorm;
-(double)variablePrevNorm;
-(double)variableNonLinChange;
-(double)variableSteadyChange;
-(int)variableNonLinConverged;
-(int)variableSteadyConverged;
-(int)variableNonLinIter;
-(double)variableValues:(int)i;
-(double)variablePrevValues:(int)i: (int)j;
-(double)variableNonLinValues:(int)i;

-(void)setVariableName:(char *)string;
-(void)setVariableDofs:(int)i;
-(void)setVariablePerm:(int)i: (int)n;
-(void)setVariableNorm:(double)n;
-(void)setVariableNonLinChange:(double)n;
-(void)setVariableSteadyChange:(double)n;
-(void)setVariableNonLinConverged:(int)n;
-(void)setVariableSteadyConverged:(int)n;
-(void)setVariableNonLinIter:(int)n;
-(void)setVariableValues:(int)i: (double)n;
-(void)setVariablePrevValues:(int)i: (int)j: (double)n;
-(void)setVariableNonLinValues:(int)i: (double)n;


#pragma mark Sizes

-(int)sizeOfBoundaryReorder;
-(int)variableSizeOfPerm;
-(int)variableSizeOfValues;
-(int)variableSizeOfNonLinValues;

-(void)setSizeOfBoundaryReorder:(int)i;
-(void)setVariableSizeOfPerm:(int)i;
-(void)setVariableSizeOfValues:(int)i;
-(void)setVariableSizeOfNonLinValues: (int)i;

#pragma Test associativity

-(BOOL)isAssociatedVariablePrevValues;
-(BOOL)isAssociatedVariableNonLinValues;
-(BOOL)isAssociatedVariableSteadyValues;


#pragma mark Methods returning pointers

-(int *)returnPointerToDefDofs;

-(double *)variableReturnPointerToValues;
-(double *)variableReturnPointerToNonLinValues;
-(double *)variableReturnPointerToSteadyValues;

-(Variable_t *)meshReturnPointerToVariables;

#pragma mark Initializations

-(void)initializeILU:(int)ilun;
-(void)initializeCILU:(int)ilun;
-(void)ilutWorkspaceCheck:(int)i: (int)n;
-(void)ilutComplexWorkspaceCheck:(int)i: (int)n;

#pragma mark Allocations

-(void)allocateVariableNonLinValues:(int)n;
-(void)freeVariableNonLinValues:(int)n;



@end
