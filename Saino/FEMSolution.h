//
//  FEMSolution.h
//  Saino
//
//  Created by Hakime Seddik on 18/01/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMMesh.h"
#import "FEMMatrix.h"
#import "FEMVariable.h"

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
    
    FEMMatrix *_matrix;
    FEMVariable *_variable;
    FEMMesh *_mesh;
    NSMutableDictionary *_exportedVariable;         // Mutable dictionary holding FEMVariable classes
    
    ValueList_t values;
    
}

@property(nonatomic, strong) FEMMatrix *matrix;
@property(nonatomic, strong) FEMVariable *variable;
@property(nonatomic, strong) FEMMesh *mesh;
@property(nonatomic, strong) NSMutableDictionary *exportedVariables;

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

#pragma mark Sizes

-(int)sizeOfBoundaryReorder;
-(void)setSizeOfBoundaryReorder:(int)i;


#pragma mark Methods returning pointers

-(int *)returnPointerToDefDofs;

#pragma mark Initializations

-(void)initializeILU:(int)ilun;
-(void)initializeCILU:(int)ilun;
-(void)ilutWorkspaceCheck:(int)i: (int)n;
-(void)ilutComplexWorkspaceCheck:(int)i: (int)n;



@end
