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
    int timeOrder, doneTime, order, nofEigenValues;
    int lineBeforeProc, linAfterProc;
    double alpha, beta, dt;
    int solutionComputeWhen;
    int numberOfActiveElements;
   
    FEMMatrix *_matrix;
    FEMVariable *_variable;
    FEMMesh *_mesh;
    NSMutableDictionary *_exportedVariables;         // Mutable dictionary holding FEMVariable classes
    
    solutionArraysContainer *_containers;
    ValueList_t values;
    
}

@property(nonatomic, strong) FEMMatrix *matrix;
@property(nonatomic, strong) FEMVariable *variable;
@property(nonatomic, strong) FEMMesh *mesh;
@property(nonatomic, strong) NSMutableDictionary *exportedVariables;

-(solutionArraysContainer *)getContainers;

#pragma mark Setters and getters for FEMSolution variables

-(id)solutionInfoForKey:(NSString *)key;
-(NSString *)normalTangentialName;
-(int)simulationID;
-(int)normalTangentialNOFNodes;
-(int)coordinateSystemDimension;
-(int)timeOrder;
-(int)doneTime;
-(int)order;
-(int)nofEigenValues;
-(double)alpha;
-(double)beta;
-(double)dt;

-(void)setNormalTangentialName:(NSString *)string;
-(void)setSimulationID:(int)i;
-(void)setNormalTangentialNOFNodes:(int)i;
-(void)setCoordinateSystemDimension:(int)i;
-(void)setTimeOrder:(int)n;
-(void)setDoneTime:(int)n;
-(void)setOrder:(int)n;
-(void)setAlpha:(double)n;
-(void)setBeta:(double)n;
-(void)setDt:(double)n;


#pragma mark Initializations

-(void)initializeILU:(int)ilun;
-(void)initializeCILU:(int)ilun;
-(void)ilutWorkspaceCheck:(int)i: (int)n;
-(void)ilutComplexWorkspaceCheck:(int)i: (int)n;



@end
