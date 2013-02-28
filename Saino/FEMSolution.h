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
#import "FEMValueList.h"

#import "Constructors.h"
#import "memory.h"
#import "Utils.h"

@interface FEMSolution : NSObject {
    
    int _simulationID;
    int _coordinateSystemDimension;
    int _timeOrder;
    int _doneTime;
    int _order;
    int _nOfEigenValues;
    int _solutionSolveWhen;
    int _solutionMode;
    int _numberOfActiveElements;
    int _multigridLevel;
    int _multiGridTotal;
    int _multiGridSweep;
    double _alpha;
    double _beta;
    double _dt;
    BOOL _multigridSolution;
    BOOL _multigridEqualPlit;
    id _plugInPrincipalClassInstance;                // The instance of the principal class if a plug-in is provided
    FEMMatrix *_matrix;
    FEMVariable *_variable;
    FEMMesh *_mesh;
    NSMutableDictionary *_solutionInfo;
    NSString *_normalTangentialName;
    NSMutableDictionary *_exportedVariables;         // Mutable dictionary holding FEMVariable classes for exported variables
    NSMutableArray *_valuesList;                     // Array of FEMValueList objects
    
    solutionArraysContainer *_containers;
}

@property(nonatomic, assign) int simulationID;
@property(nonatomic, assign) int coordinateSystemDimension;
@property(nonatomic, assign) int timeOrder;
@property(nonatomic, assign) int doneTime;
@property(nonatomic, assign) int order;
@property(nonatomic, assign) int nOfEigenValues;
@property(nonatomic, assign) int solutionSolveWhen;
@property(nonatomic, assign) int solutionMode;
@property(nonatomic, assign) int numberOfActiveElements;
@property(nonatomic, assign) int multigridLevel;
@property(nonatomic, assign) int multiGridTotal;
@property(nonatomic, assign) int multiGridSweep;
@property(nonatomic, assign) double alpha;
@property(nonatomic, assign) double beta;
@property(nonatomic, assign) double dt;
@property(nonatomic, assign, getter = isMultigridSolution) BOOL multigridSolution;
@property(nonatomic, assign, getter = isMultigridEqualPlit) BOOL multigridEqualPlit;
@property(nonatomic, strong) id plugInPrincipalClassInstance;
@property(nonatomic, strong) FEMMatrix *matrix;
@property(nonatomic, strong) FEMVariable *variable;
@property(nonatomic, strong) FEMMesh *mesh;
@property(nonatomic, strong) NSMutableDictionary *solutionInfo;
@property(nonatomic, strong) NSString *normalTangentialName;
@property(nonatomic, strong) NSMutableDictionary *exportedVariables;
@property(nonatomic, strong) NSMutableArray *valuesList;


-(void)deallocation;
-(solutionArraysContainer *)getContainers;

#pragma mark Initializations

-(void)initializeILU:(int)ilun;
-(void)initializeCILU:(int)ilun;
-(void)ilutWorkspaceCheckAtIndex:(int)i numberOfRows:(int)n;
-(void)ilutComplexWorkspaceCheckAtIndex:(int)i numberOfRows:(int)n;



@end
