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
    BOOL _hasBuiltInSolution;
    id _builtInSolution;                             // If a built-in class provides the solution we will just instanciate it,
    id _plugInPrincipalClassInstance;                // otherwise the instance of the principal class if a plug-in is provided
    NSString *_plugInName;                           // Name of the plug-in if present
    FEMMatrix *_matrix;
    FEMVariable *_variable;
    FEMMesh *_mesh;
    NSMutableDictionary *_solutionInfo;
    NSString *_normalTangentialName;
    NSMutableDictionary *_exportedVariables;         // Exported variables
    NSMutableArray *_valuesList;
    
    solutionArraysContainer * __nonnull _containers;
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
@property(nonatomic, assign) BOOL hasBuiltInSolution;
@property(nonatomic, strong, nullable) id builtInSolution;
@property(nonatomic, strong, nullable) id plugInPrincipalClassInstance;
@property(nonatomic, strong, nullable) NSString *plugInName;
@property(nonatomic, strong, nullable) FEMMatrix *matrix;
@property(nonatomic, strong, nullable) FEMVariable *variable;
@property(nonatomic, strong, nullable) FEMMesh *mesh;
@property(nonatomic, strong, nonnull) NSMutableDictionary *solutionInfo;
@property(nonatomic, strong, nullable) NSString *normalTangentialName;
@property(nonatomic, strong, nullable) NSMutableDictionary <NSString *, FEMVariable *> *exportedVariables;
@property(nonatomic, strong, nonnull) NSMutableArray <FEMValueList *> *valuesList;

-(void)deallocation;
-(solutionArraysContainer * __nonnull)getContainers;

// Instantiate the plug-in principal class
-(BOOL)instantiatePrincipalClassFromPlugIn:(NSBundle * __nonnull)bundle;

@end
