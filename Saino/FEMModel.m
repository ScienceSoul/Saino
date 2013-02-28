//
//  FEMModel.m
//  Saino
//
//  Created by Hakime Seddik on 18/01/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMModel.h"

#import "FEMKernel.h"
#import "FEMListUtilities.h"
#import "FEMElementUtils.h"
#import "FEMUtilities.h"
#import "FEMMeshUtils.h"

@interface FEMModel ()
-(void)FEMModel_setCoordinateSystem;
-(void)FEMModel_localMatrix:(double **)stiff force:(double*)force mesh:(FEMMesh *)mesh element:(Element_t *)element numberOfNodes:(int)numberOfNodes power:(double)power noWeight:(BOOL)noWeight elemMin:(double *)elemMin elemMax:(double *)elemMax;
-(void)FEMModel_getNodalElementSize:(double)expo weight:(BOOL)weight nodal:(double *)h sizeNodal:(int)sizeNodal;
@end

@implementation FEMModel

@synthesize dimension = _dimension;
@synthesize numberOfNodes = _numberOfNodes;
@synthesize numberOfBulkElements = _numberOfBulkElements;
@synthesize numberOfBoundaryElements = _numberOfBoundaryElements;
@synthesize numberOfBodies = _numberOfBodies;
@synthesize numberOfBodyForces = _numberOfBodyForces;
@synthesize numberOfBoundaryConditions = _numberOfBoundaryConditions;
@synthesize numberOfBoundaries = _numberOfBoundaries;
@synthesize numberOfInitialConditions = _numberOfInitialConditions;
@synthesize numberOfSolutions = _numberOfSolutions;
@synthesize numberOfEquations = _numberOfEquations;
@synthesize numberOfMaterials = _numberOfMaterials;
@synthesize coordinates = _coordinates;
@synthesize totalMatrixElements = _totalMatrixElements;
@synthesize maxElementNodes = _maxElementNodes;
@synthesize mesh = _mesh;
@synthesize solution = _solution;
@synthesize outputPath = _outputPath;
@synthesize boundaryID = _boundaryID;
@synthesize solutions = _solutions;
@synthesize meshes = _meshes;
@synthesize bodyForces = _bodyForces;
@synthesize boundaryConditions = _boundaryConditions;
@synthesize boundaries = _boundaries;
@synthesize equations = _equations;
@synthesize simulation = _simulation;
@synthesize variables = _variables;

#pragma mark Private methods

-(void)FEMModel_setCoordinateSystem {
    
    int i, j;
    double x, y, z;
    BOOL found, c[3], any;
    NSString *csys;
    Nodes_t *nodes;
    FEMMesh *mesh;
    FEMListUtilities *listUtil;
    
    listUtil = [[FEMListUtilities alloc] init];
    
    csys = [listUtil listGetString:self inArray:self.simulation.valuesList forVariable:@"coordinate system" info:&found];
    if (found == NO) csys = @"cartesian";
    
    if ([csys isEqualToString:@"cartesian"] || [csys isEqualToString:@"polar"]) {
        mesh = self.meshes[0];
        nodes = mesh.getNodes;
        x = nodes->x[0];
        y = nodes->y[0];
        z = nodes->z[0];
        for (i=0; i<3; i++) {
            c[i] = NO;
        }
        for (FEMMesh *mesh in self.meshes) {
            nodes = mesh.getNodes;
            any = NO;
            for (i=0; i<mesh.numberOfNodes; i++) {
                if (nodes->x[i] != x) {
                    any = YES;
                    break;
                }
            }
            c[0] = (c[0] == YES || any == YES) ? YES : NO;
            
            any = NO;
            for (i=0; i<mesh.numberOfNodes; i++) {
                if (nodes->y[i] != y) {
                    any = YES;
                    break;
                }
            }
            c[1] = (c[1] == YES || any == YES) ? YES : NO;
            
            any = NO;
            for (i=0; i<mesh.numberOfNodes; i++) {
                if (nodes->z[i] != z) {
                    any = YES;
                    break;
                }
            }
            c[2] = (c[2] == YES || any == YES) ? YES : NO;
        }
        j = 0;
        for (i=0; i<3; i++) {
            if (c[i] == YES) j++;
        }
        self.dimension = j;
    }
    
    if ([csys isEqualToString:@"cartesian"]) {
        self.coordinates = cartesian;
    } else if ([csys isEqualToString:@"cartesian 1d"]) {
        self.dimension = 1;
        self.coordinates = cartesian;
    } else if ([csys isEqualToString:@"cartesian 2d"]) {
        self.dimension = 2;
        self.coordinates = cartesian;
    } else if ([csys isEqualToString:@"cartesian 3d"]) {
        self.dimension = 3;
        self.coordinates = cartesian;
    } else if ([csys isEqualToString:@"axis symmetric"]) {
        self.dimension = 2;
        self.coordinates = axis_symmetric;
    } else if ([csys isEqualToString:@"cylindric symmetric"]) {
        self.dimension = 2;
        self.coordinates = cylindric_symmetric;
    } else if ([csys isEqualToString:@"cylindrical"]) {
        self.dimension = 3;
        self.coordinates = cylindric;
    } else if ([csys isEqualToString:@"polar"]) {
        self.coordinates = polar;
    } else if ([csys isEqualToString:@"polar 2d"]) {
        self.dimension = 2;
        self.coordinates = polar;
    } else if ([csys isEqualToString:@"polar 3d"]) {
        self.dimension = 3;
        self.coordinates = polar;
    } else {
        warnfunct("FEMModel_setCoordinateSystem", "Unknown global coordinate system:");
        printf("%s\n", [csys UTF8String]);
        errorfunct("FEMModel_setCoordinateSystem", "Program terminating now...");
    }
}

-(void)FEMModel_localMatrix:(double **)stiff force:(double*)force mesh:(FEMMesh *)mesh element:(Element_t *)element numberOfNodes:(int)numberOfNodes power:(double)power noWeight:(BOOL)noWeight elemMin:(double *)elemMin elemMax:(double *)elemMax {
    
    int i, j, t;
    double detJ, loadAtIP, weight;
    GaussIntegrationPoints *IP;
    Nodes_t nodes, *meshNodes;
    FEMNumericIntegration *integration;
    
    integration = [[FEMNumericIntegration alloc] init];
    if ([integration allocation:mesh] == NO) errorfunct("FEMModel_localMatrix", "Allocation error in FEMNumericIntegration!");
    
    nodes.x = doublevec(0, numberOfNodes-1);
    nodes.y = doublevec(0, numberOfNodes-1);
    nodes.z = doublevec(0, numberOfNodes-1);
    
    meshNodes = mesh.getNodes;
    
    for (i=0; i<numberOfNodes; i++) {
        nodes.x[i] = meshNodes->x[element->NodeIndexes[i]];
        nodes.y[i] = meshNodes->y[element->NodeIndexes[i]];
        nodes.z[i] = meshNodes->z[element->NodeIndexes[i]];
    }
    
    // Numerical integration
    IP = GaussQuadrature(element, NULL, NULL);
    
    for (t=0; t<IP->n; t++) {
        // Basis function values & derivatives at the integration point
        [integration setBasisForElement:element elementNodes:&nodes inMesh:mesh firstEvaluationPoint:IP->u[t] secondEvaluationPoint:IP->v[t] thirdEvaluationPoint:IP->w[t] withBubbles:NO basisDegree:NULL];
        [integration setMetricDeterminantForElement:element elementNodes:&nodes inMesh:mesh firstEvaluationPoint:IP->u[t] secondEvaluationPoint:IP->v[t] thirdEvaluationPoint:IP->w[t]];
        
        detJ = integration.metricDeterminant;
        
        // The source term at the integration point
        loadAtIP = pow(detJ, power);
        if (noWeight == YES) {
            weight = IP->s[t];
        } else {
            weight = IP->s[t] * detJ;
        }
        
        *elemMin = min(*elemMin, loadAtIP);
        *elemMax = max(*elemMax, loadAtIP);
        
        // Finally the elemental matrix & vector
        for (i=0; i<numberOfNodes; i++) {
            for (j=0; j<numberOfNodes; j++) {
                stiff[i][j] = stiff[i][j] + weight * integration.basis[i] * integration.basis[j];
            }
            force[i] = force[i] + weight * integration.basis[i] * loadAtIP;
        }
    }
    
    free_dvector(nodes.x, 0, numberOfNodes-1);
    free_dvector(nodes.y, 0, numberOfNodes-1);
    free_dvector(nodes.z, 0, numberOfNodes-1);
    
    free_dvector(IP->u, 0, MAX_INTEGRATION_POINTS-1);
    free_dvector(IP->v, 0, MAX_INTEGRATION_POINTS-1);
    free_dvector(IP->w, 0, MAX_INTEGRATION_POINTS-1);
    free_dvector(IP->s, 0, MAX_INTEGRATION_POINTS-1);
    free(IP);
    
    [integration deallocation:mesh];
}

-(void)FEMModel_getNodalElementSize:(double)expo weight:(BOOL)weight nodal:(double *)h sizeNodal:(int)sizeNodal {
    
    int i, j, n, t, active;
    int *cperm;
    double power, elemMin, elemMax;
    double **stiff, *force;
    BOOL onlySearch, found;
    Element_t *elements;
    FEMSolution *solution;
    FEMMatrix *matrix;
    FEMListUtilities *listUtil;
    FEMElementUtils *elementUtils;
    FEMUtilities *utils;
    FEMMatrixCRS *crsMatrix;
    FEMKernel *kernel;
    matrixArraysContainer *matrixContainers = NULL;
    variableArraysContainer *bufferContainers = NULL;
    
    solution = [[FEMSolution alloc] init];
    solution.mesh = self.mesh;
    
    [listUtil addStringInClassList:solution theVariable:@"linear system iterative method" withValue:@"cg"];
    [listUtil addLogicalInClassList:solution theVariable:@"linear system symmetric" withValue:YES];
    [listUtil addIntegerInClassList:solution theVariable:@"linear system maximum iterations" withValue:5000];
    [listUtil addStringInClassList:solution theVariable:@"linear system preconditioning" withValue:@"ilu0"];
    [listUtil addIntegerInClassList:solution theVariable:@"linear system residual output" withValue:100];
    [listUtil addConstRealInClassList:solution theVariable:@"linear system convergence tolerance" withValue:1.0E-9 string:nil];
    
    cperm = intvec(0, solution.mesh.numberOfNodes-1);
    
    elementUtils = [[FEMElementUtils alloc] init];
    matrix = [elementUtils createMatrixInModel:self forSolution:solution mesh:solution.mesh dofs:1 permutation:cperm sizeOfPermutation:solution.mesh.numberOfNodes matrixFormat:MATRIX_CRS optimizeBandwidth:NO equationName:nil discontinuousGalerkinSolution:NULL globalBubbles:NULL];
    
    // TODO: Assign here matrix field (Comm, ParMatrix) for parallel run
    
    solution.matrix = matrix;
    self.solution = solution;
    
    matrixContainers = matrix.getContainers;
    matrixContainers->RHS = doublevec(0, solution.mesh.numberOfNodes-1);
    matrixContainers->sizeRHS = solution.mesh.numberOfNodes;
    
    solution.timeOrder = 0;
    
    utils = [[FEMUtilities alloc] init];
    bufferContainers = (variableArraysContainer*)malloc(sizeof(variableArraysContainer));
    bufferContainers->Values = h;
    bufferContainers->Perm = cperm;
    bufferContainers->sizeValues = sizeNodal;
    bufferContainers->sizePerm = solution.mesh.numberOfNodes;
    [utils addVariableTo:solution.mesh.variables mesh:solution.mesh solution:solution name:@"nodal h" dofs:1 container:bufferContainers component:NO ifOutput:NULL ifSecondary:NULL type:NULL];
    onlySearch = YES;
    solution.variable = [utils getVariableFrom:solution.mesh.variables model:self name:@"nodal h" onlySearch:&onlySearch maskName:NULL info:&found];
    
    // TODO: Add parallel run support here
    
    //Allocate some storage
    force = doublevec(0, solution.mesh.maxElementNodes-1);  // Just big enough for elemental arrays
    stiff = doublematrix(0, solution.mesh.maxElementNodes-1, 0, solution.mesh.maxElementNodes-1);
    for (i=0; i<solution.mesh.maxElementNodes; i++) {
        for (j=0; j<solution.mesh.maxElementNodes; j++) {
            stiff[i][j] = 0.0;
        }
        force[i] = 0.0;
    }
    
    elemMin = HUGE_VALL;
    elemMin = -HUGE_VALL;
    
    // Initialize the system and do the assembly
    power = 1.0 / expo;
    active = solution.mesh.numberOfBulkElements;
    
    memset( matrixContainers->RHS, 0.0, (matrixContainers->sizeRHS*sizeof(matrixContainers->RHS)) );
    memset( matrixContainers->Values, 0.0, (matrixContainers->sizeValues*sizeof(matrixContainers->Values)) );
    
    crsMatrix = [[FEMMatrixCRS alloc] init];
    
    elements = solution.mesh.getElements;
    for (t=0; t<active; t++) {
        n = elements[t].Type.NumberOfNodes;
        
        // Get element local matrix and rhs vector
        [self FEMModel_localMatrix:stiff force:force mesh:solution.mesh element:&elements[t] numberOfNodes:n power:power noWeight:weight elemMin:&elemMin elemMax:&elemMax];
        
        // Update global matrix and rhs  vector from local matrix and vector
        [crsMatrix glueLocalMatrixInMatrix:matrix localMatrix:stiff numberOfNodes:n dofs:1 indexes:elements[t].NodeIndexes];
        for (i=0; i<elements[t].Type.NumberOfNodes; i++) {
            matrixContainers->RHS[elements[t].NodeIndexes[i]] = matrixContainers->RHS[elements[t].NodeIndexes[i]] + force[i];
        }
    }
    
    memset( h, 0.0, (sizeNodal*sizeof(h)) );
    
    kernel = [[FEMKernel alloc] init];
    // TODO: Add support for parallel run here
    // ....
    [kernel iterativeSolve:solution];
    
    NSLog(@"FEMModel_getNodalElementSize: minimum element size: %f %f\n", elemMin, min_array(h,sizeNodal));
    NSLog(@"FEMModel_getNodalElementSize: maximum element size: %f %f\n", elemMax, max_array(h,sizeNodal));
    NSLog(@"FEMModel_getNodalElementSize: element size ratio: %f %f\n", elemMax/elemMin, max_array(h,sizeNodal) / min_array(h,sizeNodal));
    
    self.solution = nil;
    solution.mesh.variables = nil;
    
    [matrix deallocation];
    matrix = nil;
    
    free_dvector(force, 0, solution.mesh.maxElementNodes-1);
    free_dmatrix(stiff, 0, solution.mesh.maxElementNodes-1, 0, solution.mesh.maxElementNodes-1);
    free_ivector(cperm, 0, solution.mesh.numberOfNodes-1);
    
    bufferContainers->Values = NULL;
    bufferContainers->Perm = NULL;
    free(bufferContainers);
    
    [kernel deallocation];
}


#pragma mark Public methods

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
        _dimension = 0;
        _numberOfNodes = 0;
        _numberOfBulkElements = 0;
        _numberOfBoundaryElements = 0;
        _numberOfBodies = 0;
        _numberOfBodyForces = 0;
        _numberOfBoundaryConditions = 0;
        _numberOfBoundaries = 0;
        _numberOfInitialConditions = 0;
        _numberOfSolutions = 0;
        _numberOfMaterials = 0;
        _coordinates = cartesian;
        _totalMatrixElements = 0;
        _maxElementNodes = 0;
        _mesh = nil;
        
        _outputPath = [NSMutableString stringWithString:@" "];
        _boundaryID = [[NSArray alloc] init];
        _solutions = [[NSArray alloc] init];
        _meshes = [[NSMutableArray alloc] init];
        _bodies = [[NSArray alloc] init];
        _bodyForces = [[NSArray alloc] init];
        _boundaryConditions = [[NSArray alloc] init];
        _boundaries = [[NSArray alloc] init];
        _equations = [[NSArray alloc] init];
        _variables = [[NSMutableArray alloc] init];
        _simulation = [[FEMSimulation alloc] init];
        
        _elements = NULL;
        _currentElement = NULL;
            
        _containers = (modelArraysContainer*)malloc(sizeof(modelArraysContainer));
        _containers->freeSurfaceNodes = NULL;
        _containers->rowNonZeros = NULL;
        _containers->boundaryCurvatures = NULL;
    }
    
    return self;
}

-(void)deallocation {
    free(_containers);
}

-(void)loadModelName:(NSString *)name boundariesOnly:(BOOL)bd dummy:(int *)d1 dummy:(int *)d2 {
    
    int i, j, l, nlen, meshKeep, meshLevels, sizeNodal, numberPartitions, partitionID;
    int defDofs[6];
    double meshPower, *h, maxvalue;
    BOOL transient, gotMesh, oneMeshName, meshGrading, found, single;
    unsigned long loc;
    NSString *elementDef, *meshString, *aString;
    NSMutableString *meshDir, *meshName;
    FEMMesh *mesh, *newMesh, *oldMesh, *modelMesh;
    FEMListUtilities *listUtils;
    FEMMeshUtils *meshUtils;
    solutionArraysContainer *solContainers = NULL;
    NSRange substr1;
    
    listUtils = [[FEMListUtilities alloc] init];
    
    //TODO: Here comes the Model Description File (MDF) parser
    
    transient = ([[listUtils listGetString:self inArray:self.simulation.valuesList forVariable:@"simulation type" info:&found] isEqualToString:@"transient"]) ? YES : NO;
    
    memset( defDofs, -1, sizeof(defDofs) );
    defDofs[0] = 1;
    
    for (FEMSolution *solution in self.solutions) {
        
        // TODO: Here Elmer calls a procedure post-fixed with "_Init0" for a given solver.
        // Do we need to do that?
        
        gotMesh = [listUtils listCheckPresentVariable:@"mesh" inArray:solution.valuesList];
        
        solContainers = solution.getContainers;
        if (solContainers->defDofs == NULL) {
            solContainers->defDofs = intmatrix(0, self.numberOfBodies-1, 0, 5);
            solContainers->size1DefDofs = self.numberOfBodies;
            solContainers->size2DefDofs = 6;
            for (i=0; i<self.numberOfBodies; i++) {
                for (j=0; j<6; j++) {
                    solContainers->defDofs[i][j] = -1;
                }
            }
            for (i=0; i<self.numberOfBodies; i++) {
                solContainers->defDofs[i][0] = 1;
            }
        }
        
        // Define what kind of element we are working with this solver
        if ((solution.solutionInfo)[@"element"] == nil) {
            if ([(solution.solutionInfo)[@"discontinuous galerkin"] boolValue] == YES) {
                for (i=0; i<self.numberOfBodies; i++) {
                    solContainers->defDofs[i][3] = 0;
                }
                if (gotMesh == NO) defDofs[3] = max(defDofs[3], 0);
                continue;
            } else {
                elementDef = @"n:1";
            }
        } else {
            elementDef = (solution.solutionInfo)[@"element"];
        }
        
        substr1 = [elementDef rangeOfString:@"n:"];
        if (substr1.location != NSNotFound) {
            l = [[NSString stringWithFormat:@"%c", [elementDef characterAtIndex:(substr1.location+2)]] intValue];
            for (i=0; i<self.numberOfBodies; i++) {
                solContainers->defDofs[i][0] = l;
            }
            if (gotMesh == NO) defDofs[0] = max(defDofs[0], l);
        }
        
        substr1 = [elementDef rangeOfString:@"e:"];
        if (substr1.location != NSNotFound) {
            l = [[NSString stringWithFormat:@"%c", [elementDef characterAtIndex:(substr1.location+2)]] intValue];
            for (i=0; i<self.numberOfBodies; i++) {
                solContainers->defDofs[i][1] = l;
            }
            if (gotMesh == NO) defDofs[1] = max(defDofs[1], l);
        }
        
        substr1 = [elementDef rangeOfString:@"f:"];
        if (substr1.location != NSNotFound) {
            l = [[NSString stringWithFormat:@"%c", [elementDef characterAtIndex:(substr1.location+2)]] intValue];
            for (i=0; i<self.numberOfBodies; i++) {
                solContainers->defDofs[i][2] = l;
            }
            if (gotMesh == NO) defDofs[2] = max(defDofs[2], l);
        }
        
        substr1 = [elementDef rangeOfString:@"d:"];
        if (substr1.location != NSNotFound) {
            l = [[NSString stringWithFormat:@"%c", [elementDef characterAtIndex:(substr1.location+2)]] intValue];
            for (i=0; i<self.numberOfBodies; i++) {
                solContainers->defDofs[i][3] = l;
            }
            if (gotMesh == NO) defDofs[2] = max(defDofs[3], l);
        } else {
            if ([(solution.solutionInfo)[@"discontinuous galerkin"] boolValue] == YES) {
                for (i=0; i<self.numberOfBodies; i++) {
                    solContainers->defDofs[i][3] = 0;
                }
                if (gotMesh == NO)  defDofs[3] = max(defDofs[3], 0);
            }
        }
        
        substr1 = [elementDef rangeOfString:@"b:"];
        if (substr1.location != NSNotFound) {
            l = [[NSString stringWithFormat:@"%c", [elementDef characterAtIndex:(substr1.location+2)]] intValue];
            for (i=0; i<self.numberOfBodies; i++) {
                solContainers->defDofs[i][4] = l;
            }
            if (gotMesh == NO) defDofs[4] = max(defDofs[4], l);
        }
        
        substr1 = [elementDef rangeOfString:@"p:"];
        if (substr1.location != NSNotFound) {
            if ([elementDef characterAtIndex:(substr1.location+2)] == '%') {
                for (i=0; i<self.numberOfBodies; i++) {
                    solContainers->defDofs[i][5] = 0;
                }
            } else {
                l = [[NSString stringWithFormat:@"%c", [elementDef characterAtIndex:(substr1.location+2)]] intValue];
                for (i=0; i<self.numberOfBodies; i++) {
                    solContainers->defDofs[i][5] = l;
                }
                if (gotMesh == NO) defDofs[5] = max(defDofs[5], l);
            }
        }
    }
    
    // Check the mesh
    meshString = [listUtils listGetString:self inArray:self.simulation.valuesList forVariable:@"mesh" info:&found];
    
    oneMeshName = NO;
    if (found == YES) {
        substr1 = [meshString rangeOfString:@" "];
        if (substr1.location == NSNotFound) {
            substr1.location = [meshString length];
            meshName = [NSMutableString stringWithString:[meshString substringToIndex:substr1.location]];
        } else {
            meshDir = [NSMutableString stringWithString:[meshString substringToIndex:substr1.location]];
            meshName = [NSMutableString stringWithString:[meshString substringToIndex:substr1.location]];
        }
        
        loc = substr1.location;
        while (loc < [meshString length] && [meshString characterAtIndex:loc] == ' ') {
            loc++;
        }
        
        if (loc < [meshString length]) {
            [meshName appendString:@"/"];
            [meshName appendString:[meshString substringFromIndex:loc]];
        } else {
            oneMeshName = YES;
            meshDir = [NSMutableString stringWithString:@"."];
        }
    }
    
    meshUtils = [[FEMMeshUtils alloc] init];
    
    if ([meshDir characterAtIndex:0] != ' ') {
        mesh = [[FEMMesh alloc] init];
        [mesh loadMeshForModel:self meshDirectory:meshDir meshName:meshName boundariesOnly:bd numberOfPartitions:NULL partitionID:NULL definitions:defDofs];
        [self.meshes addObject:mesh];
        
        [self FEMModel_setCoordinateSystem];
        
        meshLevels = [listUtils listGetInteger:self inArray:self.simulation.valuesList forVariable:@"mesh levels" info:&found minValue:NULL maxValue:NULL];
        if (found == NO) meshLevels = 1;
        
        meshKeep = [listUtils listGetInteger:self inArray:self.simulation.valuesList forVariable:@"mesh keep" info:&found minValue:NULL maxValue:NULL];
        if (found == NO) meshKeep = meshLevels;
        
        meshPower = [listUtils listGetConstReal:self inArray:self.simulation.valuesList forVariable:@"mesh grading power" info:&found minValue:NULL maxValue:NULL];
        meshGrading = [listUtils listGetLogical:self inArray:self.simulation.valuesList forVariable:@"mesh keep grading" info:&found];
        
        for (i=2; i<=meshLevels; i++) {
            oldMesh = self.meshes[0];
            
            if (meshGrading == YES) {
                h = doublevec(0, oldMesh.numberOfNodes-1);
                self.mesh = oldMesh;
                sizeNodal = oldMesh.numberOfNodes;
                [self FEMModel_getNodalElementSize:meshPower weight:NO nodal:h sizeNodal:sizeNodal];
                 newMesh = [meshUtils splitMeshEqual:oldMesh model:self nodal:h sizeNodal:&sizeNodal];
                free_dvector(h, 0, oldMesh.numberOfNodes-1);
            } else {
                newMesh = [meshUtils splitMeshEqual:oldMesh model:self nodal:NULL sizeNodal:NULL];
            }
            
            if (i>meshLevels-meshKeep+1) {
                [newMesh.next addObject:oldMesh];
                [newMesh.parent addObject:oldMesh];
                [oldMesh.child addObject:newMesh];
                [newMesh.name setString:oldMesh.name];
                newMesh.outputActive = YES;
                oldMesh.outputActive = NO;
            } else {
                [oldMesh deallocation];
                oldMesh = nil;
            }
            self.meshes[0] = newMesh;
        }
        
        if (oneMeshName == YES) {
            i = 0;
        } else {
            i = (int)[meshName length]-1;
            while (i >= 0 && [meshName characterAtIndex:i] != '/') {
                i--;
            }
            i++;
        }
        
        modelMesh = self.meshes[0];
        modelMesh.name = [NSString stringWithString:[meshName substringFromIndex:i]];
        
        for (FEMSolution *solution in self.solutions) {
            solution.mesh = self.meshes[0];
        }
    }
    
    for (FEMSolution *solution in self.solution) {
        
        if ((solution.solutionInfo)[@"mesh"] != nil) {
            aString = (solution.solutionInfo)[@"mesh"];
            single = NO;
            if ([[aString substringToIndex:8] isEqualToString:@"-single "] == YES) {
                single = YES;
                aString = [NSString stringWithString:[aString substringFromIndex:8]];
            }
            oneMeshName = NO;
            nlen = (int)[aString length];
            substr1 = [aString rangeOfString:@" "];
            if (substr1.location == NSNotFound) {
                substr1.location = [aString length];
                meshName = [NSMutableString stringWithString:[aString substringToIndex:substr1.location]];
            } else {
                meshDir = [NSMutableString stringWithString:[aString substringToIndex:substr1.location]];
                meshName = [NSMutableString stringWithString:[aString substringToIndex:substr1.location]];
            }
            
            loc = substr1.location;
            while (loc < nlen && [aString characterAtIndex:loc] == ' ') {
                loc++;
            }
            if (loc < nlen) {
                [meshName appendString:@"/"];
                [meshName appendString:[aString substringFromIndex:loc]];
            } else {
                oneMeshName = YES;
                 meshDir = [NSMutableString stringWithString:@"."];
            }
            
            if (oneMeshName == YES) {
                i = 0;
            } else {
                i = (int)[meshName length]-1;
                while (i >= 0 && [meshName characterAtIndex:i] != '/') {
                    i--;
                }
                i++;
            }
            
            found = NO;
            for (FEMMesh *mesh in self.meshes) {
                found = YES;
                if ([mesh.name isEqualToString:[meshName substringFromIndex:i]] != YES) found = NO;
                if ([mesh.name length] != [[meshName substringFromIndex:i] length]) found = NO;
                if (found == YES) {
                    solution.mesh = mesh;
                    break;
                }
            }
            
            if (found) continue;
            
            solContainers = solution.getContainers;
            
            for (i=0; i<6; i++) {
                maxvalue = -HUGE_VAL;
                for (j=0; j<solContainers->size1DefDofs; j++) {
                    if (solContainers->defDofs[j][i]>maxvalue) {
                        maxvalue = solContainers->defDofs[j][i];
                    }
                }
                defDofs[i] = maxvalue;
            }
            
            if (single == YES) {
                numberPartitions = 1;
                partitionID = 0;
                solution.mesh = [[FEMMesh alloc] init];
                [solution.mesh loadMeshForModel:self meshDirectory:meshDir meshName:meshName boundariesOnly:bd numberOfPartitions:&numberPartitions partitionID:&partitionID definitions:defDofs];
            } else {
                // TODO: add support for parallel run
                // We are not supporting parallel runs yet, so we should never reach here
                solution.mesh = [[FEMMesh alloc] init];
                [solution.mesh loadMeshForModel:self meshDirectory:meshDir meshName:meshName boundariesOnly:bd numberOfPartitions:NULL partitionID:NULL definitions:defDofs];
            }
            
            if ((solution.solutionInfo)[@"mesh levels"] != nil) {
                meshLevels = [(solution.solutionInfo)[@"mesh levels"] intValue];
            } else meshLevels = 1;
            
            if ((solution.solutionInfo)[@"mesh keep"] != nil) {
                meshKeep = [(solution.solutionInfo)[@"mesh keep"] intValue];
            } else meshKeep = meshLevels;
            
            meshPower = [listUtils listGetConstReal:self inArray:self.simulation.valuesList forVariable:@"mesh grading power" info:&found minValue:NULL maxValue:NULL];
            meshGrading = [listUtils listGetLogical:self inArray:self.simulation.valuesList forVariable:@"mesh keep grading" info:&found];
            
            for (i=2; i<=meshLevels; i++) {
                oldMesh = solution.mesh;
                
                if (meshGrading == YES) {
                    h = doublevec(0, oldMesh.numberOfNodes-1);
                    self.mesh = oldMesh;
                    sizeNodal = oldMesh.numberOfNodes;
                    [self FEMModel_getNodalElementSize:meshPower weight:NO nodal:h sizeNodal:sizeNodal];
                    newMesh = [meshUtils splitMeshEqual:oldMesh model:self nodal:h sizeNodal:&sizeNodal];
                    free_dvector(h, 0, oldMesh.numberOfNodes-1);
                } else {
                    newMesh = [meshUtils splitMeshEqual:oldMesh model:self nodal:NULL sizeNodal:NULL];
                }
                
                if (i>meshLevels-meshKeep+1) {
                    [newMesh.next addObject:oldMesh];
                    [newMesh.parent addObject:oldMesh];
                    [oldMesh.child addObject:newMesh];
                    [newMesh.name setString:oldMesh.name];
                    newMesh.outputActive = YES;
                    oldMesh.outputActive = NO;
                } else {
                    [oldMesh deallocation];
                    oldMesh = nil;
                }
                solution.mesh = newMesh;
            }
            
            if (oneMeshName == YES) {
                i = 0;
            } else {
                i = (int)[meshName length]-1;
                while (i >= 0 && [meshName characterAtIndex:i] != '/') {
                    i--;
                }
                i++;
            }
            
            solution.mesh.name = [NSMutableString stringWithString:[meshName substringFromIndex:i]];
            
            // Just add this solver mesh to the global model meshes table
            [self.meshes addObject:solution.mesh];
        }
    }
    
    [self FEMModel_setCoordinateSystem];
    
    [self.outputPath setString:meshDir];

    for (FEMMesh *mesh in self.meshes) {
        [meshUtils SetStabilizationParametersInMesh:mesh model:self];
    }
}

#pragma mark Elements getter

-(Element_t *)getElements {
    
    return _elements;
}

-(Element_t *)getCurrentElement {
    
    return _currentElement;
}

-(modelArraysContainer*)getContainers {
    
    return _containers;
}


@end
