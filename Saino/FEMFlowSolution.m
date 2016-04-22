//
//  FEMFlowSolution.m
//  Saino
//
//  Created by Seddik hakime on 15/03/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMFlowSolution.h"
#import "FEMCore.h"
#import "FEMUtilities.h"
#import "FEMListUtilities.h"
#import "FEMMesh.h"
#import "FEMTimeIntegration.h"
#import "FEMNavierStokes.h"
#import "FEMNavierStokesCylindrical.h"
#import "FEMNavierStokesGeneral.h"
#import "FEMFreeSurface.h"
#import "Utils.h"
#import "TimeProfile.h"
#import "GaussIntegration.h"
#import "OpenCLUtils.h"

@interface FEMFlowSolution ()
-(void)FEMFlowSolution_nullify;
-(void)FEMFlowSolution_checkCircleBoundaryModel:(FEMModel * __nonnull)model;
-(void)FEMFlowSolution_doAssembly:(FEMSolution * __nonnull)solution
                             core:(FEMCore * __nonnull)core
                            model:(FEMModel * __nonnull)model
                    listUtilities:(FEMListUtilities * __nonnull)listUtilities
               coordinatesSystems:(FEMCoordinateSystems * __nonnull)coordinatesSystems
                     navierStokes:(FEMNavierStokes * __nonnull)navierStokes
          navierStokesCylindrical:(FEMNavierStokesCylindrical * __nonnull)navierStokesCylindrical
              navierStokesGeneral:(FEMNavierStokesGeneral *)navierStokesGeneral
                      integration:(FEMNumericIntegration * __nonnull)integration
                    differentials:(FEMDifferentials * __nonnull)differentials
                   materialModels:(FEMMaterialModels * __nonnull)materialModels
                        utilities:(FEMUtilities * __nonnull)utilities
               elementDescription:(FEMElementDescription * __nonnull)elementDescription
                  timeIntegration:(FEMTimeIntegration *)timeIntegration
                        crsMatrix:(FEMMatrixCRS *)crsMatrix
                       bandMatrix:(FEMMatrixBand *)bandMatrix
                             bfID:(int)bf_id
                           bodyID:(int)body_id
                             eqID:(int)eq_id
                            matID:(int)mat_id
                           vector:(listBuffer * __nonnull)vector
                           matrix:(listBuffer * __nonnull)matrix
                             pwrk:(listBuffer * __nonnull)pwrk
                      meshVeloSol:(FEMVariable * __nonnull)meshVeloSol
                     meshVeloPerm:(int * __nullable)meshVeloPerm
                     meshVelocity:(double * __nullable)meshVelocity
                          tempSol:(FEMVariable * __nullable)tempSol
                         tempPerm:(int * __nullable)tempPerm
                      temperature:(double * __nullable)temperature
                         tempPrev:(double * __nullable)tempPrev
                      ifTransient:(BOOL)ifTransient
                       convection:(BOOL)convection
                   flowContainers:(variableArraysContainer * __nonnull)flowContainers
                       densitySol:(FEMVariable * __nullable)densitySol
                densityContainers:(variableArraysContainer * __nullable)densityContainers
                               dt:(int)dt
                          gravity:(double *)gravity
                         timeStep:(int)timeStep
              newtonLinearization:(BOOL)newtonLinearization
                        transient:(BOOL)transient
                    stabilizeFlag:(NSString * __nonnull)stabilizeFlag
                divDiscretization:(BOOL)divDiscretization
              gradPDiscretization:(BOOL)gradPDiscretization
                        stabilize:(BOOL)stabilize
                          bubbles:(BOOL)bubbles
                 getActiveElement:(Element_t* (* __nonnull)(id, SEL, int, FEMSolution*, FEMModel*))getActiveElementIMP
          getEquationIDForElement:(int (* __nonnull)(id, SEL, Element_t*, FEMModel*))getEquationIDForElementIMP
         getBodyForceIDForElement:(int (* __nonnull)(id, SEL, Element_t*, FEMModel*))getBodyForceIDForElementIMP
          getMaterialIDForElement:(int (* __nonnull)(id, SEL, Element_t*, FEMModel*))getMaterialIDForElementIMP
                    listGetString:(NSString* (* __nonnull)(id, SEL, FEMModel*, NSArray*, NSString*, BOOL*))listGetStringIMP
                   listGetLogical:(BOOL (* __nonnull)(id, SEL, FEMModel*, NSArray*, NSString*, BOOL*))listGetLogicalIMP
            listGetConstRealArray:(BOOL (* __nonnull)(id, SEL, FEMModel*, NSArray*, NSString*, listBuffer*))listGetConstRealArrayIMP
     getNumberOfBubbleDofsElement:(int (* __nonnull)(id, SEL, Element_t*, FEMSolution*))getNumberOfBubbleDofsElementIMP
           getElementDofsSolution:(int (* __nonnull)(id, SEL, FEMSolution*, FEMModel*, Element_t*, int*, BOOL*))getElementDofsSolutionIMP
                         getNodes:(void (* __nonnull)(id, SEL, FEMSolution*, FEMModel*, Element_t*, Nodes_t*, int*, FEMMesh*))getNodesIMP
                          getReal:(BOOL (* __nonnull)(id, SEL, FEMModel*, Element_t*, NSArray*, NSString*, listBuffer*, FEMListUtilities*))getRealIMP
                 listGetConstReal:(double (* __nonnull)(id, SEL, FEMModel*, NSArray*, NSString*, BOOL*, double*, double*))listGetConstRealIMP
                 listGetRealArray:(BOOL (* __nonnull)(id, SEL, FEMModel*, NSArray*, NSString*, int, int*, listBuffer*))listGetRealArrayIMP
    navierStokesComposeMassMatrix:(void (* __nonnull)(id, SEL, double**, double**, double*, double**, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, BOOL, NSString*, int, BOOL, double*, double*, BOOL, double**, BOOL, double*, double*, BOOL, BOOL, double*, BOOL, BOOL, BOOL, BOOL, Element_t*, int, int, int, Nodes_t*, FEMSolution*, FEMCore*, FEMMesh*, FEMModel*, FEMNumericIntegration*, FEMMaterial*, FEMElementDescription*, FEMCoordinateSystems*, FEMMaterialModels*, FEMDifferentials*, FEMListUtilities*, FEMUtilities*))navierStokesComposeMassMatrixIMP
navierStokesCylindricalComposeMassMatrix:(void (* __nonnull)(id, SEL, double**, double**, double*, double**, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, BOOL, NSString*, BOOL, BOOL, double*, double*, BOOL, double**, BOOL, double*, double*, BOOL, BOOL, BOOL, BOOL, Element_t*, int, int, int, Nodes_t*, FEMCore*, FEMMesh*, FEMModel*, FEMNumericIntegration*, FEMElementDescription*, FEMCoordinateSystems*, FEMMaterialModels*, FEMDifferentials*, FEMListUtilities*, FEMUtilities*))navierStokesCylindricalComposeMassMatrixIMP
navierStokesGeneralComposeMassMatrix:(void (* __nonnull)(id, SEL, double**, double**, double*, double**, double*, double*, double*, double*, double*, double*, double*, double*, BOOL, BOOL, Element_t*, int, Nodes_t*, FEMCore*, FEMMesh*, FEMModel*, FEMNumericIntegration*, FEMElementDescription*, FEMCoordinateSystems*, FEMMaterialModels*, FEMDifferentials*, FEMListUtilities*, FEMUtilities*))navierStokesGeneralComposeMassMatrixIMP
            defaultFirstOrderTime:(void (* __nonnull)(id, SEL, FEMModel*, FEMSolution*, Element_t *, double**, double**, double*, int*, int*, FEMTimeIntegration*, FEMUtilities*))defaultFirstOrderTimeIMP
                nsCondensateStiff:(void (* __nonnull)(id, SEL, double**, double*, int, int, int, double*))nsCondensateStiffIMP
           defaultUpdateEquations:(void (* __nonnull)(id, SEL, FEMModel*, FEMSolution*, Element_t *, double**, double*, int*, int*, FEMMatrixCRS*, FEMMatrixBand*))defaultUpdateEquationsIMP;
@end

@implementation FEMFlowSolution {
    
    BOOL _allocationDone;
    int _cols;
    int _localNodes;
    int _nsdofs;
    int _rows;
    int _sizePDensity0;
    int _sizePDensity1;
    int _sizePseudoPressure;
    int * __nullable _flowPerm;
    int * __nullable _indexes;
    double * __nullable _alpha;
    double * __nullable _beta;
    double * __nullable _density;
    double * __nullable * __nullable _drag;
    double * __nullable _extPressure;
    double * __nullable _flowSolution;
    double * __nullable _force;
    double * __nullable _gasConstant;
    double * __nullable _heatExpansionCoeff;
    double * __nullable _heatCapacity;
    double * __nullable _layerThickness;
    double * __nullable _localTempPrev;
    double * __nullable _localTemperature;
    double * __nullable * __nullable _loadVector;
    double * __nullable * __nullable _mass;
    double * __nullable _mu;
    double * __nullable _mv;
    double * __nullable _mw;
    double * __nullable _mx;
    double * __nullable _my;
    double * __nullable _mz;
    double * __nullable _pDensity0;
    double * __nullable _pDensity1;
    double * __nullable _permeability;
    double * __nullable _potentialCoefficient;
    double *__nullable _potentialField;
    double *__nullable _pressure;
    double * __nullable _prevDensity;
    double * __nullable _prevPressure;
    double * __nullable _pseudoCompressibility;
    double * __nullable _pseudoPressure;
    double * __nullable _pSolution;
    double * __nullable _referenceTemperature;
    double * __nullable * __nullable _slipCoeff;
    double * __nullable * __nullable _stiff;
    double * __nullable _timeForce;
    double * __nullable _surfaceRoughness;
    double * __nullable _u;
    double * __nullable _v;
    double * __nullable _viscosity;
    double * __nullable _w;
    Nodes_t * __nullable _elementNodes;
    
    // The following definitions are used the
    // for assembly on the GPU
    cl_device_id _device;
    cl_context _context;
    cl_command_queue _cmd_queue;
    cl_program _program;
    cl_int _err;
    char *_kernel_source;
    size_t _src_len;
    BOOL _initializeGPU;
}

-(void)FEMFlowSolution_nullify {
    
    _indexes = NULL;
    _alpha = NULL;
    _beta = NULL;
    _density = NULL;
    _drag = NULL;
    _extPressure = NULL;
    _force = NULL;
    _gasConstant = NULL;
    _heatExpansionCoeff = NULL;
    _heatCapacity = NULL;
    _layerThickness = NULL;
    _localTempPrev = NULL;
    _localTemperature = NULL;
    _loadVector = NULL;
    _mass = NULL;
    _mu = NULL;
    _mv = NULL;
    _mw = NULL;
    _mx = NULL;
    _my = NULL;
    _mz = NULL;
    _permeability = NULL;
    _potentialCoefficient = NULL;
    _potentialField = NULL;
    _prevDensity = NULL;
    _pressure = NULL;
    _prevPressure = NULL;
    _pseudoCompressibility = NULL;
    _pSolution = NULL;
    _referenceTemperature = NULL;
    _slipCoeff = NULL;
    _stiff = NULL;
    _timeForce = NULL;
    _surfaceRoughness = NULL;
    _u = NULL;
    _v = NULL;
    _viscosity = NULL;
    _w = NULL;
    if (_elementNodes != NULL) {
        _elementNodes->x = NULL;
        _elementNodes->y = NULL;
        _elementNodes->z = NULL;
    }
    _elementNodes = NULL;
}

-(void)FEMFlowSolution_checkCircleBoundaryModel:(FEMModel * __nonnull)model {
    
    int l=0;
    double phi, r, x, y, x0, y0;
    BOOL found;
    Element_t *elements;
    Nodes_t *nodes;
    FEMListUtilities *listUtilities = [FEMListUtilities sharedListUtilities];
    
    FEMMesh *mesh = (FEMMesh *)model.mesh;
    elements = mesh.getElements;
    nodes = mesh.getNodes;
    
    for (FEMBoundaryCondition *boundaryCondition in model.boundaryConditions) {
        if ([listUtilities listGetLogical:model inArray:boundaryCondition.valuesList forVariable:@"circle boundary" info:&found] == NO) continue;
        
        x0 = [listUtilities listGetConstReal:model inArray:boundaryCondition.valuesList forVariable:@"circle x" info:&found minValue:NULL maxValue:NULL];
        if (found == NO) x0 = 0.0;
        
        y0 = [listUtilities listGetConstReal:model inArray:boundaryCondition.valuesList forVariable:@"circle y" info:&found minValue:NULL maxValue:NULL];
        if (found == NO) y0 = 0.0;
        
        r = [listUtilities listGetConstReal:model inArray:boundaryCondition.valuesList forVariable:@"circle r" info:&found minValue:NULL maxValue:NULL];
        if (found == NO) r = 0.0;
        
        for (int i=mesh.numberOfBulkElements; i<mesh.numberOfBulkElements+mesh.numberOfBoundaryElements; i++) {
            if (elements[i].BoundaryInfo->Constraint != boundaryCondition.tag) continue;
            
            for (int k=0; k<elements[i].Type.NumberOfNodes; k++) {
                x = nodes->x[elements[i].NodeIndexes[k]] - x0;
                y = nodes->y[elements[i].NodeIndexes[k]] - y0;
                
                phi = atan2(y, x);
                x = r * cos(phi);
                y = r * sin(phi);
                
                nodes->x[elements[i].NodeIndexes[k]] = x + x0;
                nodes->y[elements[i].NodeIndexes[k]] = y + y0;
            }
            l++;
        }
    }
    
    if (l > 0) {
        NSLog(@"FEMFlowSolution:FEMFlowSolution_checkCircleBoundaryModel: number of elements on circle: %d.\n", l);
    }
}

-(void)FEMFlowSolution_doAssembly:(FEMSolution * __nonnull)solution
                             core:(FEMCore * __nonnull)core
                            model:(FEMModel * __nonnull)model
                    listUtilities:(FEMListUtilities * __nonnull)listUtilities
               coordinatesSystems:(FEMCoordinateSystems * __nonnull)coordinatesSystems
                     navierStokes:(FEMNavierStokes * __nonnull)navierStokes
          navierStokesCylindrical:(FEMNavierStokesCylindrical * __nonnull)navierStokesCylindrical
              navierStokesGeneral:(FEMNavierStokesGeneral *)navierStokesGeneral
                      integration:(FEMNumericIntegration * __nonnull)integration
                    differentials:(FEMDifferentials * __nonnull)differentials
                   materialModels:(FEMMaterialModels * __nonnull)materialModels
                        utilities:(FEMUtilities * __nonnull)utilities
               elementDescription:(FEMElementDescription * __nonnull)elementDescription
                  timeIntegration:(FEMTimeIntegration *)timeIntegration
                        crsMatrix:(FEMMatrixCRS *)crsMatrix
                       bandMatrix:(FEMMatrixBand *)bandMatrix
                             bfID:(int)bf_id
                           bodyID:(int)body_id
                             eqID:(int)eq_id
                            matID:(int)mat_id
                           vector:(listBuffer * __nonnull)vector
                           matrix:(listBuffer * __nonnull)matrix
                             pwrk:(listBuffer * __nonnull)pwrk
                      meshVeloSol:(FEMVariable * __nonnull)meshVeloSol
                     meshVeloPerm:(int * __nullable)meshVeloPerm
                     meshVelocity:(double * __nullable)meshVelocity
                          tempSol:(FEMVariable * __nullable)tempSol
                         tempPerm:(int * __nullable)tempPerm
                      temperature:(double * __nullable)temperature
                         tempPrev:(double * __nullable)tempPrev
                      ifTransient:(BOOL)ifTransient
                       convection:(BOOL)convection
                   flowContainers:(variableArraysContainer * __nonnull)flowContainers
                       densitySol:(FEMVariable * __nullable)densitySol
                densityContainers:(variableArraysContainer * __nullable)densityContainers
                               dt:(int)dt
                          gravity:(double *)gravity
                         timeStep:(int)timeStep
              newtonLinearization:(BOOL)newtonLinearization
                        transient:(BOOL)transient
                    stabilizeFlag:(NSString * __nonnull)stabilizeFlag
                divDiscretization:(BOOL)divDiscretization
              gradPDiscretization:(BOOL)gradPDiscretization
                        stabilize:(BOOL)stabilize
                          bubbles:(BOOL)bubbles
                 getActiveElement:(Element_t* (* __nonnull)(id, SEL, int, FEMSolution*, FEMModel*))getActiveElementIMP
          getEquationIDForElement:(int (* __nonnull)(id, SEL, Element_t*, FEMModel*))getEquationIDForElementIMP
         getBodyForceIDForElement:(int (* __nonnull)(id, SEL, Element_t*, FEMModel*))getBodyForceIDForElementIMP
          getMaterialIDForElement:(int (* __nonnull)(id, SEL, Element_t*, FEMModel*))getMaterialIDForElementIMP
                    listGetString:(NSString* (* __nonnull)(id, SEL, FEMModel*, NSArray*, NSString*, BOOL*))listGetStringIMP
                   listGetLogical:(BOOL (* __nonnull)(id, SEL, FEMModel*, NSArray*, NSString*, BOOL*))listGetLogicalIMP
            listGetConstRealArray:(BOOL (* __nonnull)(id, SEL, FEMModel*, NSArray*, NSString*, listBuffer*))listGetConstRealArrayIMP
     getNumberOfBubbleDofsElement:(int (* __nonnull)(id, SEL, Element_t*, FEMSolution*))getNumberOfBubbleDofsElementIMP
           getElementDofsSolution:(int (* __nonnull)(id, SEL, FEMSolution*, FEMModel*, Element_t*, int*, BOOL*))getElementDofsSolutionIMP
                         getNodes:(void (* __nonnull)(id, SEL, FEMSolution*, FEMModel*, Element_t*, Nodes_t*, int*, FEMMesh*))getNodesIMP
                          getReal:(BOOL (* __nonnull)(id, SEL, FEMModel*, Element_t*, NSArray*, NSString*, listBuffer*, FEMListUtilities*))getRealIMP
                 listGetConstReal:(double (* __nonnull)(id, SEL, FEMModel*, NSArray*, NSString*, BOOL*, double*, double*))listGetConstRealIMP
                 listGetRealArray:(BOOL (* __nonnull)(id, SEL, FEMModel*, NSArray*, NSString*, int, int*, listBuffer*))listGetRealArrayIMP
    navierStokesComposeMassMatrix:(void (* __nonnull)(id, SEL, double**, double**, double*, double**, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, BOOL, NSString*, int, BOOL, double*, double*, BOOL, double**, BOOL, double*, double*, BOOL, BOOL, double*, BOOL, BOOL, BOOL, BOOL, Element_t*, int, int, int, Nodes_t*, FEMSolution*, FEMCore*, FEMMesh*, FEMModel*, FEMNumericIntegration*, FEMMaterial*, FEMElementDescription*, FEMCoordinateSystems*, FEMMaterialModels*, FEMDifferentials*, FEMListUtilities*, FEMUtilities*))navierStokesComposeMassMatrixIMP
navierStokesCylindricalComposeMassMatrix:(void (* __nonnull)(id, SEL, double**, double**, double*, double**, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, BOOL, NSString*, BOOL, BOOL, double*, double*, BOOL, double**, BOOL, double*, double*, BOOL, BOOL, BOOL, BOOL, Element_t*, int, int, int, Nodes_t*, FEMCore*, FEMMesh*, FEMModel*, FEMNumericIntegration*, FEMElementDescription*, FEMCoordinateSystems*, FEMMaterialModels*, FEMDifferentials*, FEMListUtilities*, FEMUtilities*))navierStokesCylindricalComposeMassMatrixIMP
navierStokesGeneralComposeMassMatrix:(void (* __nonnull)(id, SEL, double**, double**, double*, double**, double*, double*, double*, double*, double*, double*, double*, double*, BOOL, BOOL, Element_t*, int, Nodes_t*, FEMCore*, FEMMesh*, FEMModel*, FEMNumericIntegration*, FEMElementDescription*, FEMCoordinateSystems*, FEMMaterialModels*, FEMDifferentials*, FEMListUtilities*, FEMUtilities*))navierStokesGeneralComposeMassMatrixIMP
            defaultFirstOrderTime:(void (* __nonnull)(id, SEL, FEMModel*, FEMSolution*, Element_t *, double**, double**, double*, int*, int*, FEMTimeIntegration*, FEMUtilities*))defaultFirstOrderTimeIMP
                nsCondensateStiff:(void (* __nonnull)(id, SEL, double**, double*, int, int, int, double*))nsCondensateStiffIMP
           defaultUpdateEquations:(void (* __nonnull)(id, SEL, FEMModel*, FEMSolution*, Element_t *, double**, double*, int*, int*, FEMMatrixCRS*, FEMMatrixBand*))defaultUpdateEquationsIMP {
    
    int i, j, k, compressibilityModel=-1, n, nb, nd;
    double angularVelocity[3], referencePressure, specificHeatRatio, tDiff;
    NSString *compressibilityFlag, *flowModel;
    Element_t * element = NULL;
    FEMEquation *equationAtID = nil;
    FEMBodyForce *bodyForceAtID = nil;
    FEMMaterial *materialAtID = nil;
    BOOL convect, found, hydrostatic = NO, magneticForce = NO, porous, potentialForce, pseudoCompressible, rotating;
    
    convect = convection;
    
    for (int t=0; t<solution.numberOfActiveElements; t++) {
        
        advanceOutput(t, solution.numberOfActiveElements, NULL, NULL);
        
        element = getActiveElementIMP(core, @selector(getActiveElement:solution:model:), t, solution, model);
        pseudoCompressible = NO;
        rotating = NO;
        if (element->BodyID != body_id) {
            body_id = element->BodyID;
            
            eq_id = getEquationIDForElementIMP(core, @selector(getEquationIDForElement:model:), element, model);
            equationAtID = (model.equations)[eq_id-1];
            
            bf_id = getBodyForceIDForElementIMP(core, @selector(getBodyForceIDForElement:model:), element, model);
            bodyForceAtID = (model.bodyForces)[bf_id-1];
            
            if ([flowModel isEqualToString:@"full"] == YES) {
                convect = listGetLogicalIMP(listUtilities, @selector(listGetLogical:inArray:forVariable:info:), model, equationAtID.valuesList, @"ns convect", &found);
                if (found == NO) convect = YES;
            }
            
            mat_id = getMaterialIDForElementIMP(core, @selector(getMaterialIDForElement:model:), element, model);
            materialAtID = (model.materials)[mat_id-1];
            
            compressibilityFlag = listGetStringIMP(listUtilities, @selector(listGetString:inArray:forVariable:info:), model, materialAtID.valuesList, @"compressibility model", &found);
            if (found == NO) compressibilityFlag = incompressible;
            
            if ([compressibilityFlag isEqualToString:@"incompressible"] == YES) {
                compressibilityModel = incompressible;
            }
            else if ([compressibilityFlag isEqualToString:@"perfect gas"] == YES || [compressibilityFlag isEqualToString:@"perfect gas equation 1"] == YES) {
                compressibilityModel = perfect_gas1;
            }
            else if ([compressibilityFlag isEqualToString:@"thermal"] == YES) {
                compressibilityModel = thermal;
            }
            else if ([compressibilityFlag isEqualToString:@"user defined"] == YES) {
                compressibilityModel = user_defined1;
            }
            else if ([compressibilityFlag isEqualToString:@"pressure dependent"] == YES) {
                compressibilityModel = user_defined2;
            }
            else if ([compressibilityFlag isEqualToString:@"artificial compressible"] == YES) {
                compressibilityModel = incompressible;
                pseudoCompressible = YES;
            }
            else {
                compressibilityModel = incompressible;
            }
            
            if (bf_id > 0) {
                magneticForce = listGetLogicalIMP(listUtilities, @selector(listGetLogical:inArray:forVariable:info:), model, bodyForceAtID.valuesList, @"lorentz force", &found);
                hydrostatic = listGetLogicalIMP(listUtilities, @selector(listGetLogical:inArray:forVariable:info:), model, bodyForceAtID.valuesList, @"hydrostatic pressure", &found);
            }
            if (found == NO) {
                hydrostatic = listGetLogicalIMP(listUtilities, @selector(listGetLogical:inArray:forVariable:info:), model, equationAtID.valuesList, @"hydrostatic pressure", &found);
            }
            
            if (bf_id > 0) {
                found = listGetConstRealArrayIMP(listUtilities, @selector(listGetConstRealArray:inArray:forVariable:buffer:), model, bodyForceAtID.valuesList, @"angular velocity", matrix);
                if (found == YES) {
                    if (coordinatesSystems.coordinates == cartesian) {
                        for (i=0; i<3; i++) {
                            angularVelocity[i] = matrix->matrix[i][0];
                        }
                        rotating = YES;
                    }
                } else {
                    memset(angularVelocity, 0.0, sizeof(angularVelocity) );
                }
            }
        }
        
        n = element->Type.NumberOfNodes;
        nb = getNumberOfBubbleDofsElementIMP(core, @selector(getNumberOfBubbleDofsElement:solution:), element, solution);
        nd = getElementDofsSolutionIMP(core, @selector(getElementDofsSolution:model:forElement:atIndexes:disableDiscontinuousGalerkin:), solution, model, element, _indexes, NULL);
        
        getNodesIMP(core, @selector(getNodes:model:inElement:resultNodes:numberOfNodes:mesh:), solution, model, element, _elementNodes, NULL, nil);
        
        switch (_nsdofs) {
            case 3:
                for (i=0; i<nd; i++) {
                    _u[i] = _flowSolution[_nsdofs*_flowPerm[_indexes[i]]];
                    _v[i] = _flowSolution[_nsdofs*_flowPerm[_indexes[i]]+1];
                    _w[i] = 0.0;
                }
                break;
            case 4:
                for (i=0; i<nd; i++) {
                    _u[i] = _flowSolution[_nsdofs*_flowPerm[_indexes[i]]];
                    _v[i] = _flowSolution[_nsdofs*_flowPerm[_indexes[i]]+1];
                    _w[i] = _flowSolution[_nsdofs*_flowPerm[_indexes[i]]+2];
                }
                break;
        }
        
        memset(_mu, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
        memset(_mv, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
        memset(_mw, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
        if (meshVeloSol != nil) {
            BOOL all;
            switch (meshVeloSol.dofs) {
                case 2:
                    all = YES;
                    for (i=0; i<nd; i++) {
                        if (meshVeloPerm[_indexes[i]] < 0) {
                            all = NO;
                            break;
                        }
                    }
                    if (all == YES) {
                        for (i=0; i<nd; i++) {
                            _mu[i] = meshVelocity[2*meshVeloPerm[_indexes[i]]];
                            _mv[i] = meshVelocity[2*meshVeloPerm[_indexes[i]]+1];
                        }
                    }
                    break;
                case 3:
                    all = YES;
                    for (i=0; i<n; i++) {
                        if (meshVeloPerm[element->NodeIndexes[i]] < 0) {
                            all = NO;
                            break;
                        }
                    }
                    if (all == YES) {
                        for (i=0; i<nd; i++) {
                            _mu[i] = meshVelocity[3*meshVeloPerm[_indexes[i]]];
                            _mv[i] = meshVelocity[3*meshVeloPerm[_indexes[i]]+1];
                            _mw[i] = meshVelocity[3*meshVeloPerm[_indexes[i]]+2];
                        }
                    }
                    break;
            }
        }
        
        memset(_localTemperature, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
        memset(_localTempPrev, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
        if (tempSol != nil) {
            BOOL all = YES;
            for (i=0; i<n; i++) {
                if (tempPerm[element->NodeIndexes[i]] < 0) {
                    all = NO;
                    break;
                }
            }
            if (all == YES) {
                for (i=0; i<nd; i++) {
                    _localTemperature[i] = temperature[tempPerm[_indexes[i]]];
                    if (ifTransient == YES && compressibilityModel != incompressible) {
                        _localTempPrev[i] = tempPrev[tempPerm[_indexes[i]]];
                    }
                }
            }
        }
        
        referencePressure = 0.0;
        memset(_prevDensity, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
        memset(_density, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
        switch (compressibilityModel) {
            case incompressible:
                for (i=0; i<nd; i++) {
                    if (_nsdofs == 3) {
                        _pressure[i] = _flowSolution[_nsdofs*_flowPerm[_indexes[i]]+2];
                    } else if (_nsdofs == 4) {
                        _pressure[i] = _flowSolution[_nsdofs*_flowPerm[_indexes[i]]+3];
                    }
                }
                found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, materialAtID.valuesList, @"density", vector, listUtilities);
                if (found == YES) memcpy(_density, vector->vector, n*sizeof(double));
                
                if (pseudoCompressible == YES) {
                    found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, materialAtID.valuesList, @"artificial pressure", vector, listUtilities);
                    if (found == NO) {
                        for (i=0; i<nd; i++) {
                            _pressure[i] = _pseudoPressure[_flowPerm[_indexes[i]]];
                        }
                    } else {
                        memcpy(_pressure, vector->vector, n*sizeof(double));
                        for (i=n; i<nd; i++) {
                            _pressure[i] = 0.0;
                        }
                    }
                }
                break;
                
            case perfect_gas1:
                // Use ReferenceTemperature in the MDF file for the fixed temperature
                // field. At the moment, can not have both fixed T ideal gas and
                // Boussinesq force:
                if (tempSol == nil) {
                    found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, materialAtID.valuesList, @"reference temperature", vector, listUtilities);
                    if (found == YES) memcpy(_localTemperature, vector->vector, n*sizeof(double));
                    memcpy(_localTempPrev, _localTemperature, n*sizeof(double));
                }
                found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, materialAtID.valuesList, @"heat capacity", vector, listUtilities);
                if (found == YES) memcpy(_heatCapacity, vector->vector, n*sizeof(double));
                
                // Read specific heat ratio
                specificHeatRatio = listGetConstRealIMP(listUtilities, @selector(listGetConstReal:inArray:forVariable:info:minValue:maxValue:), model, materialAtID.valuesList, @"specific heat ratio", &found, NULL, NULL);
                if (found == NO) specificHeatRatio = 5.0/3.0;
                
                // For an ideal gas, \gamma, c_p and R are really a constant.
                // GasConstant is an array only since HeatCapacity formally is
                for (i=0; i<n; i++) {
                    _gasConstant[i] = (specificHeatRatio - 1.0) * _heatCapacity[i] / specificHeatRatio;
                }
                
                // For ideal gases, take pressure deviation p_d as the dependent variable:
                // p = p_0 + p_d
                // Read p_0
                referencePressure = listGetConstRealIMP(listUtilities, @selector(listGetConstReal:inArray:forVariable:info:minValue:maxValue:), model, materialAtID.valuesList, @"reference pressure", &found, NULL, NULL);
                if (found == NO) referencePressure = 0.0;
                
                for (i=0; i<nd; i++) {
                    if (_nsdofs == 3) {
                        _pressure[i] = _flowSolution[_nsdofs*_flowPerm[_indexes[i]]+2];
                    } else if (_nsdofs == 4) {
                        _pressure[i] = _flowSolution[_nsdofs*_flowPerm[_indexes[i]]+3];
                    }
                }
                if (ifTransient == YES) {
                    for (i=0; i<nd; i++) {
                        if (_nsdofs == 3) {
                            _prevPressure[i] = flowContainers->PrevValues[_nsdofs*_flowPerm[_indexes[i]]+2][0];
                        } else if (_nsdofs == 4) {
                            _prevPressure[i] = flowContainers->PrevValues[_nsdofs*_flowPerm[_indexes[i]]+3][0];
                        }
                    }
                }
                for (i=0; i<n; i++) {
                    _density[i] = (_pressure[i] + referencePressure) / (_gasConstant[i] * _localTemperature[i]);
                }
                break;
                
            case user_defined1:
                for (i=0; i<nd; i++) {
                    if (_nsdofs == 3) {
                        _pressure[i] = _flowSolution[_nsdofs*_flowPerm[_indexes[i]]+2];
                    } else if (_nsdofs == 4) {
                        _pressure[i] = _flowSolution[_nsdofs*_flowPerm[_indexes[i]]+3];
                    }
                }
                if (densitySol != nil) {
                    for (i=0; i<nd; i++) {
                        _density[i] = densityContainers->Values[densityContainers->Perm[_indexes[i]]];
                    }
                    if (ifTransient == YES) {
                        for (i=0; i<nd; i++) {
                            _prevDensity[i] = densityContainers->PrevValues[densityContainers->Perm[_indexes[i]]][0];
                        }
                    }
                } else {
                    found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, materialAtID.valuesList, @"density", vector, listUtilities);
                    if (found == YES) memcpy(_density, vector->vector, n*sizeof(double));
                    if (ifTransient == YES) {
                        if (_pDensity0 == NULL) {
                            _pDensity0 = doublevec(0, _localNodes-1);
                            _pDensity1 = doublevec(0, _localNodes-1);
                            _sizePDensity0 = _localNodes;
                            _sizePDensity1 = _localNodes;
                        }
                        if (dt == 1) {
                            for (i=0; i<n; i++) {
                                _pDensity0[_indexes[i]] = _density[i];
                                _pDensity1[_indexes[i]] = _density[i];
                                _prevDensity[i] = _pDensity0[_indexes[i]];
                            }
                        }
                    }
                }
                break;
                
            case user_defined2:
                found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, materialAtID.valuesList, @"density", vector, listUtilities);
                if (found == YES) memcpy(_density, vector->vector, n*sizeof(double));
                
                for (i=0; i<nd; i++) {
                    if (_nsdofs == 3) {
                        _pressure[i] = _flowSolution[_nsdofs*_flowPerm[_indexes[i]]+2];
                    } else if (_nsdofs == 4) {
                        _pressure[i] = _flowSolution[_nsdofs*_flowPerm[_indexes[i]]+3];
                    }
                }
                break;
                
            case thermal:
                for (i=0; i<nd; i++) {
                    if (_nsdofs == 3) {
                        _pressure[i] = _flowSolution[_nsdofs*_flowPerm[_indexes[i]]+2];
                    } else if (_nsdofs == 4) {
                        _pressure[i] = _flowSolution[_nsdofs*_flowPerm[_indexes[i]]+3];
                    }
                }
                
                found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, materialAtID.valuesList, @"heat expansion coefficient", vector, listUtilities);
                if (found == YES) memcpy(_heatExpansionCoeff, vector->vector, n*sizeof(double));
                
                found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, materialAtID.valuesList, @"reference temperature", vector, listUtilities);
                if (found == YES) memcpy(_referenceTemperature, vector->vector, n*sizeof(double));
                
                found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, materialAtID.valuesList, @"density", vector, listUtilities);
                if (found == YES) memcpy(_density, vector->vector, n*sizeof(double));
                
                if (ifTransient) {
                    for (i=0; i<n; i++) {
                        _prevDensity[i] = _density[i] * ( 1.0 - _heatExpansionCoeff[i] * (_localTempPrev[i] - _referenceTemperature[i]) );
                    }
                }
                for (i=0; i<n; i++) {
                    _density[i] = _density[i] * ( 1.0 - _heatExpansionCoeff[i] * (_localTemperature[i] - _referenceTemperature[i]) );
                }
                break;
        }
        
        // Read in porous media defs
        porous = listGetLogicalIMP(listUtilities, @selector(listGetLogical:inArray:forVariable:info:), model, materialAtID.valuesList, @"porous media", &found);
        if (porous == YES) {
            found = listGetRealArrayIMP(listUtilities, @selector(listGetRealArray:inArray:forVariable:numberOfNodes:indexes:buffer:), model, materialAtID.valuesList, @"porous resistivity", n, element->NodeIndexes, pwrk);
            if (found == NO) {
                found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, materialAtID.valuesList, @"porous resistivity 1", vector, listUtilities);
                if (found == YES) {
                    for (i=0; i<n; i++) {
                        _drag[0][i] = vector->vector[i];
                    }
                }
                found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, materialAtID.valuesList, @"porous resistivity 2", vector, listUtilities);
                if (found == YES) {
                    for (i=0; i<n; i++) {
                        _drag[1][i] = vector->vector[i];
                    }
                }
                if (_nsdofs-1 > 2) {
                    found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, materialAtID.valuesList, @"porous resistivity 3", vector, listUtilities);
                    if (found == YES) {
                        for (i=0; i<n; i++) {
                            _drag[2][i] = vector->vector[i];
                        }
                    }
                }
            } else if (pwrk->m == 1) {
                for (i=0; i<_nsdofs-1; i++) {
                    for (j=0; j<n; j++) {
                        _drag[i][j] = pwrk->tensor[0][0][j];
                    }
                }
            } else {
                for (i=0; i<min(_nsdofs, pwrk->m); i++) {
                    for (j=0; j<n; j++) {
                        _drag[i][j] = pwrk->tensor[i][0][j];
                    }
                }
            }
        }
        
        // Viscosity = Laminar viscosity
        found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, materialAtID.valuesList, @"viscosity", vector, listUtilities);
        if (found == YES) memcpy(_viscosity, vector->vector, n*sizeof(double));
        
        // Set body forces if any
        memset(*_loadVector, 0.0, (4*solution.mesh.maxElementDofs)*sizeof(double) );
        
        potentialForce = NO;
        if (bf_id > 0) {
            memset(_heatExpansionCoeff, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
            memset(_referenceTemperature, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
            
            // Boussinesq body force and gravity
            if (listGetLogicalIMP(listUtilities, @selector(listGetLogical:inArray:forVariable:info:), model, bodyForceAtID.valuesList, @"boussinesq", &found) == YES) {
                found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, materialAtID.valuesList, @"heat expansion coefficient", vector, listUtilities);
                if (found == YES) memcpy(_heatExpansionCoeff, vector->vector, n*sizeof(double));
                
                found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, materialAtID.valuesList, @"reference temperature", vector, listUtilities);
                if (found == YES) memcpy(_referenceTemperature, vector->vector, n*sizeof(double));
                
                for (i=0; i<n; i++) {
                    k = tempPerm[element->NodeIndexes[i]];
                    if (k >= 0) {
                        if (hydrostatic == YES) {
                            tDiff = 1.0 - _heatExpansionCoeff[i] * (temperature[k] - _referenceTemperature[i]);
                            
                            if (tDiff <= 0.0) {
                                NSLog(@"FEMFlowSolution:solutionComputer: zero or negative density.\n");
                            }
                        } else {
                            tDiff = -_heatExpansionCoeff[i] * (temperature[k] - _referenceTemperature[i]);
                        }
                        
                        _loadVector[0][i] = gravity[0] * tDiff;
                        _loadVector[1][i] = gravity[1] * tDiff;
                        if (_nsdofs > 3) _loadVector[2][i] = gravity[2] * tDiff;
                    }
                }
            } else if (hydrostatic == YES) {
                for (i=0; i<n; i++) {
                    _loadVector[0][i] = gravity[0];
                    _loadVector[1][i] = gravity[1];
                    if (_nsdofs > 3) _loadVector[2][i] = gravity[2];
                }
            }
            
            found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, bodyForceAtID.valuesList, @"flow bodyforce 1", vector, listUtilities);
            if (found == YES) {
                for (i=0; i<n; i++) {
                    _loadVector[0][i] = _loadVector[0][i] + vector->vector[i];
                }
            }
            found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, bodyForceAtID.valuesList, @"flow bodyforce 2", vector, listUtilities);
            if (found == YES) {
                for (i=0; i<n; i++) {
                    _loadVector[1][i] = _loadVector[1][i] + vector->vector[i];
                }
            }
            if (_nsdofs > 3) {
                found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, bodyForceAtID.valuesList, @"flow bodyforce 3", vector, listUtilities);
                if (found == YES) {
                    for (i=0; i<n; i++) {
                        _loadVector[2][i] = _loadVector[2][i] + vector->vector[i];
                    }
                }
            }
            
            potentialForce = listGetLogicalIMP(listUtilities, @selector(listGetLogical:inArray:forVariable:info:), model, bodyForceAtID.valuesList, @"potential force", &found);
            if (potentialForce == YES) {
                found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, bodyForceAtID.valuesList, @"potential field", vector, listUtilities);
                if (found == YES) memcpy(_potentialField, vector->vector, n*sizeof(double));
                
                found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, bodyForceAtID.valuesList, @"potential coefficient", vector, listUtilities);
                if (found == YES) memcpy(_potentialCoefficient, vector->vector, n*sizeof(double));
            }
        }
        
        // Note: LaoadVector is multiplied by density inside *Navier* routines
        if (ifTransient == YES) {
            switch (compressibilityModel) {
                case perfect_gas1:
                    if (tempSol != nil) {
                        for (i=0; i<n; i++) {
                            k = tempPerm[element->NodeIndexes[i]];
                            if (k >= 0) {
                                _loadVector[_nsdofs-1][i] = _loadVector[_nsdofs-1][i] + (temperature[k] - tempPrev[k]) / timeStep;
                            }
                        }
                    }
                    break;
                case user_defined1:
                case thermal:
                    for (i=0; i<n; i++) {
                        _loadVector[_nsdofs-1][i] = _loadVector[_nsdofs-1][i] - (_density[i] - _prevDensity[i]) / (_density[i]*timeStep);
                    }
                    break;
            }
        }
        
        // Get element local stiffness and mass matrices
        memset( *_stiff, 0.0, ((2*_nsdofs*solution.mesh.maxElementDofs)*(2*_nsdofs*solution.mesh.maxElementDofs))*sizeof(double) );
        memset( *_mass, 0.0, ((2*_nsdofs*solution.mesh.maxElementDofs)*(2*_nsdofs*solution.mesh.maxElementDofs))*sizeof(double) );
        memset( _force, 0.0, (2*_nsdofs*solution.mesh.maxElementDofs)*sizeof(double) );
        double nodalPressure[n];
        switch (coordinatesSystems.coordinates) {
            case cartesian:
                switch (compressibilityModel) {
                    case incompressible:
                    case perfect_gas1:
                    case user_defined1:
                    case user_defined2:
                    case thermal:
                        for (i=0; i<n; i++) {
                            nodalPressure[i] = referencePressure + _pressure[i];
                        }
                        navierStokesComposeMassMatrixIMP(navierStokes,
                                                         @selector(navierStokesComposeMassMatrix:stiffMatrix:forceVector:loadVector:nodalViscosity:nodalDensity:velocityX:velocityY:velocityZ:meshVelocityX:meshVelocityY:meshVelocityZ:nodalPressure:nodalTemperature:isConvect:stabilizeFlag:compressibilityModel:isPseudoCompressible:nodalCompressibility:nodalGasConstant:isPorous:nodalDrag:isPotentialForce:potentialField:potentialCoefficient:isMagneticForce:isRotating:omega:isDivDiscretization:isGradPDriscretization:isNewtonLinearization:isTransient:element:numberOfNodes:rows:cols:nodes:solution:core:mesh:model:integration:material:elementDescription:coordinateSystems:materialModels:differentials:listUtilities:utilities:),
                                                         _mass, _stiff, _force, _loadVector, _viscosity, _density, _u, _v, _w, _mu, _mv, _mw, nodalPressure, _localTemperature, convect, stabilizeFlag, compressibilityModel, pseudoCompressible, _pseudoCompressibility, _gasConstant, porous, _drag, potentialForce, _potentialField, _potentialCoefficient, magneticForce, rotating, angularVelocity, divDiscretization, gradPDiscretization, newtonLinearization, transient, element, n, _rows, _cols, _elementNodes, solution, core, solution.mesh, model, integration, materialAtID, elementDescription, coordinatesSystems,  materialModels, differentials, listUtilities, utilities);
                        break;
                }
                break;
                
            case cylindric:
            case cylindric_symmetric:
            case axis_symmetric:
                switch (compressibilityModel) {
                    case incompressible:
                    case perfect_gas1:
                        for (i=0; i<n; i++) {
                            nodalPressure[i] = referencePressure + _pressure[i];
                        }
                        navierStokesCylindricalComposeMassMatrixIMP(navierStokesCylindrical,
                                                                    @selector(navierStokesCylindricalComposeMassMatrix:stiffMatrix:forceVector:loadVector:nodalViscosity:nodalDensity:velocityX:velocityY:velocityZ:meshVelocityX:meshVelocityY:meshVelocityZ:nodalPressure:nodalTemperature:isConvect:stabilizeFlag:isCompressible:isPseudoCompressible:nodalCompressibility:nodalGasConstant:isPorous:nodalDrag:isPotentialForce:potentialField:potentialCoefficient:isMagneticForce:isDivDiscretization:isGradPDriscretization:isNewtonLinearization:element:numberOfNodes:rows:cols:nodes:core:mesh:model:integration:elementDescription:coordinateSystems:materialModels:differentials:listUtilities:utilities:),
                                                                    _mass, _stiff, _force, _loadVector, _viscosity, _density, _u, _v, _w, _mu, _mv, _mw, nodalPressure, _localTemperature, convect, stabilizeFlag, (compressibilityModel != incompressible) ? YES : NO, pseudoCompressible, _pseudoCompressibility, _gasConstant, porous, _drag, potentialForce, _potentialField, _potentialCoefficient, magneticForce, divDiscretization, gradPDiscretization, newtonLinearization, element, n, _rows, _cols, _elementNodes, core, solution.mesh, model, integration, elementDescription, coordinatesSystems, materialModels, differentials, listUtilities, utilities);
                        break;
                }
                break;
                
            default:
                
                switch (compressibilityModel) {
                    case incompressible:
                    case perfect_gas1:
                        navierStokesGeneralComposeMassMatrixIMP(navierStokesGeneral,
                                                                @selector(navierStokesGeneralComposeMassMatrix:stiffMatrix:forceVector:loadVector:nodalViscosity:nodalDensity:velocityX:velocityY:velocityZ:meshVelocityX:meshVelocityY:meshVelocityZ:isStabilize:isNewtonLinearization:element:numberOfNodes:nodes:core:mesh:model:integration:elementDescription:coordinateSystems:materialModels:differentials:listUtilities:utilities:),
                                                                _mass, _stiff, _force, _loadVector, _viscosity, _density, _u, _v, _w, _mu, _mv, _mw, stabilize, newtonLinearization, element, n, _elementNodes, core, solution.mesh, model, integration, elementDescription, coordinatesSystems, materialModels, differentials, listUtilities, utilities);
                        break;
                        
                    default:
                        fatal("FEMFlowSolution:solutionComputer", "Missing compresibility model in general coordinates.");
                        break;
                }
                break;
        }
        
        // If time dependent simulation, add mass matrix to global
        // matrix and global RHS vector
        if (compressibilityModel != incompressible && [stabilizeFlag isEqualToString:@"stabilized"] == YES) {
            bubbles = YES;
            stabilizeFlag = @"bubbles";
        }
        if (element->Type.BasisFunctionDegree <= 1 && [stabilizeFlag isEqualToString:@"p2/p1"] == YES) {
            bubbles = YES;
            stabilizeFlag = @"bubbles";
        }
        if (nb == 0 && bubbles == YES) nb = n;
        
        memset(_timeForce, 0.0, (2*_nsdofs*solution.mesh.maxElementDofs)*sizeof(double) );
        if (ifTransient) {
            // Note: The following will replace STIFF and FORCE with the combined information
            defaultFirstOrderTimeIMP(core, @selector(defaultFirstOrderTime:inSolution:forElement:realMass:realStiff:realForce:stiffRows:stiffCols:timeIntegration:utilities:), model, solution, element, _mass, _stiff, _force, &_rows, &_cols, timeIntegration, utilities);
        }
        
        if (nb > 0) {
            nsCondensateStiffIMP(core, @selector(nsCondensateStiff:force:numberOfNodes:numberOfBubbles:dimension:force1:), _stiff, _force, nd, nb, _nsdofs-1, _timeForce);
        }
        
        // Add local stiffness matrix and force vector to global matrix and vector
        defaultUpdateEquationsIMP(core, @selector(defaultUpdateEquations:inSolution:forElement:realStiff:realForce:stiffRows:stiffCols:crsMatrix:bandMatrix:), model, solution, element, _stiff, _force, &_rows, &_cols, crsMatrix, bandMatrix);
    }
}

- (id)init
{
    self = [super init];
    if (self) {
        _allocationDone = NO;
        _cols = 0;
        _rows = 0;
        
        _flowPerm = NULL;
        _flowSolution = NULL;
        
        [self FEMFlowSolution_nullify];
        
        _pseudoPressure = NULL;
        _pDensity0 = NULL;
        _pDensity1 = NULL;
        
        _device = NULL;
        _initializeGPU = NO;
}
    
    return self;
}

-(void)deallocation:(FEMSolution * __nonnull)solution {
    
    int n = solution.mesh.maxElementDofs;
    variableArraysContainer *flowContainers = solution.variable.getContainers;
    
    if (_u != NULL) free_dvector(_u, 0, n-1);
    if (_v != NULL) free_dvector(_v, 0, n-1);
    if (_w != NULL) free_dvector(_w, 0, n-1);
    if (_mu != NULL) free_dvector(_mu, 0, n-1);
    if (_mv != NULL) free_dvector(_mv, 0, n-1);
    if (_mw != NULL) free_dvector(_mw, 0, n-1);
    if (_indexes != NULL) free_ivector(_indexes, 0, n-1);
    if (_pressure != NULL) free_dvector(_pressure, 0, n-1);
    if (_prevPressure != NULL) free_dvector(_prevPressure, 0, n-1);
    if (_pseudoCompressibility != NULL) free_dvector(_pseudoCompressibility, 0, n-1);
    if (_prevDensity != NULL) free_dvector(_prevDensity, 0, n-1);
    if (_density != NULL) free_dvector(_density, 0, n-1);
    if (_layerThickness != NULL) free_dvector(_layerThickness, 0, n-1);
    if (_surfaceRoughness != NULL) free_dvector(_surfaceRoughness, 0, n-1);
    if (_permeability != NULL) free_dvector(_permeability, 0, n-1);
    if (_mx != NULL) free_dvector(_mx, 0, n-1);
    if (_my != NULL) free_dvector(_my, 0, n-1);
    if (_mz != NULL) free_dvector(_mz, 0, n-1);
    if (_slipCoeff != NULL) free_dmatrix(_slipCoeff, 0, 2, 0, n-1);
    if (_drag != NULL) free_dmatrix(_drag, 0, 2, 0, n-1);
    if (_timeForce != NULL) free_dvector(_timeForce, 0, (2*_nsdofs*n)-1);
    if (_force != NULL) free_dvector(_force, 0, (2*_nsdofs*n)-1);
    if (_viscosity != NULL) free_dvector(_viscosity, 0, n-1);
    if (_mass != NULL) free_dmatrix(_mass, 0, (2*_nsdofs*n)-1, 0, (2*_nsdofs*n)-1);
    if (_stiff != NULL) free_dmatrix(_stiff, 0, (2*_nsdofs*n)-1, 0, (2*_nsdofs*n)-1);
    if (_heatExpansionCoeff != NULL) free_dvector(_heatExpansionCoeff, 0, n-1);
    if (_gasConstant != NULL) free_dvector(_gasConstant, 0, n-1);
    if (_heatCapacity != NULL) free_dvector(_heatCapacity, 0, n-1);
    if (_referenceTemperature != NULL) free_dvector(_referenceTemperature, 0, n-1);
    if (_localTempPrev != NULL) free_dvector(_localTempPrev, 0, n-1);
    if (_localTemperature != NULL) free_dvector(_localTemperature, 0, n-1);
    if (_pSolution != NULL) free_dvector(_pSolution, 0, flowContainers->sizeValues-1);
    if (_potentialField != NULL) free_dvector(_potentialField, 0, n-1);
    if (_potentialCoefficient != NULL) free_dvector(_potentialCoefficient, 0, n-1);
    if (_loadVector != NULL) free_dmatrix(_loadVector, 0, 3, 0, n-1);
    if (_alpha != NULL) free_dvector(_alpha, 0, n-1);
    if (_beta != NULL) free_dvector(_beta, 0, n-1);
    if (_extPressure != NULL) free_dvector(_extPressure, 0, n-1);
    if (_elementNodes->x != NULL) free_dvector(_elementNodes->x, 0, n-1);
    if (_elementNodes->y != NULL) free_dvector(_elementNodes->y, 0, n-1);
    if (_elementNodes->z != NULL) free_dvector(_elementNodes->z, 0, n-1);
    if (_elementNodes != NULL) free(_elementNodes);
    
    if (_pseudoPressure != NULL) free_dvector(_pseudoPressure, 0, _sizePseudoPressure-1);
    if (_pDensity0 != NULL) free_dvector(_pDensity0, 0, _sizePDensity0);
    if (_pDensity1 != NULL) free_dvector(_pDensity1, 0, _sizePDensity1);
}

-(void)solutionComputer:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model timeStep:(int)timeStep transientSimulation:(BOOL)transient {
    
    int i, j, k, n, t, bf_id, body_id, eq_id, mat_id, compressibilityModel=-1, dim, freeSIter, iter, modelCoords=0, modelDim=0, newtonIter, nonLinearIter;
    static int dt, saveTimeStep=-1;
    int *tempPerm = NULL, *meshVeloPerm = NULL;
    double *temperature = NULL, *tempPrev = NULL, *meshVelocity = NULL;
    double at, at0, at1, freeSTol, gravity[3], nonLinearRelax, nonLinearTol, newtonTol, pseudoCompressibilityScale, relativeChange, relaxation, st, sum, totat, totst, uNorm;
    NSString *compressibilityFlag, *flowModel, *localCoords, *stabilizeFlag, *varName;
    BOOL bubbles, convect, computeFree = NO, divDiscretization, found, freeSurfaceFlag, gradPDiscretization, gotForceBC, ifTransient, mbFlag, newtonLinearization = NO, normalTangential, pseudoPressureExists=NO, pseudoPressureUpdate=NO, relaxBefore, stabilize, useLocalCoords = NO;
    NSArray *bc = nil;
    Element_t * element = NULL, *parent = NULL;
    matrixArraysContainer *matContainers = NULL;
    variableArraysContainer *densityContainers = NULL, *flowContainers = NULL, *tempSolContainers = NULL, *meshVeloSolContainers = NULL,
                            *timeVarContainers = NULL;
    listBuffer *vector = NULL;
    listBuffer *matrix = NULL;
    listBuffer *pwrk = NULL;
    FEMMesh *mesh;
    FEMVariable *densitySol, *tempSol, *meshVeloSol, *timeVar;
    FEMCore *core = [FEMCore sharedCore];
    FEMListUtilities *listUtilities = [FEMListUtilities sharedListUtilities];
    FEMUtilities *utilities = [[FEMUtilities alloc] init];
    FEMMatrixCRS *crsMatrix = [[FEMMatrixCRS alloc] init];
    FEMMatrixBand *bandMatrix = [[FEMMatrixBand alloc] init];
    FEMCoordinateSystems *coordinatesSystems = [[FEMCoordinateSystems alloc] init];
    FEMElementDescription *elementDescription = [FEMElementDescription sharedElementDescription];
    FEMDifferentials *differentials = [[FEMDifferentials alloc] init];
    FEMNavierStokes *navierStokes = [[FEMNavierStokes alloc] init];
    FEMNavierStokesCylindrical *navierStokesCylindrical = [[FEMNavierStokesCylindrical alloc] init];
    FEMNavierStokesGeneral *navierStokesGeneral = [[FEMNavierStokesGeneral alloc] init];
    FEMMaterialModels *materialModels = [[FEMMaterialModels alloc] init];
    FEMElementUtils *elementUtils = [[FEMElementUtils alloc] init];
    FEMTimeIntegration *timeIntegration;
    
    static Element_t* (*getActiveElementIMP)(id, SEL, int, FEMSolution*, FEMModel*) = nil;
    static int (*getEquationIDForElementIMP)(id, SEL, Element_t*, FEMModel*) = nil;
    static int (*getBodyForceIDForElementIMP)(id, SEL, Element_t*, FEMModel*) = nil;
    static BOOL (*listGetLogicalIMP)(id, SEL, FEMModel*, NSArray*, NSString*, BOOL*) = nil;
    static int (*getMaterialIDForElementIMP)(id, SEL, Element_t*, FEMModel*) = nil;
    static NSString* (*listGetStringIMP)(id, SEL, FEMModel*, NSArray*, NSString*, BOOL*) = nil;
    static BOOL (*listGetConstRealArrayIMP)(id, SEL, FEMModel*, NSArray*, NSString*, listBuffer*) = nil;
    static int (*getNumberOfBubbleDofsElementIMP)(id, SEL, Element_t*, FEMSolution*) = nil;
    static int (*getElementDofsSolutionIMP)(id, SEL, FEMSolution*, FEMModel*, Element_t*, int*, BOOL*) = nil;
    static void (*getNodesIMP)(id, SEL, FEMSolution*, FEMModel*, Element_t*, Nodes_t*, int*, FEMMesh*) = nil;
    static BOOL (*getRealIMP)(id, SEL, FEMModel*, Element_t*, NSArray*, NSString*, listBuffer*, FEMListUtilities*) = nil;
    static double (*listGetConstRealIMP)(id, SEL, FEMModel*, NSArray*, NSString*, BOOL*, double*, double*) = nil;
    static BOOL (*listGetRealArrayIMP)(id, SEL, FEMModel*, NSArray*, NSString*, int, int*, listBuffer*) = nil;
    static void (*defaultFirstOrderTimeIMP)(id, SEL, FEMModel*, FEMSolution*, Element_t *, double**, double**, double*, int*, int*, FEMTimeIntegration*, FEMUtilities*) = nil;
    static void (*nsCondensateStiffIMP)(id, SEL, double**, double*, int, int, int, double*) = nil;
    static void (*defaultUpdateEquationsIMP)(id, SEL, FEMModel*, FEMSolution*, Element_t *, double**, double*, int*, int*, FEMMatrixCRS*, FEMMatrixBand*) = nil;
    static void (*navierStokesComposeMassMatrixIMP)(id, SEL, double**, double**, double*, double**, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, BOOL, NSString*, int, BOOL, double*, double*, BOOL, double**, BOOL, double*, double*, BOOL, BOOL, double*, BOOL, BOOL, BOOL, BOOL, Element_t*, int, int, int, Nodes_t*, FEMSolution*, FEMCore*, FEMMesh*, FEMModel*, FEMNumericIntegration*, FEMMaterial*, FEMElementDescription*, FEMCoordinateSystems*, FEMMaterialModels*, FEMDifferentials*, FEMListUtilities*, FEMUtilities*) = nil;
    static void (*navierStokesBoundaryIMP)(id, SEL, double**, double*, double**, double*, double*, double*, double**, BOOL, Element_t*, int, Nodes_t*, FEMMesh*, FEMModel*, FEMNumericIntegration*, FEMElementDescription*, FEMElementUtils*) = nil;
    static void (*navierStokesCylindricalComposeMassMatrixIMP)(id, SEL, double**, double**, double*, double**, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, BOOL, NSString*, BOOL, BOOL, double*, double*, BOOL, double**, BOOL, double*, double*, BOOL, BOOL, BOOL, BOOL, Element_t*, int, int, int, Nodes_t*, FEMCore*, FEMMesh*, FEMModel*, FEMNumericIntegration*, FEMElementDescription*, FEMCoordinateSystems*, FEMMaterialModels*, FEMDifferentials*, FEMListUtilities*, FEMUtilities*) = nil;
     static void (*navierStokesCylindricalBoundaryIMP)(id, SEL, double**, double*, double**, double*, double*, double*, double**, BOOL, Element_t*, int, Nodes_t*, FEMMesh*, FEMModel*, FEMNumericIntegration*, FEMElementDescription*, FEMElementUtils*, FEMCoordinateSystems*) = nil;
    static void (*navierStokesGeneralComposeMassMatrixIMP)(id, SEL, double**, double**, double*, double**, double*, double*, double*, double*, double*, double*, double*, double*, BOOL, BOOL, Element_t*, int, Nodes_t*, FEMCore*, FEMMesh*, FEMModel*, FEMNumericIntegration*, FEMElementDescription*, FEMCoordinateSystems*, FEMMaterialModels*, FEMDifferentials*, FEMListUtilities*, FEMUtilities*) = nil;
    static void (*navierStokesGeneralBoundaryIMP)(id, SEL, double**, double*, double**, double*, double*, double*, double**, Element_t*, int, Nodes_t*, FEMMesh*, FEMModel*, FEMNumericIntegration*, FEMElementDescription*, FEMCoordinateSystems*) = nil;
    
    static dispatch_once_t onceToken;
    dispatch_once(&onceToken, ^{
        if (!getActiveElementIMP) {
            getActiveElementIMP = (Element_t* (*)(id, SEL, int, FEMSolution*, FEMModel*))
            [core methodForSelector: @selector(getActiveElement:solution:model:)];
        }
        if (!getEquationIDForElementIMP) {
            getEquationIDForElementIMP = (int (*)(id, SEL, Element_t*, FEMModel*))
            [core methodForSelector: @selector(getEquationIDForElement:model:)];
        }
        if (!getBodyForceIDForElementIMP) {
            getBodyForceIDForElementIMP = (int (*)(id, SEL, Element_t*, FEMModel*))
            [core methodForSelector: @selector(getBodyForceIDForElement:model:)];
        }
        if (!listGetLogicalIMP) {
            listGetLogicalIMP = (BOOL (*)(id, SEL, FEMModel*, NSArray*, NSString*, BOOL*))
            [listUtilities methodForSelector: @selector(listGetLogical:inArray:forVariable:info:)];
        }
        if (!getMaterialIDForElementIMP) {
            getMaterialIDForElementIMP = (int (*)(id, SEL, Element_t*, FEMModel*))
            [core methodForSelector: @selector(getMaterialIDForElement:model:)];
        }
        if (!listGetStringIMP) {
            listGetStringIMP = (NSString* (*)(id, SEL, FEMModel*, NSArray*, NSString*, BOOL*))
            [listUtilities methodForSelector: @selector(listGetString:inArray:forVariable:info:)];
        }
        if (!listGetConstRealArrayIMP) {
            listGetConstRealArrayIMP = (BOOL (*)(id, SEL, FEMModel*, NSArray*, NSString*, listBuffer*))
            [listUtilities methodForSelector: @selector(listGetConstRealArray:inArray:forVariable:buffer:)];
        }
        if (!getNumberOfBubbleDofsElementIMP) {
            getNumberOfBubbleDofsElementIMP = (int (*)(id, SEL, Element_t*, FEMSolution*))
            [core methodForSelector: @selector(getNumberOfBubbleDofsElement:solution:)];
        }
        if (!getElementDofsSolutionIMP) {
            getElementDofsSolutionIMP = (int (*)(id, SEL, FEMSolution*, FEMModel*, Element_t*, int*, BOOL*))
            [core methodForSelector: @selector(getElementDofsSolution:model:forElement:atIndexes:disableDiscontinuousGalerkin:)];
        }
        if (!getNodesIMP) {
            getNodesIMP = (void (*)(id, SEL, FEMSolution*, FEMModel*, Element_t*, Nodes_t*, int*, FEMMesh*))
            [core methodForSelector: @selector(getNodes:model:inElement:resultNodes:numberOfNodes:mesh:)];
        }
        if (!getRealIMP) {
            getRealIMP = (BOOL (*)(id, SEL, FEMModel*, Element_t*, NSArray*, NSString*, listBuffer*, FEMListUtilities*))
            [core methodForSelector: @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:)];
        }
        if (!listGetConstRealIMP) {
            listGetConstRealIMP = (double (*)(id, SEL, FEMModel*, NSArray*, NSString*, BOOL*, double*, double*))
            [listUtilities methodForSelector: @selector(listGetConstReal:inArray:forVariable:info:minValue:maxValue:)];
        }
        if (!listGetRealArrayIMP) {
            listGetRealArrayIMP = (BOOL (*)(id, SEL, FEMModel*, NSArray*, NSString*, int, int*, listBuffer*))
            [listUtilities methodForSelector: @selector(listGetRealArray:inArray:forVariable:numberOfNodes:indexes:buffer:)];
        }
        if (!defaultFirstOrderTimeIMP) {
            defaultFirstOrderTimeIMP = (void (*)(id, SEL, FEMModel*, FEMSolution*, Element_t *, double**, double**, double*, int*, int*, FEMTimeIntegration*, FEMUtilities*))
            [core methodForSelector: @selector(defaultFirstOrderTime:inSolution:forElement:realMass:realStiff:realForce:stiffRows:stiffCols:timeIntegration:utilities:)];
        }
        if (!nsCondensateStiffIMP) {
            nsCondensateStiffIMP = (void (*)(id, SEL, double**, double*, int, int, int, double*))
            [core methodForSelector: @selector(nsCondensateStiff:force:numberOfNodes:numberOfBubbles:dimension:force1:)];
        }
        if (!defaultUpdateEquationsIMP) {
            defaultUpdateEquationsIMP = (void (*)(id, SEL, FEMModel*, FEMSolution*, Element_t*, double**, double*, int*, int*, FEMMatrixCRS*, FEMMatrixBand*))
            [core methodForSelector: @selector(defaultUpdateEquations:inSolution:forElement:realStiff:realForce:stiffRows:stiffCols:crsMatrix:bandMatrix:)];
        }
        if (!navierStokesComposeMassMatrixIMP) {
            navierStokesComposeMassMatrixIMP = (void (*)(id, SEL, double**, double**, double*, double**, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, BOOL, NSString*, int, BOOL, double*, double*, BOOL, double**, BOOL, double*, double*, BOOL, BOOL, double*, BOOL, BOOL, BOOL, BOOL, Element_t*, int, int, int, Nodes_t*, FEMSolution*, FEMCore*, FEMMesh*, FEMModel*, FEMNumericIntegration*, FEMMaterial*, FEMElementDescription*, FEMCoordinateSystems*, FEMMaterialModels*, FEMDifferentials*, FEMListUtilities*, FEMUtilities*))
            [navierStokes methodForSelector:@selector(navierStokesComposeMassMatrix:stiffMatrix:forceVector:loadVector:nodalViscosity:nodalDensity:velocityX:velocityY:velocityZ:meshVelocityX:meshVelocityY:meshVelocityZ:nodalPressure:nodalTemperature:isConvect:stabilizeFlag:compressibilityModel:isPseudoCompressible:nodalCompressibility:nodalGasConstant:isPorous:nodalDrag:isPotentialForce:potentialField:potentialCoefficient:isMagneticForce:isRotating:omega:isDivDiscretization:isGradPDriscretization:isNewtonLinearization:isTransient:element:numberOfNodes:rows:cols:nodes:solution:core:mesh:model:integration:material:elementDescription:coordinateSystems:materialModels:differentials:listUtilities:utilities:)];
        }
        if (!navierStokesBoundaryIMP) {
            navierStokesBoundaryIMP = (void (*)(id, SEL, double**, double*, double**, double*, double*, double*, double**, BOOL, Element_t*, int, Nodes_t*, FEMMesh*, FEMModel*, FEMNumericIntegration*, FEMElementDescription*, FEMElementUtils*))
            [navierStokes methodForSelector:@selector(navierStokesBoundary:boundaryVector:loadVector:nodalAlpha:nodalBeta:nodalExtPressure:nodalSlipCoefficient:isNormalTangential:element:numberOfNodes:nodes:mesh:model:integration:elementDescription:elementUtils:)];
        }
        if (!navierStokesCylindricalComposeMassMatrixIMP) {
            navierStokesCylindricalComposeMassMatrixIMP = (void (*)(id, SEL, double**, double**, double*, double**, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, BOOL, NSString*, BOOL, BOOL, double*, double*, BOOL, double**, BOOL, double*, double*, BOOL, BOOL, BOOL, BOOL, Element_t*, int, int, int, Nodes_t*, FEMCore*, FEMMesh*, FEMModel*, FEMNumericIntegration*, FEMElementDescription*, FEMCoordinateSystems*, FEMMaterialModels*, FEMDifferentials*, FEMListUtilities*, FEMUtilities*))
            [navierStokesCylindrical methodForSelector:@selector(navierStokesCylindricalComposeMassMatrix:stiffMatrix:forceVector:loadVector:nodalViscosity:nodalDensity:velocityX:velocityY:velocityZ:meshVelocityX:meshVelocityY:meshVelocityZ:nodalPressure:nodalTemperature:isConvect:stabilizeFlag:isCompressible:isPseudoCompressible:nodalCompressibility:nodalGasConstant:isPorous:nodalDrag:isPotentialForce:potentialField:potentialCoefficient:isMagneticForce:isDivDiscretization:isGradPDriscretization:isNewtonLinearization:element:numberOfNodes:rows:cols:nodes:core:mesh:model:integration:elementDescription:coordinateSystems:materialModels:differentials:listUtilities:utilities:)];
        }
        if (!navierStokesCylindricalBoundaryIMP) {
            navierStokesCylindricalBoundaryIMP = (void (*)(id, SEL, double**, double*, double**, double*, double*, double*, double**, BOOL, Element_t*, int, Nodes_t*, FEMMesh*, FEMModel*, FEMNumericIntegration*, FEMElementDescription*, FEMElementUtils*, FEMCoordinateSystems*))
            [navierStokesCylindrical methodForSelector:@selector(navierStokesCylindricalBoundary:boundaryVector:loadVector:nodalAlpha:nodalBeta:nodalExtPressure:nodalSlipCoefficient:isNormalTangential:element:numberOfNodes:nodes:mesh:model:integration:elementDescription:elementUtils:coordinateSystems:)];
        }
        if (!navierStokesGeneralComposeMassMatrixIMP) {
            navierStokesGeneralComposeMassMatrixIMP = (void (*)(id, SEL, double**, double**, double*, double**, double*, double*, double*, double*, double*, double*, double*, double*, BOOL, BOOL, Element_t*, int, Nodes_t*, FEMCore*, FEMMesh*, FEMModel*, FEMNumericIntegration*, FEMElementDescription*, FEMCoordinateSystems*, FEMMaterialModels*, FEMDifferentials*, FEMListUtilities*, FEMUtilities*))
            [navierStokesGeneral methodForSelector:@selector(navierStokesGeneralComposeMassMatrix:stiffMatrix:forceVector:loadVector:nodalViscosity:nodalDensity:velocityX:velocityY:velocityZ:meshVelocityX:meshVelocityY:meshVelocityZ:isStabilize:isNewtonLinearization:element:numberOfNodes:nodes:core:mesh:model:integration:elementDescription:coordinateSystems:materialModels:differentials:listUtilities:utilities:)];
        }
        if (!navierStokesGeneralBoundaryIMP) {
            navierStokesGeneralBoundaryIMP = (void (*)(id, SEL, double**, double*, double**, double*, double*, double*, double**, Element_t*, int, Nodes_t*, FEMMesh*, FEMModel*, FEMNumericIntegration*, FEMElementDescription*, FEMCoordinateSystems*))
            [navierStokesGeneral methodForSelector:@selector(navierStokesGeneralBoundary:boundaryVector:loadVector:nodalAlpha:nodalBeta:nodalExtPressure:nodalSlipCoefficient:element:numberOfNodes:nodes:mesh:model:integration:elementDescription:coordinateSystems:)];
        }
    });
    
    if (solution.matrix == nil) return;
    
    NSLog(@"FEMFlowSolution:solutionComputer: solving the Navier-Stokes equations.\n");
    
    if ((solution.solutionInfo)[@"solver coordinate system"] != nil) {
        localCoords = (solution.solutionInfo)[@"solver coordinate system"];
        useLocalCoords = YES;
    }
    if (useLocalCoords == YES) {
        modelCoords = coordinatesSystems.coordinates;
        modelDim = model.dimension;
        if ([localCoords isEqualToString:@"cartesian 2d"] == YES) {
            coordinatesSystems.coordinates = 1;
            model.dimension = 2;
            NSLog(@"FEMFlowSolution:solutionComputer: Solution coordinate system is cartesian 2D.\n");
        } else if ([localCoords isEqualToString:@"cartesian 3d"] == YES) {
            coordinatesSystems.coordinates = 1;
            model.dimension = 3;
            NSLog(@"FEMFlowSolution:solutionComputer: Solution coordinate system is cartesian 3D.\n");
        } else if ([localCoords isEqualToString:@"axi symmetric"] == YES) {
            coordinatesSystems.coordinates = 4;
            model.dimension = 2;
            NSLog(@"FEMFlowSolution:solutionComputer: Solution coordinate system is axi symmetric.\n");
        } else if ([localCoords isEqualToString:@"cylindric symmetric"] == YES) {
            coordinatesSystems.coordinates = 3;
            model.dimension = 3;
            NSLog(@"FEMFlowSolution:solutionComputer: Solution coordinate system is cylindric symmetric.\n");
        } else {
            NSLog(@"FEMFlowSolution:solutionComputer: Solution coordinate system not recognized, using original.\n");
        }
    }
    
    // Check for flow model. one of 'full', 'no convection', 'stokes'
    ifTransient = transient;
    convect = YES;
    if ((solution.solutionInfo)[@"flow model"] != nil) {
        flowModel = (solution.solutionInfo)[@"flow model"];
    }
    
    if ([flowModel isEqualToString:@"no convection"] == YES) {
        convect = NO;
    } else if ([flowModel isEqualToString:@"stokes"] == YES) {
        convect = NO;
        ifTransient = NO;
    } else {
        flowModel = @"full";
    }
    
    if (ifTransient == YES) {
        timeIntegration = [[FEMTimeIntegration alloc] init];
    }
    
    mesh = (FEMMesh *)model.mesh;
    dim = model.dimension;
    
    flowContainers = solution.variable.getContainers;
    _flowPerm = flowContainers->Perm;
    _flowSolution = flowContainers->Values;
    _nsdofs = solution.variable.dofs;
    varName = [solution.variable canonicalizeName];
    
    _localNodes = 0;
    for (i=0; i<flowContainers->sizePerm; i++) {
        if (_flowPerm[i] >= 0) _localNodes++;
    }
    if (_localNodes <= 0) return;
    
    tempSol = [utilities getVariableFrom:solution.mesh.variables model:model name:@"temperature" onlySearch:NULL maskName:nil info:&found];
    if (tempSol != nil) {
        tempSolContainers = tempSol.getContainers;
        tempPerm = tempSolContainers->Perm;
        temperature = tempSolContainers->Values;
        if (ifTransient == YES) {
            tempPrev = doublevec(0, tempSolContainers->size1PrevValues-1);
            for (i=0; i<tempSolContainers->size1PrevValues; i++) {
                tempPrev[i] = tempSolContainers->PrevValues[i][0];
            }
        }
    }
    
    meshVeloSol = [utilities getVariableFrom:solution.mesh.variables model:model name:@"mesh velocity" onlySearch:NULL maskName:nil info:&found];
    if (meshVeloSol != nil) {
        meshVeloSolContainers = meshVeloSol.getContainers;
        meshVeloPerm = meshVeloSolContainers->Perm;
        meshVelocity = meshVeloSolContainers->Values;
    }
    
    densitySol = [utilities getVariableFrom:solution.mesh.variables model:model name:@"density" onlySearch:NULL maskName:nil info:&found];
    if (densitySol != nil) densityContainers = densitySol.getContainers;
    
    matContainers = solution.matrix.getContainers;
    uNorm = solution.variable.norm;
    
    // Allocate some permanent storage, this is done first time only
    if (_allocationDone == NO || solution.mesh.changed == YES) {
        n = solution.mesh.maxElementDofs;
        if (_allocationDone == YES) {
            free_dvector(_u, 0, n-1);
            free_dvector(_v, 0, n-1);
            free_dvector(_w, 0, n-1);
            free_dvector(_mu, 0, n-1);
            free_dvector(_mv, 0, n-1);
            free_dvector(_mw, 0, n-1);
            free_ivector(_indexes, 0, n-1);
            free_dvector(_pressure, 0, n-1);
            free_dvector(_prevPressure, 0, n-1);
            free_dvector(_pseudoCompressibility, 0, n-1);
            free_dvector(_prevDensity, 0, n-1);
            free_dvector(_density, 0, n-1);
            free_dvector(_layerThickness, 0, n-1);
            free_dvector(_surfaceRoughness, 0, n-1);
            free_dvector(_permeability, 0, n-1);
            free_dvector(_mx, 0, n-1);
            free_dvector(_my, 0, n-1);
            free_dvector(_mz, 0, n-1);
            free_dmatrix(_slipCoeff, 0, 2, 0, n-1);
            free_dmatrix(_drag, 0, 2, 0, n-1);
            free_dvector(_timeForce, 0, (2*_nsdofs*n)-1);
            free_dvector(_force, 0, (2*_nsdofs*n)-1);
            free_dvector(_viscosity, 0, n-1);
            free_dmatrix(_mass, 0, (2*_nsdofs*n)-1, 0, (2*_nsdofs*n)-1);
            free_dmatrix(_stiff, 0, (2*_nsdofs*n)-1, 0, (2*_nsdofs*n)-1);
            free_dvector(_heatExpansionCoeff, 0, n-1);
            free_dvector(_gasConstant, 0, n-1);
            free_dvector(_heatCapacity, 0, n-1);
            free_dvector(_referenceTemperature, 0, n-1);
            free_dvector(_localTempPrev, 0, n-1);
            free_dvector(_localTemperature, 0, n-1);
            free_dvector(_pSolution, 0, flowContainers->sizeValues-1);
            free_dvector(_potentialField, 0, n-1);
            free_dvector(_potentialCoefficient, 0, n-1);
            free_dmatrix(_loadVector, 0, 3, 0, n-1);
            free_dvector(_alpha, 0, n-1);
            free_dvector(_beta, 0, n-1);
            free_dvector(_extPressure, 0, n-1);
            free_dvector(_elementNodes->x, 0, n-1);
            free_dvector(_elementNodes->y, 0, n-1);
            free_dvector(_elementNodes->z, 0, n-1);
            free(_elementNodes);
            [self FEMFlowSolution_nullify];
        }
        _u = doublevec(0, n-1);
        _v = doublevec(0, n-1);
        _w = doublevec(0, n-1);
        _mu = doublevec(0, n-1);
        _mv = doublevec(0, n-1);
        _mw = doublevec(0, n-1);
        _indexes = intvec(0, n-1);
        _pressure = doublevec(0, n-1);
        _prevPressure = doublevec(0, n-1);
        _pseudoCompressibility = doublevec(0, n-1);
        _prevDensity = doublevec(0, n-1);
        _density = doublevec(0, n-1);
        _layerThickness = doublevec(0, n-1);
        _surfaceRoughness = doublevec(0, n-1);
        _permeability = doublevec(0, n-1);
        _mx = doublevec(0, n-1);
        _my = doublevec(0, n-1);
        _mz = doublevec(0, n-1);
        _slipCoeff = doublematrix(0, 2, 0, n-1);
        _drag = doublematrix(0, 2, 0, n-1);
        _timeForce = doublevec(0, (2*_nsdofs*n)-1);
        _force = doublevec(0, (2*_nsdofs*n)-1);
        _viscosity = doublevec(0, n-1);
        _mass = doublematrix(0, (2*_nsdofs*n)-1, 0, (2*_nsdofs*n)-1);
        _stiff = doublematrix(0, (2*_nsdofs*n)-1, 0, (2*_nsdofs*n)-1);
        _heatExpansionCoeff = doublevec(0, n-1);
        _gasConstant = doublevec(0, n-1);
        _heatCapacity = doublevec(0, n-1);
        _referenceTemperature = doublevec(0, n-1);
        _localTempPrev = doublevec(0, n-1);
        _localTemperature = doublevec(0, n-1);
        _pSolution = doublevec(0, flowContainers->sizeValues-1);
        _potentialField = doublevec(0, n-1);
        _potentialCoefficient = doublevec(0, n-1);
        _loadVector = doublematrix(0, 3, 0, n-1);
        _alpha = doublevec(0, n-1);
        _beta = doublevec(0, n-1);
        _extPressure = doublevec(0, n-1);
        _elementNodes = (Nodes_t*)malloc(sizeof(Nodes_t));
        initNodes(_elementNodes);
        _elementNodes->x = doublevec(0, n-1);
        _elementNodes->y = doublevec(0, n-1);
        _elementNodes->z = doublevec(0, n-1);
        _elementNodes->numberOfNodes = n;
        if (_u == NULL || _v == NULL || _w == NULL || _mu == NULL || _mv == NULL || _mw == NULL || _indexes == NULL || _pressure == NULL ||
            _prevPressure == NULL || _pseudoCompressibility == NULL || _prevDensity == NULL || _density == NULL || _layerThickness == NULL ||
            _surfaceRoughness == NULL || _permeability == NULL || _mx == NULL || _my == NULL || _mz == NULL || _slipCoeff == NULL ||
            _drag == NULL || _timeForce == NULL || _force == NULL || _viscosity == NULL || _mass == NULL || _stiff == NULL ||
            _heatExpansionCoeff == NULL || _gasConstant == NULL || _heatCapacity == NULL || _referenceTemperature == NULL || _localTempPrev == NULL ||
            _localTemperature == NULL || _pSolution == NULL || _potentialField == NULL || _potentialCoefficient == NULL || _loadVector == NULL ||
            _alpha == NULL || _beta == NULL || _extPressure == NULL || _elementNodes->x == NULL || _elementNodes->y == NULL || _elementNodes->z == NULL) {
            fatal("FEMFlowSolution:solutionComputer", "Memory allocation error.");
        }
        
        memset(*_drag, 0.0, (3*n)*sizeof(double) );
        
        for (FEMMaterial *material in model.materials) {
            compressibilityFlag = [listUtilities listGetString:model inArray:material.valuesList forVariable:@"compressibility model" info:&found];
            if (found == YES && [compressibilityFlag isEqualToString:@"artificial compressible"] == YES) pseudoPressureExists = YES;
        }
        
        if (pseudoPressureExists == YES) {
            if (_allocationDone == YES) {
                free_dvector(_pseudoPressure, 0, _sizePseudoPressure-1);
            }
            _sizePseudoPressure = flowContainers->sizeValues / _nsdofs;
            _pseudoPressure = doublevec(0, _sizePseudoPressure-1);
        }
        
        _rows = 2*_nsdofs*n;
        _cols = 2*_nsdofs*n;
        
        _allocationDone = YES;
    }
    
    timeVar = [utilities getVariableFrom:solution.mesh.variables model:model name:@"time step" onlySearch:NULL maskName:nil info:&found];
    timeVarContainers = timeVar.getContainers;
    dt = round(timeVarContainers->Values[0]);
    if (saveTimeStep != dt) {
        if (_pDensity0 != NULL) memcpy(_pDensity0, _pDensity1, _sizePDensity0*sizeof(double));
        saveTimeStep = dt;
    }
    
    vector = (listBuffer*)malloc(sizeof(listBuffer));
    matrix = (listBuffer*)malloc(sizeof(listBuffer));
    pwrk = (listBuffer*)malloc(sizeof(listBuffer));
    *vector = (listBuffer){.ivector=NULL, .vector=NULL, .matrix=NULL, .m=0, .n=0, .p=0};
    *matrix = (listBuffer){.ivector=NULL, .vector=NULL, .matrix=NULL, .m=0, .n=0, .p=0};
    *pwrk = (listBuffer){.ivector=NULL, .vector=NULL, .matrix=NULL, .m=0, .n=0, .p=0};
    
    // Additional initialization
    found = [listUtilities listGetConstRealArray:model inArray:model.constants.valuesList forVariable:@"gravity" buffer:matrix];
    if (found == YES) {
        for (i=0; i<3; i++) {
            gravity[i] = matrix->matrix[i][0]*matrix->matrix[3][0];
        }
    } else {
        memset(gravity, 0.0, sizeof(gravity) );
        gravity[1] = -9.81;
    }
    
    bubbles = [(solution.solutionInfo)[@"bubbles"] boolValue];
    stabilize = [(solution.solutionInfo[@"stabilize"]) boolValue];
    
    if ((solution.solutionInfo)[@"stabilization method"] != nil) {
        stabilizeFlag = (solution.solutionInfo)[@"stabilization method"];
    } else {
        if (stabilize == YES) {
            stabilizeFlag = @"stabilized";
        } else if (bubbles == YES) {
            stabilizeFlag = @"bubbles";
        } else {
            stabilizeFlag = @"stabilized";
        }
    }
    
    if ([stabilizeFlag isEqualToString:@"bubbles"] == YES) bubbles = YES;
    
    divDiscretization = [(solution.solutionInfo)[@"div dicretization"] boolValue];
    gradPDiscretization = [(solution.solutionInfo)[@"gradp discretization"] boolValue];
    nonLinearTol = [(solution.solutionInfo)[@"nonlinear system convergence tolerance"] doubleValue];
    if (nonLinearTol < 0.0) nonLinearTol = 0.0;
    newtonTol = [(solution.solutionInfo)[@"nonlinear system newton after tolerance"] boolValue];
    if (newtonTol < 0.0) newtonTol = 0.0;
    
    newtonIter = [(solution.solutionInfo)[@"nonlinear system newton after iterations"] intValue];
    if (newtonIter == 0) newtonLinearization = YES;
    
    if ([(solution.solutionInfo)[@"nonlinear system reset newton"] boolValue] == YES) newtonLinearization = NO;
    
    nonLinearIter = [(solution.solutionInfo)[@"nonlinear system maximum iterations"] intValue];
    if (nonLinearIter < 0) nonLinearIter = 0;
    
    if ((solution.solutionInfo)[@"nonlinear system norm dofs"] == nil) {
        [solution.solutionInfo setObject:@(_nsdofs-1) forKey:@"nonlinear system norm dofs"];
    }
    
    if ((solution.solutionInfo)[@"free surface after tolerance"] != nil) {
        freeSTol = [(solution.solutionInfo)[@"free surface after tolerance"] doubleValue];
    } else {
        freeSTol = DBL_MAX;
    }
    
    if ((solution.solutionInfo)[@"free surface after iterations"] != nil) {
        freeSIter = [(solution.solutionInfo)[@"free surface after iterations"] intValue];
    } else {
        freeSIter = 0;
    }
    
    // We do our own relaxation
    if ((solution.solutionInfo)[@"nonlinear system relaxation factor"] != nil) {
        nonLinearRelax = [solution.solutionInfo[@"nonlinear system relaxation factor"] doubleValue];
    } else {
        nonLinearRelax = 1.0;
    }
    [solution.solutionInfo setObject:@1 forKey:@"nonlinear system relaxation factor"];
    
    if (nonLinearRelax != 1.0) {
        [solution.solutionInfo setObject:@YES forKey:@"skip compute nonlinear change"];
    }
    
    // Check if free surfaces present
    freeSurfaceFlag = NO;
    for (FEMBoundaryCondition *boundaryCondition in model.boundaryConditions) {
        freeSurfaceFlag = (freeSurfaceFlag == YES || [listUtilities listGetLogical:model inArray:boundaryCondition.valuesList forVariable:@"free surface" info:&found] == YES) ? YES : NO;
        if (freeSurfaceFlag == YES) break;
    }
    
    [self FEMFlowSolution_checkCircleBoundaryModel:model];
    
    totat = 0.0;
    totst = 0.0;
    
    // Initialize the pressure to be used in artificial compressibility
    if (pseudoPressureExists == YES) {
        for (i=_nsdofs-1; i<flowContainers->sizeValues; i+=_nsdofs) {
            _pseudoPressure[i] = flowContainers->Values[i];
        }
        vDSP_sveD(_pseudoPressure, 1, &sum, _sizePseudoPressure);
        NSLog(@"FEMFlowSolution:solutionComputer: pseudoPressure mean: %f.\n", sum/_sizePseudoPressure);
        
        pseudoCompressibilityScale = [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:@"artificial compressibility scaling" info:&found minValue:NULL maxValue:NULL];
        if (found == NO) pseudoCompressibilityScale = 1.0;
        
        if (ifTransient == YES) pseudoCompressibilityScale = pseudoCompressibilityScale / timeStep;
        
        pseudoPressureUpdate = [listUtilities listGetLogical:model inArray:model.simulation.valuesList forVariable:@"pseudo pressure update" info:&found];
        if (found == NO) pseudoPressureUpdate = NO;
    }
    
    // Check if we do the assembly in parallel on the GPU
    BOOL parallelAssembly = NO;
    if (solution.solutionInfo[@"parallel assembly"] != nil) {
        parallelAssembly = [solution.solutionInfo[@"parallel assembly"] boolValue];
    }
    
    FEMNumericIntegration *integration;
    if (parallelAssembly == NO) {
        integration = [[FEMNumericIntegration alloc] init];
        if ([integration allocation:mesh] == NO) fatal("FEMFlowSolution:solutionComputer", "Allocation error in FEMNumericIntegration.");
    }
    
    // If we do the matrix assembly on the GPU, prepare everything here
    if (parallelAssembly == YES) {
        if (_initializeGPU == NO) {
            // Get the GPU device
            _device = find_single_device();
            device_stats(_device);
            
            // Create the context of the command queue
            _context = clCreateContext(0, 1, &_device, NULL, NULL, &_err);
            _cmd_queue = clCreateCommandQueue(_context, _device, 0, NULL);
            
            NSString *kernelfile;
            if (solution.solutionInfo[@"gpu kernel source file"] != nil) {
                NSString *source = solution.solutionInfo[@"gpu kernel source file"];
                if ([source containsString:@".cl"] == YES) {
                    kernelfile = [NSString stringWithString:source];
                } else {
                    kernelfile = [source stringByAppendingPathExtension:@"cl"];
                }
            } else {
                fatal("FEMFlowSolution:solutionComputer", "Parallel assembly requires a kernel source file.");
            }
            
            int success = LoadFileIntoString((char *)[kernelfile UTF8String], &_kernel_source, &_src_len);
            if (success < 0) {
                fatal("FEMFlowSolution:solutionComputer", "Can't load kernel source.");
            }
            
            // Allocate memory for program and kernels
            // Create the program .cl file
            _program = clCreateProgramWithSource(_context, 1, (const char**)&_kernel_source, NULL, &_err);
            if (_err) {
                fatal("FEMFlowSolution:solutionComputer", "Can't create program. Error number: ", _err);
            }
            
            // Build the program (compile it)
            const char options[] = "-DKERNEL_FP_32 -DMESH_DIMENSION_3 -DELEMENT_DIMENSION_3 -DMODEL_DIMENSION_3";
            _err = clBuildProgram(_program, 1, &_device, options, NULL, &_err);
            if (_err < 0) {
                size_t log_size;
                clGetProgramBuildInfo(_program, _device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
                char *program_log = (char *)malloc(log_size+1);
                program_log[log_size] = '\0';
                clGetProgramBuildInfo(_program, _device, CL_PROGRAM_BUILD_LOG, log_size+1, program_log, NULL);
                NSLog(@"FEMFlowSolution:solutionComputer: %s\n", program_log);
                free(program_log);
                fatal("FEMFlowSolution:solutionComputer", "Saino will abort the simulation now...");
            }
            
            
            
            _initializeGPU = YES;
        }
    }
    
    for (iter=1; iter<=nonLinearIter; iter++) {
        
        if (pseudoPressureExists == YES && pseudoPressureUpdate == YES) {
            for (i=_nsdofs-1; i<flowContainers->sizeValues; i+=_nsdofs) {
                _pseudoPressure[i] = flowContainers->Values[i];
            }
        }
        
        at = cputime();
        at0 = realtime();
        at1 = realtime();
        
        NSLog(@"FEMFlowSolution:solutionComputer:\n");
        NSLog(@"FEMFlowSolution:solutionComputer:\n");
        NSLog(@"FEMFlowSolution:solutionComputer: -----------------------------------------------------------\n");
        NSLog(@"FEMFlowSolution:solutionComputer: NAVIER-STOKES ITERATION %d.\n", iter);
        NSLog(@"FEMFlowSolution:solutionComputer: -----------------------------------------------------------\n");
        NSLog(@"FEMFlowSolution:solutionComputer:\n");
        NSLog(@"FEMFlowSolution:solutionComputer: Starting Assembly...\n");

        [core initializeToZeroMatrix:solution.matrix forceVector: matContainers->RHS sizeForceVector:matContainers->sizeRHS model:model solution:solution];
        
        bf_id = -1;
        body_id = -1;
        
        startAdvanceOutput((char *)[@"FEMFlowSolution" UTF8String], (char *)[@"Assembly:" UTF8String]);
        
        // Bulk elements
        if (parallelAssembly == YES) {
            
            GaussIntegrationPoints *IP = NULL;
            int elementCode, gaussPoints;
            
            // Figure out the element code and the number of Gauss points for the active elements
            element = getActiveElementIMP(core, @selector(getActiveElement:solution:model:), 0, solution, model);
            elementCode = element->Type.ElementCode;
            gaussPoints = element->Type.GaussPoints;
            IP = GaussQuadratureGPU(elementCode, gaussPoints);
            
            // Get the element basis functions
            
        } else {
            [self FEMFlowSolution_doAssembly:solution
                                        core:core
                                       model:model
                               listUtilities:listUtilities
                          coordinatesSystems:coordinatesSystems
                                navierStokes:navierStokes
                     navierStokesCylindrical:navierStokesCylindrical
                         navierStokesGeneral:navierStokesGeneral
                                 integration:integration
                               differentials:differentials
                              materialModels:materialModels
                                   utilities:utilities
                          elementDescription:elementDescription
                             timeIntegration:timeIntegration
                                   crsMatrix:crsMatrix
                                  bandMatrix:bandMatrix
                                        bfID:bf_id
                                      bodyID:body_id
                                        eqID:eq_id
                                       matID:mat_id
                                      vector:vector
                                      matrix:matrix
                                        pwrk:pwrk
                                 meshVeloSol:meshVeloSol
                                meshVeloPerm:meshVeloPerm
                                meshVelocity:meshVelocity
                                     tempSol:tempSol
                                    tempPerm:tempPerm
                                 temperature:temperature
                                    tempPrev:tempPrev
                                 ifTransient:ifTransient
                                  convection:convect
                              flowContainers:flowContainers
                                  densitySol:densitySol
                           densityContainers:densityContainers
                                          dt:dt
                                     gravity:gravity
                                    timeStep:timeStep
                         newtonLinearization:newtonLinearization
                                   transient:transient
                               stabilizeFlag:stabilizeFlag
                           divDiscretization:divDiscretization
                         gradPDiscretization:gradPDiscretization
                                   stabilize:stabilize
                                     bubbles:bubbles
                            getActiveElement:getActiveElementIMP
                     getEquationIDForElement:getEquationIDForElementIMP
                    getBodyForceIDForElement:getBodyForceIDForElementIMP
                     getMaterialIDForElement:getMaterialIDForElementIMP
                               listGetString:listGetStringIMP
                              listGetLogical:listGetLogicalIMP
                       listGetConstRealArray:listGetConstRealArrayIMP
                getNumberOfBubbleDofsElement:getNumberOfBubbleDofsElementIMP
                      getElementDofsSolution:getElementDofsSolutionIMP
                                    getNodes:getNodesIMP
                                     getReal:getRealIMP
                            listGetConstReal:listGetConstRealIMP
                            listGetRealArray:listGetRealArrayIMP
               navierStokesComposeMassMatrix:navierStokesComposeMassMatrixIMP
    navierStokesCylindricalComposeMassMatrix:navierStokesCylindricalComposeMassMatrixIMP
        navierStokesGeneralComposeMassMatrix:navierStokesGeneralComposeMassMatrixIMP
                       defaultFirstOrderTime:defaultFirstOrderTimeIMP
                           nsCondensateStiff:nsCondensateStiffIMP
                      defaultUpdateEquations:defaultUpdateEquationsIMP];
        }
        
        [core defaultFinishBulkAssemblySolution:solution bulkUpdate:NULL];
        NSLog(@"FEMFlowSolution:solutionComputer: Assembly done.\n");
        
        // Newmann and Newton boundary conditions
        NSString *normalTangentialName = [@"normal-tangential " stringByAppendingString:[solution.variable canonicalizeName]];
        for (t=0; t<solution.mesh.numberOfBoundaryElements; t++) {
            element = [core getBoundaryElement:solution atIndex:t];
            if ([core isActiveBoundaryElement:element inSolution:solution model:model] == NO) continue;
            
            n = element->Type.NumberOfNodes;
            
            // The element type 101 (point element) can only be used
            // to set Dirichlet BCs, so skip them at this stage
            if ([core isFluxElement:element mesh:mesh] == NO) continue;
            
            [core getNodes:solution model:model inElement:element resultNodes:_elementNodes numberOfNodes:NULL mesh:nil];
            
            bc = [core getBoundaryCondition:model forElement:element];
            if (bc == nil) continue;
            
            gotForceBC = listGetLogicalIMP(listUtilities, @selector(listGetLogical:inArray:forVariable:info:), model, bc, @"flow force bc", &found);
            if (found == NO) gotForceBC = YES;
            
            if (gotForceBC == YES) {
                memset(*_loadVector, 0.0, (4*solution.mesh.maxElementDofs)*sizeof(double) );
                memset(_alpha, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
                memset(_extPressure, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
                memset(_beta, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
                memset(*_slipCoeff, 0.0, (3*solution.mesh.maxElementDofs)*sizeof(double) );
                memset(*_stiff, 0.0, ((2*_nsdofs*solution.mesh.maxElementDofs)*(2*_nsdofs*solution.mesh.maxElementDofs))*sizeof(double) );
                memset(_force, 0.0, (2*_nsdofs*solution.mesh.maxElementDofs)*sizeof(double) );
                
                // (at thr moment, the following is done...)
                // BC: \tau \cdot n = \alpha n + @\beta/@\t + R_k u_k + F
                
                // normal force BC: \tau \cdot n = \alpha n
                if (listGetLogicalIMP(listUtilities, @selector(listGetLogical:inArray:forVariable:info:), model, bc, @"free surface", &found) == YES) {
                    found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, bc, @"surface tension coefficient", vector, listUtilities);
                    if (found == YES) memcpy(_alpha, vector->vector, n*sizeof(double));
                }
                
                found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, bc, @"external pressure", vector, listUtilities);
                if (found == YES) {
                    memcpy(_extPressure, vector->vector, n*sizeof(double));
                } else {
                    found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, bc, @"normal pressure", vector, listUtilities);
                    if (found == YES) memcpy(_extPressure, vector->vector, n*sizeof(double));
                }
                
                // Tangential force BC:
                // \tau \cdot n = @\beta/@t (tangential derivative of something)
                if (tempSol != nil) {
                    found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, bc, @"surface tension expansion coefficient", vector, listUtilities);
                    if (found == YES) {
                         memcpy(_beta, vector->vector, n*sizeof(double));
                        for (j=0; j<n; j++) {
                            k = tempPerm[element->NodeIndexes[j]];
                            if (k >= 0) _beta[j] = 1.0 - _beta[j] * temperature[k];
                        }
                        found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, bc, @"surface tension coefficient", vector, listUtilities);
                        if (found == YES) {
                            for (j=0; j<n; j++) {
                                _beta[j] = _beta[j] * vector->vector[j];
                            }
                        }
                    } else {
                        found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, bc, @"surface tension coefficient", vector, listUtilities);
                        if (found == YES) memcpy(_beta, vector->vector, n*sizeof(double));
                    }
                }
                
                // Force in given direction BC: \tau \cdot n = F
                found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, bc, @"pressure 1", vector, listUtilities);
                if (found == YES) {
                    for (i=0; i<n; i++) {
                        _loadVector[0][i] = vector->vector[i];
                    }
                }
                found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, bc, @"pressure 2", vector, listUtilities);
                if (found == YES) {
                    for (i=0; i<n; i++) {
                        _loadVector[1][i] = vector->vector[i];
                    }
                }
                found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, bc, @"pressure 3", vector, listUtilities);
                if (found == YES) {
                    for (i=0; i<n; i++) {
                        _loadVector[2][i] = vector->vector[i];
                    }
                }
                found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, bc, @"mass flux", vector, listUtilities);
                if (found == YES) {
                    for (i=0; i<n; i++) {
                        _loadVector[3][i] = vector->vector[i];
                    }
                }
                
                // Slip boundary condition BC: \tau \cdot n = R_k u_k
                found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, bc, @"slip coefficient 1", vector, listUtilities);
                if (found == YES) {
                    for (i=0; i<n; i++) {
                        _slipCoeff[0][i] = vector->vector[i];
                    }
                }
                found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, bc, @"slip coefficient 2", vector, listUtilities);
                if (found == YES) {
                    for (i=0; i<n; i++) {
                        _slipCoeff[1][i] = vector->vector[i];
                    }
                }
                found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, bc, @"slip coefficient 3", vector, listUtilities);
                if (found == YES) {
                    for (i=0; i<n; i++) {
                        _slipCoeff[2][i] = vector->vector[i];
                    }
                }
                
                normalTangential = listGetLogicalIMP(listUtilities, @selector(listGetLogical:inArray:forVariable:info:), model, bc, @"normal-tangential velocity", &found);
                if (found == NO) {
                    normalTangential = listGetLogicalIMP(listUtilities, @selector(listGetLogical:inArray:forVariable:info:), model, bc, normalTangentialName, &found);
                }
                
                switch (coordinatesSystems.coordinates) {
                    case cartesian:
                        navierStokesBoundaryIMP(navierStokes, @selector(navierStokesBoundary:boundaryVector:loadVector:nodalAlpha:nodalBeta:nodalExtPressure:nodalSlipCoefficient:isNormalTangential:element:numberOfNodes:nodes:mesh:model:integration:elementDescription:elementUtils:), _stiff, _force, _loadVector, _alpha, _beta, _extPressure, _slipCoeff, normalTangential, element, n, _elementNodes, solution.mesh, model, integration, elementDescription, elementUtils);
                        break;
                        
                    case cylindric:
                    case cylindric_symmetric:
                    case axis_symmetric:
                        navierStokesCylindricalBoundaryIMP(navierStokesCylindrical, @selector(navierStokesCylindricalBoundary:boundaryVector:loadVector:nodalAlpha:nodalBeta:nodalExtPressure:nodalSlipCoefficient:isNormalTangential:element:numberOfNodes:nodes:mesh:model:integration:elementDescription:elementUtils:coordinateSystems:), _stiff, _force, _loadVector, _alpha, _beta, _extPressure, _slipCoeff, normalTangential, element, n, _elementNodes, solution.mesh, model, integration, elementDescription, elementUtils, coordinatesSystems);
                        break;
                        
                    default:
                        navierStokesGeneralBoundaryIMP(navierStokesGeneral, @selector(navierStokesGeneralBoundary:boundaryVector:loadVector:nodalAlpha:nodalBeta:nodalExtPressure:nodalSlipCoefficient:element:numberOfNodes:nodes:mesh:model:integration:elementDescription:coordinateSystems:), _stiff, _force, _loadVector, _alpha, _beta, _extPressure, _slipCoeff, element, n, _elementNodes, solution.mesh, model, integration, elementDescription, coordinatesSystems);
                        break;
                }
                
                if (listGetLogicalIMP(listUtilities, @selector(listGetLogical:inArray:forVariable:info:), model, bc, @"wall law", &found) == YES ||
                    listGetLogicalIMP(listUtilities, @selector(listGetLogical:inArray:forVariable:info:), model, bc, @"vms wall", &found) == YES) {
                    
                    found = [core getParentMaterialProperty:@"density" forElement:element parentElement:parent model:model listUtilities:listUtilities buffer:vector];
                    if (found == YES) memcpy(_density, vector->vector, n*sizeof(double));
                    found = [core getParentMaterialProperty:@"viscosity" forElement:element parentElement:parent model:model listUtilities:listUtilities buffer:vector];
                    if (found == YES) memcpy(_viscosity, vector->vector, n*sizeof(double));
                    
                    found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, bc, @"boundary layer thickness", vector, listUtilities);
                    if (found == YES) memcpy(_layerThickness, vector->vector, n*sizeof(double));
                    found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, bc, @"surface roughness", vector, listUtilities);
                    if (found == YES) memcpy(_surfaceRoughness, vector->vector, n*sizeof(double));
                    
                    switch (_nsdofs) {
                        case 3:
                            for (i=0; i<n; i++) {
                                _u[i] = _flowSolution[_nsdofs*_flowPerm[element->NodeIndexes[i]]];
                                _v[i] = _flowSolution[_nsdofs*_flowPerm[element->NodeIndexes[i]]+1];
                                _w[i] = 0.0;
                            }
                            break;
                        case 4:
                            for (i=0; i<n; i++) {
                                _u[i] = _flowSolution[_nsdofs*_flowPerm[element->NodeIndexes[i]]];
                                _v[i] = _flowSolution[_nsdofs*_flowPerm[element->NodeIndexes[i]]+1];
                                _w[i] = _flowSolution[_nsdofs*_flowPerm[element->NodeIndexes[i]]+2];
                            }
                            break;
                    }

                }
                if (listGetLogicalIMP(listUtilities, @selector(listGetLogical:inArray:forVariable:info:), model, bc, @"wall law", &found) == YES) {
                    [navierStokes navierStokesWallLawBoundary:_stiff boundaryVector:_force layerThickness:_layerThickness surfaceRoughness:_surfaceRoughness nodalViscosity:_viscosity nodalDensity:_density velocityX:_u velocityY:_v velocityZ:_w element:element numberOfNodes:n nodes:_elementNodes mesh:solution.mesh model:model integration:integration elementDescription:elementDescription elementUtils:elementUtils];
                } else if (listGetLogicalIMP(listUtilities, @selector(listGetLogical:inArray:forVariable:info:), model, bc, @"vms wall", &found) == YES) {
                    [navierStokes vmsWallsBoundary:_stiff boundaryVector:_force layerThickness:_layerThickness surfaceRoughness:_surfaceRoughness nodalViscosity:_viscosity nodalDensity:_density velocityX:_u velocityY:_v velocityZ:_w element:element numberOfNodes:n nodes:_elementNodes solution:solution mesh:solution.mesh model:model integration:integration elementDescription:elementDescription elementUtils:elementUtils];
                }
                
                if (ifTransient == YES) {
                    memset(*_mass, 0.0, ((2*_nsdofs*solution.mesh.maxElementDofs)*(2*_nsdofs*solution.mesh.maxElementDofs))*sizeof(double) );
                    defaultFirstOrderTimeIMP(core, @selector(defaultFirstOrderTime:inSolution:forElement:realMass:realStiff:realForce:stiffRows:stiffCols:timeIntegration:utilities:), model, solution, element, _mass, _stiff, _force, &_rows, &_cols, timeIntegration, utilities);
                }
                
                // Add local stiffness matrix and force vector to global matrix and vector
                defaultUpdateEquationsIMP(core, @selector(defaultUpdateEquations:inSolution:forElement:realStiff:realForce:stiffRows:stiffCols:crsMatrix:bandMatrix:), model, solution, element, _stiff, _force, &_rows, &_cols, crsMatrix, bandMatrix);
            }
        }
        
        // Implement no-slip wall BC code
        for (FEMBoundaryCondition *boundaryCondition in model.boundaryConditions) {
            if (listGetLogicalIMP(listUtilities, @selector(listGetLogical:inArray:forVariable:info:), model, boundaryCondition.valuesList, @"no-slip wall bc", &found) == YES) {
                if ([varName isEqualToString:@"flow solution"] == YES) {
                    double value = 0.0;
                    [listUtilities addConstRealInClassList:boundaryCondition.valuesList theVariable:@"velocity 1" withValue:&value orUsingBlock:nil string:nil];
                    [listUtilities addConstRealInClassList:boundaryCondition.valuesList theVariable:@"velocity 2" withValue:&value orUsingBlock:nil string:nil];
                    if (_nsdofs > 3) [listUtilities addConstRealInClassList:boundaryCondition.valuesList theVariable:@"velocity 3" withValue:&value orUsingBlock:nil string:nil];
                } else {
                    double value = 0.0;
                    for (j=1; j<=_nsdofs-1; j++) {
                        NSString *string = [utilities appendNameFromString:solution.variable.name component:&j];
                        [listUtilities addConstRealInClassList:boundaryCondition.valuesList theVariable:string withValue:&value orUsingBlock:nil string:nil];
                    }
                }
            }
        }
        
        [core defaultFinishAssemblySolution:solution model:model timeIntegration:timeIntegration utilities:utilities];
        
        // Dirichlet boundary conditions
        [core dirichletBoundaryConditions:model inSolution:solution usingOffset:NULL offDiaginalMatrix:NULL];
        NSLog(@"FEMFlowSolution:solutionComputer: Dirichlet conditions done.\n");
        
        // Solve the system and check for convergence
        at = cputime() - at;
        st = cputime();
    
        if (nonLinearRelax != 1.0) memcpy(_pSolution, _flowSolution, flowContainers->sizeValues*sizeof(double));
        uNorm = [core findSolution:solution model:model backRorateNT:NULL];
        
        st = cputime() -  st;
        totat = totat + at;
        totst = totst + st;
        NSLog(@"FEMFlowSolution:solutionComputer: iter: %d, Assembly (s): %f %f.\n", iter, at, totat);
        NSLog(@"FEMFlowSolution:solutionComputer: iter: %d, Solve (s): %f %f.\n", iter, st, totst);
        
        n = _nsdofs * _localNodes;
        
        // This hack is needed because of the fluctuating pressure levels
        if (nonLinearRelax != 1.0) {
            double s = 0.0, ps = 0.0;
            if (compressibilityModel == incompressible) {
                s = _flowSolution[_nsdofs-1];
                ps = _pSolution[_nsdofs-1];
                for (i=_nsdofs-1; i<n; i+=_nsdofs) {
                    _flowSolution[i] = _flowSolution[i] - s;
                    _pSolution[i] = _pSolution[i] - ps;
                }
            }
            
            for (i=0; i<n; i++) {
                _flowSolution[i] = (1.0-nonLinearRelax)*_pSolution[i] + nonLinearRelax * _flowSolution[i];
            }
            
            if (compressibilityModel == incompressible) {
                for (i=_nsdofs-1; i<n; i+=_nsdofs) {
                    _flowSolution[i] = _flowSolution[i] + s;
                }
            }
            
            found = NO;
            relaxBefore = NO;
            if ((solution.solutionInfo)[@"nonlinear system relaxation before"] != nil) {
                relaxBefore = [(solution.solutionInfo)[@"nonlinear system relaxation before"] boolValue];
                found = YES;
            }
            if (found == NO || relaxBefore == YES) {
                [solution.solutionInfo setObject:@NO forKey:@"skip compute nonlinear change"];
                [core computeChange:solution model:model isSteadyState:NO nsize:&n values:_flowSolution values0:_pSolution sizeValues0:&flowContainers->sizeValues];
            }
        }
        
        relativeChange = solution.variable.nonLinChange;
        NSLog(@"FEMFlowSolution:solutionComputer: result norm: %e.\n", solution.variable.norm);
        NSLog(@"FEMFlowSolution:solutionComputer: relative change: %e.\n", relativeChange);
        
        if (relativeChange < newtonTol || iter > newtonIter) newtonLinearization = YES;
        if (relativeChange < nonLinearTol && iter < nonLinearIter) break;
        
        // If free surface in model, this will move the nodal points
        if (freeSurfaceFlag == YES) {
            FEMFreeSurface *freeSurface = [[FEMFreeSurface alloc] init];
            if (relativeChange < freeSTol || iter > freeSIter) computeFree = YES;
            if (computeFree == YES) {
                if ((solution.solutionInfo)[@"free surface relaxation factor"] != nil) {
                    relaxation = [(solution.solutionInfo)[@"free surface relaxation factor"] doubleValue];
                } else {
                    relaxation = 1.0;
                }
                found = NO;
                mbFlag = NO;
                if ((solution.solutionInfo)[@"internal move boundary"] != nil) {
                    mbFlag = [(solution.solutionInfo)[@"internal move boundary"] boolValue];
                    found = YES;
                }
                if (mbFlag == YES || found == NO) [freeSurface moveBoundaryModel:model integration:integration relax:relaxation];
            }
        }
    }
    
    [solution.solutionInfo setObject:@(nonLinearRelax) forKey:@"nonlinear system relaxation factor"];
    
    if ([(solution.solutionInfo)[@"adaptive mesh refinement"] boolValue] == YES) {
        // TODO: implement mesh refinement
    }
    
    [self FEMFlowSolution_checkCircleBoundaryModel:model];
    
    if (useLocalCoords == YES) {
        coordinatesSystems.coordinates = modelCoords;
        model.dimension = modelDim;
    }
    
    if (parallelAssembly == NO) {
        [integration deallocation:mesh];
    }
    
    if (ifTransient == YES) {
        if (tempPrev != NULL) free_dvector(tempPrev, 0, tempSolContainers->size1PrevValues-1);
    }
    
    if (vector->vector != NULL) {
        free_dvector(vector->vector, 0, vector->m-1);
    }
    if (matrix->matrix != NULL) {
        free_dmatrix(matrix->matrix, 0, matrix->m-1, 0, matrix->n-1);
    }
    if (pwrk->tensor != NULL) {
        free_d3tensor(pwrk->tensor, 0, pwrk->m-1, 0, pwrk->n-1, 0, pwrk->p-1);
    }
    free(vector);
    free(matrix);
    free(pwrk);
}

@end
