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
#import "GPUData.h"
#import "GPUUtils.h"

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
                        flowModel:(NSString * __nonnull)flowModel
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
-(void)FEMFlowSolution_initElementBasisFunctions:(FEMCore * __nonnull)core solution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model mesh:(FEMMesh * __nonnull)mesh precisionMode:(NSString * __nonnull)precisionMode getActiveElement:(Element_t* (* __nonnull)(id, SEL, int, FEMSolution*, FEMModel*))getActiveElementIMP;
-(void)FEMFlowSolution_iniNodalData:(FEMCore * __nonnull)core solution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model mesh:(FEMMesh * __nonnull)mesh precisionMode:(NSString * __nonnull)precisionMode getActiveElement:(Element_t* (* __nonnull)(id, SEL, int, FEMSolution*, FEMModel*))getActiveElementIMP getNodes:(void (* __nonnull)(id, SEL, FEMSolution*, FEMModel*, Element_t*, Nodes_t*, int*, FEMMesh*))getNodesIMP getElementDofsSolution:(int (* __nonnull)(id, SEL, FEMSolution*, FEMModel*, Element_t*, int*, BOOL*))getElementDofsSolutionIMP;
-(void)FEMFlowSolution_createMapGPUBuffersMesh:(FEMMesh * __nonnull)mesh solution:(FEMSolution * __nonnull)solution precisionMode:(NSString * __nonnull)precisionMode;
-(void)FEMFlowSolution_setKernelArgumentsCore:(FEMCore * __nonnull)core  model:(FEMModel * __nonnull)model solution:(FEMSolution * __nonnull)solution precisionMode:(NSString * __nonnull)precisionMode getActiveElement:(Element_t* (* __nonnull)(id, SEL, int, FEMSolution*, FEMModel*))getActiveElementIMP;
-(void)FEMFlowSolution_setGPUCore:(FEMCore * __nonnull)core model:(FEMModel * __nonnull)model solution:(FEMSolution * __nonnull)solution mesh:(FEMMesh * __nonnull)mesh precisionMode:(NSString * __nonnull)precisionMode transientSimulation:(BOOL)transient getActiveElement:(Element_t* (* __nonnull)(id, SEL, int, FEMSolution*, FEMModel*))getActiveElementIMP getNodes:(void (* __nonnull)(id, SEL, FEMSolution*, FEMModel*, Element_t*, Nodes_t*, int*, FEMMesh*))getNodesIMP getElementDofsSolution:(int (* __nonnull)(id, SEL, FEMSolution*, FEMModel*, Element_t*, int*, BOOL*))getElementDofsSolutionIMP;
@end

@implementation FEMFlowSolution {
    
    BOOL _allocationDone;
    BOOL _newtonLinearization;
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
    
    // --------------------------------------------------
    // The following definitions are used for
    // the assembly on the GPU
    // --------------------------------------------------
    cl_device_id _device;
    cl_context _context;
    cl_command_queue _cmd_queue;
    cl_program _program;
    cl_kernel _kernel_assembly;
    cl_kernel _kernel_basis_dbasisdx;
    size_t _src_len;
    char * __nullable _kernel_source;
    BOOL _initializeGPU;
    NSString * __nullable _precision;
    cl_uint _positionInColorMapping;
    
    // Ice flow data on the GPU
    ice_flow_gpu *_gpuData;
    
    basis_functions_f *_basis_functions_f;
    basis_functions_d *_basis_functions_d;
    
    element_basis_f *_element_basis_f;
    element_basis_d *_element_basis_d;
    
    element_dBasisdx_f *_element_dBasisdx_f;
    element_dBasisdx_d *_element_dBasisdx_d;
    
    nodal_data_f *_nodal_data_f;
    nodal_data_d *_nodal_data_d;
    
    // GPU buffers
    cl_mem _basis_functions;
    cl_mem _matDiag;
    cl_mem _matRows;
    cl_mem _matCols;
    cl_mem _matValues;
    cl_mem _matRHS;
    cl_mem _colorMapping;
    cl_mem _elementNodeIndexesStore;
    cl_mem _nodesX;
    cl_mem _nodesY;
    cl_mem _nodesZ;
    cl_mem _varSolution;
    cl_mem _varPermutation;
    cl_mem _gpuNewtonLinear;
    cl_mem _element_basis;
    cl_mem _element_dbasisdx;
    cl_mem _nodal_data;
    void * __nullable _mapped_matValues;
    void * __nullable _mapped_matRHS;
    void * __nullable _mapped_varSolution;
    void * __nullable _mapped_gpuNewton;
    void * __nullable _mapped_nodal_data;
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
        fprintf(stdout, "FEMFlowSolution:FEMFlowSolution_checkCircleBoundaryModel: number of elements on circle: %d.\n", l);
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
                        flowModel:(NSString * __nonnull)flowModel
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
    NSString *compressibilityFlag;
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
                                fprintf(stdout, "FEMFlowSolution:solutionComputer: zero or negative density.\n");
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

-(void)FEMFlowSolution_initElementBasisFunctions:(FEMCore * __nonnull)core solution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model mesh:(FEMMesh * __nonnull)mesh precisionMode:(NSString * __nonnull)precisionMode getActiveElement:(Element_t* (* __nonnull)(id, SEL, int, FEMSolution*, FEMModel*))getActiveElementIMP {
    
    Element_t *element = getActiveElementIMP(core, @selector(getActiveElement:solution:model:), 0, solution, model);
    
    if ([precisionMode isEqualToString:@"single"] == YES) {
        
        _basis_functions_f = (basis_functions_f *)malloc(sizeof(basis_functions_f)*8);
        init_basis_functions(_basis_functions_f);
        for(int n=0; n<element->Type.NumberOfNodes; n++) {
            for(int i=0; i<element->Type.BasisFunctions[n].n; i++) {
                _basis_functions_f[n].p[i] = element->Type.BasisFunctions[n].p[i];
                _basis_functions_f[n].q[i] = element->Type.BasisFunctions[n].q[i];
                _basis_functions_f[n].r[i] = element->Type.BasisFunctions[n].r[i];
                _basis_functions_f[n].coeff[i] = (float)element->Type.BasisFunctions[n].coeff[i];
            }
        }
        
    } else if ([precisionMode isEqualToString:@"double"] == YES) {
        
        _basis_functions_d = (basis_functions_d *)malloc(sizeof(basis_functions_d)*8);
        init_basis_functions(_basis_functions_d);
        for(int n=0; n<element->Type.NumberOfNodes; n++) {
            for(int i=0; i<element->Type.BasisFunctions[n].n; i++) {
                _basis_functions_d[n].p[i] = element->Type.BasisFunctions[n].p[i];
                _basis_functions_d[n].q[i] = element->Type.BasisFunctions[n].q[i];
                _basis_functions_d[n].r[i] = element->Type.BasisFunctions[n].r[i];
                _basis_functions_d[n].coeff[i] = element->Type.BasisFunctions[n].coeff[i];
            }
        }
    }
}

-(void)FEMFlowSolution_iniNodalData:(FEMCore * __nonnull)core solution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model mesh:(FEMMesh * __nonnull)mesh precisionMode:(NSString * __nonnull)precisionMode getActiveElement:(Element_t* (* __nonnull)(id, SEL, int, FEMSolution*, FEMModel*))getActiveElementIMP getNodes:(void (* __nonnull)(id, SEL, FEMSolution*, FEMModel*, Element_t*, Nodes_t*, int*, FEMMesh*))getNodesIMP getElementDofsSolution:(int (* __nonnull)(id, SEL, FEMSolution*, FEMModel*, Element_t*, int*, BOOL*))getElementDofsSolutionIMP {
    
    int err, n, nd;
    Element_t *element = NULL;
    const char *type= "double";
    
    int *indexes = intvec(0, solution.mesh.maxElementDofs-1);
    
    Nodes_t *nodes = (Nodes_t*)malloc(sizeof(Nodes_t));
    nodes->x = doublevec(0, mesh.maxElementNodes-1);
    nodes->y = doublevec(0, mesh.maxElementNodes-1);
    nodes->z = doublevec(0, mesh.maxElementNodes-1);
    nodes->numberOfNodes = mesh.maxElementNodes;
    
    if ([precisionMode isEqualToString:@"single"] == YES) {
        _nodal_data_f = (nodal_data_f *)malloc(sizeof(nodal_data_f)*mesh.numberOfBulkElements);
        init_nodal_data(_nodal_data_f, mesh.numberOfBulkElements);
    } else if ([precisionMode isEqualToString:@"double"] == YES) {
        _nodal_data_d = (nodal_data_d *)malloc(sizeof(nodal_data_d)*mesh.numberOfBulkElements);
        init_nodal_data(_nodal_data_d, mesh.numberOfBulkElements);
    }

    for (int t=0; t<mesh.numberOfBulkElements; t++) {
        element = getActiveElementIMP(core, @selector(getActiveElement:solution:model:), t, solution, model);
        getNodesIMP(core, @selector(getNodes:model:inElement:resultNodes:numberOfNodes:mesh:), solution, model, element, nodes, NULL, nil);
        n = element->Type.NumberOfNodes;
        nd = getElementDofsSolutionIMP(core, @selector(getElementDofsSolution:model:forElement:atIndexes:disableDiscontinuousGalerkin:), solution, model, element, indexes, NULL);
        
        if ([precisionMode isEqualToString:@"double"] == YES) {
            err = init_data(precision(), type, 9, nodes->x, _nodal_data_d[t].nodes_x, (size_t)n, nodes->y, _nodal_data_d[t].nodes_y, (size_t)n,
                            nodes->z, _nodal_data_d[t].nodes_z, (size_t)n);
           switch (_nsdofs) {
               case 3:
                   for (int i=0; i<nd; i++) {
                       _nodal_data_d[t].vx[i] = _flowSolution[_nsdofs*_flowPerm[indexes[i]]];
                       _nodal_data_d[t].vy[i] = _flowSolution[_nsdofs*_flowPerm[indexes[i]]+1];
                   }
               case 4:
                   for (int i=0; i<nd; i++) {
                       _nodal_data_d[t].vx[i] = _flowSolution[_nsdofs*_flowPerm[indexes[i]]];
                       _nodal_data_d[t].vy[i] = _flowSolution[_nsdofs*_flowPerm[indexes[i]]+1];
                       _nodal_data_d[t].vz[i] = _flowSolution[_nsdofs*_flowPerm[indexes[i]]+2];
                   }
           }
            for (int i=0; i<nd; i++) {
                _nodal_data_d[t].perm[i] = _flowPerm[indexes[i]];
            }
            
        } else if ([precisionMode isEqualToString:@"single"] == YES) {
            err = init_data(precision(), type, 9, nodes->x, _nodal_data_f[t].nodes_x, (size_t)n, nodes->y, _nodal_data_f[t].nodes_y, (size_t)n,
                            nodes->z, _nodal_data_f[t].nodes_z, (size_t)n);
            switch (_nsdofs) {
                case 3:
                    for (int i=0; i<nd; i++) {
                        _nodal_data_f[t].vx[i] = (float)_flowSolution[_nsdofs*_flowPerm[indexes[i]]];
                        _nodal_data_f[t].vy[i] = (float)_flowSolution[_nsdofs*_flowPerm[indexes[i]]+1];
                    }
                case 4:
                    for (int i=0; i<nd; i++) {
                        _nodal_data_f[t].vx[i] = (float)_flowSolution[_nsdofs*_flowPerm[indexes[i]]];
                        _nodal_data_f[t].vy[i] = (float)_flowSolution[_nsdofs*_flowPerm[indexes[i]]+1];
                        _nodal_data_f[t].vz[i] = (float)_flowSolution[_nsdofs*_flowPerm[indexes[i]]+2];
                    }
            }
            for (int i=0; i<nd; i++) {
                _nodal_data_f[t].perm[i] = _flowPerm[indexes[i]];
            }
        }
     }
    
    
    free_ivector(indexes, 0, solution.mesh.maxElementDofs-1);
    
    free_dvector(nodes->x, 0, mesh.maxElementNodes-1);
    free_dvector(nodes->y, 0, mesh.maxElementNodes-1);
    free_dvector(nodes->z, 0, mesh.maxElementNodes-1);
    free(nodes);
}

-(void)FEMFlowSolution_createMapGPUBuffersMesh:(FEMMesh * __nonnull)mesh solution:(FEMSolution * __nonnull)solution precisionMode:(NSString * __nonnull)precisionMode {
    
    int *colorMapping = NULL, *elementNodeIndexesStore = NULL;
    matrixArraysContainer *matContainers = NULL;
    variableArraysContainer *flowContainers = NULL;
    Nodes_t *meshNodes = NULL;
    cl_int err;
    
    flowContainers = solution.variable.getContainers;
    matContainers = solution.matrix.getContainers;
    meshNodes = mesh.getNodes;
    colorMapping = mesh.getColorMapping;
    elementNodeIndexesStore = mesh.getElementNodeIndexesStore;
    
    const char *type= "double";
    
    // ------------------------ Create the buffers on the device
    
    if ([solution.solutionInfo[@"use global basis functions coefficients"] boolValue] == YES) {
        if ([precisionMode isEqualToString:@"single"] == YES) {
            _basis_functions = clCreateBuffer(_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(basis_functions_f)*8, _basis_functions_f, &err);
        } else if ([precisionMode isEqualToString:@"double"] == YES) {
            _basis_functions = clCreateBuffer(_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(basis_functions_d)*8, _basis_functions_d, &err);
        }
    }
    
    _matDiag = clCreateBuffer(_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_int)*matContainers->sizeDiag, matContainers->Diag, &err);
    if (err < 0) {
        fatal("FEMFlowSolution:FEMFlowSolution_createMapGPUBuffersMesh", "Couldn't create the buffer object _matDiag.");
    }
    _matRows = clCreateBuffer(_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_int)*matContainers->sizeRows, matContainers->Rows, &err);
    if (err < 0) {
        fatal("FEMFlowSolution:FEMFlowSolution_createMapGPUBuffersMesh", "Couldn't create the buffer object _matRows.");
    }
    _matCols = clCreateBuffer(_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_int)*matContainers->sizeCols, matContainers->Cols, &err);
    if (err < 0) {
        fatal("FEMFlowSolution:FEMFlowSolution_createMapGPUBuffersMesh", "Couldn't create the buffer object _matCols.");
    }
    
    // Those buffers are initialized later inside the non-linearity loop
    _gpuData->matValues_v = alloc_mem(precision(), (size_t)matContainers->sizeValues);
    _gpuData->matRHS_v = alloc_mem(precision(), (size_t)matContainers->sizeRHS);
    
    _matValues = createDeviceBuffer(precision(), _context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, matContainers->sizeValues, _gpuData->matValues_v, err);
    if (err < 0) {
        fatal("FEMFlowSolution:FEMFlowSolution_createMapGPUBuffersMesh", "Couldn't create the buffer object _matValues.");
    }
    
    _matRHS = createDeviceBuffer(precision(), _context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, matContainers->sizeRHS, _gpuData->matRHS_v, err);
    if (err < 0) {
        fatal("FEMFlowSolution:FEMFlowSolution_createMapGPUBuffersMesh", "Couldn't create the buffer object _matRHS.");
    }
    
    // If we use the GPU local memory, we need to make sure that our work group size is a multiple
    // of the number of elemements in each color set. Adjust the number of elements in the color set
    // accordingly
    if ([solution.solutionInfo[@"use gpu local memory"] boolValue] == YES) {
        int adjust = 0;
        if (solution.solutionInfo[@"adjust global work size to be a multiple of"] != nil) {
            adjust = [solution.solutionInfo[@"adjust global work size to be a multiple of"] intValue];
        } else {
            fatal("FEMFlowSolution:FEMFlowSolution_createMapGPUBuffersMesh", "Missing parameter to adjust color sets.");
        }
        
        int nbElements, totElements;
        totElements = 0;
        for (NSMutableArray *color in mesh.colors) {
            nbElements = [color[0] intValue];
            nbElements = nbElements + (adjust - (nbElements & (adjust-1)));
            totElements += nbElements;
        }
        Element_t *elements = mesh.getElements;
        int *adjusted_colorMapping = calloc(totElements, sizeof(int));
        int indx = 0, add;
        for (NSMutableArray *color in mesh.colors) {
            // First the real elements for a given color
            for (int i=0; i<mesh.numberOfBulkElements; i++) {
                if (elements[i].color.colorIndex-1 == [color[1] intValue]) {
                    adjusted_colorMapping[indx] = elements[i].ElementIndex-1;
                    indx++;
                }
            }
            // Add the gost elements (element index=-1) for a given color
            nbElements = [color[0] intValue];
            add = (nbElements + (adjust - (nbElements & (adjust-1)))) - nbElements;
            for (int i=0; i<add; i++) {
                adjusted_colorMapping[indx] = -1;
                indx++;
            }
            color[0] = @(nbElements + (adjust - (nbElements & (adjust-1))));
        }
        fprintf(stdout, "FEMFlowSolution:FEMFlowSolution_createMapGPUBuffersMesh: number of elements after adjustement: \n");
        for (NSMutableArray *color in mesh.colors) {
            fprintf(stdout, "color index: %d: %d\n", [color[1] intValue], [color[0] intValue]);
        }
        _colorMapping = clCreateBuffer(_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_int)*totElements, adjusted_colorMapping, &err);
        free(adjusted_colorMapping);
        
    } else {
        _colorMapping = clCreateBuffer(_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_int)*mesh.numberOfBulkElements, colorMapping, &err);
    }
    if (err < 0) {
        fatal("FEMFlowSolution:FEMFlowSolution_createMapGPUBuffersMesh", "Couldn't create the buffer object _colorMapping.");
    }
    
    if ([solution.solutionInfo[@"use element nodal data"] boolValue] == YES) {
        if ([precisionMode isEqualToString:@"single"] == YES) {
            _nodal_data = clCreateBuffer(_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(nodal_data_f)*mesh.numberOfBulkElements, _nodal_data_f, &err);
        } else if ([precisionMode isEqualToString:@"double"] == YES) {
            _nodal_data = clCreateBuffer(_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(nodal_data_d)*mesh.numberOfBulkElements, _nodal_data_d, &err);
        }
    } else {
        _elementNodeIndexesStore = clCreateBuffer(_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_int)*(mesh.numberOfBulkElements*mesh.maxElementDofs), elementNodeIndexesStore, &err);
        if (err < 0) {
            fatal("FEMFlowSolution:FEMFlowSolution_createMapGPUBuffersMesh", "Couldn't create the buffer object _elementNodeIndexesStore.");
        }
        
        _gpuData->nodesX_v = alloc_mem(precision(), (size_t)mesh.numberOfNodes);
        _gpuData->nodesY_v = alloc_mem(precision(), (size_t)mesh.numberOfNodes);
        _gpuData->nodesZ_v = alloc_mem(precision(), (size_t)mesh.numberOfNodes);
        _gpuData->varSol_v = alloc_mem(precision(), (size_t)flowContainers->sizeValues);
        if ([precisionMode isEqualToString:@"single"] == YES) {
            
            _gpuData->nodesX_sp = _gpuData->nodesX_v;
            _gpuData->nodesY_sp = _gpuData->nodesY_v;
            _gpuData->nodesZ_sp = _gpuData->nodesZ_v;
            _gpuData->varSol_sp = _gpuData->varSol_v;
            
            size_t size = (size_t)mesh.numberOfNodes;
            err = init_data(precision(), type, 9, meshNodes->x, _gpuData->nodesX_sp, size, meshNodes->y,  _gpuData->nodesY_sp, size,
                            meshNodes->z, _gpuData->nodesZ_sp, size);
            
            size = (size_t)flowContainers->sizeValues;
            err |= init_data(precision(), type, 3, _flowSolution, _gpuData->varSol_sp, size);
            
        } else if ([precisionMode isEqualToString:@"double"] == YES) {
            
            _gpuData->nodesX_dp = _gpuData->nodesX_v;
            _gpuData->nodesY_dp = _gpuData->nodesY_v;
            _gpuData->nodesZ_dp = _gpuData->nodesZ_v;
            _gpuData->varSol_dp = _gpuData->varSol_v;
            
            size_t size = (size_t)mesh.numberOfNodes;
            err = init_data(precision(), type, 9, meshNodes->x, _gpuData->nodesX_dp, size, meshNodes->y,  _gpuData->nodesY_dp, size,
                            meshNodes->z, _gpuData->nodesZ_dp, size);
            
            size = (size_t)flowContainers->sizeValues;
            err |= init_data(precision(), type, 3, _flowSolution, _gpuData->varSol_dp, size);
            
        }
        if (err < 0) {
            fatal("FEMFlowSolution:FEMFlowSolution_createMapGPUBuffersMesh", "Error in setting data.");
        }
        
        _nodesX = createDeviceBuffer(precision(), _context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, mesh.numberOfNodes, _gpuData->nodesX_v, err);
        if (err < 0) {
            fatal("FEMFlowSolution:FEMFlowSolution_createMapGPUBuffersMesh", "Couldn't create the buffer object _nodesX.");
        }
        
        _nodesY = createDeviceBuffer(precision(), _context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, mesh.numberOfNodes, _gpuData->nodesY_v, err);
        if (err < 0) {
            fatal("FEMFlowSolution:FEMFlowSolution_createMapGPUBuffersMesh", "Couldn't create the buffer object _nodesY.");
        }
        
        _nodesZ = createDeviceBuffer(precision(), _context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, mesh.numberOfNodes, _gpuData->nodesZ_v, err);
        if (err < 0) {
            fatal("FEMFlowSolution:FEMFlowSolution_createMapGPUBuffersMesh", "Couldn't create the buffer object _nodesZ.");
        }
        
        _varSolution = createDeviceBuffer(precision(), _context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, flowContainers->sizeValues, _gpuData->varSol_v, err);
        if (err < 0) {
            fatal("FEMFlowSolution:FEMFlowSolution_createMapGPUBuffersMesh", "Couldn't create the buffer object _varSolution.");
        }
        
        _varPermutation = clCreateBuffer(_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_int)*flowContainers->sizePerm, _flowPerm, &err);
        if (err < 0) {
            fatal("FEMFlowSolution:FEMFlowSolution_createMapGPUBuffersMesh", "Couldn't create the buffer object _varPermutation.");
        }
    }
    
    _gpuNewtonLinear = clCreateBuffer(_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_char), &_newtonLinearization, &err);
    
    if ([solution.solutionInfo[@"compute basis and basis derivatives in separate kernel"] boolValue] == YES) {
        if ([precisionMode isEqualToString:@"single"] == YES) {
            _element_basis = clCreateBuffer(_context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(element_basis_f), _element_basis_f, &err);
            _element_dbasisdx = clCreateBuffer(_context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(element_dBasisdx_f)*mesh.numberOfBulkElements, _element_dBasisdx_f, &err);
        } else if ([precisionMode isEqualToString:@"double"] == YES) {
            _element_basis = clCreateBuffer(_context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(element_basis_d), _element_basis_d, &err);
            _element_dbasisdx = clCreateBuffer(_context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(element_dBasisdx_d)*mesh.numberOfBulkElements, _element_dBasisdx_d, &err);
        }
    }
    
    GlobalMemoryAllocationSize_t *globalAllocation = (GlobalMemoryAllocationSize_t *)malloc(sizeof(GlobalMemoryAllocationSize_t));
    initGlobalMemoryAllocation(globalAllocation);
    
    int nbOfReals = 0;
    if ([solution.solutionInfo[@"use element nodal data"] boolValue] == YES) {
        nbOfReals = 256 + 24 + 32 + matContainers->sizeValues + matContainers->sizeRHS;
    } else {
        nbOfReals = 256 + 24 + 32 + matContainers->sizeValues + matContainers->sizeRHS + mesh.numberOfNodes
        + mesh.numberOfNodes +  mesh.numberOfNodes + flowContainers->sizeValues;
    }
    
    if (precision() == 1) {
       globalAllocation->nb_float = nbOfReals;
    } else {
        globalAllocation->nb_double = nbOfReals;
    }
    
    if ([solution.solutionInfo[@"use element nodal data"] boolValue] == YES) {
        globalAllocation->nb_int = matContainers->sizeDiag + matContainers->sizeRows + matContainers->sizeCols + mesh.numberOfBulkElements;
    } else {
        globalAllocation->nb_int = matContainers->sizeDiag + matContainers->sizeRows + matContainers->sizeCols + mesh.numberOfBulkElements
        + (mesh.numberOfBulkElements*mesh.maxElementDofs) + flowContainers->sizePerm;
    }
    
    globalAllocation->nb_char = 1;
    
    size_t globalAllocationSize = computeGlobalMemoryAllocation(globalAllocation);
    
    if ([solution.solutionInfo[@"use element nodal data"] boolValue] == YES) {
        if (precision() == 1) {
            globalAllocationSize = globalAllocationSize + sizeof(nodal_data_f)*mesh.numberOfBulkElements;
        } else {
            globalAllocationSize = globalAllocationSize + sizeof(nodal_data_d)*mesh.numberOfBulkElements;
        }
    }
    
    if ([solution.solutionInfo[@"use global basis functions coefficients"] boolValue] == YES) {
        if (precision() == 1) {
            globalAllocationSize = globalAllocationSize + sizeof(basis_functions_f)*8;
        } else {
            globalAllocationSize = globalAllocationSize + sizeof(basis_functions_d)*8;
        }
    }
    
    if ([solution.solutionInfo[@"compute basis and basis derivatives in separate kernel"] boolValue] == YES) {
        if (precision() == 1) {
            globalAllocationSize = globalAllocationSize + sizeof(element_basis_f) + sizeof(element_dBasisdx_f)*mesh.numberOfBulkElements;
        } else {
            globalAllocationSize = globalAllocationSize + sizeof(element_basis_d) + sizeof(element_dBasisdx_d)*mesh.numberOfBulkElements;
        }
    }
    
    fprintf(stdout, "FEMFlowSolution:FEMFlowSolution_createMapGPUBuffersMesh: %lu MBytes allocated on the GPU global address space.\n",
            globalAllocationSize/(1024*1024));
    fprintf(stdout, "FEMFlowSolution:FEMFlowSolution_createMapGPUBuffersMesh: global matrix number of rows: %d.\n", solution.matrix.numberOfRows);
    fprintf(stdout, "FEMFlowSolution:FEMFlowSolution_createMapGPUBuffersMesh: global matrix number of non-zeros: %d.\n", matContainers->sizeValues);
}

-(void)FEMFlowSolution_setKernelArgumentsCore:(FEMCore * __nonnull)core  model:(FEMModel * __nonnull)model solution:(FEMSolution * __nonnull)solution precisionMode:(NSString * __nonnull)precisionMode getActiveElement:(Element_t* (* __nonnull)(id, SEL, int, FEMSolution*, FEMModel*))getActiveElementIMP {
    
    cl_int err;
    cl_uint argIdx = 0;
    
    err  = clSetKernelArg(_kernel_assembly, argIdx, sizeof(cl_mem), &_matValues);
    err |= clSetKernelArg(_kernel_assembly, ++argIdx, sizeof(cl_mem), &_matRHS);
    err |= clSetKernelArg(_kernel_assembly, ++argIdx, sizeof(cl_mem), &_matDiag);
    err |= clSetKernelArg(_kernel_assembly, ++argIdx, sizeof(cl_mem), &_matRows);
    err |= clSetKernelArg(_kernel_assembly, ++argIdx, sizeof(cl_mem), &_matCols);
    err |= clSetKernelArg(_kernel_assembly, ++argIdx, sizeof(cl_mem), &_colorMapping);
    if ([solution.solutionInfo[@"use element nodal data"] boolValue] == YES) {
        err |= clSetKernelArg(_kernel_assembly, ++argIdx, sizeof(cl_mem), &_nodal_data);
    } else {
        err |= clSetKernelArg(_kernel_assembly, ++argIdx, sizeof(cl_mem), &_elementNodeIndexesStore);
        err |= clSetKernelArg(_kernel_assembly, ++argIdx, sizeof(cl_mem), &_nodesX);
        err |= clSetKernelArg(_kernel_assembly, ++argIdx, sizeof(cl_mem), &_nodesY);
        err |= clSetKernelArg(_kernel_assembly, ++argIdx, sizeof(cl_mem), &_nodesZ);
        err |= clSetKernelArg(_kernel_assembly, ++argIdx, sizeof(cl_mem), &_varSolution);
        err |= clSetKernelArg(_kernel_assembly, ++argIdx, sizeof(cl_mem), &_varPermutation);
    }
    err |= clSetKernelArg(_kernel_assembly, ++argIdx, sizeof(cl_mem), &_gpuNewtonLinear);
    
    if (solution.solutionInfo[@"gpu ice density"] == nil) fatal("FEMFlowSolution:FEMFlowSolution_setKernelArgumentsCore", "Missing ice density value for gpu execution.");
    if (solution.solutionInfo[@"gpu ice viscosity"] == nil) fatal("FEMFlowSolution:FEMFlowSolution_setKernelArgumentsCore", "Missing ice viscosity value for gpu execution.");
    if (solution.solutionInfo[@"gpu ice gravity"] == nil) fatal("FEMFlowSolution:FEMFlowSolution_setKernelArgumentsCore", "Missing gravity value for gpu execution.");
    
    Element_t *element = getActiveElementIMP(core, @selector(getActiveElement:solution:model:), 0, solution, model);
    double hScale;
    if (solution.solutionInfo[@"h scale"] != nil) {
        hScale = [solution.solutionInfo[@"h scale"] doubleValue];
    } else {
        hScale = 1.0;
    }
    
    _gpuData->density_v = alloc_mem(precision(), 1);
    _gpuData->viscosity_v = alloc_mem(precision(), 1);
    _gpuData->gravity_v = alloc_mem(precision(), 1);
    _gpuData->hk_v = alloc_mem(precision(), 1);
    _gpuData->mk_v = alloc_mem(precision(), 1);
    
    if ([precisionMode isEqualToString:@"single"] == YES) {
        
        _gpuData->density_sp = (float)[solution.solutionInfo[@"gpu ice density"] doubleValue];
        memcpy(_gpuData->density_v, &_gpuData->density_sp, sizeof(float));
        
        _gpuData->viscosity_sp = (float)[solution.solutionInfo[@"gpu ice viscosity"] doubleValue];
        memcpy(_gpuData->viscosity_v, &_gpuData->viscosity_sp, sizeof(float));
        
        _gpuData->gravity_sp = (float)[solution.solutionInfo[@"gpu ice gravity"] doubleValue];
        memcpy(_gpuData->gravity_v, &_gpuData->gravity_sp, sizeof(float));
        
        _gpuData->hk_sp = (float)element->hK * hScale;
        memcpy(_gpuData->hk_v, &_gpuData->hk_sp, sizeof(float));
        
         _gpuData->mk_sp = (float)element->StabilizationMK;
        memcpy(_gpuData->mk_v, &_gpuData->mk_sp, sizeof(float));
        
    } else if ([precisionMode isEqualToString:@"double"] == YES) {
        
        _gpuData->density_dp = [solution.solutionInfo[@"gpu ice density"] doubleValue];
        memcpy(_gpuData->density_v, &_gpuData->density_dp, sizeof(double));
        
        _gpuData->viscosity_dp = [solution.solutionInfo[@"gpu ice viscosity"] doubleValue];
        memcpy(_gpuData->viscosity_v, &_gpuData->viscosity_dp, sizeof(double));
        
        _gpuData->gravity_dp = [solution.solutionInfo[@"gpu ice gravity"] doubleValue];
        memcpy(_gpuData->gravity_v, &_gpuData->gravity_dp, sizeof(double));
        
        _gpuData->hk_dp = element->hK * hScale;
        memcpy(_gpuData->hk_v, &_gpuData->hk_dp, sizeof(double));
        
        _gpuData->mk_dp = element->StabilizationMK;
        memcpy(_gpuData->mk_v, &_gpuData->mk_dp, sizeof(double));
        
    }
    
    err |= setDeviceKernelArg(precision(), _kernel_assembly, ++argIdx, 1, *_gpuData->density_v);
    err |= setDeviceKernelArg(precision(), _kernel_assembly, ++argIdx, 1, *_gpuData->viscosity_v);
    err |= setDeviceKernelArg(precision(), _kernel_assembly, ++argIdx, 1, *_gpuData->gravity_v);
    err |= setDeviceKernelArg(precision(), _kernel_assembly, ++argIdx, 1, *_gpuData->hk_v);
    err |= setDeviceKernelArg(precision(), _kernel_assembly, ++argIdx, 1, *_gpuData->mk_v);
    
    _positionInColorMapping = ++argIdx;
    
    int n = element->Type.NumberOfNodes;
    int nBasis = element->Type.NumberOfNodes;
    err |= clSetKernelArg(_kernel_assembly, ++argIdx, sizeof(cl_int), &n);
    err |= clSetKernelArg(_kernel_assembly, ++argIdx, sizeof(cl_int), &n);
    err |= clSetKernelArg(_kernel_assembly, ++argIdx, sizeof(cl_int), &nBasis);
    
    int varDofs = solution.variable.dofs;
    err |= clSetKernelArg(_kernel_assembly, ++argIdx, sizeof(cl_int), &varDofs);
    
    // If we use basis function from global memory instead of the definitions
    // inside the kernel
    if ([solution.solutionInfo[@"use global basis functions coefficients"] boolValue] == YES) {
        err |= clSetKernelArg(_kernel_assembly, ++argIdx, sizeof(cl_mem), &_basis_functions);
    }
    
    // If we compute the basis and basis derivatives in another kernel
    if ([solution.solutionInfo[@"compute basis and basis derivatives in separate kernel"] boolValue] == YES) {
        err |= clSetKernelArg(_kernel_assembly, ++argIdx, sizeof(cl_mem), &_element_basis);
        err |= clSetKernelArg(_kernel_assembly, ++argIdx, sizeof(cl_mem), &_element_dbasisdx);
    }
    
    // Allocate space for local memory if required
    if ([solution.solutionInfo[@"use gpu local memory"] boolValue] == YES) {
        int wrk_groupSize = 0;
        if (solution.solutionInfo[@"parallel assembly work-group size"] != nil) {
            wrk_groupSize  = [solution.solutionInfo[@"parallel assembly work-group size"] intValue];
        }
        
        if ([solution.solutionInfo[@"compute basis and basis derivatives in separate kernel"] boolValue] == YES) {
            if ([precisionMode isEqualToString:@"single"] == YES) {
                err |= clSetKernelArg(_kernel_assembly, ++argIdx, wrk_groupSize*sizeof(element_info_extended_f), NULL);
            } else if ([precisionMode isEqualToString:@"double"] == YES) {
                err |= clSetKernelArg(_kernel_assembly, ++argIdx, wrk_groupSize*sizeof(element_info_extended_d), NULL);
            }
        } else {
            if ([precisionMode isEqualToString:@"single"] == YES) {
                err |= clSetKernelArg(_kernel_assembly, ++argIdx, wrk_groupSize*sizeof(element_info_f), NULL);
            } else if ([precisionMode isEqualToString:@"double"] == YES) {
                err |= clSetKernelArg(_kernel_assembly, ++argIdx, wrk_groupSize*sizeof(element_info_d), NULL);
            }
        }
    }
    
    if (err < 0) {
        fatal("FEMFlowSolution:FEMFlowSolution_setKernelArgumentsCore", "Error in setting < assembly > kernel arguments.");
    }
    
    if ([solution.solutionInfo[@"compute basis and basis derivatives in separate kernel"] boolValue] == YES) {
        int numberOfElements = model.numberOfBulkElements;
        argIdx = 0;
        err  = clSetKernelArg(_kernel_basis_dbasisdx, argIdx, sizeof(cl_mem), &_element_basis);
        err |= clSetKernelArg(_kernel_basis_dbasisdx, ++argIdx, sizeof(cl_mem), &_element_dbasisdx);
        if ([solution.solutionInfo[@"use element nodal data"] boolValue] == YES) {
            err |= clSetKernelArg(_kernel_basis_dbasisdx, ++argIdx, sizeof(cl_mem), &_nodal_data);
        } else {
            err |= clSetKernelArg(_kernel_basis_dbasisdx, ++argIdx, sizeof(cl_mem), &_nodesX);
            err |= clSetKernelArg(_kernel_basis_dbasisdx, ++argIdx, sizeof(cl_mem), &_nodesY);
            err |= clSetKernelArg(_kernel_basis_dbasisdx, ++argIdx, sizeof(cl_mem), &_nodesZ);
            err |= clSetKernelArg(_kernel_basis_dbasisdx, ++argIdx, sizeof(cl_mem), &_elementNodeIndexesStore);
        }
        err |= clSetKernelArg(_kernel_basis_dbasisdx, ++argIdx, sizeof(cl_int), &n);
        err |= clSetKernelArg(_kernel_basis_dbasisdx, ++argIdx, sizeof(cl_int), &n);
        err |= clSetKernelArg(_kernel_basis_dbasisdx, ++argIdx, sizeof(cl_int), &numberOfElements);
        
        if ([solution.solutionInfo[@"use global basis functions coefficients"] boolValue] == YES) {
            err |= clSetKernelArg(_kernel_basis_dbasisdx, ++argIdx, sizeof(cl_mem), &_basis_functions);
        }
        
        if ([solution.solutionInfo[@"use gpu local memory"] boolValue] == YES) {
            int wrk_groupSize = 0;
            if (solution.solutionInfo[@"parallel assembly work-group size"] != nil) {
                wrk_groupSize  = [solution.solutionInfo[@"parallel assembly work-group size"] intValue];
            }
            if ([precisionMode isEqualToString:@"single"] == YES) {
                err |= clSetKernelArg(_kernel_basis_dbasisdx, ++argIdx, wrk_groupSize*sizeof(element_info_extended_f), NULL);
            } else if ([precisionMode isEqualToString:@"double"] == YES) {
                err |= clSetKernelArg(_kernel_basis_dbasisdx, ++argIdx, wrk_groupSize*sizeof(element_info_extended_d), NULL);
            }
        }
        if (err < 0) {
            fatal("FEMFlowSolution:FEMFlowSolution_setKernelArgumentsCore", "Error in setting < basis and basis derivatives >kernel arguments.");
        }
    }
}

-(void)FEMFlowSolution_setGPUCore:(FEMCore * __nonnull)core model:(FEMModel * __nonnull)model solution:(FEMSolution * __nonnull)solution mesh:(FEMMesh * __nonnull)mesh precisionMode:(NSString * __nonnull)precisionMode transientSimulation:(BOOL)transient getActiveElement:(Element_t* (* __nonnull)(id, SEL, int, FEMSolution*, FEMModel*))getActiveElementIMP getNodes:(void (* __nonnull)(id, SEL, FEMSolution*, FEMModel*, Element_t*, Nodes_t*, int*, FEMMesh*))getNodesIMP getElementDofsSolution:(int (* __nonnull)(id, SEL, FEMSolution*, FEMModel*, Element_t*, int*, BOOL*))getElementDofsSolutionIMP {
    
    cl_int err;
    
    // Get the GPU device
    _device = find_single_device();
    fprintf(stdout, "FEMFlowSolution:FEMFlowSolution_setGPUCore: parrallel assembly computed on GPU info: \n");
    device_info(_device);
    
    // Create the context of the command queue
    _context = clCreateContext(0, 1, &_device, NULL, NULL, &err);
    if (err < 0) {
        fatal("FEMFlowSolution:FEMFlowSolution_setGPUCore", "Can't create context for device. Error: ", err);
    }
    _cmd_queue = clCreateCommandQueue(_context, _device, CL_QUEUE_PROFILING_ENABLE, NULL);
    
    NSString *kernelfile;
    if (solution.solutionInfo[@"gpu kernel source file"] != nil) {
        NSString *source = solution.solutionInfo[@"gpu kernel source file"];
        if ([source containsString:@".cl"] == YES) {
            kernelfile = [NSString stringWithString:source];
        } else {
            kernelfile = [source stringByAppendingPathExtension:@"cl"];
        }
    } else {
        fatal("FEMFlowSolution:FEMFlowSolution_setGPUCore", "Parallel assembly requires a kernel source file.");
    }
    
    int success = LoadFileIntoString((char *)[kernelfile UTF8String], &_kernel_source, &_src_len);
    if (success < 0) {
        fatal("FEMFlowSolution:FEMFlowSolution_setGPUCore", "Can't load kernel source.");
    }
    
    // Allocate memory for program and kernels
    // Create the program .cl file
    _program = clCreateProgramWithSource(_context, 1, (const char**)&_kernel_source, NULL, &err);
    if (err < 0) {
        fatal("FEMFlowSolution:FEMFlowSolution_setGPUCore", "Can't create program. Error: ", err);
    }
    
    // Build the program (compile it)
    
    // So far we only supports 3D problems
    if (mesh.dimension < 3 || model.dimension < 3) {
        fatal("FEMFlowSolution:FEMFlowSolution_setGPUCore: unsupported mesh or model dimension for GPU assembly.");
    }
    
    // Options passed to the compilation of the kernel
    NSMutableString *listOfOptions;
    listOfOptions = [NSMutableString stringWithString:@"-DKERNEL_FP_64 -DMESH_DIMENSION_3 -DELEMENT_DIMENSION_3 -DMODEL_DIMENSION_3"];
    if ([solution.solutionInfo[@"compute basis and basis derivatives in separate kernel"] boolValue] == YES) {
        [listOfOptions appendString:@" -DKERNEL_BASIS_DBASISDX"];
    }
    if ([solution.solutionInfo[@"use gpu local memory"] boolValue] == YES) {
        [listOfOptions appendString:@" -DUSE_GPU_LOCAL_MEM"];
    }
    if ([solution.solutionInfo[@"use global basis functions coefficients"] boolValue] == YES) {
        [listOfOptions appendString:@" -DGLOBAL_BASIS_FUNCTIONS"];
    }
    if ([solution.solutionInfo[@"use element nodal data"] boolValue] == YES) {
        [listOfOptions appendString:@" -DNODAL_DATA"];
    }
    if ([precisionMode isEqualToString:@"single"] == YES) {
        [listOfOptions replaceOccurrencesOfString:@"-DKERNEL_FP_64" withString:@"-DKERNEL_FP_32" options:NSBackwardsSearch range:NSMakeRange(0, [listOfOptions length])];
    }
    
    // Some extra flags
    
    // Do we use DEBUG mode in the kernel
    if ([solution.solutionInfo[@"enable GPU debug mode"] boolValue] == YES) {
        [listOfOptions appendString:@" -DDEBUG"];
    }
    
    // Do we have a transient run
    if (transient == YES) {
        [listOfOptions appendString:@" -DTRANSIENT"];
    }
    
    // Do we use the Newton Linearization
    if ([solution.solutionInfo[@"enable newton linearization"] boolValue] == YES) {
        [listOfOptions appendString:@" -DNEWTONLINEAR"];
    }
    
    // Do we have negative permutation indexes
    variableArraysContainer *flowContainers = solution.variable.getContainers;
    BOOL anyNegativePermIndx = NO;
    for (int i=0; i<flowContainers->sizePerm; i++) {
        if (_flowPerm[i] < 0) {
            anyNegativePermIndx = YES;
            break;
        }
    }
    if (anyNegativePermIndx == YES) {
        [listOfOptions appendString:@" -DCHECK_NEG_PERM"];
    }
    
    // Use MAD
    if ([solution.solutionInfo[@"enable gpu multiply-and-add operations"] boolValue] == YES) {
        [listOfOptions appendString:@" -cl-mad-enable"];
    }
    // Flush denorms to zero
    if ([solution.solutionInfo[@"disable processing of denormalized numbers"] boolValue] == YES) {
        [listOfOptions appendString:@" -cl-denorms-are-zero"];
    }
    // Relax IEEE compliance
    if ([solution.solutionInfo[@"relax IEEE compliance"] boolValue] == YES) {
        [listOfOptions appendString:@" -cl-fast-relaxed-math"];
    }
    const char *options = [listOfOptions UTF8String];
    
    fprintf(stdout, "FEMFlowSolution:FEMFlowSolution_setGPUCore: build GPU program...\n");
    err = clBuildProgram(_program, 1, &_device, options, NULL, &err);
    if (err < 0) {
        size_t log_size;
        clGetProgramBuildInfo(_program, _device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
        char *program_log = (char *)malloc(log_size+1);
        program_log[log_size] = '\0';
        clGetProgramBuildInfo(_program, _device, CL_PROGRAM_BUILD_LOG, log_size+1, program_log, NULL);
        fprintf(stderr, "FEMFlowSolution:FEMFlowSolution_setGPUCore: error when buiding the GPU kernel.\n");
        fprintf(stderr, "FEMFlowSolution:FEMFlowSolution_setGPUCore: log: %s\n", program_log);
        free(program_log);
        fatal("FEMFlowSolution:FEMFlowSolution_setGPUCore");
    }
    fprintf(stdout, "FEMFlowSolution:FEMFlowSolution_setGPUCore: done.\n");
    
    // Create the kernel(s)
    _kernel_assembly = clCreateKernel(_program, "NavierStokesCompose", &err);
    if (err < 0) {
        fatal("FEMFlowSolution:FEMFlowSolution_setGPUCore", "Can't create kernel. Error: ", err);
    }
    if ([solution.solutionInfo[@"compute basis and basis derivatives in separate kernel"] boolValue] == YES) {
        _kernel_basis_dbasisdx = clCreateKernel(_program, "ComputeBasisDBasisdx", &err);
        if (err < 0) {
            fatal("FEMFlowSolution:FEMFlowSolution_setGPUCore", "Can't create kernel. Error: ", err);
        }
    }
    
    // Create data structures needed by the GPU
    
    _gpuData = (ice_flow_gpu *)malloc(sizeof(ice_flow_gpu));
    init_gpu_data(_gpuData);
    
    if ([solution.solutionInfo[@"use global basis functions coefficients"] boolValue] == YES) {
        [self FEMFlowSolution_initElementBasisFunctions:core solution:solution model:model mesh:mesh precisionMode:precisionMode getActiveElement:getActiveElementIMP];
    }
    
    if ([solution.solutionInfo[@"compute basis and basis derivatives in separate kernel"] boolValue] == YES) {
        if (precision() == 1) {
            _element_basis_f = (element_basis_f *)malloc(sizeof(element_basis_f));
            _element_dBasisdx_f = (element_dBasisdx_f *)malloc(sizeof(element_dBasisdx_f)*mesh.numberOfBulkElements);
            
        } else {
            _element_basis_d = (element_basis_d *)malloc(sizeof(element_basis_d));
            _element_dBasisdx_d = (element_dBasisdx_d *)malloc(sizeof(element_dBasisdx_d)*mesh.numberOfBulkElements);
        }
    }
    
    if ([solution.solutionInfo[@"use element nodal data"] doubleValue] == YES) {
        [self FEMFlowSolution_iniNodalData:core solution:solution model:model mesh:mesh precisionMode:precisionMode getActiveElement:getActiveElementIMP getNodes:getNodesIMP getElementDofsSolution:getElementDofsSolutionIMP];
    }
    
    // Allocate the GPU buffers and map the written buffers to the host
    [self FEMFlowSolution_createMapGPUBuffersMesh:mesh solution:solution precisionMode:precisionMode];
    
    // Set kernel arguments
    [self FEMFlowSolution_setKernelArgumentsCore:core model:model solution:solution precisionMode:precisionMode getActiveElement:getActiveElementIMP];
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
        _context = NULL;
        _cmd_queue = NULL;
        _program = NULL;
        _kernel_assembly = NULL;
        _kernel_basis_dbasisdx = NULL;
        _kernel_source = NULL;
        _initializeGPU = NO;
        
        _gpuData = NULL;
        _basis_functions_f = NULL;
        _basis_functions_d = NULL;
        
        _element_basis_f = NULL;
        _element_basis_d = NULL;
        
        _element_dBasisdx_f = NULL;
        _element_dBasisdx_d = NULL;
        
        _nodal_data_f = NULL;
        _nodal_data_d = NULL;
        
        _basis_functions = NULL;
        
        _matDiag = NULL;
        _matRows = NULL;
        _matCols = NULL;
        _matValues = NULL;
        _matRHS = NULL;
        _colorMapping = NULL;
        _elementNodeIndexesStore = NULL;
        _nodesX = NULL;
        _nodesY = NULL;
        _nodesZ = NULL;
        _varSolution = NULL;
        _varPermutation = NULL;
        _gpuNewtonLinear = NULL;
        _element_basis = NULL;
        _element_dbasisdx = NULL;
        _nodal_data = NULL;
        
        _mapped_matValues = NULL;
        _mapped_matRHS = NULL;
        _mapped_varSolution = NULL;
        _mapped_gpuNewton = NULL;
        _mapped_nodal_data = NULL;
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
    
    if (solution.solutionInfo[@"parallel assembly"] != nil) {
        if ([solution.solutionInfo[@"parallel assembly"] boolValue] == YES) {
                
            if (_kernel_assembly != NULL) clReleaseKernel(_kernel_assembly);
            if (_kernel_basis_dbasisdx != NULL) clReleaseKernel(_kernel_basis_dbasisdx);
            if (_context != NULL) clReleaseContext(_context);
            if (_cmd_queue != NULL) clReleaseCommandQueue(_cmd_queue);
            if (_program != NULL) clReleaseProgram(_program);
            
            if (_basis_functions != NULL) clReleaseMemObject(_basis_functions);
            if (_matDiag != NULL) clReleaseMemObject(_matDiag);
            if (_matRows != NULL) clReleaseMemObject(_matRows);
            if (_matCols != NULL) clReleaseMemObject(_matCols);
            if (_matValues != NULL) clReleaseMemObject(_matValues);
            if (_matRHS != NULL) clReleaseMemObject(_matRHS);
            if (_colorMapping != NULL) clReleaseMemObject(_colorMapping);
            if (_elementNodeIndexesStore != NULL) clReleaseMemObject(_elementNodeIndexesStore);
            if (_nodesX != NULL) clReleaseMemObject(_nodesX);
            if (_nodesY != NULL) clReleaseMemObject(_nodesY);
            if (_nodesZ != NULL) clReleaseMemObject(_nodesZ);
            if (_varSolution != NULL) clReleaseMemObject(_varSolution);
            if (_varPermutation != NULL) clReleaseMemObject(_varPermutation);
            if (_gpuNewtonLinear != NULL) clReleaseMemObject(_gpuNewtonLinear);
            if (_element_basis != NULL) clReleaseMemObject(_element_basis);
            if (_element_dbasisdx != NULL) clReleaseMemObject(_element_dbasisdx);
            
            if (_gpuData != NULL) {
                if (_gpuData->nodesX_v != NULL) free(_gpuData->nodesX_v);
                if (_gpuData->nodesY_v != NULL) free(_gpuData->nodesY_v);
                if (_gpuData->nodesZ_v != NULL) free(_gpuData->nodesZ_v);
                if (_gpuData->varSol_v != NULL) free(_gpuData->varSol_v);
                if (_gpuData->matValues_v != NULL) free(_gpuData->matValues_v);
                if (_gpuData->matRHS_v != NULL) free(_gpuData->matRHS_v);
                if (_gpuData->density_v != NULL) free(_gpuData->density_v);
                if (_gpuData->viscosity_v != NULL) free(_gpuData->viscosity_v);
                if (_gpuData->gravity_v != NULL) free(_gpuData->gravity_v);
                if (_gpuData->hk_v != NULL) free(_gpuData->hk_v);
                if (_gpuData->mk_v != NULL) free(_gpuData->mk_v);
                free(_gpuData);
            }
            if (_basis_functions_f != NULL) {
                free(_basis_functions_f);
            }
            if (_basis_functions_d != NULL) {
                free(_basis_functions_d);
            }
            
            if (_element_basis_f != NULL) {
                free(_element_basis_f);
            }
            if (_element_dBasisdx_f != NULL) {
                free(_element_dBasisdx_f);
            }
            if (_element_basis_d != NULL) {
                free(_element_basis_d);
            }
            if(_element_dBasisdx_d != NULL) {
                free(_element_dBasisdx_d);
            }
            
            if (_nodal_data_f != NULL) free(_nodal_data_f);
            if (_nodal_data_d != NULL) free(_nodal_data_d);
            
            if (_kernel_source != NULL) free(_kernel_source);
        }
    }
}

-(void)solutionComputer:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model timeStep:(int)timeStep transientSimulation:(BOOL)transient {
    
    int i, j, k, n, t, bf_id, body_id, eq_id, mat_id, compressibilityModel=-1, dim, freeSIter, iter, modelCoords=0, modelDim=0, newtonIter, nonLinearIter;
    static int dt, saveTimeStep=-1;
    int *tempPerm = NULL, *meshVeloPerm = NULL;
    double *temperature = NULL, *tempPrev = NULL, *meshVelocity = NULL;
    double at, at0, at1, freeSTol, gravity[3], nonLinearRelax, nonLinearTol, newtonTol, pseudoCompressibilityScale, relativeChange, relaxation, st, sum, totat, totst, uNorm;
    NSString *compressibilityFlag, *flowModel, *localCoords, *stabilizeFlag, *varName;
    BOOL bubbles, convect, computeFree = NO, divDiscretization, found, freeSurfaceFlag, gradPDiscretization, gotForceBC, ifTransient, mbFlag, normalTangential, pseudoPressureExists=NO, pseudoPressureUpdate=NO, relaxBefore, stabilize, useLocalCoords = NO;
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
    
    fprintf(stdout, "FEMFlowSolution:solutionComputer: solving the Navier-Stokes equations.\n");
    
    _newtonLinearization = NO;
    
    if (solution.solutionInfo[@"solver coordinate system"] != nil) {
        localCoords = solution.solutionInfo[@"solver coordinate system"];
        useLocalCoords = YES;
    }
    if (useLocalCoords == YES) {
        modelCoords = coordinatesSystems.coordinates;
        modelDim = model.dimension;
        if ([localCoords isEqualToString:@"cartesian 2d"] == YES) {
            coordinatesSystems.coordinates = 1;
            model.dimension = 2;
            fprintf(stdout, "FEMFlowSolution:solutionComputer: Solution coordinate system is cartesian 2D.\n");
        } else if ([localCoords isEqualToString:@"cartesian 3d"] == YES) {
            coordinatesSystems.coordinates = 1;
            model.dimension = 3;
            fprintf(stdout, "FEMFlowSolution:solutionComputer: Solution coordinate system is cartesian 3D.\n");
        } else if ([localCoords isEqualToString:@"axi symmetric"] == YES) {
            coordinatesSystems.coordinates = 4;
            model.dimension = 2;
            fprintf(stdout, "FEMFlowSolution:solutionComputer: Solution coordinate system is axi symmetric.\n");
        } else if ([localCoords isEqualToString:@"cylindric symmetric"] == YES) {
            coordinatesSystems.coordinates = 3;
            model.dimension = 3;
            fprintf(stdout, "FEMFlowSolution:solutionComputer: Solution coordinate system is cylindric symmetric.\n");
        } else {
            fprintf(stdout, "FEMFlowSolution:solutionComputer: Solution coordinate system not recognized, using original.\n");
        }
    }
    
    // Check for flow model. one of 'full', 'no convection', 'stokes'
    ifTransient = transient;
    convect = YES;
    if (solution.solutionInfo[@"flow model"] != nil) {
        flowModel = solution.solutionInfo[@"flow model"];
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
    
    bubbles = [solution.solutionInfo[@"bubbles"] boolValue];
    stabilize = [(solution.solutionInfo[@"stabilize"]) boolValue];
    
    if (solution.solutionInfo[@"stabilization method"] != nil) {
        stabilizeFlag = solution.solutionInfo[@"stabilization method"];
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
    
    divDiscretization = [solution.solutionInfo[@"div dicretization"] boolValue];
    gradPDiscretization = [solution.solutionInfo[@"gradp discretization"] boolValue];
    nonLinearTol = [solution.solutionInfo[@"nonlinear system convergence tolerance"] doubleValue];
    if (nonLinearTol < 0.0) nonLinearTol = 0.0;
    newtonTol = [solution.solutionInfo[@"nonlinear system newton after tolerance"] doubleValue];
    if (newtonTol < 0.0) newtonTol = 0.0;
    
    newtonIter = [solution.solutionInfo[@"nonlinear system newton after iterations"] intValue];
    if (newtonIter == 0) _newtonLinearization = YES;
    
    if ([solution.solutionInfo[@"nonlinear system reset newton"] boolValue] == YES) _newtonLinearization = NO;
    
    nonLinearIter = [solution.solutionInfo[@"nonlinear system maximum iterations"] intValue];
    if (nonLinearIter < 0) nonLinearIter = 0;
    
    if (solution.solutionInfo[@"nonlinear system norm dofs"] == nil) {
        [solution.solutionInfo setObject:@(_nsdofs-1) forKey:@"nonlinear system norm dofs"];
    }
    
    if (solution.solutionInfo[@"free surface after tolerance"] != nil) {
        freeSTol = [solution.solutionInfo[@"free surface after tolerance"] doubleValue];
    } else {
        freeSTol = DBL_MAX;
    }
    
    if (solution.solutionInfo[@"free surface after iterations"] != nil) {
        freeSIter = [solution.solutionInfo[@"free surface after iterations"] intValue];
    } else {
        freeSIter = 0;
    }
    
    // We do our own relaxation
    if (solution.solutionInfo[@"nonlinear system relaxation factor"] != nil) {
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
        fprintf(stdout, "FEMFlowSolution:solutionComputer: pseudoPressure mean: %f.\n", sum/_sizePseudoPressure);
        
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
    
    FEMNumericIntegration *integration = [[FEMNumericIntegration alloc] init];
    if ([integration allocation:mesh] == NO) fatal("FEMFlowSolution:solutionComputer", "Allocation error in FEMNumericIntegration.");
    
    // If we do the matrix assembly on the GPU, set-up everything now (only once)
    if (parallelAssembly == YES) {
        if (_initializeGPU == NO) {
            // Check whether we support the type of flow we want to compute on GPU
            if (solution.solutionInfo[@"gpu flow type"] != nil) {
                if ([solution.solutionInfo[@"gpu flow type"] isEqualToString:@"ice flow"] == NO) {
                    fatal("FEMFlowSolution:solutionComputer", "Currently only the value <ice flow > for the < gpu flow type> parameter is supported for parallel assembly.");
                }
            } else fatal("FEMFlowSolution:solutionComputer", "Parameter < gpu flow type > not given.");
            
            if ([solution.solutionInfo[@"gpu floating-point precision"] isEqualToString:@"single"] == YES) {
                setPrecision(true);
                fprintf(stdout, "FEMFlowSolution:solutionComputer: single precision mode used in GPU solver.\n");
                _precision = @"single";
            } else if ([solution.solutionInfo[@"gpu floating-point precision"] isEqualToString:@"double"] == YES) {
                fprintf(stdout, "FEMFlowSolution:solutionComputer: double precision mode used in GPU solver.\n");
                _precision = @"double";
            } else {
                fatal("FEMFlowSolution:solutionComputer", "Unknown GPU precision mode.");
            }
            
            [self FEMFlowSolution_setGPUCore:core model:model solution:solution mesh:mesh precisionMode:_precision transientSimulation:transient getActiveElement:getActiveElementIMP getNodes:getNodesIMP getElementDofsSolution:getElementDofsSolutionIMP];
            _initializeGPU = YES;
        }
    }
    if (parallelAssembly == YES && transient == YES) {
        cl_int err;
        char *nonLinNewton = (char *)&_newtonLinearization;
        _mapped_gpuNewton = clEnqueueMapBuffer(_cmd_queue, _gpuNewtonLinear, CL_TRUE, CL_MAP_WRITE, 0, sizeof(cl_char), 0, NULL, NULL, &err);
        if (err < 0) {
            fatal("FEMFlowSolution:FEMFlowSolution_createMapGPUBuffersMesh", "Couldn't map _mapped_gpuNewton.");
        }
        memcpy(_mapped_gpuNewton, nonLinNewton, sizeof(cl_char));
        clEnqueueUnmapMemObject(_cmd_queue, _gpuNewtonLinear, _mapped_gpuNewton, 0, NULL, NULL);
    }
    
    if (parallelAssembly == YES) {
        cl_int err;
        size_t wg_size, wg_multiple;
        cl_ulong private_usage, local_usage;
        err = clGetKernelWorkGroupInfo(_kernel_assembly, _device, CL_KERNEL_WORK_GROUP_SIZE,
                                       sizeof(wg_size), &wg_size, NULL);
        err |= clGetKernelWorkGroupInfo(_kernel_assembly, _device,
                                        CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
                                        sizeof(wg_multiple), &wg_multiple, NULL);
        err |= clGetKernelWorkGroupInfo(_kernel_assembly, _device, CL_KERNEL_LOCAL_MEM_SIZE,
                                        sizeof(local_usage), &local_usage, NULL);
        err |= clGetKernelWorkGroupInfo(_kernel_assembly, _device, CL_KERNEL_PRIVATE_MEM_SIZE,
                                        sizeof(private_usage), &private_usage, NULL);
        if(err < 0) {
            fatal("FEMFlowSolution:solutionComputer", "Error in getting kernel work-group size information.");
        };
        fprintf(stdout, "FEMFlowSolution:solutionComputer: < assembly > kernel maximum work group size: %zu\n", wg_size);
        fprintf(stdout, "FEMFlowSolution:solutionComputer: < assembly > kernel work group size best multiple: %zu\n", wg_multiple);
        fprintf(stdout, "FEMFlowSolution:solutionComputer: local memory used by the < assembly > kernel (KB): %llu\n", local_usage/1024);
        fprintf(stdout, "FEMFlowSolution:solutionComputer: private memory used by the < assembly > kernel (KB): %llu\n", private_usage/1024);
        
        if ([solution.solutionInfo[@"compute basis and basis derivatives in separate kernel"] boolValue] == YES) {
            err = clGetKernelWorkGroupInfo(_kernel_basis_dbasisdx, _device, CL_KERNEL_WORK_GROUP_SIZE,
                                           sizeof(wg_size), &wg_size, NULL);
            err |= clGetKernelWorkGroupInfo(_kernel_basis_dbasisdx, _device,
                                            CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
                                            sizeof(wg_multiple), &wg_multiple, NULL);
            err |= clGetKernelWorkGroupInfo(_kernel_basis_dbasisdx, _device, CL_KERNEL_LOCAL_MEM_SIZE,
                                            sizeof(local_usage), &local_usage, NULL);
            err |= clGetKernelWorkGroupInfo(_kernel_basis_dbasisdx, _device, CL_KERNEL_PRIVATE_MEM_SIZE,
                                            sizeof(private_usage), &private_usage, NULL);
            if(err < 0) {
                fatal("FEMFlowSolution:solutionComputer", "Error in getting kernel work-group size information.");
            };
            fprintf(stdout, "FEMFlowSolution:solutionComputer: < basis and basis derivatives > kernel maximum work group size: %zu\n", wg_size);
            fprintf(stdout, "FEMFlowSolution:solutionComputer: < basis and basis derivatives > kernel work group size best multiple: %zu\n", wg_multiple);
            fprintf(stdout, "FEMFlowSolution:solutionComputer: local memory used by the < basis and basis derivatives > kernel (KB): %llu\n", local_usage/1024);
            fprintf(stdout, "FEMFlowSolution:solutionComputer: private memory used by the < basis and basis derivatives > kernel (KB): %llu\n", private_usage/1024);
        }
    }
    
    if ([solution.solutionInfo[@"compute basis and basis derivatives in separate kernel"] boolValue] == YES) {
        cl_int err;
        cl_event basis_dBasis_event;
        cl_ulong kernel_sart, kernel_end, kernel_profile;
        size_t *local_work_size = NULL, val;
        double kernel_profile_fp;
        
        size_t global_work_size = mesh.numberOfBulkElements;
        
        if ([solution.solutionInfo[@"use gpu local memory"] boolValue] == YES) {
            if (solution.solutionInfo[@"parallel assembly work-group size"] != nil) {
                val = [solution.solutionInfo[@"parallel assembly work-group size"] intValue];
                local_work_size = &val;
            } else {
                fatal("FEMFlowSolution:solutionComputer", "Missing work-group size for the < basis and basis derivatives > kernel.");
            }
            
            int adjust = 0;
            if (solution.solutionInfo[@"adjust global work size to be a multiple of"] != nil) {
                adjust = [solution.solutionInfo[@"adjust global work size to be a multiple of"] intValue];
            }

            global_work_size = global_work_size + (adjust - (global_work_size & (adjust-1)));
        }
        
        // Queue up the basis and basis derivatives kernel
        err = clEnqueueNDRangeKernel(_cmd_queue, _kernel_basis_dbasisdx, 1, NULL, &global_work_size, local_work_size, 0, NULL, &basis_dBasis_event);
        if (err < 0) {
            fatal("FEMFlowSolution:solutionComputer", "Can't enqueue kernel. Error: ", err);
        }
        clFinish(_cmd_queue);
        
        // Do some profiling
        clGetEventProfilingInfo(basis_dBasis_event, CL_PROFILING_COMMAND_START, sizeof(kernel_sart), &kernel_sart, NULL);
        clGetEventProfilingInfo(basis_dBasis_event, CL_PROFILING_COMMAND_END, sizeof(kernel_end), &kernel_end, NULL);
        kernel_profile = (kernel_end - kernel_sart);
        kernel_profile_fp = (double)kernel_profile;
        fprintf(stdout, "FEMFlowSolution:solutionComputer: total time for < basis and basis derivatives > kernel execution on GPU (s): %f\n", kernel_profile_fp*1.0e-9);
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
        
        fprintf(stdout, "FEMFlowSolution:solutionComputer:\n");
        fprintf(stdout, "FEMFlowSolution:solutionComputer:\n");
        fprintf(stdout, "FEMFlowSolution:solutionComputer: -----------------------------------------------------------\n");
        fprintf(stdout, "FEMFlowSolution:solutionComputer: NAVIER-STOKES ITERATION %d.\n", iter);
        fprintf(stdout, "FEMFlowSolution:solutionComputer: -----------------------------------------------------------\n");
        fprintf(stdout, "FEMFlowSolution:solutionComputer:\n");
        fprintf(stdout, "FEMFlowSolution:solutionComputer: Starting Assembly...\n");

        [core initializeToZeroMatrix:solution.matrix forceVector:matContainers->RHS sizeForceVector:matContainers->sizeRHS model:model solution:solution];
        if (parallelAssembly == YES) {
            cl_int err;
            _mapped_matValues = enqueueDeviceMapBuffer(precision(), _cmd_queue, _matValues, CL_TRUE, CL_MAP_WRITE, 0, matContainers->sizeValues, 0, NULL, NULL, err);
            if (err < 0) {
                fatal("FEMFlowSolution:FEMFlowSolution_createMapGPUBuffersMesh", "Couldn't map _mapped_matValues.");
            }
            _mapped_matRHS = enqueueDeviceMapBuffer(precision(), _cmd_queue, _matRHS, CL_TRUE, CL_MAP_WRITE, 0, matContainers->sizeRHS, 0, NULL, NULL, err);
            if (err < 0) {
                fatal("FEMFlowSolution:FEMFlowSolution_createMapGPUBuffersMesh", "Couldn't map _mapped_matRHS.");
            }
            mapmemset(precision(), _mapped_matValues, matContainers->sizeValues);
            mapmemset(precision(), _mapped_matRHS, matContainers->sizeRHS);
            clEnqueueUnmapMemObject(_cmd_queue, _matValues, _mapped_matValues, 0, NULL, NULL);
            clEnqueueUnmapMemObject(_cmd_queue, _matRHS, _mapped_matRHS, 0, NULL, NULL);
        }
        
        bf_id = -1;
        body_id = -1;
        
        startAdvanceOutput((char *)[@"FEMFlowSolution" UTF8String], (char *)[@"Assembly:" UTF8String]);
        
        // Bulk elements
        if (parallelAssembly == YES) {
            int position;
            cl_int err;
            cl_event assembly_event;
            cl_ulong kernel_sart, kernel_end, kernel_profile;
            clock_t cpu_assembly_loop_start, cpu_assembly_loop_end;
            double real_assembly_loop_start, real_assembly_loop_end;
            double kernel_profile_fp;
            
            size_t *local_work_size = NULL, val;
            if ([solution.solutionInfo[@"use gpu local memory"] boolValue] == YES) {
                if (solution.solutionInfo[@"parallel assembly work-group size"] != nil) {
                    val = [solution.solutionInfo[@"parallel assembly work-group size"] intValue];
                    local_work_size = &val;
                } else {
                    fatal("FEMFlowSolution:solutionComputer", "Missing work-group size for the < assembly > kernel.");
                }
            }
            
            kernel_profile = 0;
            cpu_assembly_loop_start = clock();
            real_assembly_loop_start = mach_absolute_time();
            for (NSMutableArray *color in mesh.colors) {
                position = 0;
                for (i=0; i<[color[1] intValue]; i++) {
                    NSMutableArray *cc = mesh.colors[i];
                    position = position + [cc[0] intValue];
                }
                
                size_t global_work_size = [color[0] intValue];
                
                err = clSetKernelArg(_kernel_assembly, _positionInColorMapping, sizeof(cl_int), &position);
                if (err < 0) {
                    fatal("FEMFlowSolution:solutionComputer", "Error in setting kernel arguments.");
                }
                
                // Queue up the assembly kernel
                err = clEnqueueNDRangeKernel(_cmd_queue, _kernel_assembly, 1, NULL, &global_work_size, local_work_size, 0, NULL, &assembly_event);
                if (err < 0) {
                    fatal("FEMFlowSolution:solutionComputer", "Can't enqueue kernel. Error: ", err);
                }
                
                // Finish the calculation
                clFinish(_cmd_queue);
                if (err < 0) {
                    fatal("FEMFlowSolution:solutionComputer", "Can't finish kernel. Error: ", err);
                }
                
                // Do some profiling
                clGetEventProfilingInfo(assembly_event, CL_PROFILING_COMMAND_START, sizeof(kernel_sart), &kernel_sart, NULL);
                clGetEventProfilingInfo(assembly_event, CL_PROFILING_COMMAND_END, sizeof(kernel_end), &kernel_end, NULL);
                kernel_profile += (kernel_end - kernel_sart);
            }
            cpu_assembly_loop_end = clock();
            real_assembly_loop_end = mach_absolute_time();
            kernel_profile_fp = (double)kernel_profile;
            
            fprintf(stdout, "FEMFlowSolution:solutionComputer: total time for < assembly > kernel execution on GPU (s): %f\n", kernel_profile_fp*1.0e-9);
            fprintf(stdout, "FEMFlowSolution:solutionComputer: assembly execution on GPU: CPU time (s): %f\n", (cpu_assembly_loop_end-cpu_assembly_loop_start)/(double)CLOCKS_PER_SEC);
            fprintf(stdout, "FEMFlowSolution:solutionComputer: assembly execution on GPU: Real time (s): %f\n", machcore(real_assembly_loop_end,real_assembly_loop_start));
            
            // Read data from the device
            _mapped_matValues = enqueueDeviceMapBuffer(precision(), _cmd_queue, _matValues, CL_TRUE, CL_MAP_READ, 0, matContainers->sizeValues, 0, NULL, NULL, err);
            if (err < 0) {
                fatal("FEMFlowSolution:solutionComputer", "Couldn't map _mapped_matValues.");
            }
            _mapped_matRHS = enqueueDeviceMapBuffer(precision(), _cmd_queue, _matRHS, CL_TRUE, CL_MAP_READ, 0, matContainers->sizeRHS, 0, NULL, NULL, err);
            if (err < 0) {
                fatal("FEMFlowSolution:solutionComputer", "Couldn't map _mapped_matRHS.");
            }
            if (precision() == 0) {
                memcpy(matContainers->Values, _mapped_matValues, matContainers->sizeValues*sizeof(cl_double));
                memcpy(matContainers->RHS, _mapped_matRHS, matContainers->sizeRHS*sizeof(cl_double));
            } else {
                float *buff1 = _mapped_matValues;
                float *buff2 = _mapped_matRHS;
                for (i=0; i<matContainers->sizeValues; i++) {
                    matContainers->Values[i] = (double)buff1[i];
                }
                for (i=0; i<matContainers->sizeRHS; i++) {
                    matContainers->RHS[i] = (double)buff2[i];
                }
            }
            clEnqueueUnmapMemObject(_cmd_queue, _matValues, _mapped_matValues, 0, NULL, NULL);
            clEnqueueUnmapMemObject(_cmd_queue, _matRHS, _mapped_matRHS, 0, NULL, NULL);
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
                         newtonLinearization:_newtonLinearization
                                   transient:transient
                               stabilizeFlag:stabilizeFlag
                           divDiscretization:divDiscretization
                         gradPDiscretization:gradPDiscretization
                                   stabilize:stabilize
                                     bubbles:bubbles
                                   flowModel:flowModel
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
        fprintf(stdout, "FEMFlowSolution:solutionComputer: Assembly done.\n");
        
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
        fprintf(stdout, "FEMFlowSolution:solutionComputer: Dirichlet conditions done.\n");
        
        // Solve the system and check for convergence
        at = cputime() - at;
        st = cputime();
    
        if (nonLinearRelax != 1.0) memcpy(_pSolution, _flowSolution, flowContainers->sizeValues*sizeof(double));
        uNorm = [core findSolution:solution model:model backRorateNT:NULL];
        
        st = cputime() -  st;
        totat = totat + at;
        totst = totst + st;
        fprintf(stdout, "FEMFlowSolution:solutionComputer: iter: %d, Assembly (s): %f %f.\n", iter, at, totat);
        fprintf(stdout, "FEMFlowSolution:solutionComputer: iter: %d, Solve (s): %f %f.\n", iter, st, totst);
        
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
            if (solution.solutionInfo[@"nonlinear system relaxation before"] != nil) {
                relaxBefore = [solution.solutionInfo[@"nonlinear system relaxation before"] boolValue];
                found = YES;
            }
            if (found == NO || relaxBefore == YES) {
                [solution.solutionInfo setObject:@NO forKey:@"skip compute nonlinear change"];
                [core computeChange:solution model:model isSteadyState:NO nsize:&n values:_flowSolution values0:_pSolution sizeValues0:&flowContainers->sizeValues];
            }
        }
        if (parallelAssembly == YES) {
            cl_int err;
            if ([solution.solutionInfo[@"use element nodal data"] boolValue] == YES) {
                
                int nd;
                Element_t *element;
                int *indexes = intvec(0, solution.mesh.maxElementDofs-1);
                
                if (precision() == 0) {
                    _mapped_nodal_data = clEnqueueMapBuffer(_cmd_queue, _nodal_data, CL_TRUE, CL_MAP_WRITE, 0, sizeof(nodal_data_d)*mesh.numberOfBulkElements, 0, NULL, NULL, &err);
                    if (err < 0) {
                        fatal("FEMFlowSolution:FEMFlowSolution_createMapGPUBuffersMesh", "Couldn't map _mapped_nodal_data.");
                    }
                    nodal_data_d *ref = _mapped_nodal_data;
                    for (int t=0; t<mesh.numberOfBulkElements; t++) {
                        element = getActiveElementIMP(core, @selector(getActiveElement:solution:model:), t, solution, model);
                        nd = getElementDofsSolutionIMP(core, @selector(getElementDofsSolution:model:forElement:atIndexes:disableDiscontinuousGalerkin:), solution, model, element, indexes, NULL);
                        switch (_nsdofs) {
                            case 3:
                                for (int i=0; i<nd; i++) {
                                    ref[t].vx[i] = _flowSolution[_nsdofs*_flowPerm[indexes[i]]];
                                    ref[t].vy[i] = _flowSolution[_nsdofs*_flowPerm[indexes[i]]+1];
                                }
                            case 4:
                                for (int i=0; i<nd; i++) {
                                    ref[t].vx[i] = _flowSolution[_nsdofs*_flowPerm[indexes[i]]];
                                    ref[t].vy[i] = _flowSolution[_nsdofs*_flowPerm[indexes[i]]+1];
                                    ref[t].vz[i] = _flowSolution[_nsdofs*_flowPerm[indexes[i]]+2];
                                }
                        }
                    }
                } else {
                    _mapped_nodal_data = clEnqueueMapBuffer(_cmd_queue, _nodal_data, CL_TRUE, CL_MAP_WRITE, 0, sizeof(nodal_data_f)*mesh.numberOfBulkElements, 0, NULL, NULL, &err);
                    if (err < 0) {
                        fatal("FEMFlowSolution:FEMFlowSolution_createMapGPUBuffersMesh", "Couldn't map _mapped_nodal_data.");
                    }
                    nodal_data_f *ref = _mapped_nodal_data;
                    for (int t=0; t<mesh.numberOfBulkElements; t++) {
                        element = getActiveElementIMP(core, @selector(getActiveElement:solution:model:), t, solution, model);
                        nd = getElementDofsSolutionIMP(core, @selector(getElementDofsSolution:model:forElement:atIndexes:disableDiscontinuousGalerkin:), solution, model, element, indexes, NULL);
                        switch (_nsdofs) {
                            case 3:
                                for (int i=0; i<nd; i++) {
                                    ref[t].vx[i] = (float)_flowSolution[_nsdofs*_flowPerm[indexes[i]]];
                                    ref[t].vy[i] = (float)_flowSolution[_nsdofs*_flowPerm[indexes[i]]+1];
                                }
                            case 4:
                                for (int i=0; i<nd; i++) {
                                    ref[t].vx[i] = (float)_flowSolution[_nsdofs*_flowPerm[indexes[i]]];
                                    ref[t].vy[i] = (float)_flowSolution[_nsdofs*_flowPerm[indexes[i]]+1];
                                    ref[t].vz[i] = (float)_flowSolution[_nsdofs*_flowPerm[indexes[i]]+2];
                                }
                        }
                    }
                }
                
                clEnqueueUnmapMemObject(_cmd_queue, _nodal_data, _mapped_nodal_data, 0, NULL, NULL);
                free_ivector(indexes, 0, solution.mesh.maxElementDofs-1);
                
            } else {
                _mapped_varSolution = enqueueDeviceMapBuffer(precision(), _cmd_queue, _varSolution, CL_TRUE, CL_MAP_WRITE, 0, flowContainers->sizeValues, 0, NULL, NULL, err);
                if (err < 0) {
                    fatal("FEMFlowSolution:FEMFlowSolution_createMapGPUBuffersMesh", "Couldn't map _mapped_matValues.");
                }
                if (precision() == 0) {
                    memcpy(_mapped_varSolution, _flowSolution, flowContainers->sizeValues*sizeof(cl_double));
                } else {
                    float *buff = _mapped_varSolution;
                    for (i=0; i<flowContainers->sizeValues; i++) {
                        buff[i] = (float)_flowSolution[i];
                    }
                }
                clEnqueueUnmapMemObject(_cmd_queue, _varSolution, _mapped_varSolution, 0, NULL, NULL);
            }
        }
        
        relativeChange = solution.variable.nonLinChange;
        fprintf(stdout, "FEMFlowSolution:solutionComputer: result norm: %e.\n", solution.variable.norm);
        fprintf(stdout, "FEMFlowSolution:solutionComputer: relative change: %e.\n", relativeChange);
        
        if (relativeChange < newtonTol || iter > newtonIter) {
            _newtonLinearization = YES;
            if (parallelAssembly) {
                cl_int err;
                char *nonLinNewton = (char *)&_newtonLinearization;
                _mapped_gpuNewton = clEnqueueMapBuffer(_cmd_queue, _gpuNewtonLinear, CL_TRUE, CL_MAP_WRITE, 0, sizeof(cl_char), 0, NULL, NULL, &err);
                if (err < 0) {
                    fatal("FEMFlowSolution:FEMFlowSolution_createMapGPUBuffersMesh", "Couldn't map _mapped_gpuNewton.");
                }
                memcpy(_mapped_gpuNewton, nonLinNewton, sizeof(cl_char));
                clEnqueueUnmapMemObject(_cmd_queue, _gpuNewtonLinear, _mapped_gpuNewton, 0, NULL, NULL);
            }
        }
        if (relativeChange < nonLinearTol && iter < nonLinearIter) break;
        
        // If free surface in model, this will move the nodal points
        if (freeSurfaceFlag == YES) {
            FEMFreeSurface *freeSurface = [[FEMFreeSurface alloc] init];
            if (relativeChange < freeSTol || iter > freeSIter) computeFree = YES;
            if (computeFree == YES) {
                if (solution.solutionInfo[@"free surface relaxation factor"] != nil) {
                    relaxation = [solution.solutionInfo[@"free surface relaxation factor"] doubleValue];
                } else {
                    relaxation = 1.0;
                }
                found = NO;
                mbFlag = NO;
                if (solution.solutionInfo[@"internal move boundary"] != nil) {
                    mbFlag = [solution.solutionInfo[@"internal move boundary"] boolValue];
                    found = YES;
                }
                if (mbFlag == YES || found == NO) [freeSurface moveBoundaryModel:model integration:integration relax:relaxation];
            }
        }
    }
    
    [solution.solutionInfo setObject:@(nonLinearRelax) forKey:@"nonlinear system relaxation factor"];
    
    if ([solution.solutionInfo[@"adaptive mesh refinement"] boolValue] == YES) {
        // TODO: implement mesh refinement
    }
    
    [self FEMFlowSolution_checkCircleBoundaryModel:model];
    
    if (useLocalCoords == YES) {
        coordinatesSystems.coordinates = modelCoords;
        model.dimension = modelDim;
    }
    
    [integration deallocation:mesh];
    
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
