//
//  FEMCore.m
//  Saino
//
//  Created by Seddik hakime on 15/04/11.
//  Copyright 2011 ScienceSoul. All rights reserved.
//

#import "FEMCore.h"

#include <math.h>
#import <complex.h>

#import "FEMHUTIter.h"
#import "FEMIterativeMethods.h"
#import "FEMPrecondition.h"
#import "FEMParallelMPI.h"
#import "FEMElementDescription.h"
#import "FEMElementUtils.h"
#import "FEMNumericIntegration.h"
#import "FEMMeshUtils.h"
#import "FEMLinearAlgebra.h"
#import "FEMFlowSolution.h"
#import "FEMMagneticInductionSolution.h"
#import "FEMStressAnalysisSolution.h"
#import "FEMMeshUpdateSolution.h"
#import "FEMHeatSolution.h"
#import "FEMHeatSolution_OpenCL.h"
#import "GaussIntegration.h"
#import "Utils.h"

#ifdef TEST
    #import "FEMTest.h"
#endif

static const int ITER_BICGSTAB     =  320;
static const int ITER_TFQMR        =  330;
static const int ITER_CG           =  340;
static const int ITER_CGS          =  350;
static const int ITER_GMRES        =  360;
static const int ITER_BICGSTAB2    =  370;
static const int ITER_SGS          =  380;
static const int ITER_JACOBI       =  390;
static const int ITER_RICHARDSON   =  391;
static const int ITER_BICGSTABL    =  400;
static const int ITER_GCR          =  410;

static const int PRECOND_NONE      =  500;
static const int PRECOND_DIAGONAL  =  510;
static const int PRECOND_ILUN      =  520;
static const int PRECOND_ILUT      =  530;
static const int PRECOND_MG        =  540;
static const int PRECOND_BILUN     =  550;
static const int PRECOND_VANKA     =  560;

@interface FEMCore ()

// Solving linear systems and compute norms
-(void)FEMCore_iterCallType:(int)iterType solution:(FEMSolution * __nonnull)solution matrix:(FEMMatrix * __nonnull)matrix result:(double * __nonnull)x rhs:(double * __nonnull)b ipar:(int * __nonnull)ipar dpar:(double *)dpar work:(double * __nonnull * __nonnull)work pcondlMethod:(SEL __nonnull)pcondlMethod pcondrMethod:(SEL __nonnull)pcondrMethod matvecMethod:(SEL __nonnull)matvecMethod mstopMethod:(SEL __nonnull)mstopMethod;
-(void)FEMCore_rotateNTSystem:(double * __nonnull)vec nodeNumber:(int)nodeNumber model:(FEMModel * __nonnull)model;
-(void)FEMCore_backRotateNTSystem:(double * __nonnull)solution permutation:(int * __nonnull)perm ndofs:(int)ndofs model:(FEMModel * __nonnull)model;
-(double)FEMCore_computeNormInSolution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model size:(int)n values:(double * __nullable)values;
-(void)FEMCore_solveLinearSystemMatrix:(FEMMatrix * __nonnull)matrix rhs:(double * __nonnull)b result:(double * __nonnull)x norm:(double * __nonnull)norm dofs:(int)dofs solution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model bulkMatrix:(FEMMatrix * __nullable)bulkMatrix;

// Check passive element
-(BOOL)FEMCore_checkPassiveElement:(Element_t * __nonnull)element model:(FEMModel * __nonnull)model solution:(FEMSolution * __nonnull)solution;

// Check normal / tangential vector boundary conditions
-(void)FEMCore_checkNormalTangentialBoundaryModel:(FEMModel * __nonnull)model variableName:(NSString * __nonnull)variableName dimension:(int)dimension;

// Average boundary normals
-(void)FEMCore_averageBoundaryNormalsModel:(FEMModel * __nonnull)model variableName:(NSString * __nonnull)variableName dimension:(int)dimension;

// Rotate matrix
-(void)FEMCore_rotateMatrix:(double * __nonnull * __nonnull)matrix solution:(FEMSolution * __nonnull)solution vector:(double * __nonnull)vector size:(int)n dimension:(int)dim dofs:(int)dofs nodeIndexes:(int * __nonnull)nodeIndexes;

// Update global force
-(void)FEMCore_updateGlobalForceModel:(FEMModel * __nonnull)model solution:(FEMSolution * __nonnull)solution element:(Element_t * __nonnull)element forceVector:(double * __nonnull * __nonnull)forceVector forceVectorUpdateAtIndex:(int)index localForce:(double * __nonnull)localForce size:(int)n dofs:(int)dofs nodeIndexes:(int * __nonnull)nodeIndexes rotateNT:(BOOL * __nullable)rotateNT;

// Update force vector
-(void)FEMCore_finishAssemblyModel:(FEMModel * __nonnull)model solution:(FEMSolution * __nonnull)solution forceVector:(double * __nonnull)forceVector sizeForceVector:(int)n;

// Loads
-(void)FEMCore_setElementLoadsModel:(FEMModel * __nonnull)model solution:(FEMSolution * __nonnull)solution element:(Element_t * __nonnull)element values:(NSArray * __nonnull)values name:(NSString * __nonnull)name indexes:(int * __nonnull)indexes doneLoad:(BOOL * __nonnull)doneLoad size:(int)n dof:(int)dof ndofs:(int)ndofs diagonalScaling:(double * __nonnull)diagScaling;
-(void)FEMCore_setPointLoadsModel:(FEMModel * __nonnull)model solution:(FEMSolution * __nonnull)solution element:(Element_t * __nonnull)element values:(NSArray * __nonnull)values name:(NSString * __nonnull)name indexes:(int * __nonnull)indexes size:(int)n dof:(int)dof ndofs:(int)ndofs diagonalScaling:(double * __nonnull)diagScaling;

// Element and point values
-(void)FEMCore_setElementValues:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution forElementNumber:(int)elno numberOfNodes:(int)n atIndexes:(int * __nonnull)indexes forValues:(NSArray * __nonnull)values variableName:(NSMutableString * __nonnull)name orderOfDofs:(int)dof activeCondition:(BOOL)conditional conditionName:(NSString * __nonnull)condName permutationOffset:(int)offset diaginalScaling:(double * __nonnull)diagScaling offDiaginal:(BOOL)offDiaginal;
-(void)FEMCore_setPointValues:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution numberofNodes:(int)n atIndexes:(int * __nonnull)indexes forValues:(NSArray * __nonnull)values variableName:(NSMutableString * __nonnull)name orderOfDofs:(int)dof activeCondition:(BOOL)conditional conditionName:(NSString * __nonnull)condName permutationOffset:(int)offset diaginalScaling:(double * __nonnull)diagScaling offDiaginal:(BOOL)offDiaginal;

// Periodic
-(void)FEMCore_setPeriodicBoundariesPass1Model:(FEMModel * __nonnull)model solution:(FEMSolution * __nonnull)solution name:(NSMutableString * __nonnull)name dof:(int)dof this:(int)this done:(BOOL * __nonnull)done diaginalScaling:(double * __nonnull)diagScaling;
-(void)FEMCore_setPeriodicBoundariesPass2Model:(FEMModel * __nonnull)model solution:(FEMSolution * __nonnull)solution name:(NSMutableString * __nonnull)name dof:(int)dof this:(int)this done:(BOOL * __nonnull)done diaginalScaling:(double * __nonnull)diagScaling;

// Limiter
-(void)FEMCore_determineSoftLimiterInSolution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model;
-(void)FEMCore_setLimiterValues:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution numberofNodes:(int)n atIndexes:(int * __nonnull)indexes forValues:(NSArray * __nonnull)values  orderOfDofs:(int)dof limitActive:(bool * __nonnull)limitActive conditionName:(NSString * __nonnull)condName permutationOffset:(int)offset diaginalScaling:(double * __nonnull)diagScaling offDiaginal:(BOOL)offDiaginal;

// Rows equilibration
-(void)FEMCore_rowEquilibrationMatrix:(FEMMatrix * __nonnull)matrix vector:(double * __nonnull)f parallel:(BOOL)parallel;
-(void)FEMCore_reverseRowEquilibrationMatrix:(FEMMatrix * __nonnull)matrix vector:(double * __nonnull)f;

// Method for coupled solution used by solveEquationsModel:timeStep:transientSimulation:coupledMinIteration:coupleMaxIteration:steadyStateReached:realTimeStep:
-(void)FEMCore_solveCoupledModel:(FEMModel * __nonnull)model timeStep:(double)dt coupledMinIteration:(int)coupledMinIter coupleMaxIteration:(int)coupleMaxIter transientSimulation:(BOOL)transient scanning:(BOOL)scanning doneThis:(BOOL * __nonnull)doneThis afterConverged:(BOOL * __nonnull)afterConverged steadyIt:(double * __nonnull)steadyIt;

// Single solution
-(void)FEMCore_singleSolution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model timeStep:(double)dt transientSimulation:(BOOL)transient;

//
-(void)FEMCore_updateExportedVariablesSolution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model;

@end

@implementation FEMCore {
    
    double * __nullable * __nullable _kernStiff;
    double * __nullable _kernWork;                         // kernStiff(maxElementDofs)(maxElementDofs), kernWork(maxElementDofs)
    int * __nullable _g_Ind;
    int * __nullable _l_Ind;                               // g_Ind(maxElementDofs), l_Ind(maxElementDofs)
    int _size1kernStiff;
    int _size2kernStiff;
    int _sizekernWork;
    int _size_g_Ind;
    int _size_l_Ind;
    
    int * __nonnull * __nonnull _lineEM;
    int * __nonnull * __nonnull _triangleEM;
    int * __nonnull * __nonnull _quadEM;
    int * __nonnull * __nonnull _tetraEM;
    int * __nonnull * __nonnull _prismEM;
    int * __nonnull * __nonnull _wedgeEM;
    int * __nonnull * __nonnull _brickEM;
    
    BOOL _initialized[8];
    
    int _k1, _n1;
    double * __nullable _saveValues;
    double * __nullable * __nullable _damp;
    double * __nullable * __nullable _stiff;
    double * __nullable * __nullable _mass;
    double * __nullable * __nullable _x;

    id <SainoSolutionsComputer> __nullable _instance;
}

#pragma mark Private methods...

#pragma mark Solve linear systems and norms

-(void)FEMCore_iterCallType:(int)iterType solution:(FEMSolution * __nonnull)solution matrix:(FEMMatrix * __nonnull)matrix result:(double * __nonnull)x rhs:(double * __nonnull)b ipar:(int * __nonnull)ipar dpar:(double *)dpar work:(double * __nonnull * __nonnull)work pcondlMethod:(SEL __nonnull)pcondlMethod pcondrMethod:(SEL __nonnull)pcondrMethod matvecMethod:(SEL __nonnull)matvecMethod mstopMethod:(SEL __nonnull)mstopMethod {

    if (pcondrMethod == 0) {
        pcondrMethod = @selector(CRSPCondDummyMatrix:afterPrecondition:rightHandSide:info:);
    }
    
    if (matrix.isComplexMatrix == NO) {
        // Solvers from the FEMHUTiter class:
        // ----------------------------------
        if (iterType == ITER_BICGSTAB) { // Solve with BI-CGSTAB
            
            FEMHUTIter *iterSolver = [[FEMHUTIter alloc] init];
            [iterSolver dbicgstabSolveInSolution:solution matrix:matrix ndim:ipar[2] result:x rhs:b ipar:ipar dpar:dpar pcondlMethod:pcondlMethod pcondrMethod:pcondrMethod matvecMethod:matvecMethod mstopMethod:mstopMethod];
        }
        else if (iterType == ITER_BICGSTAB2) { // Solve with BI-CGSTAB2
            
            FEMHUTIter *iterSolver = [[FEMHUTIter alloc] init];
            [iterSolver dbicgstab2SolveInSolution:solution matrix:matrix ndim:ipar[2] result:x rhs:b ipar:ipar dpar:dpar pcondlMethod:pcondlMethod pcondrMethod:pcondrMethod matvecMethod:matvecMethod mstopMethod:mstopMethod];
        }
        else if (iterType == ITER_TFQMR) { // Solve with TFQMR
            
            FEMHUTIter *iterSolver = [[FEMHUTIter alloc] init];
            [iterSolver dtfqmrSolveInSolution:solution matrix:matrix ndim:ipar[2] result:x rhs:b ipar:ipar dpar:dpar pcondlMethod:pcondlMethod pcondrMethod:pcondrMethod matvecMethod:matvecMethod mstopMethod:mstopMethod];
        }
        else if (iterType == ITER_CG) { // Solve with CG
            
            FEMHUTIter *iterSolver = [[FEMHUTIter alloc] init];
            [iterSolver dcgSolveInSolution:solution matrix:matrix ndim:ipar[2] result:x rhs:b ipar:ipar dpar:dpar pcondlMethod:pcondlMethod pcondrMethod:pcondrMethod matvecMethod:matvecMethod mstopMethod:mstopMethod];
        }
        else if (iterType == ITER_CGS) { // Solve with CGS
            
            FEMHUTIter *iterSolver = [[FEMHUTIter alloc] init];
            [iterSolver dcgsSolveInSolution:solution matrix:matrix ndim:ipar[2] result:x rhs:b ipar:ipar dpar:dpar pcondlMethod:pcondlMethod pcondrMethod:pcondrMethod matvecMethod:matvecMethod mstopMethod:mstopMethod];
        }
        else if (iterType == ITER_GMRES) { // Solve with GMRES
            
            FEMHUTIter *iterSolver = [[FEMHUTIter alloc] init];
            [iterSolver dgmresSolveInSolution:solution matrix:matrix ndim:ipar[2] wrkdim:ipar[3] result:x rhs:b ipar:ipar dpar:dpar work:work pcondlMethod:pcondlMethod pcondrMethod:pcondrMethod matvecMethod:matvecMethod mstopMethod:mstopMethod];
            
        }
        // Solvers from the FEMIterativeMethods class:
        // -------------------------------------------
        else if (iterType == ITER_SGS) { // Solve with SGS
            
            FEMIterativeMethods *iterSolver = [[FEMIterativeMethods alloc] init];
            [iterSolver dsgsSolveMatrix:matrix ndim:ipar[2] result:x rhs:b ipar:ipar dpar:dpar pcondlMethod:pcondlMethod pcondrMethod:pcondrMethod matvecMethod:matvecMethod mstopMethod:mstopMethod];
        }
        else if (iterType == ITER_JACOBI) { // Solve with Jacobi
            
            FEMIterativeMethods *iterSolver = [[FEMIterativeMethods alloc] init];
            [iterSolver djacobiSolveMatrix:matrix ndim:ipar[2] result:x rhs:b ipar:ipar dpar:dpar pcondlMethod:pcondlMethod pcondrMethod:pcondrMethod matvecMethod:matvecMethod mstopMethod:mstopMethod];
        }
        else if (iterType == ITER_BICGSTABL) { // Solve with BI-CGSTAB(l)
            
            FEMIterativeMethods *iterSolver = [[FEMIterativeMethods alloc] init];
            [iterSolver dbicgstablSolveMatrix:matrix ndim:ipar[2] result:x rhs:b ipar:ipar dpar:dpar pcondlMethod:pcondlMethod pcondrMethod:pcondrMethod matvecMethod:matvecMethod mstopMethod:mstopMethod];
        }
        else if (iterType == ITER_GCR) { // Solve with GCR
            
            FEMIterativeMethods *iterSolver = [[FEMIterativeMethods alloc] init];
            [iterSolver dgcrSolveMatrix:matrix ndim:ipar[2] result:x rhs:b ipar:ipar dpar:dpar pcondlMethod:pcondlMethod pcondrMethod:pcondrMethod matvecMethod:matvecMethod mstopMethod:mstopMethod];
        }
         else if (iterType == ITER_RICHARDSON) { // Solve with Richardson
             
             FEMIterativeMethods *iterSolver = [[FEMIterativeMethods alloc] init];
             [iterSolver drichardsonSolveMatrix:matrix ndim:ipar[2] result:x rhs:b ipar:ipar dpar:dpar pcondlMethod:pcondlMethod pcondrMethod:pcondrMethod matvecMethod:matvecMethod mstopMethod:mstopMethod];
        }
    } else {
        //TODO: implements solvers for complex matrix
        fatal("FEMCore:FEMCore_iterCallType:", "linear solvers for complex matrix not implemented yet.");
    }
}

-(void)FEMCore_rotateNTSystem:(double * __nonnull)vec nodeNumber:(int)nodeNumber model:(FEMModel * __nonnull)model {
    
    int i, k, dim;
    double bu, bv, bw, rm[3][3];
    
    if (self.normalTangentialNumberOfNodes <= 0) return;
    
    dim = model.dimension;
    
    k = self.boundaryReorder[nodeNumber];
    if (k < 0) return;
    
    if (dim < 3) {
        bu = vec[0];
        bv = vec[1];
        vec[0] = self.boundaryNormals[k][0]*bu + self.boundaryNormals[k][1]*bv;
        vec[1] = -self.boundaryNormals[k][1]*bu + self.boundaryNormals[k][0]*bv;
    } else {
        
        bu = vec[0];
        bv = vec[1];
        bw = vec[2];
        
        for (i=0; i<3; i++) {
            rm[i][0] = self.boundaryNormals[k][i];
            rm[i][1] = self.boundaryTangent1[k][i];
            rm[i][2] = self.boundaryTangent2[k][i];
        }
        
        vec[0] = rm[0][0]*bu + rm[1][0]*bv + rm[2][0]*bw;
        vec[1] = rm[0][1]*bu + rm[1][1]*bv + rm[2][1]*bw;
        vec[2] = rm[0][2]*bu + rm[1][2]*bv + rm[2][2]*bw;
    }
}

-(void)FEMCore_backRotateNTSystem:(double * __nonnull)solution permutation:(int * __nonnull)perm ndofs:(int)ndofs model:(FEMModel * __nonnull)model {
    
    int i, j, k, l, dim;
    double bu, bv, bw, rm[3][3];
    
    dim = model.dimension;
    
    if (self.normalTangentialNumberOfNodes <= 0 || ndofs < dim) return;
    
    for (i=0; i<self.sizeBoundaryReorder; i++) {
        k = self.boundaryReorder[i];
        if (k < 0) continue;
        j = perm[i];
        if (j < 0) continue;
        
        if (dim < 3) {
            
            bu = solution[(ndofs*j)];
            bv = solution[(ndofs*j+1)];
            
            solution[(ndofs*j)] = self.boundaryNormals[k][0]*bu - self.boundaryNormals[k][1]*bv;
            solution[(ndofs*j+1)] = self.boundaryNormals[k][1]*bu + self.boundaryNormals[k][0]*bv;
        } else {
                        
            bu = solution[(ndofs*j)];
            bv = solution[(ndofs*j+1)];
            bw = solution[(ndofs*j+2)];
            
            for (l=0; l<3; l++) {
                rm[0][l] = self.boundaryNormals[k][l];
                rm[1][l] = self.boundaryTangent1[k][l];
                rm[2][l] = self.boundaryTangent2[k][l];
            }
            
            solution[(ndofs*j)] = rm[0][0]*bu + rm[1][0]*bv + rm[2][0]*bw;
            solution[(ndofs*j+1)] = rm[0][1]*bu + rm[1][1]*bv + rm[2][1]*bw;
            solution[(ndofs*j+2)] = rm[0][2]*bu + rm[1][2]*bv + rm[2][2]*bw;
        }
    }
}

/*******************************************************
    Computes the norm related to a solution vector.
*******************************************************/
-(double)FEMCore_computeNormInSolution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model size:(int)n values:(double * __nullable)values {
    
    int normDim, normDofs, dofs, i, j, k, l, nn, permStart, totn;
    int *iPerm;
    double norm, nscale, sum;
    double *x = NULL, *y;
    FEMParallelMPI *parallelUtil;
    variableArraysContainer *varContainers = NULL;
    listBuffer normComponents = { NULL, NULL, NULL, NULL, 0, 0, 0};
    BOOL componentsAllocated, found;
    
    parallelUtil = [[FEMParallelMPI alloc] init];
    varContainers = solution.variable.getContainers;
    
    if (values != NULL) {
        x = values;
    } else {
        x = varContainers->Values;
    }
    
    if ((solution.solutionInfo)[@"nonlinear system norm degree"] != nil) {
        normDim = [(solution.solutionInfo)[@"nonlinear system norm degree"] intValue];
    } else {
        normDim = 2;
    }
    
    dofs = solution.variable.dofs;
    
    componentsAllocated = NO;
    FEMListUtilities *listUtilities = [FEMListUtilities sharedListUtilities];
    found = [listUtilities listGetIntegerArray:model inArray:solution.valuesList forVariable:@"nonlinear system norm components" buffer:&normComponents]; //Those are indexes
    if (found == YES) {
        normDofs = normComponents.m;
        componentsAllocated = YES;
    } else {
        if ((solution.solutionInfo)[@"nonlinear system norm dofs"] != nil) {
            normDofs = [(solution.solutionInfo)[@"nonlinear system norm dofs"] intValue];
            normComponents.ivector = intvec(0, normDofs-1);
            componentsAllocated = YES;
            for (i=0; i<normDofs; i++) {
                normComponents.ivector[i] = i;
            }
        } else {
            normDofs = dofs;
        }
    }
    
    nn = n;
    
    if ((solution.solutionInfo)[@"norm permutation"] != nil) { // This is an index
        permStart = [(solution.solutionInfo)[@"norm permutation"] intValue];
        iPerm = intvec(0, varContainers->sizePerm-1);
        memset( iPerm, -1, varContainers->sizePerm*sizeof(int) );
        nn = 0;
        for (i=permStart; i<varContainers->sizePerm; i++) {
            if (varContainers->Perm[i] >= 0) {
                iPerm[nn] = varContainers->Perm[i];
                nn++;
            }
        }
        y = doublevec(0, nn-1);
        for (i=0; i<nn; i++) {
            y[i] = x[iPerm[i]];
        }
        x = y;
        free_ivector(iPerm, 0, varContainers->sizePerm-1);
    }
    
    
    totn = [parallelUtil parallelReductionOfValue:(1.0*n) operArg:NULL];
    nscale = normDofs * totn/(1.0*dofs);
    
    if (normDofs < dofs) {
        norm = 0.0;
        double maxVal;
        double buffer[nn/dofs];
        switch (normDim) {
            case 0:
                for (i=0; i<normDofs; i++) {
                    j = normComponents.ivector[i];
                    memset( buffer, 0.0, sizeof(buffer) );
                    vDSP_vabsD(x+j, dofs, buffer, 1, nn/dofs);
                    vDSP_maxvD(buffer, 1, &maxVal, nn/dofs);
                    norm = max(norm, maxVal);
                }
                l = 2;
                norm = [parallelUtil parallelReductionOfValue:norm operArg:&l];
                break;
            case 1:
                for (i=0; i<normDofs; i++) {
                    j = normComponents.ivector[i];
                    memset( buffer, 0.0, sizeof(buffer) );
                    vDSP_vabsD(x+j, dofs, buffer, 1, nn/dofs);
                    vDSP_sveD(buffer, 1, &sum, nn/dofs);
                    norm = norm + sum;
                }
                norm = [parallelUtil parallelReductionOfValue:norm operArg:NULL] / nscale;
                break;
            case 2:
                for (i=0; i<normDofs; i++) {
                    j = normComponents.ivector[i];
                    vDSP_svesqD(x+j, dofs, &sum, nn/dofs);
                    norm = norm + sum;
                }
                norm = sqrt([parallelUtil parallelReductionOfValue:norm operArg:NULL] / nscale);
                break;
            default:
                for (i=0; i<normDofs; i++) {
                    j = normComponents.ivector[i];
                    memset( buffer, 0.0, sizeof(buffer) );
                    l = 0.0;
                    for (k=j; k<nn; k+=dofs) {
                        buffer[l] = pow(x[k], normDim);
                        l++;
                    }
                    vDSP_sveD(buffer, 1, &sum, nn/dofs);
                    norm = norm + sum;
                }
                norm = pow( ([parallelUtil parallelReductionOfValue:norm operArg:NULL] / nscale), (1.0/normDim) );
                break;
        }
    } else {
        double buffer[nn];
        memset( buffer, 0.0, sizeof(buffer) );
        switch (normDim) {
            case 0:
                vDSP_vabsD(x, 1, buffer, 1, nn);
                vDSP_maxvD(buffer, 1, &norm, nn);
                l = 2;
                norm = [parallelUtil parallelReductionOfValue:norm operArg:&l];
                break;
            case 1:
                vDSP_vabsD(x, 1, buffer, 1, nn);
                vDSP_sveD(buffer, 1, &sum, nn);
                norm = [parallelUtil parallelReductionOfValue:sum operArg:NULL] / nscale;
                break;
            case 2:
                vDSP_svesqD(x, 1, &sum, nn);
                norm = sqrt([parallelUtil parallelReductionOfValue:sum operArg:NULL] / nscale);
                break;
            default:
                for (i=0; i<nn; i++) {
                    buffer[i] = pow(x[i], normDim);
                }
                vDSP_sveD(buffer, 1, &sum, nn);
                norm = pow( ([parallelUtil parallelReductionOfValue:sum operArg:NULL] / nscale), (1.0/normDim));
                break;
        }
    }
    
    if (componentsAllocated == YES) {
        free_ivector(normComponents.ivector, 0, normDofs-1);
    }
    
    return norm;
}

/*********************************************************************************
    Solves a linear system and also calls the necessary preconditioning methods
*********************************************************************************/
-(void)FEMCore_solveLinearSystemMatrix:(FEMMatrix * __nonnull)matrix rhs:(double * __nonnull)b result:(double * __nonnull)x norm:(double * __nonnull)norm dofs:(int)dofs solution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model bulkMatrix:(FEMMatrix * __nullable)bulkMatrix {
    
    int i, j, k, l, ii, m, n, dof;
    double bnorm, sum, energy, *tempVector, *saveValues;
    NSString *method;
    NSMutableString *name;
    FEMVariable *nodalLoads;
    FEMMatrix *matrixAid, *projector;
    matrixArraysContainer *matContainers = NULL, *matAidContainers = NULL,
                          *projectorContainers = NULL;
    variableArraysContainer *varContainers = NULL, *nodalContainers = NULL;
    FEMUtilities *utilities;
    FEMListUtilities *listUtilities;
    FEMParallelMPI *parallelUtil;
    BOOL scaleSystem, eigenAnalysis, harmonicAnalysis, backRotation, applyLimiter, skipZeroRhs, applyRowEquilibration=NO, parallel,
         rhsScaling, found;
    
    parallelUtil = [[FEMParallelMPI alloc] init];
    
    n = matrix.numberOfRows;
    matContainers = matrix.getContainers;
    varContainers = solution.variable.getContainers;
    
    if ((solution.solutionInfo)[@"back rotate n-t solution"] != nil) {
        backRotation = [(solution.solutionInfo)[@"back rotate n-t solution"] boolValue];
    } else {
        backRotation = YES;
    }
    
    if (matrix.isLumped == YES && solution.timeOrder == 1)  {
        if ((solution.solutionInfo)[@"time stepping method"] != nil) {
            method = [NSString stringWithString:(solution.solutionInfo)[@"time stepping method"]];
            if ([method isEqualToString:@"runge-kutta"] == YES || [method isEqualToString:@"explicit euler"] == YES) {
                for (i=0; i<n; i++) {
                    if (fabs(matContainers->Values[matContainers->Diag[i]]) > 0.0)
                        x[i] = b[i] / matContainers->Values[matContainers->Diag[i]];
                }
                if (backRotation == YES) [self FEMCore_backRotateNTSystem:x permutation:varContainers->Perm ndofs:dofs model:model];
                *norm = [self FEMCore_computeNormInSolution:solution model:model size:n values:x];
                return;
            }
        }
    }
    
    // These definitions are needed if changing the iterative solver on the fly
    solution.multigridSolution = ([(solution.solutionInfo)[@"linear system solver"] isEqualToString:@"multigrid"]) ? YES : NO;
    solution.multiGridTotal = max(solution.multiGridTotal, [(solution.solutionInfo)[@"mg levels"] intValue]);
    solution.multiGridTotal = max(solution.multiGridTotal, [(solution.solutionInfo)[@"multigrid levels"] intValue]);
    solution.multigridLevel = solution.multiGridTotal;
    
    if ((solution.solutionInfo)[@"linear system scaling"] != nil) {
        scaleSystem = [(solution.solutionInfo)[@"linear system scaling"] boolValue];
    } else {
        scaleSystem = YES;
    }
    
    eigenAnalysis = (solution.nOfEigenValues  > 0 && [(solution.solutionInfo)[@"eigen analysis"] boolValue] == YES) ? YES: NO;
    harmonicAnalysis = (solution.nOfEigenValues  > 0 && [(solution.solutionInfo)[@"harmonic analysis"] boolValue] == YES) ? YES: NO;

    applyLimiter = NO;
    if ((solution.solutionInfo)[@"apply limiter"] != nil) {
        applyLimiter = [(solution.solutionInfo)[@"apply limiter"] boolValue];
    }
    
    skipZeroRhs = NO;
    if ((solution.solutionInfo)[@"skip zero rhs test"] != nil) {
        skipZeroRhs = [(solution.solutionInfo)[@"skip zero rhs test"] boolValue];
    }
    
    if ( !(harmonicAnalysis || eigenAnalysis || applyLimiter || skipZeroRhs)) {
        vDSP_svesqD(b, 1, &sum, n);
        bnorm = [parallelUtil parallelReductionOfValue:sqrt(sum) operArg:NULL];
        if (bnorm <= DBL_MIN) {
            NSLog(@"FEMCore:FEMCore_solveLinearSystemMatrix: solution trivially zero.\n");
            memset( x, 0.0, varContainers->sizeValues*sizeof(double) );
            return;
        }
    }
    
    if (solution.multigridLevel == -1) return;
        
    // If solving harmonic anaysis go there:
    if (harmonicAnalysis == YES) {
        
        // TODO: implement the harmonic pass
        
    }
    
    // If solving eigen system go there:
    if (eigenAnalysis) {
        
        // TODO:implement the eigen pass

    }
    
    // If whether b=0 sice then equation Ax=b has only the trivial solution, x=0.
    // In case of a limiter one still may need to check the limiter for contact
    vDSP_svesqD(b, 1, &sum, n);
    bnorm = [parallelUtil parallelReductionOfValue:sqrt(sum) operArg:NULL];
    if (bnorm <= DBL_MIN && skipZeroRhs == NO) {
        NSLog(@"FEMCore:FEMCore_solveLinearSystemMatrix: solution trivially zero.\n");
        memset( x, 0.0, varContainers->sizeValues*sizeof(double) );
        if (applyLimiter == YES) {
            [self FEMCore_determineSoftLimiterInSolution:solution model:model];
        }
        return;
    }
    
    // Convert rhs & initial value to the scaled system
    if (scaleSystem == YES) {
        if ((solution.solutionInfo)[@"linear system row equilibration"] != nil) {
            applyRowEquilibration = [(solution.solutionInfo)[@"linear system row equilibration"] boolValue];
        } else applyRowEquilibration = NO;
        
        if (applyRowEquilibration == YES) {
            // TODO: Test wheter we are in a parallel run, for now it's never the case
            parallel = NO;
            [self FEMCore_rowEquilibrationMatrix:matrix vector:b parallel:parallel];
        } else {
            rhsScaling = (bnorm != 0.0) ? YES : NO;
            [self scaleLinearSystem:solution matrix:matrix rhs:b result:x diagScaling:NULL applyScaling:NULL rhsScaling:&rhsScaling];
        }
    }
    
    if ((solution.solutionInfo)[@"linear system solver"] != nil) {
        method = [NSString stringWithString:(solution.solutionInfo)[@"linear system solver"]];
    } else {
        method = @"iterative";
    }

    if ([method isEqualToString:@"iterative"] == YES) {
        
        [self iterativeSolveMatrix:matrix result:x rhs:b dimensions:NULL solution:solution];
    } 
    else if ([method isEqualToString:@"multigrid"] == YES) {
        // TODO: Need to be implemented
    }
    else if ([method isEqualToString:@"feti"] == YES) {
        // TODO: Need to be implemented
    }
    else if ([method isEqualToString:@"block"]) {
        // TODO: Need to be implemented
    }
    else if ([method isEqualToString:@"mortar"]) {
        // TODO: Need to be implemented
    }
    else { // Direct solver
           // TODO: Need to be implemented
    }
    
    // TODO: Add suport for parallel run in which case parallel linear solvers must
    // be called
    
    if (scaleSystem == YES) {
        if (applyRowEquilibration == YES) {
            [self FEMCore_reverseRowEquilibrationMatrix:matrix vector:b];
        } else {
            [self backScaleLinearSystem:solution matrix:matrix rhs:b result:x diagScaling:NULL sizeOFDiagScaling:NULL];
        }
    }
    
    utilities = [[FEMUtilities alloc] init];
    nodalLoads = [utilities getVariableFrom:solution.mesh.variables model:model name:[[solution.variable canonicalizeName] stringByAppendingString:@" loads"] onlySearch:NULL maskName:NULL info:&found];
    
    matrixAid = solution.matrix;
    // If bulkMatrix is given then it's not nil
    if (bulkMatrix != nil) {
        matrixAid = bulkMatrix;
    }
    
    matAidContainers = matrixAid.getContainers;
    if (nodalLoads != nil && matAidContainers->BulkValues != NULL) {
        tempVector = doublevec(0, matrixAid.numberOfRows-1);
        saveValues = matAidContainers->Values;
        matAidContainers->Values = matAidContainers->BulkValues;
        
        //TODO: Add support for parallel run
        [self matrixVectorMultplyInMatrix:matrixAid multiplyVector:x resultVector:tempVector];
        
        if ([(solution.solutionInfo)[@"calculate energy norm"] boolValue] == YES) {
            energy = 0.0;
            for (i=0; i<matrixAid.numberOfRows; i++) {
                //TODO: Add suppor for pararrel run
                energy = energy + x[i]*tempVector[i];
            }
            energy = [parallelUtil parallelReductionOfValue:energy operArg:NULL];
            [solution.solutionInfo setObject:@(energy) forKey:@"energy norm"];
            name = [NSMutableString stringWithString:@"res: "];
            [name appendString:[solution.variable canonicalizeName]];
            [name appendString:@" energy norm"];
            listUtilities = [FEMListUtilities sharedListUtilities];
            [listUtilities addConstRealInClassList:model.simulation theVariable:name withValue:&energy orUsingBlock:nil string:nil];
            
            NSLog(@"FEMCore:FEMCore_solveLinearSystemMatrix: energy norm: %f.\n", energy);
        }
        
        matAidContainers->Values = saveValues;
        //TODO: add support for parallel run
        for (i=0; i<matrixAid.numberOfRows; i++) {
            tempVector[i] = tempVector[i] - matAidContainers->BulkRHS[i];
        }
        
        for (FEMBoundaryCondition *boundaryCondition in model.boundaryConditions) {
            projector = boundaryCondition.pMatrix;
            if (projector != nil) {
                projectorContainers = projector.getContainers;
                for (dof=0; dof<dofs; dof++) {
                    for (i=0; i<projector.numberOfRows; i++) {
                        ii = projectorContainers->InvPerm[i];
                        k = varContainers->Perm[ii];
                        if (k < 0) continue;
                        k = dofs * k + dof;
                        tempVector[k] = 0.0;
                        
                        for (l=projectorContainers->Rows[i]; l<=projectorContainers->Rows[i+1]-1; l++) {
                            if (projectorContainers->Cols[l] < 0) continue;
                            m = varContainers->Perm[projectorContainers->Cols[l]];
                            if (m >= 0) {
                                m = dofs * m + dof;
                                tempVector[k] = tempVector[k] + projectorContainers->Values[l]*tempVector[m];
                            }
                        }
                    }
                }
            }
        }
        
        nodalContainers = nodalLoads.getContainers;
        for (i=0; i<nodalContainers->sizePerm; i++) {
            if (nodalContainers->Perm[i] >= 0 && varContainers->Perm[i] >= 0) {
                for (j=0; j<dofs; j++) {
                    nodalContainers->Values[dofs*nodalContainers->Perm[i]+j] = tempVector[dofs*varContainers->Perm[i]+j];
                }
            }
        }
        free_dvector(tempVector, 0, matrixAid.numberOfRows-1);
        
        if (backRotation == YES) [self FEMCore_backRotateNTSystem:nodalContainers->Values permutation:nodalContainers->Perm ndofs:dofs model:model];
        
    }
    
    if (backRotation == YES) [self FEMCore_backRotateNTSystem:x permutation:varContainers->Perm ndofs:dofs model:model];
    
    // Compute the change of the solution with different methods
    [self computeChange:solution model:model isSteadyState:NO nsize:&n values:x values0:NULL sizeValues0:NULL];
    *norm = solution.variable.norm;
    
    // Create soft limiters to be later applied by the dirichlet conditions in the next round.
    // Within apply a hard limiter after the set is determined.
    if (applyLimiter == YES) {
        [self FEMCore_determineSoftLimiterInSolution:solution model:model];
    }
    
    solution.variable.primaryMesh = solution.mesh;
    [self invalidateVariableInTopMesh:model.meshes primaryMesh:solution.mesh name:[solution.variable canonicalizeName] model:model];
    
    if (nodalLoads != nil) {
        nodalLoads.primaryMesh = solution.mesh;
        [self invalidateVariableInTopMesh:model.meshes primaryMesh:solution.mesh name:[nodalLoads canonicalizeName] model:model];
    }
    
    // In order to ba able to change the preconditioners or solutions, the old
    //  matrix structures must be deallocated on request
    if ([(solution.solutionInfo)[@"linear system preconditioning deallocate"] boolValue] == YES) {
        // ILU predonditioning
        if (matContainers->ILUValues != NULL) {
            if (matContainers->sizeILUValues != matContainers->sizeValues) {
                free_ivector(matContainers->ILUCols, 0, matContainers->sizeILUCols-1);
                free_ivector(matContainers->ILURows, 0, matContainers->sizeILURows-1);
                free_ivector(matContainers->ILUDiag, 0, matContainers->sizeILUDiag-1);
                free_dvector(matContainers->ILUValues, 0, matContainers->sizeILUValues-1);
                matContainers->ILUCols = NULL; matContainers->ILURows = NULL;
                matContainers->ILUDiag = NULL; matContainers->ILUValues = NULL;
            }
        }
        
        // Multigrid solver / preconditioner
        if (solution.multigridLevel > 0) {
            matrixAid = solution.matrix;
            if (matrixAid.parent != nil) {
                while (matrixAid.parent != nil) {
                    matrixAid = matrixAid.parent;
                }
                while (matrixAid.child != nil) {
                    matrixAid = matrixAid.child;
                    if (matrixAid.parent != nil) matrixAid.parent = nil;
                    if (matrixAid.ematrix != nil) matrixAid.ematrix = nil;
                }
            }
        }
    }
}

#pragma mark Element info

-(BOOL)FEMCore_checkPassiveElement:(Element_t * __nonnull)element model:(FEMModel * __nonnull)model solution:(FEMSolution * __nonnull)solution {
    
    int body_id, bf_id, nbNodes;
    listBuffer passive = { NULL, NULL, NULL, NULL, 0, 0, 0};
    NSMutableString *passName;
    NSString *str;
    BOOL isPassive, found;
    
    FEMListUtilities *listUtil;
    FEMBodyForce *bodyForceAtBodyID;
    
    isPassive = NO;
    
    body_id = element->BodyID-1;
    if (body_id < 0) return isPassive; // body_id == -1 for boundary elements
    
    listUtil = [FEMListUtilities sharedListUtilities];
    
    nbNodes = element->Type.NumberOfNodes;
    
    if ((model.bodies)[body_id] == nil) {
        return isPassive;
    } else {
        bf_id = [(model.bodies)[body_id][@"body force"] intValue];
    }
    
    passName = [NSMutableString stringWithString:[solution.variable canonicalizeName]];
    str = @" passive";
    [passName appendString:str];
    
    // Returns a bodyForce object at index bf_id
    bodyForceAtBodyID = (model.bodyForces)[bf_id-1];
    
    found = [listUtil listGetReal:model inArray:bodyForceAtBodyID.valuesList forVariable:passName numberOfNodes:nbNodes indexes:element->NodeIndexes buffer:&passive minValue:NULL maxValue:NULL];
    if (found == YES) {
        if ( count(passive.vector, '>', 0.0, nbNodes) > count(passive.vector, '<', 0.0, nbNodes) ) {
            isPassive = YES;
        }
    }
    
    if (passive.vector != NULL) free_dvector(passive.vector, 0, passive.m-1);
    
    return isPassive;
}

#pragma mark Check if normal / tangential

/*************************************************************************************
    Check if normal / tangential vector boundary conditions present and allocate
    space normals and if in 3D two tangent direction vectors.
*************************************************************************************/
-(void)FEMCore_checkNormalTangentialBoundaryModel:(FEMModel * __nonnull)model variableName:(NSString * __nonnull)variableName dimension:(int)dimension {
    
    int j, k, n, t;
    BOOL conditional, found, gotIt;
    FEMMesh *mesh;
    FEMListUtilities *listUtilities;
    Element_t *element = NULL;
    listBuffer condition = { NULL, NULL, NULL, NULL, 0, 0, 0};
    
    listUtilities = [FEMListUtilities sharedListUtilities];
    
    // Need an early initialization to average normals across partitions
    // TODO: add support for parallel runs
    
    self.normalTangentialNumberOfNodes = 0;
    
    found = NO;
    for (FEMBoundaryCondition *boundaryCondition in model.boundaryConditions) {
        if ([listUtilities listGetLogical:model inArray:boundaryCondition.valuesList forVariable:variableName info:&gotIt] == YES) {
            found = [listUtilities listGetLogical:model inArray:boundaryCondition.valuesList
                                      forVariable:[variableName stringByAppendingString:@" rotate"]
                                             info:&gotIt];
            if (gotIt == NO) found = YES;
            if (found == YES) break;
        }
    }
    if (found == NO) return;
    
    mesh = (FEMMesh *)model.mesh;
    n = mesh.numberOfNodes;
    element = mesh.getElements;
    
    if (self.boundaryReorder == NULL) {
        self.boundaryReorder = intvec(0, n-1);
        self.sizeBoundaryReorder = n;
    } else if (self.sizeBoundaryReorder < n) {
        free_ivector(self.boundaryReorder, 0, self.sizeBoundaryReorder-1);
        self.boundaryReorder = intvec(0, n-1);
        self.sizeBoundaryReorder = n;
    }
    memset( self.boundaryReorder, -1, n*sizeof(int) );
    
    for (t=mesh.numberOfBulkElements; t<mesh.numberOfBulkElements+mesh.numberOfBoundaryElements; t++) {
        if (element[t].Type.ElementCode == 101) continue;
        
        n = element[t].Type.NumberOfNodes;
        for (FEMBoundaryCondition *boundaryCondition in model.boundaryConditions) {
            if (element[t].BoundaryInfo->Constraint == boundaryCondition.tag) {
                if ([listUtilities listGetLogical:model inArray:boundaryCondition.valuesList forVariable:variableName info:&gotIt] == YES) {
                    found = [listUtilities listGetLogical:model inArray:boundaryCondition.valuesList forVariable:[variableName stringByAppendingString:@" rotate"] info:&gotIt];
                    if (found == YES || gotIt == NO) {
                        conditional = [listUtilities listGetReal:model inArray:boundaryCondition.valuesList forVariable:[variableName stringByAppendingString:@" condition"] numberOfNodes:n indexes:element[t].NodeIndexes buffer:&condition minValue:NULL maxValue:NULL];
                        for (j=0; j<n; j++) {
                            if (conditional == YES && condition.vector[j] < 0.0) continue;
                            
                            k = element[t].NodeIndexes[j];
                            if (self.boundaryReorder[k] < 0) {
                                self.boundaryReorder[k] = self.normalTangentialNumberOfNodes;
                                self.normalTangentialNumberOfNodes++;
                            }
                        }
                    }
                }
            }
        }
    }
    if (condition.vector != NULL) free_dvector(condition.vector, 0, condition.m-1);
    
    // TODO: Add support for parallel runs
    
    if (self.normalTangentialNumberOfNodes == 0) {
        // No ops
    } else {
        if (self.boundaryNormals != NULL) {
            free_dmatrix(self.boundaryNormals, 0, self.size1boundaryNormals-1, 0, self.size2boundaryNormals-1);
            free_dmatrix(self.boundaryTangent1, 0, self.size1boundaryTangent1-1, 0, self.size2boundaryTangent1-1);
            free_dmatrix(self.boundaryTangent2, 0, self.size1boundaryTangent2-1, 0, self.size2boundaryTangent2-1);
        }
        self.ntElement = intmatrix(0, self.normalTangentialNumberOfNodes-1, 0, 2);
        self.size1NtElement = self.normalTangentialNumberOfNodes;
        self.size2NtElement = 3;
        
        self.ntZeroingDone = boolmatrix(0, self.normalTangentialNumberOfNodes-1, 0, 2);
        self.size1NtZeroingDone = self.normalTangentialNumberOfNodes;
        self.size2NtZeroingDone = 3;
        
        self.boundaryNormals = doublematrix(0, self.normalTangentialNumberOfNodes-1, 0, 2);
        self.size1boundaryNormals = self.normalTangentialNumberOfNodes;
        self.size2boundaryNormals = 3;
        
        self.boundaryTangent1 = doublematrix(0, self.normalTangentialNumberOfNodes-1, 0, 2);
        self.size1boundaryTangent1 = self.normalTangentialNumberOfNodes;
        self.size2boundaryTangent1 = 3;
        
        self.boundaryTangent2 = doublematrix(0, self.normalTangentialNumberOfNodes-1, 0, 2);
        self.size1boundaryTangent2 = self.normalTangentialNumberOfNodes;
        self.size2boundaryTangent2 = 3;
        
        memset( *self.boundaryNormals, 0.0, (self.normalTangentialNumberOfNodes*3)*sizeof(double) );
        memset( *self.boundaryTangent1, 0.0, (self.normalTangentialNumberOfNodes*3)*sizeof(double) );
        memset( *self.boundaryTangent2, 0.0, (self.normalTangentialNumberOfNodes*3)*sizeof(double) );
    }
}


#pragma mark Average normals

/**************************************************************************************
    Average boundary normals for nodes. The average boundary normals may be beneficial
    as they provide more continuous definition of normal curved boundaries
**************************************************************************************/
-(void)FEMCore_averageBoundaryNormalsModel:(FEMModel * __nonnull)model variableName:(NSString * __nonnull)variableName dimension:(int)dimension {
    
    int i, j, k, l, m, n, t, tt;
    double bu, bv, detJ, *lrnm, *nrm, s, sum;
    double *buffer1, *buffer2, *buffer3;
    bool *done;
    BOOL any, check, conditional, found, gotIt, massConsistent, stat;
    Element_t *elements;
    Nodes_t *elementNodes, *nodes;
    FEMMesh *mesh;
    FEMVariable *nrmVar = nil;
    FEMUtilities *utilities;
    FEMListUtilities *listUtilities;
    FEMElementDescription *elementDescription;
    FEMElementUtils *elementUtils;
    FEMNumericIntegration *integration;
    variableArraysContainer *nrmContainers = NULL;
    matrixArraysContainer *projectorContainers = NULL;
    listBuffer condition = { NULL, NULL, NULL, NULL, 0, 0, 0};
    listBuffer rot = { NULL, NULL, NULL, NULL, 0, 0, 0};
    GaussIntegrationPoints *IP = NULL;
    
    utilities = [[FEMUtilities alloc] init];
    listUtilities = [FEMListUtilities sharedListUtilities];
    
    elementNodes = (Nodes_t*)malloc(sizeof(Nodes_t));
    elementNodes->x = doublevec(0, model.maxElementNodes-1);
    elementNodes->y = doublevec(0, model.maxElementNodes-1);
    elementNodes->z = doublevec(0, model.maxElementNodes-1);
    
    mesh = (FEMMesh *)model.mesh;
    nrmVar = [utilities getVariableFrom:mesh.variables model:model name:@"normals" onlySearch:NULL maskName:nil info:&found];
    
    if (nrmVar != nil) {
        nrmContainers = nrmVar.getContainers;
        if (self.normalTangentialNumberOfNodes > 0) {
            memset( *self.boundaryNormals, 0.0, (self.size1boundaryNormals*self.size2boundaryNormals)*sizeof(double) );
            for (i=0; i<model.numberOfNodes; i++) {
                k = self.boundaryReorder[i];
                if (k >= 0) {
                    for (l=0; l<nrmVar.dofs; l++) {
                        self.boundaryNormals[k][l] = nrmContainers->Values[nrmVar.dofs*nrmContainers->Perm[i]+l];
                    }
                }
            }
        }
    } else {
        // Compute sum of element-wise normals for nodes on boundaries
        elementDescription = [FEMElementDescription sharedElementDescription];
        integration = [[FEMNumericIntegration alloc] init];
        if ([integration allocation:mesh] == NO) fatal("FEMCore:FEMCore_localMatrix", "Allocation error in FEMNumericIntegration.");
        elements = mesh.getElements;
        nodes = mesh.getNodes;
        nrm = doublevec(0, 2);
        lrnm = doublevec(0, 2);
        
        if (_normalTangentialNumberOfNodes > 0) {
            memset( *self.boundaryNormals, 0.0, (self.size1boundaryNormals*self.size2boundaryNormals)*sizeof(double) );
            for (t=model.numberOfBulkElements; t<model.numberOfBulkElements+model.numberOfBoundaryElements; t++) {
                if (elements[t].Type.ElementCode < 200) continue;
                
                n = elements[t].Type.NumberOfNodes;
                for (i=0; i<n; i++) {
                    elementNodes->x[i] = nodes->x[elements[t].NodeIndexes[i]];
                    elementNodes->y[i] = nodes->y[elements[t].NodeIndexes[i]];
                    elementNodes->z[i] = nodes->z[elements[t].NodeIndexes[i]];
                }
                
                IP = GaussQuadrature(&elements[t], NULL, NULL);
                for (FEMBoundaryCondition *boundaryCondition in model.boundaryConditions) {
                    if (elements[t].BoundaryInfo->Constraint == boundaryCondition.tag) {
                        if ([listUtilities listGetLogical:model inArray:boundaryCondition.valuesList forVariable:variableName info:&gotIt] == YES) {
                            found = [listUtilities listGetLogical:model inArray:boundaryCondition.valuesList forVariable:[variableName stringByAppendingString:@" rotate"] info:&gotIt];
                            if (found == YES || gotIt == NO) {
                                massConsistent = [listUtilities listGetLogical:model inArray:boundaryCondition.valuesList forVariable:@"mass consistent normals" info:&gotIt];
                                
                                conditional = [listUtilities listGetReal:model inArray:boundaryCondition.valuesList forVariable:[variableName stringByAppendingString:@" condition"] numberOfNodes:n indexes:elements[t].NodeIndexes buffer:&condition minValue:NULL maxValue:NULL];
                                for (j=0; j<n; j++) {
                                    if (conditional == YES && condition.vector[j] < 0.0) continue;
                                    k = self.boundaryReorder[elements[t].NodeIndexes[j]];
                                    if (k >= 0) {
                                        memset( nrm, 0.0, 3*sizeof(double) );
                                        if (massConsistent == YES) {
                                            for (tt=0; tt<IP->n; tt++) {
                                                stat = [integration setBasisForElement:&elements[t] elementNodes:elementNodes inMesh:mesh firstEvaluationPoint:IP->u[tt] secondEvaluationPoint:IP->v[tt] thirdEvaluationPoint:IP->w[tt] withBubbles:NO basisDegree:NULL];
                                                stat = [integration setMetricDeterminantForElement:&elements[t] elementNodes:elementNodes inMesh:mesh firstEvaluationPoint:IP->u[tt] secondEvaluationPoint:IP->v[tt] thirdEvaluationPoint:IP->w[tt]];
                                                detJ = integration.metricDeterminant;
                                                check = YES;
                                                [elementDescription normalVectorForBDElement:&elements[t] boundaryNodes:elementNodes mesh:mesh paraU:&IP->u[tt] paraV:&IP->v[tt] check:&check normals:lrnm];
                                                for (i=0; i<3; i++) {
                                                    nrm[i] = nrm[i] + IP->s[tt] * lrnm[i] * detJ * integration.basis[j];
                                                }
                                            }
                                        } else {
                                            bu = elements[k].Type.NodeU[j];
                                            bv = elements[k].Type.NodeV[j];
                                            check = YES;
                                            [elementDescription normalVectorForBDElement:&elements[t] boundaryNodes:elementNodes mesh:mesh paraU:&bu paraV:&bv check:&check normals:nrm];
                                        }
                                        for (i=0; i<self.size2boundaryNormals; i++) {
                                            self.boundaryNormals[k][i] = self.boundaryNormals[k][i] + nrm[i];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (condition.vector != NULL) free_dvector(condition.vector, 0, condition.m-1);
            
            for (FEMBoundaryCondition *boundaryCondition in model.boundaryConditions) {
                if (boundaryCondition.pMatrix == nil) continue;
                
                projectorContainers = boundaryCondition.pMatrix.getContainers;
                found = [listUtilities listGetConstRealArray:model inArray:boundaryCondition.valuesList forVariable:@"periodic bc rotate" buffer:&rot];
                if (found == YES && rot.matrix != NULL) {
                    any = NO;
                    for (i=0; i<rot.m; i++) {
                        for (j=0; j<rot.n; j++) {
                            if (rot.matrix[i][j] != 0.0) {
                                any = YES;
                                break;
                            }
                        }
                    }
                    if (any == YES) {
                        done = boolvec(0, self.size1boundaryNormals-1);
                        memset( done, false, self.size1boundaryNormals*sizeof(bool) );
                        for (i=0; i<boundaryCondition.pMatrix.numberOfRows; i++) {
                            k = self.boundaryReorder[projectorContainers->InvPerm[i]];
                            if (k < 0) continue;
                            for (l=projectorContainers->Rows[i]; l<=projectorContainers->Rows[i+1]-1; l++) {
                                if (projectorContainers->Cols[l] < 0) continue;
                                m = self.boundaryReorder[projectorContainers->Cols[l]];
                                if (m >= 0) {
                                    if (done[m] == false) {
                                        done[m] = true;
                                        for (j=0; j<self.size2boundaryNormals; j++) {
                                            self.boundaryNormals[m][j] = -self.boundaryNormals[m][j];
                                        }
                                    }
                                }
                            }
                        }
                        free_bvector(done, 0, self.size1boundaryNormals-1);
                        continue;
                    }
                }
                
                for (i=0; i<boundaryCondition.pMatrix.numberOfRows; i++) {
                    k = self.boundaryReorder[projectorContainers->InvPerm[i]];
                    if (k < 0) continue;
                    for (l=projectorContainers->Rows[i]; l<=projectorContainers->Rows[i+1]-1; l++) {
                        if (projectorContainers->Cols[l] < 0) continue;
                        m = self.boundaryReorder[projectorContainers->Cols[l]];
                        if (m >= 0) {
                            for (j=0; j<self.size2boundaryNormals; j++) {
                                self.boundaryNormals[m][j] = 0.0;
                            }
                        }
                    }
                }
            }
            
            for (FEMBoundaryCondition *boundaryCondition in model.boundaryConditions) {
                if (boundaryCondition.pMatrix == nil) continue;
                
                projectorContainers = boundaryCondition.pMatrix.getContainers;
                found = [listUtilities listGetConstRealArray:model inArray:boundaryCondition.valuesList forVariable:@"periodic bc rotate" buffer:&rot];
                if (found == YES && rot.matrix != NULL) {
                    any = NO;
                    for (i=0; i<rot.m; i++) {
                        for (j=0; j<rot.n; j++) {
                            if (rot.matrix[i][j] != 0.0) {
                                any = YES;
                                break;
                            }
                        }
                    }
                    if (any == YES) continue;
                }
                
                for (i=0; i<boundaryCondition.pMatrix.numberOfRows; i++) {
                    k = self.boundaryReorder[projectorContainers->InvPerm[i]];
                    if (k < 0) continue;
                    for (l=projectorContainers->Rows[i]; l<=projectorContainers->Rows[i+1]-1; l++) {
                        if (projectorContainers->Cols[l] < 0) continue;
                        m = self.boundaryReorder[projectorContainers->Cols[l]];
                        if (m >= 0) {
                            for (j=0; j<self.size2boundaryNormals; j++) {
                                self.boundaryNormals[m][j] = self.boundaryNormals[m][j] +
                                           projectorContainers->Values[l] * self.boundaryNormals[k][j];
                            }
                        }
                    }
                }
            }
            if (rot.vector != NULL) free_dvector(rot.vector, 0, rot.m-1);
        }
        
        // TODO: add support for parallel runs
        
        [integration deallocation:mesh];
        free_dvector(nrm, 0, 2);
        free_dvector(lrnm, 0, 2);
    }
    
    // Normalize
    if (self.normalTangentialNumberOfNodes > 0) {
        buffer1 = doublevec(0, self.size2boundaryNormals-1);
        buffer2 = doublevec(0, self.size2boundaryTangent1-1);
        buffer3 = doublevec(0, self.size2boundaryTangent2-1);
        elementUtils = [[FEMElementUtils alloc] init];
        for (i=0; i<model.numberOfNodes; i++) {
            k = self.boundaryReorder[i];
            if (k >= 0) {
                sum = 0.0;
                for (j=0; j<self.size2boundaryNormals; j++) {
                    sum = sum + pow(self.boundaryNormals[k][j], 2.0);
                }
                s = sqrt(sum);
                if (s != 0.0) {
                    for (j=0; j<self.size2boundaryNormals; j++) {
                        self.boundaryNormals[k][j] = self.boundaryNormals[k][j] / s;
                    }
                }
                if (model.dimension > 2) {
                    for (j=0; j<self.size2boundaryNormals; j++) {
                        buffer1[j] = self.boundaryNormals[k][j];
                    }
                    [elementUtils tangentDirectionsForNormal:buffer1 tangent1:buffer2 tangent2:buffer3];
                    for (j=0; j<self.size2boundaryNormals; j++) {
                        self.boundaryTangent1[k][j] = buffer2[j];
                        self.boundaryTangent2[k][j] = buffer3[j];
                    }
                }
            }
        }
        free_dvector(buffer1, 0, self.size2boundaryNormals-1);
        free_dvector(buffer2, 0, self.size2boundaryTangent1-1);
        free_dvector(buffer3, 0, self.size2boundaryTangent2-1);
    }
    free_dvector(elementNodes->x, 0, model.maxElementNodes-1);
    free_dvector(elementNodes->y, 0, model.maxElementNodes-1);
    free_dvector(elementNodes->z, 0, model.maxElementNodes-1);
    free(elementNodes);
}


#pragma mark Manipulate matrix

-(void)FEMCore_rotateMatrix:(double * __nonnull * __nonnull)matrix solution:(FEMSolution * __nonnull)solution vector:(double * __nonnull)vector size:(int)n dimension:(int)dim dofs:(int)dofs nodeIndexes:(int * __nonnull)nodeIndexes {
    
    int i, j, k, l;
    double s, r[n*dofs][n*dofs], q[n*dofs][n*dofs], n1[3], t1[3], t2[3];
    
    for (i=0; i<n; i++) {
        
        if (nodeIndexes[i] < 0 || nodeIndexes[i]+1 >= _size1boundaryNormals) continue;
        
        memset( *r, 0.0, ((n*dofs)*(n*dofs))*sizeof(double) );
        for (j=0; j<n*dofs; j++) {
            r[j][j] = 1.0;
        }
        
        for (j=0; j<_size2boundaryNormals; j++) {
            n1[j] = _boundaryNormals[nodeIndexes[i]][j];
        }
        
        switch (dim) {
            case 2:
                r[dofs*i][dofs*i] = n1[0];
                r[dofs*i][dofs*i+1] = n1[1];
                
                r[dofs*i+1][dofs*i] = -n1[1];
                r[dofs*i+1][dofs*i+1] = n1[0];
                break;
            case 3:
                for (j=0; j<_size2boundaryTangent1; j++) {
                    t1[j] = _boundaryTangent1[nodeIndexes[i]][j];
                }
                for (j=0; j<_size2boundaryTangent2; j++) {
                    t2[j] = _boundaryTangent2[nodeIndexes[i]][j];
                }
                
                r[dofs*i][dofs*i] = n1[0];
                r[dofs*i][dofs*i+1] = n1[1];
                r[dofs*i][dofs*i+2] = n1[2];
                
                r[dofs*i+1][dofs*i] = t1[0];
                r[dofs*i+1][dofs*i+1] = t1[1];
                r[dofs*i+1][dofs*i+2] = t1[2];
                
                r[dofs*i+2][dofs*i] = t2[0];
                r[dofs*i+2][dofs*i+1] = t2[1];
                r[dofs*i+2][dofs*i+2] = t2[2];
                break;
        }
        
        for (j=0; j<n*dofs; j++) {
            for (k=0; k<n*dofs; k++) {
                s = 0.0;
                for (l=0; l<n*dofs; l++) {
                    s = s + r[j][l] * matrix[l][k];
                }
                q[j][k] = s;
            }
        }
        
        for (j=0; j<n*dofs; j++) {
            for (k=0; k<n*dofs; k++) {
                s = 0.0;
                for (l=0; l<n*dofs; l++) {
                    s = s + q[j][l] * r[k][l];
                }
                matrix[j][k] = s;
            }
        }
        
        for (j=0; j<n*dofs; j++) {
            s = 0.0;
            for (k=0; k<n*dofs; k++) {
                s = s + r[j][k] * vector[k];
            }
            q[j][0] = s;
        }
        
        for (j=0; j<n*dofs; j++) {
            vector[j] = q[j][0];
        }
    }
}

-(void)FEMCore_updateGlobalForceModel:(FEMModel * __nonnull)model solution:(FEMSolution * __nonnull)solution element:(Element_t * __nonnull)element forceVector:(double * __nonnull * __nonnull)forceVector forceVectorUpdateAtIndex:(int)index localForce:(double * __nonnull)localForce size:(int)n dofs:(int)dofs nodeIndexes:(int * __nonnull)nodeIndexes rotateNT:(BOOL * __nullable)rotateNT {
    
    int i, j, k, dim, indexes[n];
    BOOL rotate;
    double **localStiffMatrix;
    
    // Check if this element has been defined as passive
    if ([self FEMCore_checkPassiveElement:element model:model solution:solution] == YES) return;
        
    rotate = YES;
    if (rotateNT != NULL) rotate = *rotateNT;
    
    if (rotate == YES && _normalTangentialNumberOfNodes > 0) {
        localStiffMatrix = doublematrix(0, (n*dofs)-1, 0, (n*dofs)-1);
        dim = model.dimension;
        memset( indexes, -1, sizeof(indexes) );
        for (i=0; i<element->Type.NumberOfNodes; i++) {
            indexes[i] = _boundaryReorder[element->NodeIndexes[i]];
        }
        [self FEMCore_rotateMatrix:localStiffMatrix solution:solution vector:localForce size:n dimension:dim dofs:dofs nodeIndexes:indexes];
        free_dmatrix(localStiffMatrix, 0, (n*dofs)-1, 0, (n*dofs)-1);
    }
    
    for (i=0; i<n; i++) {
        if (nodeIndexes[i] >= 0) {
            for (j=0; j<dofs; j++) {
                k = dofs*nodeIndexes[i] + j;
                forceVector[k][index] = forceVector[k][index] + localForce[dofs*i+j];
            }
        }
        
    }
}

#pragma mark Update force

/***********************************************************************************************
    Update force vector after all other assembly steps but before setting Dirichlet conditions.
    Required only for time dependent simulations.
***********************************************************************************************/
-(void)FEMCore_finishAssemblyModel:(FEMModel * __nonnull)model solution:(FEMSolution * __nonnull)solution forceVector:(double * __nonnull)forceVector sizeForceVector:(int)n {
    
    int i, order;
    NSString *method, *simulation;
    BOOL found;
    matrixArraysContainer *matContainers = NULL;
    
    if (solution.matrix.format == MATRIX_LIST) {
        // TODO: implement this case for list matrix
    }

    if ([(solution.solutionInfo)[@"linear system fct"] boolValue] == YES) {
        if (solution.variable.dofs == 1) {
            FEMMatrixCRS *matrixCRS = [[FEMMatrixCRS alloc] init];
            [matrixCRS fctlLowOrderInSolution:solution orMatrix:nil];
        } else {
            fatal("FEMCore:FEMCore_finishAssemblyModel", "FCT scheme implemented only for one dof.");
        }
    }
    
    FEMListUtilities *listUtilities = [FEMListUtilities sharedListUtilities];
    [listUtilities listGetString:model inArray:model.simulation.valuesList forVariable:@"simulation type" info:&found];
    if ([simulation isEqualToString:@"transient"] == YES) {
        method = (solution.solutionInfo)[@"time stepping method"];
        order = min(solution.doneTime, solution.order);
        
        if (order <= 0 || solution.timeOrder != 1 || [method isEqualToString:@"bdf"] == YES) return;
        
        if (solution.beta != 0.0) {
            matContainers = solution.matrix.getContainers;
            for (i=0; i<n; i++) {
                forceVector[i] = forceVector[i] + (solution.beta - 1.0) * matContainers->Force[i][0] + (1.0 - solution.beta) * matContainers->Force[i][1];
            }
        }
    }
}

#pragma mark Loads

-(void)FEMCore_setElementLoadsModel:(FEMModel * __nonnull)model solution:(FEMSolution * __nonnull)solution element:(Element_t * __nonnull)element values:(NSArray * __nonnull)values name:(NSString * __nonnull)name indexes:(int * __nonnull)indexes doneLoad:(BOOL * __nonnull)doneLoad size:(int)n dof:(int)dof ndofs:(int)ndofs diagonalScaling:(double * __nonnull)diagScaling {
    
    int j, k, l, k1;
    listBuffer work = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    listBuffer workA = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    FEMListUtilities *listUtil;
    matrixArraysContainer *matContainers = NULL;
    variableArraysContainer *varContainers = NULL;
    BOOL stat;
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
        
    listUtil = [FEMListUtilities sharedListUtilities];
    
    if (dof >= 0) {
        stat = [listUtil listGetReal:model inArray:values forVariable:name numberOfNodes:n indexes:indexes buffer:&work minValue:NULL maxValue:NULL];
        if (stat == NO) {
            NSMutableString *componentName = [NSMutableString stringWithString:name];
            [componentName appendString:@" dofs"];
            stat = [listUtil listGetReal:model inArray:values forVariable:componentName numberOfNodes:n indexes:indexes buffer:&work minValue:NULL maxValue:NULL];
        }
    } else {
        stat = [listUtil listGetRealArray:model inArray:values forVariable:name numberOfNodes:n indexes:indexes buffer:&workA];
    }
    
    if (stat == YES) {
    
        for (j=0; j<n; j++) {
            k = varContainers->Perm[indexes[j]];
            
            if (doneLoad[k] == YES) continue;
            doneLoad[k] = YES;
            
            if (k >= 0) {
                if (dof >= 0) {
                    k = ndofs * k + dof;
                    matContainers->RHS[k] = matContainers->RHS[k] + work.vector[j] * diagScaling[k];
                } else {
                    for (l=0; l<min(ndofs, workA.m); l++) {
                        k1 = ndofs * k + l;
                        matContainers->RHS[k1] = matContainers->RHS[k1] + workA.tensor[l][0][j] * diagScaling[k1];
                    }
                }
            }
        }
    
    }
    
    if (workA.tensor != NULL) {
        free_d3tensor(workA.tensor, 0,  workA.m-1, 0,  workA.n-1, 0, workA.p-1);
    }
    if (work.vector != NULL) {
        free_dvector(work.vector, 0, work.m-1);
    }
}

-(void)FEMCore_setPointLoadsModel:(FEMModel * __nonnull)model solution:(FEMSolution * __nonnull)solution element:(Element_t * __nonnull)element values:(NSArray * __nonnull)values name:(NSString * __nonnull)name indexes:(int * __nonnull)indexes size:(int)n dof:(int)dof ndofs:(int)ndofs diagonalScaling:(double * __nonnull)diagScaling {
    
    int j, k, l, k1;
    listBuffer work = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    listBuffer workA = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    FEMListUtilities *listUtil;
    matrixArraysContainer *matContainers = NULL;
    variableArraysContainer *varContainers = NULL;
    BOOL stat;
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
    
    listUtil = [FEMListUtilities sharedListUtilities];
    
    if (dof >= 0) {
        stat = [listUtil listGetReal:model inArray:values forVariable:name numberOfNodes:n indexes:indexes buffer:&work minValue:NULL maxValue:NULL];
    } else {
        stat = [listUtil listGetRealArray:model inArray:values forVariable:name numberOfNodes:n indexes:indexes buffer:&workA];
    }
    
    if (stat == YES) {
        
        for (j=0; j<n; j++) {
            
            if (indexes[j] >= varContainers->sizePerm || indexes[j] < 0) {
                NSLog(@"FEMCore:FEMCore_setPointLoadsModel: invalid node number.\n");
                continue;
            }
            
            k = varContainers->Perm[indexes[j]];
            if (k >= 0) {
                if (dof >= 0) {
                    k = ndofs * k + dof;
                    matContainers->RHS[k] = matContainers->RHS[k] + work.vector[j] * diagScaling[k];
                } else {
                    for (l=0; l<min(ndofs, workA.m); l++) {
                        k1 = ndofs * k + l;
                        matContainers->RHS[k1] = matContainers->RHS[k1] + workA.tensor[l][0][j] * diagScaling[k1];
                    }
                }
            }
        }
    }
    
    if (workA.tensor != NULL) {
        free_d3tensor(workA.tensor, 0,  workA.m-1, 0,  workA.n-1, 0, workA.p-1);
    }
    if (work.vector != NULL) {
        free_dvector(work.vector, 0, work.m-1);
    }
}

#pragma mark Element and point values

-(void)FEMCore_setElementValues:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution forElementNumber:(int)elno numberOfNodes:(int)n atIndexes:(int * __nonnull)indexes forValues:(NSArray * __nonnull)values variableName:(NSMutableString * __nonnull)name orderOfDofs:(int)dof activeCondition:(BOOL)conditional conditionName:(NSString * __nonnull)condName permutationOffset:(int)offset diaginalScaling:(double * __nonnull)diagScaling offDiaginal:(BOOL)offDiaginal {
    
    int i, j, k, l, m, dim, k1;
    BOOL checkNT, stat, all;
    double *rotvec;
    listBuffer condition = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    listBuffer work = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    listBuffer workA = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    FEMMatrixCRS *crsMatrix;
    FEMMatrixBand *bandMatrix;
    matrixArraysContainer *matContainers = NULL;
    variableArraysContainer *varContainers = NULL;
    solutionArraysContainer *solContainers = NULL;
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
    solContainers = solution.getContainers;
    
    rotvec = doublevec(0, 2);
    dim = model.dimension;
    
    NSString *str = [name stringByAppendingString:@" dofs"];
    FEMListUtilities *listUtil = [FEMListUtilities sharedListUtilities];
    if (dof >= 0) {
        stat = [listUtil listGetReal:model inArray:values forVariable:name numberOfNodes:n indexes:indexes buffer:&work minValue:NULL maxValue:NULL];
        if (stat == NO) {
            stat = [listUtil listGetReal:model inArray:values forVariable:str numberOfNodes:n indexes:indexes buffer:&work minValue:NULL maxValue:NULL];
        }
    } else {
        stat = [listUtil listGetRealArray:model inArray:values forVariable:name numberOfNodes:n indexes:indexes buffer:&workA];
    }
    
    if (stat == YES) {
        
        if (conditional == YES) {
            stat = [listUtil listGetReal:model inArray:values forVariable:condName numberOfNodes:n indexes:indexes buffer:&condition minValue:NULL maxValue:NULL];
            conditional = (conditional == YES && stat == YES) ? YES : NO;
        }
        
        // Check for nodes belonging to n-t boundary getting set by other BCs
        checkNT = NO;
        if (self.normalTangentialNumberOfNodes > 0 && dof >= 0) {
            checkNT = YES;
            all = YES;
            for (i=0; i<n; i++) {
                if (self.boundaryReorder[indexes[i]] < 0) {
                    continue;
                } else {
                    all = NO;
                    break;
                }
            }
            if (all == YES) checkNT = NO;
            if ([listUtil listGetLogical:model inArray:values forVariable:solution.normalTangentialName info:&stat] == YES) checkNT = NO;
        }
        
        for (j=0; j<n; j++) {
            
            @autoreleasepool {
                if (conditional == YES && condition.vector[j] < 0.0) continue;
                k = varContainers->Perm[indexes[j]];
                
                if (k >= 0) {
                    if (dof >= 0) {
                        m = 0;
                        if (self.normalTangentialNumberOfNodes > 0) m = self.boundaryReorder[indexes[j]];
                        if (m >= 0 && checkNT == YES) {
                            memset( rotvec, 0.0, 3*sizeof(double) );
                            rotvec[dof] = 1.0;
                            [self FEMCore_rotateNTSystem:rotvec nodeNumber:indexes[j] model:model];
                            for (k=0; k<dim; k++) {
                                if (fabs(rotvec[k]) > 1.0e-8) {
                                    if (self.ntElement[m][k] == elno) {
                                        l = solution.variable.dofs * varContainers->Perm[indexes[j]] + k;
                                        if (self.ntZeroingDone[m][k] == false) {
                                            matContainers->RHS[l] = 0.0;
                                            [self zeroTheNumberOfRows:l inSolutionMatrix:solution];
                                            if (offDiaginal == NO) [self setMatrixElementForSolution:solution atIndex:l andIndex:l value:1.0];
                                            self.ntZeroingDone[m][k] = true;
                                        }
                                        if (offDiaginal == NO) matContainers->RHS[l] = matContainers->RHS[l] + rotvec[k] * work.vector[j]/diagScaling[l];
                                    }
                                }
                            }
                        } else {
                            k = offset + solution.variable.dofs * k + dof;
                            if (solution.matrix.format == MATRIX_SBAND) {
                                bandMatrix = [[FEMMatrixBand alloc] init];
                                [bandMatrix sBand_setDirichlet:solution orderedNumber:k value:work.vector[j]];
                            } else if (solution.matrix.format == MATRIX_CRS && solution.matrix.isSymmetric == YES) {
                                crsMatrix = [[FEMMatrixCRS alloc] init];
                                [crsMatrix setSymmetricDirichletInGlobal:solution atIndex:k value:work.vector[j]/diagScaling[k]];
                            } else {
                                [self zeroTheNumberOfRows:k inSolutionMatrix:solution];
                                if (offDiaginal == NO) {
                                    [self setMatrixElementForSolution:solution atIndex:k andIndex:k value:1.0];
                                    matContainers->RHS[k] = work.vector[j]/diagScaling[k];
                                }
                            }
                        }
                    } else {
                        for (l=0; l<min(solution.variable.dofs, workA.m); l++) {
                            k1 = offset + solution.variable.dofs * k + l;
                            if (solution.matrix.format == MATRIX_SBAND) {
                                bandMatrix = [[FEMMatrixBand alloc] init];
                                [bandMatrix sBand_setDirichlet:solution orderedNumber:k1 value:workA.tensor[l][0][j]];
                            } else if (solution.matrix.format == MATRIX_CRS && solution.matrix.isSymmetric == YES) {
                                if (offDiaginal == NO) {
                                    crsMatrix = [[FEMMatrixCRS alloc] init];
                                    [crsMatrix setSymmetricDirichletInGlobal:solution atIndex:k1 value:workA.tensor[l][0][j]/diagScaling[k1]];
                                }
                            } else {
                                [self zeroTheNumberOfRows:k1 inSolutionMatrix:solution];
                                if (offDiaginal == NO) {
                                    [self setMatrixElementForSolution:solution atIndex:k1 andIndex:k1 value:1.0];
                                    matContainers->RHS[k1] = workA.tensor[l][0][j]/diagScaling[k1];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    if (condition.vector != NULL) {
        free_dvector(condition.vector, 0, condition.m-1);
    }
    if (work.vector != NULL) {
        free_dvector(work.vector, 0, work.m-1);
    }
    if (workA.tensor != NULL) {
        free_d3tensor(workA.tensor, 0,  workA.m-1, 0,  workA.n-1, 0, workA.p-1);
    }
    free_dvector(rotvec, 0, 2);
}

-(void)FEMCore_setPointValues:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution numberofNodes:(int)n atIndexes:(int * __nonnull)indexes forValues:(NSArray * __nonnull)values variableName:(NSMutableString * __nonnull)name orderOfDofs:(int)dof activeCondition:(BOOL)conditional conditionName:(NSString * __nonnull)condName permutationOffset:(int)offset diaginalScaling:(double * __nonnull)diagScaling offDiaginal:(BOOL)offDiaginal {
    
    int j, k, l, k1;
    BOOL stat;
    listBuffer condition = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    listBuffer work = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    listBuffer workA = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    FEMMatrixCRS *crsMatrix;
    FEMMatrixBand *bandMatrix;
    matrixArraysContainer *matContainers = NULL;
    variableArraysContainer *varContainers = NULL;
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
    
    FEMListUtilities *listUtil = [FEMListUtilities sharedListUtilities];
    if (dof >= 0) {
        stat = [listUtil listGetReal:model inArray:values forVariable:name numberOfNodes:n indexes:indexes buffer:&work minValue:NULL maxValue:NULL];
    } else {
        stat = [listUtil listGetRealArray:model inArray:values forVariable:name numberOfNodes:n indexes:indexes buffer:&workA];
    }
    
    if (stat == YES) {
        
        if (conditional == YES) {
            stat = [listUtil listGetReal:model inArray:values forVariable:condName numberOfNodes:n indexes:indexes buffer:&condition minValue:NULL maxValue:NULL];
            conditional = (conditional == YES && stat == YES) ? YES : NO;
        }
        
        for (j=0; j<n; j++) {
            
            @autoreleasepool {
                if (conditional == YES && condition.vector[j] < 0.0) continue;
                if (indexes[j] >= varContainers->sizePerm || indexes[j] < 0) {
                    NSLog(@"FEMCore:FEMCore_setPointValues: invalid node number.\n");
                    continue;
                }
                
                k = varContainers->Perm[indexes[j]];
                
                if (k >= 0) {
                    if (dof >= 0) {
                        k = offset + solution.variable.dofs*k + dof;
                        if (solution.matrix.format == MATRIX_SBAND) {
                            bandMatrix = [[FEMMatrixBand alloc] init];
                            [bandMatrix sBand_setDirichlet:solution orderedNumber:k value:work.vector[j]];
                        } else if (solution.matrix.format == MATRIX_CRS && solution.matrix.isSymmetric == YES) {
                            crsMatrix = [[FEMMatrixCRS alloc] init];
                            [crsMatrix setSymmetricDirichletInGlobal:solution atIndex:k value:work.vector[j]/diagScaling[k]];
                        } else {
                            [self zeroTheNumberOfRows:k inSolutionMatrix:solution];
                            if (offDiaginal == NO) {
                                [self setMatrixElementForSolution:solution atIndex:k andIndex:k value:1.0];
                                matContainers->RHS[k] = work.vector[j]/diagScaling[k];
                            }
                        }
                    } else {
                        for (l=0; l<min(solution.variable.dofs, workA.m); l++) {
                            k1 = offset + solution.variable.dofs*k + l;
                            if (solution.matrix.format== MATRIX_SBAND) {
                                bandMatrix = [[FEMMatrixBand alloc] init];
                                [bandMatrix sBand_setDirichlet:solution orderedNumber:k1 value:workA.tensor[l][0][j]];
                            } else if (solution.matrix.format == MATRIX_CRS && solution.matrix.isSymmetric == YES) {
                                if (offDiaginal == NO) {
                                    crsMatrix = [[FEMMatrixCRS alloc] init];
                                    [crsMatrix setSymmetricDirichletInGlobal:solution atIndex:k1 value:workA.tensor[l][0][j]/diagScaling[k1]];
                                }
                            } else {
                                [self zeroTheNumberOfRows:k1 inSolutionMatrix:solution];
                                if (offDiaginal == NO) {
                                    [self setMatrixElementForSolution:solution atIndex:k1 andIndex:k1 value:1.0];
                                    matContainers->RHS[k1] = workA.tensor[l][0][j]/diagScaling[k1];
                                }
                            }                        
                        }
                    }
                }
            }
        }
    }

    
    if (condition.vector != NULL) {
        free_dvector(condition.vector, 0, condition.m-1);
    }
    if (work.vector != NULL) {
        free_dvector(work.vector, 0, work.m-1);
    }
    
    if (workA.tensor != NULL) {
        free_d3tensor(workA.tensor, 0,  workA.m-1, 0,  workA.n-1, 0, workA.p-1);
    }
}

#pragma mark Periodics

-(void)FEMCore_setPeriodicBoundariesPass1Model:(FEMModel * __nonnull)model solution:(FEMSolution * __nonnull)solution name:(NSMutableString * __nonnull)name dof:(int)dof this:(int)this done:(BOOL * __nonnull)done diaginalScaling:(double * __nonnull)diagScaling {
/*****************************************************************************************************************
    A first pass sum together the rows related to the periodic dofs.
 
    FEMModel *model           -> class containing the model
    FEMSolution *solution     -> solution class containing the matrix and variable
    NSMutableString *name     -> name of the dof to be set
    int dof                   -> the order number of the dof
 
    The permutation (node reordoring info) is contained in the class solution and has been generated at the 
    beginning of the simulation by the bandwidth optimization
 
*****************************************************************************************************************/
    
    int i, j, k, l, m, n, ii, p, q, nn;
    double scale;
    BOOL stat;
    FEMVariable *var;
    FEMListUtilities *listUtil;
    FEMUtilities *util;
    FEMMatrixCRS *crsMatrix;
    FEMBoundaryCondition *boundaryConditionAtId;
    matrixArraysContainer *matContainers = NULL, *projectorContainers = NULL, *aContainers = NULL;
    variableArraysContainer *varContainers = NULL, *containers = NULL;
    
    listUtil = [FEMListUtilities sharedListUtilities];
    util = [[FEMUtilities alloc] init];
    crsMatrix = [[FEMMatrixCRS alloc] init];
    
    NSString *str1 = [@"periodic bc " stringByAppendingString:name];
    NSString *str2 = [@"anti periodic bc " stringByAppendingString:name];
    NSString *str3 = [@"periodic bc scale " stringByAppendingString:name];

    boundaryConditionAtId = (model.boundaryConditions)[this];
    
    if ([listUtil listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:str1 info:&stat] == YES) {
        scale = -1.0;
    } else if ([listUtil listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:str2 info:&stat] == YES) {
        scale = 1.0;
    } else {
        scale = [listUtil listGetConstReal:model inArray:boundaryConditionAtId.valuesList forVariable:str3 info:&stat minValue:NULL maxValue:NULL];
        if (stat == NO) return;
    }
    
    if (boundaryConditionAtId.pMatrix == NULL) return;
    
    matContainers = solution.matrix.getContainers;
    projectorContainers = boundaryConditionAtId.pMatrix.getContainers;
    varContainers = solution.variable.getContainers;
    
    // For explicit conditions just create the dependency almost like a normal Dirichlet BC, 
    // for implicit one (otherwise) do the assembly of the projector
    if ([listUtil listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:@"periodic bc explicit" info:&stat] == YES) {
        
        var = [util getVariableFrom:model.variables model:model name:name onlySearch:NULL maskName:NULL info:&stat];
        containers = var.getContainers;
        
        for (i=0; i<boundaryConditionAtId.pMatrix.numberOfRows; i++) {
            ii = projectorContainers->InvPerm[i];
            k = varContainers->Perm[ii];
            if (done[ii] == NO && k >= 0) {
                k = solution.variable.dofs * k + dof;
                [self zeroTheNumberOfRows:k inSolutionMatrix:solution];
                [self addToMatrixElementForSolution:solution atIndex:k andIndex:k value:1.0];
                matContainers->RHS[k] = 0.0;
                
                for (l=projectorContainers->Rows[i]; l<=projectorContainers->Rows[i+1]-1; l++) {
                    if (projectorContainers->Cols[l] < 0) continue;
                    m = varContainers->Perm[projectorContainers->Cols[l]];
                    if (m >= 0) {
                        m = solution.variable.dofs * m + dof;
                        matContainers->RHS[k] = matContainers->RHS[k] - scale * projectorContainers->Values[l] * containers->Values[m]/diagScaling[k];
                    }
                }
            }
        }
        containers = NULL;
        
    } else if ([listUtil listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:@"periodic bc use lagrange coefficient" info:&stat] == YES) {
        
        FEMMatrix *b = nil;
        FEMMatrix *a = solution.matrix.constraint;
        while (a != NULL) {
            b = boundaryConditionAtId.pMatrix.child;
            while (b != NULL) {
                if (a == b) break;
                b = b.child;
            }
            if (a == b) break;
            a = a.constraint;
        }
        
        if (a != b) {
            FEMMatrix *a = [[FEMMatrix alloc] init];
            a.constraint = solution.matrix.constraint;
            solution.matrix.constraint = a;
            
            a.child = boundaryConditionAtId.pMatrix.child;
            boundaryConditionAtId.pMatrix.child = a;
            
            n = solution.variable.dofs*boundaryConditionAtId.pMatrix.numberOfRows;
            a.numberOfRows = n;
            
            aContainers = a.getContainers;
            
            aContainers->RHS = doublevec(0, n-1);
            aContainers->Rows = intvec(0, (n+1)-1);
            aContainers->Cols = intvec(0, (projectorContainers->sizeCols*pow(solution.variable.dofs, 2.0)+n)-1);
            aContainers->Values = doublevec(0, (projectorContainers->sizeValues*pow(solution.variable.dofs, 2.0)+n)-1);
            
            aContainers->Rows[0] = 0;
            for (i=0; i<boundaryConditionAtId.pMatrix.numberOfRows; i++) {
                n = projectorContainers->Rows[i+1]-projectorContainers->Rows[i];
                for (j=0; j<solution.variable.dofs; j++) {
                    k = solution.variable.dofs*i+j;
                    aContainers->Rows[k+1] = aContainers->Rows[k]+solution.variable.dofs*n;
                }
            }
            
            for (i=0; i<boundaryConditionAtId.pMatrix.numberOfRows; i++) {
                k = varContainers->Perm[projectorContainers->InvPerm[i]];
                for (p=0; p<solution.variable.dofs; p++) {
                    n = aContainers->Rows[solution.variable.dofs*i+p];
                    aContainers->Cols[n] = solution.variable.dofs*k+p;
                    for (j=projectorContainers->Rows[i]; j<=projectorContainers->Rows[i+1]-1; j++) {
                        m = varContainers->Perm[projectorContainers->Cols[j]];
                        if (m >= 0) {
                            for (q=0; q<solution.variable.dofs; q++) {
                                n++;
                                aContainers->Cols[n] = solution.variable.dofs*m+q;
                            }
                        }
                    }
                    
                }
            }
            [crsMatrix sortMatrix:a alsoValues:NULL];
            
            memset( aContainers->RHS, 0.0, n*sizeof(double) );
            memset( aContainers->Values, 0.0, (projectorContainers->sizeValues*pow(solution.variable.dofs, 2.0)+n)*sizeof(double) );
        }
        
        for (i=0; i<boundaryConditionAtId.pMatrix.numberOfRows; i++) {
            k = varContainers->Perm[projectorContainers->InvPerm[i]];
            if (k < 0) continue;
            
            [util zeroTheNumberOfRows:solution.variable.dofs*i+dof inMatrix:a];
            [util setMatrixElement:a atIndex:solution.variable.dofs*i+dof andIndex:solution.variable.dofs*k+dof value:scale];
            
            for (j=projectorContainers->Rows[i]; j<=projectorContainers->Rows[i+1]-1; j++) {
                m = varContainers->Perm[projectorContainers->Cols[j]];
                if (m >= 0) [util setMatrixElement:a atIndex:solution.variable.dofs*i+dof andIndex:solution.variable.dofs*m+dof value:projectorContainers->Values[j]];
            }
        }
        
    } else {
        for (i=0; i<boundaryConditionAtId.pMatrix.numberOfRows; i++) {
            ii = projectorContainers->InvPerm[i];
            k = varContainers->Perm[ii];
            if (done[ii] == NO && k >= 0) {
                k = solution.variable.dofs*k+dof;
                for (l=projectorContainers->Rows[i]; l<=projectorContainers->Rows[i+1]-1; l++) {
                    if (projectorContainers->Cols[l] < 0 || projectorContainers->Values[l] == 0.0) continue;
                    
                    m = varContainers->Perm[projectorContainers->Cols[l]];
                    if (m >= 0 ) {
                        m = solution.variable.dofs*m+dof;
                        for (nn=matContainers->Rows[k]; nn<=matContainers->Rows[k+1]-1; nn++) {
                            [self addToMatrixElementForSolution:solution atIndex:m andIndex:matContainers->Cols[nn] value:projectorContainers->Values[l]*matContainers->Values[nn]/diagScaling[k] * diagScaling[m]];
                            if (matContainers->BulkValues != NULL) {
                                [self addToMatrixElementForSolution:solution atIndex:m andIndex:matContainers->Cols[nn] value:projectorContainers->Values[l]*matContainers->BulkValues[nn]];
                            }
                        }
                        matContainers->RHS[m] = matContainers->RHS[m] + projectorContainers->Values[l]*matContainers->RHS[k] * diagScaling[m]/diagScaling[k];
                        if (matContainers->BulkRHS != NULL) {
                            matContainers->BulkRHS[m] = matContainers->BulkRHS[m] + projectorContainers->Values[l] * matContainers->BulkRHS[k];
                        }
                    }
                }
            }
            done[ii] = YES;
        }
    }
}

-(void)FEMCore_setPeriodicBoundariesPass2Model:(FEMModel * __nonnull)model solution:(FEMSolution * __nonnull)solution name:(NSMutableString * __nonnull)name dof:(int)dof this:(int)this done:(BOOL * __nonnull)done diaginalScaling:(double * __nonnull)diagScaling {
/*****************************************************************************************************************
    At second pass add the [...1.. -1 ...] row structure that results from the equality of the pediodic dofs.
 
    FEMModel *model           -> class containing the model
    FEMSolution *solution     -> solution class containing the matrix and variable
    NSMutableString *name     -> name of the dof to be set
    int dof                   -> the order number of the dof
 
    The permutation (node reordoring info) is contained in the class solution and has been generated at the 
    beginning of the simulation by the bandwidth optimization
 
*****************************************************************************************************************/
    
    int i, k, l, m, ii;
    double scale, valueOffset;
    BOOL stat;
    FEMListUtilities *listUtil;
    FEMBoundaryCondition *boundaryConditionAtId;
    matrixArraysContainer *matContainers = NULL, *projectorContainers = NULL;
    variableArraysContainer *varContainers = NULL;
    listUtil = [FEMListUtilities sharedListUtilities];

    NSString *str1 = [@"periodic bc " stringByAppendingString:name];
    NSString *str2 = [@"anti periodic bc " stringByAppendingString:name];
    NSString *str3 = [@"periodic bc scale " stringByAppendingString:name];
    NSString *str4 = [@"periodic bc offset " stringByAppendingString:name];
    
    boundaryConditionAtId = (model.boundaryConditions)[this];
    
    if ([listUtil listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:str1 info:&stat] == YES) {
        scale = -1.0;
    } else if ([listUtil listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:str2 info:&stat] == YES) {
        scale = 1.0;
    } else {
        scale = [listUtil listGetConstReal:model inArray:boundaryConditionAtId.valuesList forVariable:str3 info:&stat minValue:NULL maxValue:NULL];
        if (stat == NO) return;
    }
    
    valueOffset = [listUtil listGetConstReal:model inArray:boundaryConditionAtId.valuesList forVariable:str4 info:&stat minValue:NULL maxValue:NULL];
    
    if (boundaryConditionAtId.pMatrix == NULL) return;
    
    matContainers = solution.matrix.getContainers;
    projectorContainers = boundaryConditionAtId.pMatrix.getContainers;
    varContainers = solution.variable.getContainers;
    
    str1 = [NSMutableString stringWithString:@"periodic bc explicit"];
    str2 = [NSMutableString stringWithString:@"periodic bc use lagrange coefficient"];
    if ([listUtil listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:str1 info:&stat] == YES) return;
    if ([listUtil listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:str2 info:&stat] == YES) return;
    
    // Do the assembly of the projector
    for (i=0; i<boundaryConditionAtId.pMatrix.numberOfRows; i++) {
        ii = projectorContainers->InvPerm[i];
        k = varContainers->Perm[ii];
        if (done[ii] == NO && k >= 0) {
            k = solution.variable.dofs*k + dof;
            [self zeroTheNumberOfRows:k inSolutionMatrix:solution];
            for (l=projectorContainers->Rows[i]; l<=projectorContainers->Rows[i+1]-1; l++) {
                if (projectorContainers->Cols[l] < 0) continue;
                m = varContainers->Perm[projectorContainers->Cols[l]];
                if (m >= 0) {
                    m = solution.variable.dofs*m + dof;
                    [self addToMatrixElementForSolution:solution atIndex:k andIndex:m value:projectorContainers->Values[l] * diagScaling[m]/diagScaling[k]];
                }
            }
            matContainers->RHS[k] = -valueOffset / diagScaling[k];
            [self addToMatrixElementForSolution:solution atIndex:k andIndex:k value:scale];
        }
        done[ii] = YES;
    }
}


#pragma mark Limiter

/************************************************************************************************
    Determine soft limiters set. This is called after solution and can therefore be active
    only on the second non-linear iteation round.
************************************************************************************************/
-(void)FEMCore_determineSoftLimiterInSolution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model {
    
    int i, j, n, t, ind, dofs, dof, bf, upper, removed, added, elemFirst, elemLast, totSize, count, comp;
    double limitSign, eps;
    BOOL found, onlySearch, newLimiter, anyDofActive, anyLimitBC, anyLimitBF, set;
    bool *limitActive, *limitDone = NULL;
    listBuffer elemLimit = { NULL, NULL, NULL, NULL, 0, 0, 0};
    NSString *limitName;
    NSMutableString *name;
    Element_t *elements;
    FEMVariable *loadVar;
    FEMBodyForce *bodyForceAtID;
    FEMUtilities *utilities;
    FEMListUtilities *listUtilities;
    variableArraysContainer *varContainers = NULL, *loadVarContainers = NULL;
    
    utilities = [[FEMUtilities alloc] init];
    listUtilities = [FEMListUtilities sharedListUtilities];
    
    name = [NSMutableString stringWithString:solution.variable.name];
    [name appendString:@" loads"];
    onlySearch = YES;
    loadVar = [utilities getVariableFrom:model.variables model:model name:name onlySearch:&onlySearch maskName:NULL info:&found];
    if (loadVar == nil) {
        NSLog(@"FEMCore:FEMCore_determineSoftLimiterInSolution: no loads associated wit variable %@.\n", [solution.variable canonicalizeName]);
        return;
    }
    
    varContainers = solution.variable.getContainers;
    loadVarContainers = loadVar.getContainers;
    
    dofs = solution.variable.dofs;
    totSize = varContainers->sizeValues;
    newLimiter = NO;
    
    // Loop though upper and lower limits
    for (upper=0; upper<=1; upper++) {
        removed = 0;
        added = 0;
        
        // If the limiter already exists then check the corresponding load to
        // determine whether the node needs to be released from the set.
        limitActive = NULL;
        
        // Check that active set vectors for limiters exist, otherwise allocate
        if (upper == 0) {
            if (varContainers->lowerLimitActive == NULL) {
                varContainers->lowerLimitActive = boolvec(0, totSize-1);
                varContainers->sizeLowerLimitActive = totSize;
                memset( varContainers->lowerLimitActive, false, totSize*sizeof(bool) );
            }
            limitActive = varContainers->lowerLimitActive;
        } else {
            if (varContainers->upperLimitActive == NULL) {
                varContainers->upperLimitActive = boolvec(0, totSize-1);
                varContainers->sizeUpperLimitActive = totSize;
                memset( varContainers->upperLimitActive, false, totSize*sizeof(bool) );
            }
            limitActive = varContainers->upperLimitActive;
        }
        
        if ((solution.solutionInfo)[@"limiter load tolerance"] != nil) {
            eps = [(solution.solutionInfo)[@"limiter load tolerance"] doubleValue];
        } else eps = DBL_EPSILON;
        
        if (limitDone == NULL) {
            n = model.maxElementNodes;
            limitDone = boolvec(0, totSize-1);
            memset( limitDone, false, totSize*sizeof(bool) );
        }
        
        // These are the default sign that come from standard formulation
        // of Laplace equation
        if (upper == 0) {
            limitSign = -1.0;
        } else {
            limitSign = 1.0;
        }
        
        // The user may want to toggle the sign for other kinds of equations
        if ((solution.solutionInfo)[@"limiter load sign negative"] != nil) {
            if ([(solution.solutionInfo)[@"limiter load sign negative"] boolValue] == YES) limitSign = -1.0 * limitSign;
        }
        
        // Go through the active set and free nodes with wrong sign
        for (i=0; i<totSize; i++) {
            if (limitActive[i] == true) {
                if (limitSign * loadVarContainers->Values[i] > limitSign * eps) {
                    removed++;
                    limitActive[i] = false;
                    limitDone[i] = true;
                }
            }
        }
        
        // Go through the field variables one dof at a time since typically
        // for vector valued field, the limiter is given component wise
        anyDofActive = NO;
        for (dof=0; dof<dofs; dof++) {
            name = [NSMutableString stringWithString:solution.variable.name];
            if (solution.variable.dofs > 1) {
                comp = dof+1;
                name = (NSMutableString *)[utilities appendNameFromString:name component:&comp];
            }
            
            if (upper == 0) {
                limitName = [name stringByAppendingString:@" lower limit"];
            } else {
                limitName = [name stringByAppendingString:@" upper limit"];
            }
            
            // Check whether limiters exist at all
            anyLimitBC = [listUtilities listCheckPresentAnyBoundaryCondition:model name:limitName];
            anyLimitBF = [listUtilities listCheckPresentAnyBodyForce:model name:limitName];
            
            if (!(anyLimitBC || anyLimitBF)) continue;
            anyDofActive = YES;
            
            if (limitDone == NULL) {
                n = model.maxElementNodes;
                limitDone = boolvec(0, totSize-1);
                memset( limitDone, false, totSize*sizeof(bool) );
            }
            
            if ((solution.solutionInfo)[@"limiter value tolerance"] != nil) {
                eps = [(solution.solutionInfo)[@"limiter value tolerance"] doubleValue];
            } else eps = DBL_EPSILON;
            
            // Define the range of elements for which the limiters are active
            elemFirst = model.numberOfBulkElements;
            elemLast = model.numberOfBulkElements;
            
            if (anyLimitBF == YES) elemFirst = 0;
            if (anyLimitBC == YES) elemLast = model.numberOfBulkElements + model.numberOfBoundaryElements;
            
            // Go through the elements
            elements = model.getElements;
            for (t=elemFirst; t<elemLast; t++) {
                n = elements[t].Type.NumberOfNodes;
                
                if (t >= model.numberOfBulkElements) {
                    found = NO;
                    for (FEMBoundaryCondition *boundaryCondition in model.boundaryConditions) {
                        if (elements[t].BoundaryInfo->Constraint == boundaryCondition.tag) {
                            found = [listUtilities listGetReal:model inArray:boundaryCondition.valuesList forVariable:limitName numberOfNodes:n indexes:elements[t].NodeIndexes buffer:&elemLimit minValue:NULL maxValue:NULL];
                            break;
                        }
                    }
                    if (found == NO) continue;
                } else {
                    if ((model.bodies)[elements[t].BodyID-1][@"body force"] != nil) {
                        bf = [(model.bodies)[elements[t].BodyID-1][@"body force"] intValue];
                    } else {
                        continue;
                    }
                    bodyForceAtID = (model.bodyForces)[bf-1];
                    found = [listUtilities listGetReal:model inArray:bodyForceAtID.valuesList forVariable:limitName numberOfNodes:n indexes:elements[t].NodeIndexes buffer:&elemLimit minValue:NULL maxValue:NULL];
                    if (found == NO) continue;
                    
                }
                
                for (i=0; i<n; i++) {
                    j = varContainers->Perm[elements[t].NodeIndexes[i]];
                    if (j < 0) continue;
                    ind = dofs * j + dof;
                    
                    if (limitDone[ind] == YES) continue;
                    
                    if (upper == 0) {
                        set = (varContainers->Values[ind] < elemLimit.vector[i] - eps) ? YES : NO;
                    } else {
                        set = (varContainers->Values[ind] > elemLimit.vector[i] + eps) ? YES : NO;
                    }
                    
                    if (set == YES) {
                        if (limitActive[ind] == NO) added++;
                        limitActive[ind] = YES;
                        limitDone[ind] = YES;
                    }
                    
                    // Enforce the value to limits because nonlinear material models
                    // may oltherwise lead to divergence of the iteration
                    if (upper == 0) {
                        varContainers->Values[ind] = max(varContainers->Values[ind], elemLimit.vector[i]);
                    } else {
                        varContainers->Values[ind] = min(varContainers->Values[ind], elemLimit.vector[i]);
                    }
                }
            }
        }
        
        if (anyDofActive == NO) continue;
        
        // Output some information before exiting
        if (upper == 0) {
            NSLog(@"FEMCore:FEMCore_determineSoftLimiterInSolution: determined lower soft limit set.");
        } else {
            NSLog(@"FEMCore:FEMCore_determineSoftLimiterInSolution: determined upper soft limit set.");
        }
        count = 0;
        for (i=0; i<totSize; i++) {
            if (limitActive[i] == YES) count++;
        }
        NSLog(@"FEMCore:FEMCore_determineSoftLimiterInSolution: number of dofs in set is %d.\n", count);
        
        if (added > 0) {
            NSLog(@"FEMCore:FEMCore_determineSoftLimiterInSolution: added %d dofs to the set.\n", added);
        }
        if (removed > 0) {
            NSLog(@"FEMCore:FEMCore_determineSoftLimiterInSolution: removed %d dofs to the set.\n", removed);
        }
    }
    
    if (limitDone != NULL) {
        free_bvector(limitDone, 0, totSize-1);
    }
    
    if (elemLimit.vector != NULL) {
        free_dvector(elemLimit.vector, 0, elemLimit.m-1);
    }
}

-(void)FEMCore_setLimiterValues:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution numberofNodes:(int)n atIndexes:(int * __nonnull)indexes forValues:(NSArray * __nonnull)values  orderOfDofs:(int)dof limitActive:(bool * __nonnull)limitActive conditionName:(NSString * __nonnull)condName permutationOffset:(int)offset diaginalScaling:(double * __nonnull)diagScaling offDiaginal:(BOOL)offDiaginal {
    
    int j, k;
    BOOL stat;
    listBuffer work = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    FEMMatrixCRS *crsMatrix;
    FEMMatrixBand *bandMatrix;
    matrixArraysContainer *matContainers = NULL;
    variableArraysContainer *varContainers = NULL;
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
    
    FEMListUtilities *listUtil = [FEMListUtilities sharedListUtilities];
    stat = [listUtil listGetReal:model inArray:values forVariable:condName numberOfNodes:n indexes:indexes buffer:&work minValue:NULL maxValue:NULL];
    
    if (stat == YES) {
        
        for (j=0; j<n; j++) {
            
            @autoreleasepool {
                k = varContainers->Perm[indexes[j]];
                if (k < 0) continue;
                if (!limitActive[solution.variable.dofs*k+dof]) continue;
                
                k = offset + solution.variable.dofs*k + dof;
                if (solution.matrix.format == MATRIX_SBAND) {
                    bandMatrix = [[FEMMatrixBand alloc] init];
                    [bandMatrix sBand_setDirichlet:solution orderedNumber:k value:work.vector[j]];
                } else if (solution.matrix.format == MATRIX_CRS && solution.matrix.isSymmetric == YES) {
                    crsMatrix = [[FEMMatrixCRS alloc] init];
                    [crsMatrix setSymmetricDirichletInGlobal:solution atIndex:k value:work.vector[j]/diagScaling[k]];
                } else {
                    [self zeroTheNumberOfRows:k inSolutionMatrix:solution];
                    if (offDiaginal == NO) {
                        [self setMatrixElementForSolution:solution atIndex:k andIndex:k value:1.0];
                        matContainers->RHS[k] = work.vector[j]/diagScaling[k];
                    }
                }
            }
        }
    }
    
    if (work.vector != NULL) {
        free_dvector(work.vector, 0, work.m-1);
    }
}

#pragma mark Rows equilibration
/****************************************************************************************
    Equilibrate the rows of the coefficient matrix in solution to minimize the condition 
    number. The associated rhs vector is also scaled
*****************************************************************************************/
-(void)FEMCore_rowEquilibrationMatrix:(FEMMatrix * __nonnull)matrix vector:(double * __nonnull)f parallel:(BOOL)parallel {
    
    int i, j, n;
    double norm, tmp;
    double complex cmpx;
    BOOL complexMatrix;
    matrixArraysContainer *matContainers = NULL;
    
    n = matrix.numberOfRows;
    complexMatrix = matrix.complexMatrix;
    
    matContainers = matrix.getContainers;
    
    if (matContainers->DiagScaling == NULL) {
        matContainers->DiagScaling = doublevec(0, n-1);
        matContainers->sizeDiagScaling = n;
    }
    
    memset( matContainers->DiagScaling, 0.0, n*sizeof(double) );
    norm = 0.0;

    // Compute 1-norm of each row
    if (complexMatrix == YES) {
        for (i=0; i<n; i+=2) {
            tmp = 0.0;
            for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j+=2) {
                cmpx = matContainers->Values[j] + (-matContainers->Values[j+1] * I);
                tmp = tmp + cabs(cmpx);
            }
            
            if (parallel == NO) {
                if (tmp > norm) norm = tmp;
            }
            
            if (tmp > 0.0) {
                matContainers->DiagScaling[i] = tmp;
                matContainers->DiagScaling[i+1] = tmp;
            }
        }
    } else {
        for (i=0; i<n; i++) {
            tmp = 0.0;
            for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
                tmp = tmp + fabs(matContainers->Values[j]);
            }
            
            if (parallel == NO) {
                if (tmp > norm) norm = tmp;
            }
            
            if (tmp > 0.0) matContainers->DiagScaling[i] = tmp;
        }
    }
    
    if (parallel == YES) {
        //TODO: Add parallel run support
    }
    
    // Now define the scalig matrix by inversion and perform
    // the actual scaling of the linear system
    if (complexMatrix == YES) {
        for (i=0; i<n; i+=2) {
            if (matContainers->DiagScaling[i] > 0.0) {
                matContainers->DiagScaling[i] = 1.0 / matContainers->DiagScaling[i];
                matContainers->DiagScaling[i+1] = 1.0 / matContainers->DiagScaling[i+1];
            } else {
                matContainers->DiagScaling[i] = 1.0;
                matContainers->DiagScaling[i+1] = 1.0;
            }
        }
    } else {
        for (i=0; i<n; i++) {
            if (matContainers->DiagScaling[i] > 0.0) {
                matContainers->DiagScaling[i] = 1.0 / matContainers->DiagScaling[i];
            } else {
                matContainers->DiagScaling[i] = 1.0;
            }
        }
    }
    
    for (i=0; i<n; i++) {
        for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
            matContainers->Values[j] = matContainers->Values[j] * matContainers->DiagScaling[i];
        }
        f[i] = matContainers->DiagScaling[i] * f[i];
    }
    
    NSLog(@"FEMCore:FEMCore_rowEquilibrationMatrix: unscaled matrix norm: %f.\n", norm);
}

/****************************************************************************************
    Scale the linear system back to original when the linear system scaling has been
    done by row equilibration.
*****************************************************************************************/
-(void)FEMCore_reverseRowEquilibrationMatrix:(FEMMatrix * __nonnull)matrix vector:(double * __nonnull)f {
    
    int i, j, n;
    matrixArraysContainer *matContainers = NULL;
    
    n = matrix.numberOfRows;
    matContainers = matrix.getContainers;
    
    if (matContainers->DiagScaling == NULL) {
        fatal("FEMCore:FEMCore_reverseRowEquilibrationMatrix", "Diag is a null pointer.");
    }
    if (matContainers->sizeDiagScaling != n) {
        fatal("FEMCore:FEMCore_reverseRowEquilibrationMatrix", "Diag of wrong size.");
    }
    
    for (i=0; i<n; i++) {
        f[i] = f[i] / matContainers->DiagScaling[i];
    }
    for (i=0; i<n; i++) {
        for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
            matContainers->Values[j] = matContainers->Values[j] / matContainers->DiagScaling[i];
        }
    }
    
    free_dvector(matContainers->DiagScaling, 0, matContainers->sizeDiagScaling-1);
    matContainers->sizeDiagScaling = 0;
    matContainers->DiagScaling = NULL;
}

-(void)FEMCore_solveCoupledModel:(FEMModel * __nonnull)model timeStep:(double)dt coupledMinIteration:(int)coupledMinIter coupleMaxIteration:(int)coupleMaxIter transientSimulation:(BOOL)transient scanning:(BOOL)scanning doneThis:(BOOL * __nonnull)doneThis afterConverged:(BOOL * __nonnull)afterConverged steadyIt:(double * __nonnull)steadyIt {
    
    int  i, j, k, n;
    double relaxation;
    BOOL found, allAfterConverged, allDoneThis, calculateDerivative = NO, needSol, testConvergence;
    NSString *when;
    FEMMesh *mesh;
    FEMListUtilities *listUtilities;
    variableArraysContainer *varContainers = NULL;
    
    for (i=1; i<=coupleMaxIter; i++) {
        if (transient == YES || scanning == YES) {
            if (coupleMaxIter > 1) {
                NSLog(@"FEMCore:FEMCore_solveCoupledModel: --------------------------------------------------------------------\n");
                NSLog(@"FEMCore:FEMCore_solveCoupledModel: coupled system iteration: %d.\n", i);
                NSLog(@"FEMCore:FEMCore_solveCoupledModel: --------------------------------------------------------------------\n");
            }
            *steadyIt = (double)i;
        }
        
        memset( doneThis, 0, model.numberOfSolutions*sizeof(BOOL) );
        
        // Initialize the mesh output flag to NO here, reactivated later
        // for meshes connected to active solutions.
        for (FEMMesh *mesh in model.meshes) {
            mesh.outputActive = NO;
        }
        
        // Go through number of solutions (heat, laminar or turbulent flow, etc...)
        k = 0;
        for (FEMSolution *solution in model.solutions) {
            if (solution.hasBuiltInSolution == NO && solution.plugInPrincipalClassInstance == nil) {
                if (!(solution.solutionMode == SOLUTION_MODE_COUPLED || solution.solutionMode == SOLUTION_MODE_ASSEMBLY
                    || solution.solutionMode == SOLUTION_MODE_BLOCK)) {
                    NSLog(@"FEMCore:SolveEquations: no routine related to solution.\n");
                    doneThis[k] = YES;
                    k++;
                    continue;
                }
            }
            
            if ((solution.solutionInfo)[@"invoke solution computer"] != nil) {
                when = (solution.solutionInfo[@"invoke solution computer"]);
                if ([when isEqualToString:@"always"] == NO) {
                    doneThis[k] = YES;
                    k++;
                    continue;
                }
            } else {
                if (solution.solutionSolveWhen != SOLUTION_SOLVE_ALWAYS) {
                    doneThis[k] = YES;
                    k++;
                    continue;
                }
            }
            
            allAfterConverged = YES;
            for (j=0; j<model.numberOfSolutions; j++) {
                if (afterConverged[j] == NO) {
                    allAfterConverged = NO;
                    break;
                }
            }
            allDoneThis = YES;
            for (j=0; j<model.numberOfSolutions; j++) {
                if (doneThis[j] == NO) {
                    allDoneThis = NO;
                    break;
                }
            }
            if (afterConverged[k] == YES && !(allAfterConverged == YES || allDoneThis == YES)) {
                k++;
                continue;
            }
            
            n = 0;
            if (solution.variable != nil) {
                varContainers = solution.variable.getContainers;
                if (varContainers->Values != NULL) n = varContainers->sizeValues;
                solution.variable.prevNorm = solution.variable.norm;
            }
            
            // There are some operations that require that the previous steady state values
            // are present. Check for these operations.
            if ((solution.solutionInfo)[@"calculate derivative"] != nil) {
                calculateDerivative = [(solution.solutionInfo)[@"calculate derivative"] boolValue];
            }
            needSol = calculateDerivative;
            if (needSol == NO) {
                if ((solution.solutionInfo)[@"steady state convergence measure"] != nil) {
                    needSol = ([(solution.solutionInfo)[@"steady state convergence measure"] isEqualToString:@"norm"] == NO) ? YES : NO;
                    needSol = (needSol == YES) ? YES : NO;
                }
            }
            if (needSol == NO) {
                if ((solution.solutionInfo)[@"steady state relaxation factor"] != nil) {
                    relaxation = [(solution.solutionInfo)[@"steady state relaxation factor"] doubleValue];
                    needSol = (relaxation != 1.0) ? YES : NO;
                }
            }
            
            if (needSol == YES && n > 0) {
                found = (varContainers->SteadyValues != NULL) ? YES : NO;
                if (found == YES && varContainers->sizeSteadyValues != n) {
                    free_dvector(varContainers->SteadyValues, 0, varContainers->sizeSteadyValues-1);
                    varContainers->SteadyValues = NULL;
                    found = NO;
                }
                if (found == NO) {
                    varContainers->SteadyValues = doublevec(0, n-1);
                    varContainers->sizeSteadyValues = n;
                }
                for (j=0; j<n; j++) {
                    varContainers->SteadyValues[j] = varContainers->Values[j];
                }
            }
            
            [self activateSolution:solution model:model timeStep:dt transientSimulation:transient];
            
            if (solution.variable != nil) {
                if (varContainers->Values != NULL) n = varContainers->sizeValues;
            }
            
            // Check for coupled system convergence
            if (scanning == YES) {
                testConvergence = (i >= coupledMinIter && i != coupleMaxIter) ? YES : NO;
                if (testConvergence == YES || calculateDerivative == YES) {
                    // TODO: add support for parallel run
                    [self computeChange:solution model:model isSteadyState:YES nsize:NULL values:NULL values0:NULL sizeValues0:NULL];
                    
                    // The ComputeChange method sets a flag to zero if not yet converged (otherwise -1/1)
                    if (testConvergence == YES) {
                        if (solution.variable.steadyConverged == 0) {
                            doneThis[k] = NO;
                        } else {
                            doneThis[k] = YES;
                        }
                    }
                }
            } else if (transient == YES) {
                testConvergence = (i >= coupledMinIter && i != coupleMaxIter) ? YES : NO;
                if (testConvergence == YES) {
                    // TODO: add support for parallel run
                    [self computeChange:solution model:model isSteadyState:YES nsize:NULL values:NULL values0:NULL sizeValues0:NULL];
                    if (solution.variable.steadyConverged == 0) {
                        doneThis[k] = NO;
                    } else {
                        doneThis[k] = YES;
                    }
                }
            } else { // Steady state
                // TODO: add support for parallel run
                [self computeChange:solution model:model isSteadyState:YES nsize:NULL values:NULL values0:NULL sizeValues0:NULL];
                if (solution.variable.steadyConverged == 0) {
                    doneThis[k] = NO;
                } else {
                    doneThis[k] = YES;
                }
            }
            
            // TODO: add support for parallel run
            
            allDoneThis = YES;
            for (j=0; j<model.numberOfSolutions; j++) {
                if (doneThis[j] == NO) {
                    allDoneThis = NO;
                    break;
                }
            }
            if (allDoneThis == YES) break;
            k++;
        }
        mesh = (FEMMesh *)model.mesh;
        mesh.changed = NO;
        
        allDoneThis = YES;
        for (j=0; j<model.numberOfSolutions; j++) {
            if (doneThis[j] == NO) {
                allDoneThis = NO;
                break;
            }
        }
        if (allDoneThis == YES) break;
    }
    
    allDoneThis = YES;
    for (j=0; j<model.numberOfSolutions; j++) {
        if (doneThis[j] == NO) {
            allDoneThis = NO;
            break;
        }
    }
    if (transient == YES && allDoneThis == NO) {
        listUtilities = [FEMListUtilities sharedListUtilities];
        if ([listUtilities listGetLogical:model inArray:model.simulation.valuesList forVariable:@"coupled system abort not converged" info:&found] == YES) {
            NSLog(@"FEMCore:SolveEquations: \n");
            NSLog(@"FEMCore:SolveEquations: coupled system iteration: %d.\n", i);
            NSLog(@"FEMCore:SolveEquations: \n");
            fatal("FEMCore:SolveEquations", "Saino will abort the simulation now...");
        }
    }
}

/*************************************************************************************************************************
    This executes the original line of solutions (legacy solutions) where each solution includes looping over elements
    and the convergence control. From generality point of view, this misses some opportunities to have control of the
    nonlinear system.
**************************************************************************************************************************/
-(void)FEMCore_singleSolution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model timeStep:(double)dt transientSimulation:(BOOL)transient {
    
    int i, maxDim;
    BOOL meActive;
    NSString *equationName;
    solutionArraysContainer *solContainers = NULL;
    Element_t *elements = NULL;
    
    if (solution.mesh.changed == YES || solution.numberOfActiveElements <= 0) {
        solution.numberOfActiveElements = 0;
        if ((solution.solutionInfo)[@"equation"] != nil) {
            equationName = (solution.solutionInfo)[@"equation"];
            solContainers = solution.getContainers;
            if (solContainers->activeElements != NULL) free_ivector(solContainers->activeElements, 0, solContainers->sizeActiveElements-1);
            solContainers->activeElements = intvec(0, (solution.mesh.numberOfBulkElements+solution.mesh.numberOfBoundaryElements)-1);
            solContainers->sizeActiveElements = solution.mesh.numberOfBulkElements+solution.mesh.numberOfBoundaryElements;
            
            elements = solution.mesh.getElements;
            maxDim = 0;
            FEMListUtilities *listUtilities = [FEMListUtilities sharedListUtilities];
            for (i=0; i<solution.mesh.numberOfBulkElements+solution.mesh.numberOfBoundaryElements; i++) {
                if ([listUtilities checkElementEquation:model forElement:&elements[i] andEquation:equationName] == YES) {
                    solContainers->activeElements[solution.numberOfActiveElements] = i;
                    maxDim = max(elements[i].Type.dimension, maxDim);
                    solution.numberOfActiveElements++;
                }
            }
            [solution.solutionInfo setObject:@(maxDim) forKey:@"active mesh dimension"];
            
            if ([(solution.solutionInfo)[@"calculate weights"] boolValue] == YES) {
                [self computeNodalWeightsInSolution:solution model:model weightBoundary:NO perm:NULL sizePerm:NULL variableName:nil];
            }
            if ([(solution.solutionInfo)[@"calculate boundary weights"] boolValue] == YES) {
                [self computeNodalWeightsInSolution:solution model:model weightBoundary:YES perm:NULL sizePerm:NULL variableName:nil];
            }
        }
    }
    
    meActive = (solution.matrix != nil) ? YES : NO;
    if (meActive == YES) meActive = (meActive == YES && (solution.matrix.numberOfRows > 0)) ? YES : NO;
    // TODO: add support for parallel run
    
    if (solution.matrix != nil) {
        // TODO: add support for parallel run
    }
    
    if (solution.hasBuiltInSolution == YES) { // Call a built-in solution computer here after allocating it only once. It will be deallocated with its asssociated solution
        if ([(solution.solutionInfo)[@"equation"] isEqualToString:@"navier-stokes"] == YES) {
            if (solution.builtInSolution == nil) {
                FEMFlowSolution *flowSolution = [[FEMFlowSolution alloc] init];
                solution.builtInSolution = flowSolution;
            }
            [solution.builtInSolution solutionComputer:solution model:model timeStep:dt transientSimulation:transient];
        } else if ([(solution.solutionInfo)[@"equation"] isEqualToString:@"magnetic induction"] == YES) {
            if (solution.builtInSolution == nil) {
                FEMMagneticInductionSolution *magneticInductionSolution = [[FEMMagneticInductionSolution alloc] init];
                solution.builtInSolution = magneticInductionSolution;
            }
            [solution.builtInSolution solutionComputer:solution model:model timeStep:dt transientSimulation:transient];
        } else if ([(solution.solutionInfo)[@"equation"] isEqualToString:@"stress analysis"] == YES) {
            if (solution.builtInSolution == nil) {
                FEMStressAnalysisSolution *stressAnalysisSolution = [[FEMStressAnalysisSolution alloc] init];
                solution.builtInSolution = stressAnalysisSolution;
            }
            [solution.builtInSolution solutionComputer:solution model:model timeStep:dt transientSimulation:transient];
        } else if ([(solution.solutionInfo)[@"equation"] isEqualToString:@"mesh update"] == YES) {
            if (solution.builtInSolution == nil) {
                FEMMeshUpdateSolution *meshUpdateSolution = [[FEMMeshUpdateSolution alloc] init];
                solution.builtInSolution = meshUpdateSolution;
            }
            [solution.builtInSolution solutionComputer:solution model:model timeStep:dt transientSimulation:transient];
        } else if ([(solution.solutionInfo)[@"equation"] isEqualToString:@"heat equation"] == YES) {
            if (solution.builtInSolution == nil) {
                if ([(solution.solutionInfo)[@"parallel assembly"] boolValue] == YES) {
                     FEMHeatSolution_OpenCL *heatSolution = [[FEMHeatSolution_OpenCL alloc] init];
                    solution.builtInSolution = heatSolution;
                } else {
                    FEMHeatSolution *heatSolution = [[FEMHeatSolution alloc] init];
                    solution.builtInSolution = heatSolution;
                }
            }
            [solution.builtInSolution solutionComputer:solution model:model timeStep:dt transientSimulation:transient];
        } else {
            NSLog(@"FEMCore:FEMCore_singleSolution: can't proceed further because there is no class implementation for %@.\n", (solution.solutionInfo)[@"equation"]);
            fatal("FEMCore:FEMCore_singleSolution", "Saino will abort the simulation now...");
        }

    } else if (solution.plugInPrincipalClassInstance != nil) { // Otherwise call the user-provided plug-in, it will also be deallocated with its asssociated solution
        _instance = solution.plugInPrincipalClassInstance;
        [_instance solutionComputer:solution model:model timeStep:dt transientSimulation:transient];
    } else {
        NSLog(@"FEMCore:FEMCore_singleSolution: can't proceed further because there is built-in computer or plug-in associated with the equation %@.\n", (solution.solutionInfo)[@"equation"]);
        fatal("FEMCore:FEMCore_singleSolution", "Saino will abort the simulation now...");
    }
}

-(void)FEMCore_updateExportedVariablesSolution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model {
    
    int i, j, l, m, dofs;
    double *localSol, *localCond;
    NSString *varName;
    FEMUtilities *utilities;
    FEMVariable *expVariable;
    variableArraysContainer *expVariableContainers = NULL;
    static BOOL allocationdDone = NO;
    BOOL *activePart, *activeCond, found, stateVariable, variableOutput, variableGlobal;
    
    utilities = [[FEMUtilities alloc] init];
    
    l = 0;
    while (1) {
        l++;
        NSString *componentName = [@"exported variable " stringByAppendingString:[@(l) stringValue]];
        if ((solution.solutionInfo)[componentName] != nil) {
            varName = (solution.solutionInfo)[componentName];
            found = YES;
        } else found = NO;
        
        if (found == NO) break;
        
        if ([varName characterAtIndex:0] == '-') {
            if ([[[varName substringFromIndex:0] substringToIndex:10] isEqualToString:@"-nooutput "] == YES) {
                variableOutput = NO;
                varName = [varName substringFromIndex:10];
            }
            if ([[[varName substringFromIndex:0] substringToIndex:8] isEqualToString:@"-global "] == YES) {
                variableGlobal = YES;
                varName = [varName substringFromIndex:8];
            }
            if ([[[varName substringFromIndex:0] substringToIndex:6] isEqualToString:@"-dofs "] == YES) {
                i = 5;
                while ([varName characterAtIndex:i] == ' ') {
                    i++;
                }
                j = i;
                while (1) {
                    if ([varName characterAtIndex:j] == ' ') break;
                    j++;
                }
                dofs = [[[varName substringFromIndex:i] substringToIndex:j-i] intValue];
                // If we get 0 then we are dealing with an invalid dof number
                if (dofs == 0) {
                    NSLog(@"FEMUtilities:addEquationBasicsToSolution: invalid dof in variable definiton.\n");
                    NSLog(@"FEMUtilities:addEquationBasicsToSolution: the incorrect value was: %@.\n", [[varName substringFromIndex:i] substringToIndex:j-i]);
                }
                varName = [varName substringFromIndex:j+1];
            }
        }
        
        expVariable = [utilities getVariableFrom:solution.mesh.variables model:model name:varName onlySearch:NULL maskName:NULL info:&found];
        if (expVariable == nil) continue;
        
        NSLog(@"FEMCore:FEMCore_updateExportedVariablesSolution: trying to set values for variable: %@.\n", varName);
        
        if (allocationdDone == NO) {
            m = model.numberOfBodyForces;
            activePart = (BOOL*)malloc(sizeof(BOOL) * m );
            activeCond = (BOOL*)malloc(sizeof(BOOL) * m );
            
            m = solution.mesh.maxElementDofs;
            localSol = doublevec(0, m-1);
            localCond = doublevec(0, m-1);
            
            allocationdDone = YES;
        }
        
        dofs = expVariable.dofs;
        expVariableContainers = expVariable.getContainers;
        
        stateVariable = (expVariableContainers->sizeValues == dofs);
        if (stateVariable == YES) {
            if (dofs > 1) {
                
            }
        }
    }
    
    //TODO: implementation incomplete, seems that there is a mistake in the Elmer implementation
}

#pragma mark Public methods...

@synthesize normalTangentialNumberOfNodes = _normalTangentialNumberOfNodes;
@synthesize size1NtZeroingDone = _size1NtZeroingDone;
@synthesize size2NtZeroingDone = _size2NtZeroingDone;
@synthesize sizeBoundaryReorder = _sizeBoundaryReorder;
@synthesize size1NtElement = _size1NtElement;
@synthesize size2NtElement = _size2NtElement;
@synthesize size1boundaryNormals = _size1boundaryNormals;
@synthesize size2boundaryNormals = _size2boundaryNormals;
@synthesize size1boundaryTangent1 = _size1boundaryTangent1;
@synthesize size2boundaryTangent1 = _size2boundaryTangent1;
@synthesize size1boundaryTangent2 = _size1boundaryTangent2;
@synthesize size2boundaryTangent2 = _size2boundaryTangent2;
@synthesize sizeIndexStore = _sizeIndexStore;
@synthesize indexStore = _indexStore;
@synthesize ntZeroingDone = _ntZeroingDone;
@synthesize boundaryReorder = _boundaryReorder;
@synthesize ntElement = _ntElement;
@synthesize boundaryNormals = _boundaryNormals;
@synthesize boundaryTangent1 = _boundaryTangent1;
@synthesize boundaryTangent2 = _boundaryTangent2;
@synthesize normalTangentialName = _normalTangentialName;
@synthesize outputLevelMask = _outputLevelMask;
@synthesize outputPrefix = _outputPrefix;
@synthesize outputCaller = _outputCaller;
@synthesize maxOutputLevel = _maxOutputLevel;
@synthesize minOutputLevel = _minOutputLevel;
@synthesize outputPE = _outputPE;

#pragma mark Singleton method

static FEMCore * __nullable sharedCore = nil;
static dispatch_once_t onceToken;

+(id __nonnull)sharedCore {
    
    dispatch_once(&onceToken, ^{
        sharedCore = [[self alloc] init];
    });
    return sharedCore;
}

+(void) selfDestruct {
    sharedCore = nil;
    onceToken = 0;
}

- (id)init
{
    int i;
    
    self = [super init];
    if (self) {
        // Initialization code here.
        _indexStore = intvec(0, 511);
        _sizeIndexStore = 512;
        memset( _indexStore, -1, _sizeIndexStore*sizeof(int) );
        
        _lineEM = intmatrix(0, 0, 0, 1);
        _triangleEM = intmatrix(0, 2, 0, 1);
        _quadEM = intmatrix(0, 3, 0, 1);
        _tetraEM = intmatrix(0, 5, 0, 1);
        _prismEM = intmatrix(0, 7, 0, 1);
        _wedgeEM = intmatrix(0, 8, 0, 1);
        _brickEM = intmatrix(0, 11, 0, 1);
        
        _ntZeroingDone = NULL;
        _boundaryReorder = NULL;
        _ntElement = NULL;
        _boundaryNormals = NULL;
        _boundaryTangent1 = NULL;
        _boundaryTangent2 = NULL;
        
        _sizeBoundaryReorder = 0;
        _size1boundaryNormals = 0;
        _size2boundaryNormals = 0;
        _size1boundaryTangent1 = 0;
        _size2boundaryTangent1 = 0;
        _size1boundaryTangent2 = 0;
        _size2boundaryTangent2 = 0;
        //TODO: Allocate boundaryReorder, boundaryNormals, boundaryTangent1, boundaryTangent2
        
        for (i=0; i<8; i++) {
            _initialized[i] = NO;
        }
        
        _outputLevelMask = [[NSMutableArray alloc] init];
        for (i=0; i<32; i++) {
            _outputLevelMask[i] = @YES;
        }
        _outputPrefix = NO;
        _outputCaller = YES;
        _maxOutputLevel = 32;
        _minOutputLevel = 0;
        _outputPE = 0;
        
        _normalTangentialName = [NSMutableString stringWithString:@""];
        
        _n1 = 0;
        _k1 = 0;
        _kernStiff = NULL;
        _kernWork = NULL;
        _g_Ind = NULL;
        _l_Ind = NULL;
        _saveValues = NULL;
        _damp = NULL;
        _stiff = NULL;
        _mass = NULL;
        _x = NULL;
        _instance = nil;
    }
    
    return self;
}

-(void)deallocation {
    free_ivector(_indexStore, 0, 511);
    free_imatrix(_lineEM, 0, 0, 0, 1);
    free_imatrix(_triangleEM, 0, 2, 0, 1);
    free_imatrix(_quadEM, 0, 3, 0, 1);
    free_imatrix(_tetraEM, 0, 5, 0, 1);
    free_imatrix(_prismEM, 0, 7, 0, 1);
    free_imatrix(_wedgeEM, 0, 8, 0, 1);
    free_imatrix(_brickEM, 0, 11, 0, 1);
    
    if (_g_Ind != NULL) {
        free_ivector(_g_Ind, 0, _size_g_Ind-1);
        _g_Ind = NULL;
    }
    
    if (_l_Ind != NULL) {
        free_ivector(_l_Ind, 0, _size_l_Ind-1);
        _l_Ind = NULL;
    }
    
    if (_kernStiff != NULL) {
        free_dmatrix(_kernStiff, 0, _size1kernStiff-1, 0, _size2kernStiff-1);
        _kernStiff = NULL;
    }
    
    if (_kernWork != NULL) {
        free_dvector(_kernWork, 0, _sizekernWork-1);
        _kernWork = NULL;
    }
    
    if (_stiff != NULL) {
        free_dmatrix(_stiff, 0, 0, 0, _n1-1);
    }
    if (_mass != NULL) {
        free_dmatrix(_mass, 0, 0, 0, _n1-1);
    }
    if (_damp != NULL) {
        free_dmatrix(_damp, 0, 0, 0, _n1-1);
    }
    if (_x != NULL) {
        free_dmatrix(_x, 0, _n1-1, 0, _k1-1);
    }
}


// Initialize solution of next time step
-(void)initializeTimeStepInSolution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model {
    
    int i, j, order;
    BOOL found;
    NSString *method;
    FEMListUtilities *listUtilities;
    variableArraysContainer *varContainers = NULL;
    matrixArraysContainer *matContainers = NULL;
    
    solution.doneTime++;
    
    listUtilities = [FEMListUtilities sharedListUtilities];
    varContainers = solution.variable.getContainers;
    
    if (solution.matrix == nil || varContainers->Values == NULL) return;
    if (solution.timeOrder <= 0) return;
    
    found = NO;
    if ((solution.solutionInfo)[@"time stepping method"] != nil) {
        method = (solution.solutionInfo)[@"time stepping method"];
        found = YES;
    }
    if ([method isEqualToString:@"none"] == YES) return;
    
    if (found == NO) {
        found = NO;
        if ((solution.solutionInfo)[@"newmark beta"] != nil) {
            solution.beta = [(solution.solutionInfo)[@"newmark beta"] doubleValue];
            found = YES;
        } else {
            solution.beta = [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:@"newmark beta" info:&found minValue:NULL maxValue:NULL];
        }
        
        if (found == NO) {
            NSLog(@"FEMCore:FEMCore_initializeTimeStepInSolution: time stepping default to IMPLICIT EULER.");
            solution.beta = 1.0;
            method = @"implicit euler";
        }
    } else {
        solution.beta = 1.0;
        if ([method isEqualToString:@"implicit euler"] == YES) {
            solution.beta = 1.0;
        } else if ([method isEqualToString:@"explicit euler"] == YES) {
            solution.beta = 0.0;
        } else if ([method isEqualToString:@"runge-kutta"] == YES) {
            solution.beta = 0.0;
        } else if ([method isEqualToString:@"crank-nicolson"] == YES) {
            solution.beta = 0.5;
        } else if ([method isEqualToString:@"fs"] == YES) {
            solution.beta = 0.5;
        } else if ([method isEqualToString:@"newmark"] == YES) {
            if ((solution.solutionInfo)[@"newmark beta"] != nil) {
                solution.beta = [(solution.solutionInfo)[@"newmark beta"] doubleValue];
            } else {
                solution.beta = [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:@"newmark beta" info:&found minValue:NULL maxValue:NULL];
            }
            
            if (solution.beta < 0 || solution.beta > 1) {
                NSLog(@"FEMCore:FEMCore_initializeTimeStepInSolution: invalid value of beta: %f.\n", solution.beta);
            }
        } else if ([method isEqualToString:@"bdf"] == YES) {
            if (solution.order < 1 || solution.order > 5) {
                NSLog(@"FEMCore:FEMCore_initializeTimeStepInSolution: invalid BDF order: %d.\n", solution.order);
            }
        } else {
            NSLog(@"FEMCore:FEMCore_initializeTimeStepInSolution: unknown time stepping method: %@.\n", method);
            fatal("FEMCore:FEMCore_initializeTimeStepInSolution", "Saino will abort the simulation now...");
        }
    }
    
    if ([method isEqualToString:@"bdf"] == NO || solution.timeOrder > 1) {
        if (solution.doneTime == 1 && solution.beta != 0.0) {
            solution.beta = 1.0;
        }
        
        switch (solution.timeOrder) {
            case 1:
                order = min(solution.doneTime, solution.order);
                for (i=order-1; i>=1; i--) {
                    for (j=0; j<varContainers->size1PrevValues; j++) {
                        varContainers->PrevValues[j][i] = varContainers->PrevValues[j][i-1];
                    }
                }
                for (i=0; i<varContainers->size1PrevValues; i++) {
                    varContainers->PrevValues[i][0] = varContainers->Values[i];
                }
                matContainers = solution.matrix.getContainers;
                for (i=0; i<matContainers->size1force; i++) {
                    matContainers->Force[i][1] = matContainers->Force[i][0];
                }
                break;
            case 2:
                found = NO;
                if ((solution.solutionInfo)[@"bossak alpha"] != nil) {
                    solution.alpha = [(solution.solutionInfo)[@"bossak alpha"] doubleValue];
                    found = YES;
                } else {
                    solution.alpha = [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:@"bossak alpha" info:&found minValue:NULL maxValue:NULL];
                }
                if (found == NO) solution.alpha = -0.05;
                
                for (i=0; i<varContainers->size1PrevValues; i++) {
                    varContainers->PrevValues[i][2] = varContainers->Values[i];
                    varContainers->PrevValues[i][3] = varContainers->PrevValues[i][0];
                    varContainers->PrevValues[i][4] = varContainers->PrevValues[i][1];
                }
                break;
        }
    } else {
        order = min(solution.doneTime, solution.order);
        for (i=order-1; i>=1; i--) {
            for (j=0; j<varContainers->size1PrevValues; j++) {
                varContainers->PrevValues[j][i] = varContainers->PrevValues[j][i-1];
            }
        }
        for (i=0; i<varContainers->size1PrevValues; i++) {
            varContainers->PrevValues[i][0] = varContainers->Values[i];
        }
    }
}

/*************************************************************************
    Initialize matrix structure and vector to zero intial value
 
    (FEMMatrix *)matrix    ->  Matrix to be initialized
    (double *)forceVector  ->  Vector to be initialized
*************************************************************************/
-(void)initializeToZeroMatrix:(FEMMatrix * __nonnull)matrix forceVector:(double * __nonnull)forceVector sizeForceVector:(int)sizeForceVector model:(FEMModel * __nonnull)model solution:(FEMSolution * __nonnull)solution {
    
    int i, dim;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = matrix.getContainers;
    
    if (matrix != nil) {
        memset( matContainers->Values, 0.0, matContainers->sizeValues*sizeof(double) );
        
        if (matContainers->MassValues != NULL) {
            memset( matContainers->MassValues, 0.0, matContainers->sizeMassValues*sizeof(double) );
        }
        if (matContainers->DampValues != NULL) {
            memset( matContainers->DampValues, 0.0, matContainers->sizeDampValues*sizeof(double) );
        }
        
        if (matContainers->Force != NULL) {
            for (i=0; i<matContainers->size1force; i++) {
                matContainers->Force[i][0] = 0.0;
            }
        }
        if (matContainers->BulkRHS != NULL) {
            memset( matContainers->BulkRHS, 0.0, matContainers->sizeBulkRHS*sizeof(double) );
        }
        if (matContainers->RHS_im != NULL) {
            memset( matContainers->RHS_im, 0.0, matContainers->sizeRHS*sizeof(double) );
        }
        if (matContainers->BulkValues != NULL) {
            memset( matContainers->BulkValues, 0.0, matContainers->sizeBulkValues*sizeof(double) );
        }
    }
    
    memset(forceVector, 0.0, sizeForceVector*sizeof(double) );
    
    self.normalTangentialNumberOfNodes = 0;
    if (solution.variable.dofs <= 1) return;
    
    [self.normalTangentialName setString:@"normal-tangential"];
    if ([[solution.variable canonicalizeName] isEqualToString:@"flow solution"] == YES) {
        [self.normalTangentialName appendString:@" velocity"];
    } else {
        [self.normalTangentialName appendString:@" "];
        [self.normalTangentialName appendString:[solution.variable canonicalizeName]];
    }
    
    dim = model.dimension;
    [self FEMCore_checkNormalTangentialBoundaryModel:model variableName:self.normalTangentialName dimension:dim];
    [self FEMCore_averageBoundaryNormalsModel:model variableName:self.normalTangentialName dimension:dim];
}

/*******************************************************************************
    Performs initialization for matrix equation related to the active solution
*******************************************************************************/
-(void)defaultInitializeSolution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model {
    
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = solution.matrix.getContainers;
    [self initializeToZeroMatrix:solution.matrix
                     forceVector:matContainers->RHS sizeForceVector:matContainers->sizeRHS model:model solution:solution];
}

/*****************************************************************************
    Returns a real by its name if found in the array structure and in the
    active element
*****************************************************************************/
-(BOOL)getReal:(FEMModel * __nonnull)model forElement:(Element_t * __nullable)element inArray:(NSArray * __nonnull)array variableName:(NSString * __nonnull)name buffer:(listBuffer * __nonnull)result listUtilities:(FEMListUtilities * __nonnull)listUtil {
    
    int n;
    int dNodes[1];
    int *nodeIndexes = NULL;
    BOOL found;
    
    if (element != NULL) {
        n = [self getNumberOfNodesForElement:element];
        nodeIndexes = element->NodeIndexes;
    } else {
        n = 1;
        nodeIndexes = dNodes;
        nodeIndexes[0] = 0;
    }
    
    found = [listUtil listGetReal:model inArray:array forVariable:name numberOfNodes:n indexes:element->NodeIndexes buffer:result minValue:NULL maxValue:NULL];
    
    return found;
}

-(int)isPElement:(Element_t * __nonnull)element {
    
    if (element->Pdefs != NULL) {
        return 1;
    } else {
        return 0;
    }
}

-(void)getNodes:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model inElement:(Element_t * __nonnull)element resultNodes:(Nodes_t * __nonnull)nodes numberOfNodes:(int * __nullable)nd mesh:(FEMMesh * __nullable)mesh {
    
    int i, n, nb, sz, sz1;
    
    FEMMesh *msh = solution.mesh;
    if (mesh != nil) msh = mesh;
    n = max(msh.maxElementNodes, msh.maxElementDofs);
    
    Nodes_t *meshNodes = msh.getNodes;
    
    if (nodes->x == NULL) {
        nodes->x = doublevec(0, n-1);
        nodes->y = doublevec(0, n-1);
        nodes->z = doublevec(0, n-1);
        nodes->numberOfNodes = n;
    } else if (nd != NULL && *nd < n) {
        free_dvector(nodes->x, 0, *nd-1);
        free_dvector(nodes->y, 0, *nd-1);
        free_dvector(nodes->z, 0, *nd-1);
        nodes->x = doublevec(0, n-1);
        nodes->y = doublevec(0, n-1);
        nodes->z = doublevec(0, n-1);
        nodes->numberOfNodes = n;
    }
    
    n = element->Type.NumberOfNodes;
    for (i=0; i<n; i++) {
        nodes->x[i] = meshNodes->x[element->NodeIndexes[i]];
        nodes->y[i] = meshNodes->y[element->NodeIndexes[i]];
        nodes->z[i] = meshNodes->z[element->NodeIndexes[i]];
    }
    
    sz = max(solution.mesh.maxElementNodes, solution.mesh.maxElementDofs); 
    if ( sz > n) {
        for (i=n; i<sz; i++) {
            nodes->x[i] = 0.0;
            nodes->y[i] = 0.0;
            nodes->z[i] = 0.0;
        }
    }
    
    sz1 = meshNodes->numberOfNodes;
    if (sz1 > solution.mesh.numberOfNodes) {
        memset( self.indexStore, -1, self.sizeIndexStore*sizeof(int) );
        nb = [self getElementDofsSolution:solution model:model forElement:element atIndexes:self.indexStore disableDiscontinuousGalerkin:NULL];
        for (i=n; i<nb; i++) {
            if (self.indexStore[i] >= 0 && self.indexStore[i] < sz1) {
                nodes->x[i] = meshNodes->x[self.indexStore[i]];
                nodes->y[i] = meshNodes->y[self.indexStore[i]];
                nodes->z[i] = meshNodes->z[self.indexStore[i]];
            }
        }
    }
}

/*****************************************************************************
    Return true if the element is a possible flux element
    Needed to skip nodal elements in 2D and 3D boundary condition setting.
*****************************************************************************/
-(BOOL)isFluxElement:(Element_t * __nonnull)element mesh:(FEMMesh * __nonnull)mesh {
    
    int family, meshDim;
    
    meshDim = mesh.dimension;
    family = [self getElementFamily:element];
    
    // This is not generic rule but happens to be true for all combinations
    // 3D: families 3 and 4
    // 2D: family 2
    // 1D: family 1
    return (meshDim <= family) ? YES : NO;
}

-(int)getElementFamily:(Element_t * __nonnull)element {
    
    return element->Type.ElementCode / 100;
}

/************************************************************
    Method corresponds to Elmer from git on October 27 2015
************************************************************/
-(int)getElementDofsSolution:(FEMSolution * __nullable)uSolution model:(FEMModel * __nonnull)model forElement:(Element_t * __nonnull)element atIndexes:(int * __nonnull)indexes disableDiscontinuousGalerkin:(BOOL * __nullable)disableDiscontinuousGalerkin {
    
    int nb, i, j, k, bid, edofs, fdofs, faceDofs, edgeDofs, bubbleDofs, ind;
    BOOL gb;
    FEMSolution *solution;
    Element_t *parent, *edges = NULL, *faces = NULL;
    solutionArraysContainer *solContainers = NULL;
    
    if (uSolution != nil) {
        solution = uSolution;
    } else {
        solution = (FEMSolution *)model.solution;
    }
    
    edges = solution.mesh.getEdges;
    faces = solution.mesh.getFaces;
    
    nb = 0;
    
    BOOL disableDG = NO;
    if (disableDiscontinuousGalerkin != NULL) disableDG = *disableDiscontinuousGalerkin;
    if (!disableDG && (solution.solutionInfo)[@"discontinuous galerkin"] != nil) {
        if ([(solution.solutionInfo)[@"discontinuous galerkin"] boolValue] == YES) {
            
            for (i=0; i<element->DGDOFs; i++) {
                indexes[nb] = element->DGIndexes[i];
                nb++;
            }
            
            if (element->BoundaryInfo != NULL) {
                if (element->BoundaryInfo->Left != NULL) {
                    for (i=0; i<element->BoundaryInfo->Left->DGDOFs; i++) {
                        indexes[nb] = element->BoundaryInfo->Left->DGIndexes[i];
                        nb++;
                    }
                }
                if (element->BoundaryInfo->Right != NULL) {
                    for (i=0; i<element->BoundaryInfo->Right->DGDOFs; i++) {
                        indexes[nb] = element->BoundaryInfo->Right->DGIndexes[i];
                        nb++;
                    }
                }
            }
            
            if (nb > 0) { 
                return nb;
            }
        }
    }
    
    bid = element->BodyID;
    if (bid == 0 && element->BoundaryInfo != NULL) {
        if (element->BoundaryInfo->Left != NULL) bid = element->BoundaryInfo->Left->BodyID;
        if (element->BoundaryInfo->Right != NULL) bid = element->BoundaryInfo->Right->BodyID;
    }
    if (bid == 0) bid = 1;
    
    solContainers = solution.getContainers;
    
    if (solContainers->defDofs[bid-1][0] > 0) {
        for (i=0; i<element->NDOFs; i++) {
            indexes[nb] = element->NodeIndexes[i];
            nb++;
        }
    }
    for (i=1; i<solContainers->size2DefDofs; i++) {
        if (solContainers->defDofs[bid-1][i] < 0) {
            return nb;
        }
    }
    
    faceDofs = solution.mesh.maxFaceDofs;
    edgeDofs = solution.mesh.maxEdgeDofs;
    bubbleDofs = solution.mesh.maxBdofs;
    
    if (element->EdgeIndexes != NULL) {
        for (j=0; j<element->Type.NumberOfEdges; j++) {
            edofs = edges[element->EdgeIndexes[j]].BDOFs;
            for (i=0; i<edofs; i++) {
                indexes[nb] = edgeDofs*(element->EdgeIndexes[j]) + i + solution.mesh.numberOfNodes;
                nb++;
            }
        }
    }
    
    if (element->FaceIndexes != NULL) {
        for (j=0; j<element->Type.NumberOfFaces; j++) {
            fdofs = faces[element->FaceIndexes[j]].BDOFs;
            for (i=0; i<fdofs; i++) {
                indexes[nb] = faceDofs*(element->FaceIndexes[j]) + i + solution.mesh.numberOfNodes +
                edgeDofs*solution.mesh.numberOfEdges;
                nb++;
            }
        }
    }
    
    if ((solution.solutionInfo)[@"bubbles in global system"] != nil) {
        gb = [(solution.solutionInfo)[@"bubbles in global system"] boolValue];
    } else {
        gb = YES;
    }
    
    if ( element->BoundaryInfo != NULL && (element->EdgeIndexes == NULL || element->FaceIndexes == NULL) ) {
        
        parent = element->BoundaryInfo->Left;
        if (parent == NULL) parent = element->BoundaryInfo->Right;
        if (parent == NULL) {
            return nb;
        }
        
        if ([self getElementFamily:element] == 2) {
            if (parent->EdgeIndexes != NULL) {
                if ([self isPElement:element] == 1) {
                    ind = element->Pdefs->LocalNumber;
                } else {
                    for (ind=0; ind<parent->Type.NumberOfEdges; ind++) {
                        k = 0;
                        for (i=0; i<edges[parent->EdgeIndexes[ind]].Type.NumberOfNodes; i++) {
                            for (j=0; j<element->Type.NumberOfNodes; j++) {
                                if (edges[parent->EdgeIndexes[ind]].NodeIndexes[i] == element->NodeIndexes[j]) k++;
                            }
                        }
                        if (k == element->Type.NumberOfNodes) break;
                    }
                }
                
                edofs = element->BDOFs;
                for (i=0; i<edofs; i++) {
                    indexes[nb] = edgeDofs*(parent->EdgeIndexes[ind]) + i + solution.mesh.numberOfNodes;
                    nb++;
                }
            }
        } 
        else if ([self getElementFamily:element] == 3 || [self getElementFamily:element] == 4) {
            if (parent->FaceIndexes != NULL) {
                if ([self isPElement:element] == 1) {
                    ind = element->Pdefs->LocalNumber;
                } else {
                    for (ind=0; ind<parent->Type.NumberOfFaces; ind++) {
                        k = 0;
                        for (i=0; i<faces[parent->FaceIndexes[ind]].Type.NumberOfNodes; i++) {
                            for (j=0; j<element->Type.NumberOfNodes; j++) {
                                if (faces[parent->FaceIndexes[ind]].NodeIndexes[i] == element->NodeIndexes[j]) k++;
                            }
                        }
                        if (k == element->Type.NumberOfNodes) break;
                    }
                }
                
                fdofs = element->BDOFs;
                for (i=0; i<fdofs; i++) {
                    indexes[nb] = faceDofs*(parent->FaceIndexes[ind]) + i + solution.mesh.numberOfNodes
                    + edgeDofs*solution.mesh.numberOfEdges;
                    nb++;
                }
            }
        }
    } else if (gb == YES) {
        if (element->BubbleIndexes != NULL) {
            for (i=0; i<element->BDOFs; i++) {
                indexes[nb] = faceDofs*solution.mesh.numberOfFaces + solution.mesh.numberOfNodes + edgeDofs*solution.mesh.numberOfEdges + element->BubbleIndexes[i];
                nb++;
            }
        }
    }
    
    return nb;
}

-(int)sgetElementDofsSolution:(FEMSolution * __nullable)uSolution model:(FEMModel * __nonnull)model forElement:(Element_t * __nonnull)element atIndexes:(int * __nonnull)indexes {
    
    int nb, i, j, edofs, fdofs, faceDofs, edgeDofs, bubbleDofs;
    BOOL gb;
    FEMSolution *solution;
    Element_t *parent, *edges = NULL, *faces = NULL;
    
    if (uSolution != nil) {
        solution = uSolution;
    } else {
        solution = (FEMSolution *)model.solution;
    }
    
    edges = solution.mesh.getEdges;
    faces = solution.mesh.getFaces;
    
    nb = 0;
    
    if ((solution.solutionInfo)[@"discontinuous galerkin"] != nil) {
        if ([(solution.solutionInfo)[@"discontinuous galerkin"] boolValue] == YES) {
            
            for (i=0; i<element->DGDOFs; i++) {
                indexes[nb] = element->DGIndexes[i];
                nb++;
            }
            
            if (element->BoundaryInfo != NULL) {
                if (element->BoundaryInfo->Left != NULL) {
                    for (i=0; i<element->BoundaryInfo->Left->DGDOFs; i++) {
                        indexes[nb] = element->BoundaryInfo->Left->DGIndexes[i];
                        nb++;
                    }
                }
                if (element->BoundaryInfo->Right != NULL) {
                    for (i=0; i<element->BoundaryInfo->Right->DGDOFs; i++) {
                        indexes[nb] = element->BoundaryInfo->Right->DGIndexes[i];
                        nb++;
                    }
                }
            }
            
            if (nb > 0) return nb;
        }
    }
    
    for (i=0; i<element->NDOFs; i++) {
        indexes[nb] = element->NodeIndexes[i];
        nb++;
    }
    
    faceDofs = solution.mesh.maxFaceDofs;
    edgeDofs = solution.mesh.maxEdgeDofs;
    bubbleDofs = solution.mesh.maxBdofs;
    
    if (element->EdgeIndexes != NULL) {
        for (j=0; j<element->Type.NumberOfEdges; j++) {
            edofs = edges[element->EdgeIndexes[j]].BDOFs;
            for (i=0; i<edofs; i++) {
                indexes[nb] = edgeDofs*element->EdgeIndexes[j] + i + solution.mesh.numberOfNodes;
                nb++;
            }
        }
    }
    
    if (element->FaceIndexes != NULL) {
        for (j=0; j<element->Type.NumberOfFaces; j++) {
            fdofs = faces[element->FaceIndexes[j]].BDOFs;
            for (i=0; i<fdofs; i++) {
                indexes[nb] = faceDofs*element->FaceIndexes[j] + i + solution.mesh.numberOfNodes +
                edgeDofs*solution.mesh.numberOfEdges;
                nb++;
            }
        }
    }
    
    if ((solution.solutionInfo)[@"bubbles in global system"] != nil) {
        gb = [(solution.solutionInfo)[@"bubbles in global system"] boolValue];
    } else {
        gb = YES;
    }
    
    if (element->BoundaryInfo != NULL) {
        
        if ([self isPElement:element] == NO) return nb;
        
        parent = element->BoundaryInfo->Left;
        if (parent == NULL) parent = element->BoundaryInfo->Right;
        if (parent == NULL) return nb;
        
        if (parent->EdgeIndexes != NULL) {
            edofs = element->BDOFs;
            for (i=0; i<edofs; i++) {
                indexes[nb] = edgeDofs*parent->EdgeIndexes[element->Pdefs->LocalNumber] + i + solution.mesh.numberOfNodes;
                nb++;
            }
        }
        
        if (parent->FaceIndexes != NULL) {
            fdofs = element->BDOFs;
            for (i=0; i<fdofs; i++) {
                indexes[nb] = faceDofs*(parent->FaceIndexes[element->Pdefs->LocalNumber]) + i + solution.mesh.numberOfNodes + edgeDofs*solution.mesh.numberOfEdges;
                nb++;
            }
        }
    } else if (gb == YES) {
        if (element->BubbleIndexes != NULL) {
            for (i=0; i<element->BDOFs; i++) {
                indexes[nb] = faceDofs*solution.mesh.numberOfFaces + solution.mesh.numberOfNodes + edgeDofs*solution.mesh.numberOfEdges + element->BubbleIndexes[i];
                nb++;
            }
        }
        
    }
    
    return nb;    
}

// Return the number of bubble degree of freedom in the active element
-(int)getNumberOfBubbleDofsElement:(Element_t * __nonnull)element solution:(FEMSolution * __nonnull)solution {
    
    int n=0;
    BOOL gb;
    
    if ((solution.solutionInfo)[@"bubbles in global system"] !=nil) {
        gb = [(solution.solutionInfo)[@"bubbles in global system"] boolValue];
    } else {
        gb = YES;
    }
    
    if (gb == NO) {
        n = element->BDOFs;
    }
    
    return n;
}

-(Element_t * __nonnull)getActiveElement:(int)t solution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model {
    
    Element_t *element = NULL, *elements = NULL, *currentElement = NULL;
    solutionArraysContainer *solutionContainers = NULL;
    
    if (t >= 0 && t < solution.numberOfActiveElements) {
        solutionContainers = solution.getContainers;
        elements = solution.mesh.getElements;
        element = &elements[solutionContainers->activeElements[t]];
        currentElement = model.getCurrentElement;
        currentElement = element;
    } else {
        NSLog(@"FEMCore:getActiveElement: invalid element number requested: %d\n", t);
        fatal("FEMCore:getActiveElement", "Saino will abort the simulation now...");
    }
    return element;
}

-(int)getNumberOfNodesForElement:(Element_t * __nonnull)element {
    
    return element->Type.NumberOfNodes;
}

-(int)getBoundaryConditionID:(FEMModel * __nonnull)model forElement:(Element_t * __nonnull)element {
    
    int bc_id = 0;
    
    for (FEMBoundaryCondition *boundaryCondition in model.boundaryConditions) {
        bc_id++;
        if (element->BoundaryInfo->Constraint == boundaryCondition.tag) break;
    }
    
    if (bc_id > model.numberOfBoundaryConditions) bc_id = 0;
    
    return bc_id;
}

-(NSArray * __nullable)getBoundaryCondition:(FEMModel * __nonnull)model forElement:(Element_t * __nonnull)element {
    
    int bc_id;
    
    FEMBoundaryCondition *boundaryConditionAtId;
    
    // Returns a boundaryCondition object at index bc_id
    bc_id = [self getBoundaryConditionID:model forElement:element];
    if (bc_id > 0) boundaryConditionAtId = (model.boundaryConditions)[bc_id-1];
    
    if (boundaryConditionAtId != nil) {
        return boundaryConditionAtId.valuesList;
    } else {
        return nil;
    }
}

-(Element_t * __nonnull)getBoundaryElement:(FEMSolution * __nonnull)solution atIndex:(int)index {
    
    Element_t *elements;
    
    elements = solution.mesh.getElements;
    
    if (index < 0 || index > solution.mesh.numberOfBoundaryElements-1) {
        NSLog(@"FEMCore:getBoundaryElement: invalid element number requested at index: %d.\n", index);
    }
    
    return &elements[solution.mesh.numberOfBulkElements+index];
}

-(void)getScalarLocalField:(double * __nonnull)field sizeField:(int)sizeField name:(NSString * __nullable)name element:(Element_t * __nonnull)element solution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model timeStep:(int * __nullable)tStep {
    
    int i, j, n;
    double *values = NULL;
    BOOL found, needsDeallocation=NO;
    FEMVariable *variable = nil;
    FEMSolution *varSolution;
    FEMUtilities *utilities;
    variableArraysContainer *varContainers = NULL;
    
    utilities = [[FEMUtilities alloc] init];
    memset( field, 0.0, sizeField*sizeof(double) );
    
    variable = solution.variable;
    if (name != nil) {
        variable = [utilities getVariableFrom:solution.mesh.variables model:model name:name onlySearch:NULL maskName:NULL info:&found];
    }
    if (variable == nil) return;
    
    memset( self.indexStore, -1, self.sizeIndexStore*sizeof(int) );
    varSolution = (FEMSolution *)variable.solution;
    if (varSolution != nil) {
        n = [self getElementDofsSolution:varSolution model:model forElement:element atIndexes:self.indexStore disableDiscontinuousGalerkin:NULL];
    } else {
        n = [self getElementDofsSolution:solution model:model forElement:element atIndexes:self.indexStore disableDiscontinuousGalerkin:NULL];
    }
    n = min(n, sizeField);
    
    varContainers = variable.getContainers;
    if (tStep != NULL) {
        if (*tStep < 0) {
            if (varContainers->PrevValues != NULL && -*tStep < varContainers->size2PrevValues) {
                values = doublevec(0, varContainers->size1PrevValues-1);
                for (i=0; i<varContainers->size1PrevValues; i++) {
                    values[i] = varContainers->PrevValues[i][-*tStep];
                }
                needsDeallocation = YES;
            }
        }
    } else values = varContainers->Values;
    
    if (varContainers->Perm != NULL) {
        for (i=0; i<n; i++) {
            j = self.indexStore[i];
            if (j >= 0 && j < varContainers->sizePerm) {
                j = varContainers->Perm[j];
                if (j >= 0) field[i] = values[j];
            }
        }
    } else {
        for (i=0; i<n; i++) {
            j = self.indexStore[i];
            if (j >= 0 && j < varContainers->sizeValues) field[i] = values[self.indexStore[i]];
        }
    }
    
    if (needsDeallocation == YES) free_dvector(values, 0, varContainers->size1PrevValues-1);
}

-(void)getVectorLocalField:(double * __nonnull * __nonnull)field size1Field:(int)size1Field size2Field:(int)size2Field name:(NSString * __nullable)name element:(Element_t * __nonnull)element solution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model timeStep:(int * __nullable)tStep {
    
    int i, j, k, n;
    double *values = NULL;
    BOOL found, needsDeallocation=NO;
    FEMVariable *variable = nil;
    FEMSolution *varSolution;
    FEMUtilities *utilities;
    variableArraysContainer *varContainers = NULL;
    
    utilities = [[FEMUtilities alloc] init];
    memset( *field, 0.0, (size1Field*size2Field)*sizeof(double) );
    
    variable = solution.variable;
    if (name != nil) {
        variable = [utilities getVariableFrom:solution.mesh.variables model:model name:name onlySearch:NULL maskName:NULL info:&found];
    }
    if (variable == nil) return;
    
    if (variable != nil) {
        memset( self.indexStore, -1, self.sizeIndexStore*sizeof(int) );
        varSolution = (FEMSolution *)variable.solution;
        if (varSolution != nil) {
            n = [self getElementDofsSolution:varSolution model:model forElement:element atIndexes:self.indexStore disableDiscontinuousGalerkin:NULL];
        } else {
            n = [self getElementDofsSolution:solution model:model forElement:element atIndexes:self.indexStore disableDiscontinuousGalerkin:NULL];
        }
        n = min(n, (size1Field*size2Field));
        
        varContainers = variable.getContainers;
        if (tStep != NULL) {
            if (*tStep < 0) {
                if (varContainers->PrevValues != NULL && -*tStep < varContainers->size2PrevValues) {
                    values = doublevec(0, varContainers->size1PrevValues-1);
                    for (i=0; i<varContainers->size1PrevValues; i++) {
                        values[i] = varContainers->PrevValues[i][-*tStep];
                    }
                    needsDeallocation = YES;
                }
            }
        } else values = varContainers->Values;
        
        for (i=0; i<variable.dofs; i++) {
            if (varContainers->Perm != NULL) {
                for (j=0; j<n; j++) {
                    k = self.indexStore[j];
                    if (k >= 0 && k < varContainers->sizePerm) {
                        k = varContainers->Perm[k];
                        if (k >= 0) field[i][j] = values[variable.dofs*k+i];
                    }
                }
            } else {
                for (j=0; j<n; j++) {
                    if (variable.dofs*self.indexStore[j]+i < varContainers->sizeValues) {
                        field[i][j] = values[variable.dofs*self.indexStore[j]+i];
                    }
                }
            }
        }
    }
    
    if (needsDeallocation == YES) free_dvector(values, 0, varContainers->size1PrevValues-1);    
}

/****************************************************
    Return the EdgeMap and its size (number of rows)
****************************************************/
-(int * __nonnull * __nonnull)getEdgeMap:(int)elementFamily mapSize:(int * __nullable)mapSize {
    
    int **edgeMap = NULL;
    
    switch (elementFamily) {
        case 2:
            edgeMap = _lineEM;
            break;
        case 3:
            edgeMap = _triangleEM;
            break;
        case 4:
            edgeMap = _quadEM;
            break;
        case 5:
            edgeMap = _tetraEM;
            break;
        case 6:
            edgeMap = _prismEM;
            break;
        case 7:
            edgeMap = _wedgeEM;
            break;
        case 8:
            edgeMap = _brickEM;
            break;
    }
    
    if (mapSize != NULL) {
        switch (elementFamily) {
            case 2:
                *mapSize = 1;
                break;
            case 3:
                *mapSize = 3;
                break;
            case 4:
                *mapSize = 4;
                break;
            case 5:
                *mapSize = 6;
                break;
            case 6:
                *mapSize = 8;
                break;
            case 7:
                *mapSize = 9;
                break;
            case 8:
                *mapSize = 12;
                break;
        }
    }
    
    if (_initialized[elementFamily-1] == NO) {
        _initialized[elementFamily-1] = YES;
        switch (elementFamily) {
            case 2:
                edgeMap[0][0] = 0;
                edgeMap[0][1] = 1;
                break;
            case 3:
                edgeMap[0][0] = 0;
                edgeMap[0][1] = 1;
                
                edgeMap[1][0] = 1;
                edgeMap[1][1] = 2;
                
                edgeMap[2][0] = 2;
                edgeMap[2][1] = 0;
                break;
            case 4:
                edgeMap[0][0] = 0;
                edgeMap[0][1] = 1;
                
                edgeMap[1][0] = 1;
                edgeMap[1][1] = 2;
                
                edgeMap[2][0] = 2;
                edgeMap[2][1] = 3;
                
                edgeMap[3][0] = 3;
                edgeMap[3][1] = 0;
                break;
            case 5:
                edgeMap[0][0] = 0;
                edgeMap[0][1] = 1;
                
                edgeMap[1][0] = 1;
                edgeMap[1][1] = 2;
                
                edgeMap[2][0] = 2;
                edgeMap[2][1] = 0;
                
                edgeMap[3][0] = 0;
                edgeMap[3][1] = 3;
                
                edgeMap[4][0] = 1;
                edgeMap[4][1] = 3;
                
                edgeMap[5][0] = 2;
                edgeMap[5][1] = 3;
                break;
            case 6:
                edgeMap[0][0] = 0;
                edgeMap[0][1] = 1;
                
                edgeMap[1][0] = 1;
                edgeMap[1][1] = 2;
                
                edgeMap[2][0] = 3;
                edgeMap[2][1] = 2;
                
                edgeMap[3][0] = 0;
                edgeMap[3][1] = 3;
                
                edgeMap[4][0] = 0;
                edgeMap[4][1] = 4;
                
                edgeMap[5][0] = 1;
                edgeMap[5][1] = 4;
                
                edgeMap[6][0] = 2;
                edgeMap[6][1] = 4;
                
                edgeMap[7][0] = 3;
                edgeMap[7][1] = 4;
                break;
            case 7:
                edgeMap[0][0] = 0;
                edgeMap[0][1] = 1;
                
                edgeMap[1][0] = 1;
                edgeMap[1][1] = 2;
                
                edgeMap[2][0] = 2;
                edgeMap[2][1] = 0;
                
                edgeMap[3][0] = 3;
                edgeMap[3][1] = 4;
                
                edgeMap[4][0] = 4;
                edgeMap[4][1] = 5;
                
                edgeMap[5][0] = 5;
                edgeMap[5][1] = 3;
                
                edgeMap[6][0] = 0;
                edgeMap[6][1] = 3;
                
                edgeMap[7][0] = 1;
                edgeMap[7][1] = 4;
                
                edgeMap[8][0] = 2;
                edgeMap[8][1] = 5;
                break;
            case 8:
                edgeMap[0][0] = 0;
                edgeMap[0][1] = 1;
                
                edgeMap[1][0] = 1;
                edgeMap[1][1] = 2;
                
                edgeMap[2][0] = 3;
                edgeMap[2][1] = 2;
                
                edgeMap[3][0] = 0;
                edgeMap[3][1] = 3;
                
                edgeMap[4][0] = 4;
                edgeMap[4][1] = 5;
                
                edgeMap[5][0] = 5;
                edgeMap[5][1] = 6;
                
                edgeMap[6][0] = 7;
                edgeMap[6][1] = 6;
                
                edgeMap[7][0] = 4;
                edgeMap[7][1] = 7;
                
                edgeMap[8][0] = 0;
                edgeMap[8][1] = 4;
                
                edgeMap[9][0] = 1;
                edgeMap[9][1] = 5;
                
                edgeMap[10][0] = 2;
                edgeMap[10][1] = 6;
                
                edgeMap[11][0] = 3;
                edgeMap[11][1] = 7;
                break;
        }
    }
    
    return edgeMap;
}

/***********************************************************************************************
    Calculate global indexes of boundary dofs for given element and its boundary
 
    Arguments:
 
        FEMMesh *mesh      -> class containing the mesh and its edges and faces of elements
        Element_t *element -> Boundary element to get indexes for
        Element_t *parent  -> Parent of boundary element
        int *indexes       -> Calculated indexes of boundary element in global system
    int resultSize     -> Size of created index vector, i.e., how many indexes were created
************************************************************************************************/
-(void)getBoundaryIndexes:(FEMMesh * __nonnull)mesh forBoundaryElement:(Element_t * __nonnull)element withParentElement:(Element_t * __nonnull)parent resultVector:(int * __nonnull)indexes sizeVector:(int)size indexSize:(int * __nonnull)indexSize {
    
    int i, j, n;
    Element_t *edges = NULL, *faces = NULL;
    
    // Clear indexes
    memset( indexes, 0, size*sizeof(int) );
    
    n = element->Type.NumberOfNodes;
    
    edges = mesh.getEdges;
    faces = mesh.getFaces;
    
    // Nodal indexes
    memcpy(indexes, element->NodeIndexes, n*sizeof(int));
    
    // Assign rest of indexes if necessary
    switch (parent->Type.dimension) {
        case 1:
            *indexSize = n;
        case 2:
            // Add index for each bubble dof in edge
            for (i=0; i<element->BDOFs; i++) {
                if (size < n) {
                    fatal("FEMCore:getBoundaryIndexes", "Not enough space reserved for indexes.");
                    return;
                }
                
                indexes[n] = mesh.numberOfNodes + parent->EdgeIndexes[element->Pdefs->LocalNumber] * mesh.maxEdgeDofs + i;
                n++;
            }
            *indexSize = n;
            break;
        case 3:
            // Add indexes of faces edges
            for (i=0; i<faces[parent->FaceIndexes[element->Pdefs->LocalNumber]].Type.NumberOfEdges; i++) {
                
                // If edge has no dofs jump to next edge
                if (edges[faces[parent->FaceIndexes[element->Pdefs->LocalNumber]].EdgeIndexes[i]].BDOFs <= 0) continue;
                
                for (j=0; j<edges[faces[parent->FaceIndexes[element->Pdefs->LocalNumber]].EdgeIndexes[i]].BDOFs; j++) {
                    
                    if (size < n) {
                        fatal("FEMCore:getBoundaryIndexes", "Not enough space reserved for indexes.");
                        return;
                    }
                    
                    indexes[n] = mesh.numberOfNodes + faces[parent->FaceIndexes[element->Pdefs->LocalNumber]].EdgeIndexes[i]*mesh.maxEdgeDofs + j;
                    n++;
                }
            }
            
            // Add indexes of faces bubbles
            for (i=0; i<faces[parent->FaceIndexes[element->Pdefs->LocalNumber]].BDOFs; i++) {
                if (size < n) {
                    fatal("FEMCore:getBoundaryIndexes", "Not enough space reserved for indexes.");
                    return;
                }
                
                indexes[n] = mesh.numberOfNodes + mesh.numberOfEdges * mesh.maxEdgeDofs + parent->FaceIndexes[element->Pdefs->LocalNumber] * mesh.maxFaceDofs + i;
                n++;
            }
            
            *indexSize = n;
            break;
        default:
            fatal("FEMCore:getBoundaryIndexes", "Unsupported dimension.");
            break;
    }
}

-(int)getBodyForceIDForElement:(Element_t * __nonnull)element model:(FEMModel * __nonnull)model {
    
    int body_id, bf_id;

    body_id = element->BodyID;
    bf_id = 1;
    if ((model.bodies)[body_id-1][@"body force"] != nil) {
        bf_id = [(model.bodies)[body_id-1][@"body force"] intValue];
        if (bf_id < 1) bf_id = 1;
        if (bf_id > model.numberOfBodyForces) bf_id = model.numberOfBodyForces;
    }
    return bf_id;
}

-(int)getMaterialIDForElement:(Element_t * __nonnull)element model:(FEMModel * __nonnull)model {
    
    int body_id, mat_id;
    
    body_id = element->BodyID;
    mat_id = 1;
    if ((model.bodies)[body_id-1][@"material"] != nil) {
        mat_id = [(model.bodies)[body_id-1][@"material"] intValue];
        if (mat_id < 1) mat_id = 1;
        if (mat_id > model.numberOfMaterials) mat_id =  model.numberOfMaterials;
    }
    return mat_id;
}

-(int)getEquationIDForElement:(Element_t * __nonnull)element model:(FEMModel * __nonnull)model {
    
    int body_id, eq_id;
    
    body_id = element->BodyID;
    eq_id = 1;
    if ((model.bodies)[body_id-1][@"equation"] != nil) {
        eq_id = [(model.bodies)[body_id-1][@"equation"] intValue];
        if (eq_id < 1) eq_id = 1;
        if (eq_id > model.numberOfEquations) eq_id =  model.numberOfEquations;
    }
    return eq_id;
}

/*******************************************************************************************
    Returns a material property from either of the parents of the current boundary element
*******************************************************************************************/
-(BOOL)getParentMaterialProperty:(NSString * __nonnull)name forElement:(Element_t * __nonnull)element parentElement:(Element_t * __nullable)parentElement model:(FEMModel * __nonnull)model listUtilities:(FEMListUtilities * __nonnull)listUtilities buffer:(listBuffer * __nonnull)result {
    
    int leftright, mat_id, n;
    BOOL gotIt, found;
    Element_t *parent = NULL;
    FEMMaterial *materialAtID;
    
    n = [self getNumberOfNodesForElement:element];
    
    gotIt = NO;
    for (leftright=1; leftright<=2; leftright++) {
        if (leftright == 1) {
            parent = element->BoundaryInfo->Left;
        } else {
            parent = element->BoundaryInfo->Right;
        }
        if (parent != NULL) {
            mat_id = [self getMaterialIDForElement:element model:model];
            materialAtID = (model.materials)[mat_id-1];
            if ([listUtilities listCheckPresentVariable:name inArray:materialAtID.valuesList] == YES) {
                found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:name numberOfNodes:n indexes:element->NodeIndexes buffer:result minValue:NULL maxValue:NULL];
                parentElement = parent;
                gotIt = YES;
                break;
            }
        }
    }
    
    if (gotIt == NO) {
        NSLog(@"FEMCore:getMaterialProperty: property %@ not found in either parents.\n", name);
    }
    return gotIt;
}

-(BOOL)isActiveBoundaryElement:(Element_t * __nonnull)element inSolution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model {
    
    BOOL active;
    int i, n;
    variableArraysContainer *varContainers = NULL;
    
    memset( self.indexStore, -1, self.sizeIndexStore*sizeof(int) );
    n = [self getElementDofsSolution:solution model:model forElement:element atIndexes:self.indexStore disableDiscontinuousGalerkin:NULL];
    if ([self isPElement:element] == YES) n = [self getNumberOfNodesForElement:element];
    varContainers = solution.variable.getContainers;
    
    active = YES;
    for (i=0; i<n; i++) {
        if (varContainers->Perm[self.indexStore[i]] < 0) {
            active = NO;
            break;
        }
    }
    return active;
}

/********************************************************************************************************************
 Check n-t node setting element
********************************************************************************************************************/
-(void)checkNormalTangential:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution forElementNumber:(int)elno numberofNodes:(int)n atIndexes:(int * __nonnull)indexes atBoundary:(int)bc variableName:(NSMutableString * __nonnull)name orderOfDofs:(int)dof activeCondition:(BOOL)conditional conditionName:(NSString * __nonnull)condName permutationOffset:(int)offset {
    
    int i, j, k, m, dim;
    listBuffer condition = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    double *rotvec;
    BOOL stat, all;
    FEMBoundaryCondition *boundaryConditionAtId;
    variableArraysContainer *varContainers = NULL;
    solutionArraysContainer *solContainers = NULL;
    
    rotvec = doublevec(0, 2);
    dim = model.dimension;
    
    boundaryConditionAtId = (model.boundaryConditions)[bc];
    
    if (dof < 0) return;
    
    solContainers = solution.getContainers;
    
    all = YES;
    for (i=0; i<n; i++) {
        if (self.boundaryReorder[indexes[i]] < 0) {
            continue;
        } else {
            all = NO;
            break;
        }
    }
    if (all == YES) return;
    
    varContainers = solution.variable.getContainers;
    
    FEMListUtilities *listUtil = [FEMListUtilities sharedListUtilities];
    if ([listUtil listCheckPresentVariable:name inArray:boundaryConditionAtId.valuesList] == NO) return;
    if ([listUtil listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:solution.normalTangentialName info:&stat] == YES) return;
    
    if (conditional == YES) {
        stat = [listUtil listGetReal:model inArray:boundaryConditionAtId.valuesList forVariable:condName numberOfNodes:n indexes:indexes buffer:&condition minValue:NULL maxValue:NULL];
        conditional = (conditional == YES && stat == YES) ? YES : NO; 
    }
    
    // Check for nodes belonging to n-t boundary getting set by other bcs
    for (j=0; j<n; j++) {
        if (conditional == YES && condition.vector[j] < 0.0) continue;
        k = varContainers->Perm[indexes[j]];
        if (k >= 0) {
            k = k + offset;
            m = self.boundaryReorder[indexes[j]];
            if (m >= 0) {
                memset( rotvec, 0.0, 3*sizeof(double) );
                rotvec[dof] = 1.0;
                [self FEMCore_rotateNTSystem:rotvec nodeNumber:indexes[j] model:model];
                for (k=0; k<dim; k++) {
                    if (fabs(rotvec[k]) > 1.0e-8) self.ntElement[m][k] = elno;
                }
            }
        }
    }
    
    if (condition.vector != NULL) {
        free_dvector(condition.vector, 0, condition.m-1);
    }
    free_dvector(rotvec, 0, 2);
}

-(void)zeroTheNumberOfRows:(int)n inSolutionMatrix:(FEMSolution * __nonnull)solution {
        
    if (solution.matrix.format == MATRIX_CRS) {
        FEMMatrixCRS *crsMatrix = [[FEMMatrixCRS alloc] init];
        [crsMatrix zeroRowInGlobal:solution numberOfRows:n];
        
    } else if (solution.matrix.format == MATRIX_LIST) {
        
        // TODO: implement the zeroRow method for list matrix.
        
    } else if (solution.matrix.format == MATRIX_BAND || solution.matrix.format == MATRIX_SBAND) {
        FEMMatrixBand *bandMatrix = [[FEMMatrixBand alloc] init];
        [bandMatrix zeroRowInGlobal:solution numberOfRows:n];
    }
}

-(void)setMatrixElementForSolution:(FEMSolution * __nonnull)solution atIndex:(int)i andIndex:(int)j value:(double)value {
        
    if (solution.matrix.format == MATRIX_CRS) {
        FEMMatrixCRS *crsMatrix = [[FEMMatrixCRS alloc] init];
        [crsMatrix setElementInGlobal:solution row:i col:j value:value];
        
    } else if (solution.matrix.format == MATRIX_LIST) {
        // TODO: implement the setMatrixElement method for list matrix.
        
    } else if (solution.matrix.format == MATRIX_BAND || solution.matrix.format == MATRIX_SBAND) {
        FEMMatrixBand *bandMatrix = [[FEMMatrixBand alloc] init];
        [bandMatrix setElementInGlobal:solution row:i col:j value:value];
    }
}

-(void)addToMatrixElementForSolution:(FEMSolution * __nonnull)solution atIndex:(int)i andIndex:(int)j value:(double)value {
        
    if (solution.matrix.format == MATRIX_CRS) {
        FEMMatrixCRS *crsMatrix = [[FEMMatrixCRS alloc] init];
        [crsMatrix addToElementInGlobal:solution row:i col:j value:value];
        
    } else if (solution.matrix.format == MATRIX_LIST) {
        // TODO: implement the setMatrixElement method for list matrix.
        
    } else if (solution.matrix.format == MATRIX_BAND || solution.matrix.format == MATRIX_SBAND) {
        FEMMatrixBand *bandMatrix = [[FEMMatrixBand alloc] init];
        [bandMatrix addToElementInGlobal:solution row:i col:j value:value];
    }
}

/***************************************************************************************************************
    Integral given function over element j
    
    Arguments:
 
        NSArray *bc           ->   boundary condition value list
        Element_t *element    ->   boundary element to get stiffness matrix to
        int nd                ->   number of boundary element nodes
        NSString *name        ->   name of boundary condition
        double *integral      ->   integral of the given function
        
***************************************************************************************************************/
-(void)localBoundaryIntegral:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution atBoundary:(NSArray * __nonnull)bc forElement:(Element_t * __nonnull)element withNumberOfNodes:(int)nd andParent:(Element_t * __nonnull)parent withNumberOfNodes:(int)np boundaryName:(NSString * __nonnull)name functionIntegral:(double * __nonnull)integral {
    
    int i, n, jj=-1, kk=-1, t, size;
    int **edgeMap = NULL;
    double s, l, sum;
    double **vLoad, g[3], *vl;
    listBuffer load = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    listBuffer buffer = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    ElementType_t savedType;
    GaussIntegrationPoints *IP = NULL;
    Nodes_t *nodes, *pNodes;
    NSString *string;
    BOOL stat;
    
    FEMListUtilities *listUtilities = [FEMListUtilities sharedListUtilities];
    FEMElementDescription *elementDescription = [FEMElementDescription sharedElementDescription];
    FEMNumericIntegration *numericIntegration = [[FEMNumericIntegration alloc] init];
    if ([numericIntegration allocation:solution.mesh] == NO) fatal("FEMCore:localBoundaryIntegral", "Allocation error in FEMNumericIntegration.");
    
    n = max(solution.mesh.maxElementNodes, solution.mesh.maxElementDofs);
    
    nodes = (Nodes_t*)malloc(sizeof(Nodes_t));
    initNodes(nodes);
    nodes->x = doublevec(0, n-1);
    nodes->y = doublevec(0, n-1);
    nodes->z = doublevec(0, n-1);
    
    pNodes = (Nodes_t*)malloc(sizeof(Nodes_t));
    initNodes(pNodes);
    pNodes->x = doublevec(0, n-1);
    pNodes->y = doublevec(0, n-1);
    pNodes->z = doublevec(0, n-1);
    
    [self getNodes:solution model:model inElement:element resultNodes:nodes numberOfNodes:&n mesh:nil];
    [self getNodes:solution model:model inElement:parent resultNodes:pNodes numberOfNodes:&n mesh:nil];
    
    vLoad = doublematrix(0, 2, 0, np-1);
    memset( *vLoad, 0.0, (3*np)*sizeof(double) );
    
    vl = doublevec(0, 2);
    memset( vl, 0.0, 3*sizeof(double) );
    
    memset( g, 0.0, sizeof(g) );
    
    stat = [self getReal:model forElement:element inArray:bc variableName:name buffer:&load listUtilities:listUtilities];
    
    edgeMap = [self getEdgeMap:[self getElementFamily:parent] mapSize:&size];
    
    for (i=0; i<size; i++) {
        jj = edgeMap[i][0];
        kk = edgeMap[i][1];
        if ( (parent->NodeIndexes[jj] == element->NodeIndexes[0]
              && parent->NodeIndexes[kk] == element->NodeIndexes[1])
              || (parent->NodeIndexes[jj] == element->NodeIndexes[1]
              && parent->NodeIndexes[kk] == element->NodeIndexes[0]) ) break;
    }
    
    string = [name stringByAppendingString:@" 1"];
    stat = [self getReal:model forElement:element inArray:bc variableName:string buffer:&buffer listUtilities:listUtilities];
    for (i=0; i<nd; i++) {
        vLoad[0][i] = buffer.vector[i];
    }
    
    string = [name stringByAppendingString:@" 2"];
    stat = [self getReal:model forElement:element inArray:bc variableName:string buffer:&buffer listUtilities:listUtilities];
    for (i=0; i<nd; i++) {
        vLoad[1][i] = buffer.vector[i];
    }
    
    string = [name stringByAppendingString:@" 3"];
    stat = [self getReal:model forElement:element inArray:bc variableName:string buffer:&buffer listUtilities:listUtilities];
    for (i=0; i<nd; i++) {
        vLoad[2][i] = buffer.vector[i];
    }
    if (buffer.vector != NULL) {
        free_dvector(buffer.vector, 0, buffer.m-1);
    }

    g[0] = pNodes->x[kk] - pNodes->x[jj];
    g[1] = pNodes->y[kk] - pNodes->y[jj];
    g[2] = pNodes->z[kk] - pNodes->z[jj];
    vDSP_svesqD(g, 1, &sum, 3);
    sum = sqrt(sum);
    vDSP_vsdivD(g, 1, &sum, g, 1, 3);
    
    savedType = element->Type;
    if ([self getElementFamily:element] == 1) element->Type = *[elementDescription getElementType:202 inMesh:solution.mesh stabilization:NULL];
    
    *integral = 0.0;
    IP = GaussQuadrature(element, NULL, NULL);
    double buff[3];
    for (t=0; t<IP->n; t++) {
        stat = [numericIntegration setMetricDeterminantForElement:element 
                                                     elementNodes:nodes 
                                                           inMesh:solution.mesh 
                                             firstEvaluationPoint:IP->u[t]
                                            secondEvaluationPoint:IP->v[t]
                                             thirdEvaluationPoint:IP->w[t]];
        s = IP->s[t] * numericIntegration.metricDeterminant;
        
        stat = [numericIntegration setBasisForElement:element 
                                         elementNodes:nodes 
                                               inMesh:solution.mesh 
                                 firstEvaluationPoint:IP->u[t]
                                secondEvaluationPoint:IP->v[t]
                                 thirdEvaluationPoint:IP->w[t]
                                          withBubbles:NO 
                                          basisDegree:NULL];
        sum = cblas_ddot(nd, load.vector, 1, numericIntegration.basis, 1);
        l = sum;
        cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, nd, 1.0, *vLoad, np, numericIntegration.basis, 1, 0.0, vl, 1);
        vDSP_vmulD(vl, 1, g, 1, buff, 1, 3);
        vDSP_sveD(buff, 1, &sum, 3);
        *integral = *integral + s * ( l + sum);
    }
    element->Type = savedType;
    
    jj = parent->NodeIndexes[jj];
    // TODO: add support for parallel run
    kk = parent->NodeIndexes[kk];
    // TODO: add support for parallel run
    
    if (kk < jj) *integral = -*integral;
        
    free_dvector(nodes->x, 0, n-1);
    free_dvector(nodes->y, 0, n-1);
    free_dvector(nodes->z, 0, n-1);
    free(nodes);
    
    free_dvector(pNodes->x, 0, n-1);
    free_dvector(pNodes->y, 0, n-1);
    free_dvector(pNodes->z, 0, n-1);
    free(pNodes);
    
    if (load.vector != NULL) {
        free_dvector(load.vector, 0, load.m-1);
    }
    free_dmatrix(vLoad, 0, 2, 0, np-1);
    free_dvector(vl, 0, 2);
    
    [numericIntegration deallocation:solution.mesh];
}

/*************************************************************************************************************************************
    Given boundary condition, element and stiffness matrix and force vector, assemble boundary problem local
    stiffness matrix and force vector
 
    Arguments:
        FEMModel *model                 ->  class containing the model
        FEMSolution *solution           ->  class solution containing the global matrix
        NSArray *bc                     ->  boundary condition value list
        Element_t *element              ->  boundary element to get stiffness matrix to
        int nd                          ->  number of degrees of freedom in boundary element
        NSMutableString *name           ->  name of boundary condition
        double **stiff, double *force   ->  Output: boundary problem stiffness matrix and force vector
 
*************************************************************************************************************************************/
-(void)localBoundaryBDOFs:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution atBoundary:(NSArray * __nonnull)bc forElement:(Element_t * __nonnull)element withNumberOfNodes:(int)nd boundaryName:(NSMutableString * __nonnull)name resultMatrix:(double * __nonnull * __nonnull)stiff resultVector:(double * __nonnull)force {
    
    int i, n, p, q, t;
    double xip, yip, zip, s, load;
    FEMNumericIntegration *numericIntegration;
    FEMListUtilities *listUtil;
    Nodes_t *nodes;
    GaussIntegrationPoints *IP = NULL;
    BOOL stat;
    
    numericIntegration = [[FEMNumericIntegration alloc] init];
    if ([numericIntegration allocation:solution.mesh] == NO) fatal("FEMCore:localBoundaryBDOFs", "Allocation error in FEMNumericIntegration.");
    
    
    listUtil = [FEMListUtilities sharedListUtilities];
    
    n = max(solution.mesh.maxElementNodes, solution.mesh.maxElementDofs);
    
    // Get nodes of boundary elements parent and gauss points for boundary
    nodes = (Nodes_t*)malloc(sizeof(Nodes_t));
    initNodes(nodes);
    nodes->x = doublevec(0, n-1);
    nodes->y = doublevec(0, n-1);
    nodes->z = doublevec(0, n-1);
    
    [self getNodes:solution model:model inElement:element resultNodes:nodes numberOfNodes:&n mesh:nil];
    IP = GaussQuadrature(element, NULL, NULL);
    
    memset( force, 0.0, nd*sizeof(double) );
    memset( *stiff, 0.0, (nd*nd)*sizeof(double) );
    
    for (t=0; t<IP->n; t++) {
        stat = [numericIntegration setMetricDeterminantForElement:element 
                                                     elementNodes:nodes 
                                                           inMesh:solution.mesh 
                                             firstEvaluationPoint:IP->u[t]
                                            secondEvaluationPoint:IP->v[t]
                                             thirdEvaluationPoint:IP->w[t]];
        
        s = IP->s[t] * numericIntegration.metricDeterminant;
        
        stat = [numericIntegration setBasisForElement:element 
                                         elementNodes:nodes 
                                               inMesh:solution.mesh 
                                 firstEvaluationPoint:IP->u[t]
                                secondEvaluationPoint:IP->v[t]
                                 thirdEvaluationPoint:IP->w[t]
                                          withBubbles:NO 
                                          basisDegree:NULL];
        
        // Get value of boundary condition
        xip = 0.0;
        yip = 0.0;
        zip = 0.0;
        for (i=0; i<nd; i++) {
            xip = xip + numericIntegration.basis[i]*nodes->x[i];
            yip = yip + numericIntegration.basis[i]*nodes->y[i];
            zip = zip + numericIntegration.basis[i]*nodes->z[i];
        }
        //TODO: xip, yip, zip are used inside listGetConstReal when a user function is used.
        // But we don't support that yet so they are not used yet
        load = [listUtil listGetConstReal:model inArray:bc forVariable:name info:&stat minValue:NULL maxValue:NULL];
        
        // Build local stiffness matrix and force vector
        for (p=0; p<nd; p++) {
            for (q=0; q<nd; q++) {
                stiff[p][q] = stiff[p][q] + s * numericIntegration.basis[p] * numericIntegration.basis[q];
            }
            force[p] = force[p] + s * load * numericIntegration.basis[p];
        }
    }
    
    free_dvector(nodes->x, 0, n-1);
    free_dvector(nodes->y, 0, n-1);
    free_dvector(nodes->z, 0, n-1);
    free(nodes);
    
    [numericIntegration deallocation:solution.mesh];
}

-(void)solveWithLapackMatrix:(double * __nonnull)a andVector:(double * __nonnull)x size:(int)n leadingDimension:(int)lda {
    
    int nhrs, info;
    int *ipiv;
    char *trans;
    
    ipiv = intvec(0, n-1);
    
    if (n <= 0) return;
    dgetrf_(&n, &n, a, &lda, ipiv, &info);
    if (info < 0 || info > 0) {
        NSLog(@"FEMCore:solveWithLapackMatrix: error in lapack routine dgetrf. Error code: %d.\n", info);
        fatal("FEMCore:solveWithLapackMatrix", "Saino will abort the simulation now...");
    }
    
    trans = "N";
    nhrs = 1;
    dgetrs_(trans, &n, &nhrs, a, &lda, ipiv, x, &n, &info);
    if (info < 0 || info > 0) {
        NSLog(@"FEMCore:solveWithLapackMatrix: error in lapack routine dgetrs. Error code: %d.\n", info);
        fatal("FEMCore:solveWithLapackMatrix", "Saino will abort the simulation now...");
    }
    
    free_ivector(ipiv, 0, n-1);
}

-(void)solveLinearSystemWithMatrix:(double * __nonnull * __nonnull)a andVector:(double * __nonnull)x size:(int)n leadingDimension:(int)lda {
    
    int i, j;
    double *at, *b;
    FEMUtilities *util;
    
    util = [[FEMUtilities alloc] init];
    
    b = doublevec(0, n-1);
    
    switch (n) {
        case 1:
            x[0] = x[0] / a[0][0];
            break;
        case 2:
            for (i=0; i<n; i++) {
                b[i] = x[i];
            }
            [util solveLinearSystem2x2:a afterSolve:x rightHandSide:b];
            break;
        case 3:
            for (i=0; i<n; i++) {
                b[i] = x[i];
            }
            [util solveLinearSystem3x3:a afterSolve:x rightHandSide:b];
            break;
        default:
            at = doublevec(0, (lda*lda)-1);
            // Transfrom the matrix for LAPACK, column-major order
            for (i=0; i<lda; i++) {
                for (j=0; j<lda; j++) {
                    at[j+lda*i] = a[j][i];
                }
            }
            [self solveWithLapackMatrix:at andVector:x size:n leadingDimension:lda];
            
            // Back to the original matrix
            for (i=0; i<lda; i++) {
                for (j=0; j<lda; j++) {
                    a[j][i] = at[j+lda*i];
                }
            }
            free_dvector(at, 0, (lda*lda)-1);
            break;
    }
    
    free_dvector(b, 0, n-1);
}

/*****************************************************************************************************************
    Set nodel load for given dof
 
    FEMModel *model           -> class containing the model
    FEMSolution *solution     -> solution class containing the matrix and variable
    NSMutableString *name     -> name of the dof to be set
    int dof                   -> the order number of the dof
 
    The permutation (node reordoring info) is contained in the class solution and has been generated at the
    beginning of the simulation by the bandwidth optimization

*****************************************************************************************************************/
-(void)setNodalLoads:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution variableName:(NSString * __nonnull)name orderOfDofs:(int)dof {
    
    int i, j, n, bc, t, bf_id, noNodes=0, noDims=0;
    int *indexes, *inNodes;
    double minDist, dist, *diagScaling = NULL;
    listBuffer nodeIndexes = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    listBuffer coordNodes = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    BOOL *activePart, *activePartAll, *doneLoad = NULL, anyActive;
    Element_t *elements = NULL;
    Nodes_t *globalNodes = NULL;
    FEMListUtilities *listUtil;
    FEMBodyForce *bodyForceAtId;
    NSMutableString *loadName, *str;
    matrixArraysContainer *matContainers = NULL;
    variableArraysContainer *varContainers = NULL;
    BOOL stat, nodesFound;
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
    
    listUtil = [FEMListUtilities sharedListUtilities];
    
    loadName = [NSMutableString stringWithString:name];
    [loadName appendString:@" load"];
    
    NSLog(@"setNodalLoads: checking load for: %@.\n", name);
    
    n = max(model.numberOfBoundaryConditions, model.numberOfBodyForces);
    activePart = (BOOL*)malloc(sizeof(BOOL) * n );
    activePartAll = (BOOL*)malloc(sizeof(BOOL) * n );
    
    elements = solution.mesh.getElements;
    globalNodes = solution.mesh.getNodes;
    
    indexes = intvec(0, solution.mesh.maxElementDofs-1);
    
    // Go through the boundaries
    
    diagScaling = matContainers->DiagScaling;
    if (diagScaling == NULL) {
        diagScaling = doublevec(0, solution.matrix.numberOfRows-1);
        for (i=0; i<solution.matrix.numberOfRows; i++) {
            diagScaling[i] = 1.0;
        }
    }
    
    memset( activePart, NO, n*sizeof(BOOL) );
    memset( activePartAll, NO, n*sizeof(BOOL) );
    
    str = [NSMutableString stringWithString:loadName];
    [str appendString:@" dofs"];
    
    bc = 0;
    for (FEMBoundaryCondition *boundaryCondition in model.boundaryConditions) {
        if ([listUtil listCheckPresentVariable:@"target boundaries" inArray:boundaryCondition.valuesList] == NO) continue;
        activePart[bc] = [listUtil listCheckPresentVariable:loadName inArray:boundaryCondition.valuesList];
        activePartAll[bc] = [listUtil listCheckPresentVariable:str inArray:boundaryCondition.valuesList];
        bc++;
    }
    
    anyActive = NO;
    for (bc=0; bc<model.numberOfBoundaryConditions; bc++) {
        if (activePart[bc] == YES || activePartAll[bc] == YES) {
            anyActive = YES;
            break;
        }
    }
    
    if (anyActive == YES) {
        NSLog(@"setNodalLoads: settings nodals on boundaries: %@.\n", name);
        doneLoad = (BOOL*)malloc(sizeof(BOOL) *  (matContainers->sizeRHS/solution.variable.dofs) );
        memset( doneLoad, NO, (matContainers->sizeRHS/solution.variable.dofs)*sizeof(BOOL) );

        bc = 0;
        for (FEMBoundaryCondition *boundaryCondition in model.boundaryConditions) {
            if (activePart[bc] == NO && activePartAll[bc] == NO) {
                bc++;
                continue;
            };
        
            for (t=model.numberOfBulkElements; t<model.numberOfBulkElements+model.numberOfBoundaryElements; t++) {
                
                if (elements[t].BoundaryInfo->Constraint != boundaryCondition.tag) continue;
                
                if (activePart[bc] == YES) {
                    n = elements[t].Type.NumberOfNodes;
                    for (i=0; i<n; i++) {
                        indexes[i] = elements[t].NodeIndexes[i];
                    }
                } else {
                    n = [self sgetElementDofsSolution:solution model:model forElement:&elements[t] atIndexes:indexes];
                }
                
                [self FEMCore_setElementLoadsModel:model solution:solution element:&elements[t] values:boundaryCondition.valuesList name:loadName indexes:indexes doneLoad:doneLoad size:n dof:dof ndofs:solution.variable.dofs diagonalScaling:diagScaling];
            }
            bc++;
        }
    }
    
    // Go though the nodal load conditions for the body force list
        
    memset( activePart, NO, n*sizeof(BOOL) );
    memset( activePartAll, NO, n*sizeof(BOOL) );
    
    bf_id = 0;
    for (FEMBodyForce *bodyForce in model.bodyForces) {
        activePart[bf_id] = [listUtil listCheckPresentVariable:loadName inArray:bodyForce.valuesList];
        activePartAll[bf_id] = [listUtil listCheckPresentVariable:str inArray:bodyForce.valuesList];
        bf_id++;
    }
    
    anyActive = NO;
    for (bf_id=0; bf_id<model.numberOfBodyForces; bf_id++) {
        if (activePart[bf_id] == YES || activePartAll[bf_id] == YES) {
            anyActive = YES;
            break;
        }
    }
    
    if (anyActive == YES) {
        NSLog(@"setNodalLoads: settings nodals on bulk: %@.\n", name);
        if (doneLoad == NULL) doneLoad = (BOOL*)malloc(sizeof(BOOL) *  (matContainers->sizeRHS/solution.variable.dofs) );
        memset( doneLoad, NO, (matContainers->sizeRHS/solution.variable.dofs)*sizeof(BOOL) );
        for (t=0; t<model.numberOfBulkElements; t++) {
            
            if ((model.bodies)[elements[t].BodyID-1][@"body force"] == nil) continue;
            bf_id = [(model.bodies)[elements[t].BodyID-1][@"body force"] intValue];
            
            if (activePart[bf_id] == NO && activePartAll[bf_id] == NO) continue;
            
            
            if (activePart[bf_id] == YES) {
                n = elements[t].Type.NumberOfNodes;
                memcpy(indexes, elements[t].NodeIndexes, n*sizeof(int));
            } else {
                n = [self sgetElementDofsSolution:solution model:model forElement:&elements[t] atIndexes:indexes];
            }
            
            bodyForceAtId = (model.bodyForces)[bf_id];
            [self FEMCore_setElementLoadsModel:model solution:solution element:&elements[t] values:bodyForceAtId.valuesList name:loadName indexes:indexes doneLoad:doneLoad size:n dof:dof ndofs:solution.variable.dofs diagonalScaling:diagScaling];
        }
    }
    
    free(activePart);
    free(activePartAll);
    if (doneLoad != NULL) free(doneLoad);
    
    // Go through the point loads which are created on the fly
    for (FEMBoundaryCondition *boundaryCondition in model.boundaryConditions) {
        if ([listUtil listCheckPresentVariable:loadName inArray:boundaryCondition.valuesList] == NO) continue;
        nodesFound = [listUtil listCheckPresentVariable:@"target nodes" inArray:boundaryCondition.valuesList];
        
        // At the first calling, the list of coordinates is transformed to list of nodes
        if (nodesFound == NO) {
            
            stat = [listUtil listGetConstRealArray:model inArray:boundaryCondition.valuesList forVariable:@"target coordinates" buffer:&coordNodes];
            
            if (stat == YES) {
                
                noNodes = coordNodes.m;
                noDims = coordNodes.n;
                
                if (noNodes > 0) {
                    inNodes = intvec(0, noNodes-1);
                    for (j=0; j<noNodes; j++) {
                        minDist = HUGE_VAL;
                        for (i=0; i<model.numberOfNodes; i++) {
                            if (varContainers->Perm[i] < 0) continue;
                            
                            dist = pow((globalNodes->x[i]-coordNodes.matrix[j][0]), 2.0);
                            if (noDims >= 2) dist = dist + pow((globalNodes->y[i]-coordNodes.matrix[j][1]), 2.0);
                            if (noDims == 3) dist = dist + pow((globalNodes->z[i]-coordNodes.matrix[j][2]), 2.0);
                            
                            if (dist<minDist) {
                                minDist = dist;
                                inNodes[j] = i;
                            }
                        }
                    }
                    // Add the found nodes to the list values
                    [listUtil addIntegerArrayInClassList:boundaryCondition theVariable:@"target nodes" withValues:inNodes size:noNodes orUsingBlock:nil];
                    free_ivector(inNodes, 0, noNodes-1);
                    nodesFound = YES;
                }
            }
        }
        
        if (coordNodes.matrix != NULL) {
            free_dmatrix(coordNodes.matrix, 0, coordNodes.m-1, 0, coordNodes.n-1);
        }
        
        if (nodesFound == YES) {
            NSLog(@"setNodalLoads: settings nodal nodals on target nodes: %@.\n", name);
            [listUtil listGetIntegerArray:model inArray:boundaryCondition.valuesList forVariable:@"target nodes" buffer:&nodeIndexes];
            n = nodeIndexes.m;
            [self FEMCore_setPointLoadsModel:model solution:solution element:elements values:boundaryCondition.valuesList name:loadName indexes:nodeIndexes.ivector size:n dof:dof ndofs:solution.variable.dofs diagonalScaling:diagScaling];
            if (nodeIndexes.ivector != NULL) {
                free_ivector(nodeIndexes.ivector, 0, nodeIndexes.m-1);
            }
        }
    }
    
    free_ivector(indexes, 0, solution.mesh.maxElementDofs-1);
    if (diagScaling != matContainers->DiagScaling) free_dvector(diagScaling, 0, solution.matrix.numberOfRows-1);
}

/*****************************************************************************************************************
    Set dirichlet boundary condition for given dof. The conditions are set based on the given name and applied 
    directly to the matrix structure so that a row is zeroed execpt for the diagonal which is set to one. Then 
    the RHS balue determines the value of the field variable in the solution of the linear system.
 
    FEMModel *model           -> class containing the model
    FEMSolution *solution     -> solution class containing the matrix and variable
    NSMutableString *name     -> name of the dof to be set
    int dof                   -> the order number of the dof
    int *offset               -> Optional: If the matrix and pernutation vectors are not in sync the offset
                                 may be used as a remedy. Needed in fully coupled systems.
    BOLL *offDiaginalMatrix   -> Optional: For block systems only the diagonal matrix should be given non-zero
                                 entries for matrix and RHS. for off-diagonal matrices just set the row to zero.

    The permutation (node reordoring info) is contained in the class solution and has been generated at the
    beginning of the simulation by the bandwidth optimization
 
*****************************************************************************************************************/
-(void)setDirichletBoundaries:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution variableName:(NSMutableString * __nonnull)name orderOfDofs:(int)dof permutationOffset:(int * __nullable)offset offDiaginalMatrix:(BOOL * __nullable)offDiaginalMatrix {
    
    int i, j, k, n, bc, t, bf_id, bndry_start, bndry_end, noNodes, noDims, permOffset, numberOfNodesFound;
    int *indexes;
    double dist, eps, *diagScaling = NULL;
    listBuffer nodeIndexes = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    listBuffer coordNodes = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    BOOL *activePart, *activePartAll, *activeCond, *donePeriodic, anyActive, offDiaginal, passive;
    Element_t *elements = NULL;
    Nodes_t *globalNodes = NULL;
    FEMListUtilities *listUtil;
    FEMBodyForce *bodyForceAtId;
    NSString *str1, *str2;
    matrixArraysContainer *matContainers = NULL;
    variableArraysContainer *varContainers = NULL;
    solutionArraysContainer *solContainers = NULL;
    BOOL stat, nodesFound, orderByBCNumbering, conditional;
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
    solContainers = solution.getContainers;
    
    listUtil = [FEMListUtilities sharedListUtilities];
    
    diagScaling = matContainers->DiagScaling;
    if (diagScaling == NULL) {
        diagScaling = doublevec(0, solution.matrix.numberOfRows-1);
        for (i=0; i<solution.matrix.numberOfRows; i++) {
            diagScaling[i] = 1.0;
        }
    }
    
    // These logical vectors are used to minimize extra effort in setting up different BCs
    int activeSize = max(model.numberOfBoundaryConditions, model.numberOfBodyForces);
    activePart = (BOOL*)malloc(sizeof(BOOL) * activeSize );
    activePartAll = (BOOL*)malloc(sizeof(BOOL) * activeSize );
    activeCond = (BOOL*)malloc(sizeof(BOOL) * activeSize );
    
    NSString *condName = [name stringByAppendingString:@" condition"];
    NSString *passName = [name stringByAppendingString:@" passive"];
    NSString *passCondName = [name stringByAppendingString:@" condition passive"];
    
    permOffset = 0;
    if (offset != NULL) permOffset = *offset;
    
    offDiaginal = NO;
    if (offDiaginalMatrix != NULL) offDiaginal = *offDiaginalMatrix;
    
    elements = solution.mesh.getElements;
    globalNodes = solution.mesh.getNodes;

    indexes = intvec(0, solution.mesh.maxElementDofs-1);
    
    // Go through the perdiodic BCs and set the linear dependence
    
    memset( activePart, NO, activeSize*sizeof(BOOL) );
    
    str1 = [@"periodic bc " stringByAppendingString:name];
    str2 = [@"anti periodic bc " stringByAppendingString:name];

    bc = 0;
    for (FEMBoundaryCondition *boundaryCondition in model.boundaryConditions) {
        if ([listUtil listGetLogical:model inArray:boundaryCondition.valuesList forVariable:str1 info:&stat] == YES) activePart[bc] = YES;
        if ([listUtil listGetLogical:model inArray:boundaryCondition.valuesList forVariable:str2 info:&stat] == YES) activePart[bc] = YES;
        bc++;
    }
    
    anyActive = NO;
    for (bc=0; bc<model.numberOfBoundaryConditions; bc++) {
        if (activePart[bc] == YES) {
            anyActive = YES;
            break;
        }
    }
    
    if (anyActive == YES) {
        
        if (permOffset > 0) fatal("setDirichletBoundaries", "Periodicity not supported with offset.");
        
        donePeriodic = (BOOL*)malloc(sizeof(BOOL) * solution.mesh.numberOfNodes );
        memset( donePeriodic, NO, solution.mesh.numberOfNodes*sizeof(BOOL) );
        
        for (bc=0; bc<model.numberOfBoundaryConditions; bc++) {
            if (activePart[bc] == YES) [self FEMCore_setPeriodicBoundariesPass1Model:model solution:solution name:name dof:dof this:bc done:donePeriodic diaginalScaling:diagScaling];
            
        }
        
        memset( donePeriodic, NO, solution.mesh.numberOfNodes*sizeof(BOOL) );
        for (bc=0; bc<model.numberOfBoundaryConditions; bc++) {
            if (activePart[bc] == YES) [self FEMCore_setPeriodicBoundariesPass2Model:model solution:solution name:name dof:dof this:bc done:donePeriodic diaginalScaling:diagScaling];
        }
        
        free(donePeriodic);
    }
    
    // Go through the normal Dirichlet BCs applied on the boundaries
    
    memset( activePart, NO, activeSize*sizeof(BOOL) );
    memset( activePartAll, NO, activeSize*sizeof(BOOL) );
    memset( activeCond, NO, activeSize*sizeof(BOOL) );

    str1 = [name stringByAppendingString:@" dofs"];
    
    bc = 0;
    for (FEMBoundaryCondition *boundaryCondition in model.boundaryConditions) {
        activePartAll[bc] = [listUtil listCheckPresentVariable:str1 inArray:boundaryCondition.valuesList];
        activePart[bc] = [listUtil listCheckPresentVariable:name inArray:boundaryCondition.valuesList];
        activeCond[bc] = [listUtil listCheckPresentVariable:condName inArray:boundaryCondition.valuesList];
        bc++;
    }
    
    orderByBCNumbering = [listUtil listGetLogical:model inArray:model.simulation.valuesList forVariable:@"set dirichlet boundaries by boundary numbering" info:&stat];
    
    bndry_start = model.numberOfBulkElements;
    bndry_end = bndry_start + model.numberOfBoundaryElements;
    
    // Check and set some flags for nodes belonging to n-t boundaries getting set by other bcs
    if (self.normalTangentialNumberOfNodes > 0) {
        if (orderByBCNumbering == YES) {
            bc = 0;
            for (FEMBoundaryCondition *boundaryCondition in model.boundaryConditions) {
                
                if (activePart[bc] == NO && activePartAll[bc] == NO) {
                    bc++;
                    continue;
                };
                conditional = activeCond[bc];
                
                for (t=bndry_start; t<bndry_end; t++) {
                    if (elements[t].BoundaryInfo->Constraint != boundaryCondition.tag) continue;
                    if (activePart[bc] == YES) {
                        n = elements[t].Type.NumberOfNodes;
                        memcpy(indexes, elements[t].NodeIndexes, n*sizeof(int));
                    } else {
                        n = [self sgetElementDofsSolution:solution model:model forElement:&elements[t] atIndexes:indexes];
                    }
                    [self checkNormalTangential:model inSolution:solution forElementNumber:t numberofNodes:n atIndexes:indexes atBoundary:bc variableName:name orderOfDofs:dof activeCondition:conditional conditionName:condName permutationOffset:permOffset];
                }
                bc++;
            }
        } else {
             for (t=bndry_start; t<bndry_end; t++) {
                 bc = 0;
                 for (FEMBoundaryCondition *boundaryCondition in model.boundaryConditions) {
                     if (activePart[bc] == NO && activePartAll[bc] == NO) {
                         bc++;
                         continue;
                     };
                     conditional = activeCond[bc];
                     
                     if (elements[t].BoundaryInfo->Constraint != boundaryCondition.tag) {
                         bc++;
                         continue;
                     };
                     
                     if (activePart[bc] == YES) {
                         n = elements[t].Type.NumberOfNodes;
                         memcpy(indexes, elements[t].NodeIndexes, n*sizeof(int));
                     } else {
                         n = [self sgetElementDofsSolution:solution model:model forElement:&elements[t] atIndexes:indexes];
                     }
                     [self checkNormalTangential:model inSolution:solution forElementNumber:t numberofNodes:n atIndexes:indexes atBoundary:bc variableName:name orderOfDofs:dof activeCondition:conditional conditionName:condName permutationOffset:permOffset];
                     bc++;
                 }
             }
        }
        
        if (dof < 0) {
            for (t=bndry_start; t<bndry_end; t++) {
                n = elements[t].Type.NumberOfNodes;
                for (j=0; j<n; j++) {
                    k = self.boundaryReorder[elements[t].NodeIndexes[j]];
                    if (k >= 0) {
                        for (i=0; i<3; i++) {
                            self.ntElement[k][i] = 0;
                            self.ntZeroingDone[k][i] = false;
                        }
                    }
                }
            }
        }
    }
    
    // Set the Dirichlet BCs from active boundary elements, if any...
    anyActive = NO;
    for (bc=0; bc<model.numberOfBoundaryConditions; bc++) {
        if (activePart[bc] == YES || activePartAll[bc] == YES) {
            anyActive = YES;
            break;
        }
    }
    
    if (anyActive == YES) {
        if (orderByBCNumbering == YES) {
            bc = 0;
            for (FEMBoundaryCondition *boundaryCondition in model.boundaryConditions) {
                
                if (activePart[bc] == NO && activePartAll[bc] == NO) {
                    bc++;
                    continue;
                };
                conditional = activeCond[bc];
                
                for (t=bndry_start; t<bndry_end; t++) {
                    if (elements[t].BoundaryInfo->Constraint != boundaryCondition.tag) continue;
                    if (activePart[bc] == YES) {
                        n = elements[t].Type.NumberOfNodes;
                        memcpy(indexes, elements[t].NodeIndexes, n*sizeof(int));
                    } else {
                        n = [self sgetElementDofsSolution:solution model:model forElement:&elements[t] atIndexes:indexes];
                    }
                    [self FEMCore_setElementValues:model inSolution:solution forElementNumber:t numberOfNodes:n atIndexes:indexes forValues:boundaryCondition.valuesList variableName:name orderOfDofs:dof activeCondition:conditional conditionName:condName permutationOffset:permOffset diaginalScaling:diagScaling offDiaginal:offDiaginal];
                }
                bc++;
            }
        } else {
            for (t=bndry_start; t<bndry_end; t++) {
                bc = 0;
                for (FEMBoundaryCondition *boundaryCondition in model.boundaryConditions) {
                    if (activePart[bc] == NO && activePartAll[bc] == NO) {
                        bc++;
                        continue;
                    };
                    conditional = activeCond[bc];
                    
                    if (elements[t].BoundaryInfo->Constraint != boundaryCondition.tag) {
                        bc++;
                        continue;
                    };
                    
                    if (activePart[bc] == YES) {
                        n = elements[t].Type.NumberOfNodes;
                        memcpy(indexes, elements[t].NodeIndexes, n*sizeof(int));
                    } else {
                        n = [self sgetElementDofsSolution:solution model:model forElement:&elements[t] atIndexes:indexes];
                    }
                    [self FEMCore_setElementValues:model inSolution:solution forElementNumber:t numberOfNodes:n atIndexes:indexes forValues:boundaryCondition.valuesList variableName:name orderOfDofs:dof activeCondition:conditional conditionName:condName permutationOffset:permOffset diaginalScaling:diagScaling offDiaginal:offDiaginal];
                    bc++;
                }
            }
        }
    }
    
    // Go through the Dirichlet conditions in the body force lists
    
    memset( activePart, NO, activeSize*sizeof(BOOL) );
    memset( activeCond, NO, activeSize*sizeof(BOOL) );
    memset( activePartAll, NO, activeSize*sizeof(BOOL) );
    passive = NO;
    
    bf_id = 0;
    for (FEMBodyForce *bodyForce in model.bodyForces) {
        activePartAll[bf_id] = [listUtil listCheckPresentVariable:str1 inArray:bodyForce.valuesList];
        activeCond[bf_id] = [listUtil listCheckPresentVariable:condName inArray:bodyForce.valuesList];
        activePart[bf_id] = [listUtil listCheckPresentVariable:name inArray:bodyForce.valuesList ];
        passive = (passive == YES || [listUtil listCheckPresentVariable:passName inArray:bodyForce.valuesList] == YES) ? YES : NO;
        bf_id++;
    }
    
    anyActive = NO;
    for (bf_id=0; bf_id<model.numberOfBodyForces; bf_id++) {
        if (activePart[bf_id] == YES || activePartAll[bf_id] == YES) {
            anyActive = YES;
            break;
        }
    }

    if (anyActive == YES) {
        int passMerm[solution.mesh.numberOfNodes];
        int nodeIndexes[1];
        for (i=1; i<=solution.mesh.numberOfPassiveBCs; i++) {
            j = solution.mesh.numberOfBulkElements+solution.mesh.numberOfBoundaryElements-i;
            for (k=0; k<elements[j].Type.NumberOfNodes; k++) {
                passMerm[elements[j].NodeIndexes[k]] = 1;
            }
        }
        
        for (t=0; t<model.numberOfBulkElements; t++) {
            
            if ((model.bodies)[elements[t].BodyID-1][@"body force"] == nil) continue;
            bf_id = [(model.bodies)[elements[t].BodyID-1][@"body force"] intValue];
            
            if (activePart[bf_id] == NO && activePartAll[bf_id] == NO) continue;
            conditional = activeCond[bf_id];
            
            if (activePart[bf_id] == YES) {
                n = elements[t].Type.NumberOfNodes;
                memcpy(indexes, elements[t].NodeIndexes, n*sizeof(int));
            } else {
                n = [self sgetElementDofsSolution:solution model:model forElement:&elements[t] atIndexes:indexes];
            }
            bodyForceAtId = (model.bodyForces)[bf_id];
            
            if ([listUtil listGetLogical:model inArray:bodyForceAtId.valuesList forVariable:passCondName info:&stat] == YES) {
                if ([self FEMCore_checkPassiveElement:&elements[t] model:model solution:solution] == NO) continue;
                for (j=0; j<n; j++) {
                    nodeIndexes[0] = indexes[j];
                    if (passMerm[nodeIndexes[0]] == 0) [self FEMCore_setPointValues:model inSolution:solution numberofNodes:1 atIndexes:nodeIndexes forValues:bodyForceAtId.valuesList variableName:name orderOfDofs:dof activeCondition:conditional conditionName:condName permutationOffset:permOffset diaginalScaling:diagScaling offDiaginal:offDiaginal];
                }
            } else {
                [self FEMCore_setElementValues:model inSolution:solution forElementNumber:t numberOfNodes:n atIndexes:indexes forValues:bodyForceAtId.valuesList variableName:name orderOfDofs:dof activeCondition:conditional conditionName:condName permutationOffset:permOffset diaginalScaling:diagScaling offDiaginal:offDiaginal];
            }
        }
    }
    free(activePart);
    free(activePartAll);
    free(activeCond);
    free_ivector(indexes, 0, solution.mesh.maxElementDofs-1);
    
    // Go through the pointwise Dirichlet BCs that are created on the fly.
    // Note that it is best that the coordinates are transformed to nodes using
    // the right variable. Otherwise, it could point to nodes that are not active

    for (FEMBoundaryCondition *boundaryCondition in model.boundaryConditions) {
        
        if ([listUtil listCheckPresentVariable:name inArray:boundaryCondition.valuesList] == NO) continue;
        nodesFound = [listUtil listCheckPresentVariable:@"target nodes" inArray:boundaryCondition.valuesList];
        
        // The coodinates are only requested for a body that has no list of nodes.
        // At the first calling, the list of coordinates is transformed to a list of nodes
        if (nodesFound == NO) {
            
            stat = [listUtil listGetConstRealArray:model inArray:boundaryCondition.valuesList forVariable:@"target coordinates" buffer:&coordNodes];
            
            if (stat == YES) {
                
                eps = [listUtil listGetConstReal:model inArray:boundaryCondition.valuesList forVariable:@"target coordinates eps" info:&stat minValue:NULL maxValue:NULL];
                if (stat == NO) {
                    eps = HUGE_VAL;
                } else {
                    eps = pow(eps, 2.0);
                }
                
                noNodes = coordNodes.m;
                noDims = coordNodes.n;
                
                if (noNodes > 0) {
                    int indNodes[noNodes];
                    int minDist[noNodes];
                    memset( indNodes, -1, sizeof(indNodes) );
                    for (j=0; j<noNodes; j++) {
                        minDist[j] = HUGE_VAL;
                    }
                    for (j=0; j<noNodes; j++) {
                        for (i=0; i<model.numberOfNodes; i++) {
                            if (varContainers->Perm[i] < 0) continue;
                            
                            dist = pow((globalNodes->x[i]-coordNodes.matrix[j][0]), 2.0);
                            if (noDims >= 2) dist = dist + pow((globalNodes->y[i]-coordNodes.matrix[j][1]), 2.0);
                            if (noDims == 3) dist = dist + pow((globalNodes->z[i]-coordNodes.matrix[j][2]), 2.0);
                            dist = sqrt(dist);
                            
                            if (dist < minDist[j] && dist <= eps) {
                                minDist[j] = dist;
                                indNodes[j] = i;
                            }
                        }
                    }
                    
                    // In parallel case, eliminate all exept the nearest node.
                    // This relies on the fact that for each node partition, the distance
                    // to nearest node is computed accurately
                    FEMParallelMPI *parallelUtil = [[FEMParallelMPI alloc] init];
                    int oper = 1;
                    double globalMinDist;
                    for (j=0; j<noNodes; j++) {
                        globalMinDist = [parallelUtil parallelReductionOfValue:minDist[j] operArg:&oper];
                        if (fabs(globalMinDist - minDist[j]) > DBL_MIN) indNodes[j] = -1;
                    }
                    
                    numberOfNodesFound = 0;
                    for (j=0; j<noNodes; j++) {
                        if (indNodes[j] >= 0) {
                            indNodes[numberOfNodesFound] = indNodes[j];
                            numberOfNodesFound++;
                        }
                    }
                    
                    // In the first time add the found nodes to the list
                    if (numberOfNodesFound > 0) {
                        [listUtil addIntegerArrayInClassList:boundaryCondition theVariable:@"target nodes" withValues:indNodes size:numberOfNodesFound orUsingBlock:nil];
                        nodesFound = YES;
                    } else {
                        // If no nodes found, add still an empty list and make sure the negative value is not
                        // treated later one. Otherwise this search would be retreated each time
                        [listUtil addIntegerArrayInClassList:boundaryCondition theVariable:@"target nodes" withValues:indNodes size:1 orUsingBlock:nil];
                    }
                }
            }
            
            if (coordNodes.matrix != NULL) {
                free_dmatrix(coordNodes.matrix, 0, coordNodes.m-1, 0, coordNodes.n-1);
            }
        }
        
        // If the target coordinates has already been assigned to an empty list,
        // cycle over it by testing the first node
        if (nodesFound == YES) {
            [listUtil listGetIntegerArray:model inArray:boundaryCondition.valuesList forVariable:@"target nodes" buffer:&nodeIndexes];
            if (nodeIndexes.ivector[0] == -1) nodesFound = NO;
        }
        
        if (nodesFound == YES) {
            conditional = [listUtil listCheckPresentVariable:condName inArray:boundaryCondition.valuesList];
            n = nodeIndexes.m;
            [self FEMCore_setPointValues:model inSolution:solution numberofNodes:n atIndexes:nodeIndexes.ivector forValues:boundaryCondition.valuesList variableName:name orderOfDofs:dof activeCondition:conditional conditionName:condName permutationOffset:permOffset diaginalScaling:diagScaling offDiaginal:offDiaginal];
        }
        if (nodeIndexes.ivector != NULL) {
            free_ivector(nodeIndexes.ivector, 0, nodeIndexes.m-1);
            nodeIndexes.ivector = NULL;
        }
    }
    
    // Go through soft upper and lower  limits
    if (dof >= 0 && [(solution.solutionInfo)[@"apply limiter"] boolValue] == YES) {
        NSString *condName;
        bool *limitActive;
        for (int upper=0; upper<=1; upper++) {
            limitActive = NULL;
            // The limiters have been implemented only component-wise
            if (upper == 0) {
                if (varContainers->lowerLimitActive != NULL) limitActive = varContainers->lowerLimitActive;
            } else {
                if (varContainers->upperLimitActive != NULL) limitActive = varContainers->upperLimitActive;
            }
            
            if (limitActive == NULL) continue;
            if (upper == 0) {
                condName = [name stringByAppendingString:@" lower limit"];
            } else {
                condName = [name stringByAppendingString:@" upper limit"];
            }
            
            // Check and set sone flags for nodes belonging to n-t boundaries
            // getting set by other bcs
            for (t=0; t<model.numberOfBulkElements+model.numberOfBoundaryConditions; t++) {
                n = elements[t].Type.NumberOfNodes;
                if (t>=model.numberOfBulkElements) {
                    for (FEMBoundaryCondition *boundaryCondition in model.boundaries) {
                        if (elements[t].BoundaryInfo->Constraint != boundaryCondition.tag) continue;
                        [self FEMCore_setLimiterValues:model inSolution:solution numberofNodes:n atIndexes:elements[t].NodeIndexes forValues:boundaryCondition.valuesList orderOfDofs:dof limitActive:limitActive conditionName:condName permutationOffset:permOffset diaginalScaling:diagScaling offDiaginal:offDiaginal];
                    }
                } else {
                    if ((model.bodies)[elements[t].BodyID-1][@"body force"] == nil) continue;
                    bf_id = [(model.bodies)[elements[t].BodyID-1][@"body force"] intValue];
                    bodyForceAtId = (model.bodyForces)[bf_id];
                    [self FEMCore_setLimiterValues:model inSolution:solution numberofNodes:n atIndexes:elements[t].NodeIndexes forValues:bodyForceAtId.valuesList orderOfDofs:dof limitActive:limitActive conditionName:condName permutationOffset:permOffset diaginalScaling:diagScaling offDiaginal:offDiaginal];
                }
            }
        }
    }
    
    // Take care of the matrix entries of passive elements
    if (passive == YES) {
        int m;
        double s;
        for (i=0; i<solution.numberOfActiveElements; i++) {
            if ([self FEMCore_checkPassiveElement:&elements[solContainers->activeElements[i]] model:model solution:solution] == YES) {
                n = [self sgetElementDofsSolution:solution model:model forElement:&elements[solContainers->activeElements[i]] atIndexes:indexes];
                for (j=0; j<n; j++) {
                    k = indexes[j];
                    if (k < 0) continue;
                    
                    k = varContainers->Perm[k];
                    if (k < 0) continue;
                    s = 0.0;
                    for (int l=0; l<solution.variable.dofs; l++) {
                        m = solution.variable.dofs * k + l;
                        s = s + fabs(matContainers->Values[matContainers->Diag[m]]);
                    }
                    if (s > DBL_EPSILON) continue;
                    
                    for (int l=0; l<solution.variable.dofs; l++) {
                        m = solution.variable.dofs * k + l;
                        [self zeroTheNumberOfRows:m inSolutionMatrix:solution];
                        matContainers->Values[matContainers->Diag[m]] = 1.0;
                        matContainers->RHS[m] = varContainers->Values[m]/diagScaling[m];
                    }
                }
            }
        }
    }
    
    if (diagScaling != matContainers->DiagScaling) free_dvector(diagScaling, 0, solution.matrix.numberOfRows-1);
}

/****************************************************************************************************************
    Scale system Ax = b as:
    (DAD) = Db, where D = 1/sqrt(Diag(A)) and y = D^-1 x
*****************************************************************************************************************/
-(void)scaleLinearSystem:(FEMSolution * __nonnull)solution matrix:(FEMMatrix * __nonnull)matrix rhs:(double * __nullable)b result:(double * __nullable)x diagScaling:(double * __nullable)diagScaling applyScaling:(BOOL * __nullable)applyScaling rhsScaling:(BOOL * __nullable)rhsScaling {
    
    int i, j, n;
    double *diag, bnorm, sum;
    double complex diagC;
    BOOL complexMatrix, doRHS;
    matrixArraysContainer *matContainers = NULL;
    FEMParallelMPI *parallelUtil;
    
    n = matrix.numberOfRows;
    matContainers = matrix.getContainers;
    
    if (diagScaling != NULL) {
        diag = diagScaling;
    } else {
        if (matContainers->DiagScaling == NULL) {
            matContainers->DiagScaling = doublevec(0, n-1);
            matContainers->sizeDiagScaling = n;
        }
        diag = matContainers->DiagScaling;
        
        complexMatrix = solution.matrix.complexMatrix;
        
        if (complexMatrix == YES) {
            for (i=0; i<n; i+=2) {
                j = matContainers->Diag[i];
                diag[i] = matContainers->Values[j];
                diag[i+1] = matContainers->Values[j+1];
            }
        } else {
            for (i=0; i<n; i++) {
                diag[i] = matContainers->Values[matContainers->Diag[i]];
            }
        }
        
        // TODO: Add support for parallel run
        
        if (complexMatrix == YES) {
            for (i=0; i<n; i+=2) {
                diagC = diag[i] + (-diag[i+1] * I);
                if (cabs(diagC) != 0.0) {
                    diag[i] = 1.0 / sqrt(cabs(diagC));
                    diag[i+1] = 1.0 / sqrt(cabs(diagC));
                } else {
                    diag[i] = 1.0;
                    diag[i+1] = 1.0;
                }
            }
        } else {
            for (i=0; i<n; i++) {
                if (fabs(diag[i]) != 0.0) {
                    diag[i] = 1.0 / sqrt(fabs(diag[i]));
                } else {
                    diag[i] = 1.0;
                }
            }
        }
    }
    
    // Optionally we may just create the diag and leave the scaling undone
    if (applyScaling != NULL) {
        if (applyScaling == NO) return;
    }
    
    for (i=0; i<n; i++) {
        for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
            matContainers->Values[j] = matContainers->Values[j] * ( diag[i] * diag[matContainers->Cols[j]]);
        }
    }
    
    if (matContainers->MassValues != NULL) {
        if (matContainers->sizeValues == matContainers->sizeMassValues) {
            for (i=0; i<n; i++) {
                for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
                    matContainers->MassValues[j] = matContainers->MassValues[j] * (diag[i] * diag[matContainers->Cols[j]]);
                }
            }
        }
    }
    
    if (matContainers->DampValues != NULL) {
        if (matContainers->sizeValues == matContainers->sizeDampValues) {
            for (i=0; i<n; i++) {
                for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
                    matContainers->DampValues[j] = matContainers->DampValues[j] * (diag[i] * diag[matContainers->Cols[j]]);
                }
            }
        }
    }
    
    // Scale r.h.s and initial guess
    if (b != NULL) {
        for (i=0; i<n; i++) {
            b[i] = b[i] * diag[i];
        }
        doRHS = YES;
        if (rhsScaling != NULL) doRHS = *rhsScaling;
        if (doRHS == YES) {
            parallelUtil = [[FEMParallelMPI alloc] init];
            vDSP_svesqD(b, 1, &sum, n);
            bnorm = [parallelUtil parallelReductionOfValue:sqrt(sum) operArg:NULL];
        } else {
            bnorm = 1.0;
        }
        matrix.rhsScaling = bnorm;
        
        for (i=0; i<n; i++) {
            diag[i] = diag[i] * bnorm;
            b[i] = b[i] / bnorm;
        }
        if (x != NULL) {
            for (i=0; i<n; i++) {
                x[i] = x[i] / diag[i];
            }
        }
    }
}

/************************************************************************************
    Scale the system back to original.
    diagscaling is optional and if given, its size sizeOFDiagScaling should also
    be given.
************************************************************************************/
-(void)backScaleLinearSystem:(FEMSolution * __nonnull)solution matrix:(FEMMatrix * __nonnull)matrix rhs:(double * __nullable)b result:(double * __nullable)x diagScaling:(double * __nullable)diagScaling sizeOFDiagScaling:(int * __nullable)sizeOfDiagScaling {
    
    int i, j, k, n;
    double *diag, bnorm;
    matrixArraysContainer *matContainers = NULL;
    variableArraysContainer *varContainers = NULL;
    
    n = matrix.numberOfRows;
    matContainers = matrix.getContainers;
    
    if (diagScaling != NULL) {
        diag = diagScaling;
    } else {
        diag = matContainers->DiagScaling;
    }
    
    if (diag == NULL) {
        fatal("FEMCore:backScaleLinearSystem", "Diag is a null pointer.");
    }
    if (diagScaling != NULL) {
        if (*sizeOfDiagScaling != n) fatal("FEMCore:backScaleLinearSystem", "Diag of wrong size.");
    } else {
        if (matContainers->sizeDiagScaling != n) fatal("FEMCore:backScaleLinearSystem", "Diag of wrong size.");
    }
    
    if (b != NULL) {
        
        // Solve x: INV(D)x = y, scale b back to original
        if (x != NULL) {
            for (i=0; i<n; i++) {
                x[i] = x[i] * diag[i];
            }
        }
        bnorm = matrix.rhsScaling;
        for (i=0; i<n; i++) {
            diag[i] = diag[i] / bnorm;
        }
        for (i=0; i<n; i++) {
            b[i] = b[i] / diag[i] * bnorm;
        }
    }
    
    varContainers = solution.variable.getContainers;
    for (i=0; i<solution.nOfEigenValues; i++) {
        
        // Solve x: INV(D)x = y
        if (solution.matrix.complexMatrix == YES) {
            k = 0;
            for (j=0; j<n/2; j++) {
                varContainers->EigenVectors[i][j] = varContainers->EigenVectors[i][j] * diag[k];
                k+=2;
            }
        } else {
            for (j=0; j<n; j++) {
                varContainers->EigenVectors[i][j] = varContainers->EigenVectors[i][j] * diag[j];
            }
            
        }
    }
    
    for (i=0; i<n; i++) {
        for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
            matContainers->Values[j] = matContainers->Values[j] / (diag[i] * diag[matContainers->Cols[j]]);
        }
    }
    
    if (matContainers->MassValues != NULL) {
        if (matContainers->sizeValues == matContainers->sizeMassValues) {
            for (i=0; i<n; i++) {
                for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
                    matContainers->MassValues[j] = matContainers->MassValues[j] / (diag[i] * diag[matContainers->Cols[j]]);
                }
            }
        }
    }
    
    if (matContainers->DampValues != NULL) {
        if (matContainers->sizeValues == matContainers->sizeDampValues) {
            for (i=0; i<n; i++) {
                for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
                    matContainers->DampValues[j] = matContainers->DampValues[j] / (diag[i] * diag[matContainers->Cols[j]]);
                }
            }
        }
    }
    
    free_dvector(matContainers->DiagScaling, 0, matContainers->sizeDiagScaling-1);
    matContainers->sizeDiagScaling = 0;
    matContainers->DiagScaling = NULL;
}

/****************************************************************************
    Matrix vector multiplication of sparse matrices.
    Operates ditectly on the solution matrix
****************************************************************************/
-(void)matrixVectorMultplyInSolution:(FEMSolution * __nonnull)solution multiplyVector:(double * __nonnull)u resultVector:(double * __nonnull)v {
    
    FEMMatrixCRS *crsMatrix;
    FEMMatrixBand *bandMatrix;
    
    switch (solution.matrix.format) {
        case MATRIX_CRS:
            crsMatrix = [[FEMMatrixCRS alloc] init];
            [crsMatrix matrixVectorMultiplyInGlobal:solution vector:u result:v];
            break;
            
        case MATRIX_BAND:
        case MATRIX_SBAND:
            bandMatrix = [[FEMMatrixBand alloc] init];
            [bandMatrix matrixVectorMultiplyInGlobal:solution vector:u result:v];
            break;
            
        case MATRIX_LIST:
            NSLog(@"FEMCore:matrixVectorMultplyInSolution: list matrix type not suppored in this method.\n");
            break;
    }
}

/****************************************************************************
    Matrix vector multiplication of sparse matrices.
****************************************************************************/
-(void)matrixVectorMultplyInMatrix:(FEMMatrix * __nonnull)matrix multiplyVector:(double * __nonnull)u resultVector:(double * __nonnull)v {
    
    FEMMatrixCRS *crsMatrix;
    FEMMatrixBand *bandMatrix;
    
    switch (matrix.format) {
        case MATRIX_CRS:
            crsMatrix = [[FEMMatrixCRS alloc] init];
            [crsMatrix matrixVectorMultiply:matrix vector:u result:v];
            break;
            
        case MATRIX_BAND:
        case MATRIX_SBAND:
            bandMatrix = [[FEMMatrixBand alloc] init];
            [bandMatrix matrixVectorMultiply:matrix vector:u result:v];
            break;
            
        case MATRIX_LIST:
            NSLog(@"FEMCore:matrixVectorMultplyInSolution: list matrix type not suppored in this method.\n");
            break;
    }
}

-(void)invalidateVariableInTopMesh:(NSArray * __nonnull)topMesh primaryMesh:(FEMMesh * __nonnull)primaryMesh name:(NSString * __nonnull)name model:(FEMModel * __nonnull)model {
    
    int i;
    BOOL onlySearch, found;
    FEMVariable *var = nil, *var1 = nil, *primVar = nil;
    FEMUtilities *utilities;
    NSString *tmpName;
    
    utilities = [[FEMUtilities alloc] init];
    onlySearch = YES;
    primVar = [utilities getVariableFrom:primaryMesh.variables model:model name:name onlySearch:&onlySearch maskName:NULL info:&found];
    if (primVar == nil) return;
    
    for (FEMMesh *mesh in topMesh) {
        if (primaryMesh != mesh) {
            var = [utilities getVariableFrom:mesh.variables model:model name:name onlySearch:&onlySearch maskName:NULL info:&found];
            if (var != nil) {
                var.valid = NO;
                var.primaryMesh = primaryMesh;
            }
            
            if (primVar.dofs > 1) {
                if (/* DISABLES CODE */ (NO)) {
                    if ([primVar.name isEqualToString:@"flow solution"] == YES) {
                        var1 = [utilities getVariableFrom:mesh.variables model:model name:@"velocity 1" onlySearch:&onlySearch maskName:NULL info:&found];
                        if (var1 != nil) {
                            var1.valid = NO;
                            var1.primaryMesh = primaryMesh;
                        }
                        var1 = [utilities getVariableFrom:mesh.variables model:model name:@"velocity 2" onlySearch:&onlySearch maskName:NULL info:&found];
                        if (var != nil) {
                            var1.valid = NO;
                            var1.primaryMesh = primaryMesh;
                        }
                        var1 = [utilities getVariableFrom:mesh.variables model:model name:@"velocity 3" onlySearch:&onlySearch maskName:NULL info:&found];
                        if (var != nil) {
                            var1.valid = NO;
                            var1.primaryMesh = primaryMesh;
                        }
                        var1 = [utilities getVariableFrom:mesh.variables model:model name:@"pressure" onlySearch:&onlySearch maskName:NULL info:&found];
                        if (var1 != nil) {
                            var1.valid = NO;
                            var1.primaryMesh = primaryMesh;
                        }
                        var1 = [utilities getVariableFrom:mesh.variables model:model name:@"surface" onlySearch:&onlySearch maskName:NULL info:&found];
                        if (var1 != nil) {
                            var1.valid = NO;
                            var1.primaryMesh = primaryMesh;
                        }
                    }
                } else {
                    for (i=1; i<=primVar.dofs; i++) {
                        tmpName = [utilities appendNameFromString:name component:&i];
                        var1 = [utilities getVariableFrom:mesh.variables model:model name:tmpName onlySearch:&onlySearch maskName:NULL info:&found];
                        if (var1 != nil) {
                            var1.valid = NO;
                            var1.primaryMesh = primaryMesh;
                        }
                    }
                }
            }
        }
    }
    
    primVar.valuesChanged = YES;
    if (primVar.dofs > 1) {
        if (/* DISABLES CODE */ (NO)) {
            if ([primVar.name isEqualToString:@"flow solution"] == YES) {
                var = [utilities getVariableFrom:primaryMesh.variables model:model name:@"surface" onlySearch:&onlySearch maskName:NULL info:&found];
                if (var != nil) var.valuesChanged = YES;
                var = [utilities getVariableFrom:primaryMesh.variables model:model name:@"pressure" onlySearch:&onlySearch maskName:NULL info:&found];
                if (var != nil) var.valuesChanged = YES;
                var = [utilities getVariableFrom:primaryMesh.variables model:model name:@"velocity 1" onlySearch:&onlySearch maskName:NULL info:&found];
                if (var != nil) var.valuesChanged = YES;
                var = [utilities getVariableFrom:primaryMesh.variables model:model name:@"velocity 2" onlySearch:&onlySearch maskName:NULL info:&found];
                if (var != nil) var.valuesChanged = YES;
                var = [utilities getVariableFrom:primaryMesh.variables model:model name:@"velocity 3" onlySearch:&onlySearch maskName:NULL info:&found];
                if (var != nil) var.valuesChanged = YES;
            }
        } else {
            for (i=1; i<=primVar.dofs; i++) {
                tmpName = [utilities appendNameFromString:name component:&i];
                var = [utilities getVariableFrom:primaryMesh.variables model:model name:tmpName onlySearch:&onlySearch maskName:NULL info:&found];
                if (var != nil) var.valuesChanged = YES;
            }
        }
    }
}

-(void)getPassiveBoundaryAtIndex:(int)bcID model:(FEMModel * __nonnull)model mesh:(FEMMesh * __nonnull)mesh solution:(FEMSolution * __nonnull)solution {
    
    int i, n=0, cnt, ind, sz;
    int *arr;
    BOOL findEdges, L1, L2;
    FEMMeshUtils *meshUtils;
    Element_t *faces = NULL, *tElements = NULL, *meshElements = NULL, *modelElements = NULL;
    
    meshUtils = [[FEMMeshUtils alloc] init];
    
    findEdges = NO;
    [meshUtils findEdgesForMesh:mesh findEdges:NO];
    switch (mesh.dimension) {
        case 2:
            faces = mesh.getEdges;
            n = mesh.numberOfEdges;
            break;
        case 3:
            faces = mesh.getFaces;
            n = mesh.numberOfFaces;
    }
    
    arr = intvec(0, n-1);
    cnt = 0;
    for (i=0; i<n; i++) {
        if (faces[i].BoundaryInfo->Right == NULL || faces[i].BoundaryInfo->Left == NULL) continue;
        
        L1 = [self FEMCore_checkPassiveElement:faces[i].BoundaryInfo->Right model:model solution:solution];
        L2 = [self FEMCore_checkPassiveElement:faces[i].BoundaryInfo->Left model:model solution:solution];
        
        if (L1 ^ L2) {
            arr[cnt] = i;
            cnt++;
        }
    }
    
    sz = mesh.numberOfBulkElements + mesh.numberOfBoundaryElements - mesh.numberOfPassiveBCs;
    if ((sz+cnt) > mesh.numberOfElements) {
        tElements = mesh.getElements;
        meshElements = mesh.getElements;
        meshElements = (Element_t*) malloc( sizeof(Element_t) * (sz+cnt) );
        modelElements = model.getElements;
        if (modelElements == tElements) modelElements = meshElements;
        
        for (i=0; i<sz; i++) {
            meshElements[i] = tElements[i];
        }
        
        // Fix boundary element parent pointers to use new array...
        for (i=0; i<mesh.numberOfBoundaryElements-mesh.numberOfPassiveBCs; i++) {
            ind = i + mesh.numberOfBulkElements;
            if (meshElements[ind].BoundaryInfo->Left != NULL) meshElements[ind].BoundaryInfo->Left = &meshElements[meshElements[ind].BoundaryInfo->Left->ElementIndex-1];
            if (meshElements[ind].BoundaryInfo->Right != NULL) meshElements[ind].BoundaryInfo->Right = &meshElements[meshElements[ind].BoundaryInfo->Right->ElementIndex-1];
        }
        
        //... likewise for faces (edges)
        for (i=0; i<n; i++) {
            if (faces[i].BoundaryInfo->Left != NULL) faces[i].BoundaryInfo->Left = &meshElements[faces[i].BoundaryInfo->Left->ElementIndex-1];
            if (faces[i].BoundaryInfo->Right != NULL) faces[i].BoundaryInfo->Right = &meshElements[faces[i].BoundaryInfo->Right->ElementIndex-1];
        }
        
        free(tElements);
    }
    
    for (i=0; i<cnt; i++) {
        meshElements[sz] = faces[arr[i]];
        meshElements[sz].copy = true;
        meshElements[sz].ElementIndex = sz+1;
        meshElements[sz].BoundaryInfo->Constraint = bcID;
        sz++;
    }
    mesh.numberOfBoundaryElements = mesh.numberOfBoundaryElements - mesh.numberOfPassiveBCs + cnt;
    mesh.numberOfPassiveBCs = cnt;
    if (modelElements == meshElements) model.numberOfBoundaryElements = mesh.numberOfBoundaryElements;
    free_ivector(arr, 0, n-1);
}

-(void)computeNodalWeightsInSolution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model weightBoundary:(BOOL)weightBoundary perm:(int * __nullable)perm sizePerm:(int * __nullable)sizePerm variableName:(NSString * __nullable)variableName {
    
    int i, e, m, n, t, elemStart, elemEnd, sizeIntPerm;
    int *intPerm = NULL, *localIndexes;
    double detJ;
    double *sol = NULL;
    BOOL any, found, stat;
    Nodes_t *elementNodes, *nodes = NULL;
    NSString *intVarName;
    FEMVariable *weightsVar = nil;
    FEMUtilities *utilities;
    FEMNumericIntegration *integration;
    Element_t *elements = NULL;
    GaussIntegrationPoints *IP = NULL;
    variableArraysContainer *variableContainers = NULL, *weightsVarContainers = NULL, *bufferContainers = NULL;
    
    if (variableName != nil) {
        intVarName = variableName;
    } else if (weightBoundary == YES) {
        intVarName = [[solution.variable canonicalizeName] stringByAppendingString:@" boundary weights"];
    } else {
        intVarName = [[solution.variable canonicalizeName] stringByAppendingString:@" weights"];
    }
    utilities = [[FEMUtilities alloc] init];
    weightsVar = [utilities getVariableFrom:solution.mesh.variables model:model name:intVarName onlySearch:NULL maskName:NULL info:&found];
    
    if (weightBoundary == YES) {
        elemStart = solution.mesh.numberOfBulkElements;
        elemEnd = solution.mesh.numberOfBulkElements + solution.mesh.numberOfBoundaryElements;
    } else {
        elemStart = 0;
        elemEnd = solution.mesh.numberOfBulkElements;
    }
    
    if (weightsVar == nil) {
        if (perm != NULL) {
            intPerm = perm;
            sizeIntPerm = *sizePerm;
        } else {
            variableContainers = solution.variable.getContainers;
            intPerm = variableContainers->Perm;
            sizeIntPerm = variableContainers->sizePerm;
        }
        if (intPerm != NULL) {
            n = max_array(intPerm, sizeIntPerm);
            sol = doublevec(0, n-1);
            memset( sol, 0.0, n*sizeof(double) );
            bufferContainers = allocateVariableContainer();
            bufferContainers->Values = sol;
            bufferContainers->sizeValues = n;
            bufferContainers->Perm = intPerm;
            bufferContainers->sizePerm = sizeIntPerm;
            [utilities addVariableTo:solution.mesh.variables mesh:solution.mesh solution:solution name:intVarName dofs:1 container:bufferContainers component:NO ifOutput:NULL ifSecondary:NULL type:NULL];
            sol = NULL;
            free(bufferContainers);
        } else {
            NSLog(@"FEMCore:computeNodalWeightsInSolution: permutation vector not present.\n");
            return;
        }
        weightsVar = [utilities getVariableFrom:solution.mesh.variables model:model name:intVarName onlySearch:NULL maskName:NULL info:&found];
    }
    
    weightsVarContainers = weightsVar.getContainers;
    if (weightsVarContainers->Values == NULL) {
        NSLog(@"FEMCore:computeNodalWeightsInSolution: solution vector not present.\n");
        return;
    }
    
    NSLog(@"FEMCore:computeNodalWeightsInSolution: computing weights for solution to variable %@.\n", intVarName);
    m = solution.mesh.maxElementNodes;
    elementNodes = (Nodes_t*)malloc(sizeof(Nodes_t));
    initNodes(elementNodes);
    elementNodes->x = doublevec(0, m-1);
    elementNodes->y = doublevec(0, m-1);
    elementNodes->z = doublevec(0, m-1);
    localIndexes = intvec(0, m-1);
    memset( weightsVarContainers->Values, 0.0, weightsVarContainers->sizeValues*sizeof(double) );
    
    integration = [[FEMNumericIntegration alloc] init];
    if ([integration allocation:solution.mesh] == NO) fatal("FEMCore:computeNodalWeightsInSolution", "Allocation error in computeNodalWeightsInSolution.");
    
    elements = solution.mesh.getElements;
    nodes = solution.mesh.getNodes;
    for (e=elemStart; e<elemEnd; e++) {
        n = elements[e].Type.NumberOfNodes;
        for (i=0; i<n; i++) {
            localIndexes[i] = intPerm[elements[e].NodeIndexes[i]];
        }
        any = NO;
        for (i=0; i<n; i++) {
            if (localIndexes[i] < 0) {
                any = YES;
                break;
            }
        }
        if (any == YES) continue;
        
        n = elements[e].Type.NumberOfNodes;
        for (i=0; i<n; i++) {
            elementNodes->x[i] = nodes->x[elements[e].NodeIndexes[i]];
            elementNodes->y[i] = nodes->y[elements[e].NodeIndexes[i]];
            elementNodes->z[i] = nodes->z[elements[e].NodeIndexes[i]];
        }
        
        IP = GaussQuadrature(&elements[e], NULL, NULL);
        for (t=0; t<IP->n; t++) {
            // Basis function values & derivatives at the integration point
            stat = [integration setBasisForElement:&elements[e] elementNodes:elementNodes  inMesh:solution.mesh firstEvaluationPoint:IP->u[t] secondEvaluationPoint:IP->v[t] thirdEvaluationPoint:IP->w[t] withBubbles:NO basisDegree:NULL];
            stat = [integration setMetricDeterminantForElement:&elements[e] elementNodes:elementNodes inMesh:solution.mesh firstEvaluationPoint:IP->u[t] secondEvaluationPoint:IP->v[t] thirdEvaluationPoint:IP->w[t]];
            detJ = integration.metricDeterminant;
            
            for (i=0; i<n; i++) {
                weightsVarContainers->Values[localIndexes[i]] = weightsVarContainers->Values[localIndexes[i]] + IP->s[t] * detJ * integration.basis[i];
            }
        }
    }
    
    [integration deallocation:solution.mesh];
    
    free_dvector( elementNodes->x, 0, m-1);
    free_dvector( elementNodes->y, 0, m-1);
    free_dvector( elementNodes->z, 0, m-1);
    free(elementNodes);
    free_ivector(localIndexes, 0, m-1);
    
    NSLog(@"FEMCore:computeNodalWeightsInSolution: all done.");
}

// Eliminates bubble degrees of freedom from a local linear system
-(void)nsCondensateStiff:(double * __nonnull * __nonnull)stiff force:(double * __nonnull)force numberOfNodes:(int)n numberOfBubbles:(int)nb dimension:(int)dim force1:(double * __nonnull)force1 {
    
    int i, j, m, p, cdofs[(dim+1)*n], bdofs[dim*nb];
    double fb[nb*dim], **kbb, kbl[nb*dim][n*(dim+1)], klb[n*(dim+1)][nb*dim], x[(dim+1)*n], y[(dim+1)*n];
    double b[nb*dim][n*(dim+1)], c[nb*dim][n*(dim+1)], cc[n*(dim+1)][n*(dim+1)];
    FEMLinearAlgebra *linearAlgebra;
    
    m = 0;
    for (p=0; p<n; p++) {
        for (i=0; i<dim+1; i++) {
            cdofs[m] = (dim+1)*p + i;
            m++;
        }
    }
    
    m = 0;
    for (p=0; p<nb; p++) {
        for (i=0; i<dim; i++) {
            bdofs[m] = (dim+1)*p + i + n*(dim+1);
            m++;
        }
    }
    
    kbb = doublematrix(0, (nb*dim)-1, 0, (nb*dim)-1);
    for (i=0; i<nb*dim; i++) {
        for (j=0; j<nb*dim; j++) {
            kbb[i][j] = stiff[bdofs[i]][bdofs[j]];
        }
        fb[i] = force[bdofs[i]];
    }

    for (i=0; i<nb*dim; i++) {
        for (j=0; j<(dim+1)*n; j++) {
             kbl[i][j] = stiff[bdofs[i]][cdofs[j]];
        }
    }
    
    for (i=0; i<(dim+1)*n; i++) {
        for (j=0; j<nb*dim; j++) {
            klb[i][j] = stiff[cdofs[i]][bdofs[j]];
        }
    }
    
    linearAlgebra = [[FEMLinearAlgebra alloc] init];
    [linearAlgebra invertMatrix:kbb ofSize:nb*dim];
    
    memset( y, 0.0, sizeof(y) );
    cblas_dgemv(CblasRowMajor, CblasNoTrans, nb*dim, nb*dim, 1.0, *kbb, nb*dim, fb, 1, 0.0, y, 1);
    memcpy(x, y, sizeof(double));
    cblas_dgemv(CblasRowMajor, CblasNoTrans, n*(dim+1), nb*dim, 1.0, (double *)klb, n*(dim+1), x, 1, 0.0, y, 1);
    for (i=0; i<(dim+1)*n; i++) {
        force[i] = force[i] - y[i];
    }
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nb*dim, n*(dim+1), nb*dim, 1.0, *kbb, nb*dim, (double *)kbl, nb*dim, 0.0,  (double *)c, nb*dim);
    memcpy(*b, *c, ((nb*dim)*(n*(dim+1)))*sizeof(double));
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n*(dim+1), n*(dim+1), nb*dim, 1.0, (double *)klb, n*(dim+1), (double *)b, nb*dim, 0.0,  (double *)cc, n*(dim+1));
    for (i=0; i<(dim+1)*n; i++) {
        for (j=0; j<(dim+1)*n; j++) {
            stiff[i][j] = stiff[i][j] - cc[i][j];
        }
    }
    
    for (i=0; i<dim*nb; i++) {
        fb[i] = force1[bdofs[i]];
    }
    memset( y, 0.0, sizeof(y) );
    cblas_dgemv(CblasRowMajor, CblasNoTrans, nb*dim, nb*dim, 1.0, *kbb, nb*dim, fb, 1, 0.0, y, 1);
    memcpy(x, y, sizeof(double));
    cblas_dgemv(CblasRowMajor, CblasNoTrans, n*(dim+1), nb*dim, 1.0, (double *)klb, n*(dim+1), x, 1, 0.0, y, 1);
    for (i=0; i<(dim+1)*n; i++) {
        force1[i] = force1[i] - y[i];
    }
    free_dmatrix(kbb, 0, (nb*dim)-1, 0, (nb*dim)-1);
}

-(void)condensateStiff:(double * __nonnull * __nonnull)stiff force:(double * __nonnull)force numberOfNodes:(int)n force1:(double * __nullable)force1 {
    
    int i, j, ldofs[n], bdofs[n];
    double fb[n], b[n][n], c[n][n], **kbb, kbl[n][n], klb[n][n], x[n], y[n];
    FEMLinearAlgebra *linearAlgebra;
    
    for (i=0; i<n; i++) {
        ldofs[i] = i;
        bdofs[i] = ldofs[i] + n;
    }
    
    kbb = doublematrix(0, n-1, 0, n-1);
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            kbb[i][j] = stiff[bdofs[i]][bdofs[j]];
            kbl[i][j] = stiff[bdofs[i]][ldofs[j]];
            klb[i][j] = stiff[ldofs[i]][bdofs[j]];
        }
        fb[i] = force[bdofs[i]];
    }
    
    linearAlgebra = [[FEMLinearAlgebra alloc] init];
    [linearAlgebra invertMatrix:kbb ofSize:n];
    
    memset( y, 0.0, sizeof(y) );
    cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, *kbb, n, fb, 1, 0.0, y, 1);
    memcpy(x, y, sizeof(double));
    cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, (double *)klb, n, x, 1, 0.0, y, 1);
    for (i=0; i<n; i++) {
        force[i] = force[i] - y[i];
    }
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, *kbb, n, (double *)kbl, n, 0.0,  (double *)c, n);
    memcpy(*b, *c, (n*n)*sizeof(double));
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, (double *)klb, n, (double *)b, n, 0.0,  (double *)c, n);
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            stiff[i][j] = stiff[i][j] - c[i][j];
        }
    }
    
    if (force1 != NULL) {
        for (i=0; i<n; i++) {
            fb[i] = force1[bdofs[i]];
        }
        memset( y, 0.0, sizeof(y) );
        cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, *kbb, n, fb, 1, 0.0, y, 1);
        memcpy(x, y, sizeof(double));
        cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, (double *)klb, n, x, 1, 0.0, y, 1);
        for (i=0; i<n; i++) {
            force1[i] = force1[i] - y[i];
        }
    }
    free_dmatrix(kbb, 0, n-1, 0, n-1);
}

/*************************************************************************************************************************
    For time dependent simulations add the time derivative coefficient terms to the matrix containing other coefficients.
*************************************************************************************************************************/
-(void)addFirstOrderTimeModel:(FEMModel * __nonnull)model solution:(FEMSolution * __nonnull)solution element:(Element_t * __nonnull)element massMatrix:(double * __nonnull * __nonnull)massMatrix stiffMatrix:(double * __nonnull * __nonnull)stiffMatrix force:(double * __nonnull)force dt:(double)dt size:(int)n dofs:(int)dofs nodeIndexes:(int * __nonnull)nodeIndexes rows:(int * __nonnull)rows cols:(int * __nonnull)cols timeIntegration:(FEMTimeIntegration * __nonnull)timeIntegration utilities:(FEMUtilities * __nonnull)utilities {
    
    int i, j, k, l, m, order;
    double s, t;
    double **prevSol, lForce[dofs*n], buffer[dofs*n];
    double dts[solution.order];
    BOOL constantDt, found;
    NSString *method;
    FEMVariable *dtVar;
    matrixArraysContainer *matContainers = NULL;
    variableArraysContainer *varContainers = NULL, *containers = NULL;
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
    
    prevSol = doublematrix(0, (dofs*n)-1, 0, solution.order-1);
    
    if (solution.matrix.isLumped == YES) {
        s = 0.0;
        t = 0.0;
        for (i=0; i<n*dofs; i++) {
            for (j=0; j<n*dofs; j++) {
                s = s + massMatrix[i][j];
                if (i != j) massMatrix[i][j] = 0.0;
            }
            t = t + massMatrix[i][i];
        }
        
        for (i=0; i<n; i++) {
            for (j=0; j<dofs; j++) {
                k = dofs * i + j;
                if (t != 0.0) massMatrix[k][k] = massMatrix[k][k] * s / t;
            }
        }
    }
    
    order = min(solution.doneTime, solution.order);
    
    for (i=0; i<n; i++) {
        for (j=0; j<dofs; j++) {
            k = dofs * i + j;
            l = dofs * nodeIndexes[i] + j;
            for (m=0; m<order; m++) {
                prevSol[k][m] = varContainers->PrevValues[l][m];
            }
        }
    }
    
    for (i=0; i<n*dofs; i++) {
        lForce[i] = force[i];
    }
    
    [self FEMCore_updateGlobalForceModel:model solution:solution element:element forceVector: matContainers->Force forceVectorUpdateAtIndex:0 localForce:lForce size:n dofs:dofs nodeIndexes:nodeIndexes rotateNT:NULL];
    
    if ((solution.solutionInfo)[@"time stepping method"] != nil) {
        method = [NSString stringWithString:(solution.solutionInfo)[@"time stepping method"]];
    }
    
    if ([method isEqualToString:@"fs"] == YES) {
        
        for (i=0; i<dofs*n; i++) {
            buffer[i] = prevSol[i][0];
        }
        [timeIntegration fractionalStepInSolution:solution numberOfNodes:n*dofs dt:dt massMatrix:massMatrix stiffMatrix:stiffMatrix force:force prevSolution:buffer rows:rows];
    }
    else if ([method isEqualToString:@"bdf"]) {
        dts[0] = dt;
        constantDt = YES;
        if (order > 1) {
            dtVar = [utilities getVariableFrom:solution.mesh.variables model:model name:@"time step size" onlySearch:NULL maskName:NULL info:&found];
            containers = dtVar.getContainers;
            for (i=1; i<order; i++) {
                dts[i] = containers->PrevValues[0][i-1];
                if ( fabs(dts[i]-dts[0]) > 1.0e-6 * dts[0] ) constantDt = NO;
            }
            containers = NULL;
        }
        if (constantDt == YES) {
            [timeIntegration bdfLocalInSolution:solution numberOfNodes:n*dofs dt:dt massMatrix:massMatrix stiffMatrix:stiffMatrix force:force prevSolution:prevSol order:order rows:rows cols:cols];
        } else {
            [timeIntegration vbdfLocalInSolution:solution numberOfNodes:n*dofs dts:dts massMatrix:massMatrix stiffMatrix:stiffMatrix force:force prevSolution:prevSol order:order rows:rows cols:cols];
        }
    }
    else {
        for (i=0; i<dofs*n; i++) {
            buffer[i] = prevSol[i][0];
        }
        [timeIntegration newMarkBetaInSolution:solution numberOfNodes:n*dofs dt:dt massMatrix:massMatrix stiffMatrix:stiffMatrix force:force prevSolution:buffer beta:solution.beta rows:rows cols:cols];
    }
    
    free_dmatrix(prevSol, 0, (dofs*n)-1, 0, solution.order-1);
}

/*************************************************************************************************************************
    For time dependent simulations add the time derivative coefficient terms to the matrix containing other coefficients.
*************************************************************************************************************************/
-(void)addSecondOrderTimeModel:(FEMModel * __nonnull)model solution:(FEMSolution * __nonnull)solution element:(Element_t * __nonnull)element massMatrix:(double * __nonnull * __nonnull)massMatrix dampMatrix:(double * __nonnull * __nonnull)dampMatrix stiffMatrix:(double * __nonnull * __nonnull)stiffMatrix force:(double * __nonnull)force dt:(double)dt size:(int)n dofs:(int)dofs nodeIndexes:(int * __nonnull)nodeIndexes rows:(int * __nonnull)rows cols:(int * __nonnull)cols timeIntegration:(FEMTimeIntegration * __nonnull)timeIntegration {
    
    int i, j, k, l;
    double s, t;
    double lForce[dofs*n];
    double x[dofs*n], v[dofs*n], a[dofs*n];;
    matrixArraysContainer *matContainers = NULL;
    variableArraysContainer *varContainers = NULL;
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
    
    if (solution.matrix.isLumped == YES) {
        s = 0.0;
        t = 0.0;
        for (i=0; i<n*dofs; i++) {
            for (j=0; j<n*dofs; j++) {
                s = s + massMatrix[i][j];
                if (i != j) massMatrix[i][j] = 0.0;
            }
            t = t + massMatrix[i][i];
        }
        
        for (i=0; i<n; i++) {
            for (j=0; j<dofs; j++) {
                k = dofs * i + j;
                if (t != 0.0) massMatrix[k][k] = massMatrix[k][k] * s / t;
            }
        }
        
        s = 0.0;
        t = 0.0;
        for (i=0; i<n*dofs; i++) {
            for (j=0; j<n*dofs; j++) {
                s = s + dampMatrix[i][j];
                if (i != j) dampMatrix[i][j] = 0.0;
            }
            t = t + dampMatrix[i][i];
        }
        
        for (i=0; i<n; i++) {
            for (j=0; j<dofs; j++) {
                k = dofs * i + j;
                if (t != 0.0) dampMatrix[k][k] = dampMatrix[k][k] * s / t;
            }
        }
    }
    
    // Get previous solution vectors and update current force
    for (i=0; i<n; i++) {
        for (j=0; j<dofs; j++) {
            k = dofs * i + j;
            if (nodeIndexes[i] >= 0) {
                l = dofs * nodeIndexes[i] + j;
                x[k] = varContainers->PrevValues[l][2];
                v[k] = varContainers->PrevValues[l][3];
                a[k] = varContainers->PrevValues[l][4];
            }
        }
    }
    
    for (i=0; i<n*dofs; i++) {
        lForce[i] = force[i];
    }
    
    [self FEMCore_updateGlobalForceModel:model solution:solution element:element forceVector: matContainers->Force forceVectorUpdateAtIndex:0 localForce:lForce size:n dofs:dofs nodeIndexes:nodeIndexes rotateNT:NULL];
    
    [timeIntegration bossakSecondOrder:solution numberOfNodes:n*dofs dt:dt massMatrix:massMatrix dampMatrix:dampMatrix stiffMatrix:stiffMatrix force:force prevSolution1:x prevSolution2:v prevSolution3:a alpha:solution.alpha rows:rows cols:cols];
}

#pragma mark First order time

/********************************************************************************************************************
    Update the mass matrix only
 
    Arguments:
 
    FEMSolution *solution       -> Class solution which contains the global matrix
    double **localStiffMatrix   -> Local matrix to be added tp the global matrix
    int n                       -> Number of nodes in element
    int dofs                    -> Number of dofs per node
    int *nodeIndexes            -> Element node to global node numbering mapping
********************************************************************************************************************/
-(void)updateMassMatrixModel:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution localMassMatrix:(double * __nonnull * __nonnull)localMassMatrix element:(Element_t * __nonnull)element numberOfNodes:(int)n dofs:(int)dofs nodeIndexes:(int * __nonnull)nodeIndexes {
    
    int i, j;
    double s, t, *saveValues;
    FEMMatrixCRS *crsMatrix;
    FEMMatrixBand *bandMatrix;
    matrixArraysContainer *matContainers = NULL;
    
    // Check first if this element has been defined passive
    if ([self FEMCore_checkPassiveElement:element model:model solution:solution] == YES) return;
    
    // Update global matrix and rhs vector
    
    if (solution.matrix.lumped == YES) {
        s = 0.0;
        t = 0.0;
        for (i=0; i<n*dofs; i++) {
            for (j=0; j<n*dofs; j++) {
                s = s + localMassMatrix[i][j];
                if (i != j) localMassMatrix[i][j] = 0.0;
            }
            t = t + localMassMatrix[i][i];
        }
        
        for (i=0; i<n*dofs; i++) {
            localMassMatrix[i][i] = localMassMatrix[i][i] * s / t;
        }
    }
    
    matContainers = solution.matrix.getContainers;
    saveValues = matContainers->Values;
    matContainers->Values = matContainers->MassValues;
    
    switch (solution.matrix.format) {
        case MATRIX_CRS:
            crsMatrix = [[FEMMatrixCRS alloc] init];
            [crsMatrix glueLocalMatrix:localMassMatrix inGlobal:solution numberOfNodes:n dofs:dofs indexes:nodeIndexes];
            break;
            
        case MATRIX_LIST:
            // TODO: implement for the List matrix case
            break;
            
        case MATRIX_BAND:
        case MATRIX_SBAND:
            bandMatrix = [[FEMMatrixBand alloc] init];
            [bandMatrix glueLocalMatrix:localMassMatrix inGlobal:solution numberOfNodes:n dofs:dofs indexes:nodeIndexes];
            break;
    }
    
    matContainers->Values = saveValues;
}

-(void)defaultUpdateMass:(double * __nonnull * __nonnull)mass element:(Element_t * __nonnull)element solution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model {
    
    int i, n;
    matrixArraysContainer *matContainers = NULL;
    variableArraysContainer *varContainers = NULL;
    
    memset( _indexStore, -1, _sizeIndexStore*sizeof(int) );
    n = [self getElementDofsSolution:solution model:model forElement:element atIndexes:_indexStore disableDiscontinuousGalerkin:NULL];
    
    // TODO: Add support for parallel run
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
    
    if (matContainers->MassValues == NULL) {
        matContainers->MassValues = doublevec(0, matContainers->sizeValues-1);
        matContainers->sizeMassValues = matContainers->sizeValues;
        memset( matContainers->MassValues, 0.0, matContainers->sizeMassValues*sizeof(double) );
    }
    
    int indexes[n];
    for (i=0; i<n; i++) {
        indexes[i] = varContainers->Perm[_indexStore[i]];
    }
    [self updateMassMatrixModel:model inSolution:solution localMassMatrix:mass element:element numberOfNodes:n dofs:solution.variable.dofs nodeIndexes:indexes];
}

-(void)defaultUpdateComplexMass:(double complex * __nonnull * __nonnull)cmass element:(Element_t * __nonnull)element solution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model {
    
    int i, j, dofs, n;
    double **mass;
    matrixArraysContainer *matContainers = NULL;
    variableArraysContainer *varContainers = NULL;
    
    dofs = solution.variable.dofs;
    memset( _indexStore, -1, _sizeIndexStore*sizeof(int) );
    n = [self getElementDofsSolution:solution model:model forElement:element atIndexes:_indexStore disableDiscontinuousGalerkin:NULL];
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
    
    // TODO: Add support for parallel run
    
    if (matContainers->MassValues == NULL) {
        matContainers->MassValues = doublevec(0, matContainers->sizeValues-1);
        matContainers->sizeMassValues = matContainers->sizeValues;
        memset( matContainers->MassValues, 0.0, matContainers->sizeMassValues*sizeof(double) );
    }
    
    mass = doublematrix(0, (dofs*n)-1, 0, (dofs*n)-1);
    for (i=0; i<n*dofs/2; i++) {
        for (j=0; j<n*dofs/2; j++) {
            mass[2*i][2*j]      = creal(cmass[i][j]);
            mass[2*i][2*j+1]    = -cimag(cmass[i][j]);
            mass[2*i+1][2*j]    = cimag(cmass[i][j]);
            mass[2*i+1][2*j+1]  = creal(cmass[i][j]);
        }
    }
    
    int indexes[n];
    for (i=0; i<n; i++) {
        indexes[i] = varContainers->Perm[_indexStore[i]];
    }
    [self updateMassMatrixModel:model inSolution:solution localMassMatrix:mass element:element numberOfNodes:n dofs:solution.variable.dofs nodeIndexes:indexes];
    free_dmatrix(mass, 0, (dofs*n)-1, 0, (dofs*n)-1);
}

-(void)defaultUpdateDamp:(double * __nonnull * __nonnull)damp element:(Element_t * __nonnull)element solution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model {
    
    int i, n;
    double *saveValues;
    matrixArraysContainer *matContainers = NULL;
    variableArraysContainer *varContainers = NULL;
    
    memset( _indexStore, -1, self.sizeIndexStore*sizeof(int) );
    n = [self getElementDofsSolution:solution model:model forElement:element atIndexes:self.indexStore disableDiscontinuousGalerkin:NULL];
    
    // TODO: Add support for parallel run
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
    
    if (matContainers->DampValues == NULL) {
        matContainers->DampValues = doublevec(0, matContainers->sizeValues-1);
        matContainers->sizeDampValues = matContainers->sizeValues;
        memset( matContainers->DampValues, 0.0, matContainers->sizeDampValues*sizeof(double) );
    }
    
    saveValues = matContainers->MassValues;
    matContainers->MassValues = matContainers->DampValues;
    int indexes[n];
    for (i=0; i<n; i++) {
        indexes[i] = varContainers->Perm[_indexStore[i]];
    }
    [self updateMassMatrixModel:model inSolution:solution localMassMatrix:damp element:element numberOfNodes:n dofs:solution.variable.dofs nodeIndexes:indexes];
    matContainers->MassValues = saveValues;
}

-(void)defaultUpdateComplexDamp:(double complex * __nonnull * __nonnull)cdamp element:(Element_t * __nonnull)element solution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model {
    
    int i, j, dofs, n;
    double **damp, *saveValues;
    matrixArraysContainer *matContainers = NULL;
    variableArraysContainer *varContainers = NULL;
    
    dofs = solution.variable.dofs;
    memset( _indexStore, -1, self.sizeIndexStore*sizeof(int) );
    n = [self getElementDofsSolution:solution model:model forElement:element atIndexes:_indexStore disableDiscontinuousGalerkin:NULL];
    
    // TODO: Add support for parallel run
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
    
    if (matContainers->DampValues == NULL) {
        matContainers->DampValues = doublevec(0, matContainers->sizeValues-1);
        matContainers->sizeDampValues = matContainers->sizeValues;
        memset( matContainers->DampValues, 0.0, matContainers->sizeDampValues*sizeof(double) );
    }
    
    damp = doublematrix(0, (dofs*n)-1, 0, (dofs*n)-1);
    for (i=0; i<n*dofs/2; i++) {
        for (j=0; j<n*dofs/2; j++) {
            damp[2*i][2*j]      = creal(cdamp[i][j]);
            damp[2*i][2*j+1]    = -cimag(cdamp[i][j]);
            damp[2*i+1][2*j]    = cimag(cdamp[i][j]);
            damp[2*i+1][2*j+1]  = creal(cdamp[i][j]);
        }
    }
    
    saveValues = matContainers->MassValues;
    matContainers->MassValues = matContainers->DampValues;
    int indexes[n];
    for (i=0; i<n; i++) {
        indexes[i] = varContainers->Perm[_indexStore[i]];
    }
    [self updateMassMatrixModel:model inSolution:solution localMassMatrix:damp element:element numberOfNodes:n dofs:solution.variable.dofs nodeIndexes:indexes];
    matContainers->MassValues = saveValues;
    free_dmatrix(damp, 0, (dofs*n)-1, 0, (dofs*n)-1);
}

-(void)defaultFirstOrderTime:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution forElement:(Element_t * __nonnull)element realMass:(double * __nonnull * __nonnull)mass realStiff:(double * __nonnull * __nonnull)stiff realForce:(double * __nonnull)force stiffRows:(int * __nonnull)rows stiffCols:(int * __nonnull)cols timeIntegration:(FEMTimeIntegration * __nonnull)timeIntegration utilities:(FEMUtilities * __nonnull)utilities { // rows and cols for stiff and mass matrix
    
    int i, n;
    double dt;
    variableArraysContainer *varContainers = NULL;
    
    if ([(solution.solutionInfo)[@"use global mass matrix"] boolValue] == YES) {
        [self defaultUpdateMass:mass element:element solution:solution model:model];
        return;
    }
    
    dt = solution.dt;
    memset( _indexStore, -1, _sizeIndexStore*sizeof(int) );
    n = [self getElementDofsSolution:solution model:model forElement:element atIndexes:_indexStore disableDiscontinuousGalerkin:NULL];
    varContainers = solution.variable.getContainers;
    
    int perm[n];
    for (i=0; i<n; i++) {
        perm[i] = varContainers->Perm[_indexStore[i]];
    }
    
    [self addFirstOrderTimeModel:model solution:solution element:element massMatrix:mass stiffMatrix:stiff force:force dt:dt size:n dofs:solution.variable.dofs nodeIndexes:perm rows:rows cols:cols timeIntegration:timeIntegration utilities:utilities];
}

-(void)defaultFirstOrderTime:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution forElement:(Element_t * __nonnull)element complexMass:(double complex * __nonnull * __nonnull)cmass complexStiff:(double complex * __nonnull * __nonnull)cstiff complexForce:(double complex * __nonnull)cforce stiffRows:(int * __nonnull)rows stiffCols:(int * __nonnull)cols timeIntegration:(FEMTimeIntegration * __nonnull)timeIntegration utilities:(FEMUtilities * __nonnull)utilities {
    
    int i, j, n, dofs;
    double dt;
    double **mass, **stiff, *force;
    variableArraysContainer *varContainers = NULL;
    
    if ([(solution.solutionInfo)[@"use global mass matrix"] boolValue] == YES) {
        [self defaultUpdateComplexMass:cmass element:element solution:solution model:model];
        return;
    }
    
    dt = solution.dt;
    dofs = solution.variable.dofs;
    memset( _indexStore, -1, _sizeIndexStore*sizeof(int) );
    n = [self getElementDofsSolution:solution model:model forElement:element atIndexes:_indexStore disableDiscontinuousGalerkin:NULL];
    varContainers = solution.variable.getContainers;
    
    mass = doublematrix(0, (n*dofs)-1, 0, (n*dofs)-1);
    stiff = doublematrix(0, (n*dofs)-1, 0, (n*dofs)-1);
    force = doublevec(0, (n*dofs)-1);
    
    for (i=0; i<n*dofs/2; i++) {
        force[2*i] = creal(cforce[i]);
        force[2*i+1] = cimag(cforce[i]);
        for (j=0; j<n*dofs; j++) {
            mass[2*i][2*j]      = creal(cmass[i][j]);
            mass[2*i][2*j+1]    = -cimag(cmass[i][j]);
            mass[2*i+1][2*j]    = cimag(cmass[i][j]);
            mass[2*i+1][2*j+1]  = creal(cmass[i][j]);
            stiff[2*i][2*j]     = creal(cstiff[i][j]);
            stiff[2*i][2*j+1]   = -cimag(cstiff[i][j]);
            stiff[2*i+1][2*j]   = cimag(cstiff[i][j]);
            stiff[2*i+1][2*j+1] = creal(cstiff[i][j]);
        }
    }
    
    int perm[n];
    for (i=0; i<n; i++) {
        perm[i] = varContainers->Perm[_indexStore[i]];
    }

    [self addFirstOrderTimeModel:model solution:solution element:element massMatrix:mass stiffMatrix:stiff force:force dt:dt size:n dofs:solution.variable.dofs nodeIndexes:perm rows:rows cols:cols timeIntegration:timeIntegration utilities:utilities];
    
    for (i=0; i<n*dofs/2; i++) {
        cforce[i] = force[2*i] + force[2*i+1]*I;
        for (j=0; j<n*dofs/2; j++) {
            cmass[i][j] = mass[2*i][2*j] + (-mass[2*i][2*j+1]*I);
            cstiff[i][j] = stiff[2*i][2*j] + (-stiff[2*i][2*j+1]*I);
        }
    }
    
    free_dmatrix(mass, 0, (n*dofs)-1, 0, (n*dofs)-1);
    free_dmatrix(stiff, 0, (n*dofs)-1, 0, (n*dofs)-1);
    free_dvector(force, 0, (n*dofs)-1);
}

-(void)defaultSecondOrderTime:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution forElement:(Element_t * __nonnull)element realMass:(double * __nonnull * __nonnull)mass realDamp:(double * __nonnull * __nonnull)damp realStiff:(double * __nonnull * __nonnull)stiff realForce:(double * __nonnull)force stiffRows:(int * __nonnull)rows stiffCols:(int * __nonnull)cols timeIntegration:(FEMTimeIntegration * __nonnull)timeIntegration {
    
    int i, n;
    double dt;
    variableArraysContainer *varContainers = NULL;
    
    if ([(solution.solutionInfo)[@"use global mass matrix"] boolValue] == YES) {
        [self defaultUpdateMass:mass element:element solution:solution model:model];
        [self defaultUpdateDamp:damp element:element solution:solution model:model];
        return;
    }
    
    dt = solution.dt;
    memset( _indexStore, -1, _sizeIndexStore*sizeof(int) );
    n = [self getElementDofsSolution:solution model:model forElement:element atIndexes:_indexStore disableDiscontinuousGalerkin:NULL];
    varContainers = solution.variable.getContainers;
    
    int perm[n];
    for (i=0; i<n; i++) {
        perm[i] = varContainers->Perm[_indexStore[i]];
    }
    
    [self addSecondOrderTimeModel:model solution:solution element:element massMatrix:mass dampMatrix:damp stiffMatrix:stiff force:force dt:dt size:n dofs:solution.variable.dofs nodeIndexes:perm rows:rows cols:cols timeIntegration:timeIntegration];
}

-(void)defaultSecondOrderTime:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution forElement:(Element_t * __nonnull)element complexMass:(double complex * __nonnull * __nonnull)cmass complexDamp:(double complex * __nonnull * __nonnull)cdamp complexStiff:(double complex * __nonnull * __nonnull)cstiff complexForce:(double complex * __nonnull)cforce stiffRows:(int * __nonnull)rows stiffCols:(int * __nonnull)cols timeIntegration:(FEMTimeIntegration * __nonnull)timeIntegration {
    
    int i, j, n, dofs;
    double dt;
    double **damp, **mass, **stiff, *force;
    variableArraysContainer *varContainers = NULL;
    
    if ([(solution.solutionInfo)[@"use global mass matrix"] boolValue] == YES) {
        [self defaultUpdateComplexMass:cmass element:element solution:solution model:model];
        [self defaultUpdateComplexDamp:cdamp element:element solution:solution model:model];
        return;
    }
    
    dt = solution.dt;
    dofs = solution.variable.dofs;
    memset( _indexStore, -1, _sizeIndexStore*sizeof(int) );
    n = [self getElementDofsSolution:solution model:model forElement:element atIndexes:_indexStore disableDiscontinuousGalerkin:NULL];
    varContainers = solution.variable.getContainers;
    
    mass = doublematrix(0, (n*dofs)-1, 0, (n*dofs)-1);
    stiff = doublematrix(0, (n*dofs)-1, 0, (n*dofs)-1);
    damp = doublematrix(0, (n*dofs)-1, 0, (n*dofs)-1);
    force = doublevec(0, (n*dofs)-1);
    
    for (i=0; i<n*dofs/2; i++) {
        force[2*i] = creal(cforce[i]);
        force[2*i+1] = cimag(cforce[i]);
        for (j=0; j<n*dofs; j++) {
            mass[2*i][2*j]      = creal(cmass[i][j]);
            mass[2*i][2*j+1]    = -cimag(cmass[i][j]);
            mass[2*i+1][2*j]    = cimag(cmass[i][j]);
            mass[2*i+1][2*j+1]  = creal(cmass[i][j]);
            damp[2*i][2*j]      = creal(cdamp[i][j]);
            damp[2*i][2*j+1]    = -cimag(cdamp[i][j]);
            damp[2*i+1][2*j]    = cimag(cdamp[i][j]);
            damp[2*i+1][2*j+1]  = creal(cdamp[i][j]);
            stiff[2*i][2*j]     = creal(cstiff[i][j]);
            stiff[2*i][2*j+1]   = -cimag(cstiff[i][j]);
            stiff[2*i+1][2*j]   = cimag(cstiff[i][j]);
            stiff[2*i+1][2*j+1] = creal(cstiff[i][j]);
        }
    }
    
    int perm[n];
    for (i=0; i<n; i++) {
        perm[i] = varContainers->Perm[_indexStore[i]];
    }
    
    [self addSecondOrderTimeModel:model solution:solution element:element massMatrix:mass dampMatrix:damp stiffMatrix:stiff force:force dt:dt size:n dofs:solution.variable.dofs nodeIndexes:perm rows:rows cols:cols timeIntegration:timeIntegration];
    
    for (i=0; i<n*dofs/2; i++) {
        cforce[i] = force[2*i] + force[2*i+1]*I;
        for (j=0; j<n*dofs/2; j++) {
            cmass[i][j] = mass[2*i][2*j] + (-mass[2*i][2*j+1]*I);
            cdamp[i][j] = damp[2*i][2*j] + (-damp[2*i][2*j+1]*I);
            cstiff[i][j] = stiff[2*i][2*j] + (-stiff[2*i][2*j+1]*I);
        }
    }
    
    free_dmatrix(mass, 0, (n*dofs)-1, 0, (n*dofs)-1);
    free_dmatrix(stiff, 0, (n*dofs)-1, 0, (n*dofs)-1);
    free_dmatrix(damp, 0, (n*dofs)-1, 0, (n*dofs)-1);
    free_dvector(force, 0, (n*dofs)-1);
}

-(void)defaultFirstOrderTimeGlobalModel:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution timeIntegration:(FEMTimeIntegration * __nonnull)timeIntegration utilities:(FEMUtilities * __nonnull)utilities {
    
    int i, j, k, n, order;
    double force[1], dts[16];
    BOOL constantDt=YES, found, hasMass;
    NSString *method;
    FEMVariable *dtVar;
    variableArraysContainer *varContainers = NULL, *dtVarContainers = NULL;
    matrixArraysContainer *matContainers = NULL;

    varContainers = solution.variable.getContainers;
    matContainers = solution.matrix.getContainers;
    if (_saveValues != varContainers->Values) {
        if (_stiff != NULL) {
            free_dmatrix(_stiff, 0, 0, 0, _n1-1);
            free_dmatrix(_mass, 0, 0, 0, _n1-1);
            free_dmatrix(_x, 0, _n1-1, 0, _k1-1);
        }
        _n1 = 0;
        for (i=0; i<solution.matrix.numberOfRows; i++) {
            _n1 = max(_n1, (matContainers->Rows[i+1]-matContainers->Rows[i]));
        }
        _k1 = varContainers->size2PrevValues;
        _stiff = doublematrix(0, 0, 0, _n1-1);
        _mass = doublematrix(0, 0, 0, _n1-1);
        _x = doublematrix(0, _n1-1, 0, _k1-1);
        
        _saveValues = varContainers->Values;
    }

    memset( *_stiff, 0.0, (1*_n1)*sizeof(double) );
    memset( *_mass, 0.0, (1*_n1)*sizeof(double) );
    memset( *_x, 0.0, (_n1*_k1)*sizeof(double) );
    
    order = min(solution.doneTime, solution.order);
    hasMass = (matContainers->MassValues != NULL) ? YES : NO;
    
    method = (solution.solutionInfo)[@"time stepping method"];
    if ([method isEqualToString:@"bdf"] == YES) {
        dts[0] = solution.dt;
        if (order > 1) {
            dtVar = [utilities getVariableFrom:solution.mesh.variables model:model name:@"time step size" onlySearch:NULL maskName:NULL info:&found];
            dtVarContainers = dtVar.getContainers;
            for (i=1; i<order; i++) {
                dts[i] = dtVarContainers->PrevValues[0][i-1];
                if (fabs(dts[i]-dts[0]) > 1.0e-6 * dts[0]) constantDt = NO;
            }
        }
    }
    
    int rows = 1;
    for (i=0; i<solution.matrix.numberOfRows; i++) {
        n = 0;
        for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
            _stiff[0][n] = matContainers->Values[j];
            if (hasMass == YES) _mass[0][n] = matContainers->MassValues[j];
            for (k=0; k<varContainers->size2PrevValues; k++) {
                _x[n][k] = varContainers->PrevValues[matContainers->Cols[j]][k];
            }
            n++;
        }
        force[0] = matContainers->RHS[i];
        matContainers->Force[i][0] = force[0];
        
        if ([method isEqualToString:@"fs"] == YES) {
            double *prevSol = doublevec(0, _n1-1);
            for (j=0; j<_n1; j++) {
                prevSol[j] = _x[j][0];
            }
            [timeIntegration fractionalStepInSolution:solution numberOfNodes:n dt:solution.dt massMatrix:_mass stiffMatrix:_stiff force:force prevSolution:prevSol rows:&rows];
            free_dvector(prevSol, 0, _n1-1);
        } else if ([method isEqualToString:@"bdf"] == YES) {
            if (constantDt == YES) {
                [timeIntegration bdfLocalInSolution:solution numberOfNodes:n dt:solution.dt massMatrix:_mass stiffMatrix:_stiff force:force prevSolution:_x order:order rows:&rows cols:&_n1];
            } else {
                [timeIntegration vbdfLocalInSolution:solution numberOfNodes:n dts:dts massMatrix:_mass stiffMatrix:_stiff force:force prevSolution:_x order:order rows:&rows cols:&_n1];
            }
        } else {
            double *prevSol = doublevec(0, _n1-1);
            for (j=0; j<_n1; j++) {
                prevSol[j] = _x[j][0];
            }
            [timeIntegration newMarkBetaInSolution:solution numberOfNodes:n dt:solution.dt massMatrix:_mass stiffMatrix:_stiff force:force prevSolution:prevSol beta:solution.beta rows:&rows cols:&_n1];
            free_dvector(prevSol, 0, _n1-1);
        }
        n = 0;
        for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
            matContainers->Values[j] = _stiff[0][n];
            n++;
        }
        matContainers->RHS[i] = force[0];
    }
}

-(void)defaultSecondOrderTimeGlobalModel:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution timeIntegration:(FEMTimeIntegration * __nonnull)timeIntegration {
    
    int i, j, k, n;
    double force[1];
    BOOL hasDamping, hasMass;
    variableArraysContainer *varContainers = NULL;
    matrixArraysContainer *matContainers = NULL;
    
    varContainers = solution.variable.getContainers;
    matContainers = solution.matrix.getContainers;
    if (_saveValues != varContainers->Values) {
        if (_stiff != NULL) {
            free_dmatrix(_stiff, 0, 0, 0, _n1-1);
            free_dmatrix(_mass, 0, 0, 0, _n1-1);
            free_dmatrix(_damp, 0, 0, 0, _n1-1);
            free_dmatrix(_x, 0, _n1-1, 0, _k1-1);
        }
        _n1 = 0;
        for (i=0; i<solution.matrix.numberOfRows; i++) {
            _n1 = max(_n1, (matContainers->Rows[i+1]-matContainers->Rows[i]));
        }
        _k1 = varContainers->size2PrevValues;
        _stiff = doublematrix(0, 0, 0, _n1-1);
        _mass = doublematrix(0, 0, 0, _n1-1);
        _damp = doublematrix(0, 0, 0, _n1-1);
        _x = doublematrix(0, _n1-1, 0, _k1-1);
        
        _saveValues = varContainers->Values;
    }
    
    memset( *_stiff, 0.0, (1*_n1)*sizeof(double) );
    memset( *_mass, 0.0, (1*_n1)*sizeof(double) );
    memset( *_damp, 0.0, (1*_n1)*sizeof(double) );
    memset( *_x, 0.0, (_n1*_k1)*sizeof(double) );
    
    hasDamping = (matContainers->DampValues != NULL) ? YES : NO;
    hasMass = (matContainers->MassValues != NULL) ? YES : NO;
    
    int rows = 1;
    for (i=0; i<solution.matrix.numberOfRows; i++) {
        n = 0;
        for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
            if (hasMass == YES) _mass[0][n] = matContainers->MassValues[j];
            if (hasDamping == YES) _damp[0][n] = matContainers->DampValues[j];
            _stiff[0][n] = matContainers->Values[j];
            for (k=0; k<varContainers->size2PrevValues; k++) {
                _x[n][k] = varContainers->PrevValues[matContainers->Cols[j]][k];
            }
            n++;
        }
        force[0] = matContainers->RHS[i];
        matContainers->Force[i][0] = force[0];
        
        double x1[n];
        double x2[n];
        double x3[n];
        for (j=0; j<n; j++) {
            x1[j] = _x[j][2];
            x2[j] = _x[j][3];
            x3[j] = _x[j][4];
        }
        [timeIntegration bossakSecondOrder:solution numberOfNodes:n dt:solution.dt massMatrix:_mass dampMatrix:_damp stiffMatrix:_stiff force:force prevSolution1:x1 prevSolution2:x2 prevSolution3:x3 alpha:solution.alpha rows:&rows cols:&_n1];
        
        n = 0;
        for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
            matContainers->Values[j] = _stiff[0][n];
            n++;
        }
        matContainers->RHS[i] = force[0];
    }
}

#pragma mark Update equations

/********************************************************************************************************************
    Add element local matrices and vectors to global matrices and vectors
 
    Arguments:
 
        FEMSolution *solution       -> Class solution which contains the global matrix and global RHS vector
        double **localStiffMatrix   -> Local matrix to be added tp the global matrix
        double *forceVector         -> The global RHS vector
        double  *loalForce          -> Element local force vector
        int n                       -> Number of nodes
        int dofs                    -> Number of dofs
        int *nodeIndexes            -> Element node to global node numbering mapping
********************************************************************************************************************/
-(void)updateGlobalEquationsModel:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution element:(Element_t * __nonnull)element localStiffMatrix:(double * __nonnull * __nonnull)localStiffMatrix forceVector:(double * __nonnull)forceVector localForce:(double * __nonnull)localForce size:(int)n dofs:(int)dofs nodeIndexes:(int * __nonnull)nodeIndexes rows:(int * __nonnull)rows cols:(int * __nonnull)cols rotateNT:(BOOL * __nullable)rotateNT crsMatrix:(FEMMatrixCRS * __nonnull)crsMatrix bandMatrix:(FEMMatrixBand * __nonnull)bandMatrix {
    
    int i, j, k, dim;
    int indexes[n];
    BOOL rotate;
    
    static BOOL (*checkPassiveElementIMP)(id, SEL, Element_t *, FEMModel *, FEMSolution *) = nil;
    static void (*rotateMatrixIMP)(id, SEL, double **, FEMSolution*, double*, int, int, int, int*) = nil;
    static void (*crsGlueLocalMatrixInGlobalIMP)(id, SEL, double **, FEMSolution *, int, int, int*) = nil;
    static void (*bandGlueLocalMatrixInGlobalIMP)(id, SEL, double **, FEMSolution *, int, int, int*) = nil;
    static dispatch_once_t onceToken;
    dispatch_once(&onceToken, ^{
        if (!checkPassiveElementIMP) {
            checkPassiveElementIMP = (BOOL (*)(id, SEL, Element_t *, FEMModel *, FEMSolution *))
            [self methodForSelector: @selector(FEMCore_checkPassiveElement:model:solution:)];
        }
        if (!rotateMatrixIMP) {
            rotateMatrixIMP = (void (*)(id, SEL, double **, FEMSolution*, double*, int, int, int, int*))
            [self methodForSelector: @selector(FEMCore_rotateMatrix:solution:vector:size:dimension:dofs:nodeIndexes:)];
        }
        if (!crsGlueLocalMatrixInGlobalIMP) {
            crsGlueLocalMatrixInGlobalIMP = (void (*)(id, SEL, double **, FEMSolution *, int, int, int*))
            [crsMatrix methodForSelector: @selector(glueLocalMatrix:inGlobal:numberOfNodes:dofs:indexes:)];
        }
        if (!bandGlueLocalMatrixInGlobalIMP) {
            bandGlueLocalMatrixInGlobalIMP = (void (*)(id, SEL, double **, FEMSolution *, int, int, int*))
            [bandMatrix methodForSelector: @selector(glueLocalMatrix:inGlobal:numberOfNodes:dofs:indexes:)];
        }
    });
    
    // Check if this element has been defined as passive
    if (checkPassiveElementIMP(self, @selector(FEMCore_checkPassiveElement:model:solution:), element, model, solution) == YES) return;
    
    rotate = YES;
    if (rotateNT != NULL) {
        rotate = *rotateNT;
    }
    
    if (rotate == YES && _normalTangentialNumberOfNodes > 0) {
        dim = model.dimension;
        memset( indexes, -1, sizeof(indexes) );
        for (i=0; i<element->Type.NumberOfNodes; i++) {
            indexes[i] = _boundaryReorder[element->NodeIndexes[i]];
        }
        rotateMatrixIMP(self, @selector(FEMCore_rotateMatrix:solution:vector:size:dimension:dofs:nodeIndexes:), localStiffMatrix, solution, localForce, n, dim, dofs, indexes);
    }
    
    switch (solution.matrix.format) {
        case MATRIX_CRS:
            crsGlueLocalMatrixInGlobalIMP(crsMatrix, @selector(glueLocalMatrix:inGlobal:numberOfNodes:dofs:indexes:), localStiffMatrix, solution, n, dofs, nodeIndexes);
            break;
            
        case MATRIX_LIST:
            // TODO: implement for the List matrix case
            break;
            
        case MATRIX_BAND:
        case MATRIX_SBAND:
            bandGlueLocalMatrixInGlobalIMP(bandMatrix, @selector(glueLocalMatrix:inGlobal:numberOfNodes:dofs:indexes:), localStiffMatrix, solution, n, dofs, nodeIndexes);
            break;
    }
    
    for (i=0; i<n; i++) {
        if (nodeIndexes[i] >= 0) {
            for (j=0; j<dofs; j++) {
                k = dofs * nodeIndexes[i] + j;
                forceVector[k] = forceVector[k] + localForce[dofs*i+j];
            }
        }
    }
}

-(void)defaultUpdateEquations:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution forElement:(Element_t * __nonnull)element realStiff:(double * __nonnull * __nonnull)stiff realForce:(double * __nonnull)force stiffRows:(int * __nonnull)rows stiffCols:(int * __nonnull)cols crsMatrix:(FEMMatrixCRS * __nonnull)crsMatrix bandMatrix:(FEMMatrixBand * __nonnull)bandMatrix {
    
    int i, n;
    matrixArraysContainer *matContainers = NULL;
    variableArraysContainer *varContainers = NULL;
    
    // TODO: add support for parallel runs
    
    static void (*updateGlobalEquationIMP)(id, SEL, FEMModel*, FEMSolution*, Element_t *, double**, double*, double*, int, int, int *, int*, int*, BOOL*, FEMMatrixCRS*, FEMMatrixBand*) = nil;
    static int (*getElementDofsIMP)(id, SEL, FEMSolution*, FEMModel*, Element_t*, int*, BOOL*) = nil;
    static dispatch_once_t onceToken;
    dispatch_once(&onceToken, ^{
        if (!updateGlobalEquationIMP) {
            updateGlobalEquationIMP = (void (*)(id, SEL, FEMModel*, FEMSolution*, Element_t *, double**, double*, double*, int, int, int *, int*, int*, BOOL*, FEMMatrixCRS*, FEMMatrixBand*))
            [self methodForSelector: @selector(updateGlobalEquationsModel:inSolution:element:localStiffMatrix:forceVector:localForce:size:dofs:nodeIndexes:rows:cols:rotateNT:crsMatrix:bandMatrix:)];
        }
        if (!getElementDofsIMP) {
            getElementDofsIMP = (int (*)(id, SEL, FEMSolution*, FEMModel*, Element_t*, int*, BOOL*))
            [self methodForSelector: @selector(getElementDofsSolution:model:forElement:atIndexes:disableDiscontinuousGalerkin:)];
        }
    });
    
    memset( _indexStore, -1, _sizeIndexStore*sizeof(int) );
    n = getElementDofsIMP(self, @selector(getElementDofsSolution:model:forElement:atIndexes:disableDiscontinuousGalerkin:), solution, model, element, _indexStore, NULL);
    
    varContainers = solution.variable.getContainers;
    int perm[n];
    for (i=0; i<n; i++) {
        perm[i] = varContainers->Perm[_indexStore[i]];
    }
    
    matContainers = solution.matrix.getContainers;
    
    updateGlobalEquationIMP(self, @selector(updateGlobalEquationsModel:inSolution:element:localStiffMatrix:forceVector:localForce:size:dofs:nodeIndexes:rows:cols:rotateNT:crsMatrix:bandMatrix:),
                            model, solution, element, stiff, matContainers->RHS, force, n, solution.variable.dofs, perm, rows, cols, NULL, crsMatrix, bandMatrix);
}

-(void)defaultUpdateEquations:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution forElement:(Element_t * __nonnull)element complexStiff:(double complex * __nonnull * __nonnull)cstiff complexForce:(double complex * __nonnull )cforce stiffRows:(int * __nonnull)rows stiffCols:(int * __nonnull)cols crsMatrix:(FEMMatrixCRS * __nonnull)crsMatrix bandMatrix:(FEMMatrixBand * __nonnull)bandMatrix {
    
    int i, j, n, dofs;
    int *perm;
    double **stiff, *force;
    matrixArraysContainer *matContainers = NULL;
    variableArraysContainer *varContainers = NULL;
    
    memset( _indexStore, -1, _sizeIndexStore*sizeof(int) );
    n = [self getElementDofsSolution:solution model:model forElement:element atIndexes:_indexStore disableDiscontinuousGalerkin:NULL];
    
    varContainers = solution.variable.getContainers;
    perm = intvec(0, n-1);
    for (i=0; i<n; i++) {
        perm[i] = varContainers->Perm[_indexStore[i]];
    }
    
    // TODO: add support for parallel runs
    
    dofs = solution.variable.dofs;
    stiff = doublematrix(0, (n*dofs)-1, 0, (n*dofs)-1);
    force = doublevec(0, (n*dofs)-1);
    
    for (i=0; i<n*dofs/2; i++) {
        
        force[2*i] = creal(cforce[i]);
        force[2*i+1] = cimag(cforce[i]);
        
        for (j=0; j<n*dofs; j++) {
            stiff[2*i][2*j]     = creal(cstiff[i][j]);
            stiff[2*i][2*j+1]   = -cimag(cstiff[i][j]);
            stiff[2*i+1][2*j]   = cimag(cstiff[i][j]);
            stiff[2*i+1][2*j+1] = creal(cstiff[i][j]);
        }
    }

     matContainers = solution.matrix.getContainers;
    
    [self updateGlobalEquationsModel:model inSolution:solution element:element localStiffMatrix:stiff forceVector:matContainers->RHS localForce:force size:n dofs:solution.variable.dofs nodeIndexes:perm rows:rows cols:cols rotateNT:NULL crsMatrix:crsMatrix bandMatrix:bandMatrix];
}

#pragma mark Finish assembly
/*******************************************************************************
    Finished the bulk assembly of the matrix equation. Optionally save the
    matrix for later use.
*******************************************************************************/
-(void)defaultFinishBulkAssemblySolution:(FEMSolution * __nonnull)solution bulkUpdate:(BOOL * __nullable)bulkUpdate {
    
    int n;
    BOOL bUpd = NO;
    matrixArraysContainer *matContainers = NULL;
    
    if (bulkUpdate != NULL) {
        bUpd = *bulkUpdate;
    } else {
        bUpd = (bUpd == YES || [(solution.solutionInfo)[@"calculate loads"] boolValue] == YES) ? YES : NO;
        bUpd = (bUpd == YES || [(solution.solutionInfo)[@"constant bulk system"] boolValue] == YES) ? YES : NO;
        bUpd = (bUpd == YES || [(solution.solutionInfo)[@"save bulk system"] boolValue] == YES) ? YES : NO;
        bUpd = (bUpd == YES || [(solution.solutionInfo)[@"constant bulk matrix"] boolValue] == YES) ? YES : NO;
    }
    
    if (bUpd == YES) {
        NSString *str = (solution.solutionInfo)[@"equation"];
        NSLog(@"FEMCore:defaultFinishBulkAssemblySolution: saving bulk values for: %@.\n", str);
        
        matContainers = solution.matrix.getContainers;
        
        n = matContainers->sizeRHS;
        if (matContainers->BulkRHS != NULL) {
            if (matContainers->sizeBulkRHS != n) {
                free_dvector(matContainers->BulkRHS, 0, matContainers->sizeBulkRHS-1);
                matContainers->BulkRHS = NULL;
            }
        }
        if (matContainers->BulkRHS == NULL) {
            matContainers->BulkRHS = doublevec(0, n-1);
            matContainers->sizeBulkRHS = n;
        }
        for (int i=0; i<n; i++) {
            matContainers->BulkRHS[i] = matContainers->RHS[i];
        }
        
        n = matContainers->sizeValues;
        if (matContainers->BulkValues != NULL) {
            if (matContainers->sizeBulkValues != n) {
                free_dvector(matContainers->BulkValues, 0, matContainers->sizeBulkValues-1);
                matContainers->BulkValues = NULL;
            }
        }
        if (matContainers->BulkValues == NULL) {
            matContainers->BulkValues = doublevec(0, n-1);
            matContainers->sizeBulkValues = n;
        }
        for (int i=0; i<n; i++) {
            matContainers->BulkValues[i] = matContainers->Values[i];
        }
    }
    
    // TODO: Add support for linear system multiply. Elmer 6100.
    // Also present in defaultFinishAssemblySolution
}

-(void)defaultFinishAssemblySolution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model timeIntegration:(FEMTimeIntegration * __nonnull)timeIntegration utilities:(FEMUtilities * __nonnull)utilities {
    
    int order;
    NSString *string;
    BOOL found;
    matrixArraysContainer *matContainers = NULL;
    
    if ([(solution.solutionInfo)[@"use global mass matrix"] boolValue] == YES) {
        FEMListUtilities *listUtilities = [FEMListUtilities sharedListUtilities];
        string = [listUtilities listGetString:model inArray:model.simulation.valuesList forVariable:@"simulation type" info:&found];
        if ([string isEqualToString:@"transient"] == YES) {
            if ((solution.solutionInfo)[@"time derivative order"] != nil) {
                order = [(solution.solutionInfo)[@"time derivative order"] intValue];
            } else order = 1;
            switch (order) {
                case 1:
                    [self defaultFirstOrderTimeGlobalModel:model inSolution:solution timeIntegration:timeIntegration utilities:utilities];
                    break;
                    
                case 2:
                    [self defaultSecondOrderTimeGlobalModel:model inSolution:solution timeIntegration:timeIntegration];
                    break;
            }
        }
    }
    matContainers = solution.matrix.getContainers;
    [self FEMCore_finishAssemblyModel:model solution:solution forceVector:matContainers->RHS sizeForceVector:matContainers->sizeRHS];
}


#pragma mark Boundary conditions
//TODO: Update this method so that it's similar to the version in Elmer 5908
-(void)dirichletBoundaryConditions:(FEMModel * __nonnull)model inSolution:(FEMSolution * __nonnull)solution usingOffset:(int * __nullable)offset offDiaginalMatrix:(BOOL * __nullable)offDiaginalMatrix {
    
    int i, j, k, kk, l, n, nb, mb, dof, numEdgeDofs, n_start, u_offset;
    double *diagScaling;
    Element_t *element, *parent, *edges, *faces;
    NSMutableString *name;
    NSString *componentName;
    NSArray *bc;
    FEMValueList *list;
    matrixArraysContainer *matContainers = NULL;
    variableArraysContainer *varContainers = NULL;
    solutionArraysContainer *solContainers = NULL;
    BOOL constantValue, scaleSystem;
    
    edges = solution.mesh.getEdges;
    faces = solution.mesh.getFaces;
    
    FEMMatrixCRS *crsMatrix = [[FEMMatrixCRS alloc] init];
    FEMListUtilities *listUtil = [FEMListUtilities sharedListUtilities];
    FEMUtilities *utilities = [[FEMUtilities alloc] init];
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
    solContainers = solution.getContainers;
    
    if ((solution.solutionInfo)[@"linear system dirichlet scaling"] != nil) {
        scaleSystem = [(solution.solutionInfo)[@"linear system dirichlet scaling"] boolValue];
    } else {
        if ((solution.solutionInfo)[@"linear system scaling"] != nil) {
            scaleSystem = [(solution.solutionInfo)[@"linear system scaling"] boolValue];
        } else scaleSystem = YES;
    }
    
    if (scaleSystem == YES) {
        BOOL rhsScaling = NO;
        [self scaleLinearSystem:solution matrix:solution.matrix rhs:matContainers->RHS result:NULL diagScaling:NULL applyScaling:NULL rhsScaling:&rhsScaling];
        diagScaling = matContainers->DiagScaling;
    } else {
        diagScaling = doublevec(0, solution.matrix.numberOfRows-1);
        for (i=0; i<solution.matrix.numberOfRows; i++) {
            diagScaling[i] = 1.0;
        }
    }
    
    u_offset = 0;
    if (offset != NULL) {
        u_offset = *offset;
    }
    
    n = solution.mesh.maxElementDofs;
    
    if (_g_Ind == NULL) {
        _g_Ind = intvec(0, n-1);
        _l_Ind = intvec(0, n-1);
        _kernStiff = doublematrix(0, n-1, 0, n-1);
        _kernWork = doublevec(0, n-1);
        if (_g_Ind == NULL || _l_Ind == NULL || _kernStiff == NULL || _kernWork == NULL)
            fatal("FEMCore:dirichletBoundaryConditions", "Memory allocation error.");
        _size1kernStiff = n;
        _size2kernStiff = n;
        _sizekernWork = n;
        _size_g_Ind = n;
        _size_l_Ind = n;
    } else if (_size_g_Ind < n) {
        free_ivector(_g_Ind, 0, _size_g_Ind);
        free_ivector(_l_Ind, 0, _size_l_Ind);
        free_dmatrix(_kernStiff, 0, _size1kernStiff, 0, _size2kernStiff);
        free_dvector(_kernWork, 0, _sizekernWork);
        
        _g_Ind = intvec(0, n-1);
        _l_Ind = intvec(0, n-1);
        _kernStiff = doublematrix(0, n-1, 0, n-1);
        _kernWork = doublevec(0, n-1);
        if (_g_Ind == NULL || _l_Ind == NULL || _kernStiff == NULL || _kernWork == NULL)
            fatal("FEMCore:dirichletBoundaryConditions", "Memory allocation error.");

        _size1kernStiff = n;
        _size2kernStiff = n;
        _sizekernWork = n;
        _size_g_Ind = n;
        _size_l_Ind = n;
    }
    
    name = [NSMutableString stringWithString:[solution.variable canonicalizeName]];
    if (solution.variable.dofs > 1) {
        [self setDirichletBoundaries:model inSolution:solution variableName:name orderOfDofs:-2 permutationOffset:NULL offDiaginalMatrix:NULL];
    }
    
    constantValue = NO;
    
    // Set Dirichlet dofs for edges and faces
    for (dof=0; dof<solution.variable.dofs; dof++) {
        
        if (solution.variable.dofs > 1) {
            int component = dof + 1;
            componentName = [utilities appendNameFromString:solution.variable.name component:&component];
        } else componentName = solution.variable.name;
        
        // Clean BC face and edge dofs
        for (i=0; i<solution.mesh.numberOfBoundaryElements; i++) {
            
            element = [self getBoundaryElement:solution atIndex:i];
            if ([self isActiveBoundaryElement:element inSolution:solution model:model] == NO) continue;
            
            // Get parent element
            parent = element->BoundaryInfo->Left;
            if (parent == NULL) parent = element->BoundaryInfo->Right;
            
            if (parent == NULL) continue;
            
            bc = [self getBoundaryCondition:model forElement:element];
            if (bc == nil) continue;
            
            list = [listUtil listFindVariable:componentName inArray:bc];
            if (list == nil) continue;
            
            constantValue = ([list type] == LIST_TYPE_CONSTANT_SCALAR) ? YES : NO;
            
            // Get indexes for boudnary and values for dofs associated to them
            n = [self getNumberOfNodesForElement:element];
            
            if (parent->Pdefs != NULL) {
                [self getBoundaryIndexes:solution.mesh forBoundaryElement:element withParentElement:parent resultVector:_g_Ind sizeVector:_size_g_Ind indexSize:&numEdgeDofs];
            } else {
                continue;
            }
            
            // Contribute this boundary to global system (i.e., solve global boundary problem)
            for (k=n; k<numEdgeDofs; k++) {
                nb = varContainers->Perm[_g_Ind[k]];
                if (nb < 0) continue;
                nb = u_offset + solution.variable.dofs*nb + dof;
                if (constantValue == YES) {
                    [crsMatrix setSymmetricDirichletInGlobal:solution atIndex:nb value:0.0];
                }else {
                    [self zeroTheNumberOfRows:nb inSolutionMatrix:solution];
                    matContainers->RHS[nb] = 0.0;
                }
            }
        }
    }
    
    // Set Dirichlet dofs for edges and faces
    for (dof=0; dof<solution.variable.dofs; dof++) {
        
        if (solution.variable.dofs > 1) {
            int component = dof + 1;
            componentName = [utilities appendNameFromString:solution.variable.name component:&component];
        } else componentName = solution.variable.name;
        
        [self setNodalLoads:model inSolution:solution variableName:componentName orderOfDofs:dof];
        [self setDirichletBoundaries:model inSolution:solution variableName:(NSMutableString *)componentName orderOfDofs:dof permutationOffset:offset offDiaginalMatrix:offDiaginalMatrix];
        
        // Dirichlet BCs for face and edge dofs
        for (i=0; i<solution.mesh.numberOfBoundaryElements; i++) {
            
            element = [self getBoundaryElement:solution atIndex:i];
            if ([self isActiveBoundaryElement:element inSolution:solution model:model] == NO) continue;
            
            bc = [self getBoundaryCondition:model forElement:element];
            if (bc == nil) continue;
            
            NSString *str1 = [componentName stringByAppendingString:@" {e}"];
            NSString *str2 = [componentName stringByAppendingString:@" {f}"];
            if ( [listUtil listCheckPresentVariable:componentName inArray:bc] == NO
                && [listUtil listCheckPresentVariable:str1 inArray:bc] == NO
                && [listUtil listCheckPresentVariable:str2 inArray:bc] == NO ) continue;
            
            // Get parent element
            parent = element->BoundaryInfo->Left;
            if (parent == NULL) {
                parent = element->BoundaryInfo->Right;
            }
            if (parent == NULL) continue;
            
            if ([listUtil listCheckPresentVariable:str1 inArray:bc] == YES) {
                if (edges != NULL) {
                    switch ([self getElementFamily:element]) {
                        case 1:
                        case 2:
                            for (j=0; j<parent->Type.NumberOfEdges; j++) {
                                n = 0;
                                for (k=0; k<element->Type.NumberOfNodes; k++) {
                                    for (l=0; l<edges[parent->EdgeIndexes[j]].Type.NumberOfNodes; l++) {
                                        if (element->NodeIndexes[k] == edges[parent->EdgeIndexes[j]].NodeIndexes[l]) n++;
                                    }
                                }
                                if (n == element->Type.NumberOfNodes) break;
                            }
                            nb = parent->Type.NumberOfNodes;
                            n = edges[parent->EdgeIndexes[j]].Type.NumberOfNodes;
                            [self localBoundaryIntegral:model inSolution:solution atBoundary:bc forElement:&edges[parent->EdgeIndexes[j]] withNumberOfNodes:n andParent:parent withNumberOfNodes:nb boundaryName:str1 functionIntegral:&_kernWork[0]];
                            
                            n = [self getElementDofsSolution:solution model:model forElement:&edges[parent->EdgeIndexes[j]] atIndexes:_g_Ind disableDiscontinuousGalerkin:NULL];
                            for (k=solContainers->defDofs[parent->BodyID-1][0]*edges[parent->EdgeIndexes[j]].NDOFs; k<n; k++) {
                                nb = varContainers->Perm[_g_Ind[k]];
                                if (nb < 0) continue;
                                nb = u_offset + solution.variable.dofs * nb + dof;
                                if (solution.matrix.isSymmetric == YES) {
                                    [crsMatrix setSymmetricDirichletInGlobal:solution atIndex:nb value:_kernWork[0]/diagScaling[nb]];
                                } else {
                                    [self zeroTheNumberOfRows:nb inSolutionMatrix:solution];
                                    matContainers->RHS[nb] = _kernWork[0]/diagScaling[nb];
                                    [self setMatrixElementForSolution:solution atIndex:nb andIndex:nb value:1.0];
                                }
                            }
                            break;
                        case 3:
                        case 4:
                            for (j=0; j<parent->Type.NumberOfFaces; j++) {
                                if (element->Type.ElementCode == faces[parent->FaceIndexes[j]].Type.ElementCode) {
                                    n = 0;
                                    for (k=0; k<element->Type.NumberOfNodes; k++) {
                                        for (l=0; l<faces[parent->FaceIndexes[j]].Type.NumberOfNodes; l++) {
                                            if (element->NodeIndexes[k] == faces[parent->FaceIndexes[j]].NodeIndexes[l]) n++;
                                        }
                                    }
                                    if (n == element->Type.NumberOfNodes) break;
                                }
                            }
                            
                            for (j=0; j<faces[parent->FaceIndexes[j]].Type.NumberOfEdges; j++) {
                                nb = edges[faces[parent->FaceIndexes[j]].EdgeIndexes[j]].Type.NumberOfNodes;
                                n = parent->Type.NumberOfNodes;
                                [self localBoundaryIntegral:model inSolution:solution atBoundary:bc forElement:&edges[faces[parent->FaceIndexes[j]].EdgeIndexes[j]] withNumberOfNodes:nb andParent:parent withNumberOfNodes:nb boundaryName:str1 functionIntegral:&_kernWork[0]];
                                
                                n = [self getElementDofsSolution:solution model:model forElement:&edges[faces[parent->FaceIndexes[j]].EdgeIndexes[j]] atIndexes:_g_Ind disableDiscontinuousGalerkin:NULL];
                                for (k=solContainers->defDofs[parent->BodyID-1][0]*edges[faces[parent->FaceIndexes[j]].EdgeIndexes[j]].NDOFs; k<n; k++) {
                                    nb = varContainers->Perm[_g_Ind[k]];
                                    if (nb < 0) continue;
                                    nb = u_offset + solution.variable.dofs * nb + dof;
                                    if (solution.matrix.isSymmetric == YES) {
                                        [crsMatrix setSymmetricDirichletInGlobal:solution atIndex:nb value:_kernWork[0]/diagScaling[nb]];
                                    } else {
                                        [self zeroTheNumberOfRows:nb inSolutionMatrix:solution];
                                        matContainers->RHS[nb] = _kernWork[0]/diagScaling[nb];
                                        [self setMatrixElementForSolution:solution atIndex:nb andIndex:nb value:1.0];
                                    }
                                }
                            }
                    }
                } // end associated edges
            } else if ([listUtil listCheckPresentVariable:str2 inArray:bc] == YES) {
                n = element->Type.NumberOfNodes;
                n = [self getElementDofsSolution:solution model:model forElement:element atIndexes:_g_Ind disableDiscontinuousGalerkin:NULL];
                for (k=0; k<n; k++) {
                    nb = varContainers->Perm[_g_Ind[k]];
                    if (nb < 0) continue;
                    nb = u_offset + solution.variable.dofs*nb + dof;
                    
                    [self zeroTheNumberOfRows:nb inSolutionMatrix:solution];
                    [self setMatrixElementForSolution:solution atIndex:nb andIndex:nb value:1.0];
                    matContainers->RHS[nb] = _kernWork[0]/diagScaling[nb];
                }
            }
            
            if (parent->Pdefs == NULL) continue;
            
            list = [listUtil listFindVariable:componentName inArray:bc];
            constantValue = ([list type] == LIST_TYPE_CONSTANT_SCALAR) ? YES : NO;
            if (constantValue == YES) continue;
            
            switch (parent->Type.dimension) {
                case 2:
                    // If no edges do not try to set boundary conditions
                    // TODO: change to break instead of continue ??
                    if (edges == NULL) continue;
                    
                    // If boundary edge has not dofs, move on to next edge
                    if (element->BDOFs <= 0) continue;
                    
                    // Number of nodes for this element
                    n = element->Type.NumberOfNodes;
                    
                    // Get indexes for boundary and values dofs associated to them
                    [self getBoundaryIndexes:solution.mesh forBoundaryElement:element withParentElement:parent resultVector:_g_Ind sizeVector:_size_g_Ind indexSize:&numEdgeDofs];
                    [self localBoundaryBDOFs:model inSolution:solution atBoundary:bc forElement:element withNumberOfNodes:numEdgeDofs boundaryName:(NSMutableString *)componentName resultMatrix:_kernStiff resultVector:_kernWork];
                    
                    if (solution.matrix.isSymmetric == YES) {
                        for (l=0; l<n; l++) {
                            nb = varContainers->Perm[_g_Ind[l]];
                            if (nb < 0) continue;
                            nb = u_offset + solution.variable.dofs*nb + dof;
                            for (k=n; k<numEdgeDofs; k++) {
                                _kernWork[k] = _kernWork[k] - _kernStiff[k][l]*matContainers->RHS[nb]*diagScaling[nb];
                            }
                        }
                        
                        for (k=n; k<numEdgeDofs; k++) {
                            for (l=n; l<numEdgeDofs; l++) {
                                _kernStiff[k-n][l-n] = _kernStiff[k][l];
                            }
                            _kernWork[k-n] = _kernWork[k];
                        }
                        l = numEdgeDofs - n;
                        if (l == 1) {
                            _kernWork[0] = _kernWork[0] / _kernStiff[0][0];
                        } else {
                            [self solveLinearSystemWithMatrix:_kernStiff andVector:_kernWork size:l leadingDimension:solution.mesh.maxElementDofs];
                        }
                        for (k=n; k<numEdgeDofs; k++) {
                            nb = varContainers->Perm[_g_Ind[k]];
                            if (nb < 0) continue;
                            nb = u_offset + solution.variable.dofs*nb + dof;
                            [crsMatrix setSymmetricDirichletInGlobal:solution atIndex:nb value:_kernWork[k-n]/diagScaling[nb]];
                        }
                    } else {
                        // Contribute this boundary to global system
                        // (i.e., solve global boundary problem)
                        for (k=n; k<numEdgeDofs; k++) {
                            nb = varContainers->Perm[_g_Ind[k]];
                            if (nb < 0) continue;
                            nb = u_offset + solution.variable.dofs*nb + dof;
                            matContainers->RHS[nb] = matContainers->RHS[nb] + _kernWork[k]/diagScaling[nb];
                            for (l=0; l<numEdgeDofs; l++) {
                                mb = varContainers->Perm[_g_Ind[l]];
                                if (mb < 0) continue;
                                mb = u_offset + solution.variable.dofs*mb + dof;
                                for (kk=matContainers->Rows[nb]+dof; kk<=matContainers->Rows[nb+1]-1; k+=solution.variable.dofs) {
                                    if (matContainers->Cols[kk] == mb) {
                                        matContainers->Values[kk] = matContainers->Values[kk] + _kernStiff[k][l] * diagScaling[mb]/diagScaling[nb];
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    break;
                case 3:
                    // If no faces present do not try to set boundary conditions
                    // TODO: change to break instead of continue ??
                    if (faces == NULL) continue;
                    
                    // Parameters of element
                    n = element->Type.NumberOfNodes;
                    
                    // Get global boundary indexes and solve dofs asscociated to them
                    [self getBoundaryIndexes:solution.mesh forBoundaryElement:element withParentElement:parent resultVector:_g_Ind sizeVector:_size_g_Ind indexSize:&numEdgeDofs];
                    
                    // If boundary face has no dofs, skip to next boundary element
                    if (numEdgeDofs == n) continue;
                    
                    // Get local solution
                    [self localBoundaryBDOFs:model inSolution:solution atBoundary:bc forElement:element withNumberOfNodes:numEdgeDofs boundaryName:(NSMutableString *)componentName resultMatrix:_kernStiff resultVector:_kernWork];
                    
                    n_start = 0;
                    if (solution.matrix.isSymmetric == YES) {
                        for (l=0; l<n; l++) {
                            nb = varContainers->Perm[_g_Ind[l]];
                            if (nb < 0) continue;
                            nb = u_offset + solution.variable.dofs*nb + dof;
                            for (k=n; k<numEdgeDofs; k++) {
                                _kernWork[k] = _kernWork[k] - _kernStiff[k][l] * matContainers->RHS[nb]*diagScaling[nb];
                            }
                        }
                        n_start++;
                    }
                    
                    // Contribute this entry to global boundary problem
                    for (k=n; k<numEdgeDofs; k++) {
                        nb = varContainers->Perm[_g_Ind[k]];
                        if (nb < 0) continue;
                        nb = u_offset + solution.variable.dofs*nb + dof;
                        matContainers->RHS[nb] = matContainers->RHS[nb] + _kernWork[k]/diagScaling[nb];
                        for (l=n_start; l<numEdgeDofs; l++) {
                            mb = varContainers->Perm[_g_Ind[l]];
                            if (mb < 0) continue;
                            mb = u_offset + solution.variable.dofs*mb + dof;
                            for (kk=matContainers->Rows[nb]+dof; kk<=matContainers->Rows[nb+1]-1; k+=solution.variable.dofs) {
                                if (matContainers->Cols[kk] == mb) {
                                    matContainers->Values[kk] = matContainers->Values[kk] + _kernStiff[k][l] * diagScaling[mb]/diagScaling[nb];
                                    break;
                                }
                            }
                        }
                    }
            }
        } //end loop over boundary elements
    } //end loop over dofs
    
    if (scaleSystem == YES) {
        [self backScaleLinearSystem:solution matrix:solution.matrix rhs:matContainers->RHS result:NULL diagScaling:NULL sizeOFDiagScaling:NULL];
    } else {
        free_dvector(diagScaling, 0, solution.matrix.numberOfRows-1);
    }
    
    NSLog(@"dirichletBoundaryConditions: dirichlet boundary conditions set.\n");
}

#pragma mark Solve

/*****************************************************************************************************
    When a new field has been computed compare it to the previous one. Different convergence measures 
    may be used. Also performs relaxation if a non-unity relaxation factor is given.
*****************************************************************************************************/
-(void)computeChange:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)aModel isSteadyState:(BOOL)steadyState nsize:(int * __nullable)nsize values:(double * __nullable)values values0:(double * __nullable)values0 sizeValues0:(int * __nullable)sizeValues0 {
    
    NSString *convergenceType, *solverName;
    int i, n, n0, relaxAfter, iterNo;
    double norm, prevNorm, bNorm, change, relaxation=1.0, maxNorm, dt, tolerance, eps;
    double *x = NULL, *r, *x0 = NULL;
    BOOL skip, convergenceAbsolute, relax, relaxBefore=NO, stat, doIt;
    FEMVariable *iterV, *timeStepVar, *veloVar;
    matrixArraysContainer *matContainers = NULL;
    variableArraysContainer *varContainers = NULL, *itervContainers = NULL;
    FEMUtilities *utilities;
    
    varContainers = solution.variable.getContainers;
    
    utilities = [[FEMUtilities alloc] init];
    
    relax = NO;
    relaxBefore = NO;
    
    if (steadyState == YES) {
        skip = [(solution.solutionInfo)[@"skip compute steady state change"] boolValue];
        if (skip == YES) return;
        
        if ((solution.solutionInfo)[@"steady state convergence measure"] != nil) {
            convergenceType = (solution.solutionInfo)[@"steady state convergence measure"];
        } else {
            convergenceType = @"norm";
        }
        
        if ((solution.solutionInfo)[@"steady state convergence absolute"] != nil) {
            convergenceAbsolute = [(solution.solutionInfo)[@"steady state convergence absolute"] boolValue];
        } else if ((solution.solutionInfo)[@"use absolute norm for convergence"] != nil) {
            convergenceAbsolute = [(solution.solutionInfo)[@"use absolute norm for convergence"] boolValue];
        } else convergenceAbsolute = NO;
        
        if ((solution.solutionInfo)[@"steady state relaxation factor"] != nil) {
            relaxation = [(solution.solutionInfo)[@"steady state relaxation factor"] doubleValue];
            relax = (relaxation != 1.0) ? YES: NO;
        }
        
        iterV = [utilities getVariableFrom:solution.mesh.variables model:aModel name:@"coupled iter" onlySearch:NULL maskName:NULL info:&stat];
        itervContainers = iterV.getContainers;
        iterNo = (int)itervContainers->Values[0];
        itervContainers = NULL;
        
        if (relax == YES) {
            if ((solution.solutionInfo)[@"steady state relaxation after"] != nil) {
                relaxAfter = [(solution.solutionInfo)[@"steady state relaxation after"] intValue];
                if (relaxAfter >= iterNo ) relax = NO;
            }
        }
        
        if (relax == YES) {
            if ((solution.solutionInfo)[@"steady state relaxation before"] != nil) {
                relaxBefore = [(solution.solutionInfo)[@"steady state relaxation before"] boolValue];
            } else {
                relaxBefore = YES;
            }
        }
    } else {
        if ((solution.solutionInfo)[@"skip compute nonlinear change"] != nil) {
            skip  = [(solution.solutionInfo)[@"skip compute nonlinear change"] boolValue];
        } else {
            skip = NO;
        }
        if (skip == YES) return;
        
        if ((solution.solutionInfo)[@"nonlinear system convergence measure"] != nil) {
            convergenceType = (solution.solutionInfo)[@"nonlinear system convergence measure"];
        } else {
            convergenceType = @"norm";
        }
        
        if ((solution.solutionInfo)[@"nonlinear system convergence absolute"] != nil) {
            convergenceAbsolute = [(solution.solutionInfo)[@"nonlinear system convergence absolute"] boolValue];
        } else if ((solution.solutionInfo)[@"use absolute norm for convergence"] != nil) {
            convergenceAbsolute = [(solution.solutionInfo)[@"use absolute norm for convergence"] boolValue];
        } else convergenceAbsolute = NO;
        
        iterV = [utilities getVariableFrom:solution.mesh.variables model:aModel name:@"nonlin iter" onlySearch:NULL maskName:NULL info:&stat];
        itervContainers = iterV.getContainers;
        iterNo = (int)itervContainers->Values[0];
        solution.variable.nonLinIter = (int)itervContainers->Values[0];
        itervContainers->Values[0] = itervContainers->Values[0]+1;
        itervContainers = NULL;
        
        if ((solution.solutionInfo)[@"nonlinear system relaxation factor"] != nil) {
            relaxation = [(solution.solutionInfo)[@"nonlinear system relaxation factor"] doubleValue];
            relax = (relaxation != 1.0) ? YES: NO;
        }
        if (relax == YES) {
            if ((solution.solutionInfo)[@"nonlinear system relaxation after"] != nil) {
                relaxAfter = [(solution.solutionInfo)[@"nonlinear system relaxation after"] intValue];
                if (relaxAfter >= solution.variable.nonLinIter) relax = NO;
            }
        }
        
        if (relax == YES) {
            if ((solution.solutionInfo)[@"nonlinear system relaxation before"] != nil) {
                relaxBefore = [(solution.solutionInfo)[@"nonlinear system relaxation before"] boolValue];
            } else relaxBefore = YES;
        }
    }
    
    if (values != NULL) {
        x = values;
    } else {
        x = varContainers->Values;
    }
    
    if (x == NULL) {
        solution.variable.norm = 0.0;
        if (steadyState == YES) {
            solution.variable.steadyChange = 0.0;
        } else {
            solution.variable.nonLinChange = 0.0;
        }
        return;
    }
    
    if (nsize != NULL && values != NULL) {
        n = *nsize;
    } else {
        n = varContainers->sizeValues;
    }
    
    stat = NO;
    n0 = 0;
    if (values0 != NULL) {
        x0 = values0;
        n0 = *sizeValues0;
        stat = YES;
    } else if (steadyState == YES) {
        if (varContainers->SteadyValues != NULL) {
            x0 = varContainers->SteadyValues;
            n0 = varContainers->sizeSteadyValues;
            stat = YES;
        }
    } else {
        if (varContainers->NonLinValues != NULL) {
            x0 = varContainers->NonLinValues;
            n0 = varContainers->sizeNonLinValues;
            stat = YES;
        }
    }
    
    if (stat == YES && n0 != n) NSLog(@"FEMCore:computeChange: possible mismatch in length of vectors.\n");
    
    if (relax == YES && relaxBefore == YES) {
        for (i=0; i<n; i++) {
            x[i] = (1.0-relaxation) * x0[i] + relaxation * x[i];
        }
    }
    
    if (steadyState == YES) {
        prevNorm = solution.variable.prevNorm;
    } else {
        prevNorm = solution.variable.norm;
    }
    
    // Compute norm here
    norm = [self FEMCore_computeNormInSolution:solution model:aModel size:n values:x];
    solution.variable.norm = norm;
    
    // The norm should be bounded in order to reach convergence
    if ((solution.solutionInfo)[@"nonlinear system max norm"] != nil) {
        maxNorm = [(solution.solutionInfo)[@"nonlinear system max norm"] doubleValue];
    } else {
        maxNorm = HUGE_VAL;
    }
    
    if (isnan(norm) != 0 || norm > maxNorm) {
        NSLog(@"FEMCore:computeChange: computed norm: %lf.\n", norm);
        fatal("FEMCore:computeChange", "Norm of solution has crossed given bounds.");
    }
    
    matContainers = solution.matrix.getContainers;
    
    if ([convergenceType isEqualToString:@"residual"] == YES) {
        // ------------------------------------------------------------------------------
        // x is solution of A(x0)x = b(x0), thus residual should be real r = b(x)-A(x)x
        // Instead we use r = b(x0)-A(x0)x0 which unfortunately is one step behind
        // ------------------------------------------------------------------------------
        r = doublevec(0, n-1);
        
        [self matrixVectorMultplyInMatrix:solution.matrix multiplyVector:x0 resultVector:r];
        for (i=0; i<n; i++) {
            r[i] = r[i] - matContainers->RHS[i];
        }
        change = [self FEMCore_computeNormInSolution:solution model:aModel size:n values:r];
        if (convergenceAbsolute == NO) {
            bNorm = [self FEMCore_computeNormInSolution:solution model:aModel size:n values:matContainers->RHS];
            if (bNorm > 0.0) {
                change = change / bNorm;
            }
        }
        free_dvector(r, 0, n-1);
    }
    else if ([convergenceType isEqualToString:@"linear system residual"] == YES) {
        // ------------------------------------------------------------------------------
        // Here the true linear system redisual r = b(x)-A(x)x is computed.
        // This option is useful for some special solvers
        // ------------------------------------------------------------------------------
        
        //TODO: add support for parallel run
        
        r = doublevec(0, n-1);
        [self matrixVectorMultplyInMatrix:solution.matrix multiplyVector:x resultVector:r];
        for (i=0; i<n; i++) {
            r[i] = r[i] - matContainers->RHS[i];
        }
        change = sqrt(cblas_ddot(n, r, 1, r, 1));
        if (convergenceAbsolute == NO) {
            bNorm = sqrt(cblas_ddot(n, matContainers->RHS, 1, matContainers->RHS, 1));
            if (bNorm > 0.0) {
                change = change / bNorm;
            }
        }
        free_dvector(r, 0, n-1);
    }
    else if ([convergenceType isEqualToString:@"solution"] == YES) {
        r = doublevec(0, n-1);
        for (i=0; i<n; i++) {
            r[i] = x[i] - x0[i];
        }
        change = [self FEMCore_computeNormInSolution:solution model:aModel size:n values:r];
        if (convergenceAbsolute == NO && (norm + prevNorm) > 0.0) {
            change = change * 2.0 / (norm+prevNorm);
        }
        free_dvector(r, 0, n-1);
    }
    else if ([convergenceType isEqualToString:@"norm"] == YES) {
        change = fabs(norm-prevNorm);
        if (convergenceAbsolute == NO && (norm+prevNorm) > 0.0) {
            change = change * 2.0 / (norm+prevNorm);
        }
    }
    else {
        NSLog(@"FEMCore:computeChange: unknown convergence measure: %@.\n", convergenceType);
        NSLog(@"FEMCore:computeChange: convergence will probably fail.\n");
        change = 1.0e6;
    }
    
    // Check for convergence: 0/1
    if (steadyState == YES) {
        solution.variable.steadyChange = change;
        if ((solution.solutionInfo)[@"steady state convergence tolerance"] != nil) {
            tolerance = [(solution.solutionInfo)[@"steady state convergence tolerance"] doubleValue];
            if (change <= tolerance) {
                solution.variable.steadyConverged = 1;
            } else {
                solution.variable.steadyConverged = 0;
            }
        }
    } else {
        solution.variable.nonLinChange = change;
        if ((solution.solutionInfo)[@"nonlinear system convergence tolerance"] != nil) {
            tolerance = [(solution.solutionInfo)[@"nonlinear system convergence tolerance"] doubleValue];
            if (change <= tolerance) {
                solution.variable.nonLinConverged = 1;
            } else {
                solution.variable.nonLinConverged = 0;
            }
        }
    }
    
    if (relax == YES && relaxBefore == NO) {
        for (i=0; i<n; i++) {
            x[i] = (1.0-relaxation)*x0[i] + relaxation*x[i];
        }
        solution.variable.norm = [self FEMCore_computeNormInSolution:solution model:aModel size:n values:x ];
    }
    
    if ((solution.solutionInfo)[@"equation"] != nil) {
        solverName = (solution.solutionInfo)[@"equation"];
    } else {
        solverName = solution.variable.name;
    }
    
    if (steadyState == YES) {
        NSLog(@"FEMCore:computeChange: SS (Iter=%d) (NRM,RELC): (%e %e) :: %@.\n", iterNo, norm, change, solverName);
    } else {
        NSLog(@"FEMCore:computeChange: NS (Iter=%d) (NRM,RELC): (%e %e) :: %@.\n", iterNo, norm, change, solverName);
    }
    
#ifdef TEST
    FEMTest *test = [FEMTest sharedTest];
    test.norm = norm;
#endif
    
    // The update of exported variables may be done internally to allow some nonlinear features
    // or in steady state level to allow coupling to other solvers.
    doIt = NO;
    if (steadyState == YES) {
        doIt = [(solution.solutionInfo)[@"update exported variable"] boolValue];
    } else {
        doIt = [(solution.solutionInfo)[@"nonlinear update exported variable"] boolValue];
    }
    if (doIt == YES) [self FEMCore_updateExportedVariablesSolution:solution model:aModel]; // TODO: This method is not fully implemented
    
    // Optional a posteriori scaling for the computed fields
    // May be usefull for some floating systems where one want to impose some intergral
    // constraints without actually using them. Then first use just one Dirichlet point
    // and then fix the level a posteriori using this condition.
    doIt = NO;
    if (steadyState == YES) {
        doIt = [(solution.solutionInfo)[@"nonlinear system set average solution"] boolValue];
    } else {
        doIt = [(solution.solutionInfo)[@"linear system set average solution"] boolValue];
    }
    if (doIt == YES) {
        // TODO: add support for parallel run
        double cTarget = [(solution.solutionInfo)[@"average solution value"] doubleValue];
        NSString *str = (solution.solutionInfo)[@"average solution weight variable"];
        if (str != nil) {
            FEMVariable *weightVar = [utilities getVariableFrom:solution.mesh.variables model:aModel name:str onlySearch:NULL maskName:NULL info:&stat];
            if (weightVar == nil) {
                NSLog(@"FEMCore:computeChange: average solution weight missing: %@.\n", str);
                fatal("FEMCore:computeChange", "Saino will abort the simulation now...");
            }
            variableArraysContainer *weightVarContainers = weightVar.getContainers;
            if (n != weightVarContainers->sizeValues) {
                NSLog(@"FEMCore:computeChange: field and weight mismatch: %@.\n", str);
                fatal("FEMCore:computeChange", "Saino will abort the simulation now...");
            }
    
            double sum;
            vDSP_sveD(weightVarContainers->Values, 1, &sum, weightVarContainers->sizeValues);
            cTarget = cTarget - cblas_ddot(weightVarContainers->sizeValues, weightVarContainers->Values, 1, x, 1) / sum;
        } else {
            double sum;
            vDSP_sveD(x, 1, &sum, n);
            cTarget = cTarget - sum / n;
        }
        for (i=0; i<n; i++) {
            x[i] = x[i] + cTarget;
        }
    }

    // Only 1st order velocity computation is implemented so far...
    if (solution.timeOrder == 1) {
        doIt = NO;
        if (steadyState == YES) {
            if ((solution.solutionInfo)[@"calculate velocity"] != nil) {
                doIt = [(solution.solutionInfo)[@"calculate velocity"] boolValue];
            } else {
                doIt = NO;
            }
        } else {
            if ((solution.solutionInfo)[@"nonlinear calculate velocity"] != nil) {
                doIt = [(solution.solutionInfo)[@"nonlinear calculate velocity"] boolValue];
            } else {
                doIt = NO;
            }
        }
        
        if (doIt == YES) {
            timeStepVar = [utilities getVariableFrom:solution.mesh.variables model:aModel name:@"time step size" onlySearch:NULL maskName:NULL info:&stat];
            variableArraysContainer *containers = timeStepVar.getContainers;
            dt = containers->Values[0];
            NSString *str1 = [solution.variable.name stringByAppendingString:@" velocity"];
            veloVar = [utilities getVariableFrom:solution.mesh.variables model:aModel name:str1 onlySearch:NULL maskName:NULL info:&stat];
            containers = veloVar.getContainers;
            for (i=0; i<containers->sizeValues; i++) {
                containers->Values[i] = (x[i] - varContainers->PrevValues[i][0]) / dt;
            }
        }
    }
    
    // Calculate derivative a.k.a sensitivity
    if (steadyState == YES) {
        if ([(solution.solutionInfo)[@"calculate derivative"]  boolValue] == YES) {
            if (iterNo > 1) {
                timeStepVar = [utilities getVariableFrom:solution.mesh.variables model:aModel name:@"derivative eps" onlySearch:NULL maskName:NULL info:&stat];
                stat = NO;
                if (timeStepVar != nil) {
                    variableArraysContainer *containers = timeStepVar.getContainers;
                    eps = containers->Values[0];
                    stat = YES;
                } else {
                    eps = [(solution.solutionInfo)[@"derivative eps"] doubleValue];
                }
                if (stat == NO) {
                    NSLog(@"FEMCore:computeChange: derivative eps not given, using one.\n");
                    eps = 1.0;
                }
                
                NSString *str1 = [solution.variable.name stringByAppendingString:@" derivative"];
                veloVar = [utilities getVariableFrom:solution.mesh.variables model:aModel name:str1 onlySearch:NULL maskName:NULL info:&stat];
                if (veloVar != nil) {
                    NSLog(@"FEMCore:computeChange: computing variable: %@.\n", str1);
                    variableArraysContainer *containers = veloVar.getContainers;
                    for (i=0; i<containers->sizeValues; i++) {
                        containers->Values[i] = (x[i] - x0[i]) / eps;
                    }
                } else {
                    NSLog(@"FEMCore:computeChange: derivative variable not present.\n");
                }
            }
        }
    }
}

/******************************************************************************************************************************
    This method decides which linear system solver to call and calls it. 
    There are two main sources of iterations whithin Saino (similar to Elmer).
    1) The FEMHUTiter class that includes the most classic iterative Krylov methods.
    2) The FEMIterativeMethods class that includes some classic iterative methods and also some more recent Krylov methods.
******************************************************************************************************************************/
-(void)iterativeSolveMatrix:(FEMMatrix * __nonnull)matrix result:(double * __nonnull)x rhs:(double * __nonnull)b dimensions:(int * __nullable)ndim solution:(FEMSolution * __nonnull)solution {
    
    NSString *str;
    int i, n, iterType, pCondType=0, *ipar, wsize, ilun=0;
    double *dpar, **work = NULL;
    double ilut_tol = 0.0;
    BOOL abortNotConverged, condition, internal, precondRecompute, refactorize;
    SEL pcondSelector=0, pcondrSelector=0, mvSelector=0;
    matrixArraysContainer *matContainers = NULL;
    variableArraysContainer *varContainers = NULL;
    
    n = matrix.numberOfRows;
    if (ndim != NULL) n = *ndim;
    
    ipar = intvec(0, 49);
    dpar = doublevec(0, 49);
    
    varContainers = solution.variable.getContainers;
    
    memset( ipar, 0, 50*sizeof(int) );
    memset( dpar, 0.0, 50*sizeof(double) );
    
    if ((solution.solutionInfo)[@"linear system iterative method"] != nil) {
        str = (solution.solutionInfo)[@"linear system iterative method"];
        NSLog(@"FEMCore:iterativeSolveMatrix: using iterative method: %@.\n", str);
    } else {
        NSLog(@"FEMCore:iterativeSolveMatrix: linear system iterative method not found, using BI-CGstab.\n");
        str = @"bi-cgstab";
    }
    
    if ([str isEqualToString:@"bi-cgstab2"] == YES)
    {
        iterType = ITER_BICGSTAB2;
    }
    else if ([str isEqualToString:@"bi-cgstab(l)"] == YES)
    {
        iterType = ITER_BICGSTABL;
    }
    else if ([str isEqualToString:@"bi-cgstab"] == YES)
    {
        iterType = ITER_BICGSTAB;
    }
    else if ([str isEqualToString:@"tfqmr"] == YES)
    {
        iterType = ITER_TFQMR;
    }
    else if ([str isEqualToString:@"cgs"] == YES)
    {
        iterType = ITER_CGS;
    }
    else if ([str isEqualToString:@"cg"] == YES)
    {
        iterType = ITER_CG;
    }
    else if ([str isEqualToString:@"gmres"] == YES)
    {
        iterType = ITER_GMRES;
    }
    else if ([str isEqualToString:@"sgs"] == YES)
    {
        iterType = ITER_SGS;
    }
    else if ([str isEqualToString:@"jacobi"] == YES)
    {
        iterType = ITER_JACOBI;
    }
    else if ([str isEqualToString:@"gcr"] == YES)
    {
        iterType = ITER_GCR;
    }
    else if ([str isEqualToString:@"richardson"] == YES)
    {
        iterType = ITER_RICHARDSON;
    }
    else {
        iterType = ITER_BICGSTAB;
    }
    
    ipar[3] = 0;
    internal = NO;
        // Solvers from the FEMHUTiter class
        // ------------------------------------------
    if (iterType == ITER_BICGSTAB) {
        ipar[3] = 8;
    }
    else if (iterType == ITER_BICGSTAB2) {
        ipar[3] = 8;
    }
    else if (iterType == ITER_TFQMR) {
        ipar[3] = 10;
    }
    else if (iterType == ITER_CG) {
        ipar[3] = 4;
    }
    else if (iterType == ITER_CGS) {
        ipar[3] = 7;
    }
    else if (iterType == ITER_GMRES) {
        if ((solution.solutionInfo)[@"linear system gmres restart"] != nil) {
            
            ipar[14] = [(solution.solutionInfo)[@"linear system gmres restart"] intValue];
        } else {
            ipar[14] = 10;
        }
        ipar[3] = 7 + ipar[14];
    }
    // Solvers from the FEMIterativeMethods class
    // ------------------------------------------
    else if (iterType == ITER_SGS) {
        ipar[3] = 1;
        if ((solution.solutionInfo)[@"sgs over relaxation factor"] != nil) {
            dpar[2] = [(solution.solutionInfo)[@"sgs over relaxation factor"] doubleValue];
            if (dpar[2] < 0.0) dpar[2] = 0.0;
            if (dpar[2] > 2.0) dpar[2] = 2.0;
        } else {
            dpar[2] = 1.8;
        }
        internal = YES;
    }
    else if (iterType == ITER_JACOBI || iterType == ITER_RICHARDSON) {
        ipar[3] = 1;
        internal = YES;
    }
    else if (iterType == ITER_GCR) {
        ipar[3] = 1;
        if ((solution.solutionInfo)[@"linear system gcr restart"] != nil) {
            ipar[16] = [(solution.solutionInfo)[@"linear system gcr restart"] intValue];
            if (ipar[16] < 1) {
                NSLog(@"FEMCore:iterativeSolveMatrix: linear system gcr restart is < 1 which is invalid. The value 1 will be used.\n");
                ipar[16] = 1;
            }
        } else {
            if ((solution.solutionInfo)[@"linear system maximum iterations"] != nil) {
                ipar[16] = [(solution.solutionInfo)[@"linear system maximum iterations"] intValue];
                if (ipar[16] < 1) {
                    NSLog(@"FEMCore:iterativeSolveMatrix: linear system gcr restart is < 1 which is invalid. The value 1 will be used.\n");
                    ipar[16] = 1;
                }
            }
        }
        internal = YES;
    }
    else if (iterType == ITER_BICGSTABL) {
        ipar[3] = 1;
        if ((solution.solutionInfo)[@"bi-cgstab(l) polynomial degree"] != nil) {
            ipar[15] = [(solution.solutionInfo)[@"bi-cgstab(l) polynomial degree"] intValue];
            if (ipar[15] < 2) ipar[15] = 2;
        } else {
            ipar[15] = 2;
        }
        internal = YES;
    }
    
    wsize = ipar[3];
    
    ipar[11] = 1;
    ipar[2] = n;
    
    if ((solution.solutionInfo)[@"linear system residual output"] != nil) {
        ipar[4] = [(solution.solutionInfo)[@"linear system residual output"] intValue];
    } else {
        ipar[4] = 1;
    }
    
    // TODO: Add support for parallel run
    
    if ((solution.solutionInfo)[@"linear system maximum iterations"] != nil) {
        ipar[9] = [(solution.solutionInfo)[@"linear system maximum iterations"] intValue];
        if (ipar[9] < 1) ipar[9] = 1;
    } else {
        ipar[9] = 1;
    }
    
    // Turned out that in our implementation, the work array is only used
    // by the gmres iterative method. So only allocate it when needed.
    if (iterType == ITER_GMRES) {
        work = doublematrix(0, n-1, 0, wsize-1);
        if (work == NULL) {
            fatal("FEMCore:iterativeSolveMatrix", "Memory allocation error.");
        }
    }
    
    if (all(x, '=', 0.0, varContainers->sizeValues) == true) {
        for (i=0; i<varContainers->sizeValues; i++) {
            x[i] = 1.0e-8;
        }
    }
    ipar[13] = 1;
    
    if ((solution.solutionInfo)[@"linear system convergence tolerance"] != nil) {
        dpar[0] = [(solution.solutionInfo)[@"linear system convergence tolerance"] doubleValue];
    }
    
    if ((solution.solutionInfo)[@"linear system divergence tolerance"] != nil) {
        dpar[1] = [(solution.solutionInfo)[@"linear system divergence tolerance"] doubleValue];
    } else {
        dpar[1] = HUGE_VAL;
    }
    
    if ((solution.solutionInfo)[@"linear system preconditioning"] != nil) {
        str = (solution.solutionInfo)[@"linear system preconditioning"];
    } else {
        str = @"none";
    }
    
    if ((solution.solutionInfo)[@"linear system symmetric ilu"] != nil) {
        matrix.cholesky = [(solution.solutionInfo)[@"linear system symmetric ilu"] boolValue];
    }
    
    if ([str isEqualToString:@"none"] == YES) {
        pCondType = PRECOND_NONE;
    }
    else if ([str isEqualToString:@"diagonal"] == YES) {
        pCondType = PRECOND_DIAGONAL;
    }
    else if ([str isEqualToString:@"ilut"] == YES) {
        if ((solution.solutionInfo)[@"linear system ilut tolerance"] != nil) {
            ilut_tol = [(solution.solutionInfo)[@"linear system ilut tolerance"] doubleValue];
        } else {
            ilut_tol = 0.0;
        }
        pCondType = PRECOND_ILUT;
    }
    else if ([[str substringToIndex:3] isEqualToString:@"ilu"] == YES) {
        if ((solution.solutionInfo)[@"linear system ilu order"] != nil) {
            ilun = [(solution.solutionInfo)[@"linear system ilu order"] intValue];
        } else ilun = [str characterAtIndex:3] - '0';
        if (ilun < 0 || ilun > 9 ) ilun = 0;        
        pCondType = PRECOND_ILUN;
    }
    else if ([[str substringToIndex:4] isEqualToString:@"bilu"] == YES) {
        ilun = [str characterAtIndex:4] - '0';
        if (ilun < 0 || ilun > 9 ) ilun = 0;
        if (solution.variable.dofs == 1) {
            NSLog(@"FEMCore:iterativeSolveMatrix: BILU for one dofs is equal to ILU.\n");
            pCondType = PRECOND_ILUN;
        } else {
            pCondType = PRECOND_BILUN;
        }
    }
    else if ([str isEqualToString:@"multigrid"] == YES) {
        pCondType = PRECOND_MG;
    }
    else if ([str isEqualToString:@"vanka"] == YES) {
        pCondType = PRECOND_VANKA;
    }
    else {
        pCondType = PRECOND_NONE;
        NSLog(@"FEMCore:iterativeSolveMatrix: unknown preconditioner type, feature disabled.\n");
    }
    
    precondRecompute = NO;
    if ((solution.solutionInfo)[@"no precondition recompute"] != nil) {
        precondRecompute = [(solution.solutionInfo)[@"no precondition recompute"] boolValue];
    }
    if (!precondRecompute) {
        n = 0;
        if ((solution.solutionInfo)[@"linear system precondition recompute"] != nil) {
            n = [(solution.solutionInfo)[@"linear system precondition recompute"] intValue];
        }
        if (n <= 0) n = 1;
        
        if ((solution.solutionInfo)[@"linear system refactorize"] != nil) {
            refactorize = [(solution.solutionInfo)[@"linear system refactorize"] boolValue];
        } else refactorize = YES;
        
        matContainers = matrix.getContainers;
        if (matContainers->ILUValues == NULL || (refactorize == YES && matrix.solveCount % n == 0)) {
            if (matrix.format == MATRIX_CRS) {
                if (matrix.complexMatrix == YES) {
                    if (pCondType == PRECOND_ILUN) {
                        FEMPrecondition *preconditioning = [[FEMPrecondition alloc] init];
                        condition = [preconditioning CRSComplexIncompleteLUMatrix:matrix fillsOrder:ilun];
                    } else if (pCondType == PRECOND_ILUT) {
                        FEMPrecondition *preconditioning = [[FEMPrecondition alloc] init];
                        condition = [preconditioning CRSComplexIlutMatrix:matrix dropTolerance:ilut_tol];
                    }
                } else {
                    FEMPrecondition *preconditioning = [[FEMPrecondition alloc] init];
                    int blocks;
                    switch (pCondType) {
                        case PRECOND_ILUN:
                            condition = [preconditioning CRSIncompleteLUMatrix:matrix fillsOrder:ilun];
                            break;
                        case PRECOND_ILUT:
                            condition = [preconditioning CRSIlutMatrix:matrix dropTolerance:ilut_tol];
                            break;
                        case PRECOND_BILUN:
                            blocks = solution.variable.dofs;
                            if (blocks <= 1) {
                                condition = [preconditioning CRSIncompleteLUMatrix:matrix fillsOrder:ilun];
                            } else {
                                if (matContainers->ILUValues == NULL) {
                                    FEMMatrix *adiag = [[FEMMatrix alloc] init];
                                    [preconditioning CRSBlockDiagonalMatrix:matrix blockDiagMatrix:adiag numberOfBlocks:blocks];
                                    condition = [preconditioning CRSIncompleteLUMatrix:adiag fillsOrder:ilun];
                                    matrixArraysContainer *aDiagContainers = adiag.getContainers;
                                    matContainers->ILURows = aDiagContainers->ILURows;
                                    matContainers->ILUCols = aDiagContainers->ILUCols;
                                    matContainers->ILUValues = aDiagContainers->ILUValues;
                                    matContainers->ILUDiag = aDiagContainers->ILUDiag;
                                    if (ilun > 0) {
                                        free_ivector(aDiagContainers->Rows, 0, aDiagContainers->sizeRows-1);
                                        free_ivector(aDiagContainers->Cols, 0, aDiagContainers->sizeCols-1);
                                        free_ivector(aDiagContainers->Diag, 0, aDiagContainers->sizeDiag-1);
                                        free_dvector(aDiagContainers->Values, 0, aDiagContainers->sizeValues-1);
                                    }
                                } else {
                                    condition = [preconditioning CRSIncompleteLUMatrix:matrix fillsOrder:ilun];
                                }
                            }
                            break;
                    case PRECOND_VANKA:
                            // No ops
                            break;
                    }
                }
            } else {
                if (pCondType == PRECOND_ILUN) {
                    NSLog(@"FEMCore:iterativeSolveMatrix: no ILU preconditioner for band matrix format.\n");
                    NSLog(@"FEMCore:iterativeSolveMatrix: using Diagonal preconditioner instead...\n");
                    pCondType = PRECOND_DIAGONAL;
                }
            }
        }
    }
    
    matrix.solveCount++;
    
    if ((solution.solutionInfo)[@"linear system abort not converged"] != nil) {
        abortNotConverged = [(solution.solutionInfo)[@"linear system abort not converged"] boolValue];
    } else {
        abortNotConverged = YES;
    }
    
    // Get the selector for the matrix-vector multiplication method we want to use
    if (matrix.isComplexMatrix == NO) {
        mvSelector= @selector(CRSMatrixVectorProduct:vector:result:info:);
    } else {
        mvSelector = @selector(CRSComplexMatrixVectorProduct:vector:result:info:);
    }
    
    // Get the selector for the preconditioning method we want to use
    if (pCondType == PRECOND_NONE) {
        if (matrix.isComplexMatrix == NO) {
            pcondSelector = @selector(CRSPCondDummyMatrix:afterPrecondition:rightHandSide:info:);
        } else {
            pcondSelector = @selector(CRSPCondDummyComplexMatrix:afterPrecondition:rightHandSide:info:);
        }
    }
    else if (pCondType == PRECOND_DIAGONAL) {
        if (matrix.isComplexMatrix == NO) {
            pcondSelector = @selector(CRSDiagPreconditionMatrix:afterPrecondition:rightHandSide:info:);
        } else {
            pcondSelector = @selector(CRSComplexDiagPreconditionMatrix:afterPrecondition:rightHandSide:info:);
        }
    }
    else if (pCondType == PRECOND_ILUN || pCondType == PRECOND_ILUT || pCondType == PRECOND_BILUN) {
        if (matrix.isComplexMatrix == NO) {
            pcondSelector = @selector(CRSLuPreconditionMatrix:afterPrecondition:rightHandSide:info:);
        } else {
            pcondSelector = @selector(CRSComplexLuPreconditionMatrix:afterPrecondition:rightHandSide:info:);
        }
    }
    
    if (matrix.isComplexMatrix == NO) {
        if (internal == YES) {
            if (ipar[4] == 0) ipar[4] = HUGE_VAL;
        }
    } else {
        ipar[2] = ipar[2] / 2;
        if (internal == YES) {
            if (ipar[4] == 0) ipar[4] = HUGE_VAL;
        }
    }
    
    // Everything is happening in this method...
    [self FEMCore_iterCallType:iterType solution:solution matrix:matrix result:x rhs:b ipar:ipar dpar:dpar work:work pcondlMethod:pcondSelector pcondrMethod:pcondrSelector matvecMethod:mvSelector mstopMethod:@selector(stopc:multiplyVector:righHandSide:ipar:)];
    
    if (matrix.isComplexMatrix == YES) ipar[2] = ipar[2] * 2;
    
    if (ipar[29] != 1 /*TODO: add support for parallel run*/) {
        if (ipar[29] == 3) {
            fatal("FEMCore:iterativeSolveMatrix", "System diverged over tolerance.");
        } else if (abortNotConverged == YES) {
            fatal("FEMCore:iterativeSolveMatrix", "Failed convergence tolerances.");
        } else {
            NSLog(@"FEMCore:iterativeSolveMatrix: Failed convergence tolerances.");
        }
    }
    
    if (work != NULL) free_dmatrix(work, 0, n-1, 0, wsize-1);
}

/******************************************************************************
    Solve a system. Various additional utilities are included and naturally a
    call to the linear system solver
 
    Arguments:
        FEMMatrix *matrix      ->   The coefficient matrix
        double *b              ->   The RHS vector
        double *x              ->   Previous solution on entry, new solution on exit (hopefully)
        double *norm           ->   L2 norm of the solution
        int dofs               ->   Number of degrees of freedom per node for this equation
        FEMSolution *solution  ->   Object for active solution
        FEMModel *model        ->   Object for active model

    TODO: Add support for parallel run
******************************************************************************/
-(void)solveSystemMatrix:(FEMMatrix * __nonnull)matrix rhs:(double * __nonnull)b result:(double * __nonnull)x norm:(double * __nonnull)norm dofs:(int)dofs solution:(FEMSolution * __nonnull)solution  model:(FEMModel * __nonnull)model {
    
    int n, i;
    double relaxation, beta, gamma;
    double t0=0.0, rt0=0.0, st, rst;
    BOOL constrainedSolve, needPrevSol, found;
    NSString *method;
    variableArraysContainer *varContainers = NULL;
    
    if ((solution.solutionInfo)[@"linear system timing"] != nil) {
        if ([(solution.solutionInfo)[@"linear system timing"] boolValue] == YES) {
            t0 = cputime();
            rt0 = realtime();
        }
    }
    
    n = matrix.numberOfRows;
    varContainers = solution.variable.getContainers;
    
    // The allocation of previous values has to be here in order to work properly
    // with the Dirichlet elimination
    needPrevSol = NO;
    if ((solution.solutionInfo)[@"nonlinear system relaxation factor"] != nil) {
        relaxation = [(solution.solutionInfo)[@"nonlinear system relaxation factor"] doubleValue];
        needPrevSol = (relaxation != 1.0) ? YES: NO;
    }
    
    if (needPrevSol == NO) {
        if ((solution.solutionInfo)[@"nonlinear system convergence measure"] != nil) {
            method = [NSString stringWithString:(solution.solutionInfo)[@"nonlinear system convergence measure"]];
            needPrevSol = ([method isEqualToString:@"residual"] == YES || [method isEqualToString:@"solution"] == YES) ? YES: NO;
        }
    }
    
    if (needPrevSol == YES) {
        found = (varContainers->NonLinValues != NULL) ? YES : NO;
        if (found == YES) {
            if (varContainers->sizeNonLinValues != n) {
                free_dvector(varContainers->NonLinValues, 0, varContainers->sizeNonLinValues-1);
                varContainers->NonLinValues = NULL;
                found = NO;
            }
        }
        if (found == NO) {
            varContainers->NonLinValues = doublevec(0, n-1);
            if (varContainers->NonLinValues == NULL) fatal("FEMCore:solveSystemMatrix", "Memory allocation error.");
            varContainers->sizeNonLinValues = n;
        }
        for (i=0; i<n; i++) {
            varContainers->NonLinValues[i] = x[i];
        }
    }
    
    
    // TODO: Add support for parallel run
    
    constrainedSolve = NO;
    if (matrix.constraint != nil) {
        if (matrix.constraint.numberOfRows >= 1) constrainedSolve = YES;
    }
    
    if (constrainedSolve == YES) {
         // TODO: Add support for constrained matrix solution
    } else {
        [self FEMCore_solveLinearSystemMatrix:matrix rhs:b result:x norm:norm dofs:dofs solution:solution model:model bulkMatrix:nil];
    }
    
    if (solution.timeOrder == 2) {
        if (varContainers->PrevValues != NULL) {
            gamma = 0.5 - solution.alpha;
            beta = ( pow((1.0 - solution.alpha), 2.0) ) / 4.0;
            for (i=0; i<n; i++) {
                varContainers->PrevValues[i][1] = (1.0/(beta*pow(solution.dt, 2.0))) * (x[i]-varContainers->PrevValues[i][2]) -
                (1.0/(beta*solution.dt))*varContainers->PrevValues[i][3] + (1.0-1.0/(2.0*beta))*varContainers->PrevValues[i][4];
                
                varContainers->PrevValues[i][0] = varContainers->PrevValues[i][3] + solution.dt*((1.0-gamma)*varContainers->PrevValues[i][4] +
                        gamma*varContainers->PrevValues[i][1]);
            }
        }
    }
    
    if ((solution.solutionInfo)[@"linear system timing"] != nil) {
        if ([(solution.solutionInfo)[@"linear system timing"] boolValue] == YES) {
            st = cputime() - t0;
            rst = realtime() - rt0;
            
            FEMListUtilities *listUtilities = [FEMListUtilities sharedListUtilities];
            [listUtilities addConstRealInClassList:model.simulation theVariable:[@"res: linsys cpu time " stringByAppendingString:[solution.variable canonicalizeName]] withValue:&st orUsingBlock:nil string:nil];
            [listUtilities addConstRealInClassList:model.simulation theVariable:[@"res: linsys real time " stringByAppendingString:[solution.variable canonicalizeName]] withValue:&rst orUsingBlock:nil string:nil];
            
            NSLog(@"FEMCore:solveSystemMatrix: linear system time (CPU, REAL) for %@: %lf %lf (s).\n", [solution.variable canonicalizeName], st, rst);
            
            if ([(solution.solutionInfo)[@"linear system timing cumulative"] boolValue] == YES) {
                double ct = [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:[@"res: cum linsys cpu time " stringByAppendingString:[solution.variable canonicalizeName]] info:&found minValue:NULL maxValue:NULL];
                st = st + ct;
                ct = [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:[@"res: cum linsys real time " stringByAppendingString:[solution.variable canonicalizeName]] info:&found minValue:NULL maxValue:NULL];
                rst = rst + ct;
                [listUtilities addConstRealInClassList:model.simulation theVariable:[@"res: cum linsys cpu time " stringByAppendingString:[solution.variable canonicalizeName]] withValue:&st orUsingBlock:nil string:nil];
                [listUtilities addConstRealInClassList:model.simulation theVariable:[@"res: cum linsys real time " stringByAppendingString:[solution.variable canonicalizeName]] withValue:&rst orUsingBlock:nil string:nil];
            }
        }
    }
}

-(double)findSolution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model backRorateNT:(BOOL * __nullable)backRorateNT {
    
    double norm;
    BOOL backBot=NO;
    variableArraysContainer *varContainers = NULL;
    matrixArraysContainer *matContainers = NULL;
    
    if ((solution.solutionInfo)[@"linear system solver disabled"] != nil) {
        if ([(solution.solutionInfo)[@"linear system solver disabled"] boolValue] == YES) return 0.0;
    }
    
    //TODO: add support for dumping the system matrix
    
    if (backRorateNT != NULL) {
        if ((solution.solutionInfo)[@"back rotate n-t solution"] != nil) {
            backBot = [(solution.solutionInfo)[@"back rotate n-t solution"] boolValue];
        } else backBot = YES;
        if (backBot ^ *backRorateNT) {
            [solution.solutionInfo setValue:@(*backRorateNT) forKey:@"back rotate n-t solution"];
        }
    }
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
    [self solveSystemMatrix:solution.matrix rhs:matContainers->RHS result:varContainers->Values norm:&norm dofs:solution.variable.dofs solution:solution model:model];
    
    if (backRorateNT != NULL) {
        if (backBot ^ *backRorateNT) {
            [solution.solutionInfo setValue:@(backBot) forKey:@"back rotate n-t solution"];
        }
    }
    solution.variable.norm = norm;
    
    return norm;
}

/*********************************************************************************
    We don't use the backward error estimate e = ||Ax-b||/||A|| ||X|| + ||b||
    as stoping criterion
*********************************************************************************/
-(double)stopc:(FEMMatrix * __nonnull)matrix multiplyVector:(double * __nonnull)x righHandSide:(double * __nonnull)b ipar:(int * __nonnull)ipar {

    int n;
    double err, *res;
    double sum1, sum2, sum3, sum4;
    FEMMatrixCRS *crsMatrix;
    matrixArraysContainer *matContainers = NULL;
    
    n = ipar[2];
    res = doublevec(0, ipar[2]-1);
    
    crsMatrix = [[FEMMatrixCRS alloc] init];
    [crsMatrix matrixVectorMultiply:matrix vector:x result:res];
    
    matContainers = matrix.getContainers;
    
    vDSP_vsubD(b, 1, res, 1, res, 1, ipar[2]);
    
    vDSP_svesqD(res, 1, &sum1, n);
    vDSP_svesqD(matContainers->Values, 1, &sum2, n);
    vDSP_svesqD(x, 1, &sum3, n);
    vDSP_svesqD(b, 1, &sum4, n);
    err = sqrt(sum1) / ( sqrt(sum2) * sqrt(sum3) + sqrt(sum4) );
    
    free_dvector(res, 0, ipar[2]-1);
    
    return err;
}

/********************************************************************************************************
    Runs a solution as defined in the MDF. There are several ways on how to skeep the execution.
    There are also three different ways how the metrices may be assembled and solved: standard (single),
    coupled amd block.
*********************************************************************************************************/
-(void)activateSolution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model timeStep:(double)dt transientSimulation:(BOOL)transient {
    
    int j, timei, timestep, passiveBCId;
    double st, dtScale = 1.0;
    double rst, rt0=0.0, tCond, t0=0.0;
    BOOL found, timing, timeDerivativeActive, isPassiveBC;
    NSString *str;
    listBuffer execIntervals = { NULL, NULL, NULL, NULL, 0, 0, 0};
    FEMVariable *timeVar = nil, *iterV = nil;
    FEMListUtilities *listUtilities;
    FEMMeshUtils *meshUtils;
    FEMUtilities *utilities;
    variableArraysContainer *timeVarContainers = NULL, *iterVContainers = NULL;
    
    listUtilities = [FEMListUtilities sharedListUtilities];
    meshUtils = [[FEMMeshUtils alloc] init];
    utilities = [[FEMUtilities alloc] init];
    
    [meshUtils setCurrentMesh:solution.mesh inModel:model];
    model.solution = solution;
    
    // The solution may be skipped by giving the start time, or stop time, or
    // the active execution intervals
    if ((solution.solutionInfo)[@"start time"] != nil) {
        st = [(solution.solutionInfo)[@"start time"] doubleValue];
        timeVar = [utilities getVariableFrom:model.variables model:model name:@"time" onlySearch:NULL maskName:NULL info:&found];
        timeVarContainers = timeVar.getContainers;
        if (timeVarContainers->Values[0] < st) return;
    }
    
    if ((solution.solutionInfo)[@"stop time"] != nil) {
        st = [(solution.solutionInfo)[@"stop time"] doubleValue];
        timeVar = [utilities getVariableFrom:model.variables model:model name:@"time" onlySearch:NULL maskName:NULL info:&found];
        timeVarContainers = timeVar.getContainers;
        if (timeVarContainers->Values[0] > st) return;
    }
    
    found = [listUtilities listGetIntegerArray:model inArray:solution.valuesList forVariable:@"exec intervals" buffer:&execIntervals];
    if (found == NO) {
        found = [listUtilities listGetIntegerArray:model inArray:solution.valuesList forVariable:@"exec interval" buffer:&execIntervals];
    }
    if (found == YES) {
        timeVar = [utilities getVariableFrom:model.variables model:model name:@"time step interval" onlySearch:NULL maskName:NULL info:&found];
        timeVarContainers = timeVar.getContainers;
        timei = round(timeVarContainers->Values[0]);
        
        timeVar = [utilities getVariableFrom:model.variables model:model name:@"time step" onlySearch:NULL maskName:NULL info:&found];
        timeVarContainers = timeVar.getContainers;
        timestep = round(timeVarContainers->Values[0]);
        
        if ((timestep-1 % execIntervals.ivector[timei]) != 0) return;
        
        free_ivector(execIntervals.ivector, 0, execIntervals.m-1);
    }
    
    // If solver timing is requested start the profiling
    timing = NO;
    if ((solution.solutionInfo)[@"solution timing"] != nil) {
        timing = [(solution.solutionInfo)[@"solution timing"] boolValue];
    }
    if (timing == YES) {
        t0 = cputime();
        rt0 = realtime();
    }
    
    solution.mesh.outputActive = YES;
    timeDerivativeActive = transient;
    
    // This is to avoid resetting of certain info that could be interesting when
    // saving data, i.e., using an auxiliary solution
    if ([(solution.solutionInfo)[@"auxiliary solution"] boolValue] == NO) {
        if ((solution.solutionInfo)[@"time step scale"] != nil) {
            dtScale = [(solution.solutionInfo)[@"time step scale"] doubleValue];
        }
        solution.dt = dtScale * dt;
        
        if (transient == YES) {
            if ((solution.solutionInfo)[@"time derivative active"] != nil) {
                timeDerivativeActive = [(solution.solutionInfo)[@"time derivative active"] boolValue];
            } else {
                timeDerivativeActive = YES;
                if ((solution.solutionInfo)[@"time derivative condition"] != nil) {
                    tCond = [(solution.solutionInfo)[@"time derivative condition"] doubleValue];
                    timeDerivativeActive = (timeDerivativeActive == YES && tCond > 0) ? YES : NO;
                }
            }
        }
        
        iterV = [utilities getVariableFrom:solution.mesh.variables model:model name:@"nonlin iter" onlySearch:NULL maskName:NULL info:&found];
        iterVContainers = iterV.getContainers;
        iterVContainers->Values[0] = 1;
        
        if ((solution.solutionInfo)[@"namespace"] != nil) {
            str = (solution.solutionInfo)[@"namespace"];
            [listUtilities listSetNameSpace:str];
        }
    }
    
    // Check for passive-active boundaries
    passiveBCId = 0;
    isPassiveBC = NO;
    j = 1;
    for (FEMBoundaryCondition *boundaryCondition in model.boundaryConditions) {
        isPassiveBC = [listUtilities listGetLogical:model inArray:boundaryCondition.valuesList forVariable:@"passive target" info:&found];
        if (passiveBCId == YES) {
            passiveBCId = j;
            break;
        }
        j++;
    }
    if (isPassiveBC == YES) {
        [self getPassiveBoundaryAtIndex:passiveBCId model:model mesh:(FEMMesh *)model.mesh solution:solution];
        NSLog(@"FEMCore:activateSolution: passive element bc no. %d assigned to bc-id no. %d.\n", j, passiveBCId);
    }
    
    // Get the correct type of solution: standard (single), coupled ot block
    // This is where everything really happens!
    if (solution.solutionMode == SOLUTION_MODE_COUPLED || solution.solutionMode == SOLUTION_MODE_ASSEMBLY ) {
        // TODO:...
    } else if (solution.solutionMode == SOLUTION_MODE_BLOCK) {
        // TODO:...
    } else {
        [self FEMCore_singleSolution:solution model:model timeStep:dtScale*dt transientSimulation:timeDerivativeActive];
    }
    
    [listUtilities listSetNameSpace:nil];
    solution.dt = dt;
    
    // After solution register the timing if requested
    if (timing == YES) {
        st = cputime() - t0;
        rst = realtime() - rt0;
        NSString *str = (solution.solutionInfo)[@"equation"];
        [listUtilities addConstRealInClassList:model.simulation theVariable:[@"res: solution cpu time " stringByAppendingString:str] withValue:&st orUsingBlock:nil string:nil];
        [listUtilities addConstRealInClassList:model.simulation theVariable:[@"res: solution real time " stringByAppendingString:str] withValue:&rst orUsingBlock:nil string:nil];
        
        NSLog(@"FEMCore:solveSystemMatrix: solution time (CPU, REAL) for %@: %lf %lf (s).\n", str, st, rst);
        
        if ([(solution.solutionInfo)[@"solution timing cumulative"] boolValue] == YES) {
            double ct = [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:[@"res: cum solution cpu time " stringByAppendingString:str] info:&found minValue:NULL maxValue:NULL];
            st = st + ct;
            ct = [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:[@"res: cum solution real time " stringByAppendingString:str] info:&found minValue:NULL maxValue:NULL];
            rst = rst + ct;
            [listUtilities addConstRealInClassList:model.simulation theVariable:[@"res: cum solution cpu time " stringByAppendingString:str] withValue:&st orUsingBlock:nil string:nil];
            [listUtilities addConstRealInClassList:model.simulation theVariable:[@"res: cum solution real time " stringByAppendingString:str] withValue:&rst orUsingBlock:nil string:nil];
        }
    }
}

// Solve the equations one-by-one
-(void)solveEquationsModel:(FEMModel * __nonnull)model timeStep:(double * __nonnull)dt transientSimulation:(BOOL)transient coupledMinIteration:(int)coupledMinIter coupleMaxIteration:(int)coupleMaxIter steadyStateReached:(BOOL * __nonnull)steadyStateReached realTimeStep:(int * __nonnull)realTimeStep {
    
    int i, j, n=0, rgOrder=0;
    double *steadyIt = NULL, prevDt=0.0;
    BOOL found, scanning, *doneThis, *afterConverged, rungeKutta;
    NSString *when;
    RungeKutta_t *rgCoeff = NULL;
    FEMVariable *iterV = nil;
    FEMSolution *solution;
    FEMListUtilities *listUtilities;
    FEMUtilities *utilities;
    variableArraysContainer *variableContainers = NULL, *varContainers = NULL;
    
    listUtilities = [FEMListUtilities sharedListUtilities];
    utilities = [[FEMUtilities alloc] init];
    
    scanning = ([[listUtilities listGetString:model inArray:model.simulation.valuesList forVariable:@"simulation type" info:&found] isEqualToString:@"scanning"] == YES) ? YES : NO;
    
    if (transient) {
        for (FEMSolution *solution in model.solutions) {
            if (solution.hasBuiltInSolution != NO || solution.plugInPrincipalClassInstance != nil) [self initializeTimeStepInSolution:solution model:model];
        }
    }
    
    if (transient == YES || scanning == YES) {
        solution = model.solutions[0];
        iterV = [utilities getVariableFrom:solution.mesh.variables model:model name:@"coupled iter" onlySearch:NULL maskName:NULL info:&found];
        varContainers = iterV.getContainers;
        steadyIt = &varContainers->Values[0];
    }
    
    for (FEMSolution *solution in model.solutions) {
        if (solution.hasBuiltInSolution == NO && solution.plugInPrincipalClassInstance == nil) continue;
        if ((solution.solutionInfo)[@"invoke solution computer"] != nil) {
            when = (solution.solutionInfo)[@"invoke solution computer"];
            if ([when isEqualToString:@"before time step"] == YES) {
                [self activateSolution:solution model:model timeStep:*dt transientSimulation:transient];
                // TODO: add support for parallel run
            }
        } else {
            if (solution.solutionSolveWhen == SOLUTION_SOLVE_AHEAD_TIME) {
                [self activateSolution:solution model:model timeStep:*dt transientSimulation:transient];
                // TODO: add support for parallel run
            }
        }
    }
    
    doneThis = (BOOL*)malloc(sizeof(BOOL) * model.numberOfSolutions );
    afterConverged = (BOOL*)malloc(sizeof(BOOL) * model.numberOfSolutions );
    
    for (i=0; i<model.numberOfSolutions; i++) {
        afterConverged[i] = [(solution.solutionInfo)[@"coupled system after others converged"] boolValue];
    }
    
    if (prevDt == 0) prevDt = *dt;
    
    if (*realTimeStep > 2) {
        i = 0;
        for (FEMSolution *solution in model.solutions) {
            if (solution.hasBuiltInSolution == NO && solution.plugInPrincipalClassInstance == nil) continue;
            
            rungeKutta = NO;
            if (transient == YES && solution.timeOrder == 1) {
                rungeKutta = ([(solution.solutionInfo)[@"time stepping method"] isEqualToString:@"runge-kutta"] == YES) ? YES : NO;
            }
            
            if (rungeKutta == NO) continue;
            if (rgCoeff == NULL) rgCoeff =  allocateRungeKutta( model.numberOfSolutions);
            
            variableContainers = solution.variable.getContainers;
            n = variableContainers->sizeValues;
            rgCoeff[i].k1 = doublevec(0, n-1);
            rgCoeff[i].k2 = doublevec(0, n-1);
            rgCoeff[i].k3 = doublevec(0, n-1);
            rgCoeff[i].k4 = doublevec(0, n-1);
            
            for (j=0; j<variableContainers->size1PrevValues; j++) {
                rgCoeff[i].k1[j] = variableContainers->PrevValues[j][0] -variableContainers->PrevValues[j][1];
            }
            for (j=0; j<n; j++) {
                rgCoeff[i].k1[j] = *dt * rgCoeff[i].k1[j] / prevDt;
            }
            rgOrder = solution.order;
            switch (rgOrder) {
                case 2:
                    for (j=0; j<variableContainers->sizeValues; j++) {
                        variableContainers->Values[j] = variableContainers->PrevValues[j][0] + rgCoeff[i].k1[j];
                    }
                    break;
                case 4:
                    for (j=0; j<variableContainers->sizeValues; j++) {
                        variableContainers->Values[j] = variableContainers->PrevValues[j][0] + rgCoeff[i].k1[j]/2.0;
                    }
                    break;
            }
            i++;
        }
        if (rgOrder == 4) *dt = *dt / 2.0;
    }
    
    [self FEMCore_solveCoupledModel:model timeStep:*dt coupledMinIteration:coupledMinIter coupleMaxIteration:coupleMaxIter transientSimulation:transient scanning:scanning doneThis:doneThis afterConverged:afterConverged steadyIt:steadyIt];
    
    if (rgCoeff != NULL) {
        i = 0;
        for (FEMSolution *solution in model.solutions) {
            if (rgCoeff[i].k1 == NULL) continue;
            variableContainers = solution.variable.getContainers;
            rgOrder = solution.order;
            switch (rgOrder) {
                case 2:
                    for (j=0; j<n; j++) {
                        rgCoeff[i].k2[j] = variableContainers->Values[j] - variableContainers->PrevValues[j][0];
                    }
                    for (j=0; j<variableContainers->sizeValues; j++) {
                        variableContainers->Values[j] = variableContainers->PrevValues[j][0] + (rgCoeff[i].k1[j] + rgCoeff[i].k2[j])/2.0;
                    }
                    solution.variable.norm = [self FEMCore_computeNormInSolution:solution model:model size:n values:variableContainers->Values];
                    break;
                case 4:
                    for (j=0; j<n; j++) {
                        rgCoeff[i].k2[j] = 2.0 * (variableContainers->Values[j] - variableContainers->PrevValues[j][0]);
                    }
                    for (j=0; j<variableContainers->sizeValues; j++) {
                        variableContainers->Values[j] = variableContainers->PrevValues[j][0] + rgCoeff[i].k2[j]/2.0;
                    }
                    break;
            }
            i++;
        }
        
        if (rgOrder > 2) {
            [self FEMCore_solveCoupledModel:model timeStep:*dt coupledMinIteration:coupledMinIter coupleMaxIteration:coupleMaxIter transientSimulation:transient scanning:scanning doneThis:doneThis afterConverged:afterConverged steadyIt:steadyIt];
            
            i = 0;
            for (FEMSolution *solution in model.solutions) {
                if (rgCoeff[i].k1 == NULL) continue;
                variableContainers = solution.variable.getContainers;
                
                for (j=0; j<n; j++) {
                    rgCoeff[i].k3[j] = 2.0 * (variableContainers->Values[j] - variableContainers->PrevValues[j][0]);
                }
                for (j=0; j<variableContainers->sizeValues; j++) {
                    variableContainers->Values[j] = variableContainers->PrevValues[j][0] + rgCoeff[i].k3[j];
                }
                i++;
            }
            *dt = 2.0 * (*dt);
            [self FEMCore_solveCoupledModel:model timeStep:*dt coupledMinIteration:coupledMinIter coupleMaxIteration:coupleMaxIter transientSimulation:transient scanning:scanning doneThis:doneThis afterConverged:afterConverged steadyIt:steadyIt];
            
            i = 0;
            for (FEMSolution *solution in model.solutions) {
                if (rgCoeff[i].k1 == NULL) continue;
                variableContainers = solution.variable.getContainers;
                
                for (j=0; j<n; j++) {
                    rgCoeff[i].k4[j] = variableContainers->Values[j] - variableContainers->PrevValues[j][0];
                }
                for (j=0; j<variableContainers->sizeValues; j++) {
                    variableContainers->Values[j] = variableContainers->PrevValues[j][0]
                       + (rgCoeff[i].k1[j] + 2.0*rgCoeff[i].k2[j] + 2.0*rgCoeff[i].k3[j] + rgCoeff[i].k4[j]) / 6.0;
                }
                solution.variable.norm = [self FEMCore_computeNormInSolution:solution model:model size:n values:variableContainers->Values];
                i++;
            }
        }
        
        for (i=0; i<model.numberOfSolutions; i++) {
            if (rgCoeff[i].k1 != NULL) {
                free_dvector(rgCoeff[i].k1, 0, n-1);
                free_dvector(rgCoeff[i].k2, 0, n-1);
                free_dvector(rgCoeff[i].k3, 0, n-1);
                free_dvector(rgCoeff[i].k4, 0, n-1);
            }
        }
        free(rgCoeff);
    }
    
    prevDt = *dt;
    
    for (FEMSolution *solution in model.solutions) {
        if (solution.hasBuiltInSolution == NO && solution.plugInPrincipalClassInstance == nil) continue;
        
        if ((solution.solutionInfo)[@"invoke solution computer"] != nil) {
            when = (solution.solutionInfo)[@"invoke solution computer"];
            if ([when isEqualToString:@"after time step"] == YES) {
                [self activateSolution:solution model:model timeStep:*dt transientSimulation:transient];
                // TODO: add support for parallel run
            }
        } else {
            if (solution.solutionSolveWhen == SOLUTION_SOLVE_AFTER_TIME) {
                [self activateSolution:solution model:model timeStep:*dt transientSimulation:transient];
                // TODO: add support for parallel run
            }
        }
    }
    
    if (transient == NO) {
        BOOL all = YES;
        for (i=0; i<model.numberOfSolutions; i++) {
            if (doneThis[i] == NO) {
                all = NO;
                break;
            }
        }
        *steadyStateReached = (all == YES) ? YES : NO;
    }
    
    free(doneThis);
    free(afterConverged);
}

-(dispatch_queue_t __nonnull)getDispatchQueueAndInfoForDeviceType:(NSString * __nonnull)deviceType {
    
    dispatch_queue_t queue = NULL;
    cl_int error;
    cl_device_local_mem_type localMemType;
    cl_ulong globalMemSize, globalMemCacheSize, localMemSize, maxMemAllocSize;
    cl_uint clockFrequency, vectorWidth, maxComputeUnits;
    size_t maxWorkItemDims, maxWorkGroupSize, maxWorkItemSizes[3];
    char name[128];
    char vendor[128];
    char driver[128];
    char profile[128];
    char extensions[1024];
    
    if ([deviceType isEqualToString:@"gpu"] == YES) {
        queue = gcl_create_dispatch_queue(CL_DEVICE_TYPE_GPU, NULL);
        // If no OpenCL compatible GPU found, try to use the CPU
        if (queue == NULL) {
            queue = gcl_create_dispatch_queue(CL_DEVICE_TYPE_CPU, NULL);
        }
    } else if ([deviceType isEqualToString:@"cpu"] == YES) {
        queue = gcl_create_dispatch_queue(CL_DEVICE_TYPE_CPU, NULL);
    } else {
        fatal("FEMCore:getDispatchQueueAndInfoForDeviceType", "OpenCL device not supported. Saino will abort the simulation now...");
    }
    
    cl_uint vectorTypes[] = {CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR, CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT, CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT, CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG, CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT, CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE};
	char *vectorTypeNames[] = {"char", "short", "int", "long", "float", "double"};
    
    // Print some informations about the device
    cl_device_id device = gcl_get_device_id_with_dispatch_queue(queue);
    error = clGetDeviceInfo(device, CL_DEVICE_NAME, 128, name, NULL);
    error |= clGetDeviceInfo(device, CL_DEVICE_VENDOR, 128, vendor, NULL);
    error |= clGetDeviceInfo(device, CL_DRIVER_VERSION, 128, driver, NULL);
    error |= clGetDeviceInfo(device, CL_DEVICE_PROFILE, 128, profile, NULL);
	error |= clGetDeviceInfo(device, CL_DEVICE_EXTENSIONS, 1024, extensions, NULL);
	error |= clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_TYPE, sizeof(localMemType), &localMemType, NULL);
	
	error |= clGetDeviceInfo(device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(globalMemSize), &globalMemSize, NULL);
	error |= clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(localMemSize), &localMemSize, NULL);
	error |= clGetDeviceInfo(device, CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE, sizeof(globalMemCacheSize), &globalMemCacheSize, NULL);
	error |= clGetDeviceInfo(device, CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(maxMemAllocSize), &maxMemAllocSize, NULL);
	
	error |= clGetDeviceInfo(device, CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(clockFrequency), &clockFrequency, NULL);
    error |= clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(maxComputeUnits), &maxComputeUnits, NULL);
	
	error |= clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(maxWorkGroupSize), &maxWorkGroupSize, NULL);
	error |= clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(maxWorkItemDims), &maxWorkItemDims, NULL);
	error |= clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(maxWorkItemSizes), &maxWorkItemSizes, NULL);
    
    NSLog(@"FEMCore:getDispatchQueueAndInfoForDeviceType: Device Name: %s.\n", name);
    NSLog(@"FEMCore:getDispatchQueueAndInfoForDeviceType: Vendor: %s.\n", vendor);
    NSLog(@"FEMCore:getDispatchQueueAndInfoForDeviceType: Driver: %s.\n", driver);
	NSLog(@"FEMCore:getDispatchQueueAndInfoForDeviceType: Profile: %s.\n", profile);
	NSLog(@"FEMCore:getDispatchQueueAndInfoForDeviceType: Supported Extensions: %s.\n", extensions);
	
	NSLog(@"FEMCore:getDispatchQueueAndInfoForDeviceType: Local Mem Type (Local=1, Global=2): %i.\n", (int)localMemType);
	NSLog(@"FEMCore:getDispatchQueueAndInfoForDeviceType: Global Mem Size (MB): %i.\n", (int)globalMemSize/(1024*1024));
	NSLog(@"FEMCore:getDispatchQueueAndInfoForDeviceType: Global Mem Cache Size (Bytes): %i.\n", (int)globalMemCacheSize);
	NSLog(@"FEMCore:getDispatchQueueAndInfoForDeviceType: Local Mem Size (Bytes): %i.\n", (int)localMemSize);
	NSLog(@"FEMCore:getDispatchQueueAndInfoForDeviceType: Max Mem Alloc Size (MB): %ld.\n", (long int)maxMemAllocSize/(1024*1024));
	
	NSLog(@"FEMCore:getDispatchQueueAndInfoForDeviceType: Clock Frequency (MHz): %i.\n", clockFrequency);
    NSLog(@"FEMCore:getDispatchQueueAndInfoForDeviceType: Max Compute Units: %i.\n", maxComputeUnits);
	
	for(int i=0;i<6;i++) {
		error |= clGetDeviceInfo(device, vectorTypes[i], sizeof(vectorWidth), &vectorWidth, NULL);
		NSLog(@"FEMCore:getDispatchQueueAndInfoForDeviceType: Vector type width for: %s = %i.\n", vectorTypeNames[i], vectorWidth);
	}
	
	NSLog(@"FEMCore:getDispatchQueueAndInfoForDeviceType: Max Work Group Size: %lu.\n", maxWorkGroupSize);
	NSLog(@"FEMCore:getDispatchQueueAndInfoForDeviceType: Max Work Item Dimensions: %i.\n", (int)maxWorkItemDims);
    NSLog(@"FEMCore:getDispatchQueueAndInfoForDeviceType: Max Work Item Sizes: %i %i %i.\n", (int)maxWorkItemSizes[0], (int)maxWorkItemSizes[1], (int)maxWorkItemSizes[2]);

    return queue;
}


@end
