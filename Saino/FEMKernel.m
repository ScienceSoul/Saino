//
//  FEMKernel.m
//  Saino
//
//  Created by Seddik hakime on 15/04/11.
//  Copyright 2011 ScienceSoul. All rights reserved.
//

#import "FEMKernel.h"

#import "math.h"
#import <complex.h>

#import "FEMValueList.h"
#import "FEMBodyForce.h"
#import "FEMBoundaryCondition.h"

#import "memory.h"
#import "Utils.h"
#import "Lists.h"

static int ITER_BICGSTAB     =  320;
static int ITER_TFQMR        =  330;
static int ITER_CG           =  340;
static int ITER_CGS          =  350;
static int ITER_GMRES        =  360;
static int ITER_BICGSTAB2    =  370;
static int ITER_SGS          =  380;
static int ITER_JACOBI       =  390;
static int ITER_BICGSTABL    =  400;
static int ITER_GCR          =  410;

static int PRECOND_NONE      =  500;
static int PRECOND_DIAGONAL  =  510;
static int PRECOND_ILUN      =  520;
static int PRECOND_ILUT      =  530;
static int PRECOND_MG        =  540;
static int PRECOND_BILUN     =  550;
static int PRECOND_VANKA     =  560;

@interface FEMKernel ()

// Solving linear stsems and compute norms
-(void)iterCall:(int)iterType: (FEMSolution *)solution: (int *)ipar: (double *)dpar: (double **)work: (SEL)pcondlMethod: (SEL)pcondrMethod: (SEL)matvecMethod: (SEL)mstopMethod;
-(void)FEMKernel_iterSolver:(FEMSolution *)solution;
-(void)FEMKernel_backRotateNTSystem:(FEMSolution *)solution;
-(double)FEMKernel_computeNorm:(FEMSolution *)solution: (int)n: (double *)values;
-(void)FEMKernel_computeChange:(FEMSolution *)solution: (BOOL)steadyState: (int*)nsize: (double *)values: (double *)values0;
-(void)FEMKernel_solveLinearSystem:(FEMSolution *)solution;
-(void)FEMKernel_solveSystem:(FEMSolution *)solution;

// Check pasisve element
-(BOOL)FEMKernel_checkPassiveElement:(FEMModel *)model: (FEMSolution *)solution: (Element_t *)element;

// Rotate matrix
-(void)FEMKernel_rotateMatrix:(FEMSolution *)solution: (double **)matrix: (double *)vector: (int)n: (int)dim: (int)dofs: (int *)nodeIndexes;

// Update global force
-(void)FEMKernel_updateGlobalForce:(FEMModel *)model: (FEMSolution *)solution: (Element_t *)element: (double *)forceVector: (double *)localForce: (int)n: (int)dofs: (int *)nodeIndexes: (BOOL *)rotateNT;

// Manipulate matrix coefficients for time dependent simulations
-(void)FEMKernel_addFirstOrderTime:(FEMModel *)model: (FEMSolution *)solution: (Element_t *)element: (double **)massMatrix: (double **)stiffMatrix: (double *)force: (double)dt: (int)n: (int)dofs: (int *)nodeIndexes: (int *)rows: (int *)cols;

// Update global equations
-(void)FEMKernel_updateGlobalEquations:(FEMModel *)model: (FEMSolution *)solution: (Element_t *)element: (double **)stiffMatrix: (double *)force: (int)n: (int)dofs: (int *)nodeIndexes: (int *)rows: (int *)cols: (BOOL *)rotateNT: (BOOL *)bulkUpdate;

@end

@implementation FEMKernel

#pragma mark Private methods...

#pragma mark Solve linear systems and norms

-(void)iterCall:(int)iterType :(FEMSolution *)solution :(int *)ipar :(double *)dpar :(double **)work :(SEL)pcondlMethod :(SEL)pcondrMethod :(SEL)matvecMethod :(SEL)mstopMethod {
    
    FEMHUTIter *hutisolver;
    
    if (pcondrMethod == 0) {
        pcondrMethod = @selector(CRS_pcond_dummy::::);
    }
    
    hutisolver = [[FEMHUTIter alloc] init];
    
    if (iterType == ITER_BICGSTAB) { // Solve with BI-CGSTAB
        
        [hutisolver dbicgstabSolve:solution :ipar[2] :ipar[3] :ipar :dpar :work :pcondlMethod :pcondrMethod :matvecMethod :mstopMethod];
    }
    else if (iterType == ITER_BICGSTAB2) { // Solve with BI-CGSTAB2
        
        [hutisolver dbicgstab2Solve:solution :ipar[2] :ipar[3] :ipar :dpar :work :pcondlMethod :pcondrMethod :matvecMethod :mstopMethod];
    }
    else if (iterType == ITER_TFQMR) { // Solve with TFQMR
        
        [hutisolver dtfqmrSolve:solution :ipar[2] :ipar[3] :ipar :dpar :work :pcondlMethod :pcondrMethod :matvecMethod :mstopMethod];
    }
    else if (iterType == ITER_CG) { // Solve with CG
        
        [hutisolver dcgSolve:solution :ipar[2] :ipar[3] :ipar :dpar :work :pcondlMethod :pcondrMethod :matvecMethod :mstopMethod];
    }
    else if (iterType == ITER_CGS) { // Solve with CGS
        
        [hutisolver dcgsSolve:solution :ipar[2] :ipar[3] :ipar :dpar :work :pcondlMethod :pcondrMethod :matvecMethod :mstopMethod];
    }
    else if (iterType == ITER_GCR) { // Solve with GMRES
        
        [hutisolver dgmresSolve:solution :ipar[2] :ipar[3] :ipar :dpar :work :pcondlMethod :pcondrMethod :matvecMethod :mstopMethod];
    }
    else if (iterType == ITER_SGS) { // Solve with SGS
    
        [hutisolver dsgsSolve:solution :ipar[2] :ipar[3] :ipar :dpar :work :pcondlMethod :pcondrMethod :matvecMethod :mstopMethod];
    }
    else if (iterType == ITER_JACOBI) { // Solve with Jacobi
        
        [hutisolver djacobiSolve:solution :ipar[2] :ipar[3] :ipar :dpar :work :pcondlMethod :pcondrMethod :matvecMethod :mstopMethod];
    }
    else if (iterType == ITER_BICGSTABL) { // Solve with BI-CGSTAB(l)
        
        [hutisolver dbicgstablSolve:solution :ipar[2] :ipar[3] :ipar :dpar :work :pcondlMethod :pcondrMethod :matvecMethod :mstopMethod];
    }
    else if (iterType == ITER_GCR) { // Solve with GCR
        
        [hutisolver dgcrSolve:solution :ipar[2] :ipar[3] :ipar :dpar :work :pcondlMethod :pcondrMethod :matvecMethod :mstopMethod];
    }
    
}

-(void)FEMKernel_iterSolver:(FEMSolution *)solution {
    
    NSString *str;
    int i, n, iterType, pCondType=0, *ipar, wsize, ilun;
    double *dpar, **work;
    double ilut_tol;
    
    BOOL abortNotConverged;
    
    SEL pcondSelector=0, pcondrSelector=0, mvSelector=0;
    
    n = [solution matrixNumberOfRows];
    ipar = intvec(0, 49);
    dpar = doublevec(0, 49);
    
    for (i=0; i<50; i++) {
        ipar[i] = 0;
        dpar[i] = 0.0;
    }
    
    str = [NSString stringWithString:[solution solutionInfoForKey:@"Linear system iteration method"]];
    
    if ([str isEqualToString:@"Bi-CGSTAB2"] == YES) 
    {
        iterType = ITER_BICGSTAB2;
    } 
    else if ([str isEqualToString:@"Bi-CGSTAB(l)"] == YES) 
    {
        iterType = ITER_BICGSTABL;
    } 
    else if ([str isEqualToString:@"Bi-CGSTAB"] == YES) 
    {
        iterType = ITER_BICGSTAB;
    } 
    else if ([str isEqualToString:@"TFQMR"] == YES) 
    {
        iterType = ITER_TFQMR;
    } 
    else if ([str isEqualToString:@"CGS"] == YES) 
    {
        iterType = ITER_CGS;
    } 
    else if ([str isEqualToString:@"CG"] == YES) 
    {
        iterType = ITER_CG;
    } 
    else if ([str isEqualToString:@"GMRES"] == YES) 
    {
        iterType = ITER_GMRES;
    } 
    else if ([str isEqualToString:@"SGS"] == YES) 
    {
        iterType = ITER_SGS;
    } 
    else if ([str isEqualToString:@"JACOBI"] == YES) 
    {
        iterType = ITER_JACOBI;
    } 
    else if ([str isEqualToString:@"GCR"] == YES) 
    {
        iterType = ITER_GCR;
    } 
    else 
    {
        iterType = ITER_BICGSTAB;
    }
    
    ipar[3] = 0;
    wsize = ipar[3];
    
    if (iterType == ITER_BICGSTAB) 
    {
        ipar[3] = 8;
        wsize = ipar[3];
    }
    else if (iterType == ITER_BICGSTAB2) 
    {
        ipar[3] = 8;
        wsize = ipar[3];
    }
    else if (iterType == ITER_TFQMR)
    {
        ipar[3] = 10;
        wsize = ipar[3];
    }
    else if (iterType == ITER_CG)
    {
        ipar[3] = 4;
        wsize = ipar[3];
    }
    else if (iterType == ITER_CGS)
    {
        ipar[3] = 7;
        wsize = ipar[3];
    }
    else if (iterType == ITER_GMRES)
    {
        if ([solution solutionInfoForKey:@"Linear system GMRES restart"] != nil) {
            
            ipar[14] = [[solution solutionInfoForKey:@"Linear system GMRES restart"] intValue];
        } else {
            ipar[14] = 10;
        }
        ipar[3] = 7 + ipar[14];
        wsize = ipar[3];
    }
    else if (iterType == ITER_SGS)
    {
        ipar[3] = 1;
        wsize =ipar[3];
        if ([solution solutionInfoForKey:@"SGS over relaxation factor"] != nil) {
            dpar[2] = [[solution solutionInfoForKey:@"SGS over relaxation factor"] doubleValue];
            if (dpar[2] < 0.0 || dpar[2] > 2.0) {
                errorfunct("FEMKernel_iterSolver", "Value for SGS over relaxation factor out of bounds (min:0; max:2).");
            }
        } else {
            dpar[2] = 1.8;
        }
    }
    else if (iterType == ITER_JACOBI)
    {
        ipar[3] = 1;
        wsize = ipar[3];
    }
    else if (iterType == ITER_GCR) 
    {
        ipar[3] = 1;
        wsize = ipar[3];
        if ([solution solutionInfoForKey:@"Linear system GCR restart"] != nil) {
            ipar[16] = [[solution solutionInfoForKey:@"Linear system GCR restart"] intValue];
        } else {
            if ([solution solutionInfoForKey:@"Linear system maximum iterations"] != nil) {
                ipar[16] = [[solution solutionInfoForKey:@"Linear system maximum iterations"] intValue];
                if (ipar[16] < 1) errorfunct("FEMKernel_iterSolver", "Parameter for GCR should be equal or greater than 1.");
            } else {
                ipar[16] = 1;
            }
        }
    }
    else if (iterType == ITER_BICGSTABL)
    {
        ipar[3] = 1;
        wsize = ipar[3];
        if ([solution solutionInfoForKey:@"Bi-CGSTAB(l) polynomial degree"] != nil) {
            ipar[15] = [[solution solutionInfoForKey:@"Bi-CGSTAB(l) polynomial degree"] intValue];
            if (ipar[15] < 2) errorfunct("FEMKernel_iterSolver", "Polynomial degree for Bi-CGSTAB(l) should be equal or greater than 2");
        } else {
            ipar[15] = 2;
        }
    }
    
    ipar[11] = 1;
    ipar[2] = n;
    
    if ([solution solutionInfoForKey:@"Linear system residual output"] != nil) {
        ipar[4] = [[solution solutionInfoForKey:@"Linear system residual output"] intValue];
    } else {
        ipar[4] = 1;
    }
    
    if ([solution solutionInfoForKey:@"Linear system maximum iterations"] != nil) {
        ipar[9] = [[solution solutionInfoForKey:@"Linear system maximum iterations"] intValue];
    } else {
        ipar[9] = 1;
    }
    
    work = doublematrix(0, n-1, 0, wsize-1);
    if (work == NULL) {
        errorfunct("FEMKernel_iterSolver", "Memory allocation error");
    }
    
    if (all([solution variableReturnPointerToValues], '=', 0.0, [solution variableSizeOfValues]) == 1) {
        for (i=0; i<[solution variableSizeOfValues]; i++) {
            [solution setVariableValues:i :1.0e-8];
        }
    }
    ipar[13] = 1;
    
    if ([solution solutionInfoForKey:@"Linear system convergence tolerance"] != nil) {
        dpar[0] = [[solution solutionInfoForKey:@"Linear system convergence tolerence"] doubleValue];
    }
    
    if ([solution solutionInfoForKey:@"Linear system divergence tolerance"] != nil) {
        dpar[1] = [[solution solutionInfoForKey:@"Linear system divergence tolerance"] doubleValue];
    } else {
        dpar[1] = HUGE_VAL;
    }
    
    if ([solution solutionInfoForKey:@"Linear system preconditioning"] != nil) {
        str = [NSString stringWithString:[solution solutionInfoForKey:@"Linear system preconditioning"]];
    } else {
        str = [NSString stringWithString:@"none"];
    }
    
    if ([str isEqualToString:@"none"] == YES) 
    {
        pCondType = PRECOND_NONE;
    }
    else if ([str isEqualToString:@"Diagonal"] == YES) 
    {
        pCondType = PRECOND_DIAGONAL;
    }
    else if ([str isEqualToString:@"ILUT"] == YES)
    {
        if ([solution solutionInfoForKey:@"Linear system ILUT tolerance"] != nil) {
            ilut_tol = [[solution solutionInfoForKey:@"Linear system ILUT tolerance"] doubleValue];
        } else {
            errorfunct("FEMKernel_iterSolver", "Linear system ILUT tolerance not found.");
        }
        pCondType = PRECOND_ILUT;
    }
    else if ([str isEqualToString:@"ILU"] == YES) 
    {
        if ([solution solutionInfoForKey:@"Linear system ILU order"] != nil) {
            ilun = [[solution solutionInfoForKey:@"Linear system ILU order"] intValue];
        } else {
            ilun = [[NSString stringWithFormat:@"%c", [str characterAtIndex:3]] intValue] -  0;
            if (ilun < 0 || ilun > 9 ) ilun = 0;
            pCondType = PRECOND_ILUN;
            
        }
    }
    else if ([str isEqualToString:@"BILU"] == YES)
    {
        ilun = [[NSString stringWithFormat:@"%c", [str characterAtIndex:4]] intValue] -  0;
        if (ilun < 0 || ilun > 9 ) ilun = 0;
        if ([solution variableDofs] == 1) {
            warnfunct("FEMKernel_iterSolver", "BILU for one dofs is equal to ILU!");
            pCondType = PRECOND_ILUN;
        } else {
            pCondType = PRECOND_BILUN;
        }
    }
    else if ([str isEqualToString:@"Multigrid"] == YES)
    {
        pCondType = PRECOND_MG;
    }
    else if ([str isEqualToString:@"Vanka"] == YES) 
    {
        pCondType = PRECOND_VANKA;
    }
    else {
        pCondType = PRECOND_NONE;
        warnfunct("FEMKernel_iterSolver", "Unknown preconditioner type, feature disabled.");
    }
    
    if ([solution solutionInfoForKey:@"No precondition recompute"] != nil) {
        if ([[solution solutionInfoForKey:@"No precondition recompute"] boolValue] == NO) {
            // To do's: Implement the code if we choose not to compute the preconditioner
            // for each method call.
        }
    }
    
    [solution setMatrixSolveCount:[solution matrixSolveCount]+1];
    
    if ([solution solutionInfoForKey:@"Linear system abort not converged"] != nil) {
        abortNotConverged = [[solution solutionInfoForKey:@"Linear system abort not converged"] boolValue];
    } else {
        abortNotConverged = YES;
    }
    
    // Get the selector for the matrix-vector multiplication method we want to use
    if ([solution matrixComplex] == 1) {
        mvSelector = @selector(CRS_ComplexMatrixVectorProd::::);
    } else {
        mvSelector= @selector(CRS_MatrixVectorProd::::);
    }
    
    // Get the selector for the preconditioning method we want to use
    if (pCondType == PRECOND_NONE) {
        
        pcondSelector = @selector(CRS_pcond_dummy::::);
    } 
    else if (pCondType == PRECOND_DIAGONAL) {
        if ([solution matrixComplex] == 1) {
            pcondSelector = @selector(CRS_ComplexDiagPrecondition::::);
        } else {
            pcondSelector = @selector(CRS_DiagPrecondition::::);
        }
    } 
    else if (pCondType == PRECOND_ILUN == pCondType == PRECOND_ILUT || pCondType == PRECOND_BILUN) {
        if ([solution matrixComplex] == 1) {
            pcondSelector = @selector(CRS_ComplexLUPrecondition::::);
        } else {
            pcondSelector = @selector(CRS_LUPrecondition::::);
        }
    }
    
    if ([solution matrixComplex] == 0) {
        
        if (iterType == ITER_SGS || iterType == ITER_JACOBI || iterType == ITER_GCR || iterType == ITER_BICGSTABL) {
            if (ipar[4] == 0) ipar[4] = HUGE_VAL;
        }
    } else {
        
        ipar[2] = ipar[2] / 2;
        if (iterType == ITER_GCR || iterType == ITER_BICGSTABL) {
            if (ipar[4] == 0) ipar[4] = HUGE_VAL;
        }
    }
    
    // Everything is happening in this method...
    [self iterCall:iterType :solution :ipar :dpar :work :pcondSelector :pcondrSelector :mvSelector :@selector(stopc:::::)];
    
    if ([solution matrixComplex] == 1) ipar[2] = ipar[2] * 2;
    
    if (ipar[29] != 1) {
        if (ipar[29] == 3) {
            errorfunct("FEMKernel_iterSolver", "System diverged over tolerance.");
        } else if (abortNotConverged == YES) {
            errorfunct("FEMKernel_iterSolver", "Failed convergence tolerances.");
        } else {
            errorfunct("FEMKernel_iterSolver", "Failed convergence tolerances.");
        }
    }
        
    free_dmatrix(work, 0, n-1, 0, wsize-1);
    work = NULL;
    
}

-(void)FEMKernel_backRotateNTSystem:(FEMSolution *)solution {
    
    int i, j, k, l, dim, ndofs;
    double bu, bv, bw, **rm;
    
    if ([solution normalTangentialNOFNodes] <= 0) return;
    
    dim = [solution coordinateSystemDimension];
    ndofs = [solution variableDofs];
    
    for (i=0; i<[solution sizeOfBoundaryReorder]; i++) {
        k = [solution boundaryReorder:i];
        if (k < 0) continue;
        j = [solution variablePerm:i];
        if (j < 0) continue;
        
        if (dim < 3) {
            
            bu = [solution variableValues:(ndofs*(j-1)+1)];
            bv = [solution variableValues:(ndofs*(j-1)+2)];
            
            [solution setVariableValues:(ndofs*(j-1)+1) :[solution boundaryNormals:k :0]*bu - [solution boundaryNormals:k :1]*bv];
            [solution setVariableValues:(ndofs*(j-1)+2) :[solution boundaryNormals:k :1]*bu + [solution boundaryNormals:k :0]*bv];
        } else {
            
            rm = doublematrix(0, 2, 0, 2);
            
            bu = [solution variableValues:(ndofs*(j-1)+1)];
            bv = [solution variableValues:(ndofs*(j-1)+2)];
            bw = [solution variableValues:(ndofs*(j-1)+3)];
            
            for (l=0; l<2; l++) {
                rm[0][l] = [solution boundaryNormals:k :l];
                rm[1][l] = [solution boundaryTangent1:k :l];
                rm[2][l] = [solution boundaryTangent2:k :l];
            }
            
            [solution setVariableValues:(ndofs*(j-1)+1) :rm[0][0]*bu + rm[1][0]*bv + rm[2][0]*bw];
            [solution setVariableValues:(ndofs*(j-1)+2) :rm[0][1]*bu + rm[1][1]*bv + rm[2][1]*bw];
            [solution setVariableValues:(ndofs*(j-1)+3) :rm[0][2]*bu + rm[1][2]*bv + rm[2][2]*bw];
            
            free_dmatrix(rm, 0, 2, 0, 2);
            
        }
    }
    
}

-(double)FEMKernel_computeNorm:(FEMSolution *)solution: (int)n: (double *)values {
    
    int normDim, normDofs, dofs, i, j, k, l, totn;
    double norm, nscale, sum;
    double *x, *buffer;
    
    FEMParallelMPI *parallelUtil;
    
    parallelUtil = [[FEMParallelMPI alloc] init];
    
    if (values != NULL) {
        x = values;
    } else {
        x = [solution variableReturnPointerToValues];
    }
    
    if ([solution solutionInfoForKey:@"Non linear system norm degree"] != nil) {
        normDim = [[solution solutionInfoForKey:@"Non linear system norm degree"] intValue];
    } else {
        normDim = 2;
    }
    
    dofs = [solution variableDofs];
    if ([solution solutionInfoForKey:@"Non linear system norm dofs"] != nil) {
        normDofs = [[solution solutionInfoForKey:@"Non linear system norm dofs"] intValue];
    } else {
        normDofs = dofs;
    }
    
    totn = [parallelUtil parallelReduction:(1.0*norm) :NULL];
    nscale = normDofs * totn/(1.0*dofs);
    
    // Zero is used instead of infinity as the sign for the max norm
    if (normDofs < dofs) {
        switch (normDim) {
            case 0:
                norm = 0.0;
                for (i=0; i<n/dofs; i++) {
                    for (j=0; j<normDofs; j++) {
                        k = dofs*(i) + j;
                        norm = max(norm, fabs(x[k]));
                    }
                }
                l = 2;
                norm = [parallelUtil parallelReduction:norm :&l];
                break;
            case 1:
                norm = 0.0;
                for (i=0; i<n/dofs; i++) {
                    for (j=0; j<normDofs; j++) {
                        k = dofs*(i) + j;
                        norm = norm + fabs(x[k]);
                    }
                }
                norm = [parallelUtil parallelReduction:norm :NULL] / nscale;
                break;
            case 2:
                norm = 0.0;
                for (i=0; i<n/dofs; i++) {
                    for (j=0; j<normDofs; j++) {
                        k = dofs*(i) + j;
                        norm = norm + pow(x[k], 2.0);
                    }
                }
                norm = sqrt([parallelUtil parallelReduction:norm :NULL] / nscale);
                break;
            default:
                norm = 0.0;
                for (i=0; i<n/dofs; i++) {
                    for (j=0; j<normDofs; j++) {
                         k = dofs*(i) + j;
                        norm = norm + pow(x[k], normDim);
                    }
                }
                norm = pow( ([parallelUtil parallelReduction:norm :NULL] / nscale), (1.0/normDim) );
                break;
        }
    } else {
        switch (normDim) {
            case 0:
                buffer = doublevec(0, n-1); 
                for (i=0; i<n; i++) {
                    buffer[i] = fabs(x[i]);
                }
                sum = max_array(buffer, n);
                l = 2;
                norm = [parallelUtil parallelReduction:sum :&l];
                free_dvector(buffer, 0, n-1);
                break;
            case 1:
                sum = 0.0;
                for (i=0; i<n; i++) {
                    sum = sum + fabs(x[i]);
                }
                norm = [parallelUtil parallelReduction:sum :NULL] / nscale;
                break;
            case 2:
                sum = 0.0;
                for (i=0; i<n; i++) {
                    sum = sum + pow(x[i], 2.0);
                }
                norm = sqrt([parallelUtil parallelReduction:sum :NULL] / nscale);
                break;
            default:
                sum = 0.0;
                for (i=0; i<n; i++) {
                    sum = sum + pow(x[i], normDim);
                }
                norm = pow( ([parallelUtil parallelReduction:sum :NULL] / nscale), (1.0/normDim));
                break;
        }
    }
    
    return norm;
    
}


-(void)FEMKernel_computeChange:(FEMSolution *)solution :(BOOL)steadyState :(int *)nsize :(double *)values :(double *)values0 {
    
    NSString *convergenceType, *solverName;
    int i, n, relaxAfter, iterNo;
    double norm, prevNorm, bNorm, change, relaxation, maxNorm, dt, tolerance, eps;
    double *rhsvec, *x, *r, *x0;
    BOOL skip, convergenceAbsolute, relax, relaxBefore, stat, doIt;
    Variable_t *iterV, *timeStepVar, *veloVar;
    char *name, *str1, *str2;
    
    FEMPrecondition *preconditioning;
    
    preconditioning = [[FEMPrecondition alloc] init];
    
    relax = NO;
    
    if (steadyState == YES) {
        skip = [[solution solutionInfoForKey:@"Skip compute steady state change"] boolValue];
        if (skip == YES) return;
        
        if ([solution solutionInfoForKey:@"Steady state convergence measure"] != nil) {
            convergenceType = [NSString stringWithString:[solution solutionInfoForKey:@"Steady state convergence measure"]];
        } else {
            convergenceType = [NSString stringWithString:@"norm"];
        }
        
        if ([solution solutionInfoForKey:@"Steady state convergence absolute"] != nil) {
            convergenceAbsolute = [[solution solutionInfoForKey:@"Steady state convergence absolute"] boolValue];
        } else if ([solution solutionInfoForKey:@"Use absolute norm for convergence"] != nil) {
            convergenceAbsolute = [[solution solutionInfoForKey:@"Use absolute norm for convergence"] boolValue];
        }
        
        if ([solution solutionInfoForKey:@"Steady state relaxation factor"] != nil) {
            relaxation = [[solution solutionInfoForKey:@"Steady state relaxation factor"] doubleValue];
            relax = (relaxation != 1.0) ? YES: NO;
        }
        
        name = "coupled iter";
        iterV = getVariable([solution meshReturnPointerToVariables], name);
        iterNo = (int)iterV->Values[0];
        
        if (relax == YES) {
            if ([solution solutionInfoForKey:@"Steady state relaxation after"] != nil) {
                relaxAfter = [[solution solutionInfoForKey:@"Steady state relaxation before"] intValue];
                if (relaxAfter >= iterNo ) relax = NO;
            }    
        }
        
        if (relax == YES) {
            if ([solution solutionInfoForKey:@"Steady state relaxation before"] != nil) {
                relaxBefore = [[solution solutionInfoForKey:@"Steady state relaxation before"] boolValue];
            } else {
                relaxBefore = YES;
            }
        }
    } else {
    
        if ([solution solutionInfoForKey:@"Skip compute non linear change"] != nil) {
            skip  = [[solution solutionInfoForKey:@"Skip compute non linear change"] boolValue];
        } else {
            skip = NO;
        }
        if (skip == YES) return;
        
        if ([solution solutionInfoForKey:@"Non linear system convergence measure"] != nil) {
            convergenceType = [solution solutionInfoForKey:@"Non linear system convergence measure"];
        } else {
            convergenceType = @"norm";
        }
        
        if ([solution solutionInfoForKey:@"Non linear system convergence absolute"] != nil) {
            convergenceAbsolute = [[solution solutionInfoForKey:@"Non linear system convergence absolute"] boolValue];
        } else if ([solution solutionInfoForKey:@"Use absolute norm for convergence"] != nil) {
            convergenceAbsolute = [[solution solutionInfoForKey:@"Use absolute norm for convergence"] boolValue];
        } else {
            convergenceAbsolute = NO;
        }
        
        name = "nonlin iter";
        iterV = getVariable([solution meshReturnPointerToVariables], name);
        iterNo = (int)iterV->Values[0];
        [solution setVariableNonLinIter:(int)iterV->Values[0]];
        iterV->Values[0] = iterV->Values[0]+1;
        
        if ([solution solutionInfoForKey:@"Non linear system relaxation factor"] != nil) {
            relaxation = [[solution solutionInfoForKey:@"Non linear system relaxation factor"] doubleValue];
            relax = (relaxation != 1.0) ? YES: NO;
        }
        if (relax == YES) {
            if ([solution solutionInfoForKey:@"Non linear system relaxation after"] != nil) {
                relaxAfter = [[solution solutionInfoForKey:@"Non linear system relaxation after"] boolValue];
                if (relaxAfter >= [solution variableNonLinIter]) relax = NO;
            }
        }
        
        if (relax == YES) {
            if ([solution solutionInfoForKey:@"Non linear system relaxation before"] != nil) {
                relaxAfter = [[solution solutionInfoForKey:@"Non linear system relaxation before"] boolValue];
                if (relaxAfter >= [solution variableNonLinIter]) relax = NO;
            }
        }
        
    }
    
    if (values != NULL) {
        x = values;
    } else {
        x = [solution variableReturnPointerToValues];
    }
    
    if (values == NULL) {
        [solution setVariableNorm:0.0];
        if (steadyState == YES) {
            [solution setVariableSteadyChange:0.0];
        } else {
            [solution setVariableNonLinChange:0.0];
        }
        return;
    }
    
    if (nsize != NULL && values != NULL) {
        n = *nsize;
    } else if ([solution variableReturnPointerToValues] != NULL){
        n = [solution variableSizeOfValues];
    } 
    
    stat = NO;
    if (values0 != NULL) {
        x0 = values0;
        stat = YES;
    } else if (steadyState == YES) {
        if ([solution isAssociatedVariableSteadyValues] == YES) {
            x0 = [solution variableReturnPointerToSteadyValues];
            stat = YES;
        }
    } else {
        if ([solution isAssociatedVariableNonLinValues] == YES) {
            x0 = [solution variableReturnPointerToNonLinValues];
            stat = YES;
        }
    }
    
    // Possibly check here if length mismatch between x0 and x
    //.....
    
    if (relax == YES && relaxBefore == YES) {
        for (i=0; i<n; i++) {
            x[i] = (1-relaxation) * x0[i] + relaxation*x[i];
        }
    }
    
    if (steadyState == YES) {
        prevNorm = [solution variablePrevNorm];
    } else {
        prevNorm = [solution variableNorm];
    }
    
    // Compute norm here
    norm = [self FEMKernel_computeNorm:solution :n :x];
    [solution setVariableNorm:norm];
    
    // The norm should be bounded in order to reach convergence
    if ([solution solutionInfoForKey:@"Non linear system max norm"] != nil) {
        maxNorm = [[solution solutionInfoForKey:@"Non linear system max norm"] doubleValue];
    } else {
        maxNorm = HUGE_VALF;
    }
    
    if (isnan(norm) != 0 && norm > maxNorm) {
        warnfunct("FEMKernel_computeChange", "Computed norm:");
        printf("%lf\n", norm);
        errorfunct("FEMKernel_computeChange", "Norm of solution has crossed given bounds.");
    }
    
    if ([convergenceType isEqualToString:@"residual"] == YES) {
        // ------------------------------------------------------------------------------
        // x is solution of A(x0)x = b(x0), thus residual should be real r = b(x)-A(x)x
        // Instead we use r = b(x0)-A(x0)x0 which unfortunately is one step behind
        // ------------------------------------------------------------------------------
        rhsvec = [solution matrixReturnPointerToRHS];
        r = doublevec(0, n-1);
        
        [preconditioning CRS_MatrixVectorMultiply:solution :x0 :r];
        for (i=0; i<n; i++) {
            r[i] = r[i] - rhsvec[i];
        }
        change = [self FEMKernel_computeNorm:solution :n :r];
        if (convergenceAbsolute == NO) {
            bNorm = [self FEMKernel_computeNorm:solution :n :rhsvec];
            if (bNorm > 0.0) {
                change = change / bNorm;
            }
        }
        free_dvector(r, 0, n-1);
        rhsvec = NULL;
    }
    else if ([convergenceType isEqualToString:@"linear system residual"] == YES) {
        // ------------------------------------------------------------------------------
        // Here the true linear system redisual r = b(x)-A(x)x is computed.
        // This option is useful for some special solvers
        // ------------------------------------------------------------------------------
        rhsvec = [solution matrixReturnPointerToRHS];
        r = doublevec(0, n-1);
        
        [preconditioning CRS_MatrixVectorMultiply:solution :x :r];
        for (i=0; i<n; i++) {
            r[i] = r[i] - rhsvec[i];
        }
        change = [self FEMKernel_computeNorm:solution :n :r];
        if (convergenceAbsolute == NO) {
            bNorm = [self FEMKernel_computeNorm:solution :n :rhsvec];
            if (bNorm > 0.0) {
                change = change / bNorm;
            }
        }
        free_dvector(r, 0, n-1);
        rhsvec = NULL;
    }
    else if ([convergenceType isEqualToString:@"solution"] == YES) {
        r = doublevec(0, n-1);
        
        for (i=0; i<n; i++) {
            r[i] = x[i] - x0[i];
        }
        change = [self FEMKernel_computeNorm:solution :n :r];
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
        warnfunct("FEMKernel_computeChange", "Unknown convergence measure:");
        printf("%s\n", [convergenceType UTF8String]);
    }
    
    // Check for convergence: 0/1
    if (steadyState == YES) {
        [solution setVariableSteadyChange:change];
        if ([solution solutionInfoForKey:@"Steady state convergence tolerance"] != nil) {
            if (change <= tolerance) {
                [solution setVariableSteadyConverged:1];
            } else {
                [solution setVariableSteadyConverged:0];
            }
        }
    } else {
        [solution setVariableNonLinChange:change];
        if ([solution solutionInfoForKey:@"Non linear system convergence tolerance"] != nil) {
            if (change <= tolerance) {
                [solution setVariableNonLinConverged:1];
            } else {
                [solution setVariableNonLinConverged:0];
            }
        }
    }
    
    if (relax == YES && relaxBefore == NO) {
        for (i=0; i<n; i++) {
            x[i] = (1-relaxation)*x0[i] + relaxation*x[i];
        }
        [solution setVariableNorm:[self FEMKernel_computeNorm:solution :n :x]];
    }
    
    if ([solution solutionInfoForKey:@"Equation"] != nil) {
        solverName = [NSString stringWithString:[solution solutionInfoForKey:@"Equation"]];
    } else {
        solverName = [[NSString alloc] initWithCString:[solution variableName] encoding:NSMacOSRomanStringEncoding];
    }
    
    if (steadyState == YES) {
        NSLog(@"FEMKernel_computeChange:SS (Iter=%d) (NRM,RELC): (%e %e):: %@\n", iterNo, norm, change, solverName);
    } else {
         NSLog(@"FEMKernel_computeChange:NS (Iter=%d) (NRM,RELC): (%e %e):: %@\n", iterNo, norm, change, solverName);
    }
    
    // Only 1st order velocity computation is implemented so far...
    if ([solution timeOrder] == 1) {
        doIt = YES;
        
        if (steadyState == YES) {
            if ([solution solutionInfoForKey:@"Calculate velocity"] != nil) {
                doIt = [[solution solutionInfoForKey:@"Calculate velocity"] boolValue];
            } else {
                doIt = NO;
            }
        } else {
            if ([solution solutionInfoForKey:@"Non Linear calculate velocity"] != nil) {
                doIt = [[solution solutionInfoForKey:@"Non Linear calculate velocity"] boolValue];
            } else {
                doIt = NO;
            }
        }
        
        if (doIt == YES) {
            name = "timestep size";
            timeStepVar = getVariable([solution meshReturnPointerToVariables], name);
            dt = timeStepVar->Values[0];
            str2 = " velocity";
            str1 = strcat([solution variableName], str2);
            veloVar = getVariable([solution meshReturnPointerToVariables], str1);
            for (i=0; i<n; i++) {
                veloVar->Values[i] = (x[i] - [solution variablePrevValues:i :0]) / dt;
            }
        }
        
    }
    
    // Calculate derivative a.k.a sensitivity
    if (steadyState == YES) {
    
        stat = NO;
        
        if ([solution solutionInfoForKey:@"Calculate derivative"] != nil) {
            if ([[solution solutionInfoForKey:@"Calculate derivative"] boolValue] == YES) {
                
                if (iterNo > 1) {
                    name = "derivative eps";
                    timeStepVar = getVariable([solution meshReturnPointerToVariables], name);
                    if (timeStepVar != NULL) {
                        eps = timeStepVar->Values[0];
                        stat = YES;
                    } else {
                        eps = [[solution solutionInfoForKey:@"derivative eps"] doubleValue];
                    }
                    if (stat == NO) {
                        warnfunct("FEMKernel_computeChange", "Derivative eps not given, using one.");
                    }
                    
                    str2 = " derivative";
                    str1 = strcat([solution variableName], str2);
                    veloVar = getVariable([solution meshReturnPointerToVariables], str1);
                    if (veloVar != NULL) {
                        NSLog(@"FEMKernel_computeChange: Computing variable: %s\n", str1);
                        for (i=0; i<n; i++) {
                            veloVar->Values[i] = (x[i] - x0[i]) / eps;
                        }
                    } else {
                        warnfunct("FEMKernel_computeChange", "Derivative variable not present.");
                    }
                }
            }
        }
    }
    
    x = NULL;
    x0 = NULL;
    
}

-(void)FEMKernel_solveLinearSystem:(FEMSolution *)solution {
    
    int i, j, n;
    
    double *diag, *x;
    double diagReal, diagImag, norm, bnorm, sum;
    double complex cmpx;
    NSString *method;
    
    FEMParallelMPI *parallelUtil;
    
    BOOL scaleSystem, eigenAnalysis, harmonicAnalysis;
    
    parallelUtil = [[FEMParallelMPI alloc] init];
    
    n = [solution matrixNumberOfRows];
    
    if ([solution matrixLumped] == 1 && [solution timeOrder] == 1)  {
        
        x = [solution variableReturnPointerToValues];
        
        if ([solution solutionInfoForKey:@"Time stepping Method"] != nil) {
            method = [NSString stringWithString:[solution solutionInfoForKey:@"Time stepping Method"]];
            if ([method isEqualToString:@"runge-kutta"] == YES || [method isEqualToString:@"explicit euler"] == YES) {
                for (i=0; i<n; i++) {
                    if (fabs([solution matrixValues:[solution matrixDiag:i]]) > 0.0) [solution setVariableValues:i :[solution matrixRHS:i] / [solution matrixValues:[solution matrixDiag:i]]];
                }
                [self FEMKernel_backRotateNTSystem:solution];
                [self FEMKernel_computeNorm:solution :n :x];
                x = NULL;
                return;
            }
        }
        
    }
    
    if ([solution solutionInfoForKey:@"Linear system scaling"] != nil) {
        scaleSystem = [[solution solutionInfoForKey:@"Linear system scaling"] boolValue];
    } else {
        scaleSystem = YES;
    }
    
    eigenAnalysis = ([solution nofEigenValues] > 0 && [[solution solutionInfoForKey:@"Eigen analysis"] boolValue] == YES) ? YES: NO;
    harmonicAnalysis = ([solution nofEigenValues] > 0 && [[solution solutionInfoForKey:@"Harmonic analysis"] boolValue] == YES) ? YES: NO;
    
    if ( harmonicAnalysis == NO && eigenAnalysis == NO ) {
        sum = 0.0;
        for (i=0; i<n; i++) {
            sum = sum + pow([solution matrixRHS:i], 2.0);
        }
        bnorm = [parallelUtil parallelReduction:sqrt(sum) :NULL];
        if (bnorm == 0.0) {
            warnfunct("SolveSystem", "Solution trivially zero");
            for (i=0; i<n; i++) {
                [solution setVariableValues:i :0.0];
            }
            return;
        }
    }
    
    if ( scaleSystem == YES && harmonicAnalysis == NO ) {
        
        // Scale System Ax = b as:
        // (DAD)y = Db, where D = 1/sqrt(diag(A)) and y = D^-1 x
        diag = doublevec(0, n-1);
        
        if ([solution matrixComplex] == 1) {
            for (i=0; i<n; i+=2) {
                j = [solution matrixDiag:i];
                diag[i] = [solution matrixValues:j];
                diag[i+1] = [solution matrixValues:j+1];
            }
        } else {
            for (i=0; i<n; i++) {
                diag[i] = [solution matrixValues:[solution matrixDiag:i]];
            }
        }
        
        if ([solution matrixComplex] == 1) {
            for (i=0; i<n; i+=2) {
                diagReal = diag[i];
                diagImag = -diag[i+1];
                cmpx = diagReal + diagImag * I;
                if (cabs(cmpx) != 0.0) {
                    diag[i] = 1.0 / (sqrt(cabs(cmpx)));
                    diag[i+1] = 1.0 / (sqrt(cabs(cmpx)));
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
        
        for (i=0; i<n; i++) {
            for (j=[solution matrixRows:i]; j<=[solution matrixRows:i+1]-1; j++) {
                [solution setMatrixValues:j :[solution matrixValues:j] * (diag[i] * diag[[solution matrixCols:j]])];
            }
        }
        
        if ([solution isAssociatedMatrixMassValues] == YES) {
            if ([solution matrixSizeOfValues] == [solution matrixSizeOfMassValues]) {
                for (i=0; i<n; i++) {
                    for (j=[solution matrixRows:i]; j<=[solution matrixRows:i+1]-1; j++) {
                        [solution setMatrixMassValues:j :[solution matrixMassValues:j] * (diag[i] * diag[[solution matrixCols:j]])];
                    }
                }
            }
        }
        
        if ([solution isAssociatedMatrixDampValues] == YES) {
            if ([solution matrixSizeOfValues] == [solution matrixSizeOfDampValues]) {
                for (i=0; i<n; i++) {
                    for (j=[solution matrixRows:i]; j<=[solution matrixRows:i+1]-1; j++) {
                        [solution setMatrixDampValues:j :[solution matrixDampValues:j] * (diag[i] * diag[[solution matrixCols:j]])];
                    }
                }
            }
        }
    }
        
    // If solving harmonic anaysis go there:
    if (harmonicAnalysis == YES) {
        
        // To do's
        
    }
    
    // If solving eigen system go there:
    if (eigenAnalysis) {
        
        // To do's
        
    }
    
    // Convert rhs & initial value to the scaled system
    if (scaleSystem == YES) {
        for (i=0; i<n; i++) {
            [solution setMatrixRHS:i :[solution matrixRHS:i] * diag[i]];
        }
    }
    sum = 0.0;
    for (i=0; i<n; i++) {
        sum = sum + pow([solution matrixRHS:i], 2.0);
    }
    bnorm = [parallelUtil parallelReduction:sqrt(sum) :NULL];
    
    if (scaleSystem == YES) {
        for (i=0; i<n; i++) {
            [solution setMatrixRHS:i :[solution matrixRHS:i] / bnorm];
            [solution setVariableValues:i :[solution variableValues:i] / diag[i] / bnorm];
        }
    }
    
    if ([solution solutionInfoForKey:@"Linear system solver"] != nil) {
        method = [NSString stringWithString:[solution solutionInfoForKey:@"Linear system solver"]];
    } else {
        method = [NSString stringWithString:@"iterative"];
    }
   
    if ([method isEqualToString:@"iterative"] == YES) {
        
        [self FEMKernel_iterSolver:solution];
    } 
    else if ([method isEqualToString:@"multigrid"] == YES) {
        // Need to be implemented
    }
    else { // By default we use the direct solver
        // Need to be implemented
    }
    
    if (scaleSystem == YES) {
        
        // Solve x: INV(D)x = y, scale b back to original
        for (i=0; i<n; i++) {
            [solution setVariableValues:i :[solution variableValues:i] * diag[i] * bnorm];
            [solution setMatrixRHS:i :[solution matrixRHS:i] / diag[i] * bnorm];
        }
        
        // Scale the system back to original
        for (i=0; i<n; i++) {
            for (j=[solution matrixRows:i]; j<=[solution matrixRows:i+1]-1; j++) {
                [solution setMatrixValues:j :[solution matrixValues:j] / (diag[i] * diag[[solution matrixCols:j]])];
                
            }
        }
        
        if ([solution isAssociatedMatrixMassValues] == YES) {
            if ([solution matrixSizeOfValues] == [solution matrixSizeOfMassValues]) {
                for (i=0; i<n; i++) {
                    for (j=[solution matrixRows:i]; j<=[solution matrixRows:i+1]-1; j++) {
                        [solution setMatrixMassValues:j :[solution matrixMassValues:j] / (diag[i] * diag[[solution matrixCols:j]])];
                    }
                }
            }
        }
        
        if ([solution isAssociatedMatrixDampValues] == YES) {
            if ([solution matrixSizeOfValues] == [solution matrixSizeOfDampValues]) {
                for (i=0; i<n; i++) {
                    for (j=[solution matrixRows:i]; j<=[solution matrixRows:i+1]-1; j++) {
                        [solution setMatrixDampValues:j :[solution matrixDampValues:j] / (diag[i] * diag[[solution matrixCols:j]])];
                    }
                }
            }
        }
    
        free_dvector(diag, 0, n-1);
        
    }
    
    // To do's: Add support for the treatment of variable loads and constrained matrix here
    // ....
    
    [self FEMKernel_backRotateNTSystem:solution];
    
    // Compute the change of the solution with different methods
    [self FEMKernel_computeChange:solution :NO :&n :[solution variableReturnPointerToValues] :NULL];
    norm = [solution variableNorm];
    
}

-(void)FEMKernel_solveSystem:(FEMSolution *)solution {
    
    int n, i;
    double relaxation, beta, gamma;
    double t0, rt0, st, rst;
    BOOL needPrevSol, stat;
    NSString *method;
    
    if ([solution solutionInfoForKey:@"Linear system timing"] != nil) {
        if ([[solution solutionInfoForKey:@"Linear system timing"] boolValue] == YES) {
            t0 = cputime();
            rt0 = realtime();
        }
    }
    
    n = [solution matrixNumberOfRows];

    // The allocation of previous values has to be here in order to work properly
    // with the Dirichlet elimination
    if ([solution solutionInfoForKey:@"Non linear system relaxation factor"] != nil) {
        relaxation = [[solution solutionInfoForKey:@"Non linear system relaxation factor"] doubleValue];
        needPrevSol = (relaxation != 0.0) ? YES: NO;
    }
    
    if (needPrevSol == NO) {
        if ([solution solutionInfoForKey:@"Non linear system convergence measure"] != nil) {
            method = [NSString stringWithString:[solution solutionInfoForKey:@"Non linear system convergence measure"]];
            needPrevSol = ([method isEqualToString:@"residual"] == YES || [method isEqualToString:@"solution"] == YES) ? YES: NO;
        }
    }
    
    if (needPrevSol == YES) {
        if ([solution isAssociatedVariableNonLinValues] == YES) {
            stat = YES;
            if ([solution variableSizeOfNonLinValues] != n) {
                [solution freeVariableNonLinValues:[solution variableSizeOfNonLinValues]];
                stat = NO;
            }
        }
        if (stat == NO) {
            [solution allocateVariableNonLinValues:n];
            if ([solution isAssociatedVariableNonLinValues] == NO) errorfunct("FEMKernel_solveSystem", "Memory allocation error.");
            [solution setVariableSizeOfNonLinValues:n];
        }
        for (i=0; i<n; i++) {
            [solution setVariableNonLinValues:i :[solution variableValues:i]];
        }
    }
    
    
    // To do's: Add support for constrained matrix
    
    [self FEMKernel_solveLinearSystem:solution];
    
    if ([solution timeOrder] == 2) {
        if ([solution isAssociatedVariablePrevValues] == YES) {
            
            gamma = 0.5 - [solution alpha];
            beta = ( pow((1.0 - [solution alpha]), 2.0) ) / 4.0;
            for (i = 0; i<n; i++) {
                [solution setVariablePrevValues:i :1 :(1.0/(beta*pow([solution dt], 2.0))) * ([solution variableValues:i]-[solution variablePrevValues:i :2]) - 
                 (1.0/(beta*[solution dt]))*[solution variablePrevValues:i :3] +(1.0-1.0/(2.0*beta))*[solution variablePrevValues:i :4]];
                
                [solution setVariablePrevValues:i :0 :[solution variablePrevValues:i :3] + [solution dt]*((1.0-gamma)*[solution variablePrevValues:i :4] + 
                                                                                                          gamma*[solution variablePrevValues:i :1])];
            }
            
        }
    }
    
    if ([solution solutionInfoForKey:@"Linear system timing"] != nil) {
        if ([[solution solutionInfoForKey:@"Linear system timing"] boolValue] == YES) {
            st = cputime() - t0;
            rst = realtime() - rt0;
            
            NSLog(@"Solution time for %s: %lf %lf (s)\n", [solution variableName], st, rst);
            
            // // To do's: Add support for cumulative timing
        
        }
    }
    
}

#pragma mark Element info

-(BOOL)FEMKernel_checkPassiveElement:(FEMModel *)model: (FEMSolution *)solution: (Element_t *)element {
    
    int body_id, bf_id, nbNodes;
    double *passive;
    NSArray *bodiesValues, *bodyForces;
    char *passName, *str;
    NSString *aString;
    BOOL isPassive, found;
    
    FEMValueList *valueList;
    FEMBodyForce *bodyForceAtBodyID;
    
    isPassive = NO;
    
    body_id = element->BodyID-1;
    if (body_id < 0) return isPassive; // body_id == -1 for boundary elements
    
    valueList = [[FEMValueList alloc] init];
    
    nbNodes = element->Type.NumberOfNodes;
    passive = doublevec(0, nbNodes-1);
    
    bodiesValues = [model bodiesValuesForKey:@"Body force"];
    if ([bodiesValues objectAtIndex:body_id] == nil) {
        return isPassive;
    } else {
        bf_id = [[bodiesValues objectAtIndex:body_id] intValue];
    }
    
    passName = [solution variableName];
    str = " passive";
    strcat(passName, str);
    
    bodyForces = [model returnBodyForces];
    // Returns a bodyForce object at index bf_id
    bodyForceAtBodyID = [bodyForces objectAtIndex:bf_id-1];
    aString = [NSString stringWithCString:str encoding:NSUTF8StringEncoding];
    
    found = [valueList listGetReal:model :[bodyForceAtBodyID returnValuesList] :aString :nbNodes :element->NodeIndexes :passive];
    if (found == YES) {
        if ( count(passive, '>', 0.0, nbNodes) > count(passive, '<', 0.0, nbNodes) ) {
            isPassive = YES;
        }
    
    }
    
    free_dvector(passive, 0, nbNodes-1);
    
    return isPassive;
    
}

#pragma mark Manipulate matrix

-(void)FEMKernel_rotateMatrix:(FEMSolution *)solution: (double **)matrix: (double *)vector: (int)n: (int)dim: (int)dofs: (int *)nodeIndexes {
    
    int i, j, k, l;
    double s, **r, **q, *n1, *t1, *t2;
    
    r = doublematrix(0, (n*dofs)-1, 0, (n*dofs)-1);
    q = doublematrix(0, (n*dofs)-1, 0, (n*dofs)-1);
    n1 = doublevec(0, 2);
    t1 = doublevec(0, 2);
    t2 = doublevec(0, 2);
    
    for (i=0; i<n; i++) {
        
        if (nodeIndexes[i] < 0 || nodeIndexes[i]+1 > [solution size1boundaryNormals]) continue;
        
        for (j=0; j<n*dofs; j++) {
            for (k=0; k<n*dofs; k++) {
                r[j][k] = 0.0;
            }
        }
        
        for (j=0; j<n*dofs; j++) {
            r[j][j] = 1.0;
        }
        
        for (j=0; j<2; j++) {
            n1[j] = [solution boundaryNormals:nodeIndexes[i] :j];
        }
        
        switch (dim) {
            case 2:
                r[dofs*(i-1)+1][dofs*(i-1)+1] = n1[0];
                r[dofs*(i-1)+1][dofs*(i-1)+2] = n1[1];
                
                r[dofs*(i-1)+2][dofs*(i-1)+1] = -n1[1];
                r[dofs*(i-1)+2][dofs*(i-1)+2] = n1[0];
                break;
            case 3:
                for (j=0; j<2; j++) {
                    t1[j] = [solution boundaryTangent1:nodeIndexes[i] :j];
                    t2[j] = [solution boundaryTangent2:nodeIndexes[i] :j];
                }
                
                r[dofs*(i-1)+1][dofs*(i-1)+1] = n1[0];
                r[dofs*(i-1)+1][dofs*(i-1)+2] = n1[1];
                r[dofs*(i-1)+1][dofs*(i-1)+3] = n1[2];
                
                r[dofs*(i-1)+2][dofs*(i-1)+1] = t1[0];
                r[dofs*(i-1)+2][dofs*(i-1)+2] = t1[1];
                r[dofs+(i-1)+2][dofs*(i-1)+3] = t1[2];
                
                r[dofs*(i-1)+3][dofs*(i-1)+1] = t2[0];
                r[dofs*(i-1)+3][dofs*(i-1)+2] = t2[1];
                r[dofs*(i-1)+3][dofs*(i-1)+3] = t2[2];
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
        
        for (j=0 ; j<n*dofs; j++) {
            vector[j] = q[j][0];
        }
    }
    
    free_dmatrix(r, 0, (n*dofs)-1, 0, (n*dofs)-1);
    free_dmatrix(q, 0, (n*dofs)-1, 0, (n*dofs)-1);
    free_dvector(n1, 0, 2);
    free_dvector(t1, 0, 2);
    free_dvector(t2, 0, 2);
    
}

-(void)FEMKernel_updateGlobalForce:(FEMModel *)model: (FEMSolution *)solution: (Element_t *)element: (double *)forceVector: (double *)localForce: (int)n: (int)dofs: (int *)nodeIndexes: (BOOL *)rotateNT {
    
    int i, j, k, dim, *indexes;
    BOOL rotate;
    double **localStiffMatrix;
    
    // Check if this element has been defined as passive
    if ([self FEMKernel_checkPassiveElement:model :solution :element] == YES) return;
    
    localStiffMatrix = doublematrix(0, (n*dofs)-1, 0, (n*dofs)-1);
    indexes = intvec(0, n-1);
    
    rotate = YES;
    if (rotateNT != NULL) rotate = *rotateNT;
    
    if (rotate == YES && [solution normalTangentialNOFNodes] > 0) {
        
        dim = [model dimension];
         memset( indexes, 0.0, (n*sizeof(indexes)) );
        for (i=0; i<element->Type.NumberOfNodes; i++) {
            indexes[i] = [solution boundaryReorder:element->NodeIndexes[i]];
        }
        [self FEMKernel_rotateMatrix:solution :localStiffMatrix :localForce :n :dim :dofs :indexes];
    
    }
    
    for (i=0; i<n; i++) {
        if (nodeIndexes[i] >= 0) {
            for (j=0; j<dofs; j++) {
                k = dofs*nodeIndexes[i] + j;
                forceVector[k] = forceVector[k] + localForce[dofs*(i-1)+j];
            }
        }
        
    }
    
    free_dmatrix(localStiffMatrix, 0, (n*dofs)-1, 0, (n*dofs)-1);
    free_ivector(indexes, 0, n-1);
    
}

-(void)FEMKernel_addFirstOrderTime:(FEMModel *)model: (FEMSolution *)solution: (Element_t *)element: (double **)massMatrix: (double **)stiffMatrix: (double *)force: (double)dt: (int)n: (int)dofs: (int *)nodeIndexes: (int *)rows: (int *)cols {
/********************************************************************************************************************
    For time dependent simulations add the time derivative coefficient terms to the matrix
    containing other coefficients.
********************************************************************************************************************/
    
    int i, j, k, l, m, order;
    double s, t;
    double **prevSol, *matrixForce, *lForce, *buffer;
    double *dts;
    BOOL constantDt;
    FEMTimeIntegration *timeIntegration;
    NSString *method;
    char *name;
    Variable_t *dtVar;
    
    prevSol = doublematrix(0, (dofs*n)-1, 0, [solution order]-1);
    lForce = doublevec(0, (dofs*n)-1);
    dts = doublevec(0, [solution order]-1);
    matrixForce = doublevec(0, [solution matrixSize1OfForce]-1);
    
    buffer = doublevec(0, (dofs*n)-1);
    
    timeIntegration = [[FEMTimeIntegration alloc] init];
    
    if ([solution matrixLumped] == 1) {
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
    
    order = min([solution doneTime], [solution order]);
    
    for (i=0; i<n; i++) {
        for (j=0; j<dofs; j++) {
            k = dofs * i + j;
            l = dofs * nodeIndexes[i] + j;
            for (m=0; n<order; m++) {
                prevSol[k][m] = [solution variablePrevValues:l :m];
            }
        }
    }
    
    for (i=0; i<n*dofs; i++) {
        lForce[i] = force[i];
    }
    
    for (i=0; i<[solution matrixSize1OfForce]; i++) {
        matrixForce[i] = [solution matrixForce:i :0];
    }
    
    [self FEMKernel_updateGlobalForce:model :solution :element :matrixForce :lForce :n :dofs :nodeIndexes :NULL];
    
    if ([solution solutionInfoForKey:@"Time stepping method"] != nil) {
        method = [NSString stringWithString:[solution solutionInfoForKey:@"Time stepping method"]];
    } else {
        method = [NSString stringWithString:@"bdf"];
    }

    if ([method isEqualToString:@"fs"] == YES) {
        
        for (i=0; i<dofs*n; i++) {
            buffer[i] = prevSol[i][0];
        }
        [timeIntegration fractionStep:solution :n*dofs :dt :massMatrix :stiffMatrix :force :buffer :rows];
    } 
    else if ([method isEqualToString:@"bfs"]) {
        dts[0] = dt;
        constantDt = NO;
        if (order > 1) {
            name = "time step size";
            dtVar = getVariable([solution meshReturnPointerToVariables], name);
            for (i=1; i<order; i++) {
                dts[i] = dtVar->PrevValues[0][i-1];
                if ( fabs(dts[i]-dts[0]) > 1.0e-6 * dts[0] ) constantDt = NO;
            }
        }
        if (constantDt == YES) {
        
            [timeIntegration bdfLocal:solution :n*dofs :dt :massMatrix :stiffMatrix :force :prevSol :order :rows :cols];
        
        } else {
            
            [timeIntegration vbdfLocal:solution :n*dofs :dts :massMatrix :stiffMatrix :force :prevSol :order :rows :cols];
            
        }
    }
    else {
        
        for (i=0; i<dofs*n; i++) {
            buffer[i] = prevSol[i][0];
        }
        [timeIntegration newMarkBeta:solution :n*dofs :dt :massMatrix :stiffMatrix :force :buffer :[solution beta] :rows];
        
    }
    
    free_dmatrix(prevSol, 0, (dofs*n)-1, 0, [solution order]-1);
    free_dvector(lForce, 0, (dofs*n)-1);
    free_dvector(dts, 0, [solution order]-1);
    free_dvector(matrixForce, 0, [solution matrixSize1OfForce]-1);
    free_dvector(buffer, 0, (dofs*n)-1);
    
}

-(void)FEMKernel_updateGlobalEquations:(FEMModel *)model: (FEMSolution *)solution: (Element_t *)element: (double **)stiffMatrix: (double *)force: (int)n: (int)dofs: (int *)nodeIndexes: (int *)rows: (int *)cols: (BOOL *)rotateNT: (BOOL *)bulkUpdate {
/********************************************************************************************************************
    Add element local matrices and vectors to global matrices and vectors
    
    Arguments:
 
        FEMSolution *solution     -> Class solution which contains the global matrix and global RHS vector
        double **stiffMatrix      -> Local matrix to be added tp the global matrix
        double  *force            -> Element local force vector
        int n                     -> Number of nodes
        int dofs                  -> Number of dofs
        int *nodeIndexes          -> Element node to global node numbering mapping
 
********************************************************************************************************************/
    
    int i, j, k, dim;
    int *indexes;
    BOOL rotate;
    FEMMatrixCRS *crsMatrix;
    FEMMatrixBand *bandMatrix;
    
    // Check if this element has been defined as passive
    if ([self FEMKernel_checkPassiveElement:model :solution :element] == YES) return;
    
    indexes = intvec(0, n-1);
    
    rotate = YES;
    if (rotateNT != NULL) {
        rotate = *rotateNT;
    }
    
    if (rotate == YES && [solution normalTangentialNOFNodes] > 0) {
        dim = [model dimension];
        memset( indexes, 0.0, (n*sizeof(indexes)) );
        for (i=0; i<element->Type.NumberOfNodes; i++) {
            indexes[i] = [solution boundaryReorder:element->NodeIndexes[i]];
        }
        [self FEMKernel_rotateMatrix:solution :stiffMatrix :force :n :dim :dofs :indexes];
    }
    
    if ([solution matrixFormat] == MATRIX_CRS) {
        
        crsMatrix = [[FEMMatrixCRS alloc] init];
        [crsMatrix CRS_glueLocalMatrix:solution :stiffMatrix :n :dofs :nodeIndexes];
        
    }
    else if ([solution matrixFormat] == MATRIX_BAND || [solution matrixFormat] == MATRIX_SBAND) {
        
        bandMatrix = [[FEMMatrixBand alloc] init];
        [bandMatrix BAND_glueLocalMatrix:solution :stiffMatrix :n :dofs :nodeIndexes];
        
    }
    
    
    if (bulkUpdate == NULL || *bulkUpdate == NO) {
        for (i=0; i<n; i++) {
            if (nodeIndexes[i] >= 0) {
                for (j=0; j<dofs; j++) {
                    k = dofs * nodeIndexes[i] + j;
                    [solution setMatrixRHS:k :[solution matrixRHS:k]+force[dofs*i+j]];
                }
            }
        }
    } else if (*bulkUpdate == YES ) {
        for (i=0; i<n; i++) {
            if (nodeIndexes[i] >= 0) {
                for (j=0; j<dofs; j++) {
                    k = dofs * nodeIndexes[i] + j;
                    [solution setMatrixBulkRHS:k :[solution matrixBulkRHS:k]+force[dofs*i+j]];
                }
            }
        }
    }
    
}


#pragma mark Public methods...

- (id)init
{
    self = [super init];
    if (self) {
        // Initialization code here.
        indexStore = intvec(0, 511);
    }
    
    return self;
}

- (void)dealloc
{
    free_ivector(indexStore, 0, 511);
}

-(int)isPElement:(Element_t *)element {
    
    if (element->PDefs != NULL) {
        return 1;
    } else {
        return 0;
    }
}

-(int)getElementFamily:(Element_t *)element {
    
    return element->Type.ElementCode / 100;
}

-(int)getElementDofs:(FEMSolution *)solution forElement:(Element_t *)element atIndexes:(int *)indexes {
    
    int nb, i, j, k, edofs, fdofs, faceDofs, edgeDofs, bubbleDofs, ind;
    BOOL gb;
    Element_t *parent;
    FEMMesh *mesh;
    
    mesh = [solution returnPointerToMesh];
    
    nb = 0;
    
    if ([solution solutionInfoForKey:@"Discontinuous galerkin"] != nil) {
        if ([[solution solutionInfoForKey:@"Discontinuous galerkin"] boolValue] == YES) {
            
            for (i=0; i<element->DGDOFs; i++) {
                indexes[nb] = element->DGIndexes[i];
                nb++;
            }
            
            if (element->BoundaryInfo != NULL) {
                if (element->BoundaryInfo->Left != NULL) {
                    for (i=0; element->BoundaryInfo->Left->DGDOFs; i++) {
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
    
    if ([solution defDofs:0] > 0) {
        for (i=0; i>element->NDOFs; i++) {
            indexes[nb] = element->NodeIndexes[i];
            nb++;
        }
    }
    if (all_in_range([solution returnPointerToDefDofs], '<', 0, 1, [solution sizeOfDefDofs]) == 1) return nb;
    
    faceDofs = mesh.maxFaceDofs;
    edgeDofs = mesh.maxEdgeDofs;
    bubbleDofs = mesh.maxBdofs;
    
    if (element->EdgeIndexes != NULL) {
        for (j=0; j<element->Type.NumberOfEdges; j++) {
            edofs = [mesh Edges_BDOFs:element->EdgeIndexes[j]];
            for (i=0; i<edofs; i++) {
                indexes[nb] = edgeDofs*(element->EdgeIndexes[j]-1) + i + mesh.numberOfNodes;
                nb++;
            }
        }
    }
    
    if (element->FaceIndexes != NULL) {
        for (j=0; j<element->Type.NumberOfFaces; j++) {
            fdofs = [mesh Faces_BDOFs:element->FaceIndexes[j]];
            for (i=0; i<fdofs; i++) {
                indexes[nb] = faceDofs*(element->FaceIndexes[j]-1) + i + mesh.numberOfNodes + 
                edgeDofs*mesh.numberOfEdges;
                nb++;
            }
        }
    }
    
    if ([solution solutionInfoForKey:@"Bubbles in global system"] != nil) {
        gb = [[solution solutionInfoForKey:@"Bubbles in global system"] boolValue];
    } else {
        gb = YES;
    }
    
    if ( element->BoundaryInfo != NULL && (element->EdgeIndexes == NULL || element->FaceIndexes == NULL) ) {
        
        parent = element->BoundaryInfo->Left;
        if (parent == NULL) parent = element->BoundaryInfo->Right;
        if (parent == NULL) return nb;
        
        if ([self getElementFamily:element] == 2) {
            if (parent->EdgeIndexes != NULL) {
                if ([self isPElement:element] == 1) {
                    ind = element->PDefs->LocalNumber;
                } else {
                    for (ind=0; ind<parent->Type.NumberOfEdges; ind++) {
                        k = 0;
                        for (i=0; i<[mesh Edges_Type_NumberOfNodes:parent->EdgeIndexes[ind]]; i++) {
                            for (j=0; element->Type.NumberOfNodes; j++) {
                                if ([mesh Edges_NodeIndexes:parent->EdgeIndexes[ind] :i] == element->NodeIndexes[j]) k++;
                            }
                        }
                        if (k == element->Type.NumberOfNodes) break;
                    }
                }
                
                edofs = element->BDOFs;
                for (i=0; i<edofs; i++) {
                    indexes[nb] = edgeDofs*(parent->EdgeIndexes[ind]-1) + i + mesh.numberOfNodes; 
                    nb++;
                }
            }
        } 
        else if ([self getElementFamily:element] == 3 || [self getElementFamily:element] == 4) {
            if (parent->FaceIndexes != NULL) {
                if ([self isPElement:element] == 1) {
                    ind = element->PDefs->LocalNumber;
                } else {
                    for (ind=0; ind<parent->Type.NumberOfFaces; ind++) {
                        k = 0;
                        for (i=0; i<[mesh Faces_Type_NumberOfNodes:parent->FaceIndexes[ind]]; i++) {
                            for (j=0; j<element->Type.NumberOfNodes; j++) {
                                if ([mesh Faces_NodeIndexes:parent->FaceIndexes[ind] :i] == element->NodeIndexes[j]) k++;
                            }
                        }
                        if (k == element->Type.NumberOfNodes) break;
                    }
                }
                
                fdofs = element->BDOFs;
                for (i=0; i<fdofs; i++) {
                    indexes[nb] = faceDofs*(parent->FaceIndexes[ind]-1) + i + mesh.numberOfNodes
                    + edgeDofs*mesh.numberOfEdges;
                    nb++;
                }
            }
        }
    } else if (gb == YES) {
        if (element->BubbleIndexes != NULL) {
            for (i=0; i<element->BDOFs; i++) {
                indexes[nb] = faceDofs*mesh.numberOfFaces + mesh.numberOfNodes + edgeDofs*mesh.numberOfEdges + element->BubbleIndexes[i];
                nb++;
            }
        }
        
    }
    
    mesh = NULL;
    
    return nb;
    
}

-(int)getNumberOfNodesForElement:(Element_t *)element {
    
    return element->Type.NumberOfNodes;
}

-(int)getBoundaryConditionID:(FEMModel *)model forElement:(Element_t *)element {
    
    int i, bc_id;
    NSArray *boundaryConditions;
    FEMBoundaryCondition *boundaryConditionAtId;
    
    boundaryConditions = [model returnBoundaryConditions];
    
    for (i=0; i<[model numberOfBoundaryConditions]; i++) {
        bc_id++;
        @autoreleasepool {
            boundaryConditionAtId = [boundaryConditions objectAtIndex:i];
            if (element->BoundaryInfo->Constraint == [boundaryConditionAtId tag]) break;
        }
    }
    if (bc_id > [model numberOfBoundaryConditions]) bc_id = 0;
    
    return bc_id;
}

-(NSArray *)getBoundaryCondition:(FEMModel *)model forElement:(Element_t *)element {
    
    int bc_id;
    
    NSArray *boundaryConditions;
    FEMBoundaryCondition *boundaryConditionAtId;
    
    boundaryConditions = [model returnBoundaryConditions];
    // Returns a boundaryCondition object at index bc_id
    bc_id = [self getBoundaryConditionID:model forElement:element];
    boundaryConditionAtId = [boundaryConditions objectAtIndex:bc_id-1];
    
    return [boundaryConditionAtId returnValuesList];
    
}

#pragma mark First order time

-(void)defaultFirstOrderTime:(FEMModel *)model inSolution:(FEMSolution *)solution forElement:(Element_t *)element realMass:(double **)mass realStiff:(double **)stiff realForce:(double *)force stiffRows:(int *)rows stiffCols:(int *)cols { // rows and cols for stiff matrix
    
    int i, n;
    double dt;
    int *perm;
    
    dt = [solution dt];
    n = [self getElementDofs:solution forElement:element atIndexes:indexStore];
    
    perm = intvec(0, n-1);
    for (i=0; i<n; i++) {
        perm[i] = [solution variablePerm:indexStore[i]];
    }
    
    [self FEMKernel_addFirstOrderTime:model :solution :element :mass :stiff :force :dt :n :[solution variableDofs] :perm :rows :cols];
    
    free_ivector(perm, 0, n-1);
    
}

-(void)defaultFirstOrderTime:(FEMModel *)model inSolution:(FEMSolution *)solution forElement:(Element_t *)element complexMass:(_Complex double **)cmass complexStiff:(_Complex double **)cstiff complexForce:(_Complex double *)cforce stiffRows:(int *)rows stiffCols:(int *)cols {
    
    int i, j, n, dofs;
    double dt;
    double **mass, **stiff, *force;
    int *perm;
    
    dt = [solution dt];
    dofs = [solution variableDofs];
    n = [self getElementDofs:solution forElement:element atIndexes:indexStore];
    
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
    
    perm = intvec(0, n-1);
    for (i=0; i<n; i++) {
        perm[i] = [solution variablePerm:indexStore[i]];
    }

    [self FEMKernel_addFirstOrderTime:model :solution :element :mass :stiff :force :dt :n :[solution variableDofs] :perm :rows :cols];
    
    for (i=0; i<n*dofs/2; i++) {
        cforce[i] = force[2*i] + force[2*i+1]*I;
        for (j=0; j<n*dofs/2; j++) {
            cmass[i][j] = mass[2*i][2*j] + (-mass[2*i][2*j+1]*I);
            cstiff[i][j] = cstiff[2*i][2*j] + (-cstiff[2*i][2*j+1]*I);
        }
    }
    
    free_dmatrix(mass, 0, (n*dofs)-1, 0, (n*dofs)-1);
    free_dmatrix(stiff, 0, (n*dofs)-1, 0, (n*dofs)-1);
    free_dvector(force, 0, (n*dofs)-1);
    
    free_ivector(perm, 0, n-1);
    
}

#pragma mark Update equations

-(void)defaultUpdateEquations:(FEMModel *)model inSolution:(FEMSolution *)solution forElement:(Element_t *)element realStiff:(double **)stiff realForce:(double *)force stiffRows:(int *)rows stiffCols:(int *)cols requestBulkUpdate:(BOOL *)bulkUpdate {
    
    int i, n;
    int *perm;
    BOOL rotateNT, bupd;
    double *saveValues;
    
    n = [self getElementDofs:solution forElement:element atIndexes:indexStore];
    
    perm = intvec(0, n-1);
    
    for (i=0; i<n; i++) {
        perm[i] = [solution variablePerm:indexStore[i]];
    }
    
    [self FEMKernel_updateGlobalEquations:model :solution :element :stiff :force :n :[solution variableDofs] :perm :rows :cols :NULL: NULL];
    
    bupd = FALSE;
    if (bulkUpdate != NULL) {
        bupd = *bulkUpdate;
        if (bupd == FALSE) return;
    } else {
        if (element->BoundaryInfo != NULL) return;
    }
    if ([solution solutionInfoForKey:@"Calculate loads"] != nil) {
        bupd = (bupd == YES || [[solution solutionInfoForKey:@"Calculate loads"] boolValue] == YES) ? YES : NO;
    }
    
    if (bupd ==YES) {
        
        if ([solution matrixReturnPointerToBulkRHS] == NULL) {
            [solution allocateMatrixBulkRHS:[solution matrixSizeOfRHS]];
            [solution setMatrixSizeOfBulkRHS:[solution matrixSizeOfRHS]];
            for (i=0; i<[solution matrixSizeOfBulkRHS]; i++) {
                [solution setMatrixBulkRHS:i :0.0];
            }
        }
        
        if ([solution matrixReturnPointerToBulkValues] == NULL) {
            [solution allocateMatrixBulkValues:[solution matrixSizeOfValues]];
            [solution setMatrixSizeOfBulkValues:[solution matrixSizeOfValues]];
            for (i=0; i<[solution matrixSizeOfBulkValues]; i++) {
                [solution setMatrixBulkValues:i :0.0];
            }
        }
        saveValues = doublevec(0, [solution matrixSizeOfValues]-1);
        for (i=0; i<[solution matrixSizeOfValues]; i++) {
            saveValues[i] = [solution matrixValues:i];
        }
        for (i=0; i<[solution matrixSizeOfValues]; i++) {
            [solution setMatrixValues:i :[solution matrixBulkValues:i]];
        }
        rotateNT = NO;
        [self FEMKernel_updateGlobalEquations:model :solution :element :stiff :force :n :[solution variableDofs] :perm :rows :cols :&rotateNT :&bupd];
        for (i=0; i<[solution matrixSizeOfValues]; i++) {
            [solution setMatrixValues:i :saveValues[i]];
        }
        free_dvector(saveValues, 0, [solution matrixSizeOfValues]-1);
        
    }
    
    free_ivector(perm, 0, n-1);
    
}

-(void)defaultUpdateEquations:(FEMModel *)model inSolution:(FEMSolution *)solution forElement:(Element_t *)element complexStiff:(double **)cstiff complexForce:(double *)cforce stiffRows:(int *)rows stiffCols:(int *)cols requestBulkUpdate:(BOOL *)bulkUpdate {
    
    int i, j, n, dofs;
    int *perm;
    double **stiff, *force;
    BOOL rotateNT, bupd;
    double *saveValues;
    
    n = [self getElementDofs:solution forElement:element atIndexes:indexStore];
    
    perm = intvec(0, n-1);
    
    for (i=0; i<n; i++) {
        perm[i] = [solution variablePerm:indexStore[i]];
    }
    
    dofs = [solution variableDofs];
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

    [self FEMKernel_updateGlobalEquations:model :solution :element :stiff :force :n :[solution variableDofs] :perm :rows :cols :NULL: NULL];
    
    bupd = FALSE;
    if (bulkUpdate != NULL) {
        bupd = *bulkUpdate;
        if (bupd == FALSE) return;
    } else {
        if (element->BoundaryInfo != NULL) return;
    }
    if ([solution solutionInfoForKey:@"Calculate loads"] != nil) {
        bupd = (bupd == YES || [[solution solutionInfoForKey:@"Calculate loads"] boolValue] == YES) ? YES : NO;
    }
    
    if (bupd ==YES) {
        
        if ([solution matrixReturnPointerToBulkRHS] == NULL) {
            [solution allocateMatrixBulkRHS:[solution matrixSizeOfRHS]];
            [solution setMatrixSizeOfBulkRHS:[solution matrixSizeOfRHS]];
            for (i=0; i<[solution matrixSizeOfBulkRHS]; i++) {
                [solution setMatrixBulkRHS:i :0.0];
            }
        }
        
        if ([solution matrixReturnPointerToBulkValues] == NULL) {
            [solution allocateMatrixBulkValues:[solution matrixSizeOfValues]];
            [solution setMatrixSizeOfBulkValues:[solution matrixSizeOfValues]];
            for (i=0; i<[solution matrixSizeOfBulkValues]; i++) {
                [solution setMatrixBulkValues:i :0.0];
            }
        }
        saveValues = doublevec(0, [solution matrixSizeOfValues]-1);
        for (i=0; i<[solution matrixSizeOfValues]; i++) {
            saveValues[i] = [solution matrixValues:i];
        }
        for (i=0; i<[solution matrixSizeOfValues]; i++) {
            [solution setMatrixValues:i :[solution matrixBulkValues:i]];
        }
        rotateNT = NO;
        [self FEMKernel_updateGlobalEquations:model :solution :element :stiff :force :n :[solution variableDofs] :perm :rows :cols :&rotateNT :&bupd];
        for (i=0; i<[solution matrixSizeOfValues]; i++) {
            [solution setMatrixValues:i :saveValues[i]];
        }
        free_dvector(saveValues, 0, [solution matrixSizeOfValues]-1);
        
    }
    
    free_dmatrix(stiff, 0, (n*dofs)-1, 0, (n*dofs)-1);
    free_dvector(force, 0, (n*dofs)-1);    
    free_ivector(perm, 0, n-1);
    
}


-(double)findSolution:(FEMSolution *)solution {
    
    double norm;
    
    if ([solution solutionInfoForKey:@"Linear system solver disabled"] != nil) {
        if ([[solution solutionInfoForKey:@"Linear system solver disabled"] boolValue] == YES) return 0.0;
    }
    
    
    [self FEMKernel_solveSystem:solution];
    norm = [solution variableNorm];
    
    return norm;
}

-(double)stopc:(FEMSolution *)solution :(double *)x :(double *)b :(double *)r :(int *)ipar {
/*********************************************************************************
    We don't use the backward error estimate e = ||Ax-b||/||A|| ||X|| + ||b||
    as stoping criterion
*********************************************************************************/ 
    int i, n;
    double err, *res;
    double sum1, sum2, sum3, sum4;
    FEMPrecondition *preconditioning;
    
    n = ipar[2];
    res = doublevec(0, ipar[2]-1);
    
    preconditioning = [[FEMPrecondition alloc] init];
    [preconditioning CRS_MatrixVectorMultiply:solution :x :res];
    
    for (i=0; i<ipar[2]; i++) {
        res[i] = res[i] - b[i];
    }
    
    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;
    sum4 = 0.0;
    
    for (i=0; i<n; i++) {
        sum1 = sum1 + pow(res[i], 2.0);
    }
    for (i=0; i<n; i++) {
        sum2 = sum2 + pow([solution matrixValues:i], 2.0);
    }
    for (i=0; i<n; i++) {
        sum3 = sum3 + pow(x[i], 2.0);
    }
    for (i=0; i<n; i++) {
        sum4 = sum4 + pow(b[i], 2.0);
    }
    err = sqrt(sum1) / ( sqrt(sum2) * sqrt(sum3) + sqrt(sum4) );
    
    free_dvector(res, 0, ipar[2]-1);
    
    return err;
    
}


#pragma mark memory allocation utility methods

-(BOOL)allocateIntegersVector:(int *)Vector from:(long)nl to:(long)nh {
    
    Vector = intvec(nl, nh);
    
    // Return a bool NO if allocation failes but don't terminate. 
    // The user calling this method will have to decide what to do. 
    // Do the same for all allcations methods
    if (Vector == NULL) {
        return NO;
    }
    else {
        return YES;
    }
}

-(BOOL)allocateFloatsVector:(float *)Vector from:(long)nl to:(long)nh {
    
    Vector = floatvec(nl, nh);
    if (Vector == NULL) {
        return NO;
    } 
    else {
        return YES;
    }
}

-(BOOL)allocateDoublesVector:(double *)Vector from:(long)nl to:(long)nh {
    
    Vector = doublevec(nl, nh);
    if (Vector == NULL) {
        return NO;
    } 
    else {
        return YES;
    }
}

-(BOOL)allocateIntegersMatrix:(int **)Matrix fromRow:(long)nrl toRow:(long)nrh fromColumn:(long)ncl toColumn:(long)nch {
    
    Matrix = intmatrix(nrl, nrh, ncl, nch);
    if (Matrix == NULL) {
        return NO;
    } 
    else {
        return YES;
    }
}

-(BOOL)allocateFloatsMatrix:(float **)Matrix fromRow:(long)nrl toRow:(long)nrh fromColumn:(long)ncl toColumn:(long)nch {
    
    Matrix = floatmatrix(nrl, nrh, ncl, nch);
    if (Matrix == NULL) {
        return NO;
    } 
    else {
        return YES;
    }
}

-(BOOL)allocateDoublesMatrix:(double **)Matrix fromRow:(long)nrl toRow:(long)nrh fromColumn:(long)ncl toColumn:(long)nch {
    
    Matrix = doublematrix(nrl, nrh, ncl, nch);
    if (Matrix == NULL) {
        return NO;
    } 
    else {
        return YES;
    }
}

-(BOOL)allocateIntegersTensor:(int ***)Tensor fromRow:(long)nrl toRow:(long)nrh fromColumn:(long)ncl toColumn:(long)nch fromDepth:(long)ndl toDepth:(long)ndh {
    
    Tensor = i3tensor(nrl, nrh, ncl, nch, ndl, ndh);
    if (Tensor == NULL) {
        return NO;
    } 
    else {
        return YES;
    }
}

-(BOOL)allocateFloatsTensor:(float ***)Tensor fromRow:(long)nrl toRow:(long)nrh fromColumn:(long)ncl toColumn:(long)nch fromDepth:(long)ndl toDepth:(long)ndh {
    
    Tensor = f3tensor(nrl, nrh, ncl, nch, ndl, ndh);
    if (Tensor == NULL) {
        return NO;
    } 
    else {
        return YES;
    }
}

-(BOOL)allocateDoublesTensor:(double ***)Tensor fromRow:(long)nrl toRow:(long)nrh fromColumn:(long)ncl toColumn:(long)nch fromDepth:(long)ndl toDepth:(long)ndh {
    
    Tensor = d3tensor(nrl, nrh, ncl, nch, ndl, ndh);
    if (Tensor == NULL) {
        return NO;
    } 
    else {
        return YES;
    }
}

-(void)freeIntegersVector:(int *)Vector from:(long)nl to:(long)nh {
    
    free_ivector(Vector, nl, nh);
}

-(void)freeFloatsVector:(float *)Vector from:(long)nl to:(long)nh {
    
    free_fvector(Vector, nl, nh);
}

-(void)freeDoublesVector:(double *)Vector from:(long)nl to:(long)nh {
    
    free_dvector(Vector, nl, nh);
}

-(void)freeIntegersMatrix:(int **)Matrix fromRow:(long)nrl toRow:(long)nrh fromColumn:(long)ncl toColumn:(long)nch {
    
    free_imatrix(Matrix, nrl, nrh, ncl, nch);
}

-(void)freeFloatsMatrix:(float **)Matrix fromRow:(long)nrl toRow:(long)nrh fromColumn:(long)ncl toColumn:(long)nch {
    
    free_fmatrix(Matrix, nrl, nrh, ncl, nch);
}

-(void)freeDoublesMatrix:(double **)Matrix fromRow:(long)nrl toRow:(long)nrh fromColumn:(long)ncl toColumn:(long)nch {
    
    free_dmatrix(Matrix, nrl, nrh, ncl, nch);
}

-(void)freeIntergersTensor:(int ***)Tensor fromRow:(long)nrl toRow:(long)nrh fromColumn:(long)ncl toColumn:(long)nch fromDepth:(long)ndl toDepth:(long)ndh {
    
    free_i3tensor(Tensor, nrl, nrh, ncl, nch, ndl, ndh);
}

-(void)freeFloatsTensor:(float ***)Tensor fromRow:(long)nrl toRow:(long)nrh fromColumn:(long)ncl toColumn:(long)nch fromDepth:(long)ndl toDepth:(long)ndh {
    
    free_f3tensor(Tensor, nrl, nrh, ncl, nch, ndl, ndh);
}

-(void)freeDoublesTensor:(double ***)Tensor fromRow:(long)nrl toRow:(long)nrh fromColumn:(long)ncl toColumn:(long)nch fromDepth:(long)ndl toDepth:(long)ndh {
    
    free_d3tensor(Tensor, nrl, nrh, ncl, nch, ndl, ndh);
}

@end
