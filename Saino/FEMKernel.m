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
#import "FEMListUtilities.h"
#import "FEMBodyForce.h"
#import "FEMBoundaryCondition.h"
#import "FEMSimulation.h"

#import "memory.h"
#import "Utils.h"
#import "Numerics.h"
#import "GaussIntegration.h"

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

// Solving linear systems and compute norms
-(void)iterCall:(int)iterType: (FEMSolution *)solution: (int *)ipar: (double *)dpar: (double **)work: (SEL)pcondlMethod: (SEL)pcondrMethod: (SEL)matvecMethod: (SEL)mstopMethod;
-(void)FEMKernel_iterSolver:(FEMSolution *)solution;
-(void)FEMKernel_rotateNTSystem:(FEMSolution *)solution: (double *)vec: (int)nodeNumber;
-(void)FEMKernel_backRotateNTSystem:(FEMSolution *)solution;
-(double)FEMKernel_computeNorm:(FEMSolution *)solution: (int)n: (double *)values;
-(void)FEMKernel_computeChange:(FEMSolution *)solution model:(FEMModel *)aModel isSteadyState:(BOOL)steadyState nsize:(int*)nsize values:(double *)values values0:(double *)values0;
-(void)FEMKernel_solveLinearSystem:(FEMSolution *)solution model:(FEMModel *)aModel;
-(void)FEMKernel_solveSystem:(FEMSolution *)solution model:(FEMModel *)aModel;

// Check passive element
-(BOOL)FEMKernel_checkPassiveElement:(FEMModel *)model: (FEMSolution *)solution: (Element_t *)element;

// Rotate matrix
-(void)FEMKernel_rotateMatrix:(FEMSolution *)solution: (double **)matrix: (double *)vector: (int)n: (int)dim: (int)dofs: (int *)nodeIndexes;

// Update global force
-(void)FEMKernel_updateGlobalForce:(FEMModel *)model: (FEMSolution *)solution: (Element_t *)element: (double *)forceVector: (double *)localForce: (int)n: (int)dofs: (int *)nodeIndexes: (BOOL *)rotateNT;

// Manipulate matrix coefficients for time dependent simulations
-(void)FEMKernel_addFirstOrderTime:(FEMModel *)model: (FEMSolution *)solution: (Element_t *)element: (double **)massMatrix: (double **)stiffMatrix: (double *)force: (double)dt: (int)n: (int)dofs: (int *)nodeIndexes: (int *)rows: (int *)cols;

// Update global equations
-(void)FEMKernel_updateGlobalEquations:(FEMModel *)model: (FEMSolution *)solution: (Element_t *)element: (double **)stiffMatrix: (double *)force: (int)n: (int)dofs: (int *)nodeIndexes: (int *)rows: (int *)cols: (BOOL *)rotateNT: (BOOL *)bulkUpdate;

// Loads
-(void)FEMKernel_setElementLoads:(FEMModel *)model: (FEMSolution *)solution: (Element_t *)element: (NSArray *)values: (NSMutableString *)name: (int*)indexes: (BOOL *)doneLoad: (int)n: (int)dof: (int)ndofs;
-(void)FEMKernel_setPointLoads:(FEMModel *)model: (FEMSolution *)solution: (Element_t *)element: (NSArray *)values: (NSMutableString *)name: (int*)indexes: (int)n: (int)dof: (int)ndofs;

// Element and point values
-(void)FEMKernel_setElementValues:(FEMModel *)model inSolution:(FEMSolution *)solution forElementNumber:(int)elno numberOfNodes:(int)n atIndexes:(int *)indexes forValues:(NSArray *)values variableName:(NSMutableString *)name orderOfDofs:(int)dof activeCondition:(BOOL)conditional conditionName:(NSString *)condName permutationOffset:(int)offset;
-(void)FEMKernel_setPointValues:(FEMModel *)model inSolution:(FEMSolution *)solution numberofNodes:(int)n atIndexes:(int *)indexes forValues:(NSArray *)values variableName:(NSMutableString *)name orderOfDofs:(int)dof activeCondition:(BOOL)conditional conditionName:(NSString *)condName permutationOffset:(int)offset;

//Periodic
-(void)FEMKernel_setPeriodicBoundariesPass1:(FEMModel *)model: (FEMSolution *)solution: (NSMutableString *)name: (int)dof: (int)this: (BOOL *)done;
-(void)FEMKernel_setPeriodicBoundariesPass2:(FEMModel *)model: (FEMSolution *)solution: (NSMutableString *)name: (int)dof: (int)this: (BOOL *)done;


@end

@implementation FEMKernel {
    
    int *_indexStore;
    double **_kernStiff, *_kernWork;           // kernStiff(maxElementDofs)(maxElementDofs), kernWork(maxElementDofs)
    int *_g_Ind, *_l_Ind;                      // g_Ind(maxElementDofs), l_Ind(maxElementDofs)
    int _sizeIndexStore;
    int _size1kernStiff;
    int _size2kernStiff;
    int _sizekernWork;
    int _size_g_Ind;
    int _size_l_Ind;
    
    int **_lineEM;
    int **_triangleEM;
    int **_quadEM;
    int **_tetraEM;
    int **_prismEM;
    int **_wedgeEM;
    int **_brickEM;
    
    BOOL _initialized[8];
}

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
    variableArraysContainer *varContainers;
    BOOL abortNotConverged;
    SEL pcondSelector=0, pcondrSelector=0, mvSelector=0;
    
    n = solution.matrix.numberOfRows;
    ipar = intvec(0, 49);
    dpar = doublevec(0, 49);
    
    varContainers = solution.variable.getContainers;
    
    for (i=0; i<50; i++) {
        ipar[i] = 0;
        dpar[i] = 0.0;
    }
    
    str = [NSString stringWithString:(solution.solutionInfo)[@"Linear system iteration method"]];
    
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
        if ((solution.solutionInfo)[@"Linear system GMRES restart"] != nil) {
            
            ipar[14] = [(solution.solutionInfo)[@"Linear system GMRES restart"] intValue];
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
        if ((solution.solutionInfo)[@"SGS over relaxation factor"] != nil) {
            dpar[2] = [(solution.solutionInfo)[@"SGS over relaxation factor"] doubleValue];
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
        if ((solution.solutionInfo)[@"Linear system GCR restart"] != nil) {
            ipar[16] = [(solution.solutionInfo)[@"Linear system GCR restart"] intValue];
        } else {
            if ((solution.solutionInfo)[@"Linear system maximum iterations"] != nil) {
                ipar[16] = [(solution.solutionInfo)[@"Linear system maximum iterations"] intValue];
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
        if ((solution.solutionInfo)[@"Bi-CGSTAB(l) polynomial degree"] != nil) {
            ipar[15] = [(solution.solutionInfo)[@"Bi-CGSTAB(l) polynomial degree"] intValue];
            if (ipar[15] < 2) errorfunct("FEMKernel_iterSolver", "Polynomial degree for Bi-CGSTAB(l) should be equal or greater than 2");
        } else {
            ipar[15] = 2;
        }
    }
    
    ipar[11] = 1;
    ipar[2] = n;
    
    if ((solution.solutionInfo)[@"Linear system residual output"] != nil) {
        ipar[4] = [(solution.solutionInfo)[@"Linear system residual output"] intValue];
    } else {
        ipar[4] = 1;
    }
    
    if ((solution.solutionInfo)[@"Linear system maximum iterations"] != nil) {
        ipar[9] = [(solution.solutionInfo)[@"Linear system maximum iterations"] intValue];
    } else {
        ipar[9] = 1;
    }
    
    work = doublematrix(0, n-1, 0, wsize-1);
    if (work == NULL) {
        errorfunct("FEMKernel_iterSolver", "Memory allocation error");
    }
    
    if (all(varContainers->Values, '=', 0.0, varContainers->sizeValues) == 1) {
        for (i=0; i<varContainers->sizeValues; i++) {
            varContainers->Values[i] = 1.0e-8;
        }
    }
    ipar[13] = 1;
    
    if ((solution.solutionInfo)[@"Linear system convergence tolerance"] != nil) {
        dpar[0] = [(solution.solutionInfo)[@"Linear system convergence tolerance"] doubleValue];
    }
    
    if ((solution.solutionInfo)[@"Linear system divergence tolerance"] != nil) {
        dpar[1] = [(solution.solutionInfo)[@"Linear system divergence tolerance"] doubleValue];
    } else {
        dpar[1] = HUGE_VAL;
    }
    
    if ((solution.solutionInfo)[@"Linear system preconditioning"] != nil) {
        str = [NSString stringWithString:(solution.solutionInfo)[@"Linear system preconditioning"]];
    } else {
        str = @"none";
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
        if ((solution.solutionInfo)[@"Linear system ILUT tolerance"] != nil) {
            ilut_tol = [(solution.solutionInfo)[@"Linear system ILUT tolerance"] doubleValue];
        } else {
            errorfunct("FEMKernel_iterSolver", "Linear system ILUT tolerance not found.");
        }
        pCondType = PRECOND_ILUT;
    }
    else if ([str isEqualToString:@"ILU"] == YES) 
    {
        if ((solution.solutionInfo)[@"Linear system ILU order"] != nil) {
            ilun = [(solution.solutionInfo)[@"Linear system ILU order"] intValue];
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
        if (solution.variable.dofs == 1) {
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
    
    if ((solution.solutionInfo)[@"No precondition recompute"] != nil) {
        if ([(solution.solutionInfo)[@"No precondition recompute"] boolValue] == NO) {
            // TODO: Implement the code if we choose not to compute the preconditioner
            // for each method call.
        }
    }
    
    solution.matrix.solveCount = solution.matrix.solveCount + 1;
    
    if ((solution.solutionInfo)[@"Linear system abort not converged"] != nil) {
        abortNotConverged = [(solution.solutionInfo)[@"Linear system abort not converged"] boolValue];
    } else {
        abortNotConverged = YES;
    }
    
    // Get the selector for the matrix-vector multiplication method we want to use
    if (solution.matrix.isComplexMatrix == YES) {
        mvSelector = @selector(CRS_ComplexMatrixVectorProd::::);
    } else {
        mvSelector= @selector(CRS_MatrixVectorProd::::);
    }
    
    // Get the selector for the preconditioning method we want to use
    if (pCondType == PRECOND_NONE) {
        
        pcondSelector = @selector(CRS_pcond_dummy::::);
    } 
    else if (pCondType == PRECOND_DIAGONAL) {
        if (solution.matrix.isComplexMatrix == YES) {
            pcondSelector = @selector(CRS_ComplexDiagPrecondition::::);
        } else {
            pcondSelector = @selector(CRS_DiagPrecondition::::);
        }
    } 
    else if (pCondType == PRECOND_ILUN == pCondType == PRECOND_ILUT || pCondType == PRECOND_BILUN) {
        if (solution.matrix.isComplexMatrix == YES) {
            pcondSelector = @selector(CRS_ComplexLUPrecondition::::);
        } else {
            pcondSelector = @selector(CRS_LUPrecondition::::);
        }
    }
    
    if (solution.matrix.isComplexMatrix == NO) {
        
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
    
    if (solution.matrix.isComplexMatrix == YES) ipar[2] = ipar[2] * 2;
    
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
    varContainers = NULL;
}

-(void)FEMKernel_rotateNTSystem:(FEMSolution *)solution: (double *)vec: (int)nodeNumber {
    
    int i, k, dim;
    double bu, bv, bw, **rm;
    solutionArraysContainer *solContainers;
    
    if ([solution normalTangentialNOFNodes] <= 0) return;
    
    solContainers = solution.getContainers;
    
    dim = [solution coordinateSystemDimension];
    
    k = solContainers->boundaryReorder[nodeNumber];
    if (k < 0) return;
    
    if (dim < 3) {
        bu = vec[0];
        bv = vec[1];
        vec[0] = solContainers->boundaryNormals[k][0]*bu + solContainers->boundaryNormals[k][1]*bv;
        vec[1] = -solContainers->boundaryNormals[k][1]*bu + solContainers->boundaryNormals[k][0]*bv;
    } else {
        rm = doublematrix(0, 2, 0, 2);
        
        bu = vec[0];
        bv = vec[1];
        bw = vec[2];
        
        for (i=0; i<2; i++) {
            rm[i][0] = solContainers->boundaryNormals[k][i];
            rm[i][1] = solContainers->boundaryTangent1[k][i];
            rm[i][2] = solContainers->boundaryTangent2[k][i];
        }
        
        vec[0] = rm[0][0]*bu + rm[1][0]*bv + rm[2][0]*bw;
        vec[1] = rm[0][1]*bu + rm[1][1]*bv + rm[2][1]*bw;
        vec[2] = rm[0][2]*bu + rm[1][2]*bv + rm[2][2]*bw;

        free_dmatrix(rm, 0, 2, 0, 2);
    }
    
    solContainers = NULL;
}

-(void)FEMKernel_backRotateNTSystem:(FEMSolution *)solution {
    
    int i, j, k, l, dim, ndofs;
    double bu, bv, bw, **rm;
    variableArraysContainer *varContainers;
    solutionArraysContainer *solContainers;
    
    if ([solution normalTangentialNOFNodes] <= 0) return;
    
    varContainers = solution.variable.getContainers;
    solContainers = solution.getContainers;
    
    dim = [solution coordinateSystemDimension];
    ndofs = solution.variable.dofs;
    
    for (i=0; i<solContainers->sizeBoundaryReorder; i++) {
        k = solContainers->boundaryReorder[i];
        if (k < 0) continue;
        j = varContainers->Perm[i];
        if (j < 0) continue;
        
        if (dim < 3) {
            
            bu = varContainers->Values[(ndofs*(j-1)+1)];
            bv = varContainers->Values[(ndofs*(j-1)+2)];
            
            varContainers->Values[(ndofs*(j-1)+1)] = solContainers->boundaryNormals[k][0]*bu - solContainers->boundaryNormals[k][1]*bv;
            varContainers->Values[(ndofs*(j-1)+2)] = solContainers->boundaryNormals[k][1]*bu + solContainers->boundaryNormals[k][0]*bv;
        } else {
            
            rm = doublematrix(0, 2, 0, 2);
            
            bu = varContainers->Values[(ndofs*(j-1)+1)];
            bv = varContainers->Values[(ndofs*(j-1)+2)];
            bw = varContainers->Values[(ndofs*(j-1)+3)];
            
            for (l=0; l<2; l++) {
                rm[0][l] = solContainers->boundaryNormals[k][l];
                rm[1][l] = solContainers->boundaryTangent1[k][l];
                rm[2][l] = solContainers->boundaryTangent2[k][l];
            }

            varContainers->Values[(ndofs*(j-1)+1)] = rm[0][0]*bu + rm[1][0]*bv + rm[2][0]*bw;
            varContainers->Values[(ndofs*(j-1)+2)] = rm[0][1]*bu + rm[1][1]*bv + rm[2][1]*bw;
            varContainers->Values[(ndofs*(j-1)+3)] = rm[0][2]*bu + rm[1][2]*bv + rm[2][2]*bw;
            
            free_dmatrix(rm, 0, 2, 0, 2);
            
        }
    }
    
    varContainers = NULL;
    solContainers = NULL;
}

-(double)FEMKernel_computeNorm:(FEMSolution *)solution: (int)n: (double *)values {
    
    int normDim, normDofs, dofs, i, j, k, l, totn;
    double norm, nscale, sum;
    double *x, *buffer;
    FEMParallelMPI *parallelUtil;
    variableArraysContainer *varContainers;
    
    parallelUtil = [[FEMParallelMPI alloc] init];
    varContainers = solution.variable.getContainers;
    
    if (values != NULL) {
        x = values;
    } else {
        x = varContainers->Values;
    }
    
    if ((solution.solutionInfo)[@"Non linear system norm degree"] != nil) {
        normDim = [(solution.solutionInfo)[@"Non linear system norm degree"] intValue];
    } else {
        normDim = 2;
    }
    
    dofs = solution.variable.dofs;
    if ((solution.solutionInfo)[@"Non linear system norm dofs"] != nil) {
        normDofs = [(solution.solutionInfo)[@"Non linear system norm dofs"] intValue];
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
    
    x = NULL;
    varContainers = NULL;
    return norm;
}


-(void)FEMKernel_computeChange:(FEMSolution *)solution model:(FEMModel *)aModel isSteadyState:(BOOL)steadyState nsize:(int*)nsize values:(double *)values values0:(double *)values0 {
    
    NSString *convergenceType, *solverName;
    int i, n, relaxAfter, iterNo;
    double norm, prevNorm, bNorm, change, relaxation, maxNorm, dt, tolerance, eps;
    double *x, *r, *x0;
    BOOL skip, convergenceAbsolute, relax, relaxBefore, stat, doIt;
    FEMVariable *iterV, *timeStepVar, *veloVar;
    NSMutableString *str1;
    NSString *str2;
    matrixArraysContainer *matContainers;
    variableArraysContainer *varContainers, *itervContainers, *containers;
    FEMPrecondition *preconditioning;
    FEMUtilities *utilities;
    
    preconditioning = [[FEMPrecondition alloc] init];
    varContainers = solution.variable.getContainers;
    
    utilities = [[FEMUtilities alloc] init];
    
    relax = NO;
    
    if (steadyState == YES) {
        skip = [(solution.solutionInfo)[@"Skip compute steady state change"] boolValue];
        if (skip == YES) return;
        
        if ((solution.solutionInfo)[@"Steady state convergence measure"] != nil) {
            convergenceType = [NSString stringWithString:(solution.solutionInfo)[@"Steady state convergence measure"]];
        } else {
            convergenceType = @"norm";
        }
        
        if ((solution.solutionInfo)[@"Steady state convergence absolute"] != nil) {
            convergenceAbsolute = [(solution.solutionInfo)[@"Steady state convergence absolute"] boolValue];
        } else if ((solution.solutionInfo)[@"Use absolute norm for convergence"] != nil) {
            convergenceAbsolute = [(solution.solutionInfo)[@"Use absolute norm for convergence"] boolValue];
        }
        
        if ((solution.solutionInfo)[@"Steady state relaxation factor"] != nil) {
            relaxation = [(solution.solutionInfo)[@"Steady state relaxation factor"] doubleValue];
            relax = (relaxation != 1.0) ? YES: NO;
        }
        
        iterV = [utilities getVariableFrom:solution.mesh.variables model:aModel name:@"Coupled iter" onlySearch:NULL maskName:NULL info:&stat];
        itervContainers = iterV.getContainers;
        iterNo = (int)itervContainers->Values[0];
        itervContainers = NULL;
        
        if (relax == YES) {
            if ((solution.solutionInfo)[@"Steady state relaxation after"] != nil) {
                relaxAfter = [(solution.solutionInfo)[@"Steady state relaxation after"] intValue];
                if (relaxAfter >= iterNo ) relax = NO;
            }    
        }
        
        if (relax == YES) {
            if ((solution.solutionInfo)[@"Steady state relaxation before"] != nil) {
                relaxBefore = [(solution.solutionInfo)[@"Steady state relaxation before"] boolValue];
            } else {
                relaxBefore = YES;
            }
        }
    } else {
    
        if ((solution.solutionInfo)[@"Skip compute non linear change"] != nil) {
            skip  = [(solution.solutionInfo)[@"Skip compute non linear change"] boolValue];
        } else {
            skip = NO;
        }
        if (skip == YES) return;
        
        if ((solution.solutionInfo)[@"Non linear system convergence measure"] != nil) {
            convergenceType = [NSString stringWithString:(solution.solutionInfo)[@"Non linear system convergence measure"]];
        } else {
            convergenceType = @"norm";
        }
        
        if ((solution.solutionInfo)[@"Non linear system convergence absolute"] != nil) {
            convergenceAbsolute = [(solution.solutionInfo)[@"Non linear system convergence absolute"] boolValue];
        } else if ((solution.solutionInfo)[@"Use absolute norm for convergence"] != nil) {
            convergenceAbsolute = [(solution.solutionInfo)[@"Use absolute norm for convergence"] boolValue];
        } else {
            convergenceAbsolute = NO;
        }
        
        iterV = [utilities getVariableFrom:solution.mesh.variables model:aModel name:@"Nonlin iter" onlySearch:NULL maskName:NULL info:&stat];
        itervContainers = iterV.getContainers;
        iterNo = (int)itervContainers->Values[0];
        solution.variable.nonLinIter = (int)itervContainers->Values[0];
        itervContainers->Values[0] = itervContainers->Values[0]+1;
        itervContainers = NULL;
        
        if ((solution.solutionInfo)[@"Non linear system relaxation factor"] != nil) {
            relaxation = [(solution.solutionInfo)[@"Non linear system relaxation factor"] doubleValue];
            relax = (relaxation != 1.0) ? YES: NO;
        }
        if (relax == YES) {
            if ((solution.solutionInfo)[@"Non linear system relaxation after"] != nil) {
                relaxAfter = [(solution.solutionInfo)[@"Non linear system relaxation after"] boolValue];
                if (relaxAfter >= solution.variable.nonLinIter) relax = NO;
            }
        }
        
        if (relax == YES) {
            if ((solution.solutionInfo)[@"Non linear system relaxation before"] != nil) {
                relaxBefore = [(solution.solutionInfo)[@"Non linear system relaxation before"] boolValue];
            } else {
                relaxBefore = YES;
            }
        }
        
    }
    
    if (values != NULL) {
        x = values;
    } else {
        x = varContainers->Values;
    }
    
    if (values == NULL) {
        solution.variable.norm = 0.0;
        if (steadyState == YES) {
            solution.variable.steadyChange = 0.0;
        } else {
            solution.variable.nonLinChange = 0.0;
        }
        x = NULL;
        varContainers = NULL;
        return;
    }
    
    if (nsize != NULL && values != NULL) {
        n = *nsize;
    } else if (varContainers->Values != NULL){
        n = varContainers->sizeValues;
    } 
    
    stat = NO;
    if (values0 != NULL) {
        x0 = values0;
        stat = YES;
    } else if (steadyState == YES) {
        if (varContainers->SteadyValues != NULL) {
            x0 = varContainers->SteadyValues;
            stat = YES;
        }
    } else {
        if (varContainers->NonLinValues != NULL) {
            x0 = varContainers->NonLinValues;
            stat = YES;
        }
    }
    
    // TODO: Possibly check here if length mismatch between x0 and x
    //.....
    
    if (relax == YES && relaxBefore == YES) {
        for (i=0; i<n; i++) {
            x[i] = (1-relaxation) * x0[i] + relaxation*x[i];
        }
    }
    
    if (steadyState == YES) {
        prevNorm = solution.variable.prevNorm;
    } else {
        prevNorm = solution.variable.norm;
    }
    
    // Compute norm here
    norm = [self FEMKernel_computeNorm:solution :n :x];
    solution.variable.norm = norm;
    
    // The norm should be bounded in order to reach convergence
    if ((solution.solutionInfo)[@"Non linear system max norm"] != nil) {
        maxNorm = [(solution.solutionInfo)[@"Non linear system max norm"] doubleValue];
    } else {
        maxNorm = HUGE_VAL;
    }
    
    if (isnan(norm) != 0 && norm > maxNorm) {
        warnfunct("FEMKernel_computeChange", "Computed norm:");
        printf("%lf\n", norm);
        errorfunct("FEMKernel_computeChange", "Norm of solution has crossed given bounds.");
    }
    
     matContainers = solution.matrix.getContainers;
    
    if ([convergenceType isEqualToString:@"residual"] == YES) {
        // ------------------------------------------------------------------------------
        // x is solution of A(x0)x = b(x0), thus residual should be real r = b(x)-A(x)x
        // Instead we use r = b(x0)-A(x0)x0 which unfortunately is one step behind
        // ------------------------------------------------------------------------------
        r = doublevec(0, n-1);
        
        [preconditioning CRS_MatrixVectorMultiply:solution :x0 :r];
        for (i=0; i<n; i++) {
            r[i] = r[i] - matContainers->RHS[i];
        }
        change = [self FEMKernel_computeNorm:solution :n :r];
        if (convergenceAbsolute == NO) {
            bNorm = [self FEMKernel_computeNorm:solution :n :matContainers->RHS];
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
        r = doublevec(0, n-1);
        
        [preconditioning CRS_MatrixVectorMultiply:solution :x :r];
        for (i=0; i<n; i++) {
            r[i] = r[i] - matContainers->RHS[i];
        }
        change = [self FEMKernel_computeNorm:solution :n :r];
        if (convergenceAbsolute == NO) {
            bNorm = [self FEMKernel_computeNorm:solution :n :matContainers->RHS];
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
        solution.variable.steadyChange = change;
        if ((solution.solutionInfo)[@"Steady state convergence tolerance"] != nil) {
            tolerance = [(solution.solutionInfo)[@"Steady state convergence tolerance"] doubleValue];
            if (change <= tolerance) {
                solution.variable.steadyConverged = 1;
            } else {
                solution.variable.steadyConverged = 0;
            }
        }
    } else {
        solution.variable.nonLinChange = change;
        if ((solution.solutionInfo)[@"Non linear system convergence tolerance"] != nil) {
            tolerance = [(solution.solutionInfo)[@"Non linear system convergence tolerance"] doubleValue];
            if (change <= tolerance) {
                solution.variable.nonLinConverged = 1;
            } else {
                solution.variable.nonLinConverged = 0;
            }
        }
    }
    
    if (relax == YES && relaxBefore == NO) {
        for (i=0; i<n; i++) {
            x[i] = (1-relaxation)*x0[i] + relaxation*x[i];
        }
        solution.variable.norm = [self FEMKernel_computeNorm:solution :n :x];
    }
    
    if ((solution.solutionInfo)[@"Equation"] != nil) {
        solverName = [NSString stringWithString:(solution.solutionInfo)[@"Equation"]];
    } else {
        solverName = solution.variable.name;
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
            if ((solution.solutionInfo)[@"Calculate velocity"] != nil) {
                doIt = [(solution.solutionInfo)[@"Calculate velocity"] boolValue];
            } else {
                doIt = NO;
            }
        } else {
            if ((solution.solutionInfo)[@"Non Linear calculate velocity"] != nil) {
                doIt = [(solution.solutionInfo)[@"Non Linear calculate velocity"] boolValue];
            } else {
                doIt = NO;
            }
        }
        
        if (doIt == YES) {
            timeStepVar = [utilities getVariableFrom:solution.mesh.variables model:aModel name:@"Time step size" onlySearch:NULL maskName:NULL info:&stat];
            containers = timeStepVar.getContainers;
            dt = containers->Values[0];
            str2 = @" velocity";
            str1 = [NSMutableString stringWithString:solution.variable.name];
            [str1 appendString:str2];
            veloVar = [utilities getVariableFrom:solution.mesh.variables model:aModel name:str1 onlySearch:NULL maskName:NULL info:&stat];
            containers = veloVar.getContainers;
            for (i=0; i<n; i++) {
                containers->Values[i] = (x[i] - varContainers->PrevValues[i][0]) / dt;
            }
            containers = NULL;
        }
        
    }
    
    // Calculate derivative a.k.a sensitivity
    if (steadyState == YES) {
    
        stat = NO;
        
        if ((solution.solutionInfo)[@"Calculate derivative"] != nil) {
            if ([(solution.solutionInfo)[@"Calculate derivative"]  boolValue] == YES) {
                
                if (iterNo > 1) {
                    timeStepVar = [utilities getVariableFrom:solution.mesh.variables model:aModel name:@"Derivative eps" onlySearch:NULL maskName:NULL info:&stat];
                    if (timeStepVar != nil) {
                        containers = timeStepVar.getContainers;
                        eps = containers->Values[0];
                        stat = YES;
                        containers = NULL;
                    } else {
                        eps = [(solution.solutionInfo)[@"derivative eps"] doubleValue];
                    }
                    if (stat == NO) {
                        warnfunct("FEMKernel_computeChange", "Derivative eps not given, using one.");
                    }
                    
                    str2 = @" derivative";
                    str1 = [NSMutableString stringWithString:solution.variable.name];
                    [str1 appendString:str2];
                    veloVar = [utilities getVariableFrom:solution.mesh.variables model:aModel name:str1 onlySearch:NULL maskName:NULL info:&stat];
                    if (veloVar != nil) {
                        NSLog(@"FEMKernel_computeChange: Computing variable: %@\n", str1);
                        containers = veloVar.getContainers;
                        for (i=0; i<n; i++) {
                            containers->Values[i] = (x[i] - x0[i]) / eps;
                        }
                        containers = NULL;
                    } else {
                        warnfunct("FEMKernel_computeChange", "Derivative variable not present.");
                    }
                }
            }
        }
    }
    
    x = NULL;
    x0 = NULL;
    matContainers = NULL;
    varContainers = NULL;
}

-(void)FEMKernel_solveLinearSystem:(FEMSolution *)solution model:(FEMModel *)aModel {
    
    int i, j, n;
    
    double *diag, *x;
    double diagReal, diagImag, norm, bnorm, sum;
    double complex cmpx;
    NSString *method;
    matrixArraysContainer *matContainers;
    variableArraysContainer *varContainers;
    FEMParallelMPI *parallelUtil;
    BOOL scaleSystem, eigenAnalysis, harmonicAnalysis;
    
    parallelUtil = [[FEMParallelMPI alloc] init];
    matContainers = solution.matrix.getContainers;
    
    n = solution.matrix.numberOfRows;
    
    varContainers = solution.variable.getContainers;
    
    if (solution.matrix.isLumped == YES && [solution timeOrder] == 1)  {
        
        x = varContainers->Values;
        
        if ((solution.solutionInfo)[@"Time stepping Method"] != nil) {
            method = [NSString stringWithString:(solution.solutionInfo)[@"Time stepping Method"]];
            if ([method isEqualToString:@"runge-kutta"] == YES || [method isEqualToString:@"explicit euler"] == YES) {
                for (i=0; i<n; i++) {
                    if (fabs(matContainers->Values[matContainers->Diag[i]]) > 0.0)
                        varContainers->Values[i] = matContainers->RHS[i] / matContainers->Values[matContainers->Diag[i]];
                }
                [self FEMKernel_backRotateNTSystem:solution];
                [self FEMKernel_computeNorm:solution :n :x];
                x = NULL;
                matContainers = NULL;
                varContainers = NULL;
                return;
            }
        }
        
    }
    
    if ((solution.solutionInfo)[@"Linear system scaling"] != nil) {
        scaleSystem = [(solution.solutionInfo)[@"Linear system scaling"] boolValue];
    } else {
        scaleSystem = YES;
    }
    
    eigenAnalysis = (solution.nOfEigenValues  > 0 && [(solution.solutionInfo)[@"Eigen analysis"] boolValue] == YES) ? YES: NO;
    harmonicAnalysis = (solution.nOfEigenValues  > 0 && [(solution.solutionInfo)[@"Harmonic analysis"] boolValue] == YES) ? YES: NO;
    
    if ( harmonicAnalysis == NO && eigenAnalysis == NO ) {
        sum = 0.0;
        for (i=0; i<n; i++) {
            sum = sum + pow(matContainers->RHS[i], 2.0);
        }
        bnorm = [parallelUtil parallelReduction:sqrt(sum) :NULL];
        if (bnorm == 0.0) {
            warnfunct("SolveSystem", "Solution trivially zero");
            for (i=0; i<n; i++) {
                varContainers->Values[i] = 0.0;
            }
            matContainers = NULL;
            varContainers = NULL;
            return;
        }
    }
    
    if ( scaleSystem == YES && harmonicAnalysis == NO ) {
        
        // Scale System Ax = b as:
        // (DAD)y = Db, where D = 1/sqrt(diag(A)) and y = D^-1 x
        diag = doublevec(0, n-1);
        
        if (solution.matrix.isComplexMatrix == YES) {
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
        
        if (solution.matrix.isComplexMatrix == YES) {
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
            for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
                matContainers->Values[j] = matContainers->Values[j] * (diag[i] * diag[matContainers->Cols[j]]);
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
    }
        
    // If solving harmonic anaysis go there:
    if (harmonicAnalysis == YES) {
        
        // TODO:
        
    }
    
    // If solving eigen system go there:
    if (eigenAnalysis) {
        
        // TODO:
        
    }
    
    // Convert rhs & initial value to the scaled system
    if (scaleSystem == YES) {
        for (i=0; i<n; i++) {
            matContainers->RHS[i] = matContainers->RHS[i] * diag[i];
        }
    }
    sum = 0.0;
    for (i=0; i<n; i++) {
        sum = sum + pow(matContainers->RHS[i], 2.0);
    }
    bnorm = [parallelUtil parallelReduction:sqrt(sum) :NULL];
    
    if (scaleSystem == YES) {
        for (i=0; i<n; i++) {
            matContainers->RHS[i] = matContainers->RHS[i] / bnorm;
            varContainers->Values[i] = varContainers->Values[i] / diag[i] / bnorm;
        }
    }
    
    if ((solution.solutionInfo)[@"Linear system solver"] != nil) {
        method = [NSString stringWithString:(solution.solutionInfo)[@"Linear system solver"]];
    } else {
        method = @"iterative";
    }
   
    if ([method isEqualToString:@"iterative"] == YES) {
        
        [self FEMKernel_iterSolver:solution];
    } 
    else if ([method isEqualToString:@"multigrid"] == YES) {
        // TODO: Need to be implemented
    }
    else { // By default we use the direct solver
        // TODO: Need to be implemented
    }
    
    if (scaleSystem == YES) {
        
        // Solve x: INV(D)x = y, scale b back to original
        for (i=0; i<n; i++) {
            varContainers->Values[i] = varContainers->Values[i]* diag[i] * bnorm;
            matContainers->RHS[i] = matContainers->RHS[i] / diag[i] * bnorm;
        }
        
        // Scale the system back to original
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
            if (matContainers->sizeValues != matContainers->sizeDampValues) {
                for (i=0; i<n; i++) {
                    for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
                        matContainers->DampValues[j] = matContainers->DampValues[j] / (diag[i] * diag[matContainers->Cols[j]]);
                    }
                }
            }
        }
    
        free_dvector(diag, 0, n-1);
        
    }
    
    // TODO: Add support for the treatment of variable loads and constrained matrix here
    // ....
    
    [self FEMKernel_backRotateNTSystem:solution];
    
    // Compute the change of the solution with different methods
    [self FEMKernel_computeChange:solution model:aModel isSteadyState:NO nsize:&n values:varContainers->Values values0:NULL];
    norm = solution.variable.norm;
    
    matContainers = NULL;
    varContainers = NULL;
}

-(void)FEMKernel_solveSystem:(FEMSolution *)solution model:(FEMModel *)aModel {
    
    int n, i;
    double relaxation, beta, gamma;
    double t0, rt0, st, rst;
    BOOL needPrevSol, stat;
    NSString *method;
    variableArraysContainer *varContainers;
    
    if ((solution.solutionInfo)[@"Linear system timing"] != nil) {
        if ([(solution.solutionInfo)[@"Linear system timing"] boolValue] == YES) {
            t0 = cputime();
            rt0 = realtime();
        }
    }
    
    n = solution.matrix.numberOfRows;
    
    varContainers = solution.variable.getContainers;

    // The allocation of previous values has to be here in order to work properly
    // with the Dirichlet elimination
    if ((solution.solutionInfo)[@"Non linear system relaxation factor"] != nil) {
        relaxation = [(solution.solutionInfo)[@"Non linear system relaxation factor"] doubleValue];
        needPrevSol = (relaxation != 0.0) ? YES: NO;
    }
    
    if (needPrevSol == NO) {
        if ((solution.solutionInfo)[@"Non linear system convergence measure"] != nil) {
            method = [NSString stringWithString:(solution.solutionInfo)[@"Non linear system convergence measure"]];
            needPrevSol = ([method isEqualToString:@"residual"] == YES || [method isEqualToString:@"solution"] == YES) ? YES: NO;
        }
    }
    
    if (needPrevSol == YES) {
        if (varContainers->NonLinValues != NULL) {
            stat = YES;
            if (varContainers->sizeNonLinValues != n) {
                free_dvector(varContainers->NonLinValues, 0, varContainers->sizeNonLinValues-1);
                varContainers->NonLinValues = NULL;
                stat = NO;
            }
        }
        if (stat == NO) {
            varContainers->NonLinValues = doublevec(0, n-1);
            if (varContainers->NonLinValues == NULL) errorfunct("FEMKernel_solveSystem", "Memory allocation error.");
            varContainers->sizeNonLinValues = n;
        }
        for (i=0; i<n; i++) {
            varContainers->NonLinValues[i] = varContainers->Values[i];
        }
    }
    
    
    // TODO: Add support for constrained matrix
    
    [self FEMKernel_solveLinearSystem:solution model:aModel];
    
    if ([solution timeOrder] == 2) {
        if (varContainers->PrevValues != NULL) {
            
            gamma = 0.5 - [solution alpha];
            beta = ( pow((1.0 - [solution alpha]), 2.0) ) / 4.0;
            for (i = 0; i<n; i++) {
                varContainers->PrevValues[i][1] = (1.0/(beta*pow([solution dt], 2.0))) * (varContainers->Values[i]-varContainers->PrevValues[i][2]) -
                                            (1.0/(beta*[solution dt]))*varContainers->PrevValues[i][3] + (1.0-1.0/(2.0*beta))*varContainers->PrevValues[i][4];
                varContainers->PrevValues[i][0] = varContainers->PrevValues[i][3] + [solution dt]*((1.0-gamma)*varContainers->PrevValues[i][4]+
                                            gamma*varContainers->PrevValues[i][1]);
            }
        }
    }
    
    if ((solution.solutionInfo)[@"Linear system timing"] != nil) {
        if ([(solution.solutionInfo)[@"Linear system timing"] boolValue] == YES) {
            st = cputime() - t0;
            rst = realtime() - rt0;
            
            NSLog(@"Solution time for %@: %lf %lf (s)\n", solution.variable.name, st, rst);
            
            // TODO: Add support for cumulative timing
        
        }
    }
    
    varContainers = NULL;
}

#pragma mark Element info

-(BOOL)FEMKernel_checkPassiveElement:(FEMModel *)model: (FEMSolution *)solution: (Element_t *)element {
    
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
    
    listUtil = [[FEMListUtilities alloc] init];
    
    nbNodes = element->Type.NumberOfNodes;
    
    if ((model.bodies)[body_id] == nil) {
        return isPassive;
    } else {
        bf_id = [(model.bodies)[body_id][@"Body force"] intValue];
    }
    
    passName = [NSMutableString stringWithString:solution.variable.name];
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
    
    free_dvector(passive.vector, 0, passive.m-1);
    passive.vector = NULL;
    
    return isPassive;
}

#pragma mark Manipulate matrix

-(void)FEMKernel_rotateMatrix:(FEMSolution *)solution: (double **)matrix: (double *)vector: (int)n: (int)dim: (int)dofs: (int *)nodeIndexes {
    
    int i, j, k, l;
    double s, **r, **q, *n1, *t1, *t2;
    solutionArraysContainer *solContainers;
    
    r = doublematrix(0, (n*dofs)-1, 0, (n*dofs)-1);
    q = doublematrix(0, (n*dofs)-1, 0, (n*dofs)-1);
    n1 = doublevec(0, 2);
    t1 = doublevec(0, 2);
    t2 = doublevec(0, 2);
    
    solContainers = solution.getContainers;
    
    for (i=0; i<n; i++) {
        
        if (nodeIndexes[i] < 0 || nodeIndexes[i]+1 > solContainers->size1boundaryNormals) continue;
        
        for (j=0; j<n*dofs; j++) {
            for (k=0; k<n*dofs; k++) {
                r[j][k] = 0.0;
            }
        }
        
        for (j=0; j<n*dofs; j++) {
            r[j][j] = 1.0;
        }
        
        for (j=0; j<2; j++) {
            n1[j] = solContainers->boundaryNormals[nodeIndexes[i]][j];
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
                    t1[j] = solContainers->boundaryTangent1[nodeIndexes[i]][j];
                    t2[j] = solContainers->boundaryTangent2[nodeIndexes[i]][j];
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
    solContainers = NULL;
}

-(void)FEMKernel_updateGlobalForce:(FEMModel *)model: (FEMSolution *)solution: (Element_t *)element: (double *)forceVector: (double *)localForce: (int)n: (int)dofs: (int *)nodeIndexes: (BOOL *)rotateNT {
    
    int i, j, k, dim, *indexes;
    BOOL rotate;
    double **localStiffMatrix;
    solutionArraysContainer *solContainers;
    
    // Check if this element has been defined as passive
    if ([self FEMKernel_checkPassiveElement:model :solution :element] == YES) return;
    
    localStiffMatrix = doublematrix(0, (n*dofs)-1, 0, (n*dofs)-1);
    indexes = intvec(0, n-1);
    
    solContainers = solution.getContainers;
    
    rotate = YES;
    if (rotateNT != NULL) rotate = *rotateNT;
    
    if (rotate == YES && [solution normalTangentialNOFNodes] > 0) {
        
        dim = [model dimension];
         memset( indexes, 0.0, (n*sizeof(indexes)) );
        for (i=0; i<element->Type.NumberOfNodes; i++) {
            indexes[i] = solContainers->boundaryReorder[element->NodeIndexes[i]];
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
    solContainers = NULL;
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
    BOOL constantDt, found;
    FEMTimeIntegration *timeIntegration;
    FEMUtilities *utilities;
    NSString *method;
    FEMVariable *dtVar;
    matrixArraysContainer *matContainers;
    variableArraysContainer *varContainers, *containers;
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
    
    prevSol = doublematrix(0, (dofs*n)-1, 0, [solution order]-1);
    lForce = doublevec(0, (dofs*n)-1);
    dts = doublevec(0, [solution order]-1);
    matrixForce = doublevec(0, matContainers->size1force-1);
    
    buffer = doublevec(0, (dofs*n)-1);
    
    timeIntegration = [[FEMTimeIntegration alloc] init];
    utilities = [[FEMUtilities alloc] init];
    
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
    
    order = min([solution doneTime], [solution order]);
    
    for (i=0; i<n; i++) {
        for (j=0; j<dofs; j++) {
            k = dofs * i + j;
            l = dofs * nodeIndexes[i] + j;
            for (m=0; n<order; m++) {
                prevSol[k][m] = varContainers->PrevValues[l][m];
            }
        }
    }
    
    for (i=0; i<n*dofs; i++) {
        lForce[i] = force[i];
    }
    
    for (i=0; i<matContainers->size1force; i++) {
        matrixForce[i] = matContainers->Force[i][0];
    }
    
    [self FEMKernel_updateGlobalForce:model :solution :element :matrixForce :lForce :n :dofs :nodeIndexes :NULL];
    
    if ((solution.solutionInfo)[@"Time stepping method"] != nil) {
        method = [NSString stringWithString:(solution.solutionInfo)[@"Time stepping method"]];
    } else {
        method = @"bdf";
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
            dtVar = [utilities getVariableFrom:solution.mesh.variables model:model name:@"Time step size" onlySearch:NULL maskName:NULL info:&found];
            containers = dtVar.getContainers;
            for (i=1; i<order; i++) {
                dts[i] = containers->PrevValues[0][i-1];
                if ( fabs(dts[i]-dts[0]) > 1.0e-6 * dts[0] ) constantDt = NO;
            }
            containers = NULL;
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
    free_dvector(matrixForce, 0, matContainers->size1force-1);
    free_dvector(buffer, 0, (dofs*n)-1);
    
    matContainers = NULL;
    varContainers = NULL;
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
    matrixArraysContainer *matContainers;
    solutionArraysContainer *solContainers;
    
    // Check if this element has been defined as passive
    if ([self FEMKernel_checkPassiveElement:model :solution :element] == YES) return;
    
    matContainers = solution.matrix.getContainers;
    solContainers = solution.getContainers;
    
    indexes = intvec(0, n-1);
    
    rotate = YES;
    if (rotateNT != NULL) {
        rotate = *rotateNT;
    }
    
    if (rotate == YES && [solution normalTangentialNOFNodes] > 0) {
        dim = [model dimension];
        memset( indexes, 0.0, (n*sizeof(indexes)) );
        for (i=0; i<element->Type.NumberOfNodes; i++) {
            indexes[i] = solContainers->boundaryReorder[element->NodeIndexes[i]];
        }
        [self FEMKernel_rotateMatrix:solution :stiffMatrix :force :n :dim :dofs :indexes];
    }
    
    if (solution.matrix.format == MATRIX_CRS) {
        
        crsMatrix = [[FEMMatrixCRS alloc] init];
        [crsMatrix glueLocalMatrixInGlobal:solution :stiffMatrix :n :dofs :nodeIndexes];
        
    }
    else if (solution.matrix.format == MATRIX_BAND || solution.matrix.format == MATRIX_SBAND) {
        
        bandMatrix = [[FEMMatrixBand alloc] init];
        [bandMatrix glueLocalMatrixInGlobal:solution :stiffMatrix :n :dofs :nodeIndexes];
        
    }
    
    
    if (bulkUpdate == NULL || *bulkUpdate == NO) {
        for (i=0; i<n; i++) {
            if (nodeIndexes[i] >= 0) {
                for (j=0; j<dofs; j++) {
                    k = dofs * nodeIndexes[i] + j;
                    matContainers->RHS[k] = matContainers->RHS[k] + force[dofs*i+j];
                }
            }
        }
    } else if (*bulkUpdate == YES ) {
        for (i=0; i<n; i++) {
            if (nodeIndexes[i] >= 0) {
                for (j=0; j<dofs; j++) {
                    k = dofs * nodeIndexes[i] + j;
                    matContainers->BulkRHS[k] = matContainers->BulkRHS[k] + force[dofs*i+j];
                }
            }
        }
    }
    
    matContainers = NULL;
    solContainers = NULL;
}

#pragma mark Loads

-(void)FEMKernel_setElementLoads:(FEMModel *)model: (FEMSolution *)solution: (Element_t *)element: (NSArray *)values: (NSMutableString *)name: (int*)indexes: (BOOL *)doneLoad: (int)n: (int)dof: (int)ndofs {
    
    int j, k, l, k1;
    listBuffer work = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    listBuffer workA = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    FEMListUtilities *listUtil;
    NSString *str;
    matrixArraysContainer *matContainers;
    variableArraysContainer *varContainers;
    BOOL stat;
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
        
    listUtil = [[FEMListUtilities alloc] init];
    
    if (dof >= 0) {
        
        stat = [listUtil listGetReal:model inArray:values forVariable:name numberOfNodes:n indexes:indexes buffer:&work minValue:NULL maxValue:NULL];
        if (stat == NO) {
            [name appendString:@" dofs"];
            str = [NSString stringWithString:name];
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
                    matContainers->RHS[k] = matContainers->RHS[k] + work.vector[j];
                } else {
                    for (l=0; l<min(ndofs, workA.m); l++) {
                        k1 = ndofs * k + l;
                        matContainers->RHS[k1] = matContainers->RHS[k1] + workA.tensor[l][0][j];
                    }
                }
            }
        }
    
    }
    
    if (workA.tensor != NULL) {
        free_d3tensor(workA.tensor, 0,  workA.m-1, 0,  workA.n-1, 0, workA.p-1);
        workA.tensor = NULL;
    }
    if (work.vector != NULL) {
        free_dvector(work.vector, 0, work.m-1);
        work.vector = NULL;
    }
    
    matContainers = NULL;
    varContainers = NULL;
}

-(void)FEMKernel_setPointLoads:(FEMModel *)model: (FEMSolution *)solution: (Element_t *)element: (NSArray *)values: (NSMutableString *)name: (int*)indexes: (int)n: (int)dof: (int)ndofs {
    
    int j, k, l, k1;
    listBuffer work = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    listBuffer workA = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    FEMListUtilities *listUtil;
    matrixArraysContainer *matContainers;
    variableArraysContainer *varContainers;
    BOOL stat;
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
    
    listUtil = [[FEMListUtilities alloc] init];
    
    if (dof >= 0) {
        
        stat = [listUtil listGetReal:model inArray:values forVariable:name numberOfNodes:n indexes:indexes buffer:&work minValue:NULL maxValue:NULL];
        
    } else {
        
        stat = [listUtil listGetRealArray:model inArray:values forVariable:name numberOfNodes:n indexes:indexes buffer:&workA];
    }
    
    if (stat == YES) {
        
        for (j=0; j<n; j++) {
            
            if (indexes[j] > varContainers->sizePerm || indexes[j] < 0) {
                warnfunct("FEMKernel_setPointLoads", "Invalid node number");
                continue;
            }
            
            k = varContainers->Perm[indexes[j]];
            if (k >= 0) {
                if (dof >= 0) {
                    k = ndofs * k + dof;
                    matContainers->RHS[k] = matContainers->RHS[k] + work.vector[j];
                } else {
                    for (l=0; l<min(ndofs, workA.m); l++) {
                        k1 = ndofs * k + l;
                        matContainers->RHS[k1] = matContainers->RHS[k1] + workA.tensor[l][0][j];
                    }
                }
            }
        }
        
    }
    
    if (workA.tensor != NULL) {
        free_d3tensor(workA.tensor, 0,  workA.m-1, 0,  workA.n-1, 0, workA.p-1);
        workA.tensor = NULL;
    }
    if (work.vector != NULL) {
        free_dvector(work.vector, 0, work.m-1);
        work.vector = NULL;
    }
    matContainers = NULL;
    varContainers = NULL;
}

#pragma mark Element and point values

-(void)FEMKernel_setElementValues:(FEMModel *)model inSolution:(FEMSolution *)solution forElementNumber:(int)elno numberOfNodes:(int)n atIndexes:(int *)indexes forValues:(NSArray *)values variableName:(NSMutableString *)name orderOfDofs:(int)dof activeCondition:(BOOL)conditional conditionName:(NSString *)condName permutationOffset:(int)offset {
    
    int i, j, k, l, m, dim, k1;
    BOOL checkNT, stat, all;
    double *rotvec;
    listBuffer condition = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    listBuffer work = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    listBuffer workA = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    NSMutableString *str;
    FEMListUtilities *listUtil;
    FEMMatrixCRS *crsMatrix;
    FEMMatrixBand *bandMatrix;
    matrixArraysContainer *matContainers;
    variableArraysContainer *varContainers;
    solutionArraysContainer *solContainers;
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
    solContainers = solution.getContainers;
    
    listUtil = [[FEMListUtilities alloc] init];
    
    rotvec = doublevec(0, 2);
    
    dim = [model dimension];
    
    str = [NSMutableString stringWithString:name];
    [str appendString:@" dofs"];
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
        if ([solution normalTangentialNOFNodes] > 0 && dof >= 0) {
            checkNT = YES;
            all = YES;
            for (i=0; i<n; i++) {
                k = indexes[i];
                if (solContainers->boundaryReorder[k] < 0) {
                    continue;
                } else {
                    all = NO;
                    break;
                }
            }
            if (all == YES) checkNT = NO;
            if ([listUtil listGetLogical:model inArray:values forVariable:[solution normalTangentialName] info:&stat] == YES) checkNT = NO;
        }
        
        for (j=0; j<n; j++) {
            
            @autoreleasepool {
                if (conditional == YES && condition.vector[j] < 0.0) continue;
                k = varContainers->Perm[indexes[j]];
                
                if (k >= 0) {
                    if (dof >= 0) {
                        m = 0;
                        if ([solution normalTangentialNOFNodes] > 0) m = solContainers->boundaryReorder[indexes[j]];
                        if (m >= 0 && checkNT == YES) {
                            memset( rotvec, 0.0, (3*sizeof(rotvec)) );
                            rotvec[dof] = 1.0;
                            [self FEMKernel_rotateNTSystem:solution :rotvec :indexes[j]];
                            for (k=0; k<dim; k++) {
                                if (fabs(rotvec[k]) > 1.0e-8) {
                                    if (solContainers->ntElement[m][k] == elno) {
                                        l = solution.variable.dofs * varContainers->Perm[indexes[j]] + k;
                                        if (solContainers->ntZeroingDone[m][k] == false) {
                                            matContainers->RHS[l] = 0.0;
                                            [self zeroTheNumberOfRows:l inSolutionMatrix:solution];
                                            [self setMatrixElement:solution :l :l :1.0];
                                            solContainers->ntZeroingDone[m][k] = true;
                                        }
                                        matContainers->RHS[l] = matContainers->RHS[l] + rotvec[k]+work.vector[j];
                                    }
                                }
                            }
                        } else {
                            k = offset + solution.variable.dofs * k + dof;
                            if (solution.matrix.format == MATRIX_SBAND) {
                                bandMatrix = [[FEMMatrixBand alloc] init];
                                [bandMatrix sBand_setDirichlet:solution :k :work.vector[j]];
                            } else if (solution.matrix.format == MATRIX_CRS && solution.matrix.isSymmetric == YES) {
                                crsMatrix = [[FEMMatrixCRS alloc] init];
                                [crsMatrix setSymmetricDirichletInGlobal:solution :k :work.vector[j]];
                            } else {
                                matContainers->RHS[k] = work.vector[j];
                                [self zeroTheNumberOfRows:k inSolutionMatrix:solution];
                                [self setMatrixElement:solution :k :k :1.0];
                            }
                        }
                    } else {
                        for (l=0; l<min(solution.variable.dofs, workA.m); l++) {
                            k1 = offset + solution.variable.dofs*k + l;
                            if (solution.matrix.format == MATRIX_SBAND) {
                                bandMatrix = [[FEMMatrixBand alloc] init];
                                [bandMatrix sBand_setDirichlet:solution :k1 :workA.tensor[l][0][j]];
                            } else if (solution.matrix.format == MATRIX_CRS && solution.matrix.isSymmetric == YES) {
                                crsMatrix = [[FEMMatrixCRS alloc] init];
                                [crsMatrix setSymmetricDirichletInGlobal:solution :k1 :workA.tensor[l][0][j]];
                            } else {
                                matContainers->RHS[k1] = workA.tensor[l][0][j];
                                [self zeroTheNumberOfRows:k1 inSolutionMatrix:solution];
                                [self setMatrixElement:solution :k1 :k1 :1.0];
                            }
                        }
                    }
                }
            }
            
        }
    }
    
    if (condition.vector != NULL) {
        free_dvector(condition.vector, 0, condition.m-1);
        condition.vector = NULL;
    }
    if (work.vector != NULL) {
        free_dvector(work.vector, 0, work.m-1);
        work.vector = NULL;
    }
    if (workA.tensor != NULL) {
        free_d3tensor(workA.tensor, 0,  workA.m-1, 0,  workA.n-1, 0, workA.p-1);
        workA.tensor = NULL;
    }
    free_dvector(rotvec, 0, 2);
    matContainers = NULL;
    varContainers = NULL;
    solContainers = NULL;
}

-(void)FEMKernel_setPointValues:(FEMModel *)model inSolution:(FEMSolution *)solution numberofNodes:(int)n atIndexes:(int *)indexes forValues:(NSArray *)values variableName:(NSMutableString *)name orderOfDofs:(int)dof activeCondition:(BOOL)conditional conditionName:(NSString *)condName permutationOffset:(int)offset {
    
    int j, k, l, k1;
    BOOL stat;
    listBuffer condition = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    listBuffer work = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    listBuffer workA = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    FEMListUtilities *listUtil;
    FEMMatrixCRS *crsMatrix;
    FEMMatrixBand *bandMatrix;
    matrixArraysContainer *matContainers;
    variableArraysContainer *varContainers;
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
    
    listUtil = [[FEMListUtilities alloc] init];
    
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
                if (indexes[j] > varContainers->sizePerm || indexes[j] < 0) {
                    warnfunct("FEMKernel_setPointValues", "Invalid node number");
                    continue;
                }
                
                k = varContainers->Perm[indexes[j]];
                
                if (k >= 0) {
                    if (dof >= 0) {
                        k = offset + solution.variable.dofs*k + dof;
                        if (solution.matrix.format == MATRIX_SBAND) {
                            bandMatrix = [[FEMMatrixBand alloc] init];
                            [bandMatrix sBand_setDirichlet:solution :k :work.vector[j]];
                        } else if (solution.matrix.format == MATRIX_CRS && solution.matrix.isSymmetric == YES) {
                            crsMatrix = [[FEMMatrixCRS alloc] init];
                            [crsMatrix setSymmetricDirichletInGlobal:solution :k :work.vector[j]];
                        } else {
                            matContainers->RHS[k] = work.vector[j];
                            [self zeroTheNumberOfRows:k inSolutionMatrix:solution];
                            [self setMatrixElement:solution :k :k :1.0];
                        }
                    } else {
                        for (l=0; l<min(solution.variable.dofs, workA.m); l++) {
                            k1 = offset + solution.variable.dofs*k + l;
                            if (solution.matrix.format== MATRIX_SBAND) {
                                bandMatrix = [[FEMMatrixBand alloc] init];
                                [bandMatrix sBand_setDirichlet:solution :k1 :workA.tensor[l][0][j]];
                            } else if (solution.matrix.format == MATRIX_CRS && solution.matrix.isSymmetric == YES) {
                                crsMatrix = [[FEMMatrixCRS alloc] init];
                                [crsMatrix setSymmetricDirichletInGlobal:solution :k1 :workA.tensor[l][0][j]];
                            } else {
                                matContainers->RHS[k1] = workA.tensor[l][0][j];
                                [self zeroTheNumberOfRows:k1 inSolutionMatrix:solution];
                                [self setMatrixElement:solution :k1 :k1 :1.0];
                            }                        
                        }
                    }
                }
            }
            
        }
    }

    
    if (condition.vector != NULL) {
        free_dvector(condition.vector, 0, condition.m-1);
        condition.vector = NULL;
    }
    if (work.vector != NULL) {
        free_dvector(work.vector, 0, work.m-1);
        work.vector = NULL;
    }
    
    if (workA.tensor != NULL) {
        free_d3tensor(workA.tensor, 0,  workA.m, 0,  workA.n-1, 0, workA.p-1);
        workA.tensor = NULL;
    }
    matContainers = NULL;
    varContainers = NULL;
}

#pragma mark Periodics

-(void)FEMKernel_setPeriodicBoundariesPass1:(FEMModel *)model: (FEMSolution *)solution: (NSMutableString *)name: (int)dof: (int)this: (BOOL *)done {
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
    FEMMatrix *a, *b;
    FEMVariable *var;
    FEMListUtilities *listUtil;
    FEMUtilities *util;
    FEMMatrixCRS *crsMatrix;
    FEMBoundaryCondition *boundaryConditionAtId;
    NSMutableString *str1, *str2, *str3;
    matrixArraysContainer *matContainers, *projectorContainers, *aContainers;
    variableArraysContainer *varContainers, *containers;
    
    listUtil = [[FEMListUtilities alloc] init];
    util = [[FEMUtilities alloc] init];
    crsMatrix = [[FEMMatrixCRS alloc] init];
    
    str1 = [NSMutableString stringWithString:@"Periodic BC "];
    [str1 appendString:name];
    
    str2 = [NSMutableString stringWithString:@"Anti periodic BC "];
    [str2 appendString:name];
    
    str3 = [NSMutableString stringWithString:@"Periodic BC scale "];
    [str3 appendString:name];

    boundaryConditionAtId = (model.boundaries)[this];
    
    if ([listUtil listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:str1 info:&stat] == YES) {
        scale = -1.0;
    } else if ([listUtil listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:str2 info:&stat] == YES) {
        scale = 1.0;
    } else {
        scale = [listUtil listGetConstReal:model inArray:boundaryConditionAtId.valuesList forVariable:str3 info:&stat minValue:NULL maxValue:NULL];
        if (stat == NO) return;
    }
    
    if (boundaryConditionAtId.pMatrix != NULL) return;
    
    matContainers = solution.matrix.getContainers;
    projectorContainers = boundaryConditionAtId.pMatrix.getContainers;
    varContainers = solution.variable.getContainers;
    
    // For explicit conditions just create the dependency almost like a normal Dirichlet BC, 
    // for implicit one (otherwise) do the assembly of the projector
    str1 = [NSMutableString stringWithString:@"Periodic BC explicit"];
    str2 = [NSMutableString stringWithString:@"Periodic BC use Lagrange coefficient"];
    if ([listUtil listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:str1 info:&stat] == YES) {
        
        var = [util getVariableFrom:model.variables model:model name:name onlySearch:NULL maskName:NULL info:&stat];
        containers = var.getContainers;
        
        for (i=0; i<boundaryConditionAtId.pMatrix.numberOfRows; i++) {
            ii = projectorContainers->InvPerm[i];
            k = varContainers->Perm[ii];
            if (done[ii] == NO && k >= 0) {
                k = solution.variable.dofs * k + dof;
                [self zeroTheNumberOfRows:k inSolutionMatrix:solution];
                [self addToMatrixElement:solution :k :k :1.0];
                matContainers->RHS[k] = 0.0;
                
                for (l=projectorContainers->Rows[i]; l<=projectorContainers->Rows[i+1]-1; l++) {
                    if (projectorContainers->Cols[l] < 0) continue;
                    m = varContainers->Perm[projectorContainers->Cols[l]];
                    if (m >= 0) {
                        m = solution.variable.dofs * m + dof;
                        matContainers->RHS[k] = matContainers->RHS[k] - scale*projectorContainers->Values[l]*containers->Values[m];
                    }
                }
            }
        }
        containers = NULL;
        
    } else if ([listUtil listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:str2 info:&stat] == YES) {
        
        b = nil;
        a = solution.matrix.constraint;
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
            a = [util allocateMatrix];
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
            [crsMatrix sortInMatrix:a :NULL];
            
            memset( aContainers->RHS, 0.0, (n*sizeof(aContainers->RHS)) );
            memset( aContainers->Values, 0.0, ((projectorContainers->sizeValues*pow(solution.variable.dofs, 2.0)+n)*sizeof(aContainers->Values)) );
        }
        
        for (i=0; i<boundaryConditionAtId.pMatrix.numberOfRows; i++) {
            k = varContainers->Perm[projectorContainers->InvPerm[i]];
            if (k < 0) continue;
            
            [util zeroTheNumberOfRows:solution.variable.dofs*i+dof inMatrix:a];
            [util setMatrixElement:a :solution.variable.dofs*i+dof :solution.variable.dofs*k+dof :scale];
            
            for (j=projectorContainers->Rows[i]; j<=projectorContainers->Rows[i+1]-1; j++) {
                m = varContainers->Perm[projectorContainers->Cols[j]];
                if (m >= 0) [util setMatrixElement:a :solution.variable.dofs*i+dof :solution.variable.dofs*m+dof :projectorContainers->Values[j]];                
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
                            [self addToMatrixElement:solution :m :matContainers->Cols[nn] :projectorContainers->Values[l]*matContainers->Values[nn]];
                        }
                        matContainers->RHS[m] = matContainers->RHS[m] + projectorContainers->Values[l]*matContainers->RHS[k];
                    }
                }
            }
            done[ii] = YES;
        }
    }
    
    matContainers = NULL;
    projectorContainers = NULL;
    aContainers = NULL;
    varContainers = NULL;
}

-(void)FEMKernel_setPeriodicBoundariesPass2:(FEMModel *)model: (FEMSolution *)solution: (NSMutableString *)name: (int)dof: (int)this: (BOOL *)done {
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
    double scale;
    BOOL stat;
    FEMListUtilities *listUtil;
    FEMBoundaryCondition *boundaryConditionAtId;
    NSMutableString *str1, *str2, *str3;
    matrixArraysContainer *matContainers, *projectorContainers;
    variableArraysContainer *varContainers;
    listUtil = [[FEMListUtilities alloc] init];

    str1 = [NSMutableString stringWithString:@"Periodic BC "];
    [str1 appendString:name];
    
    str2 = [NSMutableString stringWithString:@"Anti periodic BC "];
    [str2 appendString:name];
    
    str3 = [NSMutableString stringWithString:@"Periodic BC scale "];
    [str3 appendString:name];
    
    boundaryConditionAtId = (model.boundaries)[this];
    
    if ([listUtil listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:str1 info:&stat] == YES) {
        scale = -1.0;
    } else if ([listUtil listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:str2 info:&stat] == YES) {
        scale = 1.0;
    } else {
        scale = [listUtil listGetConstReal:model inArray:boundaryConditionAtId.valuesList forVariable:str3 info:&stat minValue:NULL maxValue:NULL];
        if (stat == NO) return;
    }
    
    if (boundaryConditionAtId.pMatrix != NULL) return;
    
    matContainers = solution.matrix.getContainers;
    projectorContainers = boundaryConditionAtId.pMatrix.getContainers;
    varContainers = solution.variable.getContainers;
    
    str1 = [NSMutableString stringWithString:@"Periodic BC explicit"];
    str2 = [NSMutableString stringWithString:@"Periodic BC use Lagrange coefficient"];
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
                    [self addToMatrixElement:solution :k :m :projectorContainers->Values[l]];
                }
            }
            matContainers->RHS[k] = 0.0;
            [self addToMatrixElement:solution :k :k :scale];
        }
        done[ii] = YES;
    }
    
    matContainers = NULL;
    projectorContainers = NULL;
    varContainers = NULL;
}

#pragma mark Public methods...

- (id)init
{
    int i;
    
    self = [super init];
    if (self) {
        // Initialization code here.
        _indexStore = intvec(0, 511);
        _sizeIndexStore = 512;
        memset( _indexStore, 0.0, (_sizeIndexStore*sizeof(_indexStore)) );
        
        _lineEM = intmatrix(0, 0, 0, 1);
        _triangleEM = intmatrix(0, 2, 0, 1);
        _quadEM = intmatrix(0, 3, 0, 1);
        _tetraEM = intmatrix(0, 5, 0, 1);
        _prismEM = intmatrix(0, 7, 0, 1);
        _wedgeEM = intmatrix(0, 8, 0, 1);
        _brickEM = intmatrix(0, 11, 0, 1);
        
        for (i=0; i<8; i++) {
            _initialized[i] = NO;
        }
        
    }
    
    return self;
}

- (void)dealloc
{

}

-(void)deallocation
{
    free_ivector(_indexStore, 0, 511);
    free_imatrix(_lineEM, 0, 0, 0, 1);
    free_imatrix(_triangleEM, 0, 2, 0, 1);
    free_imatrix(_quadEM, 0, 3, 0, 1);
    free_imatrix(_tetraEM, 0, 5, 0, 1);
    free_imatrix(_prismEM, 0, 7, 0, 1);
    free_imatrix(_wedgeEM, 0, 8, 0, 1);
    free_imatrix(_brickEM, 0, 11, 0, 1);
    
    free_ivector(_g_Ind, 0, _size_g_Ind-1);
    _g_Ind = NULL;
    
    free_ivector(_l_Ind, 0, _size_l_Ind-1);
    _l_Ind = NULL;
    
    free_dmatrix(_kernStiff, 0, _size1kernStiff-1, 0, _size2kernStiff-1);
    _kernStiff = NULL;
    
     
    free_dvector(_kernWork, 0, _sizekernWork-1);
    _kernWork = NULL;
    
}

-(BOOL)getReal:(FEMModel *)model forElement:(Element_t *)element inList:(NSArray *)list variableName:(NSString *)name buffer:(listBuffer *)result {
    
    int n;
    BOOL found;
    FEMListUtilities *listUtil;
    
    if (element != NULL) {
        n = [self getNumberOfNodesForElement:element];
    } else {
        return NO;
    }
    
    listUtil = [[FEMListUtilities alloc] init];
    
    found = [listUtil listGetReal:model inArray:list forVariable:name numberOfNodes:n indexes:element->NodeIndexes buffer:result minValue:NULL maxValue:NULL];
    
    return found;
    
}

-(int)isPElement:(Element_t *)element {
    
    if (element->Pdefs != NULL) {
        return 1;
    } else {
        return 0;
    }
}

-(void)getNodes:(FEMSolution *)solution inElement:(Element_t *)element resultNodes:(Nodes_t *)nodes numberOfNodes:(int)nd {
    
    int i, n, nb, sz, sz1;
    Nodes_t *globalNodes;
    
    globalNodes = solution.mesh.getNodes;
    
    n = max(solution.mesh.maxElementNodes, solution.mesh.maxElementDofs);
    
    if (nd < n) {
        errorfunct("getNodes", "Size of nodes structure smaller than the required max(maxElementNodes, maxElementDofs).");
    }
    
    n = element->Type.NumberOfNodes;
    
    for (i=0; i<n; i++) {
        nodes[i].x = globalNodes[element->NodeIndexes[i]].x;
        nodes[i].y = globalNodes[element->NodeIndexes[i]].y;
        nodes[i].z = globalNodes[element->NodeIndexes[i]].z;
    }
    
    sz = max(solution.mesh.maxElementNodes, solution.mesh.maxElementDofs); 
    if ( sz > n) {
        for (i=n; i<sz; i++) {
            nodes[i].x = 0.0;
            nodes[i].y = 0.0;
            nodes[i].z = 0.0;
        }
    }
    
    sz1 = solution.mesh.sizeOfGlobalNodes;
    if (sz1 > solution.mesh.numberOfNodes) {
        memset( _indexStore, 0.0, (_sizeIndexStore*sizeof(_indexStore)) );
        nb = [self getElementDofs:solution forElement:element atIndexes:_indexStore];
        for (i=n; i<nb; i++) {
            if (_indexStore[i] >= 0 && _indexStore[i] < sz1) {
                nodes[i].x = globalNodes[_indexStore[i]].x;
                nodes[i].y = globalNodes[_indexStore[i]].y;
                nodes[i].z = globalNodes[_indexStore[i]].z;
            }
        }
    }
    
    globalNodes = NULL;
}

-(int)getElementFamily:(Element_t *)element {
    
    return element->Type.ElementCode / 100;
}

-(int)getElementDofs:(FEMSolution *)solution forElement:(Element_t *)element atIndexes:(int *)indexes {
    
    int nb, i, j, k, bid, edofs, fdofs, faceDofs, edgeDofs, bubbleDofs, ind;
    BOOL gb;
    Element_t *parent, *edges, *faces;
    solutionArraysContainer *solContainers;
    
    edges = solution.mesh.getEdges;
    faces = solution.mesh.getFaces;
    
    nb = 0;
    
    if ((solution.solutionInfo)[@"Discontinuous galerkin"] != nil) {
        if ([(solution.solutionInfo)[@"Discontinuous galerkin"] boolValue] == YES) {
            
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
            
            if (nb > 0) { 
                edges = NULL;
                faces = NULL;
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
        for (i=0; i>element->NDOFs; i++) {
            indexes[nb] = element->NodeIndexes[i];
            nb++;
        }
    }
    for (i=1; i<solContainers->size2DefDofs; i++) {
        if (solContainers->defDofs[bid-1][i] < 0) {
            solContainers = NULL;
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
                indexes[nb] = edgeDofs*(element->EdgeIndexes[j]-1) + i + solution.mesh.numberOfNodes;
                nb++;
            }
        }
    }
    
    if (element->FaceIndexes != NULL) {
        for (j=0; j<element->Type.NumberOfFaces; j++) {
            fdofs = faces[element->FaceIndexes[j]].BDOFs;
            for (i=0; i<fdofs; i++) {
                indexes[nb] = faceDofs*(element->FaceIndexes[j]-1) + i + solution.mesh.numberOfNodes + 
                edgeDofs*solution.mesh.numberOfEdges;
                nb++;
            }
        }
    }
    
    if ((solution.solutionInfo)[@"Bubbles in global system"] != nil) {
        gb = [(solution.solutionInfo)[@"Bubbles in global system"] boolValue];
    } else {
        gb = YES;
    }
    
    if ( element->BoundaryInfo != NULL && (element->EdgeIndexes == NULL || element->FaceIndexes == NULL) ) {
        
        parent = element->BoundaryInfo->Left;
        if (parent == NULL) parent = element->BoundaryInfo->Right;
        if (parent == NULL) {
            edges = NULL;
            faces = NULL;
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
                            for (j=0; element->Type.NumberOfNodes; j++) {
                                if (edges[parent->EdgeIndexes[ind]].NodeIndexes[i] == element->NodeIndexes[j]) k++;
                            }
                        }
                        if (k == element->Type.NumberOfNodes) break;
                    }
                }
                
                edofs = element->BDOFs;
                for (i=0; i<edofs; i++) {
                    indexes[nb] = edgeDofs*(parent->EdgeIndexes[ind]-1) + i + solution.mesh.numberOfNodes; 
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
                    indexes[nb] = faceDofs*(parent->FaceIndexes[ind]-1) + i + solution.mesh.numberOfNodes
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
    
    edges = NULL;
    faces = NULL;
    solContainers = NULL;
    
    return nb;
}

-(int)sgetElementDofs:(FEMSolution *)solution forElement:(Element_t *)element atIndexes:(int *)indexes {
    
    int nb, i, j, edofs, fdofs, faceDofs, edgeDofs, bubbleDofs;
    BOOL gb;
    Element_t *parent, *edges, *faces;
    
    edges = solution.mesh.getEdges;
    faces = solution.mesh.getFaces;
    
    nb = 0;
    
    if ((solution.solutionInfo)[@"Discontinuous galerkin"] != nil) {
        if ([(solution.solutionInfo)[@"Discontinuous galerkin"] boolValue] == YES) {
            
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
            
            if (nb > 0) {
                edges = NULL;
                faces = NULL;
                return nb;
            }
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
                indexes[nb] = edgeDofs*(element->EdgeIndexes[j]-1) + i + solution.mesh.numberOfNodes;
                nb++;
            }
        }
    }
    
    if (element->FaceIndexes != NULL) {
        for (j=0; j<element->Type.NumberOfFaces; j++) {
            fdofs = faces[element->FaceIndexes[j]].BDOFs;
            for (i=0; i<fdofs; i++) {
                indexes[nb] = faceDofs*(element->FaceIndexes[j]-1) + i + solution.mesh.numberOfNodes + 
                edgeDofs*solution.mesh.numberOfEdges;
                nb++;
            }
        }
    }
    
    if ((solution.solutionInfo)[@"Bubbles in global system"] != nil) {
        gb = [(solution.solutionInfo)[@"Bubbles in global system"] boolValue];
    } else {
        gb = YES;
    }
    
    if (element->BoundaryInfo != NULL) {
        
        if ([self isPElement:element] == NO) return nb;
        
        parent = element->BoundaryInfo->Left;
        if (parent == NULL) parent = element->BoundaryInfo->Right;
        if (parent == NULL) {
            edges = NULL;
            faces = NULL;
            return nb;
        }
        
        if (parent->EdgeIndexes != NULL) {
            edofs = element->BDOFs;
            for (i=0; i<edofs; i++) {
                indexes[nb] = edgeDofs*(parent->EdgeIndexes[element->Pdefs->LocalNumber-1]-1) + i + solution.mesh.numberOfNodes;
                nb++;
            }
        }
        
        if (parent->FaceIndexes != NULL) {
            fdofs = element->BDOFs;
            for (i=0; i<fdofs; i++) {
                indexes[nb] = faceDofs*(parent->FaceIndexes[element->Pdefs->LocalNumber-1]-1) + i + solution.mesh.numberOfNodes + edgeDofs*solution.mesh.numberOfEdges;
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
    
    edges = NULL;
    faces = NULL;
    
    return nb;    
}

-(int)getNumberOfNodesForElement:(Element_t *)element {
    
    return element->Type.NumberOfNodes;
}

-(int)getBoundaryConditionID:(FEMModel *)model forElement:(Element_t *)element {
    
    int i, bc_id;
    FEMBoundaryCondition *boundaryConditionAtId;
    
    for (i=0; i<model.numberOfBoundaries; i++) {
        bc_id++;
        boundaryConditionAtId = (model.boundaries)[i];
        if (element->BoundaryInfo->Constraint == [boundaryConditionAtId tag]) break;
    }
    if (bc_id > model.numberOfBoundaries) bc_id = 0;
    
    return bc_id;
}

-(NSArray *)getBoundaryCondition:(FEMModel *)model forElement:(Element_t *)element {
    
    int bc_id;
    
    FEMBoundaryCondition *boundaryConditionAtId;
    
    // Returns a boundaryCondition object at index bc_id
    bc_id = [self getBoundaryConditionID:model forElement:element];
    boundaryConditionAtId = (model.boundaries)[bc_id-1];
    
    if (boundaryConditionAtId != nil) {
        return boundaryConditionAtId.valuesList;
    } else {
        return nil;
    }
    
}

-(Element_t *)getBoundaryElement:(FEMSolution *)solution atIndex:(int)index {
    
    Element_t *elements;
    
    elements = solution.mesh.getElements;
    
    if (index < 0 || index > solution.mesh.numberOfBoundaryElements-1) {   
        
        errorfunct("getBoundaryElement", "Invalid element number requested at index:");
        printf("%d\n", index);
    }
    
    return &elements[solution.mesh.numberOfBulkElements+index];
}

-(int **)getEdgeMap:(int)elementFamily {
    
    int **edgeMap;
    
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
                edgeMap[1][0] = 2;
                
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

-(BOOL)isActiveElement:(Element_t *)element inSolution:(FEMSolution *)solution {
    
    BOOL l;
    int i, n;
    int *perm;
    variableArraysContainer *varContainers;
    
    memset( _indexStore, 0.0, (_sizeIndexStore*sizeof(_indexStore)) );
    n = [self getElementDofs:solution forElement:element atIndexes:_indexStore];
    varContainers = solution.variable.getContainers;
    
    perm = intvec(0, n-1);
    for (i=0; i<n; i++) {
        perm[i] = varContainers->Perm[_indexStore[i]];
    }
    
    if (all(perm, '>', 0, n) == 1) {
        l = YES;
    } else {
        l = NO;
    }

    free_ivector(perm, 0, n-1);
    varContainers = NULL;
    
    return l;
    
}

-(void)checkNormalTangentiality:(FEMModel *)model inSolution:(FEMSolution *)solution forElementNumber:(int)elno numberofNodes:(int)n atIndexes:(int *)indexes atBoundary:(int)bc variableName:(NSMutableString *)name orderOfDofs:(int)dof activeCondition:(BOOL)conditional conditionName:(NSString *)condName permutationOffset:(int)offset {
    
    int i, j, k, m, dim;
    listBuffer condition = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    double *rotvec;
    BOOL stat, all;
    FEMListUtilities *listUtil;
    FEMBoundaryCondition *boundaryConditionAtId;
    variableArraysContainer *varContainers;
    solutionArraysContainer *solContainers;
    
    listUtil = [[FEMListUtilities alloc] init];
    
    rotvec = doublevec(0, 2);
    
    dim = [model dimension];
    
    boundaryConditionAtId = (model.boundaries)[bc];
    
    if (dof < 0) return;
    
    solContainers = solution.getContainers;
    
    all = YES;
    for (i=0; i<n; i++) {
        k = indexes[i];
        if (solContainers->boundaryReorder[k] < 0) {
            continue;
        } else {
            all = NO;
            break;
        }
    }
    if (all == YES) {
        solContainers = NULL;
        return;
    }
    
    varContainers = solution.variable.getContainers;
    
    if ([listUtil listCheckPresentVariable:name inArray:boundaryConditionAtId.valuesList] == NO) return;
    if ([listUtil listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:[solution normalTangentialName] info:&stat] == YES) return;
    
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
            m = solContainers->boundaryReorder[indexes[j]];
            if (m >= 0) {
                memset( rotvec, 0.0, (3*sizeof(rotvec)) );
                rotvec[dof] = 1.0;
                [self FEMKernel_rotateNTSystem:solution :rotvec :indexes[j]];
                for (k=0; k<dim; k++) {
                    if (fabs(rotvec[k]) > 1.0e-8) solContainers->ntElement[m][k] = elno;
                }
            }
        }
    }
    
    if (condition.vector != NULL) {
        free_dvector(condition.vector, 0, condition.m-1);
        condition.vector = NULL;
    }
    free_dvector(rotvec, 0, 2);
    varContainers = NULL;
    solContainers = NULL;
}

-(void)getBoundaryIndexes:(FEMMesh *)mesh forBoundaryElement:(Element_t *)element withParentElement:(Element_t *)parent resultVector:(int *)indexes resultSize:(int)indSize {
/****************************************************************************************************************************************************
    Calculate global indexes of boundary dofs for given element and its boundary
 
    Arguments:
 
        FEMMesh *mesh      -> class containing the mesh and its edges and faces of elements
        Element_t *element -> Boundary element to get indexes for
        Element_t *parent  -> Parent of boundary element
        int *indexes       -> Calculated indexes of boundary element in global system
        int resultSize     -> Size of created index vector, i.e., how many indexes were created
 
****************************************************************************************************************************************************/
    
    int i, j, n;
    Element_t *edges, *faces;
    
    // Clear indexes
    memset( indexes, 0.0, (mesh.numberOfNodes*sizeof(indexes)) );
    
    n = element->Type.NumberOfNodes;
    
    edges = [mesh getEdges];
    faces = [mesh getFaces];
    
    // Nodal indexes
    memcpy(indexes, element->NodeIndexes, n*sizeof(indexes));
    
    // Assign rest of indexes if necessary
    switch (parent->Type.dimension) {
        case 2:
            // Add index for each bubble dof in edge
            for (i=0; i<element->BDOFs; i++) {
                if (mesh.numberOfNodes < n) {
                    errorfunct("getBoundaryIndexes", "Not enough space reserved for indexes.");
                    return;
                }
                
                indexes[n] = mesh.numberOfNodes + parent->EdgeIndexes[element->Pdefs->LocalNumber-1] * mesh.maxEdgeDofs + i;
                n++;
            }
            indSize = n;
            break;
        case 3:            
            // Add indexes of faces edges
            for (i=0; i<faces[parent->FaceIndexes[element->Pdefs->LocalNumber-1]].Type.NumberOfEdges; i++) {
                
                // If edge has no dofs jump to next edge
                if (edges[faces[parent->FaceIndexes[element->Pdefs->LocalNumber-1]].EdgeIndexes[i]].BDOFs <= 0) continue;
                
                for (j=0; j<edges[faces[parent->FaceIndexes[element->Pdefs->LocalNumber-1]].EdgeIndexes[i]].BDOFs; j++) {
                    
                    if (mesh.numberOfNodes < n) {
                        errorfunct("getBoundaryIndexes", "Not enough space reserved for indexes.");
                        return;
                    }
                    
                    indexes[n] = mesh.numberOfNodes + faces[parent->FaceIndexes[element->Pdefs->LocalNumber-1]].EdgeIndexes[i]*mesh.maxEdgeDofs + j;
                    n++;
                }
            }
            
            // Add indexes of gaces bubbles
            for (i=0; i<faces[parent->FaceIndexes[element->Pdefs->LocalNumber-1]].BDOFs; i++) {                
                if (mesh.numberOfNodes < n) {
                    errorfunct("getBoundaryIndexes", "Not enough space reserved for indexes.");
                    return;
                }
                
                indexes[n] = mesh.numberOfNodes + mesh.numberOfEdges * mesh.maxEdgeDofs + parent->FaceIndexes[element->Pdefs->LocalNumber-1] * mesh.maxFaceDofs;
                n++;
            }
            
            indSize = n;
            faces = NULL;
            edges = NULL;
            break;
        default:
            errorfunct("getBoundaryIndexes", "Unsupported dimension.");
            break;
    }
    
}

-(void)zeroTheNumberOfRows:(int)n inSolutionMatrix:(FEMSolution *)solution {
    
    FEMMatrixCRS *crsMatrix;
    FEMMatrixBand *bandMatrix;
    
    if (solution.matrix.format == MATRIX_CRS) {
        crsMatrix = [[FEMMatrixCRS alloc] init];
        [crsMatrix zeroRowInGlobal:solution :n];
        
    } else if (solution.matrix.format == MATRIX_LIST) {
        
        // TODO: implement the zeroRow method for list matrix.
        
    } else if (solution.matrix.format == MATRIX_BAND || solution.matrix.format == MATRIX_SBAND) {
        bandMatrix = [[FEMMatrixBand alloc] init];
        [bandMatrix zeroRowInGlobal:solution :n];
    }
    
}

-(void)setMatrixElement:(FEMSolution *)solution: (int)i: (int)j: (double)value {
    
    FEMMatrixCRS *crsMatrix;
    FEMMatrixBand *bandMatrix;
    
    if (solution.matrix.format == MATRIX_CRS) {
         crsMatrix = [[FEMMatrixCRS alloc] init];
        [crsMatrix setMatrixElementInGlobal:solution :i :j :value];
        
    } else if (solution.matrix.format == MATRIX_LIST) {
        // TODO: implement the setMatrixElement method for list matrix.
        
    } else if (solution.matrix.format == MATRIX_BAND || solution.matrix.format == MATRIX_SBAND) {
        bandMatrix = [[FEMMatrixBand alloc] init];
        [bandMatrix setMatrixElementInGlobal:solution :i :j :value];
    }
    
}

-(void)addToMatrixElement:(FEMSolution *)solution: (int)i: (int)j: (double)value {
    
    FEMMatrixCRS *crsMatrix;
    FEMMatrixBand *bandMatrix;
    
    if (solution.matrix.format == MATRIX_CRS) {
        crsMatrix = [[FEMMatrixCRS alloc] init];
        [crsMatrix addToMatrixElementInGlobal:solution :i :j :value];
        
    } else if (solution.matrix.format == MATRIX_LIST) {
        // TODO: implement the setMatrixElement method for list matrix.
        
    } else if (solution.matrix.format == MATRIX_BAND || solution.matrix.format == MATRIX_SBAND) {
        bandMatrix = [[FEMMatrixBand alloc] init];
        [bandMatrix addToMatrixElementInGlobal:solution :i :j :value];
    }

}

-(void)localBoundaryIntegral:(FEMModel *)model inSolution:(FEMSolution *)solution atBoundary:(NSArray *)bc forElement:(Element_t *)element withNumberOfNodes:(int)nd andParent:(Element_t *)parent withNumberOfNodes:(int)np boundaryName:(NSMutableString *)name functionIntegral:(double)integral {
    
    int i, j, n, jj, kk, t, size;
    int **edgeMap;
    double s, l, sum;
    double **vLoad, *g, *vl;
    listBuffer load = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    listBuffer buffer = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    FEMElementDescription *elementDescription;
    FEMNumericIntegration *numericIntegration;
    ElementType_t savedType;
    GaussIntegrationPoints *integCompound;
    Nodes_t *nodes, *pNodes;
    NSMutableString *string;
    BOOL stat;
    
    n = max(solution.mesh.maxElementNodes, solution.mesh.maxElementDofs);
    
    elementDescription = [[FEMElementDescription alloc] init];
    
    numericIntegration = [[FEMNumericIntegration alloc] init];
    [numericIntegration allocation:solution.mesh];
    
    nodes = (Nodes_t*)malloc(sizeof(Nodes_t) * n );
    pNodes = (Nodes_t*)malloc(sizeof(Nodes_t) * n );
    
    [self getNodes:solution inElement:element resultNodes:nodes numberOfNodes:n];
    [self getNodes:solution inElement:parent resultNodes:pNodes numberOfNodes:n];
    
    vLoad = doublematrix(0, 2, 0, np-1);
    for (i=0; i<3; i++) {
        for (j=0; j<np; j++) {
            vLoad[i][j] = 0.0;
        }
    }
    
    vl = doublevec(0, 2);
    memset( vl, 0.0, (3*sizeof(vl)) );
    
    g = doublevec(0, 2);
    memset( g, 0.0, (3*sizeof(g)) );
    
    [self getReal:model forElement:element inList:bc variableName:name buffer:&load];
    
    edgeMap = [self getEdgeMap:[self getElementFamily:parent]];
    
    switch ([self getElementFamily:parent]) {
        case 2:
            size = 1;
            break;
        case 3:
            size = 3;
            break;
        case 4:
            size = 4;
            break;
        case 5:
            size = 6;
            break;
        case 6:
            size = 8;
            break;
        case 7:
            size = 9;
            break;
        case 8:
            size = 12;
            break;
    }
    
    for (i=0; i<size; i++) {
        jj = edgeMap[i][0];
        kk = edgeMap[i][1];
        if ( (parent->NodeIndexes[jj] == element->NodeIndexes[0] && parent->NodeIndexes[kk] == element->NodeIndexes[1]) ||
            (parent->NodeIndexes[jj] == element->NodeIndexes[1] && parent->NodeIndexes[kk] == element->NodeIndexes[0]) ) break;
    }
    
    string = [string initWithString:name];
    [string appendString:@" 1"];
    [self getReal:model forElement:element inList:bc variableName:string buffer:&buffer];
    for (i=0; i<nd; i++) {
        vLoad[0][i] = buffer.vector[i];
    }
    if (buffer.vector != NULL) {
        free_dvector(buffer.vector, 0, buffer.m-1);
        buffer.vector = NULL;
    }
    
    string = [string initWithString:name];
    [string appendString:@" 2"];
    [self getReal:model forElement:element inList:bc variableName:string buffer:&buffer];
    for (i=0; i<nd; i++) {
        vLoad[1][i] = buffer.vector[i];
    }
    if (buffer.vector != NULL) {
        free_dvector(buffer.vector, 0, buffer.m-1);
        buffer.vector = NULL;
    }
    
    string = [string initWithString:name];
    [string appendString:@" 3"];
    [self getReal:model forElement:element inList:bc variableName:string buffer:&buffer];
    for (i=0; i<nd; i++) {
        vLoad[2][i] = buffer.vector[i];
    }
    if (buffer.vector != NULL) {
        free_dvector(buffer.vector, 0, buffer.m-1);
        buffer.vector = NULL;
    }

    g[0] = pNodes[kk].x - pNodes[jj].x;
    g[1] = pNodes[kk].y - pNodes[jj].y;
    g[2] = pNodes[kk].z - pNodes[jj].z;
    sum = 0.0;
    for (i=0; i<3; i++) {
        sum = sum + pow(g[i], 2.0);
    }
    for (i=0; i<3; i++) {
        g[i] = g[i] / sqrt(sum);
    }
    
    savedType = element->Type;
    if ([self getElementFamily:element] == 1) element->Type = *[elementDescription getElementType:202 inMesh:solution.mesh stabilization:NULL];
    
    integral = 0.0;
    integCompound = GaussQuadrature(element);
    for (t=0; t<integCompound->n; t++) {
        stat = [numericIntegration setMetricDeterminantForElement:element 
                                                     elementNodes:nodes 
                                                           inMesh:solution.mesh 
                                             firstEvaluationPoint:integCompound->u[t]
                                            secondEvaluationPoint:integCompound->v[t]
                                             thirdEvaluationPoint:integCompound->w[t]];
        s = integCompound->s[t] * [numericIntegration metricDeterminant];
        
        stat = [numericIntegration setBasisForElement:element 
                                         elementNodes:nodes 
                                               inMesh:solution.mesh 
                                 firstEvaluationPoint:integCompound->u[t]
                                secondEvaluationPoint:integCompound->v[t]
                                 thirdEvaluationPoint:integCompound->w[t]
                                          withBubbles:NO 
                                          basisDegree:NULL];
        sum = 0.0;
        for (i=0; i<nd; i++) {
            sum = sum + load.vector[i]*numericIntegration.basis[i];
        }
        l = sum;
        cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, nd, 1.0, *vLoad, np, numericIntegration.basis, 1, 0.0, vl, 1);
        sum = 0.0;
        for (i=0; i<3; i++) {
            sum = sum + vl[i] * g[i];
        }
        integral = integral + s * (l+sum);
    }
    element->Type = savedType;
    
    jj = parent->NodeIndexes[jj];
    kk = parent->NodeIndexes[kk];
    
    if (kk < jj) integral = -integral;
    
    edgeMap = NULL;
    
    free(nodes);
    free(pNodes);
    if (load.vector != NULL) {
        free_dvector(load.vector, 0, load.m-1);
    }
    free_dmatrix(vLoad, 0, 2, 0, np-1);
    free_dvector(vl, 0, 2);
    free_dvector(g, 0, 2);
    free(integCompound);
    
    [numericIntegration deallocation:solution.mesh];
    [elementDescription deallocation];
}

-(void)localBoundaryBDOFs:(FEMModel *)model inSolution:(FEMSolution *)solution atBoundary:(NSArray *)bc forElement:(Element_t *)element withNumberOfNodes:(int)nd boundaryName:(NSMutableString *)name resultMatrix:(double **)stiff resultVector:(double *)force {
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
    
    GaussIntegrationPoints *integCompound;
    int i, j, n, p, q, t;
    double xip, yip, zip, s, load;
    FEMNumericIntegration *numericIntegration;
    FEMListUtilities *listUtil;
    Nodes_t *nodes;
    BOOL stat;
    
    numericIntegration = [[FEMNumericIntegration alloc] init];
    [numericIntegration allocation:solution.mesh];
    
    listUtil = [[FEMListUtilities alloc] init];
    
    n = max(solution.mesh.maxElementNodes, solution.mesh.maxElementDofs);
    
    // Get nodes of boundary elements parent and gauss points for boundary
    nodes = (Nodes_t*)malloc(sizeof(Nodes_t) * n );
    [self getNodes:solution inElement:element resultNodes:nodes numberOfNodes:n];
    integCompound = GaussQuadrature(element);
    
    memset( force, 0.0, (nd*sizeof(force)) );
    for (i=0; i<nd; i++) {
        for (j=0; j<nd; j++) {
            stiff[i][j] = 0.0;
        }
    }
    
    for (t=0; t<integCompound->n; t++) {
        stat = [numericIntegration setMetricDeterminantForElement:element 
                                                     elementNodes:nodes 
                                                           inMesh:solution.mesh 
                                             firstEvaluationPoint:integCompound->u[t]
                                            secondEvaluationPoint:integCompound->v[t]
                                             thirdEvaluationPoint:integCompound->w[t]];
        
        s = integCompound->s[t] * [numericIntegration metricDeterminant];
        
        stat = [numericIntegration setBasisForElement:element 
                                         elementNodes:nodes 
                                               inMesh:solution.mesh 
                                 firstEvaluationPoint:integCompound->u[t]
                                secondEvaluationPoint:integCompound->v[t]
                                 thirdEvaluationPoint:integCompound->w[t]
                                          withBubbles:NO 
                                          basisDegree:NULL];
        
        // Get value of boundary condition
        xip = 0.0;
        yip = 0.0;
        zip = 0.0;
        for (i=0; i<nd; i++) {
            xip = xip + numericIntegration.basis[i]*nodes[i].x;
            yip = yip + numericIntegration.basis[i]*nodes[i].y;
            zip = zip + numericIntegration.basis[i]*nodes[i].z;
        }
        load = [listUtil listGetConstReal:model inArray:bc forVariable:name info:&stat minValue:NULL maxValue:NULL];
        
        // Build local stiffness matrix and force vector
        for (p=0; p<nd; p++) {
            for (q=0; q<nd; q++) {
                stiff[p][q] = stiff[p][q] + s * numericIntegration.basis[p] * numericIntegration.basis[q];
            }
            force[p] = force[p] + s * load * numericIntegration.basis[p];
        }
    }
    
    free(nodes);
    free(integCompound);
    [numericIntegration deallocation:solution.mesh];
}

-(void)solveWithLapackMatrix:(double *)a andVector:(double *)x size:(int)n leadingDimension:(int)lda {
    
    int nhrs, info;
    int *ipiv;
    char *trans;
    
    ipiv = intvec(0, n-1);
    
    if (n <= 0) return;
    dgetrf_(&n, &n, a, &lda, ipiv, &info);
    if (info < 0 || info > 0) {
        warnfunct("solveWithLapackMatrix", "Error in lapack routine dgetrf. Error code:");
        printf("%d\n", info);
        errorfunct("solveWithLapackMatrix", "Program terminating now...");
    }
    
    trans = "N";
    nhrs = 1;
    dgetrs_(trans, &n, &nhrs, a, &lda, ipiv, x, &n, &info);
    if (info < 0 || info > 0) {
        warnfunct("solveWithLapackMatrix", "Error in lapack routine dgetrs. Error code:");
        printf("%d\n", info);
        errorfunct("solveWithLapackMatrix", "Program terminating now...");
    }
    
    free_ivector(ipiv, 0, n-1);
}

-(void)solveLinearSystemWithMatrix:(double **)a andVector:(double *)x size:(int)n leadingDimension:(int)lda {
    
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
            [util solveLinearSystem2x2:a :x :b];
            break;
        case 3:
            for (i=0; i<n; i++) {
                b[i] = x[i];
            }
            [util solveLinearSystem3x3:a :x :b];
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

-(void)setNodalLoads:(FEMModel *)model inSolution:(FEMSolution *)solution variableName:(NSMutableString *)name orderOfDofs:(int)dof {
/*****************************************************************************************************************
    Set nodel load for given dof
 
    FEMModel *model           -> class containing the model
    FEMSolution *solution     -> solution class containing the matrix and variable
    NSMutableString *name     -> name of the dof to be set
    int dof                   -> the order number of the dof
 
    The permutation (node reordoring info) is contained in the class solution and has been generated at the 
    beginning of the simulation by the bandwidth optimization
 
*****************************************************************************************************************/
    
    int i, j, n, bc, t, bf_id, noNodes, noDims;
    int *indexes, *inNodes;
    double minDist, dist;
    listBuffer nodeIndexes = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    listBuffer coordNodes = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    BOOL *activePort, *activePortAll, *doneLoad, anyActive;
    Element_t *elements;
    Nodes_t *globalNodes;
    FEMListUtilities *listUtil;
    FEMBoundaryCondition *boundaryConditionAtId;
    FEMBodyForce *bodyForceAtId;
    NSMutableString *loadName, *str;
    matrixArraysContainer *matContainers;
    variableArraysContainer *varContainers;
    BOOL stat, nodesFound;
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
    
    listUtil = [[FEMListUtilities alloc] init];
    
    loadName = [NSMutableString stringWithString:name];
    [loadName appendString:@" load"];
    
    n = max(model.numberOfBoundaries, model.numberOfBodyForces);
    activePort = (BOOL*)malloc(sizeof(BOOL) * n );
    activePortAll = (BOOL*)malloc(sizeof(BOOL) * n );
    
    elements = solution.mesh.getElements;
    globalNodes = solution.mesh.getNodes;
    
    indexes = intvec(0, solution.mesh.maxElementDofs-1);
    
    // Go through the boundaries
    
    memset( activePort, NO, (n*sizeof(activePort)) );
    memset( activePortAll, NO, (n*sizeof(activePortAll)) );
    
    str = [NSMutableString stringWithString:loadName];
    [str appendString:@" dofs"];
    
    for (bc=0; bc<model.numberOfBoundaries; bc++) {
        boundaryConditionAtId = (model.boundaries)[bc];
        if ([listUtil listCheckPresentVariable:@"Target boundaries" inArray:boundaryConditionAtId.valuesList] == NO) continue;
        activePort[bc] = [listUtil listCheckPresentVariable:loadName inArray:boundaryConditionAtId.valuesList];
        activePortAll[bc] = [listUtil listCheckPresentVariable:str inArray:boundaryConditionAtId.valuesList];
    }
    
    anyActive = NO;
    for (bc=0; bc<model.numberOfBoundaries; bc++) {
        if (activePort[bc] == YES || activePortAll[bc] == YES) {
            anyActive = YES;
            break;
        }
    }
    
    if (anyActive == YES) {
        doneLoad = (BOOL*)malloc(sizeof(BOOL) *  (matContainers->sizeRHS/solution.variable.dofs) );
        for (i=0; i<(matContainers->sizeRHS/solution.variable.dofs); i++) {
            doneLoad[i] = NO;
        }
        
        for (bc=0; bc<model.numberOfBoundaries; bc++) {
            if (activePort[bc] == NO && activePortAll[bc] == NO) continue;
            
              
            for (t=[model numberOfBulkElements]; t<[model numberOfBulkElements]+[model numberOfBoundaryElements]; t++) {
                
                boundaryConditionAtId = (model.boundaries)[bc];
                if (elements[t].BoundaryInfo->Constraint != [boundaryConditionAtId tag]) continue;
                
                if (activePort[bc] == YES) {
                    n = elements[t].Type.NumberOfNodes;
                    for (i=0; i<n; i++) {
                        indexes[i] = elements[t].NodeIndexes[i];
                    }
                } else {
                    n = [self sgetElementDofs:solution forElement:&elements[t] atIndexes:indexes];
                }
                
                [self FEMKernel_setElementLoads:model :solution :&elements[t] :boundaryConditionAtId.valuesList :loadName :indexes :doneLoad :n :dof :solution.variable.dofs];
            }
        }
    }
    
    // Go though the nodal load conditions for the body force list
        
    memset( activePort, NO, (n*sizeof(activePort)) );
    memset( activePortAll, NO, (n*sizeof(activePortAll)) );
    
    for (bf_id=0; bf_id<[model numberOfBodyForces]; bf_id++) {
        bodyForceAtId = (model.bodyForces)[bf_id];
        activePort[bf_id] = [listUtil listCheckPresentVariable:loadName inArray:bodyForceAtId.valuesList];
        activePortAll[bf_id] = [listUtil listCheckPresentVariable:str inArray:bodyForceAtId.valuesList];
    }
    
    anyActive = NO;
    for (bf_id=0; bf_id<[model numberOfBodyForces]; bf_id++) {
        if (activePort[bf_id] == YES || activePortAll[bf_id] == YES) {
            anyActive = YES;
            break;
        }
    }
    
        
    if (anyActive == YES) {
        if (doneLoad == NULL) doneLoad = (BOOL*)malloc(sizeof(BOOL) *  (matContainers->sizeRHS/solution.variable.dofs) );
        for (i=0; i<(matContainers->sizeRHS/solution.variable.dofs); i++) {
            doneLoad[i] = NO;
        }
        for (t=0; t<[model numberOfBulkElements]; t++) {
            
            bf_id = [(model.bodies)[elements[t].BodyID-1][@"Body force"] intValue];
            
            if ((model.bodies)[elements[t].BodyID-1][@"Body force"] == nil) continue;
            if (activePort[bf_id] == NO && activePortAll[bf_id] == NO) continue;
            
            
            if (activePort[bf_id] == YES) {
                n = elements[t].Type.NumberOfNodes;
                for (i=0; i<n; i++) {
                    indexes[i] = elements[t].NodeIndexes[i];
                }
            } else {
                n = [self sgetElementDofs:solution forElement:&elements[t] atIndexes:indexes];
            }
            
            bodyForceAtId = (model.bodyForces)[bf_id];
            [self FEMKernel_setElementLoads:model :solution :&elements[t] :bodyForceAtId.valuesList :loadName :indexes :doneLoad :n :dof :solution.variable.dofs];
        }
    }
    
    free(activePort);
    free(activePortAll);
    if (doneLoad != NULL) free(doneLoad);
    
    // Go through the point loads which are created on the fly
    for (bc=0; bc<model.numberOfBoundaries; bc++) {
        
        boundaryConditionAtId = (model.boundaries)[bc];
        if ([listUtil listCheckPresentVariable:loadName inArray:boundaryConditionAtId.valuesList] == NO) continue;
        nodesFound = [listUtil listCheckPresentVariable:@"Target nodes" inArray:boundaryConditionAtId.valuesList];
        
        // At the first calling, the list of coordinates is transformed to list of nodes
        if (nodesFound == NO) {
            
            stat = [listUtil listGetConstRealArray:model inArray:boundaryConditionAtId.valuesList forVariable:@"Target coordinates" buffer:&coordNodes];
            
            if (stat == YES) {
                
                noNodes = coordNodes.m;
                noDims = coordNodes.n;
                
                if (noNodes > 0) {
                    inNodes = intvec(0, noNodes-1);
                    for (j=0; j<noNodes; j++) {
                        minDist = HUGE_VAL;
                        for (i=0; i<[model numberOfNodes]; i++) {
                            if (varContainers->Perm[i] < 0) continue;
                            
                            dist = pow((globalNodes[i].x-coordNodes.matrix[j][0]), 2.0);
                            if (noDims >= 2) dist = dist + pow((globalNodes[i].y-coordNodes.matrix[j][1]), 2.0);
                            if (noDims == 3) dist = dist + pow((globalNodes[i].z-coordNodes.matrix[j][2]), 2.0);
                            
                            if (dist<minDist) {
                                minDist = dist;
                                inNodes[j] = i;
                            }
                        }
                    }
                    // Add the found nodes to the list values
                    [listUtil addIntegerArrayInClassList:boundaryConditionAtId theVariable:@"Target nodes" withValues:inNodes numberOfNodes:noNodes];
                    free_ivector(inNodes, 0, noNodes-1);
                    nodesFound = YES;
                }
            }
        }
        
        if (coordNodes.matrix != NULL) {
            free_dmatrix(coordNodes.matrix, 0, coordNodes.m-1, 0, coordNodes.n-1);
            coordNodes.matrix = NULL;
        }
        
        if (nodesFound == YES) {
            [listUtil listGetIntegerArray:model inArray:boundaryConditionAtId.valuesList forVariable:@"Target nodes" buffer:&nodeIndexes];
            
            [self FEMKernel_setPointLoads:model :solution :elements :boundaryConditionAtId.valuesList :loadName :nodeIndexes.ivector :n :dof :solution.variable.dofs];
            if (nodeIndexes.ivector != NULL) {
                free_ivector(nodeIndexes.ivector, 0, nodeIndexes.m-1);
                nodeIndexes.ivector = NULL;
            }
        }
    }
    
    free_ivector(indexes, 0, solution.mesh.maxElementDofs-1);
    
    elements = NULL;
    globalNodes = NULL;
    matContainers = NULL;
    varContainers = NULL;
}

-(void)setDirichletBoundaries:(FEMModel *)model inSolution:(FEMSolution *)solution variableName:(NSMutableString *)name orderOfDofs:(int)dof permutationOffset:(int *)offset {
/*****************************************************************************************************************
    Set dirichlet boundary condition for given dof
 
    FEMModel *model           -> class containing the model
    FEMSolution *solution     -> solution class containing the matrix and variable
    NSMutableString *name     -> name of the dof to be set
    int dof                   -> the order number of the dof
    int *offset               -> Optional: If the matrix and pernutation vectors are not in sync the offset
                                 may be used as a remedy. Needed in fully coupled systems.
 
    The permutation (node reordoring info) is contained in the class solution and has been generated at the 
    beginning of the simulation by the bandwidth optimization
 
*****************************************************************************************************************/
    
    int i, j, k, n, bc, t, bf_id, noNodes, noDims, permOffset, numberOfNodesFound;
    int *indexes, *inNodes, bndry_start, bndry_end;
    double minDist, dist, eps;
    listBuffer nodeIndexes = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    listBuffer coordNodes = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    BOOL *activePort, *activePortAll, *activeCond, *donePeriodic, anyActive, passive;
    Element_t *elements;
    Nodes_t *globalNodes;
    FEMListUtilities *listUtil;
    FEMBoundaryCondition *boundaryConditionAtId;
    FEMBodyForce *bodyForceAtId;
    NSMutableString *condName, *passName, *str1, *str2;
    matrixArraysContainer *matContainers;
    variableArraysContainer *varContainers;
    solutionArraysContainer *solContainers;
    valueListArraysContainer *valContainers;
    BOOL stat, nodesFound, orderByBCNumbering, conditional;
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
    solContainers = solution.getContainers;
    
    listUtil = [[FEMListUtilities alloc] init];
    
    // These logical vectors are used to minimize extra effort in setting up different BCs
    n = max(model.numberOfBoundaries, model.numberOfBodyForces);
    activePort = (BOOL*)malloc(sizeof(BOOL) * n );
    activePortAll = (BOOL*)malloc(sizeof(BOOL) * n );
    activeCond = (BOOL*)malloc(sizeof(BOOL) * n );
    
    condName = [NSMutableString stringWithString:name];
    [condName appendString:@" condition"];
    
    passName = [NSMutableString stringWithString:name];
    [passName appendString:@" passive"];
    
    permOffset = 0;
    if (offset != NULL) permOffset = *offset;
    
    elements = solution.mesh.getElements;
    globalNodes = solution.mesh.getNodes;

    indexes = intvec(0, solution.mesh.maxElementDofs-1);
    
    // Go through the perdiodic BCs and set the linear dependence
    
    memset( activePort, NO, (n*sizeof(activePort)) );
    
    str1 = [NSMutableString stringWithString:@"Periodic BC "];
    [str1 appendString:name];
    
    str2 = [NSMutableString stringWithString:@"Anti periodic BC "];
    [str2 appendString:name];

    for (bc=0; bc<model.numberOfBoundaries; bc++) {
        boundaryConditionAtId = (model.boundaries)[bc];
        if ([listUtil listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:str1 info:&stat] == NO) activePort[bc] = YES;
        if ([listUtil listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:str2 info:&stat] == NO) activePort[bc] = YES;
    }
    
    anyActive = NO;
    for (bc=0; bc<model.numberOfBoundaries; bc++) {
        if (activePort[bc] == YES) {
            anyActive = YES;
            break;
        }
    }
    
    if (anyActive == YES) {
        
        donePeriodic = (BOOL*)malloc(sizeof(BOOL) * solution.mesh.numberOfNodes );
        memset( donePeriodic, NO, (solution.mesh.numberOfNodes*sizeof(donePeriodic)) );
        
        for (bc=0; bc<model.numberOfBoundaries; bc++) {
            [self FEMKernel_setPeriodicBoundariesPass1:model :solution :name :dof :bc :donePeriodic];
            
        }
        
        memset( donePeriodic, NO, (solution.mesh.numberOfNodes*sizeof(donePeriodic)) );
        for (bc=0; bc<model.numberOfBoundaries; bc++) {
            [self FEMKernel_setPeriodicBoundariesPass2:model :solution :name :dof :bc :donePeriodic];
            
        }
        
        free(donePeriodic);
        
    }
    
    // Go through the normal Dirichlet BCs applied on the boundaries
    
    memset( activePort, NO, (n*sizeof(activePort)) );
    memset( activePortAll, NO, (n*sizeof(activePortAll)) );
    memset( activeCond, NO, (n*sizeof(activeCond)) );
    
    str1 = [NSMutableString stringWithString:name];
    [str1 appendString:@" dofs"];
    
    for (bc=0; bc<model.numberOfBoundaries; bc++) {
        boundaryConditionAtId = (model.boundaries)[bc];
        activePortAll[bc] = [listUtil listCheckPresentVariable:str1 inArray:boundaryConditionAtId.valuesList];
        activePort[bc] = [listUtil listCheckPresentVariable:name inArray:boundaryConditionAtId.valuesList];
        activeCond[bc] = [listUtil listCheckPresentVariable:condName inArray:boundaryConditionAtId.valuesList];
    }
    
    orderByBCNumbering = [listUtil listGetLogical:model inArray:model.simulation.valuesList forVariable:@"Set Dirichlet BCs by BC numbering" info:&stat];
    
    bndry_start = [model numberOfBulkElements];
    bndry_end = bndry_start + [model numberOfBoundaryElements]-1;
    
    // Check and sey some flags for nodes belonging to n-t boundaries getting set by other bcs
    if ([solution normalTangentialNOFNodes] > 0) {
        if (orderByBCNumbering == YES) {
            for (bc=0; bc<model.numberOfBoundaries; bc++) {
                
                if (activePort[bc] == NO && activePortAll[bc] == NO) continue;
                conditional = activeCond[bc];
                
                boundaryConditionAtId = (model.boundaries)[bc];
                
                for (t=bndry_start; t<bndry_end; t++) {
                    if (elements[t].BoundaryInfo->Constraint != [boundaryConditionAtId tag]) continue;
                    if (activePort[bc] == YES) {
                        n = elements[t].Type.NumberOfNodes;
                        for (i=0; i<n; i++) {
                            indexes[i] = elements[t].NodeIndexes[i];
                        }
                    } else {
                        n = [self sgetElementDofs:solution forElement:&elements[t] atIndexes:indexes];
                    }
                    [self checkNormalTangentiality:model inSolution:solution forElementNumber:t numberofNodes:n atIndexes:indexes atBoundary:bc variableName:name orderOfDofs:dof activeCondition:conditional conditionName:condName permutationOffset:permOffset];
                }
            }
        } else {
             for (t=bndry_start; t<bndry_end; t++) {
                 for (bc=0; bc<model.numberOfBoundaries; bc++) {
                     if (activePort[bc] == NO && activePortAll[bc] == NO) continue;
                     conditional = activeCond[bc];
                     
                     boundaryConditionAtId = (model.boundaries)[bc];
                     
                     if (elements[t].BoundaryInfo->Constraint != [boundaryConditionAtId tag]) continue;
                     
                     if (activePort[bc] == YES) {
                         n = elements[t].Type.NumberOfNodes;
                         for (i=0; i<n; i++) {
                             indexes[i] = elements[t].NodeIndexes[i];
                         }
                     } else {
                         n = [self sgetElementDofs:solution forElement:&elements[t] atIndexes:indexes];
                     }
                     [self checkNormalTangentiality:model inSolution:solution forElementNumber:t numberofNodes:n atIndexes:indexes atBoundary:bc variableName:name orderOfDofs:dof activeCondition:conditional conditionName:condName permutationOffset:permOffset];
                 }
             }
        }
        
        if (dof < 0) {
            for (t=bndry_start; t<bndry_end; t++) {
                n = elements[t].Type.NumberOfNodes;
                for (j=0; j<n; j++) {
                    k = solContainers->boundaryReorder[elements[t].NodeIndexes[j]];
                    if (k >= 0) {
                        for (i=0; i<3; i++) {
                            solContainers->ntZeroingDone[k][i] = false;
                        }
                    }
                }
            }
        }
    }
    
    // Set the Dirichlet BCs from active boundary elements, if any...
    anyActive = NO;
    for (bc=0; bc<model.numberOfBoundaries; bc++) {
        if (activePort[bc] == YES || activePortAll[bc] == YES) {
            anyActive = YES;
            break;
        }
    }
    
    if (anyActive == YES) {
        if (orderByBCNumbering == YES) {
            for (bc=0; bc<model.numberOfBoundaries; bc++) {
                
                if (activePort[bc] == NO && activePortAll[bc] == NO) continue;
                conditional = activeCond[bc];
                
                boundaryConditionAtId = (model.boundaries)[bc];
                
                for (t=bndry_start; t<bndry_end; t++) {
                    if (elements[t].BoundaryInfo->Constraint != [boundaryConditionAtId tag]) continue;
                    if (activePort[bc] == YES) {
                        n = elements[t].Type.NumberOfNodes;
                        for (i=0; i<n; i++) {
                            indexes[i] = elements[t].NodeIndexes[i];
                        }
                    } else {
                        n = [self sgetElementDofs:solution forElement:&elements[t] atIndexes:indexes];
                    }
                    [self FEMKernel_setElementValues:model inSolution:solution forElementNumber:t numberOfNodes:n atIndexes:indexes forValues:boundaryConditionAtId.valuesList variableName:name orderOfDofs:dof activeCondition:conditional conditionName:condName permutationOffset:permOffset];
                }
            }
        } else {
            for (t=bndry_start; t<bndry_end; t++) {
                for (bc=0; bc<model.numberOfBoundaries; bc++) {
                    if (activePort[bc] == NO && activePortAll[bc] == NO) continue;
                    conditional = activeCond[bc];
                    
                    boundaryConditionAtId = (model.boundaries)[bc];
                    
                    if (elements[t].BoundaryInfo->Constraint != [boundaryConditionAtId tag]) continue;
                    
                    if (activePort[bc] == YES) {
                        n = elements[t].Type.NumberOfNodes;
                        for (i=0; i<n; i++) {
                            indexes[i] = elements[t].NodeIndexes[i];
                        }
                    } else {
                        n = [self sgetElementDofs:solution forElement:&elements[t] atIndexes:indexes];
                    }
                    [self FEMKernel_setElementValues:model inSolution:solution forElementNumber:t numberOfNodes:n atIndexes:indexes forValues:boundaryConditionAtId.valuesList variableName:name orderOfDofs:dof activeCondition:conditional conditionName:condName permutationOffset:permOffset];
                }
            }
        }
    }
    
    // Go through the Dirichlet conditions in the body force lists
    
    memset( activePort, NO, (n*sizeof(activePort)) );
    memset( activeCond, NO, (n*sizeof(activeCond)) );
    memset( activePortAll, NO, (n*sizeof(activePortAll)) );
    passive = NO;
    
    for (bf_id=0; bf_id<[model numberOfBodyForces]; bf_id++) {
        bodyForceAtId = (model.bodyForces)[bf_id];
        activePort[bf_id] = [listUtil listCheckPresentVariable:name inArray:bodyForceAtId.valuesList ];
        activePortAll[bf_id] = [listUtil listCheckPresentVariable:str1 inArray:bodyForceAtId.valuesList];
        activeCond[bf_id] = [listUtil listCheckPresentVariable:condName inArray:bodyForceAtId.valuesList];
        
        passive = (passive == YES || [listUtil listCheckPresentVariable:passName inArray:bodyForceAtId.valuesList] == YES) ? YES : NO;
    }
    
    anyActive = NO;
    for (bf_id=0; bf_id<[model numberOfBodyForces]; bf_id++) {
        if (activePort[bf_id] == YES || activePortAll[bf_id] == YES) {
            anyActive = YES;
            break;
        }
    }

    if (anyActive == YES) {
        for (t=0; t<[model numberOfBulkElements]; t++) {
            
            bf_id = [(model.bodies)[elements[t].BodyID-1][@"Body force"] intValue];
            
            if ((model.bodies)[elements[t].BodyID-1][@"Body force"] == nil) continue;
            if (activePort[bf_id] == NO && activePortAll[bf_id] == NO) continue;
            conditional = activeCond[bf_id];
            
            if (activePort[bf_id] == YES) {
                n = elements[t].Type.NumberOfNodes;
                for (i=0; i<n; i++) {
                    indexes[i] = elements[t].NodeIndexes[i];
                }
            } else {
                n = [self sgetElementDofs:solution forElement:&elements[t] atIndexes:indexes];
            }
            bodyForceAtId = (model.bodyForces)[bf_id];
            [self FEMKernel_setElementValues:model inSolution:solution forElementNumber:t numberOfNodes:n atIndexes:indexes forValues:bodyForceAtId.valuesList variableName:name orderOfDofs:dof activeCondition:conditional conditionName:condName permutationOffset:permOffset];
        }
    }
    free(activePort);
    free(activePortAll);
    free(activeCond);
    free_ivector(indexes, 0, solution.mesh.maxElementDofs-1);
    
    // Go through the pointwise Dirichlet BCs that are created on the fly.
    // Note that it is best that the coordinates are transformed to nodes using
    // the right variable. Otherwise, it could point to nodes that are not active

    for (bc=0; bc<model.numberOfBoundaries; bc++) {
        
        boundaryConditionAtId = (model.boundaries)[bc];
        if ([listUtil listCheckPresentVariable:name inArray:boundaryConditionAtId.valuesList] == NO) continue;
        nodesFound = [listUtil listCheckPresentVariable:@"Target nodes" inArray:boundaryConditionAtId.valuesList];
        
        // The coodinates are only requested for a body that has no list of nodes.
        // At the first calling, the list of coordinates is transformed to a list of nodes
        if (nodesFound == NO) {
            
            stat = [listUtil listGetConstRealArray:model inArray:boundaryConditionAtId.valuesList forVariable:@"Target coordinates" buffer:&coordNodes];
            
            if (stat == YES) {
                
                eps = [listUtil listGetConstReal:model inArray:boundaryConditionAtId.valuesList forVariable:@"Target coordinates eps" info:&stat minValue:NULL maxValue:NULL];
                if (stat == NO) {
                    eps = HUGE_VAL;
                } else {
                    eps = pow(eps, 2.0);
                }
                
                noNodes = valContainers->sizeFValues1;
                noDims = valContainers->sizeFValues2;
                
                if (noNodes > 0) {
                    inNodes = intvec(0, noNodes-1);
                    memset( inNodes, -1, (noNodes*sizeof(inNodes)) );
                    for (j=0; j<noNodes; j++) {
                        minDist = HUGE_VAL;
                        for (i=0; i<[model numberOfNodes]; i++) {
                            if (varContainers->Perm[i] < 0) continue;
                            
                            dist = pow((globalNodes[i].x-coordNodes.matrix[j][0]), 2.0);
                            if (noDims >= 2) dist = dist + pow((globalNodes[i].y-coordNodes.matrix[j][1]), 2.0);
                            if (noDims == 3) dist = dist + pow((globalNodes[i].z-coordNodes.matrix[j][2]), 2.0);
                            
                            dist = sqrt(dist);
                            
                            if (dist < minDist && dist <= eps) {
                                minDist = dist;
                                inNodes[j] = i;
                            }
                        }
                    }
                    
                    numberOfNodesFound = 0;
                    for (j=0; j<noNodes; j++) {
                        if (inNodes[j] >= 0) {
                            inNodes[numberOfNodesFound] = inNodes[j];
                            numberOfNodesFound++;
                        }
                    }
                    
                    // In the first time add teh found nodes to the list
                    if (numberOfNodesFound > 0) {
                        [listUtil addIntegerArrayInClassList:boundaryConditionAtId theVariable:@"Target nodes" withValues:inNodes numberOfNodes:numberOfNodesFound];
                        free_ivector(inNodes, 0, noNodes-1);
                        nodesFound = YES;
                    }
                }
            }
            
            if (coordNodes.matrix != NULL) {
                free_dmatrix(coordNodes.matrix, 0, coordNodes.m-1, 0, coordNodes.n-1);
                coordNodes.matrix = NULL;
                
            }
        }
        
        if (nodesFound == YES) {
            conditional = [listUtil listCheckPresentVariable:condName inArray:boundaryConditionAtId.valuesList];
            [listUtil listGetIntegerArray:model inArray:boundaryConditionAtId.valuesList forVariable:@"Target nodes" buffer:&nodeIndexes];
            
            [self FEMKernel_setPointValues:model inSolution:solution numberofNodes:n atIndexes:nodeIndexes.ivector forValues:boundaryConditionAtId.valuesList variableName:name orderOfDofs:dof activeCondition:conditional conditionName:condName permutationOffset:permOffset];
            
            if (nodeIndexes.ivector != NULL) {
                free_ivector(nodeIndexes.ivector, 0, nodeIndexes.m-1);
                nodeIndexes.ivector = NULL;
            }
        }
    }
    
    // Take care of the matrix entries of passive elements
    if (passive == YES) {
        for (i=0; i<solution.matrix.numberOfRows; i++) {
            if (fabs(matContainers->Values[matContainers->Diag[i]]) < 1.0e-14) {
                matContainers->Values[matContainers->Diag[i]] = 1.0;
                matContainers->RHS[i] = varContainers->Values[i];
            }
        }
    }
    
    elements = NULL;
    globalNodes = NULL;
    matContainers = NULL;
    varContainers = NULL;
    solContainers = NULL;
}

#pragma mark First order time

-(void)defaultFirstOrderTime:(FEMModel *)model inSolution:(FEMSolution *)solution forElement:(Element_t *)element realMass:(double **)mass realStiff:(double **)stiff realForce:(double *)force stiffRows:(int *)rows stiffCols:(int *)cols { // rows and cols for stiff matrix
    
    int i, n;
    double dt;
    int *perm;
    variableArraysContainer *varContainers;
    
    dt = [solution dt];
    memset( _indexStore, 0.0, (_sizeIndexStore*sizeof(_indexStore)) );
    n = [self getElementDofs:solution forElement:element atIndexes:_indexStore];
    varContainers = solution.variable.getContainers;
    
    perm = intvec(0, n-1);
    for (i=0; i<n; i++) {
        perm[i] = varContainers->Perm[_indexStore[i]];
    }
    
    [self FEMKernel_addFirstOrderTime:model :solution :element :mass :stiff :force :dt :n :solution.variable.dofs :perm :rows :cols];
    
    free_ivector(perm, 0, n-1);
    varContainers = NULL;
}

-(void)defaultFirstOrderTime:(FEMModel *)model inSolution:(FEMSolution *)solution forElement:(Element_t *)element complexMass:(_Complex double **)cmass complexStiff:(_Complex double **)cstiff complexForce:(_Complex double *)cforce stiffRows:(int *)rows stiffCols:(int *)cols {
    
    int i, j, n, dofs;
    double dt;
    double **mass, **stiff, *force;
    int *perm;
    variableArraysContainer *varContainers;
    
    dt = [solution dt];
    dofs = solution.variable.dofs;
    memset( _indexStore, 0.0, (_sizeIndexStore*sizeof(_indexStore)) );
    n = [self getElementDofs:solution forElement:element atIndexes:_indexStore];
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
    
    perm = intvec(0, n-1);
    for (i=0; i<n; i++) {
        perm[i] = varContainers->Perm[_indexStore[i]];
    }

    [self FEMKernel_addFirstOrderTime:model :solution :element :mass :stiff :force :dt :n :solution.variable.dofs :perm :rows :cols];
    
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
    varContainers = NULL;
}

#pragma mark Update equations

-(void)defaultUpdateEquations:(FEMModel *)model inSolution:(FEMSolution *)solution forElement:(Element_t *)element realStiff:(double **)stiff realForce:(double *)force stiffRows:(int *)rows stiffCols:(int *)cols requestBulkUpdate:(BOOL *)bulkUpdate {
    
    int i, n;
    int *perm;
    BOOL rotateNT, bupd;
    double *saveValues;
    matrixArraysContainer *matContainers;
    variableArraysContainer *varContainers;
    
    memset( _indexStore, 0.0, (_sizeIndexStore*sizeof(_indexStore)) );
    n = [self getElementDofs:solution forElement:element atIndexes:_indexStore];
    
    varContainers = solution.variable.getContainers;
    perm = intvec(0, n-1);
    for (i=0; i<n; i++) {
        perm[i] = varContainers->Perm[_indexStore[i]];
    }
    varContainers = NULL;
    
    [self FEMKernel_updateGlobalEquations:model :solution :element :stiff :force :n :solution.variable.dofs :perm :rows :cols :NULL: NULL];
    
    bupd = NO;
    if (bulkUpdate != NULL) {
        bupd = *bulkUpdate;
        if (bupd == NO) return;
    } else {
        if (element->BoundaryInfo != NULL) return;
    }
    if ((solution.solutionInfo)[@"Calculate loads"] != nil) {
        bupd = (bupd == YES || [(solution.solutionInfo)[@"Calculate loads"] boolValue] == YES) ? YES : NO;
    }
    
    matContainers = solution.matrix.getContainers;
    
    if (bupd ==YES) {
        
        if (matContainers->BulkRHS == NULL) {
            matContainers->BulkRHS = doublevec(0, matContainers->sizeRHS-1);
            matContainers->sizeBulkRHS = matContainers->sizeRHS;
            memset(matContainers->BulkRHS, 0.0, (matContainers->sizeBulkRHS*sizeof(matContainers->BulkRHS)) );
        }
        
        if (matContainers->BulkValues == NULL) {
            matContainers->BulkValues = doublevec(0, matContainers->sizeValues);
            matContainers->sizeBulkValues = matContainers->sizeValues;
            memset(matContainers->BulkValues, 0.0, (matContainers->sizeBulkValues*sizeof(matContainers->BulkValues)) );
        }
        
        saveValues = doublevec(0, matContainers->sizeValues-1);
        for (i=0; i<matContainers->sizeValues; i++) {
            saveValues[i] = matContainers->Values[i];
        }
        for (i=0; i<matContainers->sizeValues; i++) {
            matContainers->Values[i] = matContainers->BulkValues[i];
        }
        rotateNT = NO;
        [self FEMKernel_updateGlobalEquations:model :solution :element :stiff :force :n :solution.variable.dofs :perm :rows :cols :&rotateNT :&bupd];
        for (i=0; i<matContainers->sizeValues; i++) {
            matContainers->Values[i] = saveValues[i];
        }
        free_dvector(saveValues, 0, matContainers->sizeValues-1);
        
    }
    
    free_ivector(perm, 0, n-1);
    matContainers = NULL;
}

-(void)defaultUpdateEquations:(FEMModel *)model inSolution:(FEMSolution *)solution forElement:(Element_t *)element complexStiff:(double **)cstiff complexForce:(double *)cforce stiffRows:(int *)rows stiffCols:(int *)cols requestBulkUpdate:(BOOL *)bulkUpdate {
    
    int i, j, n, dofs;
    int *perm;
    double **stiff, *force;
    BOOL rotateNT, bupd;
    double *saveValues;
    matrixArraysContainer *matContainers;
    variableArraysContainer *varContainers;
    
    memset( _indexStore, 0.0, (_sizeIndexStore*sizeof(_indexStore)) );
    n = [self getElementDofs:solution forElement:element atIndexes:_indexStore];
    
    varContainers = solution.variable.getContainers;
    perm = intvec(0, n-1);
    for (i=0; i<n; i++) {
        perm[i] = varContainers->Perm[_indexStore[i]];
    }
    varContainers = NULL;
    
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

    [self FEMKernel_updateGlobalEquations:model :solution :element :stiff :force :n :solution.variable.dofs :perm :rows :cols :NULL: NULL];
    
    bupd = NO;
    if (bulkUpdate != NULL) {
        bupd = *bulkUpdate;
        if (bupd == NO) return;
    } else {
        if (element->BoundaryInfo != NULL) return;
    }
    if ((solution.solutionInfo)[@"Calculate loads"] != nil) {
        bupd = (bupd == YES || [(solution.solutionInfo)[@"Calculate loads"] boolValue] == YES) ? YES : NO;
    }
    
    matContainers = solution.matrix.getContainers;
    
    if (bupd ==YES) {
        
        if (matContainers->BulkRHS == NULL) {
            matContainers->BulkRHS = doublevec(0, matContainers->sizeRHS-1);
            matContainers->sizeBulkRHS = matContainers->sizeRHS;
            memset(matContainers->BulkRHS, 0.0, (matContainers->sizeBulkRHS*sizeof(matContainers->BulkRHS)) );
        }
        
        if (matContainers->BulkValues == NULL) {
            matContainers->BulkValues = doublevec(0, matContainers->sizeValues);
            matContainers->sizeBulkValues = matContainers->sizeValues;
            memset(matContainers->BulkValues, 0.0, (matContainers->sizeBulkValues*sizeof(matContainers->BulkValues)) );
        }
        
        saveValues = doublevec(0, matContainers->sizeValues-1);
        for (i=0; i<matContainers->sizeValues; i++) {
            saveValues[i] = matContainers->Values[i];
        }
        for (i=0; i<matContainers->sizeValues; i++) {
            matContainers->Values[i] = matContainers->BulkValues[i];
        }
        rotateNT = NO;
        [self FEMKernel_updateGlobalEquations:model :solution :element :stiff :force :n :solution.variable.dofs :perm :rows :cols :&rotateNT :&bupd];
        for (i=0; i<matContainers->sizeValues; i++) {
            matContainers->Values[i] = saveValues[i];
        }
        free_dvector(saveValues, 0, matContainers->sizeValues-1);
        
    }
    
    free_dmatrix(stiff, 0, (n*dofs)-1, 0, (n*dofs)-1);
    free_dvector(force, 0, (n*dofs)-1);    
    free_ivector(perm, 0, n-1);
    matContainers = NULL;
}


#pragma mark Boundary conditions

-(void)dirichletBoundaryConditions:(FEMModel *)model inSolution:(FEMSolution *)solution usingOffset:(int *)offset {
    
    int i, j, k, kk, l, n, nb, mb, dof, numEdgeDofs, n_start, u_offset;
    Element_t *element, *parent, *edges, *faces;
    NSNumber *appendDof;
    NSMutableString *name, *componentName;
    NSMutableString *str1, *str2;
    NSArray *bc;    
    FEMMatrixCRS *crsMatrix;
    FEMListUtilities *listUtil; 
    FEMValueList *list;
    matrixArraysContainer *matContainers;
    variableArraysContainer *varContainers;
    solutionArraysContainer *solContainers;
    
    BOOL constantValue;
    
    edges = solution.mesh.getEdges;
    faces = solution.mesh.getFaces;
    
    crsMatrix = [[FEMMatrixCRS alloc] init];
    listUtil = [[FEMListUtilities alloc] init];
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
    solContainers = solution.getContainers;
    
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
            errorfunct("dirichletBoundaryConditions", "Memory allocation error.");
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
            errorfunct("dirichletBoundaryConditions", "Memory allocation error.");

        _size1kernStiff = n;
        _size2kernStiff = n;
        _sizekernWork = n;
        _size_g_Ind = n;
        _size_l_Ind = n;
    }
    
    name = [NSMutableString stringWithString:solution.variable.name];
        
    if (solution.variable.dofs > 1) {
        [self setDirichletBoundaries:model inSolution:solution variableName:name orderOfDofs:-1 permutationOffset:NULL];
    }
    
    constantValue = NO;
    
    // Set Dirichlet dofs for edges and faces
    for (dof=0; dof<solution.variable.dofs; dof++) {
        
        componentName = [NSMutableString stringWithString:solution.variable.name];
        if (solution.variable.dofs > 1) {
            appendDof = @(dof+1);
            [componentName appendString:@" "];
            [componentName appendString:[appendDof stringValue]];
        }
        
        // Clean BC face and edge dofs
        for (i=0; i<solution.mesh.numberOfBoundaryElements; i++) {
            
            element = [self getBoundaryElement:solution atIndex:i];
            if ([self isActiveElement:element inSolution:solution] == NO) continue;
            
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
                [self getBoundaryIndexes:solution.mesh forBoundaryElement:element withParentElement:parent resultVector:_g_Ind resultSize:numEdgeDofs];
            } else {
                continue;
            }
            
            // Contribute this boundary to global system (i.e., solve global boundary problem)
            for (k=n; k<numEdgeDofs; k++) {
                nb = varContainers->Perm[_g_Ind[k]];
                if (nb < 0) continue;
                nb = u_offset + solution.variable.dofs*nb + dof;
                if (constantValue == YES) {
                    [crsMatrix setSymmetricDirichletInGlobal:solution :nb :0.0];
                }else {
                    [self zeroTheNumberOfRows:nb inSolutionMatrix:solution];
                    matContainers->RHS[nb] = 0.0;
                }
            }
        }
    }
    
    // Set Dirichlet dofs for edges and faces
    for (dof=0; dof<solution.variable.dofs; dof++) {
        
        componentName = [NSMutableString stringWithString:solution.variable.name];
        if (solution.variable.dofs > 1) {
            appendDof = @(dof+1);
            [componentName appendString:@" "];
            [componentName appendString:[appendDof stringValue]];
        }
        
        [self setNodalLoads:model inSolution:solution variableName:componentName orderOfDofs:dof];
        
        [self setDirichletBoundaries:model inSolution:solution variableName:componentName orderOfDofs:dof permutationOffset:offset];
        
        // Dirichlet BCs for face and edge dofs
        for (i=0; i<solution.mesh.numberOfBoundaryElements; i++) {
            
            element = [self getBoundaryElement:solution atIndex:i];
            if ([self isActiveElement:element inSolution:solution] == NO) continue;
            
            bc = [self getBoundaryCondition:model forElement:element];
            if (bc == nil) continue;
            
            
            str1 = [NSString stringWithString:componentName];
            str2 = [NSString stringWithString:componentName];
            [str1 appendString:@" {e}"];
            [str2 appendString:@" {f}"];
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
                if (solution.mesh.isAssociatedEdges == YES) {
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
                            [self localBoundaryIntegral:model inSolution:solution atBoundary:bc forElement:&edges[parent->EdgeIndexes[j]] withNumberOfNodes:n andParent:parent withNumberOfNodes:nb boundaryName:str1 functionIntegral:_kernWork[0]];
                            
                            n = [self getElementDofs:solution forElement:&edges[parent->EdgeIndexes[j]] atIndexes:_g_Ind];
                            for (k=solContainers->defDofs[parent->BodyID-1][0]*edges[parent->EdgeIndexes[j]].NDOFs; k<n; k++) {
                                nb = varContainers->Perm[_g_Ind[k]];
                                if (nb < 0) continue;
                                nb = u_offset + solution.variable.dofs*nb + dof;
                                if (solution.matrix.isSymmetric == YES) {
                                    [crsMatrix setSymmetricDirichletInGlobal:solution :nb :_kernWork[0]];
                                } else {
                                    [self zeroTheNumberOfRows:nb inSolutionMatrix:solution];
                                    matContainers->RHS[nb] = _kernWork[0];
                                    [self setMatrixElement:solution :nb :nb :1.0];
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
                                [self localBoundaryIntegral:model inSolution:solution atBoundary:bc forElement:&edges[faces[parent->FaceIndexes[j]].EdgeIndexes[j]] withNumberOfNodes:nb andParent:parent withNumberOfNodes:nb boundaryName:str1 functionIntegral:_kernWork[0]];
                                
                                n = [self getElementDofs:solution forElement:&edges[faces[parent->FaceIndexes[j]].EdgeIndexes[j]] atIndexes:_g_Ind];
                                for (k=solContainers->defDofs[parent->BodyID-1][0]*edges[faces[parent->FaceIndexes[j]].EdgeIndexes[j]].NDOFs; k<n; k++) {
                                    nb = varContainers->Perm[_g_Ind[k]];
                                    if (nb < 0) continue;
                                    nb = u_offset + solution.variable.dofs*nb + dof;
                                    if (solution.matrix.isSymmetric == YES) {
                                        [crsMatrix setSymmetricDirichletInGlobal:solution :nb :_kernWork[0]];
                                    } else {
                                        [self zeroTheNumberOfRows:nb inSolutionMatrix:solution];
                                        matContainers->RHS[nb] = _kernWork[0];
                                        [self setMatrixElement:solution :nb :nb :1.0];
                                    }
                                }
                            }
                    }
                } // end associated edges
            } else if ([listUtil listCheckPresentVariable:str2 inArray:bc] == YES) {
                n = [self getElementDofs:solution forElement:element atIndexes:_g_Ind];
                for (k=0; k<n; k++) {
                    nb = varContainers->Perm[_g_Ind[k]];
                    if (nb < 0) continue;
                    nb = u_offset + solution.variable.dofs+nb + dof;
                    
                    [self zeroTheNumberOfRows:nb inSolutionMatrix:solution];
                    [self setMatrixElement:solution :nb :nb :1.0];
                    matContainers->RHS[nb] = _kernWork[0];
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
                    if (solution.mesh.isAssociatedEdges == NO) continue;
                    
                    // If boundary edge has not dofs, move on to next edge
                    if (element->BDOFs <= 0) continue;
                    
                    // Number of nodes for this element
                    n = element->Type.NumberOfNodes;
                    
                    // Get indexes for boundary and values dofs associated to them
                    [self getBoundaryIndexes:solution.mesh forBoundaryElement:element withParentElement:parent resultVector:_g_Ind resultSize:numEdgeDofs];
                    [self localBoundaryBDOFs:model inSolution:solution atBoundary:bc forElement:element withNumberOfNodes:numEdgeDofs boundaryName:componentName resultMatrix:_kernStiff resultVector:_kernWork];
                    
                    if (solution.matrix.isSymmetric == YES) {
                        
                        for (l=0; l<n; l++) {
                            nb = varContainers->Perm[_g_Ind[l]];
                            if (nb < 0) continue;
                            nb = u_offset + solution.variable.dofs*nb + dof;
                            for (k=n; k<numEdgeDofs; k++) {
                                _kernWork[k] = _kernWork[k] - _kernStiff[k][l]*matContainers->RHS[nb];
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
                            [crsMatrix setSymmetricDirichletInGlobal:solution :nb :_kernWork[k-n]];
                        }
                    } else {
                        // Contribute this boundary to global system
                        // (i.e., solve global boundary problem)
                        for (k=n; k<numEdgeDofs; k++) {
                            nb = varContainers->Perm[_g_Ind[k]];
                            if (nb < 0) continue;
                            nb = u_offset + solution.variable.dofs*nb + dof;
                            matContainers->RHS[nb] = matContainers->RHS[nb] + _kernWork[k];
                            for (l=0; l<numEdgeDofs; l++) {
                                mb = varContainers->Perm[_g_Ind[l]];
                                if (mb < 0) continue;
                                mb = u_offset + solution.variable.dofs*mb + dof;
                                for (kk=matContainers->Rows[nb]+dof; kk<=matContainers->Rows[nb+1]-1; k+=solution.variable.dofs) {
                                    if (matContainers->Cols[kk] == mb) {
                                        matContainers->Values[kk] = matContainers->Values[kk] + _kernStiff[k][l];
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
                    if (solution.mesh.isAssociatedFaces == NO) continue;
                    
                    // Parameters of element
                    n = element->Type.NumberOfNodes;
                    
                    // Get global boundary indexes and solve dofs asscociated to them
                    [self getBoundaryIndexes:solution.mesh forBoundaryElement:element withParentElement:parent resultVector:_g_Ind resultSize:numEdgeDofs];
                    
                    // If boundary face has no dofs, skip to next boundary element
                    if (numEdgeDofs == n) continue;
                    
                    // Get local solution
                    [self localBoundaryBDOFs:model inSolution:solution atBoundary:bc forElement:element withNumberOfNodes:numEdgeDofs boundaryName:componentName resultMatrix:_kernStiff resultVector:_kernWork];
                    
                    n_start = 0;
                    if (solution.matrix.isSymmetric == YES) {
                        for (l=0; l<n; l++) {
                            nb = varContainers->Perm[_g_Ind[l]];
                            if (nb < 0) continue;
                            nb = u_offset + solution.variable.dofs*nb + dof;
                            for (k=n; k<numEdgeDofs; k++) {
                                _kernWork[k] = _kernWork[k] - _kernStiff[k][l] * matContainers->RHS[nb];
                            }
                        }
                        n_start++;
                    }
                    
                    // Contribute this entry to global boundary problem
                    for (k=n; k<numEdgeDofs; k++) {
                        nb = varContainers->Perm[_g_Ind[k]];
                        if (nb < 0) continue;
                        nb = u_offset + solution.variable.dofs*nb + dof;
                        matContainers->RHS[nb] = matContainers->RHS[nb] + _kernWork[k];
                        for (l=n_start; l<numEdgeDofs; l++) {
                            mb = varContainers->Perm[_g_Ind[l]];
                            if (mb < 0) continue;
                            mb = u_offset + solution.variable.dofs*mb + dof;
                            for (kk=matContainers->Rows[nb]+dof; kk<=matContainers->Rows[nb+1]-1; k+=solution.variable.dofs) {
                                if (matContainers->Cols[kk] == mb) {
                                    matContainers->Values[kk] = matContainers->Values[kk] + _kernStiff[k][l];
                                    break;
                                }
                            }
                        }
                    }
            }
        } //end loop over boundary elements
    } //end loop over dofs
    
    matContainers = NULL;
    varContainers = NULL;
    solContainers = NULL;
}

#pragma mark Solve

-(double)findSolution:(FEMSolution *)solution model:(FEMModel *)aModel {
    
    double norm;
    
    if ((solution.solutionInfo)[@"Linear system solver disabled"] != nil) {
        if ([(solution.solutionInfo)[@"Linear system solver disabled"] boolValue] == YES) return 0.0;
    }
    
    //TODO: add support for dump system matrix and back rotated N-T solution
    
    [self FEMKernel_solveSystem:solution model:aModel];
    norm = solution.variable.norm;
    
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
    matrixArraysContainer *matContainers;
    
    n = ipar[2];
    res = doublevec(0, ipar[2]-1);
    
    preconditioning = [[FEMPrecondition alloc] init];
    [preconditioning CRS_MatrixVectorMultiply:solution :x :res];
    
     matContainers = solution.matrix.getContainers;
    
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
        sum2 = sum2 + pow(matContainers->Values[i], 2.0);
    }
    for (i=0; i<n; i++) {
        sum3 = sum3 + pow(x[i], 2.0);
    }
    for (i=0; i<n; i++) {
        sum4 = sum4 + pow(b[i], 2.0);
    }
    err = sqrt(sum1) / ( sqrt(sum2) * sqrt(sum3) + sqrt(sum4) );
    
    free_dvector(res, 0, ipar[2]-1);
    
    matContainers = NULL;
    
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
