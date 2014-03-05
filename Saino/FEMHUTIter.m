//
//  HUTIter.m
//  Saino
//
//  Created by Hakime Seddik on 16/02/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Accelerate/Accelerate.h>

#import "FEMHUTIter.h"
#import "FEMCore.h"
#import "Utils.h"

// Used in TFQRM solve. This is the magic ratio for upperb and tolerance used in upper bound convergence test
static double UPPERB_TOL_RATIO  =  10.0;

@interface FEMHUTIter ()

-(void)HUTI_dRandVector:(double *)u ipar:(int *)ipar;
-(void)HUTI_dLuSolveAt:(int)n luMatrix:(double **)lumat afterSolve:(double *)u rightHandSide:(double *)v;

@end

@implementation FEMHUTIter

#pragma mark private methods

/**********************************************************
    This method fills a vector with pseudo random numbers
**********************************************************/
-(void)HUTI_dRandVector:(double *)u ipar:(int *)ipar {
    
    for (int i=0; i<HUTI_NDIM; i++) {
        u[i] = rand() % 100;
    }
}

/*******************************************************************
    This method constructs LU decomposition of the given matrix
    and solve LUu = v
******************************************************************/
-(void)HUTI_dLuSolveAt:(int)n luMatrix:(double **)lumat afterSolve:(double *)u rightHandSide:(double *)v {
    
    int i, j, k;
    
    // This is from Saad's book, algorithm 10.4
    
    for (i=1; i<=n; i++) {
        for (k=0; k<=i-1; k++) {
            
            // Check for small pivot
            if ( fabs(lumat[k][k]) < 1.0e-16 ) {
                NSLog(@"FEMHUTIter:HUTI_dLuSolveAt: HUTI GMRES: small pivot %lf\n", lumat[k][k]);
            }
            
            // Compute a_ik = a_ik / a_kk
            
            lumat[i][k] = lumat[i][k] / lumat[k][k];
            
            for (j=k+1; j<=n; j++) {
                // Compute a_ij = a_ij - a_ik * a_kj
                lumat[i][j] = lumat[i][j] - lumat[i][k] * lumat[k][j];
            }
        }
    }
    
    // Forward solve, Lu = v
    for (i=0; i<=n; i++) {
        
        // Compute u(i) = v(i) - sum L(i,j) u(j)
        u[i] = v[i];
        for (k=0; k<=i-1; k++) {
            u[i] = u[i] - lumat[i][k] * u[k];
        }
    }
    
    // Backward solve, u = inv(U) u
    
    for (i=n; i>=0; i--) {
        
        // Compute u(i) = u(i) - sum U(i,j) u(j)
        for (k=i+1; k<=n; k++) {
            u[i] = u[i] - lumat[i][k] * u[k];
        }
        
        // Compute u(i) = u(i) / U(i,i)
        u[i] = u[i] / lumat[i][i];
    }
}

#pragma mark public methods

- (id)init
{
    self = [super init];
    if (self) {
        // Initialization code here.
        [self hutiInit];
    }
    
    return self;
}

/*
 *
 * HUTI_Init - Initialize HUTI environment
 *
 */
-(void)hutiInit {
    
    char *evname;
    huti_init_done = NO;
    
    if (huti_init_done == YES)
        return;
    
    /* Read environment variables */
    
    if ((evname = getenv(NUMBER_OF_PROCESSORS)) == NULL)
        huti_num_of_procs = 1;
    else
        if ((huti_num_of_procs = atoi(evname)) == 0) {
            
            NSLog(@"FEMHUTIter:hutiInit: environment variable NUMBER_OF_PROCESSORS has an illegal value: %s\n", evname);
        }
    
    huti_init_done = YES;
}

-(void)hutiExit {
    
}

#pragma mark BI-CGSTAB

/**********************************************************************************************************************************
 
    This method is based on Barret et al. book: "Templates for the Solution of Linear Systems: Building blocks for
    Iterative Methods", 1993.
 
    All matrix-vector operations are done externally, so we don't need to know about the matrix structure (sparse or dense).
 
    Double precision version.
 
 **********************************************************************************************************************************/
-(void)dbicgstabSolveInSolution:(FEMSolution *)solution matrix:(FEMMatrix *)matrix ndim:(int)ndim result:(double *)x rhs:(double *)b ipar:(int *)ipar dpar:(double *)dpar pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod {
    
    int i, iter_count;
    double rho, oldrho, alpha, beta, omega;
    double *rtld, *p, *t1v, *v, *s, *t2v, *t, *r;
    double residual, rhsnorm, precrhsnorm;
    
    NSMethodSignature *pCondlSignature, *pCondrSignature, *matvecSignature, *mstopSignature;
    NSInvocation *pcondlInvocation, *pcondrInvocation, *matvecInvocation, *mstopInvocation;
    
    FEMPrecondition *preconditioning;
    FEMCore *core;
    
    // Acquire signature invocations and set selector for invocations
    pCondlSignature = [FEMPrecondition instanceMethodSignatureForSelector:pcondlMethod];
    pcondlInvocation = [NSInvocation invocationWithMethodSignature:pCondlSignature];
    [pcondlInvocation setSelector:pcondlMethod];
    
    pCondrSignature = [FEMPrecondition instanceMethodSignatureForSelector:pcondrMethod];
    pcondrInvocation = [NSInvocation invocationWithMethodSignature:pCondrSignature];
    [pcondrInvocation setSelector:pcondrMethod];
    
    matvecSignature = [FEMPrecondition instanceMethodSignatureForSelector:matvecMethod];
    matvecInvocation = [NSInvocation invocationWithMethodSignature:matvecSignature];
    [matvecInvocation setSelector:matvecMethod];
    
    mstopSignature = [FEMCore instanceMethodSignatureForSelector:mstopMethod];
    mstopInvocation = [NSInvocation invocationWithMethodSignature:mstopSignature];
    [mstopInvocation setSelector:mstopMethod];
    
    // Instanciate class for the preconditioning
    preconditioning = [[FEMPrecondition alloc] init];
    
    //Instanciate class for the kernel
    core = [FEMCore sharedCore];
    
    // Targets to invocations
    [pcondlInvocation setTarget:preconditioning];
    [pcondrInvocation setTarget:preconditioning];
    [matvecInvocation setTarget:preconditioning];
    [mstopInvocation setTarget:core];
    
    // Vector allocations
    rtld = doublevec(0, ndim-1);
    p = doublevec(0, ndim-1);
    t1v = doublevec(0, ndim-1);
    v = doublevec(0, ndim-1);
    s = doublevec(0, ndim-1);
    t2v = doublevec(0, ndim-1);
    t = doublevec(0, ndim-1);
    r = doublevec(0, ndim-1);
    
    // The actual BI-CGSTAB begins here (look at the pseudo-code in the book, page 27)
    
    // First the initialization part
    iter_count = 1;
    
    // The following applies for all matrix operations in this solver
    HUTI_EXTOP_MATTYPE = HUTI_MAT_NOTTRPSED;
    
    if (HUTI_STOPC == HUTI_TRESID_SCALED_BYB || HUTI_STOPC == HUTI_PRESID_SCALED_BYB) {
        rhsnorm = cblas_dnrm2(HUTI_NDIM, b, 1);
    }
    if (HUTI_STOPC == HUTI_PRESID_SCALED_BYPRECB) {
        [pcondlInvocation setArgument:&matrix atIndex:2];
        [pcondlInvocation setArgument:&p atIndex:3];
        [pcondlInvocation setArgument:&b atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];
        
        precrhsnorm = cblas_dnrm2(HUTI_NDIM, p, 1);
    }
    
    // Generate vector xvec if needed
    if (HUTI_INITIALX == HUTI_RANDOMX) {
        [self HUTI_dRandVector:x ipar:ipar];
    } else if (HUTI_INITIALX != HUTI_USERSUPPLIEDX) {
        for (i=0; i<ndim; i++) {
            x[i] = 1;
        }
    }
    
    [matvecInvocation setArgument:&matrix atIndex:2];
    [matvecInvocation setArgument:&x atIndex:3];
    [matvecInvocation setArgument:&r atIndex:4];
    [matvecInvocation setArgument:&ipar atIndex:5];
    [matvecInvocation invoke];
    
    for (i=0; i<ndim; i++) {
        r[i] = b[i] - r[i];
        rtld[i] = r[i];
    }
    
    memset( p, 0.0, ndim*sizeof(double) );
    memset( v, 0.0, ndim*sizeof(double) );
    
    oldrho = 1;
    omega = 1;
    alpha = 0;
    
    // This is where the loop starts
    while (1) {
        
        rho = cblas_ddot(HUTI_NDIM, rtld, 1, r, 1);
        if (rho == 0.0) {
            HUTI_INFO = HUTI_BICGSTAB_RHO;
            break;
        }
        
        beta = ( rho * alpha) / ( oldrho * omega );
        for (i=0; i<ndim; i++) {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }
        
        [pcondlInvocation setArgument:&matrix atIndex:2];
        [pcondlInvocation setArgument:&v atIndex:3];
        [pcondlInvocation setArgument:&p atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];
        
        [pcondrInvocation setArgument:&matrix atIndex:2];
        [pcondrInvocation setArgument:&t1v atIndex:3];
        [pcondrInvocation setArgument:&v atIndex:4];
        [pcondrInvocation setArgument:&ipar atIndex:5];
        [pcondrInvocation invoke];
        
        [matvecInvocation setArgument:&matrix atIndex:2];
        [matvecInvocation setArgument:&t1v atIndex:3];
        [matvecInvocation setArgument:&v atIndex:4];
        [matvecInvocation setArgument:&ipar atIndex:5];
        [matvecInvocation invoke];
        
        alpha = rho / cblas_ddot(HUTI_NDIM, rtld, 1, v, 1);
        for (i=0; i<ndim; i++) {
            s[i] = r[i] - alpha * v[i];
        }
        
        residual = cblas_dnrm2(HUTI_NDIM, s, 1);
        if (residual < HUTI_EPSILON) {
            for (i=0; i<ndim; i++) {
                x[i] = x[i] + alpha * t1v[i];
                HUTI_INFO = HUTI_BICGSTAB_SNORM;
                break;
            }
        }
        
        [pcondlInvocation setArgument:&matrix atIndex:2];
        [pcondlInvocation setArgument:&t atIndex:3];
        [pcondlInvocation setArgument:&s atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];
        
        [pcondrInvocation setArgument:&matrix atIndex:2];
        [pcondrInvocation setArgument:&t2v atIndex:3];
        [pcondrInvocation setArgument:&t atIndex:4];
        [pcondrInvocation setArgument:&ipar atIndex:5];
        [pcondrInvocation invoke];
        
        [matvecInvocation setArgument:&matrix atIndex:2];
        [matvecInvocation setArgument:&t2v atIndex:3];
        [matvecInvocation setArgument:&t atIndex:4];
        [matvecInvocation setArgument:&ipar atIndex:5];
        [matvecInvocation invoke];
        
        omega = ( cblas_ddot(HUTI_NDIM, t, 1, s, 1) ) / ( cblas_ddot(HUTI_NDIM, t, 1, t, 1) );
        for (i=0; i<ndim; i++) {
            x[i] = x[i] + alpha * t1v[i] + omega * t2v[i];
            r[i] = s[i] - omega * t[i];
        }
        
        // Check the convergence against selected stopping criterion
        switch (HUTI_STOPC) {
            case HUTI_TRUERESIDUAL:
                [matvecInvocation setArgument:&matrix atIndex:2];
                [matvecInvocation setArgument:&x atIndex:3];
                [matvecInvocation setArgument:&t2v atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    t1v[i] = t2v[i] - b[i];
                }
                residual = cblas_dnrm2(HUTI_NDIM, t1v, 1);
                break;
            case HUTI_TRESID_SCALED_BYB:
                [matvecInvocation setArgument:&matrix atIndex:2];
                [matvecInvocation setArgument:&x atIndex:3];
                [matvecInvocation setArgument:&t2v atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    t1v[i] = t2v[i] - b[i];
                }
                residual = cblas_dnrm2(HUTI_NDIM, t1v, 1) / rhsnorm;
                break;
            case HUTI_PSEUDORESIDUAL:
                residual = cblas_dnrm2(HUTI_NDIM, r, 1);
                break;
            case HUTI_PRESID_SCALED_BYB:
                residual = cblas_dnrm2(HUTI_NDIM, r, 1) / rhsnorm;
                break;
            case HUTI_PRESID_SCALED_BYPRECB:
                residual = cblas_dnrm2(HUTI_NDIM, r, 1) / precrhsnorm;
                break;
            case HUTI_XDIFF_NORM:
                for (i=0; i<ndim; i++) {
                    t1v[i] = alpha * t1v[i] + omega * t2v[i];
                }
                residual = cblas_dnrm2(HUTI_NDIM, t1v, 1);
                break;
            case HUTI_USUPPLIED_STOPC:
                [mstopInvocation setArgument:&matrix atIndex:2];
                [mstopInvocation setArgument:&x atIndex:3];
                [mstopInvocation setArgument:&b atIndex:4];
                [mstopInvocation setArgument:&r atIndex:5];
                [mstopInvocation setArgument:&ipar atIndex:6];
                [mstopInvocation invoke];
                [mstopInvocation getReturnValue:&residual];
                break;
            default:
                [matvecInvocation setArgument:&matrix atIndex:2];
                [matvecInvocation setArgument:&x atIndex:3];
                [matvecInvocation setArgument:&t2v atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    t1v[i] = t2v[i] - b[i];
                }
                residual = cblas_dnrm2(HUTI_NDIM, t1v, 1);
                break;
        }
        
        // Print debugging info if required
        if (HUTI_DBUGLVL != HUTI_NO_DEBUG) {
            if ( (iter_count % HUTI_DBUGLVL) == 0 ) {
                NSLog(@"FEMHUTIter:dbicgstabSolveInSolution: %d %11.4e\n", iter_count, residual);
            }
        }
        
        if (residual < HUTI_TOLERANCE) {
            HUTI_INFO = HUTI_CONVERGENCE;
            break;
        }
        
        if (omega == 0) {
            HUTI_INFO = HUTI_BICGSTAB_OMEGA;
            break;
        }
        
        oldrho = rho;
        
        // Return back to the top of the iteration loop or exit
        iter_count = iter_count + 1;
        if (iter_count > HUTI_MAXIT) {
            HUTI_INFO = HUTI_MAXITER;
            break;
        }
        
    }
    
    // We have exited the loop after enough iterations or broke down
    if (HUTI_DBUGLVL != HUTI_NO_DEBUG) {
        NSLog(@"FEMHUTIter:dbicgstabSolveInSolution: %d %11.4e\n", iter_count, residual);
    }
    
    HUTI_ITERS = iter_count;
    
    // Release memory
    free_dvector(rtld, 0, ndim-1);
    free_dvector(p, 0, ndim-1);
    free_dvector(t1v, 0, ndim-1);
    free_dvector(v, 0, ndim-1);
    free_dvector(s, 0, ndim-1);
    free_dvector(t2v, 0, ndim-1);
    free_dvector(t, 0, ndim-1);
    free_dvector(r, 0, ndim-1);
}

#pragma mark BI-CGSTAB(2)

/**********************************************************************************************************************************
 
    This method is based on the Henk A. Van der Vorst paper: "Parallel Iretative Solution Methods for the Linear Systems
    arising from Discretized PDE's". This is the Bi-CGSTAB(2) version.
 
    All matrix-vector operations are done externally, so we don't need to know about the matrix structure (sparse or dense).
 
    Double precision version.
 
**********************************************************************************************************************************/
-(void)dbicgstab2SolveInSolution:(FEMSolution *)solution matrix:(FEMMatrix *)matrix ndim:(int)ndim result:(double *)x rhs:(double *)b ipar:(int *)ipar dpar:(double *)dpar pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod {
    
    int i, iter_count;
    double rho, oldrho, alpha, beta, omega1, omega2;
    double tau, delta, myy;
    double *rtld, *u, *t1v, *v, *s, *w, *t, *r;
    double residual, rhsnorm, precrhsnorm;
    
    NSMethodSignature *pCondlSignature, *pCondrSignature, *matvecSignature, *mstopSignature;
    NSInvocation *pcondlInvocation, *pcondrInvocation, *matvecInvocation, *mstopInvocation;
    
    FEMPrecondition *preconditioning;
    FEMCore *core;
    
    // Acquire signature invocations and set selector for invocations
    pCondlSignature = [FEMPrecondition instanceMethodSignatureForSelector:pcondlMethod];
    pcondlInvocation = [NSInvocation invocationWithMethodSignature:pCondlSignature];
    [pcondlInvocation setSelector:pcondlMethod];
    
    pCondrSignature = [FEMPrecondition instanceMethodSignatureForSelector:pcondrMethod];
    pcondrInvocation = [NSInvocation invocationWithMethodSignature:pCondrSignature];
    [pcondrInvocation setSelector:pcondrMethod];
    
    matvecSignature = [FEMPrecondition instanceMethodSignatureForSelector:matvecMethod];
    matvecInvocation = [NSInvocation invocationWithMethodSignature:matvecSignature];
    [matvecInvocation setSelector:matvecMethod];
    
    mstopSignature = [FEMCore instanceMethodSignatureForSelector:mstopMethod];
    mstopInvocation = [NSInvocation invocationWithMethodSignature:mstopSignature];
    [mstopInvocation setSelector:mstopMethod];
    
    // Instanciate class for the preconditioning
    preconditioning = [[FEMPrecondition alloc] init];
    
    //Instanciate class for the kernel
    core = [FEMCore sharedCore];
    
    // Targets to invocations
    [pcondlInvocation setTarget:preconditioning];
    [pcondrInvocation setTarget:preconditioning];
    [matvecInvocation setTarget:preconditioning];
    [mstopInvocation setTarget:core];
    
    // Vector allocations
    rtld = doublevec(0, ndim-1);
    u = doublevec(0, ndim-1);
    t1v = doublevec(0, ndim-1);
    v = doublevec(0, ndim-1);
    s = doublevec(0, ndim-1);
    w = doublevec(0, ndim-1);
    t = doublevec(0, ndim-1);
    r = doublevec(0, ndim-1);
    
    // The actual BI-CGSTAB begins here (look at the pseudo-code in the book, page 27)
    
    // First the initialization part
    iter_count = 1;

    // The following applied for all matrix operations in this solver
    HUTI_EXTOP_MATTYPE = HUTI_MAT_NOTTRPSED;
    
    // Norms of right-hand side vector are used in convergence tests
    if (HUTI_STOPC == HUTI_TRESID_SCALED_BYB || HUTI_STOPC == HUTI_PRESID_SCALED_BYB) {
        rhsnorm = cblas_dnrm2(HUTI_NDIM, b, 1);
    }
    if (HUTI_STOPC == HUTI_PRESID_SCALED_BYPRECB) {
        [pcondlInvocation setArgument:&matrix atIndex:2];
        [pcondlInvocation setArgument:&t1v atIndex:3];
        [pcondlInvocation setArgument:&b atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];
        
        precrhsnorm = cblas_dnrm2(HUTI_NDIM, t1v, 1);
    }
    
    // Generate vector xvec if needed
    if (HUTI_INITIALX == HUTI_RANDOMX) {
        [self HUTI_dRandVector:x ipar:ipar];
    } else if (HUTI_INITIALX != HUTI_USERSUPPLIEDX) {
        for (i=0; i<ndim; i++) {
            x[i] = 1;
        }
    }
    
    [pcondrInvocation setArgument:&matrix atIndex:2];
    [pcondrInvocation setArgument:&u atIndex:3];
    [pcondrInvocation setArgument:&x atIndex:4];
    [pcondrInvocation setArgument:&ipar atIndex:5];
    [pcondrInvocation invoke];
    
    [matvecInvocation setArgument:&matrix atIndex:2];
    [matvecInvocation setArgument:&u atIndex:3];
    [matvecInvocation setArgument:&r atIndex:4];
    [matvecInvocation setArgument:&ipar atIndex:5];
    [matvecInvocation invoke];
    
    for (i=0; i<ndim; i++) {
        u[i] = b[i] - r[i];
    }
    
    [pcondlInvocation setArgument:&matrix atIndex:2];
    [pcondlInvocation setArgument:&r atIndex:3];
    [pcondlInvocation setArgument:&u atIndex:4];
    [pcondlInvocation setArgument:&ipar atIndex:5];
    [pcondlInvocation invoke];
    
    for (i=0; i<ndim; i++) {
        rtld[i] = r[i];
    }
    memset( u, 0.0, ndim*sizeof(double) );
    oldrho = 1;
    omega2 = 1;
    alpha = 0;
    
    // This is where the loop starts
    while (1) {
        
        oldrho = -omega2 * oldrho;
        
        // This is the even BI-CG step
        rho = cblas_ddot(HUTI_NDIM, rtld, 1, r, 1);
        if (rho == 0.0) {
            HUTI_INFO = HUTI_BICGSTAB_2_RHO;
            break;
        }
        
        beta = ( rho * alpha ) / oldrho;
        oldrho = rho;
        for (i=0; i<ndim; i++) {
            u[i] = r[i] - beta * u[i];
        }
        
        [pcondrInvocation setArgument:&matrix atIndex:2];
        [pcondrInvocation setArgument:&v atIndex:3];
        [pcondrInvocation setArgument:&u atIndex:4];
        [pcondrInvocation setArgument:&ipar atIndex:5];
        [pcondrInvocation invoke];
        
        [matvecInvocation setArgument:&matrix atIndex:2];
        [matvecInvocation setArgument:&v atIndex:3];
        [matvecInvocation setArgument:&t1v atIndex:4];
        [matvecInvocation setArgument:&ipar atIndex:5];
        [matvecInvocation invoke];
        
        [pcondlInvocation setArgument:&matrix atIndex:2];
        [pcondlInvocation setArgument:&v atIndex:3];
        [pcondlInvocation setArgument:&t1v atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];

        alpha = oldrho / cblas_ddot(HUTI_NDIM, rtld, 1, v, 1);
        for (i=0; i<ndim; i++) {
            r[i] = r[i] - alpha * v[i];
        }
        
        [pcondrInvocation setArgument:&matrix atIndex:2];
        [pcondrInvocation setArgument:&s atIndex:3];
        [pcondrInvocation setArgument:&r atIndex:4];
        [pcondrInvocation setArgument:&ipar atIndex:5];
        [pcondrInvocation invoke];
        
        [matvecInvocation setArgument:&matrix atIndex:2];
        [matvecInvocation setArgument:&s atIndex:3];
        [matvecInvocation setArgument:&t1v atIndex:4];
        [matvecInvocation setArgument:&ipar atIndex:5];
        [matvecInvocation invoke];
        
        [pcondlInvocation setArgument:&matrix atIndex:2];
        [pcondlInvocation setArgument:&s atIndex:3];
        [pcondlInvocation setArgument:&t1v atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];
        
        for (i=0; i<ndim; i++) {
            x[i] = x[i] + alpha * u[i];
        }
        
        // This is the odd BI-CG step
        rho = cblas_ddot(HUTI_NDIM, rtld, 1, s, 1);
        if (rho == 0) {
            HUTI_INFO = HUTI_BICGSTAB_2_RHO;
            break;
        }

        beta = ( rho * alpha ) / oldrho;
        oldrho = rho;
        for (i=0; i<ndim; i++) {
            v[i] = s[i] - beta * v[i];
        }
        
        [pcondrInvocation setArgument:&matrix atIndex:2];
        [pcondrInvocation setArgument:&w atIndex:3];
        [pcondrInvocation setArgument:&v atIndex:4];
        [pcondrInvocation setArgument:&ipar atIndex:5];
        [pcondrInvocation invoke];
        
        [matvecInvocation setArgument:&matrix atIndex:2];
        [matvecInvocation setArgument:&w atIndex:3];
        [matvecInvocation setArgument:&t1v atIndex:4];
        [matvecInvocation setArgument:&ipar atIndex:5];
        [matvecInvocation invoke];
        
        [pcondlInvocation setArgument:&matrix atIndex:2];
        [pcondlInvocation setArgument:&w atIndex:3];
        [pcondlInvocation setArgument:&t1v atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];
        
        alpha = oldrho / cblas_ddot(HUTI_NDIM, rtld, 1, w, 1);
        for (i=0; i<ndim; i++) {
            u[i] = r[i] - beta * u[i];
            r[i] = r[i] - alpha * v[i];
            s[i] = s[i] - alpha * w[i];
        }
        
        [pcondrInvocation setArgument:&matrix atIndex:2];
        [pcondrInvocation setArgument:&t atIndex:3];
        [pcondrInvocation setArgument:&s atIndex:4];
        [pcondrInvocation setArgument:&ipar atIndex:5];
        [pcondrInvocation invoke];
        
        [matvecInvocation setArgument:&matrix atIndex:2];
        [matvecInvocation setArgument:&t atIndex:3];
        [matvecInvocation setArgument:&t1v atIndex:4];
        [matvecInvocation setArgument:&ipar atIndex:5];
        [matvecInvocation invoke];
        
        [pcondlInvocation setArgument:&matrix atIndex:2];
        [pcondlInvocation setArgument:&t atIndex:3];
        [pcondlInvocation setArgument:&t1v atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];
        
        // This is the GCR(2) part
        omega1 = cblas_ddot(HUTI_NDIM, r, 1, s, 1);
        myy = cblas_ddot(HUTI_NDIM, s, 1, s, 1);
        delta = cblas_ddot(HUTI_NDIM, s, 1, t, 1);
        tau = cblas_ddot(HUTI_NDIM, t, 1, t, 1);
        omega2 = cblas_ddot(HUTI_NDIM, r, 1, t, 1);
        
        tau = tau - ( delta * delta ) / myy;
        omega2 = ( omega2 - ( delta * omega1 ) / myy ) / tau;
        omega1 = ( omega1 - delta * omega2 ) / myy;
        
        for (i=0; i<ndim; i++) {
            x[i] = x[i] + omega1 * r[i] + omega2 * s[i] + alpha * u[i];
            r[i] = r[i] - omega1 * s[i] - omega2 * t[i];
        }
        
        // Check the convergence against selected stopping criterion
        switch (HUTI_STOPC) {
            case HUTI_TRUERESIDUAL:
                [pcondrInvocation setArgument:&matrix atIndex:2];
                [pcondrInvocation setArgument:&s atIndex:3];
                [pcondrInvocation setArgument:&x atIndex:4];
                [pcondrInvocation setArgument:&ipar atIndex:5];
                [pcondrInvocation invoke];
                
                [matvecInvocation setArgument:&matrix atIndex:2];
                [matvecInvocation setArgument:&s atIndex:3];
                [matvecInvocation setArgument:&t1v atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    t1v[i] = t1v[i] - b[i];
                }
                [pcondlInvocation setArgument:&matrix atIndex:2];
                [pcondlInvocation setArgument:&s atIndex:3];
                [pcondlInvocation setArgument:&t1v atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, s, 1);
                break;
            case HUTI_TRESID_SCALED_BYB:
                [pcondrInvocation setArgument:&matrix atIndex:2];
                [pcondrInvocation setArgument:&s atIndex:3];
                [pcondrInvocation setArgument:&x atIndex:4];
                [pcondrInvocation setArgument:&ipar atIndex:5];
                [pcondrInvocation invoke];
                
                [matvecInvocation setArgument:&matrix atIndex:2];
                [matvecInvocation setArgument:&s atIndex:3];
                [matvecInvocation setArgument:&t1v atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    t1v[i] = t1v[i] - b[i];
                }
                [pcondlInvocation setArgument:&matrix atIndex:2];
                [pcondlInvocation setArgument:&s atIndex:3];
                [pcondlInvocation setArgument:&t1v atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, s, 1) / rhsnorm;
                break;
            case HUTI_PSEUDORESIDUAL:
                residual = cblas_dnrm2(HUTI_NDIM, r, 1);
                break;
            case HUTI_PRESID_SCALED_BYB:
                residual = cblas_dnrm2(HUTI_NDIM, r, 1) / rhsnorm;
                break;
            case HUTI_PRESID_SCALED_BYPRECB:
                residual = cblas_dnrm2(HUTI_NDIM, r, 1) / precrhsnorm;
                break;
            case HUTI_XDIFF_NORM:
                for (i=0; i<ndim; i++) {
                    t1v[i] = omega1 * r[i] + omega2 * s[i] + alpha * u[i];
                }
                residual = cblas_dnrm2(HUTI_NDIM, t1v, 1);
                break;
            case HUTI_USUPPLIED_STOPC:
                [mstopInvocation setArgument:&matrix atIndex:2];
                [mstopInvocation setArgument:&x atIndex:3];
                [mstopInvocation setArgument:&b atIndex:4];
                [mstopInvocation setArgument:&r atIndex:5];
                [mstopInvocation setArgument:&ipar atIndex:6];
                [mstopInvocation invoke];
                [mstopInvocation getReturnValue:&residual];
                break;
            default:
                [pcondrInvocation setArgument:&matrix atIndex:2];
                [pcondrInvocation setArgument:&s atIndex:3];
                [pcondrInvocation setArgument:&x atIndex:4];
                [pcondrInvocation setArgument:&ipar atIndex:5];
                [pcondrInvocation invoke];

                [matvecInvocation setArgument:&matrix atIndex:2];
                [matvecInvocation setArgument:&s atIndex:3];
                [matvecInvocation setArgument:&t1v atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    t1v[i] = t1v[i] - b[i];
                }
                [pcondlInvocation setArgument:&matrix atIndex:2];
                [pcondlInvocation setArgument:&s atIndex:3];
                [pcondlInvocation setArgument:&t1v atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, s, 1);
                break;
        }
        
        // Print debugging info if required
        if (HUTI_DBUGLVL != HUTI_NO_DEBUG) {
            if ( (iter_count % HUTI_DBUGLVL) == 0 ) {
                NSLog(@"FEMHUTIter:dbicgstab2SolveInSolution: %d %11.4e\n", iter_count, residual);
            }
        }
        
        if (residual < HUTI_TOLERANCE) {
            HUTI_INFO = HUTI_CONVERGENCE;
            break;
        }
        
        for (i=0; i<ndim; i++) {
            u[i] = u[i] - omega1 * v[i] - omega2 * w[i];
        }
        
        // Return back to the top of the iteration loop (without initialization)
        iter_count = iter_count + 1;
        if (iter_count > HUTI_MAXIT) {
            HUTI_INFO = HUTI_MAXITER;
            break;
        }
        
    }
    
    // We have exited the loop after enough iterations or broke down
    if (HUTI_DBUGLVL != HUTI_NO_DEBUG) {
        NSLog(@"FEMHUTIter:dbicgstab2SolveInSolution: %d %11.4e\n", iter_count, residual);
    }
    
    HUTI_ITERS = iter_count;

    // Release memory
    free_dvector(rtld, 0, ndim-1);
    free_dvector(u, 0, ndim-1);
    free_dvector(t1v, 0, ndim-1);
    free_dvector(v, 0, ndim-1);
    free_dvector(s, 0, ndim-1);
    free_dvector(w, 0, ndim-1);
    free_dvector(t, 0, ndim-1);
    free_dvector(r, 0, ndim-1);
}

#pragma mark TFQMR

/**********************************************************************************************************************************
 
    This method is based on the Roland W. Freund paper: "A Transpose-Free Quasi-Minimal Residual Algorithm for Non-Hermitian
    Linear Systems" (SIAM J. Sci. Comput., March 1993).
 
    All matrix-vector operations are done externally, so we don't need to know about the matrix structure (sparse or dense).
 
    Double precision version.
 
**********************************************************************************************************************************/
-(void)dtfqmrSolveInSolution:(FEMSolution *)solution matrix:(FEMMatrix *)matrix ndim:(int)ndim result:(double *)x rhs:(double *)b ipar:(int *)ipar dpar:(double *)dpar pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod {
    
    int i, iter_count;
    double rho, oldrho=0, alpha, beta, gamma, oldgamma, eta, tau, c;
    double *v, *y, *ynew, *rtld, *t1v, *t2v, *w, *d, *r, *trv;
    double residual, upperb, rhsnorm, precrhsnorm;
    
    NSMethodSignature *pCondlSignature, *pCondrSignature, *matvecSignature, *mstopSignature;
    NSInvocation *pcondlInvocation, *pcondrInvocation, *matvecInvocation, *mstopInvocation;
    
    FEMPrecondition *preconditioning;
    FEMCore *core;
    
    // Acquire signature invocations and set selector for invocations
    pCondlSignature = [FEMPrecondition instanceMethodSignatureForSelector:pcondlMethod];
    pcondlInvocation = [NSInvocation invocationWithMethodSignature:pCondlSignature];
    [pcondlInvocation setSelector:pcondlMethod];
    
    pCondrSignature = [FEMPrecondition instanceMethodSignatureForSelector:pcondrMethod];
    pcondrInvocation = [NSInvocation invocationWithMethodSignature:pCondrSignature];
    [pcondrInvocation setSelector:pcondrMethod];
    
    matvecSignature = [FEMPrecondition instanceMethodSignatureForSelector:matvecMethod];
    matvecInvocation = [NSInvocation invocationWithMethodSignature:matvecSignature];
    [matvecInvocation setSelector:matvecMethod];
    
    mstopSignature = [FEMCore instanceMethodSignatureForSelector:mstopMethod];
    mstopInvocation = [NSInvocation invocationWithMethodSignature:mstopSignature];
    [mstopInvocation setSelector:mstopMethod];
    
    // Instanciate class for the preconditioning
    preconditioning = [[FEMPrecondition alloc] init];
    
    //Instanciate class for the kernel
    core = [FEMCore sharedCore];
    
    // Targets to invocations
    [pcondlInvocation setTarget:preconditioning];
    [pcondrInvocation setTarget:preconditioning];
    [matvecInvocation setTarget:preconditioning];
    [mstopInvocation setTarget:core];
        
    // Vector allocations
    v = doublevec(0, ndim-1);
    y = doublevec(0, ndim-1);
    ynew = doublevec(0, ndim-1);
    rtld = doublevec(0, ndim-1);
    t1v = doublevec(0, ndim-1);
    t2v = doublevec(0, ndim-1);
    w = doublevec(0, ndim-1);
    d = doublevec(0, ndim-1);
    r = doublevec(0, ndim-1);
    trv = doublevec(0, ndim-1);
    
    // The actual TFQMR begins here (look the pseudo code in the "A Transpose-Free...." paper, algorithm 5.1)
    
    // First the initialization part
    iter_count = 1;
    
    // The following applied for all matrix operations in this solver
    HUTI_EXTOP_MATTYPE = HUTI_MAT_NOTTRPSED;
    
    // Norms of right-hand side vector are used in convergence tests
    if (HUTI_STOPC == HUTI_TRESID_SCALED_BYB || HUTI_STOPC == HUTI_PRESID_SCALED_BYB || HUTI_STOPC == HUTI_UPPERB_STOPC) {
        rhsnorm = cblas_dnrm2(HUTI_NDIM, b, 1);
    }
    if (HUTI_STOPC == HUTI_PRESID_SCALED_BYPRECB) {
        [pcondlInvocation setArgument:&matrix atIndex:2];
        [pcondlInvocation setArgument:&d atIndex:3];
        [pcondlInvocation setArgument:&b atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];
        
        precrhsnorm = cblas_dnrm2(HUTI_NDIM, d, 1);
    }
    
    // Part 1A - 1C
    
    //  Generate vector xvec if needed
    if (HUTI_INITIALX == HUTI_RANDOMX) {
        [self HUTI_dRandVector:x ipar:ipar];
    } else if (HUTI_INITIALX != HUTI_USERSUPPLIEDX) {
        for (i=0; i<ndim; i++) {
            x[i] = 1;
        }
    }
    
    [pcondrInvocation setArgument:&matrix atIndex:2];
    [pcondrInvocation setArgument:&d atIndex:3];
    [pcondrInvocation setArgument:&x atIndex:4];
    [pcondrInvocation setArgument:&ipar atIndex:5];
    [pcondrInvocation invoke];
    
    [matvecInvocation setArgument:&matrix atIndex:2];
    [matvecInvocation setArgument:&d atIndex:3];
    [matvecInvocation setArgument:&r atIndex:4];
    [matvecInvocation setArgument:&ipar atIndex:5];
    [matvecInvocation invoke];
    
    for (i=0; i<ndim; i++) {
        d[i] = b[i] - r[i];
    }
    
    [pcondlInvocation setArgument:&matrix atIndex:2];
    [pcondlInvocation setArgument:&r atIndex:3];
    [pcondlInvocation setArgument:&d atIndex:4];
    [pcondlInvocation setArgument:&ipar atIndex:5];
    [pcondlInvocation invoke];
    
    for (i=0; i<ndim; i++) {
        y[i] = r[i];
        w[i] = r[i];
    }
    [pcondrInvocation setArgument:&matrix atIndex:2];
    [pcondrInvocation setArgument:&v atIndex:3];
    [pcondrInvocation setArgument:&y atIndex:4];
    [pcondrInvocation setArgument:&ipar atIndex:5];
    [pcondrInvocation invoke];
    
    [matvecInvocation setArgument:&matrix atIndex:2];
    [matvecInvocation setArgument:&v atIndex:3];
    [matvecInvocation setArgument:&d atIndex:4];
    [matvecInvocation setArgument:&ipar atIndex:5];
    [matvecInvocation invoke];
    
    [pcondlInvocation setArgument:&matrix atIndex:2];
    [pcondlInvocation setArgument:&v atIndex:3];
    [pcondlInvocation setArgument:&d atIndex:4];
    [pcondlInvocation setArgument:&ipar atIndex:5];
    [pcondlInvocation invoke];
    for (i=0; i<ndim; i++) {
        t2v[i] = v[i];
    }
    
    memset( d, 0.0, ndim*sizeof(double) );
    tau = cblas_dnrm2(HUTI_NDIM, r, 1);
    oldgamma = 0;
    gamma = 0;
    eta = 0;
    
    for (i=0; i<ndim; i++) {
        rtld[i] = r[i];
    }
    
    oldrho = cblas_ddot(HUTI_NDIM, rtld, 1, r, 1);
    if (oldrho == 0.0) {
        HUTI_INFO = HUTI_TFQMR_RHO;
        goto jump;
    }
    
    // This is where the loop starts
    
    // Part 2A
    
    while (1) {
        
        alpha = oldrho / cblas_ddot(HUTI_NDIM, rtld, 1, v, 1);
        for (i=0; i<ndim; i++) {
            ynew[i] = y[i] - alpha * v[i];
        }
        
        // Part 2B
        
        // This is the inner loop from 2n-1 to 2n
        // First the 2n-1
        // Note: We have already MATRIX * Y in t2v
        
        for (i=0; i<ndim; i++) {
            w[i] = w[i] - alpha * t2v[i];
        }
        gamma = ( cblas_dnrm2(HUTI_NDIM, w, 1) ) / tau;
        c = 1 / sqrt( 1 + gamma * gamma );
        tau = tau * gamma * c;
        
        for (i=0; i<ndim; i++) {
            d[i] = y[i] + ( ( oldgamma * oldgamma * eta ) / alpha ) * d[i];
        }
        eta = c * c * alpha;
        for (i=0; i<ndim; i++) {
            x[i] = x[i] + eta * d[i];
        }
        
        oldgamma = gamma;
        
        // Check the convergence against selected stopping criterion
        switch (HUTI_STOPC) {
            case HUTI_TRUERESIDUAL:
                [pcondrInvocation setArgument:&matrix atIndex:2];
                [pcondrInvocation setArgument:&trv atIndex:3];
                [pcondrInvocation setArgument:&x atIndex:4];
                [pcondrInvocation setArgument:&ipar atIndex:5];
                [pcondrInvocation invoke];
                
                [matvecInvocation setArgument:&matrix atIndex:2];
                [matvecInvocation setArgument:&trv atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    trv[i] = r[i] - b[i];
                }
                [pcondlInvocation setArgument:&matrix atIndex:2];
                [pcondlInvocation setArgument:&r atIndex:3];
                [pcondlInvocation setArgument:&trv atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, r, 1);
                break;
            case HUTI_TRESID_SCALED_BYB:
                [pcondrInvocation setArgument:&matrix atIndex:2];
                [pcondrInvocation setArgument:&trv atIndex:3];
                [pcondrInvocation setArgument:&x atIndex:4];
                [pcondrInvocation setArgument:&ipar atIndex:5];
                [pcondrInvocation invoke];
                
                [matvecInvocation setArgument:&matrix atIndex:2];
                [matvecInvocation setArgument:&trv atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    trv[i] = r[i] - b[i];
                }
                [pcondlInvocation setArgument:&matrix atIndex:2];
                [pcondlInvocation setArgument:&r atIndex:3];
                [pcondlInvocation setArgument:&trv atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, r, 1) / rhsnorm;
                break;
            case HUTI_PSEUDORESIDUAL:
                [matvecInvocation setArgument:&matrix atIndex:2];
                [matvecInvocation setArgument:&x atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    r[i] = r[i] - b[i];
                }
                [pcondlInvocation setArgument:&matrix atIndex:2];
                [pcondlInvocation setArgument:&trv atIndex:3];
                [pcondlInvocation setArgument:&r atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, trv, 1);
                break;
            case HUTI_PRESID_SCALED_BYB:
                [matvecInvocation setArgument:&matrix atIndex:2];
                [matvecInvocation setArgument:&x atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    r[i] = r[i] - b[i];
                }
                [pcondlInvocation setArgument:&matrix atIndex:2];
                [pcondlInvocation setArgument:&trv atIndex:3];
                [pcondlInvocation setArgument:&r atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, trv, 1) / rhsnorm;
                break;
            case HUTI_PRESID_SCALED_BYPRECB:
                [matvecInvocation setArgument:&matrix atIndex:2];
                [matvecInvocation setArgument:&x atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    r[i] = r[i] - b[i];
                }
                [pcondlInvocation setArgument:&matrix atIndex:2];
                [pcondlInvocation setArgument:&trv atIndex:3];
                [pcondlInvocation setArgument:&r atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, trv, 1) / precrhsnorm;
                break;
            case HUTI_XDIFF_NORM:
                for (i=0; i<ndim; i++) {
                    r[i] = eta * d[i];
                }
                residual = cblas_dnrm2(HUTI_NDIM, r, 1);
                break;
            case HUTI_UPPERB_STOPC:
                upperb = sqrt( 2.0 * (double) iter_count ) * tau / rhsnorm;
                if ( (upperb / HUTI_TOLERANCE) < UPPERB_TOL_RATIO ) {
                    [pcondrInvocation setArgument:&matrix atIndex:2];
                    [pcondrInvocation setArgument:&trv atIndex:3];
                    [pcondrInvocation setArgument:&x atIndex:4];
                    [pcondrInvocation setArgument:&ipar atIndex:5];
                    [pcondrInvocation invoke];
                    
                    [matvecInvocation setArgument:&matrix atIndex:2];
                    [matvecInvocation setArgument:&trv atIndex:3];
                    [matvecInvocation setArgument:&r atIndex:4];
                    [matvecInvocation setArgument:&ipar atIndex:5];
                    [matvecInvocation invoke];
                    for (i=0; i<ndim; i++) {
                        trv[i] = r[i] - b[i];
                    }
                    [pcondlInvocation setArgument:&matrix atIndex:2];
                    [pcondlInvocation setArgument:&r atIndex:3];
                    [pcondlInvocation setArgument:&trv atIndex:4];
                    [pcondlInvocation setArgument:&ipar atIndex:5];
                    [pcondlInvocation invoke];
                    residual = cblas_dnrm2(HUTI_NDIM, r, 1) / rhsnorm;
                } else {
                    residual = upperb;
                }
                break;
            case HUTI_USUPPLIED_STOPC:
                [mstopInvocation setArgument:&matrix atIndex:2];
                [mstopInvocation setArgument:&x atIndex:3];
                [mstopInvocation setArgument:&b atIndex:4];
                [mstopInvocation setArgument:&r atIndex:5];
                [mstopInvocation setArgument:&ipar atIndex:6];
                [mstopInvocation invoke];
                [mstopInvocation getReturnValue:&residual];
                break;
            default:
                [pcondrInvocation setArgument:&matrix atIndex:2];
                [pcondrInvocation setArgument:&trv atIndex:3];
                [pcondrInvocation setArgument:&x atIndex:4];
                [pcondrInvocation setArgument:&ipar atIndex:5];
                [pcondrInvocation invoke];
                
                [matvecInvocation setArgument:&matrix atIndex:2];
                [matvecInvocation setArgument:&trv atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    trv[i] = r[i] - b[i];
                }
                [pcondlInvocation setArgument:&matrix atIndex:2];
                [pcondlInvocation setArgument:&r atIndex:3];
                [pcondlInvocation setArgument:&trv atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, r, 1);
                break;
        }
        
        if (residual < HUTI_TOLERANCE) {
            HUTI_INFO
            = HUTI_CONVERGENCE;
            break;
        }
        
        // And then 2n case
        
        [pcondrInvocation setArgument:&matrix atIndex:2];
        [pcondrInvocation setArgument:&t1v atIndex:3];
        [pcondrInvocation setArgument:&ynew atIndex:4];
        [pcondrInvocation setArgument:&ipar atIndex:5];
        [pcondrInvocation invoke];
        
        [matvecInvocation setArgument:&matrix atIndex:2];
        [matvecInvocation setArgument:&t1v atIndex:3];
        [matvecInvocation setArgument:&r atIndex:4];
        [matvecInvocation setArgument:&ipar atIndex:5];
        [matvecInvocation invoke];
        
        [pcondlInvocation setArgument:&matrix atIndex:2];
        [pcondlInvocation setArgument:&t1v atIndex:3];
        [pcondlInvocation setArgument:&r atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];
        
        for (i=0; i<ndim; i++) {
            w[i] = w[i] - alpha * t1v[i];
        }
        gamma = ( cblas_dnrm2(HUTI_NDIM, w, 1) ) / tau;
        c = 1 / sqrt( 1 + gamma * gamma );
        tau = tau * gamma * c;
        
        for (i=0; i<ndim; i++) {
            d[i] = ynew[i] + ( ( oldgamma * oldgamma * eta ) / alpha ) * d[i];
        }
        eta = c * c * alpha;
        for (i=0; i<ndim; i++) {
            x[i] = x[i] + eta * d[i];
        }
        
        oldgamma = gamma;
        
        // Check the convergence against selected stopping criterion
        switch (HUTI_STOPC) {
            case HUTI_TRUERESIDUAL:
                [pcondrInvocation setArgument:&matrix atIndex:2];
                [pcondrInvocation setArgument:&trv atIndex:3];
                [pcondrInvocation setArgument:&x atIndex:4];
                [pcondrInvocation setArgument:&ipar atIndex:5];
                [pcondrInvocation invoke];
                
                [matvecInvocation setArgument:&matrix atIndex:2];
                [matvecInvocation setArgument:&trv atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    trv[i] = r[i] - b[i];
                }
                [pcondlInvocation setArgument:&matrix atIndex:2];
                [pcondlInvocation setArgument:&r atIndex:3];
                [pcondlInvocation setArgument:&trv atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, r, 1);
                break;
            case HUTI_TRESID_SCALED_BYB:
                [pcondrInvocation setArgument:&matrix atIndex:2];
                [pcondrInvocation setArgument:&trv atIndex:3];
                [pcondrInvocation setArgument:&x atIndex:4];
                [pcondrInvocation setArgument:&ipar atIndex:5];
                [pcondrInvocation invoke];
                
                [matvecInvocation setArgument:&matrix atIndex:2];
                [matvecInvocation setArgument:&trv atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    trv[i] = r[i] - b[i];
                }
                [pcondlInvocation setArgument:&matrix atIndex:2];
                [pcondlInvocation setArgument:&r atIndex:3];
                [pcondlInvocation setArgument:&trv atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, r, 1) / rhsnorm;
                break;
            case HUTI_PSEUDORESIDUAL:
                [matvecInvocation setArgument:&matrix atIndex:2];
                [matvecInvocation setArgument:&x atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    r[i] = r[i] - b[i];
                }
                [pcondlInvocation setArgument:&matrix atIndex:2];
                [pcondlInvocation setArgument:&trv atIndex:3];
                [pcondlInvocation setArgument:&r atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, trv, 1);
                break;
            case HUTI_PRESID_SCALED_BYB:
                [matvecInvocation setArgument:&matrix atIndex:2];
                [matvecInvocation setArgument:&x atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    r[i] = r[i] - b[i];
                }
                [pcondlInvocation setArgument:&matrix atIndex:2];
                [pcondlInvocation setArgument:&trv atIndex:3];
                [pcondlInvocation setArgument:&r atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, trv, 1) / rhsnorm;
                break;
            case HUTI_PRESID_SCALED_BYPRECB:
                [matvecInvocation setArgument:&matrix atIndex:2];
                [matvecInvocation setArgument:&x atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    r[i] = r[i] - b[i];
                }
                [pcondlInvocation setArgument:&matrix atIndex:2];
                [pcondlInvocation setArgument:&trv atIndex:3];
                [pcondlInvocation setArgument:&r atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, trv, 1) / precrhsnorm;
                break;
            case HUTI_XDIFF_NORM:
                for (i=0; i<ndim; i++) {
                    r[i] = eta * d[i];
                }
                residual = cblas_dnrm2(HUTI_NDIM, r, 1);
                break;
            case HUTI_UPPERB_STOPC:
                upperb = sqrt( 2.0 * (double) iter_count ) * tau / rhsnorm;
                if ( (upperb / HUTI_TOLERANCE) < UPPERB_TOL_RATIO ) {
                    [pcondrInvocation setArgument:&matrix atIndex:2];
                    [pcondrInvocation setArgument:&trv atIndex:3];
                    [pcondrInvocation setArgument:&x atIndex:4];
                    [pcondrInvocation setArgument:&ipar atIndex:5];
                    [pcondrInvocation invoke];
                    
                    [matvecInvocation setArgument:&matrix atIndex:2];
                    [matvecInvocation setArgument:&trv atIndex:3];
                    [matvecInvocation setArgument:&r atIndex:4];
                    [matvecInvocation setArgument:&ipar atIndex:5];
                    [matvecInvocation invoke];
                    for (i=0; i<ndim; i++) {
                        trv[i] = r[i] - b[i];
                    }
                    [pcondlInvocation setArgument:&matrix atIndex:2];
                    [pcondlInvocation setArgument:&r atIndex:3];
                    [pcondlInvocation setArgument:&trv atIndex:4];
                    [pcondlInvocation setArgument:&ipar atIndex:5];
                    [pcondlInvocation invoke];
                    residual = cblas_dnrm2(HUTI_NDIM, r, 1) / rhsnorm;
                } else {
                    residual = upperb;
                }
                break;
            case HUTI_USUPPLIED_STOPC:
                [mstopInvocation setArgument:&matrix atIndex:2];
                [mstopInvocation setArgument:&x atIndex:3];
                [mstopInvocation setArgument:&b atIndex:4];
                [mstopInvocation setArgument:&r atIndex:5];
                [mstopInvocation setArgument:&ipar atIndex:6];
                [mstopInvocation invoke];
                [mstopInvocation getReturnValue:&residual];
                break;
            default:
                [pcondrInvocation setArgument:&matrix atIndex:2];
                [pcondrInvocation setArgument:&trv atIndex:3];
                [pcondrInvocation setArgument:&x atIndex:4];
                [pcondrInvocation setArgument:&ipar atIndex:5];
                [pcondrInvocation invoke];
                
                [matvecInvocation setArgument:&matrix atIndex:2];
                [matvecInvocation setArgument:&trv atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    trv[i] = r[i] - b[i];
                }
                [pcondlInvocation setArgument:&matrix atIndex:2];
                [pcondlInvocation setArgument:&r atIndex:3];
                [pcondlInvocation setArgument:&trv atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, r, 1);
                break;
        }

        if (residual < HUTI_TOLERANCE) {
            HUTI_INFO = HUTI_CONVERGENCE;
            break;
        }
        
        // Print debugging info if required
        if (HUTI_DBUGLVL != HUTI_NO_DEBUG) {
            if ( (iter_count % HUTI_DBUGLVL) == 0 ) {
                if (HUTI_STOPC == HUTI_UPPERB_STOPC) {
                    NSLog(@"FEMHUTIter:dtfqmrSolveInSolution: %d %11.4e %11.4e\n", iter_count, residual, upperb);
                } else {
                    NSLog(@"FEMHUTIter:dtfqmrSolveInSolution: %d %11.4e\n", iter_count, residual);
                }
            }
        }
        
        // Part 2C
        
        rho = cblas_ddot(HUTI_NDIM, rtld, 1, w, 1);
        beta = rho / oldrho;
        for (i=0; i<ndim; i++) {
            ynew[i] = w[i] + beta * ynew[i];
        }
        
        [pcondrInvocation setArgument:&matrix atIndex:2];
        [pcondrInvocation setArgument:&t2v atIndex:3];
        [pcondrInvocation setArgument:&ynew atIndex:4];
        [pcondrInvocation setArgument:&ipar atIndex:5];
        [pcondrInvocation invoke];
        
        [matvecInvocation setArgument:&matrix atIndex:2];
        [matvecInvocation setArgument:&t2v atIndex:3];
        [matvecInvocation setArgument:&r atIndex:4];
        [matvecInvocation setArgument:&ipar atIndex:5];
        [matvecInvocation invoke];
        
        [pcondlInvocation setArgument:&matrix atIndex:2];
        [pcondlInvocation setArgument:&t2v atIndex:3];
        [pcondlInvocation setArgument:&r atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];

        // Note: we still have MATRIX * YNEW in t1v
        
        for (i=0; i<ndim; i++) {
            v[i] = t2v[i] + beta * t1v[i] + beta * beta * v[i];
            y[i] = ynew[i];
        }
        
        oldrho = rho;
        
        // Return back to the top of the iteration loop (without initialization)
        iter_count = iter_count + 1;
        if (iter_count > HUTI_MAXIT) {
            HUTI_INFO = HUTI_MAXITER;
            break;
        }
        
    }
    
 // We have exited the loop after enough iterations or broke down
jump:
    
    // Compute the unpreconditioned xvec
    [pcondrInvocation setArgument:&matrix atIndex:2];
    [pcondrInvocation setArgument:&trv atIndex:3];
    [pcondrInvocation setArgument:&x atIndex:4];
    [pcondrInvocation setArgument:&ipar atIndex:5];
    [pcondrInvocation invoke];
    for (i=0; i<ndim; i++) {
        x[i] = trv[i];
    }
    
    if (HUTI_DBUGLVL != HUTI_NO_DEBUG) {
        if (HUTI_STOPC == HUTI_UPPERB_STOPC) {
            NSLog(@"FEMHUTIter:dtfqmrSolveInSolution: %d %11.4e %11.4e\n", iter_count, residual, upperb);
        } else {
            NSLog(@"FEMHUTIter:dtfqmrSolveInSolution: %d %11.4e\n", iter_count, residual);
        }
    }
    
    HUTI_ITERS = iter_count;
    
    // Release memory
    free_dvector(v, 0, ndim-1);
    free_dvector(y, 0, ndim-1);
    free_dvector(ynew, 0, ndim-1);
    free_dvector(rtld, 0, ndim-1);
    free_dvector(t1v, 0, ndim-1);
    free_dvector(t2v, 0, ndim-1);
    free_dvector(w, 0, ndim-1);
    free_dvector(d, 0, ndim-1);
    free_dvector(r, 0, ndim-1);
    free_dvector(trv, 0, ndim-1);
}

#pragma mark CG

/**********************************************************************************************************************************
 
    This method is based on Barret et al. book: "Templates for the Solution of Linear Systems: Building blocks for
    Iterative Methods", 1993.
 
    All matrix-vector operations are done externally, so we don't need to know about the matrix structure (sparse or dense).
 
    Double precision version.
 
**********************************************************************************************************************************/
-(void)dcgSolveInSolution:(FEMSolution *)solution matrix:(FEMMatrix *)matrix ndim:(int)ndim result:(double *)x rhs:(double *)b ipar:(int *)ipar dpar:(double *)dpar pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod {
    
    int i, iter_count;
    double rho, oldrho, alpha, beta;
    double *z, *p, *q, *r;
    double residual, rhsnorm, precrhsnorm;
    
    NSMethodSignature *pCondlSignature, *pCondrSignature, *matvecSignature, *mstopSignature;
    NSInvocation *pcondlInvocation, *pcondrInvocation, *matvecInvocation, *mstopInvocation;
    
    FEMPrecondition *preconditioning;
    FEMCore *core;
    
    // Acquire signature invocations and set selector for invocations
    pCondlSignature = [FEMPrecondition instanceMethodSignatureForSelector:pcondlMethod];
    pcondlInvocation = [NSInvocation invocationWithMethodSignature:pCondlSignature];
    [pcondlInvocation setSelector:pcondlMethod];
    
    pCondrSignature = [FEMPrecondition instanceMethodSignatureForSelector:pcondrMethod];
    pcondrInvocation = [NSInvocation invocationWithMethodSignature:pCondrSignature];
    [pcondrInvocation setSelector:pcondrMethod];
    
    matvecSignature = [FEMPrecondition instanceMethodSignatureForSelector:matvecMethod];
    matvecInvocation = [NSInvocation invocationWithMethodSignature:matvecSignature];
    [matvecInvocation setSelector:matvecMethod];
    
    mstopSignature = [FEMCore instanceMethodSignatureForSelector:mstopMethod];
    mstopInvocation = [NSInvocation invocationWithMethodSignature:mstopSignature];
    [mstopInvocation setSelector:mstopMethod];
    
    // Instanciate class for the preconditioning
    preconditioning = [[FEMPrecondition alloc] init];
    
    //Instanciate class for the kernel
    core = [FEMCore sharedCore];
    
    // Targets to invocations
    [pcondlInvocation setTarget:preconditioning];
    [pcondrInvocation setTarget:preconditioning];
    [matvecInvocation setTarget:preconditioning];
    [mstopInvocation setTarget:core];
    
    // Vector allocations
    z = doublevec(0, ndim-1);
    p = doublevec(0, ndim-1);
    q = doublevec(0, ndim-1);
    r = doublevec(0, ndim-1);
    
    // The actual CG begins here (look at the pseudo-code in the book, page 27)
    
    // First the initialization part
    iter_count = 1;
    
    // Norms of right-hand side vector are used in convergence tests
    if (HUTI_STOPC == HUTI_TRESID_SCALED_BYB || HUTI_STOPC == HUTI_PRESID_SCALED_BYB) {
        rhsnorm = cblas_dnrm2(HUTI_NDIM, b, 1);
    }
    if (HUTI_STOPC == HUTI_PRESID_SCALED_BYPRECB) {
        [pcondlInvocation setArgument:&matrix atIndex:2];
        [pcondlInvocation setArgument:&p atIndex:3];
        [pcondlInvocation setArgument:&b atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];
        
        precrhsnorm = cblas_dnrm2(HUTI_NDIM, p, 1);
    }
    
    // The following applies for all matrix operations in this solver
    HUTI_EXTOP_MATTYPE = HUTI_MAT_NOTTRPSED;
    
    // Generate vector xvec if needed
    if (HUTI_INITIALX == HUTI_RANDOMX) {
        [self HUTI_dRandVector:x ipar:ipar];
    } else if (HUTI_INITIALX != HUTI_USERSUPPLIEDX) {
        for (i=0; i<ndim; i++) {
            x[i] = 1;
        }
    }
    
    [matvecInvocation setArgument:&matrix atIndex:2];
    [matvecInvocation setArgument:&x atIndex:3];
    [matvecInvocation setArgument:&r atIndex:4];
    [matvecInvocation setArgument:&ipar atIndex:5];
    [matvecInvocation invoke];
    
    for (i=0; i<ndim; i++) {
        r[i] = b[i] - r[i];
    }
    
    // This is where the loop starts
    while (1) {
        
        [pcondlInvocation setArgument:&matrix atIndex:2];
        [pcondlInvocation setArgument:&q atIndex:3];
        [pcondlInvocation setArgument:&r atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];
        
        [pcondrInvocation setArgument:&matrix atIndex:2];
        [pcondrInvocation setArgument:&z atIndex:3];
        [pcondrInvocation setArgument:&q atIndex:4];
        [pcondrInvocation setArgument:&ipar atIndex:5];
        [pcondrInvocation invoke];
        
        rho = cblas_ddot(HUTI_NDIM, r, 1, z, 1);
        if (rho == 0.0) {
            HUTI_INFO = HUTI_CG_RHO;
            break;
        }
        
        if (iter_count == 1) {
            for (i=0; i<ndim; i++) {
                p[i] = z[i];
            }
        } else {
            beta = rho / oldrho;
            for (i=0; i<ndim; i++) {
                p[i] = z[i] + beta * p[i];
            }
        }
        
        [matvecInvocation setArgument:&matrix atIndex:2];
        [matvecInvocation setArgument:&p atIndex:3];
        [matvecInvocation setArgument:&q atIndex:4];
        [matvecInvocation setArgument:&ipar atIndex:5];
        [matvecInvocation invoke];
        
        alpha = rho / cblas_ddot(HUTI_NDIM, p, 1, q, 1);
        
        for (i=0; i<ndim; i++) {
            x[i] = x[i] + alpha * p[i];
            r[i] = r[i] - alpha * q[i];
        }
        
        // Check the convergence against selected stopping criterion
        switch (HUTI_STOPC) {
            case HUTI_TRUERESIDUAL:
                [matvecInvocation setArgument:&matrix atIndex:2];
                [matvecInvocation setArgument:&x atIndex:3];
                [matvecInvocation setArgument:&z atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    z[i] = z[i] - b[i];
                }
                residual = cblas_dnrm2(HUTI_NDIM, z, 1);
                break;
            case HUTI_TRESID_SCALED_BYB:
                [matvecInvocation setArgument:&matrix atIndex:2];
                [matvecInvocation setArgument:&x atIndex:3];
                [matvecInvocation setArgument:&z atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    z[i] = z[i] - b[i];
                }
                residual = cblas_dnrm2(HUTI_NDIM, z, 1) / rhsnorm;
                break;
            case HUTI_PSEUDORESIDUAL:
                residual = cblas_dnrm2(HUTI_NDIM, r, 1);
                break;
            case HUTI_PRESID_SCALED_BYB:
                residual = cblas_dnrm2(HUTI_NDIM, r, 1) / rhsnorm;
                break;
            case HUTI_PRESID_SCALED_BYPRECB:
                residual = cblas_dnrm2(HUTI_NDIM, r, 1) / precrhsnorm;
                break;
            case HUTI_XDIFF_NORM:
                for (i=0; i<ndim; i++) {
                    z[i] = alpha * p[i];
                }
                residual = cblas_dnrm2(HUTI_NDIM, z, 1);
                break;
            case HUTI_USUPPLIED_STOPC:
                [mstopInvocation setArgument:&matrix atIndex:2];
                [mstopInvocation setArgument:&x atIndex:3];
                [mstopInvocation setArgument:&b atIndex:4];
                [mstopInvocation setArgument:&r atIndex:5];
                [mstopInvocation setArgument:&ipar atIndex:6];
                [mstopInvocation invoke];
                [mstopInvocation getReturnValue:&residual];
                break;
            default:
                [matvecInvocation setArgument:&matrix atIndex:2];
                [matvecInvocation setArgument:&x atIndex:3];
                [matvecInvocation setArgument:&z atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    z[i] = z[i] - b[i];
                }
                residual = cblas_dnrm2(HUTI_NDIM, z, 1);
                break;
        }
        
        // Print debugging info if required
        if (HUTI_DBUGLVL != HUTI_NO_DEBUG) {
            if ( (iter_count % HUTI_DBUGLVL) == 0 ) {
                NSLog(@"FEMHUTIter:dcgSolveInSolution: %d %11.4e\n", iter_count, residual);
            }
        }
        
        if (residual < HUTI_TOLERANCE) {
            HUTI_INFO = HUTI_CONVERGENCE;
            break;
        }
        
        oldrho = rho;
        
        // Return back to the top of the iteration loop (without initialization)
        iter_count = iter_count + 1;
        if (iter_count > HUTI_MAXIT) {
            HUTI_INFO = HUTI_MAXITER;
            break;
        }
        
    }
    
    // We have exited the loop after enough iterations or broke down
    if (HUTI_DBUGLVL != HUTI_NO_DEBUG) {
        NSLog(@"FEMHUTIter:dcgSolveInSolution: %d %11.4e\n", iter_count, residual);
    }
    
    HUTI_ITERS = iter_count;
    
    // Release memory
    free_dvector(z, 0, ndim-1);
    free_dvector(p, 0, ndim-1);
    free_dvector(q, 0, ndim-1);
    free_dvector(r, 0, ndim-1);
}

#pragma mark CGS

/**********************************************************************************************************************************
 
    This method is based on Barret et al. book: "Templates for the Solution of Linear Systems: Building blocks for
    Iterative Methods", 1993.
 
    All matrix-vector operations are done externally, so we don't need to know about the matrix structure (sparse or dense).
    Memory allocationf for the work array has also been done externally.
 
    Double precision version.
 
**********************************************************************************************************************************/
-(void)dcgsSolveInSolution:(FEMSolution *)solution matrix:(FEMMatrix *)matrix ndim:(int)ndim result:(double *)x rhs:(double *)b ipar:(int *)ipar dpar:(double *)dpar pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod {
    
    int i, iter_count;
    double rho, oldrho, alpha, beta;
    double *rtld, *p, *q, *u, *t1v, *t2v, *r;
    double residual, rhsnorm, precrhsnorm;
        
    NSMethodSignature *pCondlSignature, *pCondrSignature, *matvecSignature, *mstopSignature;
    NSInvocation *pcondlInvocation, *pcondrInvocation, *matvecInvocation, *mstopInvocation;
    
    FEMPrecondition *preconditioning;
    FEMCore *core;
    
    // Acquire signature invocations and set selector for invocations
    pCondlSignature = [FEMPrecondition instanceMethodSignatureForSelector:pcondlMethod];
    pcondlInvocation = [NSInvocation invocationWithMethodSignature:pCondlSignature];
    [pcondlInvocation setSelector:pcondlMethod];
    
    pCondrSignature = [FEMPrecondition instanceMethodSignatureForSelector:pcondrMethod];
    pcondrInvocation = [NSInvocation invocationWithMethodSignature:pCondrSignature];
    [pcondrInvocation setSelector:pcondrMethod];
    
    matvecSignature = [FEMPrecondition instanceMethodSignatureForSelector:matvecMethod];
    matvecInvocation = [NSInvocation invocationWithMethodSignature:matvecSignature];
    [matvecInvocation setSelector:matvecMethod];
    
    mstopSignature = [FEMCore instanceMethodSignatureForSelector:mstopMethod];
    mstopInvocation = [NSInvocation invocationWithMethodSignature:mstopSignature];
    [mstopInvocation setSelector:mstopMethod];
    
    // Instanciate class for the preconditioning
    preconditioning = [[FEMPrecondition alloc] init];
    
    //Instanciate class for the kernel
    core = [FEMCore sharedCore];
    
    // Targets to invocations
    [pcondlInvocation setTarget:preconditioning];
    [pcondrInvocation setTarget:preconditioning];
    [matvecInvocation setTarget:preconditioning];
    [mstopInvocation setTarget:core];

    // Vector allocations
    rtld = doublevec(0, ndim-1);
    p = doublevec(0, ndim-1);
    q = doublevec(0, ndim-1);
    u = doublevec(0, ndim-1);
    t1v = doublevec(0, ndim-1);
    t2v = doublevec(0, ndim-1);
    r = doublevec(0, ndim-1);
    
    // The actual CGS begins here (look at the pseudo-code in the book, page 27)
    
    // First the initialization part
    iter_count = 1;
    
    // The following applies for all matrix operations in this solver
    HUTI_EXTOP_MATTYPE = HUTI_MAT_NOTTRPSED;
    
    // Norms of right-hand side vector are used in convergence tests
    if (HUTI_STOPC == HUTI_TRESID_SCALED_BYB || HUTI_STOPC == HUTI_PRESID_SCALED_BYB) {
        rhsnorm = cblas_dnrm2(HUTI_NDIM, b, 1);
    }
    if (HUTI_STOPC == HUTI_PRESID_SCALED_BYPRECB) {
        [pcondlInvocation setArgument:&matrix atIndex:2];
        [pcondlInvocation setArgument:&p atIndex:3];
        [pcondlInvocation setArgument:&b atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];
        
        precrhsnorm = cblas_dnrm2(HUTI_NDIM, p, 1);
    }
    
    // Generate vector xvec if needed
    if (HUTI_INITIALX == HUTI_RANDOMX) {
        [self HUTI_dRandVector:x ipar:ipar];
    } else if (HUTI_INITIALX != HUTI_USERSUPPLIEDX) {
        for (i=0; i<ndim; i++) {
            x[i] = 1;
        }
    }
    
    [matvecInvocation setArgument:&matrix atIndex:2];
    [matvecInvocation setArgument:&x atIndex:3];
    [matvecInvocation setArgument:&r atIndex:4];
    [matvecInvocation setArgument:&ipar atIndex:5];
    [matvecInvocation invoke];
    
    for (i=0; i<ndim; i++) {
        r[i] = b[i] - r[i];
        rtld[i] = r[i];
    }
    
    // This is where the loop starts
    while (1) {
        
        rho = cblas_ddot(HUTI_NDIM, rtld, 1, r, 1);
        if (rho == 0.0) {
            HUTI_INFO = HUTI_CGS_RHO;
            break;
        }
        
        if (iter_count == 1) {
            for (i=0; i<ndim; i++) {
                u[i] = r[i];
                p[i] = u[i];
            }
        } else {
            beta = rho / oldrho;
            for (i=0; i<ndim; i++) {
                u[i] = r[i] + beta * q[i];
                p[i] = u[i] + beta * q[i] + beta * beta * p[i];
            }
        }
        
        [pcondlInvocation setArgument:&matrix atIndex:2];
        [pcondlInvocation setArgument:&t2v atIndex:3];
        [pcondlInvocation setArgument:&p atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];
        
        [pcondrInvocation setArgument:&matrix atIndex:2];
        [pcondrInvocation setArgument:&t1v atIndex:3];
        [pcondrInvocation setArgument:&t2v atIndex:4];
        [pcondrInvocation setArgument:&ipar atIndex:5];
        [pcondrInvocation invoke];
        
        [matvecInvocation setArgument:&matrix atIndex:2];
        [matvecInvocation setArgument:&t1v atIndex:3];
        [matvecInvocation setArgument:&t2v atIndex:4];
        [matvecInvocation setArgument:&ipar atIndex:5];
        [matvecInvocation invoke];

        alpha = rho / cblas_ddot(HUTI_NDIM, rtld, 1, t2v, 1);
        for (i=0; i<ndim; i++) {
            q[i] = u[i] - alpha * t2v[i];
        }
        for (i=0; i<ndim; i++) {
            t2v[i] = u[i] + q[i];
        }
        
        [pcondlInvocation setArgument:&matrix atIndex:2];
        [pcondlInvocation setArgument:&u atIndex:3];
        [pcondlInvocation setArgument:&t2v atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];
        
        [pcondrInvocation setArgument:&matrix atIndex:2];
        [pcondrInvocation setArgument:&t1v atIndex:3];
        [pcondrInvocation setArgument:&u atIndex:4];
        [pcondrInvocation setArgument:&ipar atIndex:5];
        [pcondrInvocation invoke];
        for (i=0; i<ndim; i++) {
            x[i] = x[i] + alpha * t1v[i];
        }
        
        [matvecInvocation setArgument:&matrix atIndex:2];
        [matvecInvocation setArgument:&t1v atIndex:3];
        [matvecInvocation setArgument:&t2v atIndex:4];
        [matvecInvocation setArgument:&ipar atIndex:5];
        [matvecInvocation invoke];
        for (i=0; i<ndim; i++) {
            r[i] = r[i] - alpha * t2v[i];
        }
        
        // Check the convergence against selected stopping criterion
        switch (HUTI_STOPC) {
            case HUTI_TRUERESIDUAL:
                [matvecInvocation setArgument:&matrix atIndex:2];
                [matvecInvocation setArgument:&x atIndex:3];
                [matvecInvocation setArgument:&t1v atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    t1v[i] = t1v[i] - b[i];
                }
                residual = cblas_dnrm2(HUTI_NDIM, t1v, 1);
                break;
            case HUTI_TRESID_SCALED_BYB:
                [matvecInvocation setArgument:&matrix atIndex:2];
                [matvecInvocation setArgument:&x atIndex:3];
                [matvecInvocation setArgument:&t1v atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    t1v[i] = t1v[i] - b[i];
                }
                residual = cblas_dnrm2(HUTI_NDIM, t1v, 1) / rhsnorm;
                break;
            case HUTI_PSEUDORESIDUAL:
                residual = cblas_dnrm2(HUTI_NDIM, r, 1);
                break;
            case HUTI_PRESID_SCALED_BYB:
                residual = cblas_dnrm2(HUTI_NDIM, r, 1) / rhsnorm;
                break;
            case HUTI_PRESID_SCALED_BYPRECB:
                residual = cblas_dnrm2(HUTI_NDIM, r, 1) / precrhsnorm;
                break;
            case HUTI_XDIFF_NORM:
                for (i=0; i<ndim; i++) {
                    t1v[i] = alpha * t1v[i];
                }
                residual = cblas_dnrm2(HUTI_NDIM, t1v, 1);
                break;
            case HUTI_USUPPLIED_STOPC:
                [mstopInvocation setArgument:&matrix atIndex:2];
                [mstopInvocation setArgument:&x atIndex:3];
                [mstopInvocation setArgument:&b atIndex:4];
                [mstopInvocation setArgument:&r atIndex:5];
                [mstopInvocation setArgument:&ipar atIndex:6];
                [mstopInvocation invoke];
                [mstopInvocation getReturnValue:&residual];
                break;
            default:
                [matvecInvocation setArgument:&matrix atIndex:2];
                [matvecInvocation setArgument:&x atIndex:3];
                [matvecInvocation setArgument:&t1v atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    t1v[i] = t1v[i] - b[i];
                }
                residual = cblas_dnrm2(HUTI_NDIM, t1v, 1);
                break;
        }
        
        // Print debugging info if required
        if (HUTI_DBUGLVL != HUTI_NO_DEBUG) {
            if ( (iter_count % HUTI_DBUGLVL) == 0 ) {
                NSLog(@"FEMHUTIter:dcgsSolveInSolution: %d %11.4e\n", iter_count, residual);
            }
        }
        
        if (residual < HUTI_TOLERANCE) {
            HUTI_INFO = HUTI_CONVERGENCE;
            break;
        }

        oldrho = rho;
        
        // Return back to the top of the iteration loop (without initialization)
        iter_count = iter_count + 1;
        if (iter_count > HUTI_MAXIT) {
            HUTI_INFO = HUTI_MAXITER;
            break;
        }
        
    }
    
    // We have exited the loop after enough iterations or broke down
    if (HUTI_DBUGLVL != HUTI_NO_DEBUG) {
        NSLog(@"FEMHUTIter:dcgsSolveInSolution: %d %11.4e\n", iter_count, residual);
    }
    
    HUTI_ITERS = iter_count;
    
    // Release memory
    free_dvector(rtld, 0, ndim-1);
    free_dvector(p, 0, ndim-1);
    free_dvector(q, 0, ndim-1);
    free_dvector(u, 0, ndim-1);
    free_dvector(t1v, 0, ndim-1);
    free_dvector(t2v, 0, ndim-1);
    free_dvector(r, 0, ndim-1);
}

#pragma mark GMRES

/**********************************************************************************************************************************
 
    This method is based on Barret et al. book: "Templates for the Solution of Linear Systems: Building blocks for
    Iterative Methods", 1993.
 
    All matrix-vector operations are done externally, so we don't need to know about the matrix structure (sparse or dense).
    Memory allocationf for the work array has also been done externally.
 
    Double precision version.
 
**********************************************************************************************************************************/
-(void)dgmresSolveInSolution:(FEMSolution *)solution matrix:(FEMMatrix *)matrix ndim:(int)ndim wrkdim:(int)wrkdim result:(double *)x rhs:(double *)b ipar:(int *)ipar dpar:(double *)dpar work:(double **)work pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod {
    
    int i, j, k, iter_count, m;
    int s_ind, vtmp_ind, v_ind;
    double bnrm, alpha, beta;
    double temp, temp2, error;
    double *w, *r, *s, *t1v, *v, *buffer;
    double **h, *cs, *sn, *y; 
    double residual, rhsnorm, precrhsnorm;
    
    NSMethodSignature *pCondlSignature, *pCondrSignature, *matvecSignature, *mstopSignature;
    NSInvocation *pcondlInvocation, *pcondrInvocation, *matvecInvocation, *mstopInvocation;
    
    FEMPrecondition *preconditioning;
    FEMCore *core;
    
    // Acquire signature invocations and set selector for invocations
    pCondlSignature = [FEMPrecondition instanceMethodSignatureForSelector:pcondlMethod];
    pcondlInvocation = [NSInvocation invocationWithMethodSignature:pCondlSignature];
    [pcondlInvocation setSelector:pcondlMethod];
    
    pCondrSignature = [FEMPrecondition instanceMethodSignatureForSelector:pcondrMethod];
    pcondrInvocation = [NSInvocation invocationWithMethodSignature:pCondrSignature];
    [pcondrInvocation setSelector:pcondrMethod];
    
    matvecSignature = [FEMPrecondition instanceMethodSignatureForSelector:matvecMethod];
    matvecInvocation = [NSInvocation invocationWithMethodSignature:matvecSignature];
    [matvecInvocation setSelector:matvecMethod];
    
    mstopSignature = [FEMCore instanceMethodSignatureForSelector:mstopMethod];
    mstopInvocation = [NSInvocation invocationWithMethodSignature:mstopSignature];
    [mstopInvocation setSelector:mstopMethod];
    
    // Instanciate class for the preconditioning
    preconditioning = [[FEMPrecondition alloc] init];
    
    //Instanciate class for the kernel
    core = [FEMCore sharedCore];
    
    // Targets to invocations
    [pcondlInvocation setTarget:preconditioning];
    [pcondrInvocation setTarget:preconditioning];
    [matvecInvocation setTarget:preconditioning];
    [mstopInvocation setTarget:core];
    
    // Vector allocations
    w = doublevec(0, ndim-1);
    r = doublevec(0, ndim-1);
    s = doublevec(0, ndim-1);
    t1v = doublevec(0, ndim-1);
    v = doublevec(0, ndim-1);
    
    buffer = doublevec(0, ndim-1);
    
    cs = doublevec(0, (HUTI_GMRES_RESTART+1)-1);
    sn = doublevec(0, (HUTI_GMRES_RESTART+1)-1);
    y = doublevec(0, (HUTI_GMRES_RESTART+1)-1);
    
    h = doublematrix(0, (HUTI_GMRES_RESTART+1)-1, 0, (HUTI_GMRES_RESTART+1)-1);
    
    // The actual GMRES begins here (look at the pseudo-code in the book, page 27)
    
    // First the initialization part
    iter_count = 1;
    
    bnrm = cblas_dnrm2(HUTI_NDIM, b, 1);
    
    // Norms of right-hand side vector are used in convergence tests
    if (HUTI_STOPC == HUTI_TRESID_SCALED_BYB || HUTI_STOPC == HUTI_PRESID_SCALED_BYB) {
        rhsnorm = bnrm;
    }
    if (HUTI_STOPC == HUTI_PRESID_SCALED_BYPRECB) {
        [pcondlInvocation setArgument:&matrix atIndex:2];
        [pcondlInvocation setArgument:&t1v atIndex:3];
        [pcondlInvocation setArgument:&b atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];
        
        precrhsnorm = cblas_dnrm2(HUTI_NDIM, t1v, 1);
    }
    
    // The following applies for all matrix operations in this solver
    HUTI_EXTOP_MATTYPE = HUTI_MAT_NOTTRPSED;

    // Generate vector xvec if needed
    if (HUTI_INITIALX == HUTI_RANDOMX) {
        [self HUTI_dRandVector:x ipar:ipar];
    } else if (HUTI_INITIALX != HUTI_USERSUPPLIEDX) {
        for (i=0; i<ndim; i++) {
            x[i] = 1;
        }
    }
    
    [pcondrInvocation setArgument:&matrix atIndex:2];
    [pcondrInvocation setArgument:&t1v atIndex:3];
    [pcondrInvocation setArgument:&x atIndex:4];
    [pcondrInvocation setArgument:&ipar atIndex:5];
    [pcondrInvocation invoke];
    
    [matvecInvocation setArgument:&matrix atIndex:2];
    [matvecInvocation setArgument:&t1v atIndex:3];
    [matvecInvocation setArgument:&r atIndex:4];
    [matvecInvocation setArgument:&ipar atIndex:5];
    [matvecInvocation invoke];
    for (i=0; i<ndim; i++) {
        t1v[i] = b[i] - r[i];
    }
    [pcondlInvocation setArgument:&matrix atIndex:2];
    [pcondlInvocation setArgument:&r atIndex:3];
    [pcondlInvocation setArgument:&t1v atIndex:4];
    [pcondlInvocation setArgument:&ipar atIndex:5];
    [pcondlInvocation invoke];
    
    m = HUTI_GMRES_RESTART;
    s_ind = 3;
    vtmp_ind = 4;
    v_ind = 6;
    for (i=(v_ind+1-1)-1; i<v_ind+m+1-1; i++) {
        for (j=0; j<ndim; j++) {
            work[j][i] = 0.0;
        }
    }
    
    memset( *h, 0.0, ((HUTI_GMRES_RESTART+1)*(HUTI_GMRES_RESTART+1))*sizeof(double) );
    memset( cs, 0.0, (HUTI_GMRES_RESTART+1)*sizeof(double) );
    memset( sn, 0.0, (HUTI_GMRES_RESTART+1)*sizeof(double) );
    work[0][vtmp_ind-1] = 1.0;
    
    // This is where the loop starts
    while (1) {
        
        [pcondrInvocation setArgument:&matrix atIndex:2];
        [pcondrInvocation setArgument:&t1v atIndex:3];
        [pcondrInvocation setArgument:&x atIndex:4];
        [pcondrInvocation setArgument:&ipar atIndex:5];
        [pcondrInvocation invoke];
        
        [matvecInvocation setArgument:&matrix atIndex:2];
        [matvecInvocation setArgument:&t1v atIndex:3];
        [matvecInvocation setArgument:&r atIndex:4];
        [matvecInvocation setArgument:&ipar atIndex:5];
        [matvecInvocation invoke];
        for (i=0; i<ndim; i++) {
            t1v[i] = b[i] - r[i];
        }
        [pcondlInvocation setArgument:&matrix atIndex:2];
        [pcondlInvocation setArgument:&r atIndex:3];
        [pcondlInvocation setArgument:&t1v atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];
        
        alpha = cblas_dnrm2(HUTI_NDIM, r, 1);
        if (alpha == 0.0) {
            HUTI_INFO = HUTI_GMRES_ALPHA;
            break;
        }

        for (i=0; i<ndim; i++) {
            work[i][(v_ind+1-1)-1] = r[i] / alpha;
            s[i] = alpha * work[i][vtmp_ind-1];
        }
        
        // Construct orthonormal
        
        for (i=0; i<m; i++) {
            
            for (j=0; j<ndim; j++) {
                buffer[j] = work[j][v_ind+i-1];
            }
            
            [pcondrInvocation setArgument:&matrix atIndex:2];
            [pcondrInvocation setArgument:&w atIndex:3];
            [pcondrInvocation setArgument:&buffer atIndex:4];
            [pcondrInvocation setArgument:&ipar atIndex:5];
            [pcondrInvocation invoke];
            
            [matvecInvocation setArgument:&matrix atIndex:2];
            [matvecInvocation setArgument:&w atIndex:3];
            [matvecInvocation setArgument:&t1v atIndex:4];
            [matvecInvocation setArgument:&ipar atIndex:5];
            [matvecInvocation invoke];
            
            [pcondlInvocation setArgument:&matrix atIndex:2];
            [pcondlInvocation setArgument:&w atIndex:3];
            [pcondlInvocation setArgument:&t1v atIndex:4];
            [pcondlInvocation setArgument:&ipar atIndex:5];
            [pcondlInvocation invoke];

            for (k=0; k<=i; k++) {
                
                for (j=0; j<ndim; j++) {
                    buffer[j] = work[j][v_ind+k-1];
                }
                h[k][i] = cblas_ddot(HUTI_NDIM, w, 1, buffer, 1);
                for (j=0; j<ndim; j++) {
                    w[j] = w[j] - h[k][i] * work[j][v_ind+k-1];
                }
            }
            
            beta = cblas_dnrm2(HUTI_NDIM, w, 1);
            if (beta == 0.0) {
                HUTI_INFO = HUTI_GMRES_BETA;
                goto jump;
            }
            
            h[i+1][i] = beta;
            for (j=0; j<ndim; j++) {
                work[j][v_ind+i+1-1] = w[j] / h[i+1][i];
            }
            
            // Compute the givens rotation
            
            for (k=0; k<=i-1; k++) {
                
                temp = cs[k] * h[k][i] + sn[k] * h[k+1][i];
                h[k+1][i] = -1.0 * sn[k] * h[k][i] + cs[k] * h[k+1][i];
                h[k][i] = temp;
            }
            
            if (h[i+1][i] == 0.0) {
                cs[i] = 1.0;
                sn[i] = 0.0;
            } else {
                if ( fabs(h[i+1][i]) > fabs(h[i][i]) ) {
                    temp2 = h[i][i] / h[i+1][i];
                    sn[i] = 1.0 / sqrt( 1.0 + (temp2 * temp2) );
                    cs[i] = temp2 * sn[i];
                } else {
                    temp2 = h[i+1][i] / h[i][i];
                    cs[i] = 1.0 / sqrt( 1.0 + (temp2 * temp2) );
                    sn[i] = temp2 * cs[i];
                }
            }

            temp = cs[i] * s[i];
            s[i+1] = -1.0 * sn[i] * s[i];
            s[i] = temp;
            
            h[i][i] = ( cs[i] * h[i][i] ) + ( sn[i] * h[i+1][i] );
            h[i+1][i] = 0.0;
            
            error = fabs(s[i+1]) / bnrm;
            if (error < HUTI_TOLERANCE) {
                
                [self HUTI_dLuSolveAt:i luMatrix:h afterSolve:y rightHandSide:s];
                
                // Matrix-vector product
                cblas_dgemv(CblasRowMajor, CblasNoTrans, matrix.numberOfRows, ((v_ind+i)-(v_ind+1-1))+1, 1.0, *work+(v_ind-1), wrkdim, y, 1, 1.0, x, 1);
                break;
            }
        }
        
        if (error < HUTI_TOLERANCE) {
            goto error_lower_than_tolerance;
        }
        
        [self HUTI_dLuSolveAt:m-1 luMatrix:h afterSolve:y rightHandSide:s];
        // Matrix-vector product
        cblas_dgemv(CblasRowMajor, CblasNoTrans, matrix.numberOfRows, ((v_ind+m-1)-(v_ind+1-1))+1, 1.0, *work+(v_ind-1), wrkdim, y, 1, 1.0, x, 1);
        
    error_lower_than_tolerance:
        
        // Check the convergence against selected stopping criterion
        switch (HUTI_STOPC) {
            case HUTI_TRUERESIDUAL:
                [pcondrInvocation setArgument:&matrix atIndex:2];
                [pcondrInvocation setArgument:&t1v atIndex:3];
                [pcondrInvocation setArgument:&x atIndex:4];
                [pcondrInvocation setArgument:&ipar atIndex:5];
                [pcondrInvocation invoke];

                [matvecInvocation setArgument:&matrix atIndex:2];
                [matvecInvocation setArgument:&t1v atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    t1v[i] = b[i] - r[i];
                }
                [pcondlInvocation setArgument:&matrix atIndex:2];
                [pcondlInvocation setArgument:&r atIndex:3];
                [pcondlInvocation setArgument:&t1v atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, r, 1);
                break;
            case HUTI_TRESID_SCALED_BYB:
                [pcondrInvocation setArgument:&matrix atIndex:2];
                [pcondrInvocation setArgument:&t1v atIndex:3];
                [pcondrInvocation setArgument:&x atIndex:4];
                [pcondrInvocation setArgument:&ipar atIndex:5];
                [pcondrInvocation invoke];

                [matvecInvocation setArgument:&matrix atIndex:2];
                [matvecInvocation setArgument:&t1v atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    t1v[i] = b[i] - r[i];
                }
                [pcondlInvocation setArgument:&matrix atIndex:2];
                [pcondlInvocation setArgument:&r atIndex:3];
                [pcondlInvocation setArgument:&t1v atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, r, 1) / rhsnorm;
                break;
            case HUTI_PSEUDORESIDUAL:
                [matvecInvocation setArgument:&matrix atIndex:2];
                [matvecInvocation setArgument:&x atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    r[i] = r[i] - b[i];
                }
                [pcondlInvocation setArgument:&matrix atIndex:2];
                [pcondlInvocation setArgument:&t1v atIndex:3];
                [pcondlInvocation setArgument:&r atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, t1v, 1);
                break;
            case HUTI_PRESID_SCALED_BYB:
                [matvecInvocation setArgument:&matrix atIndex:2];
                [matvecInvocation setArgument:&x atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    r[i] = r[i] - b[i];
                }
                [pcondlInvocation setArgument:&matrix atIndex:2];
                [pcondlInvocation setArgument:&t1v atIndex:3];
                [pcondlInvocation setArgument:&r atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, t1v, 1) / rhsnorm;
                break;
            case HUTI_PRESID_SCALED_BYPRECB:
                [matvecInvocation setArgument:&matrix atIndex:2];
                [matvecInvocation setArgument:&x atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    r[i] = r[i] - b[i];
                }
                [pcondlInvocation setArgument:&matrix atIndex:2];
                [pcondlInvocation setArgument:&t1v atIndex:3];
                [pcondlInvocation setArgument:&r atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, t1v, 1) / precrhsnorm;
                break;
            case HUTI_XDIFF_NORM:
                // Matrix-vector product
                cblas_dgemv(CblasRowMajor, CblasNoTrans, matrix.numberOfRows, ((v_ind+m-1)-(v_ind+1-1))+1, 1.0, *work+(v_ind-1), wrkdim, y, 1, 0.0, t1v, 1);
                residual = cblas_dnrm2(HUTI_NDIM, t1v, 1);
                break;
            case HUTI_USUPPLIED_STOPC:
                [mstopInvocation setArgument:&matrix atIndex:2];
                [mstopInvocation setArgument:&x atIndex:3];
                [mstopInvocation setArgument:&b atIndex:4];
                [mstopInvocation setArgument:&r atIndex:5];
                [mstopInvocation setArgument:&ipar atIndex:6];
                [mstopInvocation invoke];
                [mstopInvocation getReturnValue:&residual];
                break;
            default:
                [pcondrInvocation setArgument:&matrix atIndex:2];
                [pcondrInvocation setArgument:&t1v atIndex:3];
                [pcondrInvocation setArgument:&x atIndex:4];
                [pcondrInvocation setArgument:&ipar atIndex:5];
                [pcondrInvocation invoke];

                [matvecInvocation setArgument:&matrix atIndex:2];
                [matvecInvocation setArgument:&t1v atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    t1v[i] = r[i] - b[i];
                }
                [pcondlInvocation setArgument:&matrix atIndex:2];
                [pcondlInvocation setArgument:&r atIndex:3];
                [pcondlInvocation setArgument:&t1v atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, r, 1);
                break;
        }
        
        s[m] = cblas_dnrm2(HUTI_NDIM, r, 1);

        // Print debugging info if required
        if (HUTI_DBUGLVL != HUTI_NO_DEBUG) {
            if ( (iter_count % HUTI_DBUGLVL) == 0 ) {
                NSLog(@"FEMHUTIter:dgmresSolveInSolution: %d %11.4e\n", iter_count, residual);
            }
        }
        
        if (residual < HUTI_TOLERANCE) {
            HUTI_INFO = HUTI_CONVERGENCE;
            break;
        }
        
        // Return back to the top of the iteration loop (without initialization)
        iter_count = iter_count + 1;
        if (iter_count > HUTI_MAXIT) {
            HUTI_INFO = HUTI_MAXITER;
            break;
        }
        
    }
    
jump:
    
    // We have exited the loop after enough iterations or broke down
    if (HUTI_DBUGLVL != HUTI_NO_DEBUG) {
        NSLog(@"FEMHUTIter:dgmresSolveInSolution: %d %11.4e\n", iter_count, residual);
    }
    
    HUTI_ITERS = iter_count;
    
    // Release memory
    free_dvector(w, 0, ndim-1);
    free_dvector(r, 0, ndim-1);
    free_dvector(s, 0, ndim-1);
    free_dvector(t1v, 0, ndim-1);
    free_dvector(v, 0, ndim-1);
    free_dvector(buffer, 0, ndim-1);
    free_dvector(cs, 0, (HUTI_GMRES_RESTART+1)-1);
    free_dvector(sn, 0, (HUTI_GMRES_RESTART+1)-1);
    free_dvector(y, 0, (HUTI_GMRES_RESTART+1)-1);
    free_dmatrix(h, 0, (HUTI_GMRES_RESTART+1)-1, 0, (HUTI_GMRES_RESTART+1)-1);
}

@end
