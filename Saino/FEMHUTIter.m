//
//  HUTIter.m
//  Saino
//
//  Created by Hakime Seddik on 16/02/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMHUTIter.h"

#import "FEMKernel.h"

#import <Accelerate/Accelerate.h>

// Used in TFQRM solve. This is the magic ratio for upperb and tolerance used in upper bound convergence test
static double UPPERB_TOL_RATIO  =  10.0;

@interface FEMHUTIter ()

-(void)HUTI_dRandVector:(double *)u ipar:(int *)ipar;
-(void)HUTI_dLuSolveAt:(int)n luMatrix:(double **)lumat afterSolve:(double *)u rightHandSide:(double *)v;
-(void)HUTI_sgsNumberOfDimension:(int)n solution:(FEMSolution *)solution afterSolve:(double *)xvec rightHandSide:(double *)rhsvec ipar:(int *)ipar rounds:(double)rounds minTolerance:(double)minTolerance maxTolerance:(double)maxTolerance residual:(double *)residual converged:(BOOL *)converged diverged:(BOOL *)diverged outputInterval:(int)outputInterval omega:(double)omega matvecMethod:(SEL)matvecMethod;
-(void)HUTI_jacobiNumberOfDimension:(int)n solution:(FEMSolution *)solution afterSolve:(double *)xvec rightHandSide:(double *)rhsvec ipar:(int *)ipar rounds:(double)rounds minTolerance:(double)minTolerance maxTolerance:(double)maxTolerance residual:(double *)residual converged:(BOOL *)converged diverged:(BOOL *)diverged outputInterval:(int)outputInterval matvecMethod:(SEL)matvecMethod;
-(void)HUTI_BICGStabLNumberOfDimension:(int)n solution:(FEMSolution *)solution afterSolve:(double *)xvec rightHandSide:(double *)rhsvec ipar:(int *)ipar rounds:(double)rounds minTolerance:(double)minTolerance maxTolerance:(double)maxTolerance converged:(BOOL *)converged diverged:(BOOL *)diverged outputInterval:(int)outputInterval polynomialDegree:(int)polyDegree pcondlMethod:(SEL)pcondlMethod matvecMethod:(SEL)matvecMethod;
-(void)HUTI_gcrNumberOfDimension:(int)n solution:(FEMSolution *)solution afterSolve:(double *)xvec rightHandSide:(double *)rhsvec ipar:(int *)ipar rounds:(double)rounds minTolerance:(double)minTolerance maxTolerance:(double)maxTolerance residual:(double *)residual converged:(BOOL *)converged diverged:(BOOL *)diverged outputInterval:(int)outputInterval restart:(int)m pcondlMethod:(SEL)pcondlMethod matvecMethod:(SEL)matvecMethod;




@end

@implementation FEMHUTIter

#pragma mark private methods

-(void)HUTI_dRandVector:(double *)u ipar:(int *)ipar {
/**********************************************************
    This method fills a vector with pseudo random numbers
**********************************************************/
    int i;
    
    for (i=0; i<HUTI_NDIM; i++) {
        u[i] = rand() % 100;
    }
    
}

-(void)HUTI_dLuSolveAt:(int)n luMatrix:(double **)lumat afterSolve:(double *)u rightHandSide:(double *)v {
/*******************************************************************
    This method constructs LU decomposition of the given matrix 
    and solve LUu = v 
 ******************************************************************/   
    
    int i, j, k;
    
    // This is from Saad's book, algorithm 10.4
    
    for (i=1; i<=n; i++) {
        for (k=0; k<=i-1; k++) {
            
            // Check for small pivot
            if ( fabs(lumat[k][k]) < 1.0e-16 ) {
                NSLog(@"HUTI GMRES: small pivot %lf\n", lumat[k][k]);
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
        
        // Ccompute u(i) = u(i) / U(i,i)
        u[i] = u[i] / lumat[i][i];
    }
}

-(void)HUTI_sgsNumberOfDimension:(int)n solution:(FEMSolution *)solution afterSolve:(double *)xvec rightHandSide:(double *)rhsvec ipar:(int *)ipar rounds:(double)rounds minTolerance:(double)minTolerance maxTolerance:(double)maxTolerance residual:(double *)residual converged:(BOOL *)converged diverged:(BOOL *)diverged outputInterval:(int)outputInterval omega:(double)omega matvecMethod:(SEL)matvecMethod {
    
    int i, j, k;
    double *r;
    double bnorm, rnorm, s;
    matrixArraysContainer *matContainers;
    
    matContainers = solution.matrix.getContainers;
    
    NSMethodSignature *matvecSignature;
    NSInvocation *matvecInvocation;
    
    FEMPrecondition *preconditioning;
    
    // Acquire signature invocations and set selector for invocations
    matvecSignature = [FEMPrecondition instanceMethodSignatureForSelector:matvecMethod];
    matvecInvocation = [NSInvocation invocationWithMethodSignature:matvecSignature];
    [matvecInvocation setSelector:matvecMethod];
    
    // Instanciate class for the preconditioning
    preconditioning = [[FEMPrecondition alloc] init];
    
    // Targets to invocations
    [matvecInvocation setTarget:preconditioning];
    
    converged = NO;
    diverged = NO;
    
    r = doublevec(0, n-1);
    
    [matvecInvocation setArgument:&solution atIndex:2];
    [matvecInvocation setArgument:&xvec atIndex:3];
    [matvecInvocation setArgument:&r atIndex:4];
    [matvecInvocation setArgument:&ipar atIndex:5];
    [matvecInvocation invoke];
    
    for (i=0; i<n; i++) {
        r[i] = rhsvec[i] - r[i];
    }
    bnorm = cblas_dnrm2(n, rhsvec, 1);
    rnorm = cblas_dnrm2(n, r, 1);
    
    *residual = rnorm / bnorm;
    *converged = (*residual < minTolerance) ? YES : NO;
    *diverged = (*residual > maxTolerance) ? YES : NO;
    if (*converged == YES || *diverged == YES) {
            
        free_dvector(r, 0, n-1);

        return;
    }
    
    for (k=1; k<=rounds; k++) {
        for (i=0; i<n; i++) {
            s = 0.0;
            for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
                s = s + xvec[matContainers->Cols[j]] * matContainers->Values[j];
            }
            xvec[i] = xvec[i] + omega * (rhsvec[i]-s) / matContainers->Values[matContainers->Diag[i]];
        }
        
        for (i=n-1; i>=0; i--) {
            s = 0.0;
            for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
                s = s + xvec[matContainers->Cols[j]] * matContainers->Values[j];
            }
            xvec[i] = xvec[i] + omega * (rhsvec[i]-s) / matContainers->Values[matContainers->Diag[i]];
        }
        
        [matvecInvocation setArgument:&solution atIndex:2];
        [matvecInvocation setArgument:&xvec atIndex:3];
        [matvecInvocation setArgument:&r atIndex:4];
        [matvecInvocation setArgument:&ipar atIndex:5];
        [matvecInvocation invoke];
        
        for (i=0; i<n; i++) {
            r[i] = rhsvec[i] - r[i];
        }
        rnorm = cblas_dnrm2(n, r, 1);
        
        *residual = rnorm / bnorm;
        if ( (k % outputInterval) == 0) {
            NSLog(@"%d %lf %lf\n", k, rnorm, *residual);
        }
        
        *converged = (*residual < minTolerance) ? YES : NO;
        *diverged = (*residual > maxTolerance) ? YES : NO;
        if (*converged == YES || *diverged == YES) break;

    }
    
    free_dvector(r, 0, n-1);
    matContainers = NULL;
}

-(void)HUTI_jacobiNumberOfDimension:(int)n solution:(FEMSolution *)solution afterSolve:(double *)xvec rightHandSide:(double *)rhsvec ipar:(int *)ipar rounds:(double)rounds minTolerance:(double)minTolerance maxTolerance:(double)maxTolerance residual:(double *)residual converged:(BOOL *)converged diverged:(BOOL *)diverged outputInterval:(int)outputInterval matvecMethod:(SEL)matvecMethod {
    
    int i, k;
    double *r;
    double bnorm, rnorm;
    matrixArraysContainer *matContainers;
    
    matContainers = solution.matrix.getContainers;
    
    NSMethodSignature *matvecSignature;
    NSInvocation *matvecInvocation;
    
    FEMPrecondition *preconditioning;
    
    // Acquire signature invocations and set selector for invocations
    matvecSignature = [FEMPrecondition instanceMethodSignatureForSelector:matvecMethod];
    matvecInvocation = [NSInvocation invocationWithMethodSignature:matvecSignature];
    [matvecInvocation setSelector:matvecMethod];
    
    // Instanciate class for the preconditioning
    preconditioning = [[FEMPrecondition alloc] init];
    
    // Targets to invocations
    [matvecInvocation setTarget:preconditioning];
    
    converged = NO;
    diverged = NO;
    
    r = doublevec(0, n-1);
    
    [matvecInvocation setArgument:&solution atIndex:2];
    [matvecInvocation setArgument:&xvec atIndex:3];
    [matvecInvocation setArgument:&r atIndex:4];
    [matvecInvocation setArgument:&ipar atIndex:5];
    [matvecInvocation invoke];
    
    for (i=0; i<n; i++) {
        r[i] = rhsvec[i] - r[i];
    }
    bnorm = cblas_dnrm2(n, rhsvec, 1);
    rnorm = cblas_dnrm2(n, r, 1);

    *residual = rnorm / bnorm;
    *converged = (*residual < minTolerance) ? YES : NO;
    *diverged = (*residual > maxTolerance) ? YES : NO;
    if (*converged == YES || *diverged == YES) {
        
        free_dvector(r, 0, n-1);
    
        return;
    }
    
    for (k=1; k<=rounds; k++) {
        
        for (i=0; i<n; i++) {
            xvec[i] = xvec[i] + r[i] / matContainers->Values[matContainers->Diag[i]];
        }
        [matvecInvocation setArgument:&solution atIndex:2];
        [matvecInvocation setArgument:&xvec atIndex:3];
        [matvecInvocation setArgument:&r atIndex:4];
        [matvecInvocation setArgument:&ipar atIndex:5];
        [matvecInvocation invoke];
        
        for (i=0; i<n; i++) {
            r[i] = rhsvec[i] - r[i];
        }
        rnorm = cblas_dnrm2(n, r, 1);
        
        *residual = rnorm / bnorm;
        
        if (k % outputInterval == 0) {
            NSLog(@"%d %lf %lf\n", k, rnorm, *residual);
        }
        
        *converged = (*residual < minTolerance) ? YES : NO;
        *diverged = (*residual > maxTolerance) ? YES : NO;
        if (*converged == YES || *diverged == YES) break;
        
    }

    free_dvector(r, 0, n-1);
    matContainers = NULL;
    
}

-(void)HUTI_BICGStabLNumberOfDimension:(int)n solution:(FEMSolution *)solution afterSolve:(double *)xvec rightHandSide:(double *)rhsvec ipar:(int *)ipar rounds:(double)rounds minTolerance:(double)minTolerance maxTolerance:(double)maxTolerance converged:(BOOL *)converged diverged:(BOOL *)diverged outputInterval:(int)outputInterval polynomialDegree:(int)polyDegree pcondlMethod:(SEL)pcondlMethod matvecMethod:(SEL)matvecMethod {
/**********************************************************************************************************************************
 
 This method solves real linear system Ax = b by usinf the BICGStab(l) algorithm wirh l ≥ 2 and the right-oriented ILU(n) 
 preconditioning.
 
 The method has been written using as a starting point the work of D.R. Fokkema (subroutine zbistbl v1.1 1998). 
 
 Double precision version.
 
**********************************************************************************************************************************/
    
    double zero, one, *t, kappa0, kappal;
    double rnrm0, rnrm, mxnrmx, mxnrmr, errorind, delta = 1.0e-2, bnrm;
    int i, j, rr, r, u, xp, bp, z, zz, y0, yl, y, k, rows, cols, lda, nrhs, order, *iwork, round, info;
    double **work, **rwork, *rwork_transpose, alpha, beta, omega, rho0,rho1, sigma, varrho, hatgamma;
    double *buffer, *buffer2;
    char str;
    BOOL rcmp, xpdt, earlyExit;
    
    NSMethodSignature *pCondlSignature, *matvecSignature;
    NSInvocation *pcondlInvocation, *matvecInvocation;
    
    FEMPrecondition *preconditioning;
    
    if (polyDegree < 2) {
        errorfunct("HUTI_realBICGStabl", "Polynomial degree < 2.");
    }
    
    t = doublevec(0, n-1);
    work = doublematrix(0, n-1, 0, (3+2*(polyDegree+1))-1);
    rwork = doublematrix(0, (polyDegree+1)-1, 0, (3+2*(polyDegree+1))-1);
    rwork_transpose = doublevec(0, ( (polyDegree+1)*(3+2*(polyDegree+1)) )-1);
    iwork = intvec(0, (polyDegree-1)-1);
    
    buffer = doublevec(0, n-1);
    buffer2 = doublevec(0, n-1);
    
    // Acquire signature invocations and set selector for invocations
    pCondlSignature = [FEMPrecondition instanceMethodSignatureForSelector:pcondlMethod];
    pcondlInvocation = [NSInvocation invocationWithMethodSignature:pCondlSignature];
    [pcondlInvocation setSelector:pcondlMethod];
    
    matvecSignature = [FEMPrecondition instanceMethodSignatureForSelector:matvecMethod];
    matvecInvocation = [NSInvocation invocationWithMethodSignature:matvecSignature];
    [matvecInvocation setSelector:matvecMethod];
    
    // Instanciate class for the preconditioning
    preconditioning = [[FEMPrecondition alloc] init];
    
    // Targets to invocations
    [pcondlInvocation setTarget:preconditioning];
    [matvecInvocation setTarget:preconditioning];
    
    converged = NO;
    diverged = NO;

    if (all(xvec, '=', 0.0, n) == true) {
        for (i=0; i<n; i++) {
            xvec[i] = rhsvec[i];
        }
    }
    
    zero = 0.0;
    one = 1.0;
    
    for (i=0; i<n; i++) {
        for (j=0; j<3+2*(polyDegree+1); j++) {
            work[i][j] = 0.0;
        }
    }
    for (i=0; i<polyDegree+1; i++) {
        for (j=0; j<3+2*(polyDegree+1); j++) {
            rwork[i][j] = 0.0;
        }
    }
    
    rr = 0;
    r = rr + 1;
    u = r+(polyDegree+1);
    xp = u+(polyDegree+1);
    bp = xp+1;
    
    z = 0;
    zz = z+(polyDegree+1);
    y0 = zz+(polyDegree+1);
    yl = y0+1;
    y = yl+1;
    

    memset( buffer, 0.0, (n*sizeof(buffer)) );
    
    [matvecInvocation setArgument:&solution atIndex:2];
    [matvecInvocation setArgument:&xvec atIndex:3];
    [matvecInvocation setArgument:&buffer atIndex:4];
    [matvecInvocation setArgument:&ipar atIndex:5];
    [matvecInvocation invoke];
    // TODO: Add support for constraint matrix
    for (i=0; i<n; i++) {
        work[i][r] = buffer[i];
        work[i][r] = rhsvec[i] - work[i][r];
    }
    bnrm = cblas_dnrm2(n, rhsvec, 1);
    for (i=0; i<n; i++) {
        buffer[i] = work[i][r];
    }
    rnrm0 = cblas_dnrm2(n, buffer, 1);
    
    // Check whether the initial guess satisfies the stopping criterion
    errorind = rnrm0 / bnrm;
    *converged = (errorind < minTolerance) ? YES : NO;
    *diverged = (errorind > maxTolerance) ? YES : NO;
    
    if (*converged == YES || *diverged == YES) {
        
        free_dvector(t, 0, n-1);
        free_dmatrix(work, 0, n-1, 0, (3+2*(polyDegree+1))-1);
        free_dmatrix(rwork, 0, (polyDegree+1)-1, 0, (3+2*(polyDegree+1))-1);
        free_dvector(rwork_transpose, 0, ( (polyDegree+1)*(3+2*(polyDegree+1)) )-1);
        free_dvector(buffer, 0, n-1);
        free_dvector(buffer2, 0, n-1);
        
        return;
    }
    
    earlyExit = NO;
    
    for (i=0; i<n; i++) {
        work[i][rr] = work[i][r];
        work[i][bp] = work[i][r];
        work[i][xp] = xvec[i];
    }
    
    rnrm = rnrm0;
    mxnrmx = rnrm0;
    mxnrmr = rnrm0;
    alpha = zero;
    omega = one;
    sigma = one;
    rho0 = one;
    for (i=0; i<n; i++) {
        xvec[i] = zero;
    }
    
    free_dvector(buffer, 0, n-1);
    free_dvector(buffer2, 0, n-1);
    
    buffer = doublevec(0, (polyDegree+1)-1);
    buffer2 = doublevec(0, (polyDegree+1)-1);
    
    for (round = 1; round <= rounds; round++) {
        
        // The BiCG part
        
        rho0 = -omega*rho0;
        
        for (k=1; k<=polyDegree; k++) {
            
            for (i=0; i<n; i++) {
                buffer[i] = work[i][rr];
                buffer2[i] = work[i][r+k-1];
            }
            rho1 = cblas_ddot(n, buffer, 1, buffer2, 1);
            if (rho0 == zero) {
                errorfunct("HUTI_realBICGStabl", "Breakdown error.");
            }
            
            beta = alpha * (rho1/rho0);
            rho0 = rho1;
            for (j=0; j<=k-1; j++) {
                for (i=0; i<n; i++) {
                    work[i][u+j] = work[i][r+j] - beta * work[i][u+j];
                }
            }
            for (i=0; i<n; i++) {
                buffer[i] = work[i][u+k-1];
            }
            [pcondlInvocation setArgument:&solution atIndex:2];
            [pcondlInvocation setArgument:&t atIndex:3];
            [pcondlInvocation setArgument:&buffer atIndex:4];
            [pcondlInvocation setArgument:&ipar atIndex:5];
            [pcondlInvocation invoke];

            [matvecInvocation setArgument:&solution atIndex:2];
            [matvecInvocation setArgument:&t atIndex:3];
            [matvecInvocation setArgument:&buffer atIndex:4];
            [matvecInvocation setArgument:&ipar atIndex:5];
            [matvecInvocation invoke];
            for (i=0; i<n; i++) {
                work[i][u+k] = buffer[i];
            }
            
            for (i=0; i<n; i++) {
                buffer[i] = work[i][rr];
                buffer2[i] = work[i][u+k];
            }
            sigma = cblas_ddot(n, buffer, 1, buffer2, 1);
            if (sigma == zero) {
                errorfunct("HUTI_realBICGStabl", "Breakdown error.");
            }
            
            alpha = rho1/sigma;
            for (i=0; i<n; i++) {
                xvec[i] = xvec[i] + alpha * work[i][u];
            }
            for (j=0; j<=k-1; j++) {
                for (i=0; i<n; i++) {
                    work[i][r+j] = work[i][r+j] - alpha * work[i][u+j+1];
                }
            }
            for (i=0; i<n; i++) {
                buffer[i] = work[i][r+k-1];
            }
            [pcondlInvocation setArgument:&solution atIndex:2];
            [pcondlInvocation setArgument:&t atIndex:3];
            [pcondlInvocation setArgument:&buffer atIndex:4];
            [pcondlInvocation setArgument:&ipar atIndex:5];
            [pcondlInvocation invoke];
            
            [matvecInvocation setArgument:&solution atIndex:2];
            [matvecInvocation setArgument:&t atIndex:3];
            [matvecInvocation setArgument:&buffer atIndex:4];
            [matvecInvocation setArgument:&ipar atIndex:5];
            [matvecInvocation invoke];
            for (i=0; i<n; i++) {
                work[i][r+k] = buffer[i];
            }
            
            for (i=0; i<n; i++) {
                buffer[i] = work[i][r];
            }
            rnrm = cblas_dnrm2(n, buffer, 1);
            mxnrmx = max(mxnrmx, rnrm);
            mxnrmr = max(mxnrmr, rnrm);
            
            //  In some simple cases, a few BiCG updates may already be enough to
            //  obtain the solution. The following is for handling this special case.
            errorind = rnrm / bnrm;
            *converged = (errorind < minTolerance) ? YES : NO;
            if (*converged == YES) {
                earlyExit = YES;
                break;
            }
        }
        
        if (earlyExit == YES) break;
        
        // The convex polynomial part
        for (i=1; i<=polyDegree+1; i++) {
            for (j=1; j<=i; j++) {
                for (k=0; k<n; k++) {
                    buffer[k] = work[k][r+i-1];
                    buffer2[k] = work[k][r+j-1];
                }
                rwork[i-1][j-1] = cblas_ddot(n, buffer, 1, buffer2, 1);
                
            }
        }
        for (j=1; j<polyDegree+1; j++) {
            for (i=0; i<=j-1; i++) {
                rwork[i][j] = rwork[j][i];
            }
        }
        
        k = z;
        for (i=0; i<polyDegree+1; i++) {
            for (j=zz; j<=zz+polyDegree; j++) {
                rwork[i][zz] = rwork[i][k];
                k++;
            }
        }
        
        // Transfrom rwork for LAPACK, column-major order
        for (i=0; i<polyDegree+1; i++) {
            for (j=0; j<(3+2*(polyDegree+1)); j++) {
                rwork_transpose[j+(polyDegree+1)*i] = rwork[j][i];
            }
        }
        
        rows = polyDegree-1;
        cols = polyDegree-1;
        lda = polyDegree+1;
        dgetrf_(&rows, &cols, rwork_transpose+(((zz+2)*lda)-polyDegree), &lda, iwork, &info);
        
        // tild r0 and tild rl (small vectors)
        
        rwork_transpose[((y0+1)*lda)-(polyDegree+1)] = -one;
        for (i=2; i<=polyDegree; i++) {
            rwork_transpose[((y0+1)*lda)-i] = rwork_transpose[((z+1)*lda)-i];
        }
        
        order = polyDegree-1;
        nrhs = 1;
        str = 'N';
        lda = polyDegree+1;
        dgetrs_(&str, &order, &nrhs, rwork_transpose+(((zz+2)*lda)-polyDegree), &lda, iwork, rwork_transpose+(((y0+1)*lda)-polyDegree), &lda, &info);
        rwork[(polyDegree+1)-1][y0] = zero;
        
        rwork_transpose[((yl+1)*lda)-(polyDegree+1)] = zero;
        for (i=2; i<=polyDegree; i++) {
            rwork_transpose[((yl+1)*lda)-i] = rwork_transpose[((z+1+polyDegree)*lda)-i];
        }
        dgetrs_(&str, &order, &nrhs, rwork_transpose+(((zz+2)*lda)-polyDegree), &lda, iwork, rwork_transpose+(((yl+1)*lda)-polyDegree), &lda, &info);
        rwork_transpose[((yl+1)*lda)-1] = -one;
        
        // Back to the original rwork matrix
        for (i=0; i<polyDegree+1; i++) {
            for (j=0; j<(3+2*(polyDegree+1)); j++) {
                rwork[j][i] = rwork_transpose[j+(polyDegree+1)*i];
            }
        }

        // Convex combination
        
        memset( buffer, 0.0, ((polyDegree+1)*sizeof(buffer)) );
        memset( buffer2, 0.0, ((polyDegree+1)*sizeof(buffer2)) );
        for (i=0; i<polyDegree+1; i++) {
             buffer[i] = rwork[i][y0];
        }
        cblas_dsymv(CblasRowMajor, CblasUpper, polyDegree+1, one, *rwork, (3+2*(polyDegree+1)), buffer, 1, zero, buffer2, 1);
        kappa0 = sqrt( cblas_ddot(polyDegree+1, buffer, 1, buffer2, 1) );
        
        memset( buffer, 0.0, ((polyDegree+1)*sizeof(buffer)) );
        for (i=0; i<polyDegree+1; i++) {
            buffer[i] = rwork[i][yl];
        }
        cblas_dsymv(CblasRowMajor, CblasUpper, polyDegree+1, one, *rwork, (3+2*(polyDegree+1)), buffer, 1, zero, buffer2, 1);
        kappal = sqrt( cblas_ddot(polyDegree+1, buffer, 1, buffer2, 1) );
        
        memset( buffer, 0.0, ((polyDegree+1)*sizeof(buffer)) );
        for (i=0; i<polyDegree+1; i++) {
            buffer[i] = rwork[i][y0];
        }
        cblas_dsymv(CblasRowMajor, CblasUpper, polyDegree+1, one, *rwork, (3+2*(polyDegree+1)), buffer, 1, zero, buffer2, 1);
        memset( buffer, 0.0, ((polyDegree+1)*sizeof(buffer)) );
        for (i=0; i<polyDegree+1; i++) {
            buffer[i] = rwork[i][yl];
        }
        varrho = cblas_ddot(polyDegree+1, buffer, 1, buffer2, 1) / (kappa0*kappal);
        hatgamma = varrho / fabs(varrho) * max(fabs(varrho), 7e-1) * kappa0/kappal;
        for (i=0; i<polyDegree+1; i++) {
            rwork[i][y0] = rwork[i][y0] - hatgamma * rwork[i][yl];
        }
        
        // Update
    
        omega = rwork[(polyDegree+1)-1][y0];
        for (j=1; j<=polyDegree; j++) {
            for (i=0; i<n; i++) {
                work[i][u] = work[i][u] - rwork[(j+1)-1][y0] * work[i][u+j];
                xvec[i] = xvec[i] + rwork[(j+1)-1][y0] * work[i][r+j-1];
                work[i][r] = work[i][r] - rwork[(j+1)-1][y0] * work[i][r+j];
            }
        }
        
        memset( buffer, 0.0, ((polyDegree+1)*sizeof(buffer)) );
        for (i=0; i<polyDegree+1; i++) {
            buffer[i] = rwork[i][y0];
        }
        cblas_dsymv(CblasRowMajor, CblasUpper, polyDegree+1, one, *rwork, (3+2*(polyDegree+1)), buffer, 1, zero, buffer2, 1);
        rnrm = sqrt( cblas_ddot(polyDegree+1, buffer, 1, buffer2, 1) );
        
        // The reliable update part
        
        mxnrmx = max(mxnrmx, rnrm);
        mxnrmr = max(mxnrmr, rnrm);
        xpdt = (rnrm < delta*rnrm0 && rnrm0 < mxnrmx) ? YES: NO;
        rcmp = ((rnrm < delta*mxnrmr && rnrm0 < mxnrmr) || xpdt == YES) ? YES : NO;
        if (rcmp == YES) {
            [pcondlInvocation setArgument:&solution atIndex:2];
            [pcondlInvocation setArgument:&t atIndex:3];
            [pcondlInvocation setArgument:&xvec atIndex:4];
            [pcondlInvocation setArgument:&ipar atIndex:5];
            [pcondlInvocation invoke];
            
            memset( buffer, 0.0, (n*sizeof(buffer)) );
            [matvecInvocation setArgument:&solution atIndex:2];
            [matvecInvocation setArgument:&t atIndex:3];
            [matvecInvocation setArgument:&buffer atIndex:4];
            [matvecInvocation setArgument:&ipar atIndex:5];
            [matvecInvocation invoke];
            for (i=0; i<n; i++) {
                 work[i][r] = buffer[i];
            }
            for (i=0; i<n; i++) {
                work[i][r] = work[i][bp] - work[i][r];
            }
            mxnrmr = rnrm;
            if (xpdt == YES) {
                for (i=0; i<n; i++) {
                    work[i][xp] = work[i][xp] + t[i];
                    xvec[i] = zero;
                    work[i][bp] = work[i][r];
                    mxnrmr = rnrm;
                }
            }
        }
        
        if (rcmp == YES) {
            if (xpdt == YES) {
                for (i=0; i<n; i++) {
                    t[i] = work[i][xp];
                }
            } else {
                for (i=0; i<n; i++) {
                    t[i] = t[i] + work[i][xp];
                }
            }
        } else {
            [pcondlInvocation setArgument:&solution atIndex:2];
            [pcondlInvocation setArgument:&t atIndex:3];
            [pcondlInvocation setArgument:&xvec atIndex:4];
            [pcondlInvocation setArgument:&ipar atIndex:5];
            [pcondlInvocation invoke];
            for (i=0; i<n; i++) {
                t[i] = t[i] + work[i][xp];
            }
        }
        
        errorind = rnrm / bnrm;
        if (round % outputInterval == 0) {
            NSLog(@"%d %lf %lf\n", round, rnrm, errorind);
        }
        
        *converged = (errorind < minTolerance) ? YES: NO;
        *diverged = (errorind > maxTolerance) ? YES : NO;
        if (*converged == YES || *diverged == YES) break;
    
    } // end of rounds
    
    if (outputInterval != HUGE_VAL) {
        NSLog(@"%d %lf %lf\n", round, rnrm, errorind);
    }
    
    // We have solved z = P*x, with P the preconditioner, so finally 
    // solve the true unknown x
    for (i=0; i<n; i++) {
        t[i] = xvec[i];
    }
    [pcondlInvocation setArgument:&solution atIndex:2];
    [pcondlInvocation setArgument:&xvec atIndex:3];
    [pcondlInvocation setArgument:&t atIndex:4];
    [pcondlInvocation setArgument:&ipar atIndex:5];
    [pcondlInvocation invoke];
    for (i=0; i<n; i++) {
        xvec[i] = xvec[i] + work[i][xp];
    }
    
    free_dvector(t, 0, n-1);
    free_dmatrix(work, 0, n-1, 0, (3+2*(polyDegree+1))-1);
    free_dmatrix(rwork, 0, (polyDegree+1)-1, 0, (3+2*(polyDegree+1))-1);
    free_dvector(rwork_transpose, 0, ( (polyDegree+1)*(3+2*(polyDegree+1)) )-1);
    free_dvector(buffer, 0, (polyDegree+1)-1);
    free_dvector(buffer2, 0, (polyDegree+1)-1);
}

-(void)HUTI_gcrNumberOfDimension:(int)n solution:(FEMSolution *)solution afterSolve:(double *)xvec rightHandSide:(double *)rhsvec ipar:(int *)ipar rounds:(double)rounds minTolerance:(double)minTolerance maxTolerance:(double)maxTolerance residual:(double *)residual converged:(BOOL *)converged diverged:(BOOL *)diverged outputInterval:(int)outputInterval restart:(int)m pcondlMethod:(SEL)pcondlMethod matvecMethod:(SEL)matvecMethod {
    
    int i, j, k, l;
    double bnorm, rnorm;
    double alpha, beta;
    double *r, **s, **v, *t1, *t2;
    double *buffer;
    
    NSMethodSignature *pCondlSignature, *matvecSignature;
    NSInvocation *pcondlInvocation, *matvecInvocation;
    
    FEMPrecondition *preconditioning;
    
    converged = NO;
    diverged = NO;
    
    r = doublevec(0, n-1);
    t1 = doublevec(0, n-1);
    t2 = doublevec(0, n-1);
    
    if (m > 1) {
        v = doublematrix(0, n-1, 0, (m-1)-1);
        s = doublematrix(0, n-1, 0, (m-1)-1);
        for (i=0; i<n; i++) {
            for (j=0; j<m-1; j++) {
                v[i][j] = 0.0;
                s[i][j] = 0.0;
            }
        }
    }
    
    // Acquire signature invocations and set selector for invocations
    pCondlSignature = [FEMPrecondition instanceMethodSignatureForSelector:pcondlMethod];
    pcondlInvocation = [NSInvocation invocationWithMethodSignature:pCondlSignature];
    [pcondlInvocation setSelector:pcondlMethod];
    
    matvecSignature = [FEMPrecondition instanceMethodSignatureForSelector:matvecMethod];
    matvecInvocation = [NSInvocation invocationWithMethodSignature:matvecSignature];
    [matvecInvocation setSelector:matvecMethod];
    
    // Instanciate class for the preconditioning
    preconditioning = [[FEMPrecondition alloc] init];
    
    // Targets to invocations
    [pcondlInvocation setTarget:preconditioning];
    [matvecInvocation setTarget:preconditioning];

    
    memcpy(xvec, rhsvec, n*sizeof(rhsvec));

    [matvecInvocation setArgument:&solution atIndex:2];
    [matvecInvocation setArgument:&xvec atIndex:3];
    [matvecInvocation setArgument:&r atIndex:4];
    [matvecInvocation setArgument:&ipar atIndex:5];
    [matvecInvocation invoke];
    
    for (i=0; i<n; i++) {
        r[i] = rhsvec[i] - r[i];
    }
    
    bnorm = cblas_dnrm2(n, rhsvec, 1);
    rnorm = cblas_dnrm2(n, r, 1);
    
    *residual = rnorm / bnorm;
    *converged = (*residual < minTolerance) ? YES: NO;
    *diverged = (*residual > maxTolerance) ? YES: NO;
    if (*converged == YES || *diverged == YES) {
        
        free_dvector(r, 0, n-1);
        free_dvector(t1, 0, n-1);
        free_dvector(t2, 0, n-1);
        if (m > 1) {
            free_dmatrix(v, 0, n-1, 0, (m-1)-1);
            free_dmatrix(s, 0, n-1, 0, (m-1)-1);
        }
        
        return;
    }
    
    buffer = doublevec(0, n-1);
    
    for (k=1; k<=rounds; k++) {
        
        // Check for restarting
        if ( k % m == 0) {
            j = m;
        } else {
            j = k % m;
        }
        
        // Perform the preconditioning
        
        [pcondlInvocation setArgument:&solution atIndex:2];
        [pcondlInvocation setArgument:&t1 atIndex:3];
        [pcondlInvocation setArgument:&r atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];
        
        [matvecInvocation setArgument:&solution atIndex:2];
        [matvecInvocation setArgument:&t1 atIndex:3];
        [matvecInvocation setArgument:&t2 atIndex:4];
        [matvecInvocation setArgument:&ipar atIndex:5];
        [matvecInvocation invoke];

        // Perform the orthogonalization of the search directions
        
        for (i=0; i<j-1; i++) {
            for (l=0; l<n; l++) {
                buffer[l] = v[l][i];
            }
            beta = cblas_ddot(n, buffer, 1, t2, 1);
            for (l=0; l<n; l++) {
                t1[l] = t1[l] - beta * s[l][i];
                t2[l] = t2[l] - beta * v[l][i];
            }
        }
        
        alpha = cblas_dnrm2(n, t2, 1);
        for (i=0; i<n; i++) {
            t1[i] = 1.0 / alpha * t1[i];
            t2[i] = 1.0 / alpha *t2[i];
        }
        
        // The update of the solution and savethe search data
        
        beta = cblas_ddot(n, t2, 1, r, 1);
        for (i=0; i<n; i++) {
            xvec[i] = xvec[i] + beta * t1[i];
            r[i] = r[i] - beta * t2[i];
        }
        if (j != m) {
            for (i=0; i<n; i++) {
                s[i][j-1] = t1[i];
                v[i][j-1] = t2[i];
            }
        }
        
        // Check whether the convergence criterion is fulfilled
        
        rnorm = cblas_dnrm2(n, r, 1);
        *residual = rnorm / bnorm;
        
        if (k % outputInterval == 0) {
            NSLog(@"%d %lf\n", k, *residual);
        }
        
        *converged = (*residual < minTolerance) ? YES: NO;
        *diverged = (*residual > maxTolerance) ? YES: NO;
        if (*converged == YES || *diverged == YES) break;
        
    }
    
    free_dvector(r, 0, n-1);
    free_dvector(t1, 0, n-1);
    free_dvector(t2, 0, n-1);
    free_dvector(buffer, 0, n-1);
    if (m > 1) {
        free_dmatrix(v, 0, n-1, 0, (m-1)-1);
        free_dmatrix(s, 0, n-1, 0, (m-1)-1);
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
            
            NSLog(@"Environment variable NUMBER_OF_PROCESSORS has an illegal value: %s\n", evname);
        }
    
    huti_init_done = YES;
}

-(void)hutiExit {
    
}

#pragma mark BI-CGSTAB

-(void)dbicgstabSolveInSolution:(FEMSolution *)solution ndim:(int)ndim wrkdim:(int)wrkdim ipar:(int *)ipar dpar:(double *)dpar work:(double **)work pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod {
/**********************************************************************************************************************************
 
    This method is based on Barret et al. book: "Templates for the Solution of Linear Systems: Building blocks for 
    Iterative Methods", 1993.
 
    All matrix-vector operations are done externally, so we don't need to know about the matrix structure (sparse or dense).
    Memory allocationf for the work array has also been done externally.
    
    Double precision version.
 
**********************************************************************************************************************************/
    
    int i, iter_count;
    double rho, oldrho, alpha, beta, omega;
    double *rtld, *p, *t1v, *v, *s, *t2v, *t, *r;
    double residual, rhsnorm, precrhsnorm;
    matrixArraysContainer *matContainers;
    variableArraysContainer *varContainers;
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
    
    NSMethodSignature *pCondlSignature, *pCondrSignature, *matvecSignature, *mstopSignature;
    NSInvocation *pcondlInvocation, *pcondrInvocation, *matvecInvocation, *mstopInvocation;
    
    FEMPrecondition *preconditioning;
    FEMKernel *kernel;
    
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
    
    mstopSignature = [FEMKernel instanceMethodSignatureForSelector:mstopMethod];
    mstopInvocation = [NSInvocation invocationWithMethodSignature:mstopSignature];
    [mstopInvocation setSelector:mstopMethod];
    
    // Instanciate class for the preconditioning
    preconditioning = [[FEMPrecondition alloc] init];
    
    //Instanciate class for the kernel
    kernel = [[FEMKernel alloc] init];
    
    // Targets to invocations
    [pcondlInvocation setTarget:preconditioning];
    [pcondrInvocation setTarget:preconditioning];
    [matvecInvocation setTarget:preconditioning];
    [mstopInvocation setTarget:kernel];
    
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
        rhsnorm = cblas_dnrm2(HUTI_NDIM, matContainers->RHS, 1);
    }
    if (HUTI_STOPC == HUTI_PRESID_SCALED_BYPRECB) {
        [pcondlInvocation setArgument:&solution atIndex:2];
        [pcondlInvocation setArgument:&p atIndex:3];
        [pcondlInvocation setArgument:&matContainers->RHS atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];
        
        precrhsnorm = cblas_dnrm2(HUTI_NDIM, p, 1);
    }
    
    // Generate vector xvec if needed
    if (HUTI_INITIALX == HUTI_RANDOMX) {
        [self HUTI_dRandVector:varContainers->Values ipar:ipar];
    } else if (HUTI_INITIALX != HUTI_USERSUPPLIEDX) {
        for (i=0; i<varContainers->sizeValues; i++) {
            varContainers->Values[i] = 1;
        }
    }
    
    [matvecInvocation setArgument:&solution atIndex:2];
    [matvecInvocation setArgument:&varContainers->Values atIndex:3];
    [matvecInvocation setArgument:&r atIndex:4];
    [matvecInvocation setArgument:&ipar atIndex:5];
    [matvecInvocation invoke];
    
    for (i=0; i<ndim; i++) {
        r[i] = matContainers->RHS[i] - r[i];
        rtld[i] = r[i];
    }
    
    memset( p, 0.0, (ndim*sizeof(p)) );
    memset( v, 0.0, (ndim*sizeof(v)) );
    
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
        
        [pcondlInvocation setArgument:&solution atIndex:2];
        [pcondlInvocation setArgument:&v atIndex:3];
        [pcondlInvocation setArgument:&p atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];
        
        [pcondrInvocation setArgument:&solution atIndex:2];
        [pcondrInvocation setArgument:&t1v atIndex:3];
        [pcondrInvocation setArgument:&v atIndex:4];
        [pcondrInvocation setArgument:&ipar atIndex:5];
        [pcondrInvocation invoke];
        
        [matvecInvocation setArgument:&solution atIndex:2];
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
                varContainers->Values[i] = varContainers->Values[i] + alpha * t1v[i];
                HUTI_INFO = HUTI_BICGSTAB_SNORM;
                break;
            }
        }
        
        [pcondlInvocation setArgument:&solution atIndex:2];
        [pcondlInvocation setArgument:&t atIndex:3];
        [pcondlInvocation setArgument:&s atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];
        
        [pcondrInvocation setArgument:&solution atIndex:2];
        [pcondrInvocation setArgument:&t2v atIndex:3];
        [pcondrInvocation setArgument:&t atIndex:4];
        [pcondrInvocation setArgument:&ipar atIndex:5];
        [pcondrInvocation invoke];
        
        [matvecInvocation setArgument:&solution atIndex:2];
        [matvecInvocation setArgument:&t2v atIndex:3];
        [matvecInvocation setArgument:&t atIndex:4];
        [matvecInvocation setArgument:&ipar atIndex:5];
        [matvecInvocation invoke];
        
        omega = ( cblas_ddot(HUTI_NDIM, t, 1, s, 1) ) / ( cblas_ddot(HUTI_NDIM, t, 1, t, 1));
        for (i=0; i<ndim; i++) {
            varContainers->Values[i] = varContainers->Values[i] + alpha * t1v[i] + omega * t2v[i];
            r[i] = s[i] - omega * t[i];
        }
        
        // Check the convergence against selected stopping criterion
        switch (HUTI_STOPC) {
            case HUTI_TRUERESIDUAL:
                [matvecInvocation setArgument:&solution atIndex:2];
                [matvecInvocation setArgument:&varContainers->Values atIndex:3];
                [matvecInvocation setArgument:&t2v atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    t1v[i] = t2v[i] - matContainers->RHS[i];
                }
                residual = cblas_dnrm2(HUTI_NDIM, t1v, 1);
                break;
            case HUTI_TRESID_SCALED_BYB:
                [matvecInvocation setArgument:&solution atIndex:2];
                [matvecInvocation setArgument:&varContainers->Values atIndex:3];
                [matvecInvocation setArgument:&t2v atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    t1v[i] = t2v[i] - matContainers->RHS[i];
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
                [mstopInvocation setArgument:&solution atIndex:2];
                [mstopInvocation setArgument:&varContainers->Values atIndex:3];
                [mstopInvocation setArgument:&matContainers->RHS atIndex:4];
                [mstopInvocation setArgument:&r atIndex:5];
                [mstopInvocation setArgument:&ipar atIndex:6];
                [mstopInvocation invoke];
                [mstopInvocation getReturnValue:&residual];
                break;
            default:
                [matvecInvocation setArgument:&solution atIndex:2];
                [matvecInvocation setArgument:&varContainers->Values atIndex:3];
                [matvecInvocation setArgument:&t2v atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    t1v[i] = t2v[i] - matContainers->RHS[i];
                }
                residual = cblas_dnrm2(HUTI_NDIM, t1v, 1);
                break;
        }
        
        // Print debugging info if required
        if (HUTI_DBUGLVL != HUTI_NO_DEBUG) {
            if ( (iter_count % HUTI_DBUGLVL) == 0 ) {
                NSLog(@"%d %lf\n", iter_count, residual);
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
        NSLog(@"%d %lf\n", iter_count, residual);
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
    
    matContainers = NULL;
    varContainers = NULL;
}

#pragma mark BI-CGSTAB(2)

-(void)dbicgstab2SolveInSolution:(FEMSolution *)solution ndim:(int)ndim wrkdim:(int)wrkdim ipar:(int *)ipar dpar:(double *)dpar work:(double **)work pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod {
/**********************************************************************************************************************************
 
    This method is based on the Henk A. Van der Vorst paper: "Parallel Iretative Solution Methods for the Linear Systems
    arising from Discretized PDE's". This is the Bi-CGSTAB(2) version.
 
    All matrix-vector operations are done externally, so we don't need to know about the matrix structure (sparse or dense).
    Memory allocationf for the work array has also been done externally.
 
    Double precision version.
 
**********************************************************************************************************************************/
    
    int i, iter_count;
    double rho, oldrho, alpha, beta, omega1, omega2;
    double tau, delta, myy;
    double *rtld, *u, *t1v, *v, *s, *w, *t, *r;
    double residual, rhsnorm, precrhsnorm;
    matrixArraysContainer *matContainers;
    variableArraysContainer *varContainers;
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
    
    NSMethodSignature *pCondlSignature, *pCondrSignature, *matvecSignature, *mstopSignature;
    NSInvocation *pcondlInvocation, *pcondrInvocation, *matvecInvocation, *mstopInvocation;
    
    FEMPrecondition *preconditioning;
    FEMKernel *kernel;
    
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
    
    mstopSignature = [FEMKernel instanceMethodSignatureForSelector:mstopMethod];
    mstopInvocation = [NSInvocation invocationWithMethodSignature:mstopSignature];
    [mstopInvocation setSelector:mstopMethod];
    
    // Instanciate class for the preconditioning
    preconditioning = [[FEMPrecondition alloc] init];
    
    //Instanciate class for the kernel
    kernel = [[FEMKernel alloc] init];
    
    // Targets to invocations
    [pcondlInvocation setTarget:preconditioning];
    [pcondrInvocation setTarget:preconditioning];
    [matvecInvocation setTarget:preconditioning];
    [mstopInvocation setTarget:kernel];
    
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
        rhsnorm = cblas_dnrm2(HUTI_NDIM, matContainers->RHS, 1);
    }
    if (HUTI_STOPC == HUTI_PRESID_SCALED_BYPRECB) {
        [pcondlInvocation setArgument:&solution atIndex:2];
        [pcondlInvocation setArgument:&t1v atIndex:3];
        [pcondlInvocation setArgument:&matContainers->RHS atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];
        
        precrhsnorm = cblas_dnrm2(HUTI_NDIM, t1v, 1);
    }
    
    // Generate vector xvec if needed
    if (HUTI_INITIALX == HUTI_RANDOMX) {
        [self HUTI_dRandVector:varContainers->Values ipar:ipar];
    } else if (HUTI_INITIALX != HUTI_USERSUPPLIEDX) {
        for (i=0; i<varContainers->sizeValues; i++) {
            varContainers->Values[i] = 1;
        }
    }
    
    [pcondrInvocation setArgument:&solution atIndex:2];
    [pcondrInvocation setArgument:&u atIndex:3];
    [pcondrInvocation setArgument:&varContainers->Values atIndex:4];
    [pcondrInvocation setArgument:&ipar atIndex:5];
    [pcondrInvocation invoke];
    
    [matvecInvocation setArgument:&solution atIndex:2];
    [matvecInvocation setArgument:&u atIndex:3];
    [matvecInvocation setArgument:&r atIndex:4];
    [matvecInvocation setArgument:&ipar atIndex:5];
    [matvecInvocation invoke];
    
    for (i=0; i<ndim; i++) {
        u[i] = matContainers->RHS[i] - r[i];
    }
    
    [pcondlInvocation setArgument:&solution atIndex:2];
    [pcondlInvocation setArgument:&r atIndex:3];
    [pcondlInvocation setArgument:&u atIndex:4];
    [pcondlInvocation setArgument:&ipar atIndex:5];
    [pcondlInvocation invoke];
    
    for (i=0; i<ndim; i++) {
        rtld[i] = r[i];
    }
    memset( u, 0.0, (ndim*sizeof(u)) );
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
        
        [pcondrInvocation setArgument:&solution atIndex:2];
        [pcondrInvocation setArgument:&v atIndex:3];
        [pcondrInvocation setArgument:&u atIndex:4];
        [pcondrInvocation setArgument:&ipar atIndex:5];
        [pcondrInvocation invoke];
        
        [matvecInvocation setArgument:&solution atIndex:2];
        [matvecInvocation setArgument:&v atIndex:3];
        [matvecInvocation setArgument:&t1v atIndex:4];
        [matvecInvocation setArgument:&ipar atIndex:5];
        [matvecInvocation invoke];
        
        [pcondlInvocation setArgument:&solution atIndex:2];
        [pcondlInvocation setArgument:&v atIndex:3];
        [pcondlInvocation setArgument:&t1v atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];

        alpha = oldrho / cblas_ddot(HUTI_NDIM, rtld, 1, v, 1);
        for (i=0; i<ndim; i++) {
            r[i] = r[i] - alpha * v[i];
        }
        
        [pcondrInvocation setArgument:&solution atIndex:2];
        [pcondrInvocation setArgument:&s atIndex:3];
        [pcondrInvocation setArgument:&r atIndex:4];
        [pcondrInvocation setArgument:&ipar atIndex:5];
        [pcondrInvocation invoke];
        
        [matvecInvocation setArgument:&solution atIndex:2];
        [matvecInvocation setArgument:&s atIndex:3];
        [matvecInvocation setArgument:&t1v atIndex:4];
        [matvecInvocation setArgument:&ipar atIndex:5];
        [matvecInvocation invoke];
        
        [pcondlInvocation setArgument:&solution atIndex:2];
        [pcondlInvocation setArgument:&s atIndex:3];
        [pcondlInvocation setArgument:&t1v atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];
        
        for (i=0; i<ndim; i++) {
            varContainers->Values[i] = varContainers->Values[i] + alpha * u[i];
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
        
        [pcondrInvocation setArgument:&solution atIndex:2];
        [pcondrInvocation setArgument:&w atIndex:3];
        [pcondrInvocation setArgument:&v atIndex:4];
        [pcondrInvocation setArgument:&ipar atIndex:5];
        [pcondrInvocation invoke];
        
        [matvecInvocation setArgument:&solution atIndex:2];
        [matvecInvocation setArgument:&w atIndex:3];
        [matvecInvocation setArgument:&t1v atIndex:4];
        [matvecInvocation setArgument:&ipar atIndex:5];
        [matvecInvocation invoke];
        
        [pcondlInvocation setArgument:&solution atIndex:2];
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
        
        [pcondrInvocation setArgument:&solution atIndex:2];
        [pcondrInvocation setArgument:&t atIndex:3];
        [pcondrInvocation setArgument:&s atIndex:4];
        [pcondrInvocation setArgument:&ipar atIndex:5];
        [pcondrInvocation invoke];
        
        [matvecInvocation setArgument:&solution atIndex:2];
        [matvecInvocation setArgument:&t atIndex:3];
        [matvecInvocation setArgument:&t1v atIndex:4];
        [matvecInvocation setArgument:&ipar atIndex:5];
        [matvecInvocation invoke];
        
        [pcondlInvocation setArgument:&solution atIndex:2];
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
            varContainers->Values[i] = varContainers->Values[i] + omega1 * r[i] + omega2 * s[i] + alpha * u[i];
            r[i] = r[i] - omega1 * s[i] - omega2 * t[i];
        }
        
        // Check the convergence against selected stopping criterion
        switch (HUTI_STOPC) {
            case HUTI_TRUERESIDUAL:
                [pcondrInvocation setArgument:&solution atIndex:2];
                [pcondrInvocation setArgument:&s atIndex:3];
                [pcondrInvocation setArgument:&varContainers->Values atIndex:4];
                [pcondrInvocation setArgument:&ipar atIndex:5];
                [pcondrInvocation invoke];
                
                [matvecInvocation setArgument:&solution atIndex:2];
                [matvecInvocation setArgument:&s atIndex:3];
                [matvecInvocation setArgument:&t1v atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    t1v[i] = t1v[i] - matContainers->RHS[i];
                }
                [pcondlInvocation setArgument:&solution atIndex:2];
                [pcondlInvocation setArgument:&s atIndex:3];
                [pcondlInvocation setArgument:&t1v atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, s, 1);
                break;
            case HUTI_TRESID_SCALED_BYB:
                [pcondrInvocation setArgument:&solution atIndex:2];
                [pcondrInvocation setArgument:&s atIndex:3];
                [pcondrInvocation setArgument:&varContainers->Values atIndex:4];
                [pcondrInvocation setArgument:&ipar atIndex:5];
                [pcondrInvocation invoke];
                
                [matvecInvocation setArgument:&solution atIndex:2];
                [matvecInvocation setArgument:&s atIndex:3];
                [matvecInvocation setArgument:&t1v atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    t1v[i] = t1v[i] - matContainers->RHS[i];
                }
                [pcondlInvocation setArgument:&solution atIndex:2];
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
                [mstopInvocation setArgument:&solution atIndex:2];
                [mstopInvocation setArgument:&varContainers->Values atIndex:3];
                [mstopInvocation setArgument:&matContainers->RHS atIndex:4];
                [mstopInvocation setArgument:&r atIndex:5];
                [mstopInvocation setArgument:&ipar atIndex:6];
                [mstopInvocation invoke];
                [mstopInvocation getReturnValue:&residual];
                break;
            default:
                [pcondrInvocation setArgument:&solution atIndex:2];
                [pcondrInvocation setArgument:&s atIndex:3];
                [pcondrInvocation setArgument:&varContainers->Values atIndex:4];
                [pcondrInvocation setArgument:&ipar atIndex:5];
                [pcondrInvocation invoke];

                [matvecInvocation setArgument:&solution atIndex:2];
                [matvecInvocation setArgument:&s atIndex:3];
                [matvecInvocation setArgument:&t1v atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    t1v[i] = t1v[i] - matContainers->RHS[i];
                }
                [pcondlInvocation setArgument:&solution atIndex:2];
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
                NSLog(@"%d %lf\n", iter_count, residual);
            }
        }
        
        if (residual < HUTI_TOLERANCE) {
            HUTI_INFO = HUTI_CONVERGENCE;
            break;
        }
        
        for (i=0; i<ndim; i++) {
            u[i] = u[i] - omega1 * v[i] - omega2 * w[i];
        }
        
        // Return back to the top of the iteration loop or exit
        iter_count = iter_count + 1;
        if (iter_count > HUTI_MAXIT) {
            HUTI_INFO = HUTI_MAXITER;
            break;
        }
        
    }
    
    // We have exited the loop after enough iterations or broke down
    if (HUTI_DBUGLVL != HUTI_NO_DEBUG) {
        NSLog(@"%d %lf\n", iter_count, residual);
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

    matContainers = NULL;
    varContainers = NULL;
}

#pragma mark TFQMR

-(void)dtfqmrSolveInSolution:(FEMSolution *)solution ndim:(int)ndim wrkdim:(int)wrkdim ipar:(int *)ipar dpar:(double *)dpar work:(double **)work pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod {
/**********************************************************************************************************************************
 
    This method is based on the Roland W. Freund paper: "A Transpose-Free Quasi-Minimal Residual Algorithm for Non-Hermitian
    Linear Systems" (SIAM J. Sci. Comput., March 1993).
 
    All matrix-vector operations are done externally, so we don't need to know about the matrix structure (sparse or dense).
    Memory allocationf for the work array has also been done externally.
 
    Double precision version.
 
**********************************************************************************************************************************/
    
    int i, iter_count;
    double rho, oldrho=0, alpha, beta, gamma, oldgamma, eta, tau, c;
    double *v, *y, *ynew, *rtld, *t1v, *t2v, *w, *d, *r, *trv;
    double residual, upperb, rhsnorm, precrhsnorm;
    matrixArraysContainer *matContainers;
    variableArraysContainer *varContainers;
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
    
    NSMethodSignature *pCondlSignature, *pCondrSignature, *matvecSignature, *mstopSignature;
    NSInvocation *pcondlInvocation, *pcondrInvocation, *matvecInvocation, *mstopInvocation;
    
    FEMPrecondition *preconditioning;
    FEMKernel *kernel;
    
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
    
    mstopSignature = [FEMKernel instanceMethodSignatureForSelector:mstopMethod];
    mstopInvocation = [NSInvocation invocationWithMethodSignature:mstopSignature];
    [mstopInvocation setSelector:mstopMethod];
    
    // Instanciate class for the preconditioning
    preconditioning = [[FEMPrecondition alloc] init];
    
    //Instanciate class for the kernel
    kernel = [[FEMKernel alloc] init];
    
    // Targets to invocations
    [pcondlInvocation setTarget:preconditioning];
    [pcondrInvocation setTarget:preconditioning];
    [matvecInvocation setTarget:preconditioning];
    [mstopInvocation setTarget:kernel];
        
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
        rhsnorm = cblas_dnrm2(HUTI_NDIM, matContainers->RHS, 1);
    }
    if (HUTI_STOPC == HUTI_PRESID_SCALED_BYPRECB) {
        [pcondlInvocation setArgument:&solution atIndex:2];
        [pcondlInvocation setArgument:&d atIndex:3];
        [pcondlInvocation setArgument:&matContainers->RHS atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];
        
        precrhsnorm = cblas_dnrm2(HUTI_NDIM, d, 1);

    }
    
    // Part 1A - 1C
    
    //  Generate vector xvec if needed
    if (HUTI_INITIALX == HUTI_RANDOMX) {
        [self HUTI_dRandVector:varContainers->Values ipar:ipar];
    } else if (HUTI_INITIALX != HUTI_USERSUPPLIEDX) {
        for (i=0; i<varContainers->sizeValues; i++) {
            varContainers->Values[i] = 1;
        }
    }
    
    [pcondrInvocation setArgument:&solution atIndex:2];
    [pcondrInvocation setArgument:&d atIndex:3];
    [pcondrInvocation setArgument:&varContainers->Values atIndex:4];
    [pcondrInvocation setArgument:&ipar atIndex:5];
    [pcondrInvocation invoke];
    
    [matvecInvocation setArgument:&solution atIndex:2];
    [matvecInvocation setArgument:&d atIndex:3];
    [matvecInvocation setArgument:&r atIndex:4];
    [matvecInvocation setArgument:&ipar atIndex:5];
    [matvecInvocation invoke];
    
    for (i=0; i<ndim; i++) {
        d[i] = matContainers->RHS[i] - r[i];
    }
    
    [pcondlInvocation setArgument:&solution atIndex:2];
    [pcondlInvocation setArgument:&r atIndex:3];
    [pcondlInvocation setArgument:&d atIndex:4];
    [pcondlInvocation setArgument:&ipar atIndex:5];
    [pcondlInvocation invoke];
    
    for (i=0; i<ndim; i++) {
        y[i] = r[i];
        w[i] = r[i];
    }
    [pcondrInvocation setArgument:&solution atIndex:2];
    [pcondrInvocation setArgument:&v atIndex:3];
    [pcondrInvocation setArgument:&y atIndex:4];
    [pcondrInvocation setArgument:&ipar atIndex:5];
    [pcondrInvocation invoke];
    
    [matvecInvocation setArgument:&solution atIndex:2];
    [matvecInvocation setArgument:&v atIndex:3];
    [matvecInvocation setArgument:&d atIndex:4];
    [matvecInvocation setArgument:&ipar atIndex:5];
    [matvecInvocation invoke];
    
    [pcondlInvocation setArgument:&solution atIndex:2];
    [pcondlInvocation setArgument:&v atIndex:3];
    [pcondlInvocation setArgument:&d atIndex:4];
    [pcondlInvocation setArgument:&ipar atIndex:5];
    [pcondlInvocation invoke];
    for (i=0; i<ndim; i++) {
        t2v[i] = v[i];
    }
    
    memset( d, 0.0, (ndim*sizeof(d)) );
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
            varContainers->Values[i] = varContainers->Values[i] + eta * d[i];
        }
        
        oldgamma = gamma;
        
        // Check the convergence against selected stopping criterion
        switch (HUTI_STOPC) {
            case HUTI_TRUERESIDUAL:
                [pcondrInvocation setArgument:&solution atIndex:2];
                [pcondrInvocation setArgument:&trv atIndex:3];
                [pcondrInvocation setArgument:&varContainers->Values atIndex:4];
                [pcondrInvocation setArgument:&ipar atIndex:5];
                [pcondrInvocation invoke];
                
                [matvecInvocation setArgument:&solution atIndex:2];
                [matvecInvocation setArgument:&trv atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    trv[i] = r[i] - matContainers->RHS[i];
                }
                [pcondlInvocation setArgument:&solution atIndex:2];
                [pcondlInvocation setArgument:&r atIndex:3];
                [pcondlInvocation setArgument:&trv atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, r, 1);
                break;
            case HUTI_TRESID_SCALED_BYB:
                [pcondrInvocation setArgument:&solution atIndex:2];
                [pcondrInvocation setArgument:&trv atIndex:3];
                [pcondrInvocation setArgument:&varContainers->Values atIndex:4];
                [pcondrInvocation setArgument:&ipar atIndex:5];
                [pcondrInvocation invoke];
                
                [matvecInvocation setArgument:&solution atIndex:2];
                [matvecInvocation setArgument:&trv atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    trv[i] = r[i] - matContainers->RHS[i];
                }
                [pcondlInvocation setArgument:&solution atIndex:2];
                [pcondlInvocation setArgument:&r atIndex:3];
                [pcondlInvocation setArgument:&trv atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, r, 1) / rhsnorm;
                break;
            case HUTI_PSEUDORESIDUAL:
                [matvecInvocation setArgument:&solution atIndex:2];
                [matvecInvocation setArgument:&varContainers->Values atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    r[i] = r[i] - matContainers->RHS[i];
                }
                [pcondlInvocation setArgument:&solution atIndex:2];
                [pcondlInvocation setArgument:&trv atIndex:3];
                [pcondlInvocation setArgument:&r atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, trv, 1);
                break;
            case HUTI_PRESID_SCALED_BYB:
                [matvecInvocation setArgument:&solution atIndex:2];
                [matvecInvocation setArgument:&varContainers->Values atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    r[i] = r[i] - matContainers->RHS[i];
                }
                [pcondlInvocation setArgument:&solution atIndex:2];
                [pcondlInvocation setArgument:&trv atIndex:3];
                [pcondlInvocation setArgument:&r atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, trv, 1) / rhsnorm;
                break;
            case HUTI_PRESID_SCALED_BYPRECB:
                [matvecInvocation setArgument:&solution atIndex:2];
                [matvecInvocation setArgument:&varContainers->Values atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    r[i] = r[i] - matContainers->RHS[i];
                }
                [pcondlInvocation setArgument:&solution atIndex:2];
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
                    [pcondrInvocation setArgument:&solution atIndex:2];
                    [pcondrInvocation setArgument:&trv atIndex:3];
                    [pcondrInvocation setArgument:&varContainers->Values atIndex:4];
                    [pcondrInvocation setArgument:&ipar atIndex:5];
                    [pcondrInvocation invoke];
                    
                    [matvecInvocation setArgument:&solution atIndex:2];
                    [matvecInvocation setArgument:&trv atIndex:3];
                    [matvecInvocation setArgument:&r atIndex:4];
                    [matvecInvocation setArgument:&ipar atIndex:5];
                    [matvecInvocation invoke];
                    for (i=0; i<ndim; i++) {
                        trv[i] = r[i] - matContainers->RHS[i];
                    }
                    [pcondlInvocation setArgument:&solution atIndex:2];
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
                [mstopInvocation setArgument:&solution atIndex:2];
                [mstopInvocation setArgument:&varContainers->Values atIndex:3];
                [mstopInvocation setArgument:&matContainers->RHS atIndex:4];
                [mstopInvocation setArgument:&r atIndex:5];
                [mstopInvocation setArgument:&ipar atIndex:6];
                [mstopInvocation invoke];
                [mstopInvocation getReturnValue:&residual];
                break;
            default:
                [pcondrInvocation setArgument:&solution atIndex:2];
                [pcondrInvocation setArgument:&trv atIndex:3];
                [pcondrInvocation setArgument:&varContainers->Values atIndex:4];
                [pcondrInvocation setArgument:&ipar atIndex:5];
                [pcondrInvocation invoke];
                
                [matvecInvocation setArgument:&solution atIndex:2];
                [matvecInvocation setArgument:&trv atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    trv[i] = r[i] - matContainers->RHS[i];
                }
                [pcondlInvocation setArgument:&solution atIndex:2];
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
        
        [pcondrInvocation setArgument:&solution atIndex:2];
        [pcondrInvocation setArgument:&t1v atIndex:3];
        [pcondrInvocation setArgument:&ynew atIndex:4];
        [pcondrInvocation setArgument:&ipar atIndex:5];
        [pcondrInvocation invoke];
        
        [matvecInvocation setArgument:&solution atIndex:2];
        [matvecInvocation setArgument:&t1v atIndex:3];
        [matvecInvocation setArgument:&r atIndex:4];
        [matvecInvocation setArgument:&ipar atIndex:5];
        [matvecInvocation invoke];
        
        [pcondlInvocation setArgument:&solution atIndex:2];
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
            varContainers->Values[i] = varContainers->Values[i] + eta * d[i];
        }
        
        oldgamma = gamma;
        
        // Check the convergence against selected stopping criterion
        switch (HUTI_STOPC) {
            case HUTI_TRUERESIDUAL:
                [pcondrInvocation setArgument:&solution atIndex:2];
                [pcondrInvocation setArgument:&trv atIndex:3];
                [pcondrInvocation setArgument:&varContainers->Values atIndex:4];
                [pcondrInvocation setArgument:&ipar atIndex:5];
                [pcondrInvocation invoke];
                
                [matvecInvocation setArgument:&solution atIndex:2];
                [matvecInvocation setArgument:&trv atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    trv[i] = r[i] - matContainers->RHS[i];
                }
                [pcondlInvocation setArgument:&solution atIndex:2];
                [pcondlInvocation setArgument:&r atIndex:3];
                [pcondlInvocation setArgument:&trv atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, r, 1);
                break;
            case HUTI_TRESID_SCALED_BYB:
                [pcondrInvocation setArgument:&solution atIndex:2];
                [pcondrInvocation setArgument:&trv atIndex:3];
                [pcondrInvocation setArgument:&varContainers->Values atIndex:4];
                [pcondrInvocation setArgument:&ipar atIndex:5];
                [pcondrInvocation invoke];
                
                [matvecInvocation setArgument:&solution atIndex:2];
                [matvecInvocation setArgument:&trv atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    trv[i] = r[i] - matContainers->RHS[i];
                }
                [pcondlInvocation setArgument:&solution atIndex:2];
                [pcondlInvocation setArgument:&r atIndex:3];
                [pcondlInvocation setArgument:&trv atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, r, 1) / rhsnorm;
                break;
            case HUTI_PSEUDORESIDUAL:
                [matvecInvocation setArgument:&solution atIndex:2];
                [matvecInvocation setArgument:&varContainers->Values atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    r[i] = r[i] - matContainers->RHS[i];
                }
                [pcondlInvocation setArgument:&solution atIndex:2];
                [pcondlInvocation setArgument:&trv atIndex:3];
                [pcondlInvocation setArgument:&r atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, trv, 1);
                break;
            case HUTI_PRESID_SCALED_BYB:
                [matvecInvocation setArgument:&solution atIndex:2];
                [matvecInvocation setArgument:&varContainers->Values atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    r[i] = r[i] - matContainers->RHS[i];
                }
                [pcondlInvocation setArgument:&solution atIndex:2];
                [pcondlInvocation setArgument:&trv atIndex:3];
                [pcondlInvocation setArgument:&r atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, trv, 1) / rhsnorm;
                break;
            case HUTI_PRESID_SCALED_BYPRECB:
                [matvecInvocation setArgument:&solution atIndex:2];
                [matvecInvocation setArgument:&varContainers->Values atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    r[i] = r[i] - matContainers->RHS[i];
                }
                [pcondlInvocation setArgument:&solution atIndex:2];
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
                    [pcondrInvocation setArgument:&solution atIndex:2];
                    [pcondrInvocation setArgument:&trv atIndex:3];
                    [pcondrInvocation setArgument:&varContainers->Values atIndex:4];
                    [pcondrInvocation setArgument:&ipar atIndex:5];
                    [pcondrInvocation invoke];
                    
                    [matvecInvocation setArgument:&solution atIndex:2];
                    [matvecInvocation setArgument:&trv atIndex:3];
                    [matvecInvocation setArgument:&r atIndex:4];
                    [matvecInvocation setArgument:&ipar atIndex:5];
                    [matvecInvocation invoke];
                    for (i=0; i<ndim; i++) {
                        trv[i] = r[i] - matContainers->RHS[i];
                    }
                    [pcondlInvocation setArgument:&solution atIndex:2];
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
                [mstopInvocation setArgument:&solution atIndex:2];
                [mstopInvocation setArgument:&varContainers->Values atIndex:3];
                [mstopInvocation setArgument:&matContainers->RHS atIndex:4];
                [mstopInvocation setArgument:&r atIndex:5];
                [mstopInvocation setArgument:&ipar atIndex:6];
                [mstopInvocation invoke];
                [mstopInvocation getReturnValue:&residual];
                break;
            default:
                [pcondrInvocation setArgument:&solution atIndex:2];
                [pcondrInvocation setArgument:&trv atIndex:3];
                [pcondrInvocation setArgument:&varContainers->Values atIndex:4];
                [pcondrInvocation setArgument:&ipar atIndex:5];
                [pcondrInvocation invoke];
                
                [matvecInvocation setArgument:&solution atIndex:2];
                [matvecInvocation setArgument:&trv atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    trv[i] = r[i] - matContainers->RHS[i];
                }
                [pcondlInvocation setArgument:&solution atIndex:2];
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
                    NSLog(@"%d %lf %lf\n", iter_count, residual, upperb);
                } else {
                    NSLog(@"%d %lf\n", iter_count, residual);
                }
            }
        }
        
        // Part 2C
        
        rho = cblas_ddot(HUTI_NDIM, rtld, 1, w, 1);
        beta = rho / oldrho;
        for (i=0; i<ndim; i++) {
            ynew[i] = w[i] + beta * ynew[i];
        }
        
        [pcondrInvocation setArgument:&solution atIndex:2];
        [pcondrInvocation setArgument:&t2v atIndex:3];
        [pcondrInvocation setArgument:&ynew atIndex:4];
        [pcondrInvocation setArgument:&ipar atIndex:5];
        [pcondrInvocation invoke];
        
        [matvecInvocation setArgument:&solution atIndex:2];
        [matvecInvocation setArgument:&t2v atIndex:3];
        [matvecInvocation setArgument:&r atIndex:4];
        [matvecInvocation setArgument:&ipar atIndex:5];
        [matvecInvocation invoke];
        
        [pcondlInvocation setArgument:&solution atIndex:2];
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
        
        // Return back to the top of the iteration loop or exit
        iter_count = iter_count + 1;
        if (iter_count > HUTI_MAXIT) {
            HUTI_INFO = HUTI_MAXITER;
            break;
        }
        
    }
    
 // We have exited the loop after enough iterations or broke down
jump:
    
    // Compute the unpreconditioned xvec
    [pcondrInvocation setArgument:&solution atIndex:2];
    [pcondrInvocation setArgument:&trv atIndex:3];
    [pcondrInvocation setArgument:&varContainers->Values atIndex:4];
    [pcondrInvocation setArgument:&ipar atIndex:5];
    [pcondrInvocation invoke];
    for (i=0; i<ndim; i++) {
        varContainers->Values[i] = trv[i];
    }
    
    if (HUTI_DBUGLVL != HUTI_NO_DEBUG) {
        if (HUTI_STOPC == HUTI_UPPERB_STOPC) {
            NSLog(@"%d %lf %lf\n", iter_count, residual, upperb);
        } else {
            NSLog(@"%d %lf\n", iter_count, residual);
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
    
    matContainers = NULL;
    varContainers = NULL;
}

#pragma mark CG

-(void)dcgSolveInSolution:(FEMSolution *)solution ndim:(int)ndim wrkdim:(int)wrkdim ipar:(int *)ipar dpar:(double *)dpar work:(double **)work pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod {
/**********************************************************************************************************************************
 
    This method is based on Barret et al. book: "Templates for the Solution of Linear Systems: Building blocks for 
    Iterative Methods", 1993.
 
    All matrix-vector operations are done externally, so we don't need to know about the matrix structure (sparse or dense).
    Memory allocationf for the work array has also been done externally.
 
    Double precision version.
 
 **********************************************************************************************************************************/
    
    int i, iter_count;
    double rho, oldrho, alpha, beta;
    double *z, *p, *q, *r;
    double residual, rhsnorm, precrhsnorm;
    matrixArraysContainer *matContainers;
    variableArraysContainer *varContainers;
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
    
    NSMethodSignature *pCondlSignature, *pCondrSignature, *matvecSignature, *mstopSignature;
    NSInvocation *pcondlInvocation, *pcondrInvocation, *matvecInvocation, *mstopInvocation;
    
    FEMPrecondition *preconditioning;
    FEMKernel *kernel;
    
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
    
    mstopSignature = [FEMKernel instanceMethodSignatureForSelector:mstopMethod];
    mstopInvocation = [NSInvocation invocationWithMethodSignature:mstopSignature];
    [mstopInvocation setSelector:mstopMethod];
    
    // Instanciate class for the preconditioning
    preconditioning = [[FEMPrecondition alloc] init];
    
    //Instanciate class for the kernel
    kernel = [[FEMKernel alloc] init];
    
    // Targets to invocations
    [pcondlInvocation setTarget:preconditioning];
    [pcondrInvocation setTarget:preconditioning];
    [matvecInvocation setTarget:preconditioning];
    [mstopInvocation setTarget:kernel];
    
    // Vector allocations
    z = doublevec(0, ndim-1);
    p = doublevec(0, ndim-1);
    q = doublevec(0, ndim-1);
    r = doublevec(0, ndim-1);
    
    // The actual CG begins here (look at the pseudo-code in the book, page 27)
    
    // First the initialization part
    iter_count = 1;
    
    // The following applies for all matrix operations in this solver
    HUTI_EXTOP_MATTYPE = HUTI_MAT_NOTTRPSED;
    
    // Norms of right-hand side vector are used in convergence tests
    if (HUTI_STOPC == HUTI_TRESID_SCALED_BYB || HUTI_STOPC == HUTI_PRESID_SCALED_BYB) {
        rhsnorm = cblas_dnrm2(HUTI_NDIM, matContainers->RHS, 1);
    }
    if (HUTI_STOPC == HUTI_PRESID_SCALED_BYPRECB) {
        [pcondlInvocation setArgument:&solution atIndex:2];
        [pcondlInvocation setArgument:&p atIndex:3];
        [pcondlInvocation setArgument:&matContainers->RHS atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];
        
        precrhsnorm = cblas_dnrm2(HUTI_NDIM, p, 1);
    }
    
    // Generate vector xvec if needed
    if (HUTI_INITIALX == HUTI_RANDOMX) {
        [self HUTI_dRandVector:varContainers->Values ipar:ipar];
    } else if (HUTI_INITIALX != HUTI_USERSUPPLIEDX) {
        for (i=0; i<varContainers->sizeValues; i++) {
            varContainers->Values[i] = 1;
        }
    }
    
    [matvecInvocation setArgument:&solution atIndex:2];
    [matvecInvocation setArgument:&varContainers->Values atIndex:3];
    [matvecInvocation setArgument:&r atIndex:4];
    [matvecInvocation setArgument:&ipar atIndex:5];
    [matvecInvocation invoke];
    
    for (i=0; i<ndim; i++) {
        r[i] = matContainers->RHS[i] - r[i];
    }
    
    // This is where the loop starts
    while (1) {
        
        [pcondlInvocation setArgument:&solution atIndex:2];
        [pcondlInvocation setArgument:&q atIndex:3];
        [pcondlInvocation setArgument:&r atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];
        
        [pcondrInvocation setArgument:&solution atIndex:2];
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
        
        [matvecInvocation setArgument:&solution atIndex:2];
        [matvecInvocation setArgument:&p atIndex:3];
        [matvecInvocation setArgument:&q atIndex:4];
        [matvecInvocation setArgument:&ipar atIndex:5];
        [matvecInvocation invoke];
        
        alpha = rho / cblas_ddot(HUTI_NDIM, p, 1, q, 1);
        
        for (i=0; i<ndim; i++) {
            varContainers->Values[i] = varContainers->Values[i] + alpha * p[i];
            r[i] = r[i] - alpha * q[i];
        }
        
        // Check the convergence against selected stopping criterion
        switch (HUTI_STOPC) {
            case HUTI_TRUERESIDUAL:
                [matvecInvocation setArgument:&solution atIndex:2];
                [matvecInvocation setArgument:&varContainers->Values atIndex:3];
                [matvecInvocation setArgument:&z atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    z[i] = z[i] - matContainers->RHS[i];
                }
                residual = cblas_dnrm2(HUTI_NDIM, z, 1);
                break;
            case HUTI_TRESID_SCALED_BYB:
                [matvecInvocation setArgument:&solution atIndex:2];
                [matvecInvocation setArgument:&varContainers->Values atIndex:3];
                [matvecInvocation setArgument:&z atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    z[i] = z[i] - matContainers->RHS[i];
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
                [mstopInvocation setArgument:&solution atIndex:2];
                [mstopInvocation setArgument:&varContainers->Values atIndex:3];
                [mstopInvocation setArgument:&matContainers->RHS atIndex:4];
                [mstopInvocation setArgument:&r atIndex:5];
                [mstopInvocation setArgument:&ipar atIndex:6];
                [mstopInvocation invoke];
                [mstopInvocation getReturnValue:&residual];
                break;
            default:
                [matvecInvocation setArgument:&solution atIndex:2];
                [matvecInvocation setArgument:&varContainers->Values atIndex:3];
                [matvecInvocation setArgument:&z atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    z[i] = z[i] - matContainers->RHS[i];
                }
                residual = cblas_dnrm2(HUTI_NDIM, z, 1);
                break;
        }
        
        // Print debugging info if required
        if (HUTI_DBUGLVL != HUTI_NO_DEBUG) {
            if ( (iter_count % HUTI_DBUGLVL) == 0 ) {
                NSLog(@"%d %lf\n", iter_count, residual);
            }
        }
        
        if (residual < HUTI_TOLERANCE) {
            HUTI_INFO = HUTI_CONVERGENCE;
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
        NSLog(@"%d %lf\n", iter_count, residual);
    }
    
    HUTI_ITERS = iter_count;
    
    // Release memory
    free_dvector(z, 0, ndim-1);
    free_dvector(p, 0, ndim-1);
    free_dvector(q, 0, ndim-1);
    free_dvector(r, 0, ndim-1);
    
    matContainers = NULL;
    varContainers = NULL;
}

#pragma mark CGS

-(void)dcgsSolveInSolution:(FEMSolution *)solution ndim:(int)ndim wrkdim:(int)wrkdim ipar:(int *)ipar dpar:(double *)dpar work:(double **)work pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod {
/**********************************************************************************************************************************
 
    This method is based on Barret et al. book: "Templates for the Solution of Linear Systems: Building blocks for 
    Iterative Methods", 1993.
 
    All matrix-vector operations are done externally, so we don't need to know about the matrix structure (sparse or dense).
    Memory allocationf for the work array has also been done externally.
 
    Double precision version.
 
 **********************************************************************************************************************************/
    
    int i, iter_count;
    double rho, oldrho, alpha, beta;
    double *rtld, *p, *q, *u, *t1v, *t2v, *r;
    double residual, rhsnorm, precrhsnorm;
    matrixArraysContainer *matContainers;
    variableArraysContainer *varContainers;
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
    
    NSMethodSignature *pCondlSignature, *pCondrSignature, *matvecSignature, *mstopSignature;
    NSInvocation *pcondlInvocation, *pcondrInvocation, *matvecInvocation, *mstopInvocation;
    
    FEMPrecondition *preconditioning;
    FEMKernel *kernel;
    
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
    
    mstopSignature = [FEMKernel instanceMethodSignatureForSelector:mstopMethod];
    mstopInvocation = [NSInvocation invocationWithMethodSignature:mstopSignature];
    [mstopInvocation setSelector:mstopMethod];
    
    // Instanciate class for the preconditioning
    preconditioning = [[FEMPrecondition alloc] init];
    
    //Instanciate class for the kernel
    kernel = [[FEMKernel alloc] init];
    
    // Targets to invocations
    [pcondlInvocation setTarget:preconditioning];
    [pcondrInvocation setTarget:preconditioning];
    [matvecInvocation setTarget:preconditioning];
    [mstopInvocation setTarget:kernel];

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
        rhsnorm = cblas_dnrm2(HUTI_NDIM, matContainers->RHS, 1);
    }
    if (HUTI_STOPC == HUTI_PRESID_SCALED_BYPRECB) {
        [pcondlInvocation setArgument:&solution atIndex:2];
        [pcondlInvocation setArgument:&p atIndex:3];
        [pcondlInvocation setArgument:&matContainers->RHS atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];
        
        precrhsnorm = cblas_dnrm2(HUTI_NDIM, p, 1);
    }
    
    // Generate vector xvec if needed
    if (HUTI_INITIALX == HUTI_RANDOMX) {
        [self HUTI_dRandVector:varContainers->Values ipar:ipar];
    } else if (HUTI_INITIALX != HUTI_USERSUPPLIEDX) {
        for (i=0; i<varContainers->sizeValues; i++) {
            varContainers->Values[i] = 1;
        }
    }
    
    [matvecInvocation setArgument:&solution atIndex:2];
    [matvecInvocation setArgument:&varContainers->Values atIndex:3];
    [matvecInvocation setArgument:&r atIndex:4];
    [matvecInvocation setArgument:&ipar atIndex:5];
    [matvecInvocation invoke];
    
    for (i=0; i<ndim; i++) {
        r[i] = matContainers->RHS[i] - r[i];
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
        
        [pcondlInvocation setArgument:&solution atIndex:2];
        [pcondlInvocation setArgument:&t2v atIndex:3];
        [pcondlInvocation setArgument:&p atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];
        
        [pcondrInvocation setArgument:&solution atIndex:2];
        [pcondrInvocation setArgument:&t1v atIndex:3];
        [pcondrInvocation setArgument:&t2v atIndex:4];
        [pcondrInvocation setArgument:&ipar atIndex:5];
        [pcondrInvocation invoke];
        
        [matvecInvocation setArgument:&solution atIndex:2];
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
        
        [pcondlInvocation setArgument:&solution atIndex:2];
        [pcondlInvocation setArgument:&u atIndex:3];
        [pcondlInvocation setArgument:&t2v atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];
        
        [pcondrInvocation setArgument:&solution atIndex:2];
        [pcondrInvocation setArgument:&t1v atIndex:3];
        [pcondrInvocation setArgument:&u atIndex:4];
        [pcondrInvocation setArgument:&ipar atIndex:5];
        [pcondrInvocation invoke];
        for (i=0; i<ndim; i++) {
            varContainers->Values[i] = varContainers->Values[i] + alpha * t1v[i];
        }
        
        [matvecInvocation setArgument:&solution atIndex:2];
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
                [matvecInvocation setArgument:&solution atIndex:2];
                [matvecInvocation setArgument:&varContainers->Values atIndex:3];
                [matvecInvocation setArgument:&t1v atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    t1v[i] = t1v[i] - matContainers->RHS[i];
                }
                residual = cblas_dnrm2(HUTI_NDIM, t1v, 1);
                break;
            case HUTI_TRESID_SCALED_BYB:
                [matvecInvocation setArgument:&solution atIndex:2];
                [matvecInvocation setArgument:&varContainers->Values atIndex:3];
                [matvecInvocation setArgument:&t1v atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    t1v[i] = t1v[i] - matContainers->RHS[i];
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
                [mstopInvocation setArgument:&solution atIndex:2];
                [mstopInvocation setArgument:&varContainers->Values atIndex:3];
                [mstopInvocation setArgument:&matContainers->RHS atIndex:4];
                [mstopInvocation setArgument:&r atIndex:5];
                [mstopInvocation setArgument:&ipar atIndex:6];
                [mstopInvocation invoke];
                [mstopInvocation getReturnValue:&residual];
                break;
            default:
                [matvecInvocation setArgument:&solution atIndex:2];
                [matvecInvocation setArgument:&varContainers->Values atIndex:3];
                [matvecInvocation setArgument:&t1v atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    t1v[i] = t1v[i] - matContainers->RHS[i];
                }
                residual = cblas_dnrm2(HUTI_NDIM, t1v, 1);
                break;
        }
        
        // Print debugging info if required
        if (HUTI_DBUGLVL != HUTI_NO_DEBUG) {
            if ( (iter_count % HUTI_DBUGLVL) == 0 ) {
                NSLog(@"%d %lf\n", iter_count, residual);
            }
        }
        
        if (residual < HUTI_TOLERANCE) {
            HUTI_INFO = HUTI_CONVERGENCE;
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
        NSLog(@"%d %lf\n", iter_count, residual);
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
    
    matContainers = NULL;
    varContainers = NULL;
}

#pragma mark GMRES

-(void)dgmresSolveInSolution:(FEMSolution *)solution ndim:(int)ndim wrkdim:(int)wrkdim ipar:(int *)ipar dpar:(double *)dpar work:(double **)work pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod{
/**********************************************************************************************************************************
 
    This method is based on Barret et al. book: "Templates for the Solution of Linear Systems: Building blocks for 
    Iterative Methods", 1993.
 
    All matrix-vector operations are done externally, so we don't need to know about the matrix structure (sparse or dense).
    Memory allocationf for the work array has also been done externally.
 
    Double precision version.
 
**********************************************************************************************************************************/
    
    int i, j, k, iter_count, m;
    int s_ind, vtmp_ind, v_ind;
    double bnrm, alpha, beta;
    double temp, temp2, error;
    double *w, *r, *s, *vtmp, *t1v, *v, *buffer;
    double **h, *cs, *sn, *y; 
    double residual, rhsnorm, precrhsnorm;
    matrixArraysContainer *matContainers;
    variableArraysContainer *varContainers;
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
    
    NSMethodSignature *pCondlSignature, *pCondrSignature, *matvecSignature, *mstopSignature;
    NSInvocation *pcondlInvocation, *pcondrInvocation, *matvecInvocation, *mstopInvocation;
    
    FEMPrecondition *preconditioning;
    FEMKernel *kernel;
    
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
    
    mstopSignature = [FEMKernel instanceMethodSignatureForSelector:mstopMethod];
    mstopInvocation = [NSInvocation invocationWithMethodSignature:mstopSignature];
    [mstopInvocation setSelector:mstopMethod];
    
    // Instanciate class for the preconditioning
    preconditioning = [[FEMPrecondition alloc] init];
    
    //Instanciate class for the kernel
    kernel = [[FEMKernel alloc] init];
    
    // Targets to invocations
    [pcondlInvocation setTarget:preconditioning];
    [pcondrInvocation setTarget:preconditioning];
    [matvecInvocation setTarget:preconditioning];
    [mstopInvocation setTarget:kernel];
    
    // Vector allocations
    w = doublevec(0, ndim-1);
    r = doublevec(0, ndim-1);
    s = doublevec(0, ndim-1);
    vtmp = doublevec(0, ndim-1);
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
    
    // The following applies for all matrix operations in this solver
    HUTI_EXTOP_MATTYPE = HUTI_MAT_NOTTRPSED;
    
    bnrm = cblas_dnrm2(HUTI_NDIM, matContainers->RHS, 1);
    
    // Norms of right-hand side vector are used in convergence tests
    if (HUTI_STOPC == HUTI_TRESID_SCALED_BYB || HUTI_STOPC == HUTI_PRESID_SCALED_BYB) {
        rhsnorm = bnrm;
    }
    if (HUTI_STOPC == HUTI_PRESID_SCALED_BYPRECB) {
        [pcondlInvocation setArgument:&solution atIndex:2];
        [pcondlInvocation setArgument:&t1v atIndex:3];
        [pcondlInvocation setArgument:&matContainers->RHS atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];
        
        precrhsnorm = cblas_dnrm2(HUTI_NDIM, t1v, 1);
    }

    // Generate vector xvec if needed
    if (HUTI_INITIALX == HUTI_RANDOMX) {
        [self HUTI_dRandVector:varContainers->Values ipar:ipar];
    } else if (HUTI_INITIALX != HUTI_USERSUPPLIEDX) {
        for (i=0; i<varContainers->sizeValues; i++) {
            varContainers->Values[i] = 1;
        }
    }
    
    [pcondrInvocation setArgument:&solution atIndex:2];
    [pcondrInvocation setArgument:&t1v atIndex:3];
    [pcondrInvocation setArgument:&varContainers->Values atIndex:4];
    [pcondrInvocation setArgument:&ipar atIndex:5];
    [pcondrInvocation invoke];
    
    [matvecInvocation setArgument:&solution atIndex:2];
    [matvecInvocation setArgument:&t1v atIndex:3];
    [matvecInvocation setArgument:&r atIndex:4];
    [matvecInvocation setArgument:&ipar atIndex:5];
    [matvecInvocation invoke];
    for (i=0; i<ndim; i++) {
        t1v[i] = matContainers->RHS[i] - r[i];
    }
    [pcondlInvocation setArgument:&solution atIndex:2];
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
    
    for (i=0; i<(HUTI_GMRES_RESTART+1); i++) {
        for (j=0; j<(HUTI_GMRES_RESTART+1); j++) {
            h[i][j] = 0.0;
        }
    }
    memset( cs, 0.0, ((HUTI_GMRES_RESTART+1)*sizeof(cs)) );
    memset( sn, 0.0, ((HUTI_GMRES_RESTART+1)*sizeof(sn)) );
    memset( vtmp, 0.0, (ndim*sizeof(vtmp)) );
    work[0][vtmp_ind-1] = 1.0;
    
    // This is where the loop starts
    while (1) {
        
        [pcondrInvocation setArgument:&solution atIndex:2];
        [pcondrInvocation setArgument:&t1v atIndex:3];
        [pcondrInvocation setArgument:&varContainers->Values atIndex:4];
        [pcondrInvocation setArgument:&ipar atIndex:5];
        [pcondrInvocation invoke];
        
        [matvecInvocation setArgument:&solution atIndex:2];
        [matvecInvocation setArgument:&t1v atIndex:3];
        [matvecInvocation setArgument:&r atIndex:4];
        [matvecInvocation setArgument:&ipar atIndex:5];
        [matvecInvocation invoke];
        for (i=0; i<ndim; i++) {
            t1v[i] = matContainers->RHS[i] - r[i];
        }
        [pcondlInvocation setArgument:&solution atIndex:2];
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
            s[i] = alpha * vtmp[i];
        }
        
        // Construct orthonormal
        
        for (i=0; i<m; i++) {
            
            for (j=0; j<ndim; j++) {
                buffer[j] = work[j][v_ind+i-1];
            }
            
            [pcondrInvocation setArgument:&solution atIndex:2];
            [pcondrInvocation setArgument:&w atIndex:3];
            [pcondrInvocation setArgument:&buffer atIndex:4];
            [pcondrInvocation setArgument:&ipar atIndex:5];
            [pcondrInvocation invoke];
            
            [matvecInvocation setArgument:&solution atIndex:2];
            [matvecInvocation setArgument:&w atIndex:3];
            [matvecInvocation setArgument:&t1v atIndex:4];
            [matvecInvocation setArgument:&ipar atIndex:5];
            [matvecInvocation invoke];
            
            [pcondlInvocation setArgument:&solution atIndex:2];
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
                h[k+1][i] = -1 * sn[k] * h[k][i] + cs[k] * h[k+1][i];
                h[k][i] = temp;
            }
            
            if (h[i+1][i] == 0.0) {
                cs[i] = 1;
                sn[i] = 0.0;
            } else {
                if ( fabs(h[i+1][i]) > fabs(h[i][i]) ) {
                    temp2 = h[i][i] / h[i+1][i];
                    sn[i] = 1 / sqrt( 1 + (temp2 * temp2) );
                    cs[i] = temp2 * sn[i];
                } else {
                    temp2 = h[i+1][i] / h[i][i];
                    cs[i] = 1 / sqrt( 1 + (temp2 * temp2) );
                    sn[i] = temp2 * cs[i];
                }
            }
            
            temp = cs[i] * work[i][s_ind-1];
            work[i+1][s_ind-1] = -1 * sn[i] * work[i][s_ind-1];
            work[i][s_ind-1] = temp;
            
            h[i][i] = ( cs[i] * h[i][i] ) + ( sn[i] * h[i+1][i] );
            h[i+1][i] = 0.0;
            
            error = fabs( work[i+1][s_ind-1] ) / bnrm;
            if (error < HUTI_TOLERANCE) {
                
                for (j=0; j<ndim; j++) {
                    buffer[j] = work[j][s_ind-1];
                }
                [self HUTI_dLuSolveAt:i luMatrix:h afterSolve:y rightHandSide:buffer];
                
                // Matrix-vector product
                cblas_dgemv(CblasRowMajor, CblasNoTrans, solution.matrix.numberOfRows, ((v_ind+i)-(v_ind+1-1))+1, 1.0, *work+(v_ind+1-1), wrkdim, y, 1, 1.0, varContainers->Values, 1);
                break;
            }
        }
        
        if (error < HUTI_TOLERANCE) {
            goto error_lower_than_tolerance;
        }
        
        for (j=0; j<ndim; j++) {
            buffer[j] = work[j][s_ind-1];
        }
        [self HUTI_dLuSolveAt:m-1 luMatrix:h afterSolve:y rightHandSide:buffer];
        // Matrix-vector product
        cblas_dgemv(CblasRowMajor, CblasNoTrans, solution.matrix.numberOfRows, (v_ind+m-1)-(v_ind+1-1), 1.0, *work+(v_ind+1-1), wrkdim, y, 1, 1.0, varContainers->Values, 1);
        
    error_lower_than_tolerance:
        
        // Check the convergence against selected stopping criterion
        switch (HUTI_STOPC) {
            case HUTI_TRUERESIDUAL:
                [pcondrInvocation setArgument:&solution atIndex:2];
                [pcondrInvocation setArgument:&t1v atIndex:3];
                [pcondrInvocation setArgument:&varContainers->Values atIndex:4];
                [pcondrInvocation setArgument:&ipar atIndex:5];
                [pcondrInvocation invoke];

                [matvecInvocation setArgument:&solution atIndex:2];
                [matvecInvocation setArgument:&t1v atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    t1v[i] = matContainers->RHS[i] - r[i];
                }
                [pcondlInvocation setArgument:&solution atIndex:2];
                [pcondlInvocation setArgument:&r atIndex:3];
                [pcondlInvocation setArgument:&t1v atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, r, 1);
                break;
            case HUTI_TRESID_SCALED_BYB:
                [pcondrInvocation setArgument:&solution atIndex:2];
                [pcondrInvocation setArgument:&t1v atIndex:3];
                [pcondrInvocation setArgument:&varContainers->Values atIndex:4];
                [pcondrInvocation setArgument:&ipar atIndex:5];
                [pcondrInvocation invoke];

                [matvecInvocation setArgument:&solution atIndex:2];
                [matvecInvocation setArgument:&t1v atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    t1v[i] = matContainers->RHS[i] - r[i];
                }
                [pcondlInvocation setArgument:&solution atIndex:2];
                [pcondlInvocation setArgument:&r atIndex:3];
                [pcondlInvocation setArgument:&t1v atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, r, 1) / rhsnorm;
                break;
            case HUTI_PSEUDORESIDUAL:
                [matvecInvocation setArgument:&solution atIndex:2];
                [matvecInvocation setArgument:&varContainers->Values atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    r[i] = r[i] - matContainers->RHS[i];
                }
                [pcondlInvocation setArgument:&solution atIndex:2];
                [pcondlInvocation setArgument:&t1v atIndex:3];
                [pcondlInvocation setArgument:&r atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, t1v, 1);
                break;
            case HUTI_PRESID_SCALED_BYB:
                [matvecInvocation setArgument:&solution atIndex:2];
                [matvecInvocation setArgument:&varContainers->Values atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    r[i] = r[i] - matContainers->RHS[i];
                }
                [pcondlInvocation setArgument:&solution atIndex:2];
                [pcondlInvocation setArgument:&t1v atIndex:3];
                [pcondlInvocation setArgument:&r atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, t1v, 1) / rhsnorm;
                break;
            case HUTI_PRESID_SCALED_BYPRECB:
                [matvecInvocation setArgument:&solution atIndex:2];
                [matvecInvocation setArgument:&varContainers->Values atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    r[i] = r[i] - matContainers->RHS[i];
                }
                [pcondlInvocation setArgument:&solution atIndex:2];
                [pcondlInvocation setArgument:&t1v atIndex:3];
                [pcondlInvocation setArgument:&r atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, t1v, 1) / precrhsnorm;
                break;
            case HUTI_XDIFF_NORM:
                for (i=0; i<ndim; i++) {
                    t1v[i] = buffer[i];
                }
                residual = cblas_dnrm2(HUTI_NDIM, t1v, 1);
                break;
            case HUTI_USUPPLIED_STOPC:
                [mstopInvocation setArgument:&solution atIndex:2];
                [mstopInvocation setArgument:&varContainers->Values atIndex:3];
                [mstopInvocation setArgument:&matContainers->RHS atIndex:4];
                [mstopInvocation setArgument:&r atIndex:5];
                [mstopInvocation setArgument:&ipar atIndex:6];
                [mstopInvocation invoke];
                [mstopInvocation getReturnValue:&residual];
                break;
            default:
                [pcondrInvocation setArgument:&solution atIndex:2];
                [pcondrInvocation setArgument:&t1v atIndex:3];
                [pcondrInvocation setArgument:&varContainers->Values atIndex:4];
                [pcondrInvocation setArgument:&ipar atIndex:5];
                [pcondrInvocation invoke];

                [matvecInvocation setArgument:&solution atIndex:2];
                [matvecInvocation setArgument:&t1v atIndex:3];
                [matvecInvocation setArgument:&r atIndex:4];
                [matvecInvocation setArgument:&ipar atIndex:5];
                [matvecInvocation invoke];
                for (i=0; i<ndim; i++) {
                    t1v[i] = r[i] - matContainers->RHS[i];
                }
                [pcondlInvocation setArgument:&solution atIndex:2];
                [pcondlInvocation setArgument:&r atIndex:3];
                [pcondlInvocation setArgument:&t1v atIndex:4];
                [pcondlInvocation setArgument:&ipar atIndex:5];
                [pcondlInvocation invoke];
                residual = cblas_dnrm2(HUTI_NDIM, r, 1);
                break;
        }
        
        work[(m+1)-1][s_ind-1] = cblas_dnrm2(HUTI_NDIM, r, 1);

        // Print debugging info if required
        if (HUTI_DBUGLVL != HUTI_NO_DEBUG) {
            if ( (iter_count % HUTI_DBUGLVL) == 0 ) {
                NSLog(@"%d %lf\n", iter_count, residual);
            }
        }
        
        if (residual < HUTI_TOLERANCE) {
            HUTI_INFO = HUTI_CONVERGENCE;
            break;
        }
        
        // Return back to the top of the iteration loop or exit
        iter_count = iter_count + 1;
        if (iter_count > HUTI_MAXIT) {
            HUTI_INFO = HUTI_MAXITER;
            break;
        }
        
    }
    
jump:
    
    // We have exited the loop after enough iterations or broke down
    if (HUTI_DBUGLVL != HUTI_NO_DEBUG) {
        NSLog(@"%d %lf\n", iter_count, residual);
    }
    
    HUTI_ITERS = iter_count;
    
    // Release memory
    free_dvector(w, 0, ndim-1);
    free_dvector(r, 0, ndim-1);
    free_dvector(s, 0, ndim-1);
    free_dvector(vtmp, 0, ndim-1);
    free_dvector(t1v, 0, ndim-1);
    free_dvector(v, 0, ndim-1);
    free_dvector(buffer, 0, ndim-1);
    free_dvector(cs, 0, (HUTI_GMRES_RESTART+1)-1);
    free_dvector(sn, 0, (HUTI_GMRES_RESTART+1)-1);
    free_dvector(y, 0, (HUTI_GMRES_RESTART+1)-1);
    
    free_dmatrix(h, 0, (HUTI_GMRES_RESTART+1)-1, 0, (HUTI_GMRES_RESTART+1)-1);
    
    matContainers = NULL;
    varContainers = NULL;
}

#pragma mark SGS

-(void)dsgsSolveInsolution:(FEMSolution *)solution ndim:(int)ndim wrkdim:(int)wrkdim ipar:(int *)ipar dpar:(double *)dpar work:(double **)work pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod {
    
    int rounds, outputInterval;
    double minTol, maxTol, residual, omega;
    BOOL converged, diverged;
    matrixArraysContainer *matContainers;
    variableArraysContainer *varContainers;
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
    
    rounds = ipar[9];
    minTol = dpar[0];
    maxTol = dpar[1];
    outputInterval = ipar[4];
    omega = dpar[2];
    
    [self HUTI_sgsNumberOfDimension:ndim solution:solution afterSolve:varContainers->Values rightHandSide:matContainers->RHS ipar:ipar rounds:rounds minTolerance:minTol maxTolerance:maxTol residual:&residual converged:&converged diverged:&diverged outputInterval:outputInterval omega:omega matvecMethod:matvecMethod];
    
    if (converged == YES) ipar[29] = 1;
    if (diverged == YES) ipar[29] = 3;
    
    matContainers = NULL;
    varContainers = NULL;
}

#pragma mark JACOBI

-(void)djacobiSolveInSolution:(FEMSolution *)solution ndim:(int)ndim wrkdim:(int)wrkdim ipar:(int *)ipar dpar:(double *)dpar work:(double **)work pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod {
    
    int rounds, outputInterval;
    double minTol, maxTol, residual;
    BOOL converged, diverged;
    matrixArraysContainer *matContainers;
    variableArraysContainer *varContainers;
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;

    rounds = ipar[9];
    minTol = dpar[0];
    maxTol = dpar[1];
    outputInterval = ipar[4];
    
    [self HUTI_jacobiNumberOfDimension:ndim solution:solution afterSolve:varContainers->Values rightHandSide:matContainers->RHS ipar:ipar rounds:rounds minTolerance:minTol maxTolerance:maxTol residual:&residual converged:&converged diverged:&diverged outputInterval:outputInterval matvecMethod:matvecMethod];
    
    if (converged == YES) ipar[29] = 1;
    if (diverged == YES) ipar[29] = 3;
    
    matContainers = NULL;
    varContainers = NULL;
}

#pragma mark BI-CGSTAB(l)

-(void)dbicgstablSolveInSolution:(FEMSolution *)solution ndim:(int)ndim wrkdim:(int)wrkdim ipar:(int *)ipar dpar:(double *)dpar work:(double **)work pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod {
    
    int rounds, outputInterval, polynomialDegree;
    double minTol, maxTol;
    BOOL converged, diverged;
    matrixArraysContainer *matContainers;
    variableArraysContainer *varContainers;
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
    
    //TODO: Add support for constraint matrix
    
    rounds = ipar[9];
    minTol = dpar[0];
    maxTol = dpar[1];
    outputInterval = ipar[4];
    polynomialDegree = ipar[15];
    
    [self HUTI_BICGStabLNumberOfDimension:ndim solution:solution afterSolve:varContainers->Values rightHandSide:matContainers->RHS ipar:ipar rounds:rounds minTolerance:minTol maxTolerance:maxTol converged:&converged diverged:&diverged outputInterval:outputInterval polynomialDegree:polynomialDegree pcondlMethod:pcondlMethod matvecMethod:matvecMethod];
    
    if (converged == YES) ipar[29] = 1;
    if (diverged == YES) ipar[29] = 3;
    
    matContainers = NULL;
    varContainers = NULL;
}

#pragma mark GCR

-(void)dgcrSolveInSolution:(FEMSolution *)solution ndim:(int)ndim wrkdim:(int)wrkdim ipar:(int *)ipar dpar:(double *)dpar work:(double **)work pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod {
    
    int rounds, outputInterval, restartN;
    double minTol, maxTol, residual;
    BOOL converged, diverged;
    matrixArraysContainer *matContainers;
    variableArraysContainer *varContainers;
    
    matContainers = solution.matrix.getContainers;
    varContainers = solution.variable.getContainers;
    
    rounds = ipar[9];
    minTol = dpar[0];
    maxTol = dpar[1];
    outputInterval = ipar[4];
    restartN = ipar[16];
    
    [self HUTI_gcrNumberOfDimension:ndim solution:solution afterSolve:varContainers->Values rightHandSide:matContainers->RHS ipar:ipar rounds:rounds minTolerance:minTol maxTolerance:maxTol residual:&residual converged:&converged diverged:&diverged outputInterval:outputInterval restart:restartN pcondlMethod:pcondlMethod matvecMethod:matvecMethod];
    
    if (converged == YES) ipar[29] = 1;
    if (diverged == YES) ipar[29] = 3;
    
    matContainers = NULL;
    varContainers = NULL;
}

@end
