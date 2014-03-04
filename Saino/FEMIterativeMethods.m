//
//  FEMIterativeMethods.m
//  Saino
//
//  Created by Seddik hakime on 10/06/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import <Accelerate/Accelerate.h>

#import "FEMIterativeMethods.h"
#import "FEMPrecondition.h"
#import "Utils.h"

@interface FEMIterativeMethods ()
-(void)FEMIterativeMethods_sgsNumberOfDimensions:(int)n matrix:(FEMMatrix *)matrix afterSolve:(double *)xvec rightHandSide:(double *)rhsvec ipar:(int *)ipar rounds:(double)rounds minTolerance:(double)minTolerance maxTolerance:(double)maxTolerance residual:(double *)residual converged:(BOOL *)converged diverged:(BOOL *)diverged outputInterval:(int)outputInterval omega:(double)omega matvecMethod:(SEL)matvecMethod;
-(void)FEMIterativeMethods_jacobiNumberOfDimensions:(int)n matrix:(FEMMatrix *)matrix afterSolve:(double *)xvec rightHandSide:(double *)rhsvec ipar:(int *)ipar rounds:(double)rounds minTolerance:(double)minTolerance maxTolerance:(double)maxTolerance residual:(double *)residual converged:(BOOL *)converged diverged:(BOOL *)diverged outputInterval:(int)outputInterval matvecMethod:(SEL)matvecMethod;
-(void)FEMIterativeMethods_BICGStabLNumberOfDimensions:(int)n matrix:(FEMMatrix *)matrix afterSolve:(double *)xvec rightHandSide:(double *)rhsvec ipar:(int *)ipar rounds:(double)rounds minTolerance:(double)minTolerance maxTolerance:(double)maxTolerance converged:(BOOL *)converged diverged:(BOOL *)diverged outputInterval:(int)outputInterval polynomialDegree:(int)polyDegree constrained:(BOOL)constrained numberOfConstrains:(int)nc pcondlMethod:(SEL)pcondlMethod matvecMethod:(SEL)matvecMethod;
-(void)FEMIterativeMethods_gcrNumberOfDimensions:(int)n matrix:(FEMMatrix *)matrix afterSolve:(double *)xvec rightHandSide:(double *)rhsvec ipar:(int *)ipar rounds:(double)rounds minTolerance:(double)minTolerance maxTolerance:(double)maxTolerance residual:(double *)residual converged:(BOOL *)converged diverged:(BOOL *)diverged outputInterval:(int)outputInterval restart:(int)m pcondlMethod:(SEL)pcondlMethod matvecMethod:(SEL)matvecMethod;
-(void)FEMIterativeMethods_richardsonNumberOfDimensions:(int)n matrix:(FEMMatrix *)matrix afterSolve:(double *)xvec rightHandSide:(double *)rhsvec ipar:(int *)ipar rounds:(double)rounds minTolerance:(double)minTolerance maxTolerance:(double)maxTolerance residual:(double *)residual converged:(BOOL *)converged diverged:(BOOL *)diverged outputInterval:(int)outputInterval pcondlMethod:(SEL)pcondlMethod matvecMethod:(SEL)matvecMethod;
-(void)FEMIterativeMethods_ConstraintVectorMultiplyInMatrix:(FEMMatrix *)constraintMatrix constrained:(BOOL)constrained numberOfDimensions:(int)n multiplyVector:(double *)u resultVector:(double *)v;


@end

@implementation FEMIterativeMethods

#pragma mark Private methods

-(void)FEMIterativeMethods_sgsNumberOfDimensions:(int)n matrix:(FEMMatrix *)matrix afterSolve:(double *)xvec rightHandSide:(double *)rhsvec ipar:(int *)ipar rounds:(double)rounds minTolerance:(double)minTolerance maxTolerance:(double)maxTolerance residual:(double *)residual converged:(BOOL *)converged diverged:(BOOL *)diverged outputInterval:(int)outputInterval omega:(double)omega matvecMethod:(SEL)matvecMethod {
    
    int i, j, k;
    double *r;
    double bnorm, rnorm, s;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = matrix.getContainers;
    
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
    
    r = doublevec(0, n-1);
    
    [matvecInvocation setArgument:&matrix atIndex:2];
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
        
        [matvecInvocation setArgument:&matrix atIndex:2];
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
            NSLog(@"FEMIterativeMethods:FEMIterativeMethods_sgsNumberOfDimension: %d %11.4e %11.4e\n", k, rnorm, *residual);
        }
        
        *converged = (*residual < minTolerance) ? YES : NO;
        *diverged = (*residual > maxTolerance) ? YES : NO;
        if (*converged == YES || *diverged == YES) break;
    }
    
    free_dvector(r, 0, n-1);
}

-(void)FEMIterativeMethods_jacobiNumberOfDimensions:(int)n matrix:(FEMMatrix *)matrix afterSolve:(double *)xvec rightHandSide:(double *)rhsvec ipar:(int *)ipar rounds:(double)rounds minTolerance:(double)minTolerance maxTolerance:(double)maxTolerance residual:(double *)residual converged:(BOOL *)converged diverged:(BOOL *)diverged outputInterval:(int)outputInterval matvecMethod:(SEL)matvecMethod {
    
    int i, k;
    double *r;
    double bnorm, rnorm;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = matrix.getContainers;
    
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
    
    r = doublevec(0, n-1);
    
    [matvecInvocation setArgument:&matrix atIndex:2];
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
        [matvecInvocation setArgument:&matrix atIndex:2];
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
            NSLog(@"FEMIterativeMethods:FEMIterativeMethods_jacobiNumberOfDimension: %d %11.4e %11.4e\n", k, rnorm, *residual);
        }
        
        *converged = (*residual < minTolerance) ? YES : NO;
        *diverged = (*residual > maxTolerance) ? YES : NO;
        if (*converged == YES || *diverged == YES) break;
    }
    
    free_dvector(r, 0, n-1);
}

-(void)FEMIterativeMethods_ConstraintVectorMultiplyInMatrix:(FEMMatrix *)constraintMatrix constrained:(BOOL)constrained numberOfDimensions:(int)n multiplyVector:(double *)u resultVector:(double *)v {
    
    int nc;
    double s;
    matrixArraysContainer *matContainers = NULL;
    
    if (constrained == NO) return;
    
    matContainers = constraintMatrix.getContainers;
    
    nc = constraintMatrix.numberOfRows;
    s = 1.0;
    for (int i=0; i<nc; i++) {
        v[n-nc+i] = 0.0;
        for (int j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
            if (matContainers->DiagScaling != NULL) s = matContainers->DiagScaling[matContainers->Cols[j]];
            v[n-nc+i] = v[n-nc+i] + s*matContainers->Values[j]*u[matContainers->Cols[j]];
            v[matContainers->Cols[j]] = v[matContainers->Cols[j]] + s*matContainers->Values[j]*u[n-nc+i];
        }
    }
}

/**********************************************************************************************************************
    The subroutine has been written using as a starting point the work of D.R. Fokkema (subroutine zbistbl v1.1 1998). 
    Dr. Fokkema has given the right to distribute the derived work under GPL and hence the original more conservative
    copyright notice of the subroutine has been removed accordingly.
 
    Double precision version.
***********************************************************************************************************************/
-(void)FEMIterativeMethods_BICGStabLNumberOfDimensions:(int)n matrix:(FEMMatrix *)matrix afterSolve:(double *)xvec rightHandSide:(double *)rhsvec ipar:(int *)ipar rounds:(double)rounds minTolerance:(double)minTolerance maxTolerance:(double)maxTolerance converged:(BOOL *)converged diverged:(BOOL *)diverged outputInterval:(int)outputInterval polynomialDegree:(int)polyDegree constrained:(BOOL)constrained numberOfConstrains:(int)nc pcondlMethod:(SEL)pcondlMethod matvecMethod:(SEL)matvecMethod {
    
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
        errorfunct("FEMIterativeMethods:FEMIterativeMethods_BICGStabLNumberOfDimension", "Polynomial degree < 2.");
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
    
    if (all(xvec, '=', 0.0, n) == true) {
        for (i=0; i<n; i++) {
            xvec[i] = rhsvec[i];
        }
    }
    
    zero = 0.0;
    one = 1.0;
    
    memset( *work, 0.0, (n*(3+2*(polyDegree+1)))*sizeof(double) );
    memset( *rwork, 0.0, ((polyDegree+1)*(3+2*(polyDegree+1)))*sizeof(double) );
    
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
    
    
    memset( buffer, 0.0, n*sizeof(double) );
    
    [matvecInvocation setArgument:&matrix atIndex:2];
    [matvecInvocation setArgument:&xvec atIndex:3];
    [matvecInvocation setArgument:&buffer atIndex:4];
    [matvecInvocation setArgument:&ipar atIndex:5];
    [matvecInvocation invoke];
    if (constrained == YES)
        [self FEMIterativeMethods_ConstraintVectorMultiplyInMatrix:matrix.constraint constrained:constrained numberOfDimensions:n multiplyVector:xvec resultVector:buffer];
    for (i=0; i<n; i++) {
        work[i][r] = buffer[i];
    }
    for (i=0; i<n; i++) {
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
        free_ivector(iwork, 0, (polyDegree-1)-1);
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
                errorfunct("FEMIterativeMethods:FEMIterativeMethods_BICGStabLNumberOfDimension", "Breakdown error.");
            }
            
            beta = alpha * (rho1/rho0);
            rho0 = rho1;
            for (j=0; j<=k-1; j++) {
                for (i=0; i<n; i++) {
                    work[i][u+j] = work[i][r+j] - beta * work[i][u+j];
                }
            }
            
            for (i=0; i<n-nc; i++) {
                buffer[i] = work[i][u+k-1];
            }
            [pcondlInvocation setArgument:&matrix atIndex:2];
            [pcondlInvocation setArgument:&t atIndex:3];
            [pcondlInvocation setArgument:&buffer atIndex:4];
            [pcondlInvocation setArgument:&ipar atIndex:5];
            [pcondlInvocation invoke];
            
            [matvecInvocation setArgument:&matrix atIndex:2];
            [matvecInvocation setArgument:&t atIndex:3];
            [matvecInvocation setArgument:&buffer atIndex:4];
            [matvecInvocation setArgument:&ipar atIndex:5];
            [matvecInvocation invoke];
            if (constrained == YES) {
                for (i=n-nc; i<n; i++) {
                    t[i] = work[i][u+k-1];
                }
                [self FEMIterativeMethods_ConstraintVectorMultiplyInMatrix:matrix.constraint constrained:constrained numberOfDimensions:n multiplyVector:t resultVector:buffer];
            }
            for (i=0; i<n; i++) {
                work[i][u+k] = buffer[i];
            }
            
            for (i=0; i<n; i++) {
                buffer[i] = work[i][rr];
                buffer2[i] = work[i][u+k];
            }
            sigma = cblas_ddot(n, buffer, 1, buffer2, 1);
            if (sigma == zero) {
                errorfunct("FEMIterativeMethods:FEMIterativeMethods_BICGStabLNumberOfDimension", "Breakdown error.");
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
            [pcondlInvocation setArgument:&matrix atIndex:2];
            [pcondlInvocation setArgument:&t atIndex:3];
            [pcondlInvocation setArgument:&buffer atIndex:4];
            [pcondlInvocation setArgument:&ipar atIndex:5];
            [pcondlInvocation invoke];

            [matvecInvocation setArgument:&matrix atIndex:2];
            [matvecInvocation setArgument:&t atIndex:3];
            [matvecInvocation setArgument:&buffer atIndex:4];
            [matvecInvocation setArgument:&ipar atIndex:5];
            [matvecInvocation invoke];
            if (constrained == YES) {
                for (i=n-nc; i<n; i++) {
                    t[i] = work[i][r+k-1];
                }
                [self FEMIterativeMethods_ConstraintVectorMultiplyInMatrix:matrix.constraint constrained:constrained numberOfDimensions:n multiplyVector:t resultVector:buffer];
            }
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
        
        for (i=0; i<polyDegree+1; i++) {
            k = z;
            for (j=zz; j<=zz+polyDegree; j++) {
                rwork[i][j] = rwork[i][k];
                k++;
            }
        }
        
        // Transfrom rwork for LAPACK, column-major order
        for (i=0; i<(3+2*(polyDegree+1)); i++) {
            for (j=0; j<polyDegree+1; j++) {
                rwork_transpose[j+(polyDegree+1)*i] = rwork[j][i];
            }
        }
        
        rows = polyDegree-1;
        cols = polyDegree-1;
        lda = polyDegree+1;
        dgetrf_(&rows, &cols, rwork_transpose+((zz+1)*lda+1), &lda, iwork, &info);
        
        // tild r0 and tild rl (small vectors)
        rwork_transpose[y0*(polyDegree+1)] = -one;
        for (i=1; i<polyDegree; i++) {
            rwork_transpose[y0*(polyDegree+1)+i] = rwork_transpose[z*(polyDegree+1)+i];
        }
        
        order = polyDegree-1;
        nrhs = 1;
        str = 'N';
        lda = polyDegree+1;
        dgetrs_(&str, &order, &nrhs, rwork_transpose+((zz+1)*lda+1), &lda, iwork, rwork_transpose+(y0*(polyDegree+1)+1), &lda, &info);
        rwork_transpose[y0*(polyDegree+1)+polyDegree] = zero;
        
        rwork_transpose[yl*(polyDegree+1)] = zero;
        for (i=1; i<polyDegree; i++) {
            rwork_transpose[yl*(polyDegree+1)+i] = rwork_transpose[(z+polyDegree)*(polyDegree+1)+i];
        }
        dgetrs_(&str, &order, &nrhs, rwork_transpose+((zz+1)*lda+1), &lda, iwork, rwork_transpose+(yl*(polyDegree+1)+1), &lda, &info);
        rwork_transpose[yl*(polyDegree+1)+polyDegree] = -one;
        
        // Back to the original rwork matrix
        for (i=0; i<(3+2*(polyDegree+1)); i++) {
            for (j=0; j<polyDegree+1; j++) {
                rwork[j][i] = rwork_transpose[j+(polyDegree+1)*i];
            }
        }
        
        // Convex combination
        
        memset( buffer, 0.0, n*sizeof(double) );
        memset( buffer2, 0.0, n*sizeof(double) );
        for (i=0; i<polyDegree+1; i++) {
            buffer[i] = rwork[i][y0];
        }
        cblas_dsymv(CblasRowMajor, CblasUpper, polyDegree+1, one, *rwork, (3+2*(polyDegree+1)), buffer, 1, zero, buffer2, 1);
        kappa0 = sqrt( cblas_ddot(polyDegree+1, buffer, 1, buffer2, 1) );
        
        memset( buffer, 0.0, n*sizeof(double) );
        for (i=0; i<polyDegree+1; i++) {
            buffer[i] = rwork[i][yl];
        }
        cblas_dsymv(CblasRowMajor, CblasUpper, polyDegree+1, one, *rwork, (3+2*(polyDegree+1)), buffer, 1, zero, buffer2, 1);
        kappal = sqrt( cblas_ddot(polyDegree+1, buffer, 1, buffer2, 1) );
        
        memset( buffer, 0.0, n*sizeof(double) );
        for (i=0; i<polyDegree+1; i++) {
            buffer[i] = rwork[i][y0];
        }
        cblas_dsymv(CblasRowMajor, CblasUpper, polyDegree+1, one, *rwork, (3+2*(polyDegree+1)), buffer, 1, zero, buffer2, 1);
        memset( buffer, 0.0, n*sizeof(double) );
        for (i=0; i<polyDegree+1; i++) {
            buffer[i] = rwork[i][yl];
        }
        varrho = cblas_ddot(polyDegree+1, buffer, 1, buffer2, 1) / (kappa0*kappal);
        hatgamma = varrho / fabs(varrho) * max(fabs(varrho), 7e-1) * kappa0/kappal;
        for (i=0; i<polyDegree+1; i++) {
            rwork[i][y0] = rwork[i][y0] - hatgamma * rwork[i][yl];
        }
        
        // Update
        
        omega = rwork[polyDegree][y0];
        for (j=1; j<=polyDegree; j++) {
            for (i=0; i<n; i++) {
                work[i][u] = work[i][u] - rwork[j][y0] * work[i][u+j];
                xvec[i] = xvec[i] + rwork[j][y0] * work[i][r+j-1];
                work[i][r] = work[i][r] - rwork[j][y0] * work[i][r+j];
            }
        }
        
        memset( buffer, 0.0, n*sizeof(double) );
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
            [pcondlInvocation setArgument:&matrix atIndex:2];
            [pcondlInvocation setArgument:&t atIndex:3];
            [pcondlInvocation setArgument:&xvec atIndex:4];
            [pcondlInvocation setArgument:&ipar atIndex:5];
            [pcondlInvocation invoke];
            
            memset( buffer, 0.0, n*sizeof(double) );
            [matvecInvocation setArgument:&matrix atIndex:2];
            [matvecInvocation setArgument:&t atIndex:3];
            [matvecInvocation setArgument:&buffer atIndex:4];
            [matvecInvocation setArgument:&ipar atIndex:5];
            [matvecInvocation invoke];
            if (constrained == YES) {
                for (i=n-nc; i<n; i++) {
                    t[i] = xvec[i];
                }
                [self FEMIterativeMethods_ConstraintVectorMultiplyInMatrix:matrix.constraint constrained:constrained numberOfDimensions:n multiplyVector:t resultVector:buffer];
            }
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
            [pcondlInvocation setArgument:&matrix atIndex:2];
            [pcondlInvocation setArgument:&t atIndex:3];
            [pcondlInvocation setArgument:&xvec atIndex:4];
            [pcondlInvocation setArgument:&ipar atIndex:5];
            [pcondlInvocation invoke];
            if (constrained == YES) {
                for (i=n-nc; i<n; i++) {
                    t[i] = xvec[i];
                }
            }
            for (i=0; i<n; i++) {
                t[i] = t[i] + work[i][xp];
            }
        }
        
        errorind = rnrm / bnrm;
        if (round % outputInterval == 0) {
            NSLog(@"FEMIterativeMethods:FEMIterativeMethods_BICGStabLNumberOfDimension: %d %11.4e %11.4e\n", round, rnrm, errorind);
        }
        
        *converged = (errorind < minTolerance) ? YES: NO;
        *diverged = (errorind > maxTolerance) ? YES : NO;
        if (*converged == YES || *diverged == YES) break;
        
    } // end of rounds
    
    if (outputInterval != HUGE_VAL) {
        NSLog(@"FEMIterativeMethods:FEMIterativeMethods_BICGStabLNumberOfDimension: %d %11.4e %11.4e\n", round, rnrm, errorind);
    }
    
    // We have solved z = P*x, with P the preconditioner, so finally
    // solve the true unknown x
    for (i=0; i<n; i++) {
        t[i] = xvec[i];
    }
    [pcondlInvocation setArgument:&matrix atIndex:2];
    [pcondlInvocation setArgument:&xvec atIndex:3];
    [pcondlInvocation setArgument:&t atIndex:4];
    [pcondlInvocation setArgument:&ipar atIndex:5];
    [pcondlInvocation invoke];
    if (constrained == YES) {
        for (i=n-nc; i<n; i++) {
            xvec[i] = t[i];
        }
    }
    for (i=0; i<n; i++) {
        xvec[i] = xvec[i] + work[i][xp];
    }
    
    free_dvector(t, 0, n-1);
    free_dmatrix(work, 0, n-1, 0, (3+2*(polyDegree+1))-1);
    free_dmatrix(rwork, 0, (polyDegree+1)-1, 0, (3+2*(polyDegree+1))-1);
    free_dvector(rwork_transpose, 0, ( (polyDegree+1)*(3+2*(polyDegree+1)) )-1);
    free_dvector(buffer, 0, n-1);
    free_dvector(buffer2, 0, n-1);
    free_ivector(iwork, 0, (polyDegree-1)-1);
}

-(void)FEMIterativeMethods_gcrNumberOfDimensions:(int)n matrix:(FEMMatrix *)matrix afterSolve:(double *)xvec rightHandSide:(double *)rhsvec ipar:(int *)ipar rounds:(double)rounds minTolerance:(double)minTolerance maxTolerance:(double)maxTolerance residual:(double *)residual converged:(BOOL *)converged diverged:(BOOL *)diverged outputInterval:(int)outputInterval restart:(int)m pcondlMethod:(SEL)pcondlMethod matvecMethod:(SEL)matvecMethod {
    
    int i, j, k, l;
    double bnorm, rnorm;
    double alpha, beta;
    double *r, **s, **v, *t1, *t2;
    double *buffer;
    
    NSMethodSignature *pCondlSignature, *matvecSignature;
    NSInvocation *pcondlInvocation, *matvecInvocation;
    
    FEMPrecondition *preconditioning;
    
    r = doublevec(0, n-1);
    t1 = doublevec(0, n-1);
    t2 = doublevec(0, n-1);
    
    if (m > 1) {
        v = doublematrix(0, n-1, 0, (m-1)-1);
        s = doublematrix(0, n-1, 0, (m-1)-1);
        memset( *v, 0.0, (n*(m-1))*sizeof(double) );
        memset( *s, 0.0, (n*(m-1))*sizeof(double) );
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
    
    [matvecInvocation setArgument:&matrix atIndex:2];
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
        
        [pcondlInvocation setArgument:&matrix atIndex:2];
        [pcondlInvocation setArgument:&t1 atIndex:3];
        [pcondlInvocation setArgument:&r atIndex:4];
        [pcondlInvocation setArgument:&ipar atIndex:5];
        [pcondlInvocation invoke];
        
        [matvecInvocation setArgument:&matrix atIndex:2];
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
            NSLog(@"FEMIterativeMethods:FEMIterativeMethods_gcrNumberOfDimension: %d %11.4e\n", k, *residual);
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

-(void)FEMIterativeMethods_richardsonNumberOfDimensions:(int)n matrix:(FEMMatrix *)matrix afterSolve:(double *)xvec rightHandSide:(double *)rhsvec ipar:(int *)ipar rounds:(double)rounds minTolerance:(double)minTolerance maxTolerance:(double)maxTolerance residual:(double *)residual converged:(BOOL *)converged diverged:(BOOL *)diverged outputInterval:(int)outputInterval pcondlMethod:(SEL)pcondlMethod matvecMethod:(SEL)matvecMethod {
    
    int i, j, k;
    double bnorm, rnorm, s;
    double *r, *m;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = matrix.getContainers;
    
    NSMethodSignature *pCondlSignature, *matvecSignature;
    NSInvocation *pcondlInvocation, *matvecInvocation;

    FEMPrecondition *preconditioning;
    
    r = doublevec(0, n-1);
    m = doublevec(0, n-1);
    
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

    [matvecInvocation setArgument:&matrix atIndex:2];
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
        free_dvector(m, 0, n-1);
        return;
    }
    
    // Perform preconditioning by mass lumping
    
    for (i=0; i<n; i++) {
        s = 0.0;
        for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
            s = s + matContainers->Values[j];
        }
        m[i] = s;
    }
    
    for (k=1; k<=rounds; k++) {
        for (i=0; i<n; i++) {
            if (k == 1) {
                xvec[i] = rhsvec[i] / m[i];
            } else {
                xvec[i] = xvec[i] + r[i] / m[i];
            }
        }
        
        [matvecInvocation setArgument:&matrix atIndex:2];
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
            NSLog(@"FEMIterativeMethods:FEMIterativeMethods_richardsonNumberOfDimensions: %d %11.4e %11.4e\n", k, rnorm, *residual);
        }

        *converged = (*residual < minTolerance) ? YES : NO;
        *diverged = (*residual > maxTolerance) ? YES : NO;
        if (*converged == YES || *diverged == YES) break;
    }
    
    free_dvector(r, 0, n-1);
    free_dvector(m, 0, n-1);
}

#pragma mark Public methods

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
        
    }
    
    return self;
}

#pragma mark SGS
/***********************************************************************************
    Symmetric Gauss-Seidel iterative method for linear systems. This is not really
    of practical use but may be used for testing, for example.
***********************************************************************************/
-(void)dsgsSolveMatrix:(FEMMatrix *)matrix ndim:(int)ndim result:(double *)x rhs:(double *)b ipar:(int *)ipar dpar:(double *)dpar pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod {
    
    int rounds, outputInterval;
    double minTol, maxTol, residual, omega;
    BOOL converged = NO, diverged = NO;
    
    rounds = ipar[9];
    minTol = dpar[0];
    maxTol = dpar[1];
    outputInterval = ipar[4];
    omega = dpar[2];
    
    [self FEMIterativeMethods_sgsNumberOfDimensions:ndim matrix:matrix afterSolve:x rightHandSide:b ipar:ipar rounds:rounds minTolerance:minTol maxTolerance:maxTol residual:&residual converged:&converged diverged:&diverged outputInterval:outputInterval omega:omega matvecMethod:matvecMethod];
    
    if (converged == YES) ipar[29] = 1;
    if (diverged == YES) ipar[29] = 3;
    if ( converged == NO && diverged == NO) ipar[29] = 2;
}

#pragma mark JACOBI
/******************************************************************************
    Jacobi iterative method for linear systems. This is not really of practical
    use but may be used for testing, for example.
    Note that if the scaling is performed so that the diagonal entry is one
    the division by it is unnecessary. Hence for this method scaling is not
    needed.
*******************************************************************************/
-(void)djacobiSolveMatrix:(FEMMatrix *)matrix ndim:(int)ndim result:(double *)x rhs:(double *)b ipar:(int *)ipar dpar:(double *)dpar pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod {
    
    int rounds, outputInterval;
    double minTol, maxTol, residual;
    BOOL converged = NO, diverged = NO;
    
    rounds = ipar[9];
    minTol = dpar[0];
    maxTol = dpar[1];
    outputInterval = ipar[4];
    
    [self FEMIterativeMethods_jacobiNumberOfDimensions:ndim matrix:matrix afterSolve:x rightHandSide:b ipar:ipar rounds:rounds minTolerance:minTol maxTolerance:maxTol residual:&residual converged:&converged diverged:&diverged outputInterval:outputInterval matvecMethod:matvecMethod];
    
    if (converged == YES) ipar[29] = 1;
    if (diverged == YES) ipar[29] = 3;
    if ( converged == NO && diverged == NO) ipar[29] = 2;
}

#pragma mark BI-CGSTAB(l)
/***************************************************************************************
    This routine solves real linear systems Ax = b by using the BiCGStab(l) algorithm
    with l >= 2 and the right-oriented ILU(n) preconditioning.
***************************************************************************************/
-(void)dbicgstablSolveMatrix:(FEMMatrix *)matrix ndim:(int)ndim result:(double *)x rhs:(double *)b ipar:(int *)ipar dpar:(double *)dpar pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod {
    
    int nc = 0, rounds, outputInterval, polynomialDegree;
    double minTol, maxTol;
    double *xx, *bb;
    BOOL converged = NO, diverged = NO;
    BOOL constrained;
    
    constrained = (matrix.constraint != Nil) ? YES : NO;
    if (constrained == YES) {
        matrixArraysContainer *matContainers = matrix.constraint.getContainers;
        nc = matrix.constraint.numberOfRows;
        xx = doublevec(0, (ndim+nc)-1);
        bb = doublevec(0, (ndim+nc)-1);
        for (int i=0; i<ndim; i++) {
            xx[i] = x[i];
        }
        for (int i=ndim; i<ndim+nc; i++) {
            xx[i] = 0.0;
        }
        for (int i=0; i<ndim; i++) {
            bb[i] = b[i];
        }
        for (int i=ndim; i<ndim+nc; i++) {
            bb[i] = matContainers->RHS[i];
        }
    }
    
    rounds = ipar[9];
    minTol = dpar[0];
    maxTol = dpar[1];
    outputInterval = ipar[4];
    polynomialDegree = ipar[15];
    
    if (constrained) {
            [self FEMIterativeMethods_BICGStabLNumberOfDimensions:ndim+nc matrix:matrix afterSolve:xx rightHandSide:bb ipar:ipar rounds:rounds minTolerance:minTol maxTolerance:maxTol converged:&converged diverged:&diverged outputInterval:outputInterval polynomialDegree:polynomialDegree constrained:constrained numberOfConstrains:nc pcondlMethod:pcondlMethod matvecMethod:matvecMethod];
    } else {
         [self FEMIterativeMethods_BICGStabLNumberOfDimensions:ndim matrix:matrix afterSolve:x rightHandSide:b ipar:ipar rounds:rounds minTolerance:minTol maxTolerance:maxTol converged:&converged diverged:&diverged outputInterval:outputInterval polynomialDegree:polynomialDegree constrained:constrained numberOfConstrains:nc pcondlMethod:pcondlMethod matvecMethod:matvecMethod];
    }
    
    if (constrained) {
        for (int i=0; i<ndim; i++) {
            x[i] = xx[i];
            b[i] = bb[i];
        }
        free_dvector(xx, 0, (ndim+nc)-1);
        free_dvector(bb, 0, (ndim+nc)-1);
    }
    
    if (converged == YES) ipar[29] = 1;
    if (diverged == YES) ipar[29] = 3;
    if ( converged == NO && diverged == NO) ipar[29] = 2;
}

#pragma mark GCR
/********************************************************************************
    This routine solves real linear systems Ax = b by using the GCR algorithm
    (Generalized Conjugate Residual).
********************************************************************************/
-(void)dgcrSolveMatrix:(FEMMatrix *)matrix ndim:(int)ndim result:(double *)x rhs:(double *)b ipar:(int *)ipar dpar:(double *)dpar pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod {
    
    int rounds, outputInterval, restartN;
    double minTol, maxTol, residual;
    BOOL converged = NO, diverged = NO;
    
    rounds = ipar[9];
    minTol = dpar[0];
    maxTol = dpar[1];
    outputInterval = ipar[4];
    restartN = ipar[16];
    
    [self FEMIterativeMethods_gcrNumberOfDimensions:ndim matrix:matrix afterSolve:x rightHandSide:b ipar:ipar rounds:rounds minTolerance:minTol maxTolerance:maxTol residual:&residual converged:&converged diverged:&diverged outputInterval:outputInterval restart:restartN pcondlMethod:pcondlMethod matvecMethod:matvecMethod];
    
    if (converged == YES) ipar[29] = 1;
    if (diverged == YES) ipar[29] = 3;
    if ( converged == NO && diverged == NO) ipar[29] = 2;
}

#pragma mark Richardson
/***************************************************************************************************
    Richardson iterative method for linear systems. This may of actual use for mass matrices.
    Actually this is not the simple Richardson iteration method as it is preconditioned with 
    the lumped mass matrix. Note that if scaling is performed by the "row equilibrium" method then
    lumped mass is by construction unity (assuming all-positive entries). 
    So for this method scaling is not needed.
****************************************************************************************************/
-(void)drichardsonSolveMatrix:(FEMMatrix *)matrix ndim:(int)ndim result:(double *)x rhs:(double *)b ipar:(int *)ipar dpar:(double *)dpar pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod {
    
    int rounds, outputInterval;
    double minTol, maxTol, residual;
    BOOL converged = NO, diverged = NO;

    rounds = ipar[9];
    minTol = dpar[0];
    maxTol = dpar[1];
    outputInterval = ipar[4];
    
    [self FEMIterativeMethods_richardsonNumberOfDimensions:ndim matrix:matrix afterSolve:x rightHandSide:b ipar:ipar rounds:rounds minTolerance:minTol maxTolerance:maxTol residual:&residual converged:&converged diverged:&diverged outputInterval:outputInterval pcondlMethod:pcondlMethod matvecMethod:matvecMethod];
    
    if (converged == YES) ipar[29] = 1;
    if (diverged == YES) ipar[29] = 3;
    if ( converged == NO && diverged == NO) ipar[29] = 2;
}

@end
