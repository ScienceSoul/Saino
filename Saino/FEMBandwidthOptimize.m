//
//  FEMBandwidthOptimize.m
//  Saino
//
//  Created by Seddik hakime on 08/01/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMBandwidthOptimize.h"
#import "Utils.h"

@interface FEMBandwidthOptimize ()
-(void)FEMBandwidthOptimize_renumber:(ListMatrixEntry_t * __nonnull)current index:(int * __nonnull)index permLocal:(int * __nonnull)permLocal doneIndex:(int * __nonnull)doneIndex localNodes:(int)localNodes;
-(void)FEMBandwidthOptimize_levelSize:(ListMatrix_t * __nonnull)listMatrix maxLevel:(int * __nonnull)maxlevel localNodes:(int)localNodes doneAlready:(BOOL *)doneAlready nin:(int)nin levelin:(int)levelin;
@end

@implementation FEMBandwidthOptimize

#pragma mark Private methods

-(void)FEMBandwidthOptimize_renumber:(ListMatrixEntry_t * __nonnull)current index:(int * __nonnull)index permLocal:(int * __nonnull)permLocal doneIndex:(int * __nonnull)doneIndex localNodes:(int)localNodes {
    
    int k;
    ListMatrixEntry_t *p = NULL;
    
    p = current;
    while (p != NULL) {
        k = p->Index;
        if (k < localNodes) {
            if (doneIndex[k] < 0) {
                permLocal[*index] = k;
                doneIndex[k] = *index;
                *index = *index + 1;
            }
        }
        p = p->Next;
    }
}

-(void)FEMBandwidthOptimize_levelSize:(ListMatrix_t * __nonnull)listMatrix maxLevel:(int * __nonnull)maxlevel localNodes:(int)localNodes doneAlready:(BOOL *)doneAlready nin:(int)nin levelin:(int)levelin {
    
    int n, level, stackp;
    ListMatrixEntry_t *p = NULL;
    Stack_t *stack = NULL;
    
    n = nin;
    level = levelin;
    
    stack = (Stack_t*)malloc(sizeof(Stack_t) * localNodes );
    stackp = -1;
    
    p = listMatrix[n].Head;
    while (p != NULL) {
        stackp++;
        stack[stackp].p = p;
        
        listMatrix[n].Level = level;
        doneAlready[n] = YES;
        *maxlevel = max(*maxlevel, level);
        
        p = listMatrix[n].Head;
        
        while (1) {
            if (p != NULL) {
                n = p->Index;
                if (n < localNodes) {
                    if (doneAlready[n] == NO) {
                        level++;
                        break;
                    }
                }
            } else if (stackp >= 0) {
                p = stack[stackp].p;
                level--;
                stackp--;
            } else {
                break;
            }
            p = p->Next;
        }
    }
    free(stack);
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

/***********************************************************************************
    Method for computing the bandwidth of a sparse matrix
***********************************************************************************/
-(int)computeBandWidthInListMatrix:(ListMatrix_t * __nonnull)list size:(int)n reorder:(int * __nullable)reorder invInitialReorder:(int * __nullable)invInitialReorder {
    
    int i, j, k;
    int halfBandwidth;
    ListMatrixEntry_t *clist = NULL;
    
    halfBandwidth = 0;
    for (i=0; i<n; i++) {
        clist = list[i].Head;
        j = i;
        if (invInitialReorder != NULL) j = invInitialReorder[j];
        while (clist != NULL) {
            k = clist->Index;
            if (invInitialReorder != NULL) k = invInitialReorder[k];
            if (reorder == NULL) {
                halfBandwidth = max(halfBandwidth, abs(j-k));
            } else {
                halfBandwidth = max(halfBandwidth, abs(reorder[j]-reorder[k]));
            }
            clist = clist->Next;
        }
    }
    return halfBandwidth;
}

/***************************************************************************************************************************
    Method for reording variables for bandwidth and/or gaussian elimination filling optimization. Also computes node
    to element connections (which implies node to node connections and thus the global matrix structure).
****************************************************************************************************************************/
-(int)optimizeBandwidthInListMatrix:(ListMatrix_t * __nonnull)listMatrix permutation:(int * __nonnull)perm sizeOfPerm:(int)sizeOfPerm invInitialReorder:(int * __nonnull)invInitialReorder localNodes:(int)localNodes optimize:(BOOL)optimize useOptimized:(BOOL)useOptimized equation:(NSString * __nonnull)equation {
    
    int i, j, k, halfBandwidth, halfBandwidthBefore, halfBandwidthAfter, index;
    int minDegree, startNode, maxLevel;
    int *permLocal, *doneIndex;
    BOOL newRoot, *doneAlready;
    
    fprintf(stdout, "FEMBandwidthOptimize:optimizeBandwidthInListMatrix:---------------------------------------------------------------------------\n");
    fprintf(stdout, "FEMBandwidthOptimize:optimizeBandwidthInListMatrix: computing matrix structure for %s...\n", [equation UTF8String]);
    
    halfBandwidth = [self computeBandWidthInListMatrix:listMatrix size:localNodes reorder:NULL invInitialReorder:NULL] + 1;
    fprintf(stdout, "FEMBandwidthOptimize:optimizeBandwidthInListMatrix: done.\n");
    fprintf(stdout, "FEMBandwidthOptimize:optimizeBandwidthInListMatrix: half bandwidth without optimization: %d.\n", halfBandwidth);
    if (optimize == NO) {
        fprintf(stdout, "FEMBandwidthOptimize:optimizeBandwidthInListMatrix:---------------------------------------------------------------------------\n");
        return halfBandwidth;
    }
    
    halfBandwidthBefore = halfBandwidth;
    fprintf(stdout, "FEMBandwidthOptimize:optimizeBandwidthInListMatrix: bandwidth optimization...\n");
    
    // Search for node to start
    startNode = 0;
    minDegree = listMatrix[startNode].Degree;
    for (i=0; i<localNodes; i++) {
        if (listMatrix[i].Degree < minDegree) {
            startNode = i;
            minDegree = listMatrix[i].Degree;
        }
        listMatrix[i].Level = 0;
    }
    
    doneAlready = (BOOL*)malloc(sizeof(BOOL) * localNodes);
    maxLevel = 0;
    memset( doneAlready, 0, localNodes*sizeof(BOOL) );
    
    [self FEMBandwidthOptimize_levelSize:listMatrix maxLevel:&maxLevel localNodes:localNodes doneAlready:doneAlready nin:startNode levelin:0];
    
    newRoot = YES;
    while (newRoot == YES) {
        newRoot = NO;
        minDegree = listMatrix[startNode].Degree;
        k = startNode;
        
        for (i=0; i<localNodes; i++) {
            if (listMatrix[i].Level == maxLevel) {
                if (listMatrix[i].Degree < minDegree) {
                    k = i;
                    minDegree = listMatrix[i].Degree;
                }
            }
        }
    
        if (k != startNode) {
            j = maxLevel;
            maxLevel = 0;
            memset( doneAlready, 0, localNodes*sizeof(BOOL) );            
            [self FEMBandwidthOptimize_levelSize:listMatrix maxLevel:&maxLevel localNodes:localNodes doneAlready:doneAlready nin:k levelin:0];
            
            if (j > maxLevel) {
                newRoot = YES;
                startNode = j;
            }
        }
    }
    
    permLocal = intvec(0, sizeOfPerm-1);
    doneIndex = intvec(0, localNodes-1);
    memset( permLocal, -1, sizeOfPerm*sizeof(int) );
    memset( doneIndex, -1, localNodes*sizeof(int) );
    
    // This loop really does the thing
    index = 0;
    permLocal[index] = startNode;
    doneIndex[startNode] = index;
    index++;
    
    for (i=0; i<localNodes; i++) {
        if (permLocal[i] < 0) {
            for (j=0; j<localNodes; j++) {
                if (doneIndex[j] < 0) {
                    permLocal[index] = j;
                    doneIndex[j] = index;
                    index++;
                }
            }
        }
        [self FEMBandwidthOptimize_renumber:listMatrix[permLocal[i]].Head index:&index permLocal:permLocal doneIndex:doneIndex localNodes:localNodes];
    }
    
    // Store it the other way around for FEM and reverse order for profile optimization
    memset( doneIndex, -1, localNodes*sizeof(int) );
    for (i=0; i<localNodes; i++) {
        doneIndex[permLocal[i]] = localNodes-i-1;
    }
    
    memcpy(permLocal, perm, sizeOfPerm*sizeof(int));
    memset( perm, -1, sizeOfPerm*sizeof(int) );
    for (i=0; i<sizeOfPerm; i++) {
        k = permLocal[i];
        if (k >= 0) perm[i] = doneIndex[k];
    }
    free_ivector(doneIndex, 0, localNodes-1);

    halfBandwidthAfter = [self computeBandWidthInListMatrix:listMatrix size:localNodes reorder:perm invInitialReorder:invInitialReorder] + 1;
    fprintf(stdout, "FEMBandwidthOptimize:optimizeBandwidthInListMatrix: ...done.\n");
    fprintf(stdout, "FEMBandwidthOptimize:optimizeBandwidthInListMatrix: half bandwidth after optimization: %d.\n", halfBandwidthAfter);
    halfBandwidth = halfBandwidthAfter;
    if (halfBandwidthBefore < halfBandwidth && useOptimized == NO) {
        fprintf(stdout, "FEMBandwidthOptimize:optimizeBandwidthInListMatrix: optimization rejected, using original ordering.");
        halfBandwidth = halfBandwidthBefore;
        memcpy(perm, permLocal, sizeOfPerm*sizeof(int));
    }
    fprintf(stdout, "FEMBandwidthOptimize:optimizeBandwidthInListMatrix:---------------------------------------------------------------------------\n");

    free_ivector(permLocal, 0, sizeOfPerm-1);
    free(doneAlready);

    return halfBandwidth;
}


@end
