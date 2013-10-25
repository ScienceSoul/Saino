//
//  FEMPost.m
//  Saino
//
//  Created by Seddik hakime on 29/03/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#include <time.h>  

#import "FEMPost.h"
#import "FEMUtilities.h"
#import "FEMBoundaryCondition.h"
#import "FEMListUtilities.h"
#import "FileReader.h"
#import "Utils.h"

@interface FEMPost ()
-(void)FEMPost_writeString:(NSString *)string toFileHandle:(NSFileHandle *)fileHandle;
-(void)FEMPost_writeInteger:(int)number toFileHandle:(NSFileHandle *)fileHandle;
-(void)FEMPost_writeDouble:(double)number toFileHandle:(NSFileHandle *)fileHandle;
-(void)FEMPost_writeBytes:(const void *)bytes length:(int)length toFileHandle:(NSFileHandle *)fileHandle;
@end

@implementation FEMPost

#pragma mark Private methods

-(void)FEMPost_writeString:(NSString *)string toFileHandle:(NSFileHandle *)fileHandle {
    
    NSData *buffer;
    buffer = (NSMutableData *)[string dataUsingEncoding:NSUTF8StringEncoding];
    [fileHandle writeData:buffer];
}

-(void)FEMPost_writeInteger:(int)number toFileHandle:(NSFileHandle *)fileHandle {
    
    NSData *buffer;
    NSString *strBuffer;
    strBuffer = [NSString stringWithFormat:@"%d", number];
    buffer = (NSMutableData *)[strBuffer dataUsingEncoding:NSUTF8StringEncoding];
    [fileHandle writeData:buffer];
}

-(void)FEMPost_writeDouble:(double)number toFileHandle:(NSFileHandle *)fileHandle {
    
    NSData *buffer;
    NSString *strBuffer;
    strBuffer = [NSString stringWithFormat:@"%17.8e", number];
    buffer = (NSMutableData *)[strBuffer dataUsingEncoding:NSUTF8StringEncoding];
    [fileHandle writeData:buffer];
}

-(void)FEMPost_writeBytes:(const void *)bytes length:(int)length toFileHandle:(NSFileHandle *)fileHandle {
    
    NSData *buffer;
    buffer = [NSData dataWithBytes:bytes length:length];
    [fileHandle writeData:buffer];
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

/*********************************************************
    Writes data in ElmerPost format
*********************************************************/
-(void)writeElmerPostFile:(NSString *)postFile resultFile:(NSString *)resultFile model:(FEMModel *)model timeCount:(int)timeCount append:(BOOL *)append {
    
    int i, ii, j, jj, k, l, n, q, numberOfNodes, numberOfElements, *maskOrder, dofs, nDofs, meshDim, index, savedCount, timeStep, node, iDummy, nZeros;
    double meshScale, coord[3], time, dummy;
    BOOL found, freeSurfaceFlag, moveBoundary, maskExists=NO, all, onlySearch;
    NSRange ind, ind1;
    NSString *str, *line;
    NSMutableString *varName, *bodyName;
    NSMutableString *outputPath;
    NSArray *stringParts, *filteredArray;
    FEMMesh *mesh;
    FEMVariable *maskVar = nil, *displacement = nil, *meshUpdate = nil, *var = nil;
    FEMListUtilities *listUtiltiies;
    FEMUtilities *utilities;
    FEMBoundaryCondition *boundaryConditionAtId;
    NSFileManager *fileManager;
    NSFileHandle *postFileHandle;
    FileReader *reader;
    Element_t *elements = NULL;
    Nodes_t *nodes = NULL;
    variableArraysContainer *varContainers = NULL, *var1Containers = NULL, *maskVarContainers = NULL, *displaceVarContainers = NULL;
    
    listUtiltiies = [[FEMListUtilities alloc] init];
    utilities = [[FEMUtilities alloc] init];
    
    mesh = (FEMMesh *)model.mesh;
    if (mesh.savesDone == 0) {
        NSLog(@"FEMPost:writeElmerPostFile: saving results in ElmerPost format to file %@\n", postFile);
    }
    
    fileManager = [NSFileManager defaultManager];
    
    if ([utilities isFileNameQualified:postFile] == NO) {
        if ([model.outputPath length] > 0) {
            outputPath = [NSMutableString stringWithString:model.outputPath];
            [outputPath appendString:@"/"];
            [outputPath appendString:postFile];
            if ([fileManager fileExistsAtPath:outputPath] == NO) {
                if ([fileManager createFileAtPath:outputPath contents:nil attributes:nil] == YES) {
                    postFileHandle = [NSFileHandle fileHandleForWritingAtPath:outputPath];
                } else {
                    errorfunct("FEMPost:writeElmerPostFile", "Can't create post file.");
                }
            } else {
                if (*append == YES && mesh.savesDone != 0) {
                    postFileHandle = [NSFileHandle fileHandleForWritingAtPath:outputPath];
                    [postFileHandle seekToEndOfFile];
                } else { // File already exists, erase it and start from scratch
                    [fileManager removeItemAtPath:outputPath error:nil];
                    if ([fileManager createFileAtPath:outputPath contents:nil attributes:nil] == YES) {
                        postFileHandle = [NSFileHandle fileHandleForWritingAtPath:outputPath];
                    } else {
                        errorfunct("FEMPost:writeElmerPostFile", "Can't create post file.");
                    }
                }
            }
        } else {
            if ([fileManager fileExistsAtPath:postFile] == NO) {
                if ([fileManager createFileAtPath:postFile contents:nil attributes:nil] == YES) {
                    postFileHandle = [NSFileHandle fileHandleForWritingAtPath:postFile];
                } else {
                    errorfunct("FEMPost:writeElmerPostFile", "Can't create post file.");
                }
            } else {
                 if (*append == YES && mesh.savesDone != 0) {
                     postFileHandle = [NSFileHandle fileHandleForWritingAtPath:postFile];
                     [postFileHandle seekToEndOfFile];
                 } else { // File already exists, erase it and start from scratch
                     [fileManager removeItemAtPath:outputPath error:nil];
                     if ([fileManager createFileAtPath:outputPath contents:nil attributes:nil] == YES) {
                         postFileHandle = [NSFileHandle fileHandleForWritingAtPath:outputPath];
                     } else {
                         errorfunct("FEMPost:writeElmerPostFile", "Can't create post file.");
                     }
                 }
            }
        }
    } else {
        if ([fileManager fileExistsAtPath:postFile] == NO) {
            if ([fileManager createFileAtPath:postFile contents:nil attributes:nil] == YES) {
                postFileHandle = [NSFileHandle fileHandleForWritingAtPath:postFile];
            } else {
                errorfunct("FEMPost:writeElmerPostFile", "Can't create post file.");
            }
        } else {
            if (*append == YES && mesh.savesDone != 0) {
                postFileHandle = [NSFileHandle fileHandleForWritingAtPath:postFile];
                [postFileHandle seekToEndOfFile];
            } else { // File already exists, erase it and start from scratch
                [fileManager removeItemAtPath:outputPath error:nil];
                if ([fileManager createFileAtPath:outputPath contents:nil attributes:nil] == YES) {
                    postFileHandle = [NSFileHandle fileHandleForWritingAtPath:outputPath];
                } else {
                    errorfunct("FEMPost:writeElmerPostFile", "Can't create post file.");
                }
            }
        }
    }
    
    if (*append == NO) {
        if ([utilities isFileNameQualified:resultFile] == NO) {
            if ([model.outputPath length] > 0) {
                outputPath = [NSMutableString stringWithString:model.outputPath];
                [outputPath appendString:@"/"];
                [outputPath appendString:resultFile];
                if ([fileManager fileExistsAtPath:outputPath] == NO) {
                    errorfunct("FEMPost:writeElmerPostFile", "Result file does not exist.");
                } else {
                    reader = [[FileReader alloc] initWithFilePath:outputPath];
                }
            } else {
                if ([fileManager fileExistsAtPath:resultFile] == NO) {
                    errorfunct("FEMPost:writeElmerPostFile", "Result file does not exist.");
                } else {
                    reader = [[FileReader alloc] initWithFilePath:resultFile];
                }
            }
        } else {
            if ([fileManager fileExistsAtPath:resultFile] == NO) {
                errorfunct("FEMPost:writeElmerPostFile", "Result file does not exist.");
            } else {
                reader = [[FileReader alloc] initWithFilePath:resultFile];
            }
        }
    }
    
    freeSurfaceFlag = NO;
    moveBoundary = NO;
    for (FEMBoundaryCondition *boundaryCondition in model.boundaryConditions) {
        freeSurfaceFlag = (freeSurfaceFlag == YES || [listUtiltiies listGetLogical:model inArray:boundaryCondition.valuesList forVariable:@"free surface" info:&found] == YES) ? YES : NO;
        if (freeSurfaceFlag == YES) {
            moveBoundary = [listUtiltiies listGetLogical:model inArray:boundaryCondition.valuesList forVariable:@"internal move boundary" info:&found];
            if (found == NO) moveBoundary = YES;
            freeSurfaceFlag = (freeSurfaceFlag == YES && moveBoundary == YES) ? YES : NO;
        }
        if (freeSurfaceFlag == YES) break;
    }
    
    // Initialize stuff for masked saving
    str = [listUtiltiies listGetString:model inArray:model.simulation.valuesList forVariable:@"post file mask variable" info:&found];
    if (found == NO) {
        str = [listUtiltiies listGetString:model inArray:model.simulation.valuesList forVariable:@"elmerpost mask variable" info:&found];
    }
    
    if (found == YES) {
        maskVar = [utilities getVariableFrom:model.variables model:model name:str onlySearch:NULL maskName:NULL info:&found];
        if (maskVar != nil) {
            maskVarContainers = maskVar.getContainers;
            maskExists = (maskVarContainers->Perm != NULL) ? YES : NO;
        }
    }
    if (maskExists == YES) {
        NSLog(@"FEMPost:writeElmerPostFile: using %@ as mask variable\n", str);
        numberOfNodes = max_array(maskVarContainers->Perm, maskVarContainers->sizePerm);
        maskOrder = intvec(0, numberOfNodes);
        for (i=0; maskVarContainers->sizePerm; i++) {
            j = maskVarContainers->Perm[i];
            if (j >= 0) maskOrder[j] = i;
        }
        numberOfElements = 0;
        elements = model.getElements;
        for (i=0; i<model.numberOfBulkElements+model.numberOfBoundaryElements; i++) {
            all = YES;
            for (j=0; j<elements[i].Type.NumberOfNodes; j++) {
                if (maskVarContainers->Perm[elements[i].NodeIndexes[j]] == -1) {
                    all = NO;
                    break;
                }
            }
            if (all == YES) numberOfElements++;
        }
    } else {
        numberOfNodes = model.numberOfNodes;
        numberOfElements = model.numberOfBulkElements + model.numberOfBoundaryElements;
    }
    
    // Count degrees of freedom to be saved
    dofs = 0;
    for (FEMVariable *variable in model.variables) {
        if (variable.output == NO) continue;
        
        varContainers = variable.getContainers;
        if (varContainers->sizeValues == variable.dofs) continue;
        
        if (variable.type != VARIABLE_ON_NODES) continue;
        
        if ([variable.name isEqualToString:@"mesh update"] == YES) {
            found = NO;
            for (FEMVariable *var1 in model.variables) {
                if ([var1.name isEqualToString:@"displacement"] == YES) {
                    found = YES;
                    break;
                }
            }
            if (found == NO) dofs = dofs + 3;
        } else if ([variable.name isEqualToString:@"mesh update 1"] == YES || [variable.name isEqualToString:@"mesh update 2"] == YES ||
                   [variable.name isEqualToString:@"mesh update 3"] == YES) {
            // No ops
        } else if ([variable.name isEqualToString:@"displacement"] == YES) {
            dofs = dofs + 3;
            if (varContainers->CValues != NULL) dofs = dofs + 3;
        } else if ([variable.name isEqualToString:@"displacement 1"] == YES || [variable.name isEqualToString:@"displacement 2"] == YES ||
                   [variable.name isEqualToString:@"displacement 3"] == YES) {
            // No ops
        } else if ([variable.name isEqualToString:@"flow solution"] == YES) {
            dofs = dofs + 4;
        } else if ([variable.name isEqualToString:@"velocity 1"] == YES || [variable.name isEqualToString:@"velocity 2"] == YES ||
                   [variable.name isEqualToString:@"velocity 3"] == YES) {
            // No ops
        } else if ([variable.name isEqualToString:@"magnetic field"] == YES) {
            dofs = dofs + 3;
        } else if ([variable.name isEqualToString:@"magnetic field 1"] == YES || [variable.name isEqualToString:@"magnetic field 2"] == YES ||
                   [variable.name isEqualToString:@"magnetic field 3"] == YES) {
            // No ops
        } else {
            nDofs = 1;
            if (varContainers->CValues != NULL) nDofs = 2;
            
            if (variable.dofs == 1) {
                dofs = dofs + nDofs;
            } else {
                ind = [variable.name rangeOfString:@"["];
                if (ind.location != NSNotFound) {
                    j = 0;
                    while (1) {
                        str = [variable.name substringFromIndex:j];
                        if (j == 0) {
                            ind = [variable.name rangeOfString:@"["];
                        } else {
                            ind.location = -1;
                        }
                        ind1 = [str rangeOfString:@":"];
                        if (ind1.location == NSNotFound) errorfunct("FEMPost:writeElmerPostFile", "Missing separator ':' in variable definition using '[ ]' syntax.");
                        if (j == 0) {
                            if (ind1.location < ind.location) errorfunct("FEMPost:writeElmerPostFile", "Syntax error in variable definition.");
                        } else {
                            if (ind1.location == 0) errorfunct("FEMPost:writeElmerPostFile", "Syntax error in variable definition.");
                        }
                        i = (int)ind1.location + 1;
                        while (1) {
                            if ([str characterAtIndex:i] == ' ' || [str characterAtIndex:i] == ']') break;
                            i++;
                        }
                        k = [[[str substringFromIndex:ind1.location+1] substringToIndex:i-(ind1.location+1)] intValue];
                        dofs = dofs + nDofs*k;
                        if (k == 2) dofs = dofs + nDofs;
                        while ([str characterAtIndex:i] == ' ') {
                            i++;
                        }
                        j = i;
                    }
                }
            }
        }
    }
    
    if (freeSurfaceFlag == NO) dofs = dofs - 3;
    
    char space = ' ';
    char newLine = '\n';
    NSData *spaceBuff = [NSMutableData dataWithBytes:&space length:sizeof(space)];
    NSData *newLineBuff = [NSMutableData dataWithBytes:&newLine length:sizeof(newLine)];
    
    // Write header to output
    if (append == NO || mesh.savesDone == 0) {

        [self FEMPost_writeInteger:numberOfNodes toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
        [self FEMPost_writeInteger:numberOfElements toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
        [self FEMPost_writeInteger:dofs toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
        [self FEMPost_writeInteger:timeCount toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
        
        for (FEMVariable *variable in model.variables) {
            if (variable.output == NO) continue;
            
            varContainers = variable.getContainers;
            if (varContainers->sizeValues == variable.dofs) continue;
            
            if (variable.type != VARIABLE_ON_NODES) continue;
            
            if ([variable.name isEqualToString:@"mesh update"] == YES) {
                found = NO;
                for (FEMVariable *var1 in model.variables) {
                    if ([var1.name isEqualToString:@"displacement"] == YES) {
                        found = YES;
                        break;
                    }
                }
                if (found == NO) {
                    [self FEMPost_writeString: @"vector: Mesh.update" toFileHandle:postFileHandle];
                    displacement = variable;
                } else {
                    meshUpdate = variable;
                }
            } else if ([variable.name isEqualToString:@"mesh update 1"] == YES || [variable.name isEqualToString:@"mesh update 2"] == YES ||
                       [variable.name isEqualToString:@"mesh update 3"] == YES) {
                // No ops
            } else if ([variable.name isEqualToString:@"displacement"] == YES) {
                [self FEMPost_writeString:@" vector: Displacement" toFileHandle:postFileHandle];
                if (varContainers->CValues != NULL) {
                    [self FEMPost_writeString:@" vector: Displacement.im" toFileHandle:postFileHandle];
                }
                displacement = variable;
            } else if ([variable.name isEqualToString:@"displacement 1"] == YES || [variable.name isEqualToString:@"displacement 2"] == YES ||
                       [variable.name isEqualToString:@"displacement 3"] == YES) {
                // No ops
            } else if ([variable.name isEqualToString:@"flow solution"] == YES) {
                [self FEMPost_writeString:@" vector: Velocity scalar: Pressure" toFileHandle:postFileHandle];
            } else if ([variable.name isEqualToString:@"velocity 1"] == YES || [variable.name isEqualToString:@"velocity 2"] == YES ||
                       [variable.name isEqualToString:@"velocity 3"] == YES || [variable.name isEqualToString:@"pressure"] == YES) {
                // No ops
            } else if ([variable.name isEqualToString:@"magnetic field"] == YES) {
                [self FEMPost_writeString:@" vector: Magfield" toFileHandle:postFileHandle];
            } else if ([variable.name isEqualToString:@"magnetic field 1"] == YES || [variable.name isEqualToString:@"magnetic field 2"] == YES ||
                       [variable.name isEqualToString:@"magnetic field 3"] == YES) {
                // No ops
            } else if ([variable.name isEqualToString:@"coordinate 1"] == YES || [variable.name isEqualToString:@"coordinate 2"] == YES ||
                       [variable.name isEqualToString:@"coordinate 3"] == YES) {
                // No ops
            } else {
                if (variable.dofs == 1) {
                    varName = [NSMutableString stringWithString:variable.name];
                    [varName replaceOccurrencesOfString:@" " withString:@"." options:NSLiteralSearch range:NSMakeRange(0, [varName length])];
                    char character =  (int)[varName characterAtIndex:0] - (int)'a' + (int)'A';
                    [varName replaceCharactersInRange:NSMakeRange(0, 1) withString:[NSString stringWithUTF8String:&character]];
                    [self FEMPost_writeString:@" scalar: " toFileHandle:postFileHandle];
                    [self FEMPost_writeString:varName toFileHandle:postFileHandle];
                    if (varContainers->CValues != NULL) {
                        [self FEMPost_writeString:@" scalar: " toFileHandle:postFileHandle];
                        [self FEMPost_writeString:[varName stringByAppendingString:@".im"] toFileHandle:postFileHandle];
                    }
                } else {
                    ind = [variable.name rangeOfString:@"["];
                    if (ind.location != NSNotFound) {
                        j = 0;
                        while (1) {
                            str = [variable.name substringFromIndex:j];
                            if (j == 0) {
                                ind = [variable.name rangeOfString:@"["];
                            } else {
                                ind.location = -1;
                            }
                            ind1 = [str rangeOfString:@":"];
                            if (ind1.location == NSNotFound) errorfunct("FEMPost:writeElmerPostFile", "Missing separator ':' in variable definition using '[ ]' syntax.");
                            if (j == 0) {
                                if (ind1.location < ind.location) errorfunct("FEMPost:writeElmerPostFile", "Syntax error in variable definition.");
                            } else {
                                if (ind1.location == 0) errorfunct("FEMPost:writeElmerPostFile", "Syntax error in variable definition.");
                            }
                            i = (int)ind1.location + 1;
                            while (1) {
                                if ([str characterAtIndex:i] == ' ' || [str characterAtIndex:i] == ']') break;
                                i++;
                            }
                            k = [[[str substringFromIndex:ind1.location+1] substringToIndex:i-(ind1.location+1)] intValue];
                            dofs = dofs + k;
                            varName = (NSMutableString *)[str substringToIndex:ind1.location];
                            [varName replaceOccurrencesOfString:@" " withString:@"." options:NSLiteralSearch range:NSMakeRange(0, [varName length])];
                            char character = (int)[varName characterAtIndex:0] - (int)'a' + (int)'A';
                            [varName replaceCharactersInRange:NSMakeRange(0, 1) withString:[NSString stringWithUTF8String:&character]];
                            if (k == 1) {
                                [self FEMPost_writeString:@" scalar: " toFileHandle:postFileHandle];
                                [self FEMPost_writeString:varName toFileHandle:postFileHandle];
                                if (varContainers->CValues != NULL) {
                                    [self FEMPost_writeString:@" scalar: " toFileHandle:postFileHandle];
                                    [self FEMPost_writeString:[varName stringByAppendingString:@".im"] toFileHandle:postFileHandle];
                                }
                            } else if (k <= 3) {
                                [self FEMPost_writeString:@" vector: " toFileHandle:postFileHandle];
                                [self FEMPost_writeString:varName toFileHandle:postFileHandle];
                                if (varContainers->CValues != NULL) {
                                    [self FEMPost_writeString:@" vector: " toFileHandle:postFileHandle];
                                    [self FEMPost_writeString:[varName stringByAppendingString:@".im"] toFileHandle:postFileHandle];
                                }
                            } else {
                                for (l=1; l<=k; l++) {
                                    [self FEMPost_writeString:@" scalar: " toFileHandle:postFileHandle];
                                    NSMutableString *componentName = [NSMutableString stringWithString:varName];
                                    [componentName appendString:@"."];
                                    [self FEMPost_writeString:[componentName stringByAppendingString:[NSString stringWithFormat:@"%d",l]] toFileHandle:postFileHandle];
                                    if (varContainers->CValues != NULL) {
                                        [self FEMPost_writeString:@" scalar: " toFileHandle:postFileHandle];
                                        NSMutableString *componentName = [NSMutableString stringWithString:varName];
                                        [componentName appendString:@"."];
                                        [componentName appendString:[NSString stringWithFormat:@"%d",l]];
                                        [self FEMPost_writeString:[componentName stringByAppendingString:@".im"] toFileHandle:postFileHandle];
                                    }
                                }
                            }
                            while ([str characterAtIndex:i] == ' ') {
                                i++;
                            }
                            j = i;
                        }
                    }
                }
            }
        }
        
        if (freeSurfaceFlag == YES) {
            [self FEMPost_writeString:@" vector: Coordinates" toFileHandle:postFileHandle];
        }
        
        [postFileHandle writeData:newLineBuff];
        
        [self FEMPost_writeString:@"#File started at: " toFileHandle:postFileHandle];
        NSString *dateString = [NSString stringWithCString:dateAndTime() encoding:NSASCIIStringEncoding];
        [self FEMPost_writeString:dateString toFileHandle:postFileHandle];
        [postFileHandle writeData:newLineBuff];
        
        // Coordinates
        meshScale = 1.0;
        for (FEMSolution *solution in model.solutions) {
            if ((solution.solutionInfo)[@"displace mesh"] != nil) {
                if ([(solution.solutionInfo)[@"displace mesh"] boolValue] == NO) meshScale = 0.0;
            } else {
                if ((solution.solutionInfo)[@"output mesh deformation"] != nil) {
                    if ([(solution.solutionInfo)[@"output mesh deformation"] boolValue] == YES) meshScale = 0.0;
                } else {
                    if (solution.nOfEigenValues > 0) meshScale = 0.0;
                }
            }
        }
        
        meshDim = mesh.dimension;
        nodes = model.getNodes;
        for (ii=0; ii<numberOfNodes; ii++) {
            
            i = ii;
            if (maskExists == YES) i = maskOrder[i];
            
            coord[0] = nodes->x[i];
            coord[1] = nodes->y[i];
            coord[2] = nodes->z[i];
            
            if (displacement != nil) {
                displaceVarContainers = displacement.getContainers;
                k = displaceVarContainers->Perm[i];
                
                if (k >= 0) {
                    for (l=0; l<displacement.dofs; l++) {
                        coord[l] = coord[l] - meshScale*displaceVarContainers->Values[displacement.dofs * k + l];
                    }
                } else if (meshUpdate != nil) {
                    displaceVarContainers = meshUpdate.getContainers;
                    k = displaceVarContainers->Perm[i];
                    if (k >= 0) {
                        for (l=0; l<meshUpdate.dofs; l++) {
                            coord[l] = coord[l] - meshScale*displaceVarContainers->Values[meshUpdate.dofs * k + l];
                        }
                    }
                }
            }
            
            if (meshDim == 3) {
                for (j=0; j<meshDim; j++) {
                    [self FEMPost_writeDouble:coord[j] toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                }
            } else if (meshDim == 2) {
                for (j=0; j<meshDim; j++) {
                    [self FEMPost_writeDouble:coord[j] toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                }
                [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle];
            } else {
                [self FEMPost_writeDouble:coord[0] toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle];
            }
            [postFileHandle writeData:newLineBuff];
        }
        
        // Elements
        [self FEMPost_writeString:@"#group all" toFileHandle:postFileHandle];
        [postFileHandle writeData:newLineBuff];
        
        elements = model.getElements;
        for (i=0; i<model.numberOfBulkElements; i++) {
            
            if (maskExists == YES) {
                all = YES;
                for (j=0; j<elements[i].Type.NumberOfNodes; j++) {
                    if (maskVarContainers->Perm[elements[i].NodeIndexes[j]] < 0) {
                        all = NO;
                        break;
                    }
                }
                if (all == NO) continue;
            }
            
            k = elements[i].BodyID;
            if (k >= 1 && k <= model.numberOfBodies) {
                if ((model.bodies)[k-1][@"name"] != nil) {
                    bodyName = [NSMutableString stringWithString:(model.bodies)[k-1][@"name"]];
                    [bodyName replaceOccurrencesOfString:@" " withString:@"." options:NSLiteralSearch range:NSMakeRange(0, [bodyName length])];
                    [self FEMPost_writeString:bodyName toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                } else {
                    bodyName = [NSMutableString stringWithString:@"body"];
                    [bodyName appendString:[NSString stringWithFormat:@"%d", k]];
                    [self FEMPost_writeString:bodyName toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                }
            }
            
            [self FEMPost_writeInteger:elements[i].Type.ElementCode toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
            n = 0;
            for (j=0; j<elements[i].Type.NumberOfNodes; j+=4) {
                l = min(4, elements[i].Type.NumberOfNodes-n);
                for (k=0; k<l; k++) {
                    index = elements[i].NodeIndexes[n];
                    if (maskExists == YES) index = maskVarContainers->Perm[index];
                    [self FEMPost_writeInteger:index toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                    n++;
                }
                [postFileHandle writeData:newLineBuff];
            }
        }
        
        for (i=model.numberOfBulkElements; i<model.numberOfBulkElements+model.numberOfBoundaryElements; i++) {
            
            if (maskExists == YES) {
                all = YES;
                for (j=0; j<elements[i].Type.NumberOfNodes; j++) {
                    if (maskVarContainers->Perm[elements[i].NodeIndexes[j]] < 0) {
                        all = NO;
                        break;
                    }
                }
                if (all == NO) continue;
            }
            
            k = elements[i].BoundaryInfo->Constraint;
            if (k >= 1 && k <= model.numberOfBoundaryConditions) {
                boundaryConditionAtId = model.boundaryConditions[k-1];
                bodyName = [NSMutableString stringWithString:[listUtiltiies listGetString:model inArray:boundaryConditionAtId.valuesList forVariable:@"name" info:&found]];
            }
            if (found == YES) {
                [bodyName replaceOccurrencesOfString:@" " withString:@"." options:NSLiteralSearch range:NSMakeRange(0, [bodyName length])];
                [self FEMPost_writeString:bodyName toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
            } else {
                bodyName = [NSMutableString stringWithString:@"Constraint"];
                [bodyName appendString:[NSString stringWithFormat:@"%d", k]];
                [self FEMPost_writeString:bodyName toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
            }
            [self FEMPost_writeInteger: elements[i].Type.ElementCode toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
            for (k=0; k<elements[i].Type.NumberOfNodes; k++) {
                index = elements[i].NodeIndexes[k];
                if (maskExists == YES) index = maskVarContainers->Perm[index];
                [self FEMPost_writeInteger:index toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
            }
            [postFileHandle writeData:newLineBuff];
        }
        
        [self FEMPost_writeString:@"#endgroup all" toFileHandle:postFileHandle];
        [postFileHandle writeData:newLineBuff];
    }
    
    if (*append == YES && mesh.savesDone == 0) {
        [postFileHandle closeFile];
        return;
    }
    
    // Used to separate strings and filter them from white spaces
    NSCharacterSet *whitespaces = [NSCharacterSet whitespaceCharacterSet];
    NSPredicate *noEmptyStrings = [NSPredicate predicateWithFormat:@"SELF != ''"];
    
    line = nil;
    while (1) {
        if (*append == YES) {
            savedCount = mesh.savesDone;
            timeStep = savedCount;
            onlySearch = YES;
            var = [utilities getVariableFrom:model.variables model:model name:@"time" onlySearch:&onlySearch maskName:NULL info:&found];
            time = 1.0;
            if (var != nil) {
                varContainers = var.getContainers;
                time = varContainers->Values[0];
            }
        } else {
            // Read one time step to memory (if not already there)...
            while ((line = [reader readLine])) {
                stringParts = [line componentsSeparatedByCharactersInSet:whitespaces];
                filteredArray = [stringParts filteredArrayUsingPredicate:noEmptyStrings];
                if ([filteredArray[0] isEqualToString:@"total"] == YES && [filteredArray[1] isEqualToString:@"dofs:"] == YES) {
                    dofs = [filteredArray[2] intValue];
                }
                if ([filteredArray[0] isEqualToString:@"time:"] == YES) break;
            }
            
            if ([filteredArray[0] isEqualToString:@"time:"] == NO) break;
            savedCount = [filteredArray[1] intValue];
            timeStep = [filteredArray[2] intValue];
            time = [filteredArray[3] doubleValue];
        }
        
        [self FEMPost_writeString:@"#time " toFileHandle:postFileHandle];
        [self FEMPost_writeInteger:savedCount toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
        [self FEMPost_writeInteger:timeStep toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
        [self FEMPost_writeDouble:time toFileHandle:postFileHandle];
        [postFileHandle writeData:newLineBuff];
        
        if (*append == NO) {
            for (i=0; i<dofs; i++) {
                line = [reader readLine];
                stringParts = [line componentsSeparatedByCharactersInSet:whitespaces];
                filteredArray = [stringParts filteredArrayUsingPredicate:noEmptyStrings];
                onlySearch = YES;
                var = [utilities getVariableFrom:model.variables model:model name:filteredArray[0] onlySearch:&onlySearch maskName:NULL info:&found];
                if (var != nil) {
                    varContainers = var.getContainers;
                    for (j=0; j<numberOfNodes; j++) {
                        k = j;
                        if (maskExists == YES) k = maskVarContainers->Perm[k];
                        if (varContainers->Perm != NULL) k = varContainers->Perm[k];
                        line = [reader readLine];
                        stringParts = [line componentsSeparatedByCharactersInSet:whitespaces];
                        filteredArray = [stringParts filteredArrayUsingPredicate:noEmptyStrings];
                        if (k >= 0) {
                            node = [filteredArray[0] intValue];
                            iDummy = [filteredArray[1] intValue];
                            varContainers->Values[k] = [filteredArray[2] doubleValue];
                        } else {
                            node = [filteredArray[0] intValue];
                            iDummy = [filteredArray[1] intValue];
                            dummy = [filteredArray[2] doubleValue];
                        }
                    }
                }
            }
        }
        
        //... then save it to post file
        for (ii=0; ii<numberOfNodes; ii++) {
            
            i = ii;
            if (maskExists == YES) i = maskOrder[i];
            
            for (FEMVariable *variable in model.variables) {
                if (variable.output == NO) continue;
                varContainers = variable.getContainers;
                if (varContainers->sizeValues == variable.dofs) continue;
                if (variable.type != VARIABLE_ON_NODES) continue;
                
                if ([variable.name isEqualToString:@"mesh update"] == YES) {
                    found = NO;
                    for (FEMVariable *var1 in model.variables) {
                        if ([var1 isEqualTo:@"displacement"] == YES) {
                            found = YES;
                            break;
                        }
                    }
                    if (found == NO) {
                        k = i;
                        if (varContainers->Perm != NULL) k = varContainers->Perm[k];
                        if (k >= 0) {
                            for (j=0; j<variable.dofs; j++) {
                                [self FEMPost_writeDouble:varContainers->Values[variable.dofs*k+j] toFileHandle:postFileHandle];
                                [postFileHandle writeData:spaceBuff];
                            }
                            if (variable.dofs == 2) {
                                [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];

                            }
                        } else {
                            [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                            [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                            [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                        }
                    }
                } else if ([variable.name isEqualToString:@"mesh update 1"] == YES || [variable.name isEqualToString:@"mesh update 2"] == YES ||
                           [variable.name isEqualToString:@"mesh update 3"] == YES) {
                    // No ops
                } else if ([variable.name isEqualToString:@"displacement"] == YES) {
                    k = i;
                    if (varContainers->Perm != NULL) k = varContainers->Perm[k];
                    if (k >= 0) {
                        if (varContainers->CValues != NULL) {
                            for (j=0; j<variable.dofs; j++) {
                                [self FEMPost_writeDouble:creal(*(varContainers->CValues[variable.dofs*k+j])) toFileHandle:postFileHandle];
                                [postFileHandle writeData:spaceBuff];
                            }
                            if (variable.dofs == 2) {
                                [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                            }
                            for (j=0; j<variable.dofs; j++) {
                                [self FEMPost_writeDouble:cimag(*(varContainers->CValues[variable.dofs*k+j])) toFileHandle:postFileHandle];
                                [postFileHandle writeData:spaceBuff];
                            }
                            if (variable.dofs == 2) {
                                [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                            }
                        } else {
                            for (j=0; j<variable.dofs; j++) {
                                [self FEMPost_writeDouble:varContainers->Values[variable.dofs*k+j] toFileHandle:postFileHandle];
                                [postFileHandle writeData:spaceBuff];
                            }
                            if (variable.dofs == 2) {
                                [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                            }
                        }
                    } else {
                        found = NO;
                        for (FEMVariable *var1 in model.variables) {
                            if ([var1.name isEqualToString:@"mesh update"] == YES) {
                                var1Containers = var1.getContainers;
                                k = i;
                                if (var1Containers->Perm != NULL) k = var1Containers->Perm[k];
                                if (k >= 0) {
                                    for (j=0; j<var1.dofs; j++) {
                                        [self FEMPost_writeDouble:var1Containers->Values[var1.dofs*k+j] toFileHandle:postFileHandle];
                                        [postFileHandle writeData:spaceBuff];
                                    }
                                    if (var1.dofs == 2) {
                                        [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                                    }
                                } else {
                                    [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                                    [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                                    [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                                }
                                found = YES;
                                break;
                            }
                        }
                        if (found == NO) {
                            if (varContainers->CValues != NULL) {
                                [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                                [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                                [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                                [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                                [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                                [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                            } else {
                                [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                                [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                                [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                            }
                        }
                    }
                }
                else if ([variable.name isEqualToString:@"displacement 1"] == YES || [variable.name isEqualToString:@"displacement 2"] ||
                         [variable.name isEqualToString:@"displacement 3"] == YES) {
                    // No ops
                } else if ([variable.name isEqualToString:@"flow solution"] == YES) {
                    k = i;
                    if (varContainers->Perm != NULL) k = varContainers->Perm[k];
                    if (k >= 0) {
                        for (j=0; j<variable.dofs-1; j++) {
                            [self FEMPost_writeDouble:varContainers->Values[variable.dofs*k+j] toFileHandle:postFileHandle];
                            [postFileHandle writeData:spaceBuff];
                        }
                        if (variable.dofs < 4) {
                            [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                        }
                        [self FEMPost_writeDouble:varContainers->Values[variable.dofs*k] toFileHandle:postFileHandle];
                        [postFileHandle writeData:spaceBuff];
                    } else {
                        [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                        [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                        [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                        [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                    }
                } else if ([variable.name isEqualToString:@"velocity 1"] == YES || [variable.name isEqualToString:@"velocity 2"] == YES ||
                           [variable.name isEqualToString:@"velocity 3"] == YES) {
                    // No ops
                } else if ([variable.name isEqualToString:@"magnetic field"] == YES) {
                    k = i;
                    if (varContainers->Perm != NULL) k = varContainers->Perm[k];
                    if (k >= 0) {
                        for (j=0; j<variable.dofs; j++) {
                            [self FEMPost_writeDouble:varContainers->Values[variable.dofs*k+j] toFileHandle:postFileHandle];
                            [postFileHandle writeData:spaceBuff];
                        }
                        if (variable.dofs == 2) {
                            [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                        }
                    } else {
                        [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                        [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                        [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                    }
                } else if ([variable.name isEqualToString:@"magnetic field 1"] == YES || [variable.name isEqualToString:@"magnetic field 2"] == YES ||
                           [variable.name isEqualToString:@"magnetic field 3"] == YES) {
                    // No ops
                } else if ([variable.name isEqualToString:@"coordinate 1"] == YES || [variable.name isEqualToString:@"coordinate 2"] == YES ||
                           [variable.name isEqualToString:@"coordinate 3"] == YES) {
                    // No ops
                } else {
                    if (variable.dofs == 1) {
                        k = i;
                        if (varContainers->Perm != NULL) k = varContainers->Perm[k];
                        if (k >= 0) {
                            if (varContainers->CValues != NULL) {
                                [self FEMPost_writeDouble:creal(*(varContainers->CValues[k])) toFileHandle:postFileHandle];
                                [postFileHandle writeData:spaceBuff];
                                [self FEMPost_writeDouble:cimag(*(varContainers->CValues[k])) toFileHandle:postFileHandle];
                                [postFileHandle writeData:spaceBuff];
                            } else {
                                [self FEMPost_writeDouble:varContainers->Values[k] toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                            }
                        } else {
                            [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                        }
                    } else {
                        ind = [variable.name rangeOfString:@"["];
                        if (ind.location != NSNotFound) {
                            dofs = 0;
                            j = 0;
                            while (1) {
                                str = [variable.name substringFromIndex:j];
                                str = [variable.name substringFromIndex:j];
                                if (j == 0) {
                                    ind = [variable.name rangeOfString:@"["];
                                } else {
                                    ind.location = -1;
                                }
                                ind1 = [str rangeOfString:@":"];
                                if (ind1.location == NSNotFound) errorfunct("FEMPost:writeElmerPostFile", "Missing separator ':' in variable definition using '[ ]' syntax.");
                                if (j == 0) {
                                    if (ind1.location < ind.location) errorfunct("FEMPost:writeElmerPostFile", "Syntax error in variable definition.");
                                } else {
                                    if (ind1.location == 0) errorfunct("FEMPost:writeElmerPostFile", "Syntax error in variable definition.");
                                }
                                l = (int)ind1.location + 1;
                                while (1) {
                                    if ([str characterAtIndex:l] == ' ' || [str characterAtIndex:l] == ']') break;
                                    l++;
                                }
                                q = [[[str substringFromIndex:ind1.location+1] substringToIndex:l-(ind1.location+1)] intValue];
                                k = i;
                                if (varContainers->Perm != NULL) k = varContainers->Perm[k];
                                if (k >= 0) {
                                    if (q == 2 || q == 3) {
                                        if (varContainers->CValues != NULL) {
                                            for (jj=dofs; jj<dofs+q; jj++) {
                                                [self FEMPost_writeDouble:creal(*(varContainers->CValues[variable.dofs*k+jj])) toFileHandle:postFileHandle];
                                                [postFileHandle writeData:spaceBuff];
                                            }
                                            if (q == 2) {
                                                [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                                            }
                                            for (jj=dofs; jj<dofs+q; jj++) {
                                                [self FEMPost_writeDouble:cimag(*(varContainers->CValues[variable.dofs*k+jj])) toFileHandle:postFileHandle];
                                                [postFileHandle writeData:spaceBuff];
                                            }
                                            if (q == 2) {
                                                [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                                            }
                                        } else {
                                            for (jj=dofs; jj<dofs+q; jj++) {
                                                [self FEMPost_writeDouble:varContainers->Values[variable.dofs*k+jj] toFileHandle:postFileHandle];
                                                [postFileHandle writeData:spaceBuff];
                                            }
                                            if (q == 2) {
                                                [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];                                                
                                            }
                                        }
                                    } else {
                                        for (jj=dofs; jj<dofs+q; jj++) {
                                            if (varContainers->CValues != NULL) {
                                                [self FEMPost_writeDouble:creal(*(varContainers->CValues[variable.dofs*k+jj])) toFileHandle:postFileHandle];
                                                [postFileHandle writeData:spaceBuff];
                                                [self FEMPost_writeDouble:cimag(*(varContainers->CValues[variable.dofs*k+jj])) toFileHandle:postFileHandle];
                                                [postFileHandle writeData:spaceBuff];
                                            } else {
                                                [self FEMPost_writeDouble:varContainers->Values[variable.dofs*k+jj] toFileHandle:postFileHandle];
                                                [postFileHandle writeData:spaceBuff];
                                            }
                                        }
                                    }
                                } else {
                                    if (q == 2) {
                                        nZeros = 3;
                                    } else {
                                        nZeros = q;
                                    }
                                    if (varContainers->CValues != NULL) nZeros = 2 * nZeros;
                                    for (jj=0; jj<nZeros; jj++) {
                                        [self FEMPost_writeDouble:0.0 toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                                    }
                                }
                                dofs = dofs + q;
                                while ([str characterAtIndex:i] == ' ') {
                                    l++;
                                }
                                j = l;
                            }
                        }
                    }
                }
            }
            
            if (freeSurfaceFlag == YES) {
                var = [utilities getVariableFrom:model.variables model:model name:@"coordinate 1" onlySearch:NULL maskName:NULL info:&found];
                varContainers = var.getContainers;
                [self FEMPost_writeDouble:varContainers->Values[i] toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                
                var = [utilities getVariableFrom:model.variables model:model name:@"coordinate 2" onlySearch:NULL maskName:NULL info:&found];
                varContainers = var.getContainers;
                [self FEMPost_writeDouble:varContainers->Values[i] toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
                
                var = [utilities getVariableFrom:model.variables model:model name:@"coordinate 3" onlySearch:NULL maskName:NULL info:&found];
                varContainers = var.getContainers;
                [self FEMPost_writeDouble:varContainers->Values[i] toFileHandle:postFileHandle]; [postFileHandle writeData:spaceBuff];
            }
            
            [postFileHandle writeData:newLineBuff];
        }
        if (*append == YES) break;
    }
    
    // We are done here so close the files and deallocate
    [postFileHandle closeFile];
    if (*append == NO) [reader closeHandle];
    if (maskExists == YES) free_ivector(maskOrder, 0, numberOfNodes);
}

@end
