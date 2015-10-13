//
//  SainoCore.h
//  SainoCore
//
//  Created by Seddik hakime on 09/10/2015.
//  Copyright © 2015 Institute of Low Temperature Science. All rights reserved.
//

#import <Cocoa/Cocoa.h>

//! Project version number for SainoCore.
FOUNDATION_EXPORT double SainoCoreVersionNumber;

//! Project version string for SainoCore.
FOUNDATION_EXPORT const unsigned char SainoCoreVersionString[];

// In this header, you should import all the public headers of your framework using statements like #import <SainoCore/PublicHeader.h>

#import <SainoCore/FEMBandwidthOptimize.h>
#import <SainoCore/FEMBodyForce.h>
#import <SainoCore/FEMBoundary.h>
#import <SainoCore/FEMBoundaryCondition.h>
#import <SainoCore/FEMConstants.h>
#import <SainoCore/FEMCoordinateSystems.h>
#import <SainoCore/FEMCore.h>
#import <SainoCore/FEMDifferentials.h>
#import <SainoCore/FEMDiffuseConvectiveAnisotropic.h>
#import <SainoCore/FEMDiffuseConvectiveGeneralAnisotropic.h>
#import <SainoCore/FEMElementDescription.h>
#import <SainoCore/FEMElementsDefinition.h>
#import <SainoCore/FEMElementUtils.h>
#import <SainoCore/FEMEquation.h>
#import <SainoCore/FEMFreeSurface.h>
#import <SainoCore/FEMHUTIter.h>
#import <SainoCore/FEMInitialConditions.h>
#import <SainoCore/FEMInterpolation.h>
#import <SainoCore/FEMIterativeMethods.h>
#import <SainoCore/FEMJob.h>
#import <SainoCore/FEMLinearAlgebra.h>
#import <SainoCore/FEMListMatrix.h>
#import <SainoCore/FEMListUtilities.h>
#import <SainoCore/FEMMaterial.h>
#import <SainoCore/FEMMaterialModels.h>
#import <SainoCore/FEMMatrix.h>
#import <SainoCore/FEMMatrixBand.h>
#import <SainoCore/FEMMatrixCRS.h>
#import <SainoCore/FEMMesh.h>
#import <SainoCore/FEMMeshUtils.h>
#import <SainoCore/FEMModel.h>
#import <SainoCore/FEMNavierStokes.h>
#import <SainoCore/FEMNavierStokesCylindrical.h>
#import <SainoCore/FEMNavierStokesGeneral.h>
#import <SainoCore/FEMNumericIntegration.h>
#import <SainoCore/FEMParallelMPI.h>
#import <SainoCore/FEMPElementMaps.h>
#import <SainoCore/FEMPost.h>
#import <SainoCore/FEMPrecondition.h>
#import <SainoCore/FEMProjector.h>
#import <SainoCore/FEMRadiation.h>
#import <SainoCore/FEMSimulation.h>
#import <SainoCore/FEMSolution.h>
#import <SainoCore/FEMTimeIntegration.h>
#import <SainoCore/FEMUtilities.h>
#import <SainoCore/FEMValueList.h>
#import <SainoCore/FEMVariable.h>
#import <SainoCore/SIOInfoParallel.h>
#import <SainoCore/SIOMeshAgent.h>
#import <SainoCore/SIOMeshIO.h>
#import <SainoCore/SIOModelManager.h>
#import <SainoCore/FEMFlowSolution.h>
#import <SainoCore/FEMHeatSolution.h>
#import <SainoCore/FEMHeatSolution_OpenCL.h>
#import <SainoCore/FEMMagneticInductionSolution.h>
#import <SainoCore/FEMMeshUpdateSolution.h>
#import <SainoCore/FEMStressAnalysisSolution.h>
#import <SainoCore/FemTest.h>

#import <SainoCore/Constructors.h>
#import <SainoCore/DirectoryReader.h>
#import <SainoCore/FileReader.h>
#import <SainoCore/GaussIntegration.h>
#import <SainoCore/huti_defs.h>
#import <SainoCore/memory.h>
#import <SainoCore/NodalBasisFunctions.h>
#import <SainoCore/NodeCompare.h>
#import <SainoCore/NSDataExtensions.h>
#import <SainoCore/Numerics.h>
#import <SainoCore/OpenCLUtils.h>
#import <SainoCore/SainoSolutionsComputer.h>
#import <SainoCore/sio_config.h>
#import <SainoCore/TimeProfile.h>
#import <SainoCore/Utils.h>
#import <SainoCore/Walls.h>
