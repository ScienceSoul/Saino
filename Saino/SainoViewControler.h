//
//  SainoViewControler.h
//  Saino
//
//  Created by Seddik hakime on 24/06/2015.
//  Copyright (c) 2015 Institute of Low Temperature Science. All rights reserved.
//

#import <Cocoa/Cocoa.h>

#import <SainoCore/FEMJob.h>

@interface SainoViewControler : NSViewController {
    FEMJob * __nullable _job;
}

@property (weak, nullable) IBOutlet NSTextField *displayField;

-(IBAction)press:(id __nullable)sender;

@end
