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
    FEMJob *_job;
}

@property (weak) IBOutlet NSTextField *displayField;

-(IBAction)press:(id)sender;

@end
