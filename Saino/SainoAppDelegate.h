//
//  SainoAppDelegate.h
//  Saino
//
//  Created by Seddik hakime on 15/04/11.
//  Copyright 2011 Institute of Low Temperature Science. All rights reserved.
//

#import <Cocoa/Cocoa.h>

#import "SainoViewControler.h"

@interface SainoAppDelegate : NSObject <NSApplicationDelegate> {
@private
    NSWindow *window;
    SainoViewControler *saino_view_controller;
}

@property (strong) IBOutlet NSWindow *window;
@property (nonatomic) SainoViewControler *saino_view_controller;

@end
