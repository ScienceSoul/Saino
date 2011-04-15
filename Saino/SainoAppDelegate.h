//
//  SainoAppDelegate.h
//  Saino
//
//  Created by Seddik hakime on 15/04/11.
//  Copyright 2011 Institute of Low Temperature Science. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@interface SainoAppDelegate : NSObject <NSApplicationDelegate> {
@private
    NSWindow *window;
}

@property (assign) IBOutlet NSWindow *window;

@end
