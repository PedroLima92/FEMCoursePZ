//
//  Header.h
//  FEMCoursePZ
//
//  Created by Philippe Devloo on 7/30/15.
//
//

#ifndef FEMCoursePZ_Panic_h
#define FEMCoursePZ_Panic_h

#include <iostream>
#include <exception>

static void DebugStop()
{
    std::cout << "Your chance to put a breakpoint here\n";
    std::bad_exception myex;
    throw myex;

}
#endif
