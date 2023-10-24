/**
 * Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

// The classes defined here duplicate the data structures in the Fortran CHNLMOD module
// defined in CHNL.f. 

// This module defines data structures for 

#include "Array.h"
#include "Vector.h"

#include <string>

#ifndef CHNL_MOD_H 
#define CHNL_MOD_H 

/// @brief Channel type, used in I/O
class chnlType
{
  // Whether it is open to the screen
  bool oTS = false;

  // Whether it is open to the file
  bool oTF = false; 

  // Channel identifier
  int id;

  // File ID
  int fId;

  // File address
  std::string fName = "histor";

  // Channel tag
  std::string tag = "";

  //         Creates a new channel
  //       PROCEDURE :: new => newChnl
  //         Closes the channel
  //  PROCEDURE :: close => closeChnl
  //         To send a string to channel
   //        PROCEDURE chnlAssign
      //     GENERIC :: ASSIGNMENT(=) => chnlAssign
};

/// @brief Only to group four channels, in case I rather have them as one
/// variable
class ioType
{
  // Standard output
  chnlType o;

  // Error
  chnlType e;

  // Warning
  chnlType w;

  // Debugging
  chnlType d;

  // Status file
  chnlType s;

  //      CONTAINS
  //!        Opens all as standard channels
           //PROCEDURE :: new => newIO
  //!        Closes the channel
           //PROCEDURE :: close => closeIO
};

class ChnlMod
{
  public:
    ChnlMod();

    // Channels cases: standard output, error output, warning output,
    // debugging output
    int CHNL_O = 601;
    int CHNL_E = 602;
    int CHNL_W = 603;
    int CHNL_D = 604;

    // Whether to use color in printing outputs
    bool pClr = true;

    // A general counter for file ID
    int gFID = 314;

    // Appended path to all files that are going to be saved
    std::string appPath = "";

};

#endif

