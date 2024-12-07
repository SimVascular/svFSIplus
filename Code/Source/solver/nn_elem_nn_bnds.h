/* Copyright (c) Stanford University, The Regents of the University of California, and others.
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

/// @brief Define a map type used to set the bounds of element shape functions.
///
/// Sets   
/// \code {.cpp}
///   mesh.xib(2, com_mod.nsd);   // Bounds on Gauss integration points in parametric space
///   mesh.Nb2, mesh.eNoN);       // Bounds on shape functions
/// \endcode

using SetElementShapeBoundsMapType = std::map<ElementType, std::function<void(int, mshType&)>>;

SetElementShapeBoundsMapType set_element_shape_bounds_data = {

  {ElementType::HEX8, [](int nsd, mshType& mesh) -> void { 
      for (int i = 0; i < nsd; i++) {
        mesh.xib(0,i) = -1.0;
        mesh.xib(1,i) =  1.0;
      }
    }
  },

  {ElementType::LIN1, [](int nsd, mshType& mesh) -> void { 
    std::cout << "[set_element_shape_bounds_data] **************************" << std::endl;
    std::cout << "[set_element_shape_bounds_data] ERROR: LIN1 not supported." << std::endl;
    std::cout << "[set_element_shape_bounds_data] **************************" << std::endl;
    }
  },

  {ElementType::TET4, [](int nsd, mshType& mesh) -> void { 
      for (int i = 0; i < nsd; i++) {
        mesh.xib(0,i) = 0.0;
      }
    }
  },

  {ElementType::TRI3, [](int nsd, mshType& mesh) -> void { 
      for (int i = 0; i < nsd; i++) {
        mesh.xib(0,i) = 0.0;
      }
    }
  },

  {ElementType::WDG, [](int nsd, mshType& mesh) -> void { 
      mesh.xib(0,0) = 0.0;
      mesh.xib(0,1) = 0.0;
    }
  },

};

