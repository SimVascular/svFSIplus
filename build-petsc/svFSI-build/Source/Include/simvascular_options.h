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

/* #undef SV_USE_TIMER */
#define SV_USE_ZLIB ON 
/* #undef SV_USE_ITK */

/* #undef SV_USE_PARASOLID */
/* #undef SV_USE_PARASOLID_SHARED */
/* #undef SV_USE_OPENCASCADE */
/* #undef SV_USE_OPENCASCADE_SHARED */
/* #undef SV_USE_MESHSIM_DISCRETE_MODEL */
/* #undef SV_USE_MESHSIM_DISCRETE_MODEL_SHARED */

/* #undef SV_USE_MESHSIM */
/* #undef MESHSIM_LICENSE_IN_WIN32_REGISTRY */
/* #undef SV_USE_MESHSIM_SHARED */
/* #undef SV_USE_MESHSIM_ADAPTOR */
/* #undef SV_USE_PYTHON */
#define SV_USE_TETGEN ON
/* #undef SV_USE_MMG */
#ifdef SV_USE_TETGEN
/* #undef TETGEN150 */
/* #undef TETGEN143 */
/* #undef TETGEN151 */
#endif
/* #undef SV_USE_TET_ADAPTOR */
/* #undef SV_USE_VMTK */

/* #undef SV_USE_WIN32_REGISTRY */
#define SV_REGISTRY_TOPLEVEL "SV"
#define SV_USE_NOTIMER ON

/* #undef SV_USE_QT_GUI */
/* #undef SV_USE_QT */
/* #undef SV_USE_MITK */
/* #undef SV_NO_PYTHONQT_ALL */
