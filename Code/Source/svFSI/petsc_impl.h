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

// This file contains PETSc-dependent data structures. 

#ifndef PETSC_INTERFACE_H
#define PETSC_INTERFACE_H

#include <petscksp.h>
#include <petscao.h>
#include <unistd.h>
#include <stdbool.h>

#include "consts.h"

//--------
// LHSCtx
//---------
// PETSc lhs context. 
//
typedef struct {
    PetscBool created;  /* Whether petsc lhs is created */

    PetscInt  nNo;      /* local number of vertices */
    PetscInt  mynNo;    /* number of owned vertices */

    PetscInt *map;      /* local to local mapping, map[O2] = O1 */
    PetscInt *ltg;      /* local to global in PETSc ordering */
    PetscInt *ghostltg; /* local to global in PETSc ordering */
    PetscInt *rowPtr;   /* row pointer for adjacency info */
    PetscInt *colPtr;   /* column pointer for adjacency info */
} LHSCtx;

//-------
// LSCtx
//-------
// PETSc linear solver context. 
//
typedef struct {
    PetscBool created;  /* Whether mat and vec is created */
    const char *pre;      /* option prefix for different equations */

    PetscInt  lpPts;    /* number of dofs with lumped parameter BC */
    PetscInt *lpBC_l;   /* O2 index for dofs with lumped parameter BC */
    PetscInt *lpBC_g;   /* PETSc index for dofs with lumped parameter BC */

    PetscInt  DirPts;   /* number of dofs with Dirichlet BC */
    PetscInt *DirBC;    /* PETSc index for dofs with Dirichlet BC */

    Vec       b;        /* rhs/solution vector of owned vertices */
    Mat       A;        /* stiffness matrix */
    KSP       ksp;      /* linear solver context */

    PetscBool rcs;      /* whether rcs preconditioner is activated */
    Vec       Dr;       /* diagonal matrix from row maxabs */
    Vec       Dc;       /* diagonal matrix from col maxabs */
} LSCtx;

void petsc_initialize(const PetscInt nNo, const PetscInt mynNo, 
    const PetscInt nnz, const PetscInt nEq, const PetscInt *svFSI_ltg, 
    const PetscInt *svFSI_map, const PetscInt *svFSI_rowPtr, 
    const PetscInt *svFSI_colPtr, char *inp);

void petsc_create_linearsystem(const PetscInt dof, const PetscInt iEq, const PetscInt nEq, 
    const PetscReal *svFSI_DirBC, const PetscReal *svFSI_lpBC);

void petsc_create_linearsolver(const consts::PreconditionerType lsType, const consts::PreconditionerType pcType, 
    const PetscInt kSpace, const PetscInt maxIter, const PetscReal relTol, 
    const PetscReal absTol, const consts::EquationType phys, const PetscInt dof, 
    const PetscInt iEq, const PetscInt nEq);

void petsc_set_values(const PetscInt dof, const PetscInt iEq, const PetscReal *R, 
    const PetscReal *Val, const PetscReal *svFSI_DirBC, const PetscReal *svFSI_lpBC);

void petsc_solve(PetscReal *resNorm,  PetscReal *initNorm,  PetscReal *dB, 
    PetscReal *execTime, bool *converged, PetscInt *numIter, 
    PetscReal *R, const PetscInt maxIter, const PetscInt dof, 
    const PetscInt iEq);

void petsc_destroy_all(const PetscInt);

PetscErrorCode petsc_create_lhs(const PetscInt, const PetscInt, const PetscInt,  
                                const PetscInt *, const PetscInt *, 
                                const PetscInt *, const PetscInt *);

PetscErrorCode petsc_create_bc(const PetscInt, const PetscInt, const PetscReal *, 
                               const PetscReal *);

PetscErrorCode petsc_create_vecmat(const PetscInt, const PetscInt, const PetscInt);

PetscErrorCode petsc_set_vec(const PetscInt, const PetscInt, const PetscReal *);

PetscErrorCode petsc_set_mat(const PetscInt, const PetscInt, const PetscReal *);

PetscErrorCode petsc_set_bc(const PetscInt, const PetscReal *, const PetscReal *);

PetscErrorCode petsc_set_pcfieldsplit(const PetscInt, const PetscInt);

PetscErrorCode petsc_pc_rcs(const PetscInt, const PetscInt);


PetscErrorCode petsc_debug_save_vec(const char *, Vec);
PetscErrorCode petsc_debug_save_mat(const char *, Mat);

#endif
