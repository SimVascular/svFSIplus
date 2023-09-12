/* 
 * Author: Chi Zhu, Peking University
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

#ifndef PETSC_LINEAR_SOLVER_H
#define PETSC_LINEAR_SOLVER_H

#include <petscksp.h>
#include <petscao.h>
#include <unistd.h>
#include <stdbool.h>

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    User-defined data type
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* PETSc lhs context */
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

/* PETSc linear solver context */
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


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    Global variables
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
LHSCtx         plhs;       /* PETSc lhs */
LSCtx         *psol;       /* PETSc solver */
PetscLogStage  stages[6];  /* performance tuning. */


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
   Macro definitions
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* For linear solver and preconditioner */
#define PETSc_CG     798
#define PETSc_GMRES  797
#define PETSc_BICGS  795
    
#define PETSc_PC_FSILS  701
#define PETSc_PC_RCS    709
#define PETSc_PC        713

/* Physics solved */
#define EQ_fluid   201
#define EQ_struct  202
#define EQ_heatS   203
#define EQ_lElas   204
#define EQ_heatF   205 
#define EQ_FSI     206 
#define EQ_mesh    207
#define EQ_shell   208 
#define EQ_CMM     209 
#define EQ_CEP     210
#define EQ_ustruct 211 
#define EQ_stokes  212


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    Functions that interact with Fortran
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#ifdef __cplusplus
extern "C"{
#endif

void petsc_initialize_(const PetscInt *, const PetscInt *, const PetscInt *, \
                       const PetscInt *, const PetscInt *, const PetscInt *, \
                       const PetscInt *, const PetscInt *, char *);
void petsc_create_linearsystem_(const PetscInt *, const PetscInt *, const PetscInt *, \
                                const PetscReal *, const PetscReal *);
void petsc_create_linearsolver_(const PetscInt *, const PetscInt *, const PetscInt *, \
                                const PetscInt *, const PetscReal *, const PetscReal *, \
                                const PetscInt *, const PetscInt *, const PetscInt *, \
                                const PetscInt *);
void petsc_set_values_(const PetscInt *, const PetscInt *, const PetscReal *, \
                       const PetscReal *, const PetscReal *, const PetscReal *);
void petsc_solve_(PetscReal *, PetscReal *, PetscReal *, PetscReal *, bool *, \
                  PetscInt *,  PetscReal *, const PetscInt *, const PetscInt *, \
                  const PetscInt *);
void petsc_destroy_all_(const PetscInt *);

#ifdef __cplusplus
}
#endif


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    Private functions
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#ifdef __cplusplus
extern "C"{
#endif

PetscErrorCode petsc_create_lhs(const PetscInt, const PetscInt, const PetscInt,  \
                                const PetscInt *, const PetscInt *, \
                                const PetscInt *, const PetscInt *);
PetscErrorCode petsc_create_bc(const PetscInt, const PetscInt, const PetscReal *, \
                               const PetscReal *);
PetscErrorCode petsc_create_vecmat(const PetscInt, const PetscInt, const PetscInt);
PetscErrorCode petsc_set_vec(const PetscInt, const PetscInt, const PetscReal *);
PetscErrorCode petsc_set_mat(const PetscInt, const PetscInt, const PetscReal *);
PetscErrorCode petsc_set_bc(const PetscInt, const PetscReal *, const PetscReal *);
PetscErrorCode petsc_set_pcfieldsplit(const PetscInt, const PetscInt);
PetscErrorCode petsc_pc_rcs(const PetscInt, const PetscInt);

#ifdef __cplusplus
}
#endif


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
    Private functions for debugging
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#ifdef __cplusplus
extern "C"{
#endif

PetscErrorCode petsc_debug_save_vec(const char *, Vec);
PetscErrorCode petsc_debug_save_mat(const char *, Mat);

#ifdef __cplusplus
}
#endif

char * rm_blank(char *string);

#endif
