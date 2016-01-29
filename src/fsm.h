/*
Copyright (c) 2014-2016, Alexandr Nikitin
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of seis-eikonal nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef H_FSM
#define H_FSM

#include <omp.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>

/**  memory address for (i,j,k) index in 3d array with sizes (NI,NJ,NK).*/
#define M_MEM_ADR_3D(i,j,k,NI,NJ,NK) M_MEM_ADR_3D_RW(i,j,k,NI,NJ,NK)

/**  memory address for (i,j,k) index in row major 3d array with sizes
 * (NI,NJ,NK).*/
#define M_MEM_ADR_3D_RW(i,j,k,NI,NJ,NK) \
    ((i) * (NJ) * (NK) + (j) * (NK) + (k))

/**  memory address for (i,j,k) index in column major 3d array with sizes
 * (NI,NJ,NK).*/
#define M_MEM_ADR_3D_CW(i,j,k,NI,NJ,NK) \
    ((k) * (NJ) * (NI) + (j) * (NI) + (i))

/**
 * @brief FSM3D Fast Sweeping Method with first order upwind scheme in 3D space.
 * @param[in,out] U first arrival travel times
 * @param[in] F slowness
 * @param[in] H grid spacing
 * @param[in] NI I dimension of numerical grid in
 * @param[in] NJ J dimension of numerical grid in
 * @param[in] NK K dimension of numerical grid in
 * @param[in] SRCI I index of source (fixed U value)
 * @param[in] SRCJ J index of source (fixed U value)
 * @param[in] SRCK K index of source (fixed U value)
 * @param[in] max_iter
 * @return number of iterations completed
 * @author Alexandr Nikitin
 */
int FSM3D(double *U, double *F, double H, size_t NI, size_t NJ, size_t NK,
    size_t SRCI, size_t SRCJ, size_t SRCK, int max_iter, int *converged);

/**
 * @brief FSM3DInit Initialize FSM3D before starting computation.
 * @param[in,out] U first arrival travel times
 * @param[in] NI I dimension of numerical grid in
 * @param[in] NJ J dimension of numerical grid in
 * @param[in] NK K dimension of numerical grid in
 * @param[in] SRCI I index of source (fixed U value)
 * @param[in] SRCJ J index of source (fixed U value)
 * @param[in] SRCK K index of source (fixed U value)
 * @author Alexandr Nikitin
 */
void FSM3DInit(double *U, size_t NI, size_t NJ, size_t NK,
               size_t SRCI, size_t SRCJ, size_t SRCK);

int FSM3D_node_update(double *U, double *F, double H,
                      size_t NI, size_t NJ, size_t NK,
                      size_t SRCI, size_t SRCJ, size_t SRCK,
                      int REVI, int REVJ, int REVK,
                      size_t ir, size_t jr, size_t kr);

int FSM3D_serial(double *U, double *F, double H,
                 size_t NI, size_t NJ, size_t NK,
                 size_t SRCI, size_t SRCJ, size_t SRCK,
                 int REVI, int REVJ, int REVK,
                 size_t istart, size_t jstart, size_t kstart,
                 size_t isize, size_t jsize, size_t ksize);

#endif
