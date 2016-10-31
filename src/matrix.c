

/******************************************************************************
 * INCLUDES
 *****************************************************************************/
#include "base.h"
#include "matrix.h"
#include "util.h"
#include "timer.h"
#include <math.h>



/******************************************************************************
 * BLAS/LAPACKE FUNCTIONS
 *****************************************************************************/

/* override memory allocators */
#define LAPACK_malloc splatt_malloc
#define LAPACK_free   splatt_free


/* MKL requires mkl_lapacke.h */
#ifndef SPLATT_INC_LAPACKE
#define SPLATT_INC_LAPACKE "lapacke.h"
#endif

#ifndef SPLATT_INC_CBLAS
#define SPLATT_INC_CBLAS   "cblas.h"
#endif

/* Actual includes */
#include SPLATT_INC_LAPACKE
#include SPLATT_INC_CBLAS

#undef I /* <complex.h> has '#define I _Complex_I'. WHY WOULD THEY DO THAT? */
#undef LAPACK_malloc
#undef LAPACK_free


#if   SPLATT_VAL_TYPEWIDTH == 32
  #define SPLATT_CBLAS(x)     cblas_s ## x
  #define SPLATT_LAPACKE(x) LAPACKE_s ## x
#else
  #define SPLATT_CBLAS(x)     cblas_d ## x
  #define SPLATT_LAPACKE(x) LAPACKE_d ## x
#endif





/******************************************************************************
 * PRIVATE FUNCTIONS
 *****************************************************************************/


/**
* @brief Form the Gram matrix from A^T * A.
*
* @param[out] neq_matrix The matrix to fill.
* @param aTa The individual Gram matrices.
* @param mode Which mode we are computing for.
* @param nmodes How many total modes.
* @param reg Regularization parameter (to add to the diagonal).
*/
static void p_form_gram(
    matrix_t * neq_matrix,
    matrix_t * * aTa,
    idx_t const mode,
    idx_t const nmodes,
    val_t const reg)
{
  /* nfactors */
  int N = aTa[0]->J;

  /* form upper-triangual normal equations */
  val_t * const restrict neqs = neq_matrix->vals;
  #pragma omp parallel
  {
    /* first initialize with 1s */
    #pragma omp for schedule(static, 1)
    for(int i=0; i < N; ++i) {
      neqs[i+(i*N)] = 1. + reg;
      for(int j=0; j < N; ++j) {
        neqs[j+(i*N)] = 1.;
      }
    }

    /* now Hadamard product all (A^T * A) matrices */
    for(idx_t m=0; m < nmodes; ++m) {
      if(m == mode) {
        continue;
      }

      val_t const * const restrict mat = aTa[m]->vals;
      #pragma omp for schedule(static, 1)
      for(int i=0; i < N; ++i) {
        /* 
         * `mat` is symmetric but stored upper right triangular, so be careful
         * to only access that.
         */

        /* copy upper triangle */
        for(int j=i; j < N; ++j) {
          neqs[j+(i*N)] *= mat[j+(i*N)];
        }
      }
    } /* foreach mode */

    #pragma omp barrier

    /* copy lower triangular */
    #pragma omp for schedule(static, 1)
    for(int i=0; i < N; ++i) {
      for(int j=0; j < i; ++j) {
        neqs[j+(i*N)] = neqs[i+(j*N)];
      }
    }
  } /* omp parallel */
}



static void p_mat_2norm(
  matrix_t * const A,
  val_t * const restrict lambda,
  rank_info * const rinfo,
  thd_info * const thds)
{
  idx_t const I = A->I;
  idx_t const J = A->J;
  val_t * const restrict vals = A->vals;

  #pragma omp parallel
  {
    int const tid = splatt_omp_get_thread_num();
    val_t * const mylambda = (val_t *) thds[tid].scratch[0];
    for(idx_t j=0; j < J; ++j) {
      mylambda[j] = 0;
    }

    #pragma omp for schedule(static)
    for(idx_t i=0; i < I; ++i) {
      for(idx_t j=0; j < J; ++j) {
        mylambda[j] += vals[j + (i*J)] * vals[j + (i*J)];
      }
    }

    /* do reduction on partial sums */
    thd_reduce(thds, 0, J, REDUCE_SUM);

    #pragma omp master
    {
#ifdef SPLATT_USE_MPI
      /* now do an MPI reduction to get the global lambda */
      timer_start(&timers[TIMER_MPI_NORM]);
      timer_start(&timers[TIMER_MPI_COMM]);
      MPI_Allreduce(mylambda, lambda, J, SPLATT_MPI_VAL, MPI_SUM, rinfo->comm_3d);
      timer_stop(&timers[TIMER_MPI_COMM]);
      timer_stop(&timers[TIMER_MPI_NORM]);
#else
      memcpy(lambda, mylambda, J * sizeof(val_t));
#endif
    }

    #pragma omp barrier

    #pragma omp for schedule(static)
    for(idx_t j=0; j < J; ++j) {
      lambda[j] = sqrt(lambda[j]);
    }

    /* do the normalization */
    #pragma omp for schedule(static)
    for(idx_t i=0; i < I; ++i) {
      for(idx_t j=0; j < J; ++j) {
        vals[j+(i*J)] /= lambda[j];
      }
    }
  } /* end omp for */
}


static void p_mat_maxnorm(
  matrix_t * const A,
  val_t * const restrict lambda,
  rank_info * const rinfo,
  thd_info * const thds)
{
  idx_t const I = A->I;
  idx_t const J = A->J;
  val_t * const restrict vals = A->vals;

  #pragma omp parallel
  {
    int const tid = splatt_omp_get_thread_num();
    val_t * const mylambda = (val_t *) thds[tid].scratch[0];
    for(idx_t j=0; j < J; ++j) {
      mylambda[j] = 0;
    }

    #pragma omp for schedule(static)
    for(idx_t i=0; i < I; ++i) {
      for(idx_t j=0; j < J; ++j) {
        mylambda[j] = SS_MAX(mylambda[j], vals[j+(i*J)]);
      }
    }

    /* do reduction on partial maxes */
    thd_reduce(thds, 0, J, REDUCE_MAX);

    #pragma omp master
    {
#ifdef SPLATT_USE_MPI
      /* now do an MPI reduction to get the global lambda */
      timer_start(&timers[TIMER_MPI_NORM]);
      timer_start(&timers[TIMER_MPI_COMM]);
      MPI_Allreduce(mylambda, lambda, J, SPLATT_MPI_VAL, MPI_MAX, rinfo->comm_3d);
      timer_stop(&timers[TIMER_MPI_COMM]);
      timer_stop(&timers[TIMER_MPI_NORM]);
#else
      memcpy(lambda, mylambda, J * sizeof(val_t));
#endif

    }

    #pragma omp barrier

    #pragma omp for schedule(static)
    for(idx_t j=0; j < J; ++j) {
      lambda[j] = SS_MAX(lambda[j], 1.);
    }

    /* do the normalization */
    #pragma omp for schedule(static)
    for(idx_t i=0; i < I; ++i) {
      for(idx_t j=0; j < J; ++j) {
        vals[j+(i*J)] /= lambda[j];
      }
    }
  } /* end omp parallel */
}


/**
* @brief Solve the system LX = B.
*
* @param L The lower triangular matrix of coefficients.
* @param B The right-hand side which is overwritten with X.
*/
static void p_mat_forwardsolve(
  matrix_t const * const L,
  matrix_t * const B)
{
  /* check dimensions */
  idx_t const N = L->I;

  val_t const * const restrict lv = L->vals;
  val_t * const restrict bv = B->vals;

  /* first row of X is easy */
  for(idx_t j=0; j < N; ++j) {
    bv[j] /= lv[0];
  }

  /* now do forward substitution */
  for(idx_t i=1; i < N; ++i) {
    /* X(i,f) = B(i,f) - \sum_{j=0}^{i-1} L(i,j)X(i,j) */
    for(idx_t j=0; j < i; ++j) {
      for(idx_t f=0; f < N; ++f) {
        bv[f+(i*N)] -= lv[j+(i*N)] * bv[f+(j*N)];
      }
    }
    for(idx_t f=0; f < N; ++f) {
      bv[f+(i*N)] /= lv[i+(i*N)];
    }
  }
}


/**
* @brief Solve the system UX = B.
*
* @param U The upper triangular matrix of coefficients.
* @param B The right-hand side which is overwritten with X.
*/
static void p_mat_backwardsolve(
  matrix_t const * const U,
  matrix_t * const B)
{
  /* check dimensions */
  idx_t const N = U->I;

  val_t const * const restrict rv = U->vals;
  val_t * const restrict bv = B->vals;

  /* last row of X is easy */
  for(idx_t f=0; f < N; ++f) {
    idx_t const i = N-1;
    bv[f+(i*N)] /= rv[i+(i*N)];
  }

  /* now do backward substitution */
  for(idx_t row=2; row <= N; ++row) {
    /* operate with (N - row) to make unsigned comparisons easy */
    idx_t const i = N - row;

    /* X(i,f) = B(i,f) - \sum_{j=0}^{i-1} R(i,j)X(i,j) */
    for(idx_t j=i+1; j < N; ++j) {
      for(idx_t f=0; f < N; ++f) {
        bv[f+(i*N)] -= rv[j+(i*N)] * bv[f+(j*N)];
      }
    }
    for(idx_t f=0; f < N; ++f) {
      bv[f+(i*N)] /= rv[i+(i*N)];
    }
  }
}

/******************************************************************************
 * PUBLIC FUNCTIONS
 *****************************************************************************/


void mat_cholesky(
  matrix_t const * const A,
  matrix_t * const L)
{
  /* check dimensions */
  assert(A->I == A->J);
  assert(A->I == L->J);
  assert(L->I == L->J);

  idx_t const N = A->I;
  val_t const * const restrict av = A->vals;
  val_t * const restrict lv = L->vals;

  memset(lv, 0, N*N*sizeof(val_t));
  for (idx_t i = 0; i < N; ++i) {
    for (idx_t j = 0; j <= i; ++j) {
      val_t inner = 0;
      for (idx_t k = 0; k < j; ++k) {
        inner += lv[k+(i*N)] * lv[k+(j*N)];
      }

      if(i == j) {
        lv[j+(i*N)] = sqrt(av[i+(i*N)] - inner);
      } else {
        lv[j+(i*N)] = 1.0 / lv[j+(j*N)] * (av[j+(i*N)] - inner);
      }
    }
  }
}


void mat_aTa(
  matrix_t const * const A,
  matrix_t * const ret,
  rank_info * const rinfo,
  thd_info * const thds,
  idx_t const nthreads)
{
  timer_start(&timers[TIMER_ATA]);
  /* check matrix dimensions */
  assert(ret->I == ret->J);
  assert(ret->I == A->J);
  assert(ret->vals != NULL);
  assert(A->rowmajor);
  assert(ret->rowmajor);

  /* A^T * A */
  SPLATT_CBLAS(syrk)(
      CblasRowMajor, CblasUpper, CblasTrans,
      A->J, A->I, /* swapped due to trans */
      1.,
      A->vals, A->J,
      0.,
      ret->vals, A->J);

#ifdef SPLATT_USE_MPI
  timer_start(&timers[TIMER_MPI_ATA]);
  timer_start(&timers[TIMER_MPI_COMM]);
  MPI_Allreduce(MPI_IN_PLACE, ret->vals, F * F, SPLATT_MPI_VAL, MPI_SUM,
      rinfo->comm_3d);
  timer_stop(&timers[TIMER_MPI_COMM]);
  timer_stop(&timers[TIMER_MPI_ATA]);
#endif

  timer_stop(&timers[TIMER_ATA]);
}

void mat_matmul(
  matrix_t const * const A,
  matrix_t const * const B,
  matrix_t  * const C)
{
  timer_start(&timers[TIMER_MATMUL]);
  /* check dimensions */
  assert(A->J == B->I);
  assert(C->I * C->J <= A->I * B->J);

  /* set dimensions */
  C->I = A->I;
  C->J = B->J;

  val_t const * const restrict av = A->vals;
  val_t const * const restrict bv = B->vals;
  val_t       * const restrict cv = C->vals;

  idx_t const M  = A->I;
  idx_t const N  = B->J;
  idx_t const Na = A->J;

  /* tiled matrix multiplication */
  idx_t const TILE = 16;
  #pragma omp parallel for schedule(static)
  for(idx_t i=0; i < M; ++i) {
    for(idx_t jt=0; jt < N; jt += TILE) {
      for(idx_t kt=0; kt < Na; kt += TILE) {
        idx_t const JSTOP = SS_MIN(jt+TILE, N);
        for(idx_t j=jt; j < JSTOP; ++j) {
          val_t accum = 0;
          idx_t const KSTOP = SS_MIN(kt+TILE, Na);
          for(idx_t k=kt; k < KSTOP; ++k) {
            accum += av[k + (i*Na)] * bv[j + (k*N)];
          }
          cv[j + (i*N)] += accum;
        }
      }
    }
  }

  timer_stop(&timers[TIMER_MATMUL]);
}

void mat_normalize(
  matrix_t * const A,
  val_t * const restrict lambda,
  splatt_mat_norm const which,
  rank_info * const rinfo,
  thd_info * const thds,
  idx_t const nthreads)
{
  timer_start(&timers[TIMER_MATNORM]);

  splatt_omp_set_num_threads(nthreads);

  switch(which) {
  case MAT_NORM_2:
    p_mat_2norm(A, lambda, rinfo, thds);
    break;
  case MAT_NORM_MAX:
    p_mat_maxnorm(A, lambda, rinfo, thds);
    break;
  default:
    fprintf(stderr, "SPLATT: mat_normalize supports 2 and MAX only.\n");
    abort();
  }
  timer_stop(&timers[TIMER_MATNORM]);
}



void mat_solve_normals(
  idx_t const mode,
  idx_t const nmodes,
	matrix_t * * aTa,
  matrix_t * rhs,
  val_t const reg)
{
  timer_start(&timers[TIMER_INV]);

  p_form_gram(aTa[MAX_NMODES], aTa, mode, nmodes, reg);

  lapack_int info;
  char uplo = 'L';

  /* nfactors */
  lapack_int N = aTa[0]->J;
  lapack_int lda = N;
  lapack_int ldb = N;
  lapack_int order = N;
  lapack_int nrhs = rhs->I;

  val_t * const neqs = aTa[MAX_NMODES]->vals;

  /* Cholesky factorization */
  bool is_spd = true;
  info = SPLATT_LAPACKE(potrf)(
      LAPACK_COL_MAJOR, uplo, order,
      neqs, lda);
  if(info) {
    fprintf(stderr, "SPLATT: Gram matrix is not SPD. Trying `GELSS`.\n");
    is_spd = false;
  }

  /* Continue with Cholesky */
  if(is_spd) {
    /* Solve against rhs */
    info = SPLATT_LAPACKE(potrs)(
        LAPACK_COL_MAJOR, uplo,
        order, nrhs,
        neqs, lda,
        rhs->vals, ldb);
    if(info) {
      fprintf(stderr, "SPLATT: DPOTRS returned %d\n", info);
    }
  } else {
    /* restore gram matrix */
    p_form_gram(aTa[MAX_NMODES], aTa, mode, nmodes, reg);

    /* Use column major to avoid any internal transposition. */
    val_t * conditions = splatt_malloc(N * sizeof(*conditions));
    val_t rcond = -1.0f;
    lapack_int effective_rank;
    info = SPLATT_LAPACKE(gelss)(
        LAPACK_COL_MAJOR, N, N, nrhs,
        neqs, lda,
        rhs->vals, ldb,
        conditions, rcond, &effective_rank);
    if(info) {
      printf("SPLATT: DGELSS returned %d\n", info);
    }
    printf("SPLATT:   DGELSS effective rank: %d\n", effective_rank);
    splatt_free(conditions);
  }

  timer_stop(&timers[TIMER_INV]);
}




matrix_t * mat_alloc(
  idx_t const nrows,
  idx_t const ncols)
{
  matrix_t * mat = (matrix_t *) splatt_malloc(sizeof(matrix_t));
  mat->I = nrows;
  mat->J = ncols;
  mat->vals = (val_t *) splatt_malloc(nrows * ncols * sizeof(val_t));
  mat->rowmajor = 1;
  return mat;
}

matrix_t * mat_rand(
  idx_t const nrows,
  idx_t const ncols)
{
  matrix_t * mat = mat_alloc(nrows, ncols);
  val_t * const vals = mat->vals;

  fill_rand(vals, nrows * ncols);

  return mat;
}

void mat_free(
  matrix_t * mat)
{
  free(mat->vals);
  free(mat);
}

matrix_t * mat_mkrow(
  matrix_t const * const mat)
{
  assert(mat->rowmajor == 0);

  idx_t const I = mat->I;
  idx_t const J = mat->J;

  matrix_t * row = mat_alloc(I, J);
  val_t       * const restrict rowv = row->vals;
  val_t const * const restrict colv = mat->vals;

  for(idx_t i=0; i < I; ++i) {
    for(idx_t j=0; j < J; ++j) {
      rowv[j + (i*J)] = colv[i + (j*I)];
    }
  }

  return row;
}

matrix_t * mat_mkcol(
  matrix_t const * const mat)
{
  assert(mat->rowmajor == 1);
  idx_t const I = mat->I;
  idx_t const J = mat->J;

  matrix_t * col = mat_alloc(I, J);
  val_t       * const restrict colv = col->vals;
  val_t const * const restrict rowv = mat->vals;

  for(idx_t i=0; i < I; ++i) {
    for(idx_t j=0; j < J; ++j) {
      colv[i + (j*I)] = rowv[j + (i*J)];
    }
  }

  col->rowmajor = 0;

  return col;
}


spmatrix_t * spmat_alloc(
  idx_t const nrows,
  idx_t const ncols,
  idx_t const nnz)
{
  spmatrix_t * mat = (spmatrix_t*) splatt_malloc(sizeof(spmatrix_t));
  mat->I = nrows;
  mat->J = ncols;
  mat->nnz = nnz;
  mat->rowptr = (idx_t*) splatt_malloc((nrows+1) * sizeof(idx_t));
  mat->colind = (idx_t*) splatt_malloc(nnz * sizeof(idx_t));
  mat->vals   = (val_t*) splatt_malloc(nnz * sizeof(val_t));
  return mat;
}

void spmat_free(
  spmatrix_t * mat)
{
  free(mat->rowptr);
  free(mat->colind);
  free(mat->vals);
  free(mat);
}

