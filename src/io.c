/*
 * io.c
 *
 * Files of models of the fast method: the set-up of a large model (a high
 * degree, or many fields on one grid) takes seconds to minutes and is the
 * same on every run, so a model can be written to a file after its set-up and
 * read instead of being set up again.
 *
 * The file holds what the evaluation needs: the grid, the radial basis, the
 * window and the output map. It does not hold the coefficients, and a model
 * that was read cannot be changed. The file is in the binary format of the
 * machine (not portable between byte orders).
 *
 * Layout: 8 bytes "PSYNTHW1"; a 64-bit byte order mark; NHEAD_U unsigned and
 * NHEAD_I signed 64-bit integers; NHEAD_D doubles; the Chebyshev coefficients
 * of the radial basis; the output map (if any); the grid.
 *
 * Copyright (C) 2026 Patrick Alken and Johnathan Burchill, BSD 3-Clause License
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "potsynth_internal.h"

#define IO_MAGIC                     "PSYNTHW1"
#define IO_BOM                       UINT64_C(0x0102030405060708)
#define IO_MAXIDX                    (1L << 40) /* bound on grid row and column indices */
#define NHEAD_U                      17
#define NHEAD_I                      4
#define NHEAD_D                      15

/*
potsynth_write()
  Write a model of the fast method to a file

Inputs: filename - file name
        w        - model (fast method)

Return: success or error

Notes:
1) Radial functions of the caller are not stored: give them again to
potsynth_read()
*/

int
potsynth_write(const char * filename, const potsynth_workspace * w)
{
  if (filename == NULL || w == NULL)
    {
      POTSYNTH_ERROR ("file name and model must be given", POTSYNTH_EINVAL);
    }
  else if (w->method != POTSYNTH_METHOD_FAST)
    {
      POTSYNTH_ERROR ("only models of the fast method can be written", POTSYNTH_EUNSUP);
    }
  else
    {
      const potsynth_radial_workspace * rad = w->radial_workspace_p;
      const uint64_t bom = IO_BOM;
      const uint64_t hu[NHEAD_U] = { w->lmin, w->lmax, w->mmax, w->nfunc, w->nout, w->nrad, w->ncomp,
                                     (uint64_t) w->w, w->nf1, w->nf2, w->nrows, w->ncols, w->nch, w->nchp,
                                     rad->rho, rad->ncheb, (w->omap != NULL) };
      const int64_t hi[NHEAD_I] = { w->row_lo, w->row_hi, w->col_lo, w->col_hi };
      const double hd[NHEAD_D] = { w->a, w->rmin, w->rmax, w->th_lo, w->th_hi, w->sigma, w->beta, w->half,
                                   w->h1, w->h2, w->inv_h1, w->inv_h2, rad->uc, rad->ud, w->params.eps };
      const size_t ngrid = w->nrows * w->ncols * w->nchp;
      FILE * fp = fopen(filename, "wb");
      size_t ok;

      if (fp == NULL)
        {
          POTSYNTH_ERROR ("unable to open the model file for writing", POTSYNTH_EINVAL);
        }

      ok = (fwrite(IO_MAGIC, 1, 8, fp) == 8);
      ok = ok && (fwrite(&bom, sizeof(bom), 1, fp) == 1);
      ok = ok && (fwrite(hu, sizeof(uint64_t), NHEAD_U, fp) == NHEAD_U);
      ok = ok && (fwrite(hi, sizeof(int64_t), NHEAD_I, fp) == NHEAD_I);
      ok = ok && (fwrite(hd, sizeof(double), NHEAD_D, fp) == NHEAD_D);
      ok = ok && (fwrite(rad->basis, sizeof(double), rad->rho * rad->ncheb, fp) == rad->rho * rad->ncheb);

      if (w->omap != NULL)
        ok = ok && (fwrite(w->omap, sizeof(double), w->nout * w->nfunc, fp) == w->nout * w->nfunc);

      ok = ok && (fwrite(w->grid, sizeof(double), ngrid, fp) == ngrid);
      ok = (fclose(fp) == 0) && ok;

      if (!ok)
        {
          POTSYNTH_ERROR ("error writing the model file", POTSYNTH_EINVAL);
        }

      return POTSYNTH_SUCCESS;
    }
}

/* *p = a * b; 0 if the product overflows */
static int
io_mul(const size_t a, const size_t b, size_t * p)
{
  if (b != 0 && a > SIZE_MAX / b)
    return 0;

  *p = a * b;
  return 1;
}

/*
potsynth_read()
  Read a model written by potsynth_write()

Inputs: filename - file name
        radial   - radial functions of the caller, for a model that was set up
                   with them (the same number of functions); NULL otherwise
        nthreads - OpenMP threads of the evaluation; 0 = all available

Return: pointer to workspace, or NULL on error; evaluate it with potsynth_eval()
or potsynth_eval_rows() and free it with potsynth_free()
*/

potsynth_workspace *
potsynth_read(const char * filename, const potsynth_radial * radial, const int nthreads)
{
  FILE * fp;
  potsynth_workspace * w;
  potsynth_radial_workspace * rad;
  char magic[8];
  uint64_t bom = 0, hu[NHEAD_U];
  int64_t hi[NHEAD_I];
  double hd[NHEAD_D];
  size_t ngrid = 0, nbasis = 0, nbytes = 0;
  int ok;

  if (filename == NULL || nthreads < 0)
    {
      POTSYNTH_ERROR_NULL ("file name must be given and nthreads >= 0", POTSYNTH_EINVAL);
    }

  fp = fopen(filename, "rb");
  if (fp == NULL)
    {
      POTSYNTH_ERROR_NULL ("unable to open the model file", POTSYNTH_EINVAL);
    }

  ok = (fread(magic, 1, 8, fp) == 8) && (memcmp(magic, IO_MAGIC, 8) == 0);
  ok = ok && (fread(&bom, sizeof(bom), 1, fp) == 1) && (bom == IO_BOM);
  ok = ok && (fread(hu, sizeof(uint64_t), NHEAD_U, fp) == NHEAD_U);
  ok = ok && (fread(hi, sizeof(int64_t), NHEAD_I, fp) == NHEAD_I);
  ok = ok && (fread(hd, sizeof(double), NHEAD_D, fp) == NHEAD_D);

  if (!ok)
    {
      fclose(fp);
      POTSYNTH_ERROR_NULL ("not a model file of this version and byte order", POTSYNTH_EINVAL);
    }
  else if ((hu[5] > 0) != (radial != NULL) || (radial != NULL && radial->n != hu[5]) ||
           (radial != NULL && radial->f == NULL && radial->n != 1))
    {
      fclose(fp);
      POTSYNTH_ERROR_NULL ("the radial functions do not match those of the model file", POTSYNTH_EINVAL);
    }

  w = calloc(1, sizeof(potsynth_workspace));
  rad = calloc(1, sizeof(potsynth_radial_workspace));
  if (w == NULL || rad == NULL)
    {
      fclose(fp);
      free(w);
      free(rad);
      POTSYNTH_ERROR_NULL ("unable to allocate workspace", POTSYNTH_ENOMEM);
    }

  w->radial_workspace_p = rad;
  w->check_err = -1.0; /* not stored in the file */
  rad->err = -1.0;
  w->params = potsynth_default_parameters();
  w->method = POTSYNTH_METHOD_FAST;
  w->nthreads = potsynth_max_threads(nthreads);

  w->lmin = (size_t) hu[0];
  w->lmax = (size_t) hu[1];
  w->mmax = (size_t) hu[2];
  w->nfunc = (size_t) hu[3];
  w->nout = (size_t) hu[4];
  w->nrad = (size_t) hu[5];
  w->ncomp = (size_t) hu[6];
  w->w = (int) hu[7];
  w->nf1 = (size_t) hu[8];
  w->nf2 = (size_t) hu[9];
  w->nrows = (size_t) hu[10];
  w->ncols = (size_t) hu[11];
  w->nch = (size_t) hu[12];
  w->nchp = (size_t) hu[13];
  rad->rho = (size_t) hu[14];
  rad->ncheb = (size_t) hu[15];

  w->row_lo = (long) hi[0];
  w->row_hi = (long) hi[1];
  w->col_lo = (long) hi[2];
  w->col_hi = (long) hi[3];

  w->a = hd[0];
  w->rmin = hd[1];
  w->rmax = hd[2];
  w->th_lo = hd[3];
  w->th_hi = hd[4];
  w->sigma = hd[5];
  w->beta = hd[6];
  w->half = hd[7];
  w->h1 = hd[8];
  w->h2 = hd[9];
  w->inv_h1 = hd[10];
  w->inv_h2 = hd[11];
  rad->uc = hd[12];
  rad->ud = hd[13];
  w->params.eps = hd[14];
  w->params.nthreads = nthreads;
  w->params.method = POTSYNTH_METHOD_FAST;
  w->params.theta_min = w->th_lo;
  w->params.theta_max = w->th_hi;

  if (radial != NULL)
    {
      w->rfn = radial->f;
      w->rfn_params = radial->params;
    }

  /* sizes that the header must satisfy before anything is allocated from it: the
   * products must not overflow, the window must fit in the stored rows and columns
   * (window_origin() clamps to them), and a model with radial functions of the caller
   * has one layer per function */
  ok = (w->nfunc > 0 && rad->rho > 0 && rad->ncheb > 0 && w->lmin <= w->lmax && w->mmax <= w->lmax &&
        io_mul(rad->rho, rad->ncheb, &nbasis) && io_mul(nbasis, sizeof(double), &nbytes) &&
        io_mul(w->nrows, w->ncols, &ngrid) && io_mul(ngrid, w->nchp, &ngrid) &&
        io_mul(ngrid, sizeof(double), &nbytes) &&
        io_mul(w->nfunc, rad->rho, &nbytes) && w->nch == nbytes &&
        io_mul(w->nout, w->nfunc, &nbytes) && io_mul(nbytes, sizeof(double), &nbytes) &&
        w->nchp >= w->nch && w->nchp % POTSYNTH_CHAN_ALIGN == 0 && w->w >= 2 && w->w <= 32 &&
        labs(w->row_lo) < IO_MAXIDX && labs(w->row_hi) < IO_MAXIDX &&
        labs(w->col_lo) < IO_MAXIDX && labs(w->col_hi) < IO_MAXIDX &&
        w->row_lo <= w->row_hi && w->col_lo <= w->col_hi &&
        w->nrows == (size_t) (w->row_hi - w->row_lo + 1) && w->ncols == (size_t) (w->col_hi - w->col_lo + 1) &&
        w->nrows >= (size_t) w->w && w->ncols >= (size_t) w->w &&
        isfinite(w->inv_h1) && isfinite(w->inv_h2) && isfinite(w->half) &&
        (w->nrad == 0 || (rad->rho == w->nrad && rad->ncheb == 1)) &&
        w->nout >= w->nfunc && w->nout <= 2 * w->nfunc &&
        (hu[16] || w->nout == w->nfunc));

  if (ok)
    {
      rad->basis = potsynth_aligned_calloc(nbasis * sizeof(double));
      w->grid = potsynth_aligned_calloc(ngrid * sizeof(double));
      if (hu[16])
        w->omap = malloc(w->nout * w->nfunc * sizeof(double));

      ok = (rad->basis != NULL && w->grid != NULL && (!hu[16] || w->omap != NULL));
      ok = ok && (fread(rad->basis, sizeof(double), nbasis, fp) == nbasis);
      if (ok && hu[16])
        ok = (fread(w->omap, sizeof(double), w->nout * w->nfunc, fp) == w->nout * w->nfunc);
      ok = ok && (fread(w->grid, sizeof(double), ngrid, fp) == ngrid);
    }

  fclose(fp);

  if (!ok)
    {
      potsynth_free(w);
      POTSYNTH_ERROR_NULL ("the model file is damaged or too large for the memory", POTSYNTH_EINVAL);
    }

  return w;
}
