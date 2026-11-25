/*
 * potsynth_fft.h
 *
 * FFTs of the set-up (fast.c), from the C++ pocketfft (extern/pocketfft_cpp) through
 * src/fft.cc. Many transforms of the same length are done in one call, vectorized
 * across them; data are interleaved complex (re, im) or real, with strides in units
 * of double. Conventions (as the FFTPACK-derived C pocketfft used before):
 *
 *   forward:  F_k = fct sum_j c_j e^{-2 pi i jk/n}
 *   backward: c_j = fct sum_k F_k e^{+2 pi i jk/n}
 *   real backward, half-complex input r0, r1, i1, r2, i2, ... along the transform:
 *     x_k = fct [r0 + 2 sum_m (r_m cos - i_m sin)(2 pi mk/n) + r_{n/2} (-1)^k]
 *
 * Copyright (C) 2026 Patrick Alken and Johnathan Burchill, BSD 3-Clause License
 */

#ifndef __POTSYNTH_FFT_H__
#define __POTSYNTH_FFT_H__

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* howmany complex transforms of length n, in place: transform t starts at
   data + 2*t*dist (complex values next to each other); forward or backward */
int potsynth_fft_c2c(const size_t n, const size_t howmany, const size_t dist, const int forward,
                     double * data, const double fct);

/* howmany real backward transforms of length n from half-complex input, in place:
   element k of transform t is data[k*stride + t*dist] */
int potsynth_fft_r_backward(const size_t n, const size_t howmany, const ptrdiff_t stride,
                            const ptrdiff_t dist, double * data, const double fct);

#ifdef __cplusplus
}
#endif

#endif /* __POTSYNTH_FFT_H__ */
