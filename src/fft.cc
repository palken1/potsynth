/*
 * fft.cc
 *
 * C interface (potsynth_fft.h) to the C++ pocketfft for the set-up: many transforms
 * of one length per call, which pocketfft vectorizes across. Exceptions (memory,
 * invalid sizes) are turned into error codes.
 *
 * Copyright (C) 2026 Patrick Alken and Johnathan Burchill, BSD 3-Clause License
 */

#include <complex>
#include <cstddef>
#include <new>

#include "pocketfft_hdronly.h"
#include "potsynth_fft.h"

namespace
{
  const int FFT_SUCCESS = 0, FFT_FAILED = 2; /* POTSYNTH_SUCCESS, POTSYNTH_ENOMEM */
}

extern "C" int
potsynth_fft_c2c(const size_t n, const size_t howmany, const size_t dist, const int forward,
                 double * data, const double fct)
{
  try
    {
      std::complex<double> * c = reinterpret_cast<std::complex<double> *>(data);
      const pocketfft::shape_t shape{ howmany, n };
      const pocketfft::stride_t stride{ static_cast<ptrdiff_t>(dist * sizeof(std::complex<double>)),
                                        static_cast<ptrdiff_t>(sizeof(std::complex<double>)) };

      pocketfft::c2c<double>(shape, stride, stride, { 1 }, forward != 0, c, c, fct, 1);
      return FFT_SUCCESS;
    }
  catch (...)
    {
      return FFT_FAILED;
    }
}

extern "C" int
potsynth_fft_r_backward(const size_t n, const size_t howmany, const ptrdiff_t stride,
                        const ptrdiff_t dist, double * data, const double fct)
{
  try
    {
      const pocketfft::shape_t shape{ n, howmany };
      const pocketfft::stride_t strides{ stride * static_cast<ptrdiff_t>(sizeof(double)),
                                         dist * static_cast<ptrdiff_t>(sizeof(double)) };

      /* half-complex (FFTPACK order) to real, no sign change: as rfft_backward of the C version */
      pocketfft::r2r_fftpack<double>(shape, strides, strides, { 0 }, false, false, data, data, fct, 1);
      return FFT_SUCCESS;
    }
  catch (...)
    {
      return FFT_FAILED;
    }
}
