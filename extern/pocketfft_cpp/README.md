pocketfft (C++, header-only) by Martin Reinecke and Peter Bell, BSD-3-Clause
(LICENSE.md and the header's own notice).

Source: https://github.com/mreineck/pocketfft, branch cpp,
commit 4d865cba24359b05b37da4e75aee24498c3b78e7 (fetched 2026-11-25).
pocketfft_hdronly.h sha256 518a95c552aeaee352aa30122137cd9e135d23c91e2863d9272aaa758744c86d

Unmodified. potsynth uses it through src/fft.cc (C interface, src/potsynth_fft.h),
compiled with POCKETFFT_NO_MULTITHREADING (potsynth threads with OpenMP) and a
plan cache (POCKETFFT_CACHE_SIZE).
