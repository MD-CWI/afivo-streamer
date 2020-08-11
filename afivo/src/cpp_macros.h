#ifdef __GFORTRAN__
#define PASTE(a) a
#define CONCAT(a,b) PASTE(a)b
#else
#define PASTE(a) a ## b
#define CONCAT(a,b) PASTE(a,b)
#endif

#if NDIM == 1
#define DTIMES(TXT) TXT
#define DINDEX(TXT) TXT(1)
#define DSLICE(lo,hi) lo(1):hi(1)
#define KJI_DO(lo,hi) i = lo, hi
#define CLOSE_DO
#define IJK i
#define IJK_(s) CONCAT(i_,s)
#define DIMNAME "1d"
#elif NDIM == 2
#define DTIMES(TXT) TXT, TXT
#define DINDEX(TXT) TXT(1), TXT(2)
#define DSLICE(lo,hi) lo(1):hi(1), lo(2):hi(2)
#define KJI_DO(lo,hi) j = lo, hi; do i = lo, hi
#define CLOSE_DO end do
#define IJK i, j
#define IJK_(s) CONCAT(i_,s), CONCAT(j_,s)
#define DIMNAME "2d"
#elif NDIM == 3
#define DTIMES(TXT) TXT, TXT, TXT
#define DINDEX(TXT) TXT(1), TXT(2), TXT(3)
#define DSLICE(lo,hi) lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)
#define KJI_DO(lo,hi) k = lo, hi; do j = lo, hi; do i = lo, hi
#define CLOSE_DO end do; end do
#define IJK i, j, k
#define IJK_(s) CONCAT(i_,s), CONCAT(j_,s), CONCAT(k_,s)
#define DIMNAME "3d"
#endif
