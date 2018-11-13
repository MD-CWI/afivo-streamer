#if NDIM == 2
#define DTIMES(TXT) TXT, TXT
#define KJI_DO(lo,hi) j = lo, hi; do i = lo, hi
#define CLOSE_DO end do
#define IJK i, j
#define DIMNAME "2d"
#elif NDIM == 3
#define DTIMES(TXT) TXT, TXT, TXT
#define KJI_DO(lo,hi) k = lo, hi; do j = lo, hi; do i = lo, hi
#define CLOSE_DO end do; end do
#define IJK i, j, k
#define DIMNAME "3d"
#endif
