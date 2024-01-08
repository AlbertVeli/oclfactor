#include <stdio.h>
#include <assert.h>
/* For CHAR_BIT */
#include <limits.h>
/* Do all mpz stuff on CPU with gmp library */
#include <gmp.h>

/* 8< -----
 *
 * Everything down to ----- >8 is from mini-gmp
 * which in turn has code from longlong.h and
 * the mpn/generic subdirectory in the gmp sources.
 *
 *   https://gmplib.org/repo/gmp/file/tip/mini-gmp
 *
 * And it is all to make minigmp_div_qr_1() work.
 * which I want to run on the GPU so it has some
 * restrictions. I have rewritten it to use no
 * malloc/free. Instead it is hardcoded to use
 * two 64-bit limbs = 128 bit numbers.
 *
 * This will currently only work on 64-bit machines
 * (CL_DEVICE_ADDRESS_BITS == 64) because mp_limb_t
 * is assumed to be 64 bits and ulong is 64 bits.
 * */

/* mp_limb_t, mp_size_t, mp_ptr, mp_srcptr
 * are defined in gmp.h
 * */

#define MINIGMP_LIMB_BITS (sizeof(mp_limb_t) * CHAR_BIT)

#define MINIGMP_LIMB_MAX (~ (mp_limb_t) 0)
#define MINIGMP_LIMB_HIGHBIT ((mp_limb_t) 1 << (MINIGMP_LIMB_BITS - 1))

#define MINIGMP_HLIMB_BIT ((mp_limb_t) 1 << (MINIGMP_LIMB_BITS / 2))
#define MINIGMP_LLIMB_MASK (MINIGMP_HLIMB_BIT - 1)

#define MINIGMP_ULONG_BITS (sizeof(unsigned long) * CHAR_BIT)
#define MINIGMP_ULONG_HIGHBIT ((unsigned long) 1 << (MINIGMP_ULONG_BITS - 1))

#define MINIGMP_ABS(x) ((x) >= 0 ? (x) : -(x))
#define MINIGMP_NEG_CAST(T,x) (-((T)((x) + 1) - 1))

#define MINIGMP_MIN(a, b) ((a) < (b) ? (a) : (b))
#define MINIGMP_MAX(a, b) ((a) > (b) ? (a) : (b))

#define minigmp_umul_ppmm(w1, w0, u, v)					\
  do {									\
    mp_limb_t __x0, __x1, __x2, __x3;					\
    unsigned __ul, __vl, __uh, __vh;					\
    mp_limb_t __u = (u), __v = (v);					\
									\
    __ul = __u & MINIGMP_LLIMB_MASK;					\
    __uh = __u >> (MINIGMP_LIMB_BITS / 2);				\
    __vl = __v & MINIGMP_LLIMB_MASK;					\
    __vh = __v >> (MINIGMP_LIMB_BITS / 2);				\
									\
    __x0 = (mp_limb_t) __ul * __vl;					\
    __x1 = (mp_limb_t) __ul * __vh;					\
    __x2 = (mp_limb_t) __uh * __vl;					\
    __x3 = (mp_limb_t) __uh * __vh;					\
									\
    __x1 += __x0 >> (MINIGMP_LIMB_BITS / 2);/* this can't give carry */	\
    __x1 += __x2;		/* but this indeed can */		\
    if (__x1 < __x2)		/* did we get it? */			\
      __x3 += MINIGMP_HLIMB_BIT; /* yes, add it in the proper pos. */	\
									\
    (w1) = __x3 + (__x1 >> (MINIGMP_LIMB_BITS / 2));			\
    (w0) = (__x1 << (MINIGMP_LIMB_BITS / 2)) + (__x0 & MINIGMP_LLIMB_MASK); \
  } while (0)

#define minigmp_clz(count, x) do {					\
    mp_limb_t __clz_x = (x);						\
    unsigned __clz_c;							\
    for (__clz_c = 0;							\
	 (__clz_x & ((mp_limb_t) 0xff << (MINIGMP_LIMB_BITS - 8))) == 0; \
	 __clz_c += 8)							\
      __clz_x <<= 8;							\
    for (; (__clz_x & MINIGMP_LIMB_HIGHBIT) == 0; __clz_c++)		\
      __clz_x <<= 1;							\
    (count) = __clz_c;							\
  } while (0)

#define minigmp_ctz(count, x) do {					\
    mp_limb_t __ctz_x = (x);						\
    unsigned __ctz_c = 0;						\
    minigmp_clz (__ctz_c, __ctz_x & - __ctz_x);				\
    (count) = MINIGMP_LIMB_BITS - 1 - __ctz_c;				\
  } while (0)

#define minigmp_add_ssaaaa(sh, sl, ah, al, bh, bl) \
  do {									\
    mp_limb_t __x;							\
    __x = (al) + (bl);							\
    (sh) = (ah) + (bh) + (__x < (al));					\
    (sl) = __x;								\
  } while (0)

#define minigmp_udiv_qrnnd_preinv(q, r, nh, nl, d, di)			\
  do {									\
    mp_limb_t _qh, _ql, _r, _mask;					\
    minigmp_umul_ppmm (_qh, _ql, (nh), (di));				\
    minigmp_add_ssaaaa (_qh, _ql, _qh, _ql, (nh) + 1, (nl));		\
    _r = (nl) - _qh * (d);						\
    _mask = -(mp_limb_t) (_r > _ql); /* both > and >= are OK */		\
    _qh += _mask;							\
    _r += _mask & (d);							\
    if (_r >= (d))							\
      {									\
	_r -= (d);							\
	_qh++;								\
      }									\
									\
    (r) = _r;								\
    (q) = _qh;								\
  } while (0)

static mp_limb_t
minigmp_invert_3by2(mp_limb_t u1, mp_limb_t u0)
{
	mp_limb_t r, p, m;
	unsigned ul, uh;
	unsigned ql, qh;

	/* First, do a 2/1 inverse. */
	/* The inverse m is defined as floor( (B^2 - 1 - u1)/u1 ), so that 0 <
	 * B^2 - (B + m) u1 <= u1 */
	assert(u1 >= MINIGMP_LIMB_HIGHBIT);

	ul = u1 & MINIGMP_LLIMB_MASK;
	uh = u1 >> (MINIGMP_LIMB_BITS / 2);

	qh = ~u1 / uh;
	r = ((~u1 - (mp_limb_t) qh * uh) << (MINIGMP_LIMB_BITS / 2)) | MINIGMP_LLIMB_MASK;

	p = (mp_limb_t) qh *ul;

	/* Adjustment steps taken from udiv_qrnnd_c */
	if (r < p) {
		qh--;
		r += u1;
		if (r >= u1)	/* i.e. we didn't get carry when adding to r */
			if (r < p) {
				qh--;
				r += u1;
			}
	}
	r -= p;

	/* Do a 3/2 division (with half limb size) */
	p = (r >> (MINIGMP_LIMB_BITS / 2)) * qh + r;
	ql = (p >> (MINIGMP_LIMB_BITS / 2)) + 1;

	/* By the 3/2 method, we don't need the high half limb. */
	r = (r << (MINIGMP_LIMB_BITS / 2)) + MINIGMP_LLIMB_MASK - ql * u1;

	if (r >= (p << (MINIGMP_LIMB_BITS / 2))) {
		ql--;
		r += u1;
	}
	m = ((mp_limb_t) qh << (MINIGMP_LIMB_BITS / 2)) + ql;
	if (r >= u1) {
		m++;
		r -= u1;
	}

	if (u0 > 0) {
		mp_limb_t th, tl;

		r = ~r;
		r += u0;
		if (r < u0) {
			m--;
			if (r >= u1) {
				m--;
				r -= u1;
			}
			r -= u1;
		}
		minigmp_umul_ppmm(th, tl, u0, m);
		r += th;
		if (r < th) {
			m--;
			m -= ((r > u1) | ((r == u1) & (tl > u0)));
		}
	}

	return m;
}

#define minigmp_invert_limb(x) minigmp_invert_3by2 ((x), 0)

struct minigmp_div_inverse {
	/* Normalization shift count. */
	unsigned shift;
	/* Normalized divisor (d0 unused for mpn_div_qr_1) */
	mp_limb_t d1, d0;
	/* Inverse, for 2/1 or 3/2. */
	mp_limb_t di;
};

static void
minigmp_div_qr_1_invert(struct minigmp_div_inverse *inv, mp_limb_t d)
{
	unsigned shift;

	assert(d > 0);
	minigmp_clz(shift, d);
	inv->shift = shift;
	inv->d1 = d << shift;
	inv->di = minigmp_invert_limb(inv->d1);
}

static mp_limb_t
minigmp_div_qr_1_preinv(mp_ptr qp, mp_srcptr np, mp_size_t nn, const struct minigmp_div_inverse *inv)
{
	mp_limb_t d, di;
	mp_limb_t r;
	mp_limb_t tmp[2];
	mp_ptr tp = tmp;

	if (inv->shift > 0) {
		r = mpn_lshift(tp, np, nn, inv->shift);
		np = tp;
	} else
		r = 0;

	d = inv->d1;
	di = inv->di;
	while (--nn >= 0) {
		mp_limb_t q;

		minigmp_udiv_qrnnd_preinv(q, r, r, np[nn], d, di);
		if (qp)
			qp[nn] = q;
	}

	return r >> inv->shift;
}

/* This is the only function called from our code */
mp_limb_t
minigmp_div_qr_1(mp_ptr qp, mp_srcptr np, mp_size_t nn, mp_limb_t d)
{
	assert(d > 0);

	/* Special case for powers of two. */
	if ((d & (d - 1)) == 0) {
		mp_limb_t r = np[0] & (d - 1);

		if (qp) {
			if (d <= 1) {
				mpn_copyi(qp, np, nn);
			} else {
				unsigned shift;

				minigmp_ctz(shift, d);
				mpn_rshift(qp, np, nn, shift);
			}
		}
		return r;
	} else {
		struct minigmp_div_inverse inv;

		minigmp_div_qr_1_invert(&inv, d);
		return minigmp_div_qr_1_preinv(qp, np, nn, &inv);
	}
}

/* ----- >8
 *
 * End of copied code from mini-gmp.
 *
 * All defines and functions between scissors are
 * there to make minigmp_div_qr_1() work.
 */


/* Since we are using the mpn interface
 * we need some macros to copy 2 limbs.
 */
#define MPZ_SET_MPN_2(mpz, mpn)      \
  do {                               \
    mpz->_mp_d[1] = mpn[1];          \
    mpz->_mp_d[0] = mpn[0];          \
    if (mpn[1] == 0 && mpn[0] == 0)  \
        mpz->_mp_size = 0;           \
    else if (mpn[1] == 0)            \
	  mpz->_mp_size = 1;         \
  } while (0)

#define MPZ_TO_MPN_2(mpn, mpz)       \
  do {                               \
    mpn[1] = mpz->_mp_d[1];          \
    mpn[0] = mpz->_mp_d[0];          \
  } while (0)

#define MPN_TO_MPN_2(mpndst, mpnsrc) \
  do {                               \
    mpndst[1] = mpnsrc[1];           \
    mpndst[0] = mpnsrc[0];           \
  } while (0)

#define IS_PRIME(n) mpz_probab_prime_p(n, 2)

typedef unsigned long ulong;

int main(void)
{
	mpz_t n;
	ulong d, r;
	mp_limb_t n_limbs[2];
	mp_limb_t q_limbs[2];
	mp_size_t n_len = 2;

	/* TODO: Read input from argv */
	mpz_init_set_ui(n, (ulong)(2 * 435505820298201979));
	mpz_mul_ui(n, n, 251);

	if (n->_mp_size > 2) {
		fprintf(stderr,
			"n has %d limbs.\n"
			"Currently only up to 2 limbs are supported.\n", n->_mp_size);
		mpz_clear(n);
		return 1;
	}

	/* Copy from mpz to mpn interface */
	MPZ_TO_MPN_2(n_limbs, n);

	/* TODO: Don't try all numbers, do a sieve */
	for (d = 2; d < MINIGMP_LIMB_MAX; d++) {

		/* This should go to the GPU (in ranges, to reduce overhead) */
		r = minigmp_div_qr_1(q_limbs, n_limbs, n_len, d);

		/* Back on CPU */
		if (r == 0) {
			printf("%lu\n", d);
			/* n = q */
			MPZ_SET_MPN_2(n, q_limbs);
			MPN_TO_MPN_2(n_limbs, q_limbs);
			if (IS_PRIME(n)) {
				mpz_out_str(stdout, 10, n);
				puts("");
				break;
			}
		}
	}

	mpz_clear(n);

	return 0;
}
