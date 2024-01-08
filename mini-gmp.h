/* mini-gmp, a minimalistic implementation of a GNU GMP subset.

Copyright 2011-2015 Free Software Foundation, Inc.

This file is part of the GNU MP Library.

The GNU MP Library is free software; you can redistribute it and/or modify
it under the terms of either:

  * the GNU Lesser General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your
    option) any later version.

or

  * the GNU General Public License as published by the Free Software
    Foundation; either version 2 of the License, or (at your option) any
    later version.

or both in parallel, as here.

The GNU MP Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received copies of the GNU General Public License and the
GNU Lesser General Public License along with the GNU MP Library.  If not,
see https://www.gnu.org/licenses/.  */

/* About mini-gmp: This is a minimal implementation of a subset of the
   GMP interface. It is intended for inclusion into applications which
   have modest bignums needs, as a fallback when the real GMP library
   is not installed.

   This file defines the public interface. */

/* NOTE: This is not the original GMP-file.
 * This version was further stripped by Albert to
 * create a truly minimal standalone GMP. This version
 * lacks many features even of the original mini-gmp */


#ifndef __MINI_GMP_H__
#define __MINI_GMP_H__

/* For size_t */
# include <stddef.h>
/* CHAR_BIT */
# include <limits.h>

#if defined (__cplusplus)
extern "C" {
#endif

typedef unsigned long mp_limb_t;
typedef long mp_size_t;
typedef unsigned long mp_bitcnt_t;

typedef mp_limb_t *mp_ptr;
typedef const mp_limb_t *mp_srcptr;

typedef struct {
	int _mp_alloc;	/* Number of *limbs* allocated and pointed
			   to by the _mp_d field.  */
	int _mp_size;	/* abs(_mp_size) is the number of limbs the
			   last field points to.  If _mp_size is
			   negative this is a negative number.  */
	mp_limb_t *_mp_d;	/* Pointer to the limbs.  */
} __mpz_struct;

typedef __mpz_struct mpz_t[1];

typedef __mpz_struct *mpz_ptr;
typedef const __mpz_struct *mpz_srcptr;

/* Macros */
#define GMP_LIMB_BITS (sizeof(mp_limb_t) * CHAR_BIT)

#define GMP_LIMB_MAX (~ (mp_limb_t) 0)
#define GMP_LIMB_HIGHBIT ((mp_limb_t) 1 << (GMP_LIMB_BITS - 1))

#define GMP_HLIMB_BIT ((mp_limb_t) 1 << (GMP_LIMB_BITS / 2))
#define GMP_LLIMB_MASK (GMP_HLIMB_BIT - 1)

#define GMP_ULONG_BITS (sizeof(unsigned long) * CHAR_BIT)
#define GMP_ULONG_HIGHBIT ((unsigned long) 1 << (GMP_ULONG_BITS - 1))

#define GMP_ABS(x) ((x) >= 0 ? (x) : -(x))
#define GMP_NEG_CAST(T,x) (-((T)((x) + 1) - 1))

#define GMP_MIN(a, b) ((a) < (b) ? (a) : (b))
#define GMP_MAX(a, b) ((a) > (b) ? (a) : (b))

#define gmp_assert_nocarry(x) do { \
    mp_limb_t __cy = x;		   \
    assert (__cy == 0);		   \
  } while (0)

#define gmp_clz(count, x) do {						\
    mp_limb_t __clz_x = (x);						\
    unsigned __clz_c;							\
    for (__clz_c = 0;							\
	 (__clz_x & ((mp_limb_t) 0xff << (GMP_LIMB_BITS - 8))) == 0;	\
	 __clz_c += 8)							\
      __clz_x <<= 8;							\
    for (; (__clz_x & GMP_LIMB_HIGHBIT) == 0; __clz_c++)		\
      __clz_x <<= 1;							\
    (count) = __clz_c;							\
  } while (0)

#define gmp_ctz(count, x) do {						\
    mp_limb_t __ctz_x = (x);						\
    unsigned __ctz_c = 0;						\
    gmp_clz (__ctz_c, __ctz_x & - __ctz_x);				\
    (count) = GMP_LIMB_BITS - 1 - __ctz_c;				\
  } while (0)

#define gmp_add_ssaaaa(sh, sl, ah, al, bh, bl) \
  do {									\
    mp_limb_t __x;							\
    __x = (al) + (bl);							\
    (sh) = (ah) + (bh) + (__x < (al));					\
    (sl) = __x;								\
  } while (0)

#define gmp_sub_ddmmss(sh, sl, ah, al, bh, bl) \
  do {									\
    mp_limb_t __x;							\
    __x = (al) - (bl);							\
    (sh) = (ah) - (bh) - ((al) < (bl));					\
    (sl) = __x;								\
  } while (0)

#define gmp_umul_ppmm(w1, w0, u, v)					\
  do {									\
    mp_limb_t __x0, __x1, __x2, __x3;					\
    unsigned __ul, __vl, __uh, __vh;					\
    mp_limb_t __u = (u), __v = (v);					\
									\
    __ul = __u & GMP_LLIMB_MASK;					\
    __uh = __u >> (GMP_LIMB_BITS / 2);					\
    __vl = __v & GMP_LLIMB_MASK;					\
    __vh = __v >> (GMP_LIMB_BITS / 2);					\
									\
    __x0 = (mp_limb_t) __ul * __vl;					\
    __x1 = (mp_limb_t) __ul * __vh;					\
    __x2 = (mp_limb_t) __uh * __vl;					\
    __x3 = (mp_limb_t) __uh * __vh;					\
									\
    __x1 += __x0 >> (GMP_LIMB_BITS / 2);/* this can't give carry */	\
    __x1 += __x2;		/* but this indeed can */		\
    if (__x1 < __x2)		/* did we get it? */			\
      __x3 += GMP_HLIMB_BIT;	/* yes, add it in the proper pos. */	\
									\
    (w1) = __x3 + (__x1 >> (GMP_LIMB_BITS / 2));			\
    (w0) = (__x1 << (GMP_LIMB_BITS / 2)) + (__x0 & GMP_LLIMB_MASK);	\
  } while (0)

#define gmp_udiv_qrnnd_preinv(q, r, nh, nl, d, di)			\
  do {									\
    mp_limb_t _qh, _ql, _r, _mask;					\
    gmp_umul_ppmm (_qh, _ql, (nh), (di));				\
    gmp_add_ssaaaa (_qh, _ql, _qh, _ql, (nh) + 1, (nl));		\
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

#define gmp_udiv_qr_3by2(q, r1, r0, n2, n1, n0, d1, d0, dinv)		\
  do {									\
    mp_limb_t _q0, _t1, _t0, _mask;					\
    gmp_umul_ppmm ((q), _q0, (n2), (dinv));				\
    gmp_add_ssaaaa ((q), _q0, (q), _q0, (n2), (n1));			\
									\
    /* Compute the two most significant limbs of n - q'd */		\
    (r1) = (n1) - (d1) * (q);						\
    gmp_sub_ddmmss ((r1), (r0), (r1), (n0), (d1), (d0));		\
    gmp_umul_ppmm (_t1, _t0, (d0), (q));				\
    gmp_sub_ddmmss ((r1), (r0), (r1), (r0), _t1, _t0);			\
    (q)++;								\
									\
    /* Conditionally adjust q and the remainders */			\
    _mask = - (mp_limb_t) ((r1) >= _q0);				\
    (q) += _mask;							\
    gmp_add_ssaaaa ((r1), (r0), (r1), (r0), _mask & (d1), _mask & (d0)); \
    if ((r1) >= (d1))							\
      {									\
	if ((r1) > (d1) || (r0) >= (d0))				\
	  {								\
	    (q)++;							\
	    gmp_sub_ddmmss ((r1), (r0), (r1), (r0), (d1), (d0));	\
	  }								\
      }									\
  } while (0)

/* Swap macros. */
#define MP_LIMB_T_SWAP(x, y)						\
  do {									\
    mp_limb_t __mp_limb_t_swap__tmp = (x);				\
    (x) = (y);								\
    (y) = __mp_limb_t_swap__tmp;					\
  } while (0)
#define MP_SIZE_T_SWAP(x, y)						\
  do {									\
    mp_size_t __mp_size_t_swap__tmp = (x);				\
    (x) = (y);								\
    (y) = __mp_size_t_swap__tmp;					\
  } while (0)
#define MP_BITCNT_T_SWAP(x,y)			\
  do {						\
    mp_bitcnt_t __mp_bitcnt_t_swap__tmp = (x);	\
    (x) = (y);					\
    (y) = __mp_bitcnt_t_swap__tmp;		\
  } while (0)
#define MP_PTR_SWAP(x, y)						\
  do {									\
    mp_ptr __mp_ptr_swap__tmp = (x);					\
    (x) = (y);								\
    (y) = __mp_ptr_swap__tmp;						\
  } while (0)
#define MP_SRCPTR_SWAP(x, y)						\
  do {									\
    mp_srcptr __mp_srcptr_swap__tmp = (x);				\
    (x) = (y);								\
    (y) = __mp_srcptr_swap__tmp;					\
  } while (0)

#define MPN_PTR_SWAP(xp,xs, yp,ys)					\
  do {									\
    MP_PTR_SWAP (xp, yp);						\
    MP_SIZE_T_SWAP (xs, ys);						\
  } while(0)
#define MPN_SRCPTR_SWAP(xp,xs, yp,ys)					\
  do {									\
    MP_SRCPTR_SWAP (xp, yp);						\
    MP_SIZE_T_SWAP (xs, ys);						\
  } while(0)

#define MPZ_PTR_SWAP(x, y)						\
  do {									\
    mpz_ptr __mpz_ptr_swap__tmp = (x);					\
    (x) = (y);								\
    (y) = __mpz_ptr_swap__tmp;						\
  } while (0)
#define MPZ_SRCPTR_SWAP(x, y)						\
  do {									\
    mpz_srcptr __mpz_srcptr_swap__tmp = (x);				\
    (x) = (y);								\
    (y) = __mpz_srcptr_swap__tmp;					\
  } while (0)

/* MPZ division */
enum mpz_div_round_mode { GMP_DIV_FLOOR, GMP_DIV_CEIL, GMP_DIV_TRUNC };

/* Utlyfta av albert */
mp_limb_t mpn_div_qr_1(mp_ptr qp, mp_srcptr np, mp_size_t nn, mp_limb_t d);
unsigned long mpz_div_qr_ui(mpz_t q, mpz_t r, const mpz_t n, unsigned long d, enum mpz_div_round_mode mode);

extern const int mp_bits_per_limb;

void mpn_copyi(mp_ptr, mp_srcptr, mp_size_t);
void mpn_copyd(mp_ptr, mp_srcptr, mp_size_t);
void mpn_zero(mp_ptr, mp_size_t);

int mpn_cmp(mp_srcptr, mp_srcptr, mp_size_t);
int mpn_zero_p(mp_srcptr, mp_size_t);

mp_limb_t mpn_add_1(mp_ptr, mp_srcptr, mp_size_t, mp_limb_t);
mp_limb_t mpn_add_n(mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);
mp_limb_t mpn_add(mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t);

mp_limb_t mpn_sub_1(mp_ptr, mp_srcptr, mp_size_t, mp_limb_t);
mp_limb_t mpn_sub_n(mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);
mp_limb_t mpn_sub(mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t);

mp_limb_t mpn_mul_1(mp_ptr, mp_srcptr, mp_size_t, mp_limb_t);
mp_limb_t mpn_addmul_1(mp_ptr, mp_srcptr, mp_size_t, mp_limb_t);
mp_limb_t mpn_submul_1(mp_ptr, mp_srcptr, mp_size_t, mp_limb_t);

mp_limb_t mpn_mul(mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t);
void mpn_mul_n(mp_ptr, mp_srcptr, mp_srcptr, mp_size_t);

mp_limb_t mpn_lshift(mp_ptr, mp_srcptr, mp_size_t, unsigned int);
mp_limb_t mpn_rshift(mp_ptr, mp_srcptr, mp_size_t, unsigned int);

mp_bitcnt_t mpn_scan0(mp_srcptr, mp_bitcnt_t);
mp_bitcnt_t mpn_scan1(mp_srcptr, mp_bitcnt_t);

mp_limb_t mpn_invert_3by2(mp_limb_t, mp_limb_t);
#define mpn_invert_limb(x) mpn_invert_3by2 ((x), 0)

void mpz_init(mpz_t);
void mpz_init2(mpz_t, mp_bitcnt_t);
void mpz_clear(mpz_t);

#define mpz_odd_p(z)   (((z)->_mp_size != 0) & (int) (z)->_mp_d[0])
#define mpz_even_p(z)  (! mpz_odd_p (z))

int mpz_sgn(const mpz_t);
int mpz_cmp_si(const mpz_t, long);
int mpz_cmp_ui(const mpz_t, unsigned long);
int mpz_cmp(const mpz_t, const mpz_t);
int mpz_cmpabs_ui(const mpz_t, unsigned long);
int mpz_cmpabs(const mpz_t, const mpz_t);
int mpz_cmp_d(const mpz_t, double);
int mpz_cmpabs_d(const mpz_t, double);

void mpz_abs(mpz_t, const mpz_t);
void mpz_neg(mpz_t, const mpz_t);
void mpz_swap(mpz_t, mpz_t);

void mpz_add_ui(mpz_t, const mpz_t, unsigned long);
void mpz_add(mpz_t, const mpz_t, const mpz_t);
void mpz_sub_ui(mpz_t, const mpz_t, unsigned long);
void mpz_ui_sub(mpz_t, unsigned long, const mpz_t);
void mpz_sub(mpz_t, const mpz_t, const mpz_t);

void mpz_mul_si(mpz_t, const mpz_t, long int);
void mpz_mul_ui(mpz_t, const mpz_t, unsigned long int);
void mpz_mul(mpz_t, const mpz_t, const mpz_t);
void mpz_mul_2exp(mpz_t, const mpz_t, mp_bitcnt_t);
void mpz_addmul_ui(mpz_t, const mpz_t, unsigned long int);
void mpz_addmul(mpz_t, const mpz_t, const mpz_t);
void mpz_submul_ui(mpz_t, const mpz_t, unsigned long int);
void mpz_submul(mpz_t, const mpz_t, const mpz_t);

void mpz_cdiv_qr(mpz_t, mpz_t, const mpz_t, const mpz_t);
void mpz_fdiv_qr(mpz_t, mpz_t, const mpz_t, const mpz_t);
void mpz_tdiv_qr(mpz_t, mpz_t, const mpz_t, const mpz_t);
void mpz_cdiv_q(mpz_t, const mpz_t, const mpz_t);
void mpz_fdiv_q(mpz_t, const mpz_t, const mpz_t);
void mpz_tdiv_q(mpz_t, const mpz_t, const mpz_t);
void mpz_cdiv_r(mpz_t, const mpz_t, const mpz_t);
void mpz_fdiv_r(mpz_t, const mpz_t, const mpz_t);
void mpz_tdiv_r(mpz_t, const mpz_t, const mpz_t);

void mpz_cdiv_q_2exp(mpz_t, const mpz_t, mp_bitcnt_t);
void mpz_fdiv_q_2exp(mpz_t, const mpz_t, mp_bitcnt_t);
void mpz_tdiv_q_2exp(mpz_t, const mpz_t, mp_bitcnt_t);
void mpz_cdiv_r_2exp(mpz_t, const mpz_t, mp_bitcnt_t);
void mpz_fdiv_r_2exp(mpz_t, const mpz_t, mp_bitcnt_t);
void mpz_tdiv_r_2exp(mpz_t, const mpz_t, mp_bitcnt_t);

void mpz_mod(mpz_t, const mpz_t, const mpz_t);

void mpz_divexact(mpz_t, const mpz_t, const mpz_t);

int mpz_divisible_p(const mpz_t, const mpz_t);
int mpz_congruent_p(const mpz_t, const mpz_t, const mpz_t);

unsigned long mpz_cdiv_qr_ui(mpz_t, mpz_t, const mpz_t, unsigned long);
unsigned long mpz_fdiv_qr_ui(mpz_t, mpz_t, const mpz_t, unsigned long);
unsigned long mpz_tdiv_qr_ui(mpz_t, mpz_t, const mpz_t, unsigned long);
unsigned long mpz_cdiv_q_ui(mpz_t, const mpz_t, unsigned long);
unsigned long mpz_fdiv_q_ui(mpz_t, const mpz_t, unsigned long);
unsigned long mpz_tdiv_q_ui(mpz_t, const mpz_t, unsigned long);
unsigned long mpz_cdiv_r_ui(mpz_t, const mpz_t, unsigned long);
unsigned long mpz_fdiv_r_ui(mpz_t, const mpz_t, unsigned long);
unsigned long mpz_tdiv_r_ui(mpz_t, const mpz_t, unsigned long);
unsigned long mpz_cdiv_ui(const mpz_t, unsigned long);
unsigned long mpz_fdiv_ui(const mpz_t, unsigned long);
unsigned long mpz_tdiv_ui(const mpz_t, unsigned long);

unsigned long mpz_mod_ui(mpz_t, const mpz_t, unsigned long);

void mpz_divexact_ui(mpz_t, const mpz_t, unsigned long);

int mpz_divisible_ui_p(const mpz_t, unsigned long);

unsigned long mpz_gcd_ui(mpz_t, const mpz_t, unsigned long);
void mpz_gcd(mpz_t, const mpz_t, const mpz_t);
void mpz_gcdext(mpz_t, mpz_t, mpz_t, const mpz_t, const mpz_t);
void mpz_lcm_ui(mpz_t, const mpz_t, unsigned long);
void mpz_lcm(mpz_t, const mpz_t, const mpz_t);
int mpz_invert(mpz_t, const mpz_t, const mpz_t);

void mpz_sqrtrem(mpz_t, mpz_t, const mpz_t);
void mpz_sqrt(mpz_t, const mpz_t);
int mpz_perfect_square_p(const mpz_t);

void mpz_pow_ui(mpz_t, const mpz_t, unsigned long);
void mpz_ui_pow_ui(mpz_t, unsigned long, unsigned long);
void mpz_powm(mpz_t, const mpz_t, const mpz_t, const mpz_t);
void mpz_powm_ui(mpz_t, const mpz_t, unsigned long, const mpz_t);

void mpz_rootrem(mpz_t, mpz_t, const mpz_t, unsigned long);
int mpz_root(mpz_t, const mpz_t, unsigned long);

void mpz_fac_ui(mpz_t, unsigned long);
void mpz_bin_uiui(mpz_t, unsigned long, unsigned long);

int mpz_probab_prime_p(const mpz_t, int);

int mpz_tstbit(const mpz_t, mp_bitcnt_t);
void mpz_setbit(mpz_t, mp_bitcnt_t);
void mpz_clrbit(mpz_t, mp_bitcnt_t);
void mpz_combit(mpz_t, mp_bitcnt_t);

void mpz_com(mpz_t, const mpz_t);
void mpz_and(mpz_t, const mpz_t, const mpz_t);
void mpz_ior(mpz_t, const mpz_t, const mpz_t);
void mpz_xor(mpz_t, const mpz_t, const mpz_t);

mp_bitcnt_t mpz_popcount(const mpz_t);
mp_bitcnt_t mpz_hamdist(const mpz_t, const mpz_t);
mp_bitcnt_t mpz_scan0(const mpz_t, mp_bitcnt_t);
mp_bitcnt_t mpz_scan1(const mpz_t, mp_bitcnt_t);

int mpz_fits_slong_p(const mpz_t);
int mpz_fits_ulong_p(const mpz_t);
long int mpz_get_si(const mpz_t);
unsigned long int mpz_get_ui(const mpz_t);
double mpz_get_d(const mpz_t);
size_t mpz_size(const mpz_t);
mp_limb_t mpz_getlimbn(const mpz_t, mp_size_t);

void mpz_realloc2(mpz_t, mp_bitcnt_t);
mp_srcptr mpz_limbs_read(mpz_srcptr);
mp_ptr mpz_limbs_modify(mpz_t, mp_size_t);
mp_ptr mpz_limbs_write(mpz_t, mp_size_t);
void mpz_limbs_finish(mpz_t, mp_size_t);
mpz_srcptr mpz_roinit_n(mpz_t, mp_srcptr, mp_size_t);

#define MPZ_ROINIT_N(xp, xs) {{0, (xs),(xp) }}

void mpz_set_si(mpz_t, signed long int);
void mpz_set_ui(mpz_t, unsigned long int);
void mpz_set(mpz_t, const mpz_t);
void mpz_set_d(mpz_t, double);

void mpz_init_set_si(mpz_t, signed long int);
void mpz_init_set_ui(mpz_t, unsigned long int);
void mpz_init_set(mpz_t, const mpz_t);
void mpz_init_set_d(mpz_t, double);

size_t mpz_sizeinbase(const mpz_t, int);

#define __gmpz_add mpz_add
#define __gmpz_clear mpz_clear
#define __gmpz_cmp mpz_cmp
#define __gmpz_com mpz_com
#define __gmpz_get_si mpz_get_si
#define __gmpz_get_ui mpz_get_ui
#define __gmpz_init mpz_init
#define __gmpz_init_set mpz_init_set
#define __gmpz_init_set_si mpz_init_set_si
#define __gmpz_init_set_ui mpz_init_set_ui
#define __gmpz_mul mpz_mul
#define __gmpz_neg mpz_neg
#define __gmpz_set mpz_set
#define __gmpz_set_si mpz_set_si
#define __gmpz_set_ui mpz_set_ui
#define __gmpz_sub mpz_sub
#define __gmpz_tdiv_q mpz_tdiv_q
#define __gmpz_tdiv_r mpz_tdiv_r

#if defined (__cplusplus)
}
#endif
#endif				/* __MINI_GMP_H__ */
