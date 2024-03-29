/* mini-gmp, a minimalistic implementation of a GNU GMP subset.

   Contributed to the GNU project by Niels Möller

Copyright 1991-1997, 1999-2015 Free Software Foundation, Inc.

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

/* NOTE: This is not the original GMP-file.
 * This version was further stripped by Albert to
 * create a truly minimal standalone GMP. This version
 * lacks many features even of the original mini-gmp */

/* NOTE: All functions in this file which are not declared in
   mini-gmp.h are internal, and are not intended to be compatible
   neither with GMP nor with future versions of mini-gmp. */

/* Much of the material copied from GMP files, including: gmp-impl.h,
   longlong.h, mpn/generic/add_n.c, mpn/generic/addmul_1.c,
   mpn/generic/lshift.c, mpn/generic/mul_1.c,
   mpn/generic/mul_basecase.c, mpn/generic/rshift.c,
   mpn/generic/sbpi1_div_qr.c, mpn/generic/sub_n.c,
   mpn/generic/submul_1.c. */

# include <assert.h>
# include <ctype.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

#include "mini-gmp.h"


const int mp_bits_per_limb = GMP_LIMB_BITS;

/* Memory allocation and other helper functions. */
static void gmp_die(const char *msg)
{
	fprintf(stderr, "%s\n", msg);
	abort();
}

static void *gmp_default_alloc(size_t size)
{
	void *p;

	assert(size > 0);

	p = malloc(size);
	if (!p)
		gmp_die("gmp_default_alloc: Virtual memory exhausted.");

	return p;
}

static void *gmp_default_realloc(void *old, size_t old_size __attribute__((unused)), size_t new_size)
{
	void *p;

	p = realloc(old, new_size);

	if (!p)
		gmp_die("gmp_default_realloc: Virtual memory exhausted.");

	return p;
}

static void gmp_default_free(void *p, size_t size __attribute__((unused)))
{
	free(p);
}

static void *(*gmp_allocate_func)(size_t) = gmp_default_alloc;
static void *(*gmp_reallocate_func)(void *, size_t, size_t) = gmp_default_realloc;
static void (*gmp_free_func)(void *, size_t) = gmp_default_free;

#define gmp_xalloc(size) ((*gmp_allocate_func)((size)))
#define gmp_free(p) ((*gmp_free_func) ((p), 0))

static mp_ptr gmp_xalloc_limbs(mp_size_t size)
{
	return (mp_ptr) gmp_xalloc(size * sizeof(mp_limb_t));
}

static mp_ptr gmp_xrealloc_limbs(mp_ptr old, mp_size_t size)
{
	assert(size > 0);
	return (mp_ptr) (*gmp_reallocate_func) (old, 0, size * sizeof(mp_limb_t));
}


/* MPN interface */

void mpn_copyi(mp_ptr d, mp_srcptr s, mp_size_t n)
{
	mp_size_t i;

	for (i = 0; i < n; i++)
		d[i] = s[i];
}

void mpn_copyd(mp_ptr d, mp_srcptr s, mp_size_t n)
{
	while (--n >= 0)
		d[n] = s[n];
}

int mpn_cmp(mp_srcptr ap, mp_srcptr bp, mp_size_t n)
{
	while (--n >= 0) {
		if (ap[n] != bp[n])
			return ap[n] > bp[n] ? 1 : -1;
	}
	return 0;
}

static int mpn_cmp4(mp_srcptr ap, mp_size_t an, mp_srcptr bp, mp_size_t bn)
{
	if (an != bn)
		return an < bn ? -1 : 1;
	else
		return mpn_cmp(ap, bp, an);
}

static mp_size_t mpn_normalized_size(mp_srcptr xp, mp_size_t n)
{
	while (n > 0 && xp[n - 1] == 0)
		--n;
	return n;
}

int mpn_zero_p(mp_srcptr rp, mp_size_t n)
{
	return mpn_normalized_size(rp, n) == 0;
}

void mpn_zero(mp_ptr rp, mp_size_t n)
{
	while (--n >= 0)
		rp[n] = 0;
}

mp_limb_t mpn_add_1(mp_ptr rp, mp_srcptr ap, mp_size_t n, mp_limb_t b)
{
	mp_size_t i;

	assert(n > 0);
	i = 0;
	do {
		mp_limb_t r = ap[i] + b;

		/* Carry out */
		b = (r < b);
		rp[i] = r;
	}
	while (++i < n);

	return b;
}

mp_limb_t mpn_add_n(mp_ptr rp, mp_srcptr ap, mp_srcptr bp, mp_size_t n)
{
	mp_size_t i;
	mp_limb_t cy;

	for (i = 0, cy = 0; i < n; i++) {
		mp_limb_t a, b, r;

		a = ap[i];
		b = bp[i];
		r = a + cy;
		cy = (r < cy);
		r += b;
		cy += (r < b);
		rp[i] = r;
	}
	return cy;
}

mp_limb_t mpn_add(mp_ptr rp, mp_srcptr ap, mp_size_t an, mp_srcptr bp, mp_size_t bn)
{
	mp_limb_t cy;

	assert(an >= bn);

	cy = mpn_add_n(rp, ap, bp, bn);
	if (an > bn)
		cy = mpn_add_1(rp + bn, ap + bn, an - bn, cy);
	return cy;
}

mp_limb_t mpn_sub_1(mp_ptr rp, mp_srcptr ap, mp_size_t n, mp_limb_t b)
{
	mp_size_t i;

	assert(n > 0);

	i = 0;
	do {
		mp_limb_t a = ap[i];

		/* Carry out */
		mp_limb_t cy = a < b;;
		rp[i] = a - b;
		b = cy;
	}
	while (++i < n);

	return b;
}

mp_limb_t mpn_sub_n(mp_ptr rp, mp_srcptr ap, mp_srcptr bp, mp_size_t n)
{
	mp_size_t i;
	mp_limb_t cy;

	for (i = 0, cy = 0; i < n; i++) {
		mp_limb_t a, b;

		a = ap[i];
		b = bp[i];
		b += cy;
		cy = (b < cy);
		cy += (a < b);
		rp[i] = a - b;
	}
	return cy;
}

mp_limb_t mpn_sub(mp_ptr rp, mp_srcptr ap, mp_size_t an, mp_srcptr bp, mp_size_t bn)
{
	mp_limb_t cy;

	assert(an >= bn);

	cy = mpn_sub_n(rp, ap, bp, bn);
	if (an > bn)
		cy = mpn_sub_1(rp + bn, ap + bn, an - bn, cy);
	return cy;
}

mp_limb_t mpn_mul_1(mp_ptr rp, mp_srcptr up, mp_size_t n, mp_limb_t vl)
{
	mp_limb_t ul, cl, hpl, lpl;

	assert(n >= 1);

	cl = 0;
	do {
		ul = *up++;
		gmp_umul_ppmm(hpl, lpl, ul, vl);

		lpl += cl;
		cl = (lpl < cl) + hpl;

		*rp++ = lpl;
	}
	while (--n != 0);

	return cl;
}

mp_limb_t mpn_addmul_1(mp_ptr rp, mp_srcptr up, mp_size_t n, mp_limb_t vl)
{
	mp_limb_t ul, cl, hpl, lpl, rl;

	assert(n >= 1);

	cl = 0;
	do {
		ul = *up++;
		gmp_umul_ppmm(hpl, lpl, ul, vl);

		lpl += cl;
		cl = (lpl < cl) + hpl;

		rl = *rp;
		lpl = rl + lpl;
		cl += lpl < rl;
		*rp++ = lpl;
	}
	while (--n != 0);

	return cl;
}

mp_limb_t mpn_submul_1(mp_ptr rp, mp_srcptr up, mp_size_t n, mp_limb_t vl)
{
	mp_limb_t ul, cl, hpl, lpl, rl;

	assert(n >= 1);

	cl = 0;
	do {
		ul = *up++;
		gmp_umul_ppmm(hpl, lpl, ul, vl);

		lpl += cl;
		cl = (lpl < cl) + hpl;

		rl = *rp;
		lpl = rl - lpl;
		cl += lpl > rl;
		*rp++ = lpl;
	}
	while (--n != 0);

	return cl;
}

mp_limb_t mpn_mul(mp_ptr rp, mp_srcptr up, mp_size_t un, mp_srcptr vp, mp_size_t vn)
{
	assert(un >= vn);
	assert(vn >= 1);

	/* We first multiply by the low order limb. This result can be
	   stored, not added, to rp. We also avoid a loop for zeroing this
	   way. */

	rp[un] = mpn_mul_1(rp, up, un, vp[0]);

	/* Now accumulate the product of up[] and the next higher limb from
	   vp[]. */

	while (--vn >= 1) {
		rp += 1, vp += 1;
		rp[un] = mpn_addmul_1(rp, up, un, vp[0]);
	}
	return rp[un];
}

void mpn_mul_n(mp_ptr rp, mp_srcptr ap, mp_srcptr bp, mp_size_t n)
{
	mpn_mul(rp, ap, n, bp, n);
}

mp_limb_t mpn_lshift(mp_ptr rp, mp_srcptr up, mp_size_t n, unsigned int cnt)
{
	mp_limb_t high_limb, low_limb;
	unsigned int tnc;
	mp_limb_t retval;

	assert(n >= 1);
	assert(cnt >= 1);
	assert(cnt < GMP_LIMB_BITS);

	up += n;
	rp += n;

	tnc = GMP_LIMB_BITS - cnt;
	low_limb = *--up;
	retval = low_limb >> tnc;
	high_limb = (low_limb << cnt);

	while (--n != 0) {
		low_limb = *--up;
		*--rp = high_limb | (low_limb >> tnc);
		high_limb = (low_limb << cnt);
	}
	*--rp = high_limb;

	return retval;
}

mp_limb_t mpn_rshift(mp_ptr rp, mp_srcptr up, mp_size_t n, unsigned int cnt)
{
	mp_limb_t high_limb, low_limb;
	unsigned int tnc;
	mp_limb_t retval;

	assert(n >= 1);
	assert(cnt >= 1);
	assert(cnt < GMP_LIMB_BITS);

	tnc = GMP_LIMB_BITS - cnt;
	high_limb = *up++;
	retval = (high_limb << tnc);
	low_limb = high_limb >> cnt;

	while (--n != 0) {
		high_limb = *up++;
		*rp++ = low_limb | (high_limb << tnc);
		low_limb = high_limb >> cnt;
	}
	*rp = low_limb;

	return retval;
}

static mp_bitcnt_t mpn_common_scan(mp_limb_t limb, mp_size_t i, mp_srcptr up, mp_size_t un, mp_limb_t ux)
{
	unsigned cnt;

	assert(ux == 0 || ux == GMP_LIMB_MAX);
	assert(0 <= i && i <= un);

	while (limb == 0) {
		i++;
		if (i == un)
			return (ux == 0 ? ~(mp_bitcnt_t) 0 : un * GMP_LIMB_BITS);
		limb = ux ^ up[i];
	}
	gmp_ctz(cnt, limb);
	return (mp_bitcnt_t) i *GMP_LIMB_BITS + cnt;
}

mp_bitcnt_t mpn_scan1(mp_srcptr ptr, mp_bitcnt_t bit)
{
	mp_size_t i;

	i = bit / GMP_LIMB_BITS;

	return mpn_common_scan(ptr[i] & (GMP_LIMB_MAX << (bit % GMP_LIMB_BITS)), i, ptr, i, 0);
}

mp_bitcnt_t mpn_scan0(mp_srcptr ptr, mp_bitcnt_t bit)
{
	mp_size_t i;

	i = bit / GMP_LIMB_BITS;

	return mpn_common_scan(~ptr[i] & (GMP_LIMB_MAX << (bit % GMP_LIMB_BITS)), i, ptr, i, GMP_LIMB_MAX);
}

/* MPN division interface. */
mp_limb_t mpn_invert_3by2(mp_limb_t u1, mp_limb_t u0)
{
	mp_limb_t r, p, m;
	unsigned ul, uh;
	unsigned ql, qh;

	/* First, do a 2/1 inverse. */
	/* The inverse m is defined as floor( (B^2 - 1 - u1)/u1 ), so that 0 <
	 * B^2 - (B + m) u1 <= u1 */
	assert(u1 >= GMP_LIMB_HIGHBIT);

	ul = u1 & GMP_LLIMB_MASK;
	uh = u1 >> (GMP_LIMB_BITS / 2);

	qh = ~u1 / uh;
	r = ((~u1 - (mp_limb_t) qh * uh) << (GMP_LIMB_BITS / 2)) | GMP_LLIMB_MASK;

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
	p = (r >> (GMP_LIMB_BITS / 2)) * qh + r;
	ql = (p >> (GMP_LIMB_BITS / 2)) + 1;

	/* By the 3/2 method, we don't need the high half limb. */
	r = (r << (GMP_LIMB_BITS / 2)) + GMP_LLIMB_MASK - ql * u1;

	if (r >= (p << (GMP_LIMB_BITS / 2))) {
		ql--;
		r += u1;
	}
	m = ((mp_limb_t) qh << (GMP_LIMB_BITS / 2)) + ql;
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
		gmp_umul_ppmm(th, tl, u0, m);
		r += th;
		if (r < th) {
			m--;
			m -= ((r > u1) | ((r == u1) & (tl > u0)));
		}
	}

	return m;
}

struct gmp_div_inverse {
	/* Normalization shift count. */
	unsigned shift;
	/* Normalized divisor (d0 unused for mpn_div_qr_1) */
	mp_limb_t d1, d0;
	/* Inverse, for 2/1 or 3/2. */
	mp_limb_t di;
};

static void mpn_div_qr_1_invert(struct gmp_div_inverse *inv, mp_limb_t d)
{
	unsigned shift;

	assert(d > 0);
	gmp_clz(shift, d);
	inv->shift = shift;
	inv->d1 = d << shift;
	inv->di = mpn_invert_limb(inv->d1);
}

static void mpn_div_qr_2_invert(struct gmp_div_inverse *inv, mp_limb_t d1, mp_limb_t d0)
{
	unsigned shift;

	assert(d1 > 0);
	gmp_clz(shift, d1);
	inv->shift = shift;
	if (shift > 0) {
		d1 = (d1 << shift) | (d0 >> (GMP_LIMB_BITS - shift));
		d0 <<= shift;
	}
	inv->d1 = d1;
	inv->d0 = d0;
	inv->di = mpn_invert_3by2(d1, d0);
}

static void mpn_div_qr_invert(struct gmp_div_inverse *inv, mp_srcptr dp, mp_size_t dn)
{
	assert(dn > 0);

	if (dn == 1)
		mpn_div_qr_1_invert(inv, dp[0]);
	else if (dn == 2)
		mpn_div_qr_2_invert(inv, dp[1], dp[0]);
	else {
		unsigned shift;
		mp_limb_t d1, d0;

		d1 = dp[dn - 1];
		d0 = dp[dn - 2];
		assert(d1 > 0);
		gmp_clz(shift, d1);
		inv->shift = shift;
		if (shift > 0) {
			d1 = (d1 << shift) | (d0 >> (GMP_LIMB_BITS - shift));
			d0 = (d0 << shift) | (dp[dn - 3] >> (GMP_LIMB_BITS - shift));
		}
		inv->d1 = d1;
		inv->d0 = d0;
		inv->di = mpn_invert_3by2(d1, d0);
	}
}

/* Not matching current public gmp interface, rather corresponding to
   the sbpi1_div_* functions. */
static mp_limb_t mpn_div_qr_1_preinv(mp_ptr qp, mp_srcptr np, mp_size_t nn, const struct gmp_div_inverse *inv)
{
	mp_limb_t d, di;
	mp_limb_t r;
	mp_ptr tp = NULL;

	if (inv->shift > 0) {
		tp = gmp_xalloc_limbs(nn);
		r = mpn_lshift(tp, np, nn, inv->shift);
		np = tp;
	} else
		r = 0;

	d = inv->d1;
	di = inv->di;
	while (--nn >= 0) {
		mp_limb_t q;

		gmp_udiv_qrnnd_preinv(q, r, r, np[nn], d, di);
		if (qp)
			qp[nn] = q;
	}
	if (inv->shift > 0)
		gmp_free(tp);

	return r >> inv->shift;
}

mp_limb_t mpn_div_qr_1(mp_ptr qp, mp_srcptr np, mp_size_t nn, mp_limb_t d)
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

				gmp_ctz(shift, d);
				mpn_rshift(qp, np, nn, shift);
			}
		}
		return r;
	} else {
		struct gmp_div_inverse inv;

		mpn_div_qr_1_invert(&inv, d);
		return mpn_div_qr_1_preinv(qp, np, nn, &inv);
	}
}

static void mpn_div_qr_2_preinv(mp_ptr qp, mp_ptr rp, mp_srcptr np, mp_size_t nn, const struct gmp_div_inverse *inv)
{
	unsigned shift;
	mp_size_t i;
	mp_limb_t d1, d0, di, r1, r0;
	mp_ptr tp;

	assert(nn >= 2);
	shift = inv->shift;
	d1 = inv->d1;
	d0 = inv->d0;
	di = inv->di;

	if (shift > 0) {
		tp = gmp_xalloc_limbs(nn);
		r1 = mpn_lshift(tp, np, nn, shift);
		np = tp;
	} else
		r1 = 0;

	r0 = np[nn - 1];

	i = nn - 2;
	do {
		mp_limb_t n0, q;

		n0 = np[i];
		gmp_udiv_qr_3by2(q, r1, r0, r1, r0, n0, d1, d0, di);

		if (qp)
			qp[i] = q;
	}
	while (--i >= 0);

	if (shift > 0) {
		assert((r0 << (GMP_LIMB_BITS - shift)) == 0);
		r0 = (r0 >> shift) | (r1 << (GMP_LIMB_BITS - shift));
		r1 >>= shift;

		gmp_free(tp);
	}

	rp[1] = r1;
	rp[0] = r0;
}

static void mpn_div_qr_pi1(mp_ptr qp, mp_ptr np, mp_size_t nn, mp_limb_t n1, mp_srcptr dp, mp_size_t dn, mp_limb_t dinv)
{
	mp_size_t i;

	mp_limb_t d1, d0;
	mp_limb_t cy, cy1;
	mp_limb_t q;

	assert(dn > 2);
	assert(nn >= dn);

	d1 = dp[dn - 1];
	d0 = dp[dn - 2];

	assert((d1 & GMP_LIMB_HIGHBIT) != 0);
	/* Iteration variable is the index of the q limb.
	 *
	 * We divide <n1, np[dn-1+i], np[dn-2+i], np[dn-3+i],..., np[i]>
	 * by            <d1,          d0,        dp[dn-3],  ..., dp[0] >
	 */

	i = nn - dn;
	do {
		mp_limb_t n0 = np[dn - 1 + i];

		if (n1 == d1 && n0 == d0) {
			q = GMP_LIMB_MAX;
			mpn_submul_1(np + i, dp, dn, q);
			n1 = np[dn - 1 + i];	/* update n1, last loop's value will now be invalid */
		} else {
			gmp_udiv_qr_3by2(q, n1, n0, n1, n0, np[dn - 2 + i], d1, d0, dinv);

			cy = mpn_submul_1(np + i, dp, dn - 2, q);

			cy1 = n0 < cy;
			n0 = n0 - cy;
			cy = n1 < cy1;
			n1 = n1 - cy1;
			np[dn - 2 + i] = n0;

			if (cy != 0) {
				n1 += d1 + mpn_add_n(np + i, np + i, dp, dn - 1);
				q--;
			}
		}

		if (qp)
			qp[i] = q;
	}
	while (--i >= 0);

	np[dn - 1] = n1;
}

static void mpn_div_qr_preinv(mp_ptr qp, mp_ptr np, mp_size_t nn, mp_srcptr dp, mp_size_t dn, const struct gmp_div_inverse *inv)
{
	assert(dn > 0);
	assert(nn >= dn);

	if (dn == 1)
		np[0] = mpn_div_qr_1_preinv(qp, np, nn, inv);
	else if (dn == 2)
		mpn_div_qr_2_preinv(qp, np, np, nn, inv);
	else {
		mp_limb_t nh;
		unsigned shift;

		assert(inv->d1 == dp[dn - 1]);
		assert(inv->d0 == dp[dn - 2]);
		assert((inv->d1 & GMP_LIMB_HIGHBIT) != 0);

		shift = inv->shift;
		if (shift > 0)
			nh = mpn_lshift(np, np, nn, shift);
		else
			nh = 0;

		mpn_div_qr_pi1(qp, np, nn, nh, dp, dn, inv->di);

		if (shift > 0)
			gmp_assert_nocarry(mpn_rshift(np, np, dn, shift));
	}
}

static void mpn_div_qr(mp_ptr qp, mp_ptr np, mp_size_t nn, mp_srcptr dp, mp_size_t dn)
{
	struct gmp_div_inverse inv;
	mp_ptr tp = NULL;

	assert(dn > 0);
	assert(nn >= dn);

	mpn_div_qr_invert(&inv, dp, dn);
	if (dn > 2 && inv.shift > 0) {
		tp = gmp_xalloc_limbs(dn);
		gmp_assert_nocarry(mpn_lshift(tp, dp, dn, inv.shift));
		dp = tp;
	}
	mpn_div_qr_preinv(qp, np, nn, dp, dn, &inv);
	if (tp)
		gmp_free(tp);
}

static mp_bitcnt_t mpn_limb_size_in_base_2(mp_limb_t u)
{
	unsigned shift;

	assert(u > 0);
	gmp_clz(shift, u);
	return GMP_LIMB_BITS - shift;
}

/* MPZ interface */
void mpz_init(mpz_t r)
{
	r->_mp_alloc = 1;
	r->_mp_size = 0;
	r->_mp_d = gmp_xalloc_limbs(1);
}

/* The utility of this function is a bit limited, since many functions
   assigns the result variable using mpz_swap. */
void mpz_init2(mpz_t r, mp_bitcnt_t bits)
{
	mp_size_t rn;

	bits -= (bits != 0);	/* Round down, except if 0 */
	rn = 1 + bits / GMP_LIMB_BITS;

	r->_mp_alloc = rn;
	r->_mp_size = 0;
	r->_mp_d = gmp_xalloc_limbs(rn);
}

void mpz_clear(mpz_t r)
{
	gmp_free(r->_mp_d);
}

static mp_ptr mpz_realloc(mpz_t r, mp_size_t size)
{
	size = GMP_MAX(size, 1);

	r->_mp_d = gmp_xrealloc_limbs(r->_mp_d, size);
	r->_mp_alloc = size;

	if (GMP_ABS(r->_mp_size) > size)
		r->_mp_size = 0;

	return r->_mp_d;
}

/* Realloc for an mpz_t WHAT if it has less than NEEDED limbs.  */
#define MPZ_REALLOC(z,n) ((n) > (z)->_mp_alloc			\
			  ? mpz_realloc(z,n)			\
			  : (z)->_mp_d)

/* MPZ assignment and basic conversions. */
void mpz_set_si(mpz_t r, signed long int x)
{
	if (x >= 0)
		mpz_set_ui(r, x);
	else {			/* (x < 0) */
		r->_mp_size = -1;
		r->_mp_d[0] = GMP_NEG_CAST(unsigned long int, x);
	}
}

void mpz_set_ui(mpz_t r, unsigned long int x)
{
	if (x > 0) {
		r->_mp_size = 1;
		r->_mp_d[0] = x;
	} else
		r->_mp_size = 0;
}

void mpz_set(mpz_t r, const mpz_t x)
{
	/* Allow the NOP r == x */
	if (r != x) {
		mp_size_t n;
		mp_ptr rp;

		n = GMP_ABS(x->_mp_size);
		rp = MPZ_REALLOC(r, n);

		mpn_copyi(rp, x->_mp_d, n);
		r->_mp_size = x->_mp_size;
	}
}

void mpz_init_set_si(mpz_t r, signed long int x)
{
	mpz_init(r);
	mpz_set_si(r, x);
}

void mpz_init_set_ui(mpz_t r, unsigned long int x)
{
	mpz_init(r);
	mpz_set_ui(r, x);
}

void mpz_init_set(mpz_t r, const mpz_t x)
{
	mpz_init(r);
	mpz_set(r, x);
}

int mpz_fits_slong_p(const mpz_t u)
{
	mp_size_t us = u->_mp_size;

	if (us == 1)
		return u->_mp_d[0] < GMP_LIMB_HIGHBIT;
	else if (us == -1)
		return u->_mp_d[0] <= GMP_LIMB_HIGHBIT;
	else
		return (us == 0);
}

int mpz_fits_ulong_p(const mpz_t u)
{
	mp_size_t us = u->_mp_size;

	return (us == (us > 0));
}

long int mpz_get_si(const mpz_t u)
{
	mp_size_t us = u->_mp_size;

	if (us > 0)
		return (long)(u->_mp_d[0] & ~GMP_LIMB_HIGHBIT);
	else if (us < 0)
		return (long)(-u->_mp_d[0] | GMP_LIMB_HIGHBIT);
	else
		return 0;
}

unsigned long int mpz_get_ui(const mpz_t u)
{
	return u->_mp_size == 0 ? 0 : u->_mp_d[0];
}

size_t mpz_size(const mpz_t u)
{
	return GMP_ABS(u->_mp_size);
}

mp_limb_t mpz_getlimbn(const mpz_t u, mp_size_t n)
{
	if (n >= 0 && n < GMP_ABS(u->_mp_size))
		return u->_mp_d[n];
	else
		return 0;
}

void mpz_realloc2(mpz_t x, mp_bitcnt_t n)
{
	mpz_realloc(x, 1 + (n - (n != 0)) / GMP_LIMB_BITS);
}

mp_srcptr mpz_limbs_read(mpz_srcptr x)
{
	return x->_mp_d;;
}

mp_ptr mpz_limbs_modify(mpz_t x, mp_size_t n)
{
	assert(n > 0);
	return MPZ_REALLOC(x, n);
}

mp_ptr mpz_limbs_write(mpz_t x, mp_size_t n)
{
	return mpz_limbs_modify(x, n);
}

void mpz_limbs_finish(mpz_t x, mp_size_t xs)
{
	mp_size_t xn;

	xn = mpn_normalized_size(x->_mp_d, GMP_ABS(xs));
	x->_mp_size = xs < 0 ? -xn : xn;
}

mpz_srcptr mpz_roinit_n(mpz_t x, mp_srcptr xp, mp_size_t xs)
{
	x->_mp_alloc = 0;
	x->_mp_d = (mp_ptr) xp;
	mpz_limbs_finish(x, xs);
	return x;
}

/* Conversions and comparison to double. */
void mpz_set_d(mpz_t r, double x)
{
	int sign;
	mp_ptr rp;
	mp_size_t rn, i;
	double B;
	double Bi;
	mp_limb_t f;

	/* x != x is true when x is a NaN, and x == x * 0.5 is true when x is
	   zero or infinity. */
	if (x != x || x == x * 0.5) {
		r->_mp_size = 0;
		return;
	}

	sign = x < 0.0;
	if (sign)
		x = -x;

	if (x < 1.0) {
		r->_mp_size = 0;
		return;
	}
	B = 2.0 * (double)GMP_LIMB_HIGHBIT;
	Bi = 1.0 / B;
	for (rn = 1; x >= B; rn++)
		x *= Bi;

	rp = MPZ_REALLOC(r, rn);

	f = (mp_limb_t) x;
	x -= f;
	assert(x < 1.0);
	i = rn - 1;
	rp[i] = f;
	while (--i >= 0) {
		x = B * x;
		f = (mp_limb_t) x;
		x -= f;
		assert(x < 1.0);
		rp[i] = f;
	}

	r->_mp_size = sign ? -rn : rn;
}

void mpz_init_set_d(mpz_t r, double x)
{
	mpz_init(r);
	mpz_set_d(r, x);
}

double mpz_get_d(const mpz_t u)
{
	mp_size_t un;
	double x;
	double B = 2.0 * (double)GMP_LIMB_HIGHBIT;

	un = GMP_ABS(u->_mp_size);

	if (un == 0)
		return 0.0;

	x = u->_mp_d[--un];
	while (un > 0)
		x = B * x + u->_mp_d[--un];

	if (u->_mp_size < 0)
		x = -x;

	return x;
}

int mpz_cmpabs_d(const mpz_t x, double d)
{
	mp_size_t xn;
	double B, Bi;
	mp_size_t i;

	xn = x->_mp_size;
	d = GMP_ABS(d);

	if (xn != 0) {
		xn = GMP_ABS(xn);

		B = 2.0 * (double)GMP_LIMB_HIGHBIT;
		Bi = 1.0 / B;

		/* Scale d so it can be compared with the top limb. */
		for (i = 1; i < xn; i++)
			d *= Bi;

		if (d >= B)
			return -1;

		/* Compare floor(d) to top limb, subtract and cancel when equal. */
		for (i = xn; i-- > 0;) {
			mp_limb_t f, xl;

			f = (mp_limb_t) d;
			xl = x->_mp_d[i];
			if (xl > f)
				return 1;
			else if (xl < f)
				return -1;
			d = B * (d - f);
		}
	}
	return -(d > 0.0);
}

int mpz_cmp_d(const mpz_t x, double d)
{
	if (x->_mp_size < 0) {
		if (d >= 0.0)
			return -1;
		else
			return -mpz_cmpabs_d(x, d);
	} else {
		if (d < 0.0)
			return 1;
		else
			return mpz_cmpabs_d(x, d);
	}
}

/* MPZ comparisons and the like. */
int mpz_sgn(const mpz_t u)
{
	mp_size_t usize = u->_mp_size;

	return (usize > 0) - (usize < 0);
}

int mpz_cmp_si(const mpz_t u, long v)
{
	mp_size_t usize = u->_mp_size;

	if (usize < -1)
		return -1;
	else if (v >= 0)
		return mpz_cmp_ui(u, v);
	else if (usize >= 0)
		return 1;
	else {			/* usize == -1 */
		mp_limb_t ul = u->_mp_d[0];
		if ((mp_limb_t) GMP_NEG_CAST(unsigned long int, v) < ul)
			 return -1;

		else
			return (mp_limb_t) GMP_NEG_CAST(unsigned long int, v) > ul;
	}
}

int mpz_cmp_ui(const mpz_t u, unsigned long v)
{
	mp_size_t usize = u->_mp_size;

	if (usize > 1)
		return 1;
	else if (usize < 0)
		return -1;
	else {
		mp_limb_t ul = (usize > 0) ? u->_mp_d[0] : 0;

		return (ul > v) - (ul < v);
	}
}

int mpz_cmp(const mpz_t a, const mpz_t b)
{
	mp_size_t asize = a->_mp_size;
	mp_size_t bsize = b->_mp_size;

	if (asize != bsize)
		return (asize < bsize) ? -1 : 1;
	else if (asize >= 0)
		return mpn_cmp(a->_mp_d, b->_mp_d, asize);
	else
		return mpn_cmp(b->_mp_d, a->_mp_d, -asize);
}

int mpz_cmpabs_ui(const mpz_t u, unsigned long v)
{
	mp_size_t un = GMP_ABS(u->_mp_size);
	mp_limb_t ul;

	if (un > 1)
		return 1;

	ul = (un == 1) ? u->_mp_d[0] : 0;

	return (ul > v) - (ul < v);
}

int mpz_cmpabs(const mpz_t u, const mpz_t v)
{
	return mpn_cmp4(u->_mp_d, GMP_ABS(u->_mp_size), v->_mp_d, GMP_ABS(v->_mp_size));
}

void mpz_abs(mpz_t r, const mpz_t u)
{
	mpz_set(r, u);
	r->_mp_size = GMP_ABS(r->_mp_size);
}

void mpz_neg(mpz_t r, const mpz_t u)
{
	mpz_set(r, u);
	r->_mp_size = -r->_mp_size;
}

void mpz_swap(mpz_t u, mpz_t v)
{
	MP_SIZE_T_SWAP(u->_mp_size, v->_mp_size);
	MP_SIZE_T_SWAP(u->_mp_alloc, v->_mp_alloc);
	MP_PTR_SWAP(u->_mp_d, v->_mp_d);
}

/* MPZ addition and subtraction */

/* Adds to the absolute value. Returns new size, but doesn't store it. */
static mp_size_t mpz_abs_add_ui(mpz_t r, const mpz_t a, unsigned long b)
{
	mp_size_t an;
	mp_ptr rp;
	mp_limb_t cy;

	an = GMP_ABS(a->_mp_size);
	if (an == 0) {
		r->_mp_d[0] = b;
		return b > 0;
	}

	rp = MPZ_REALLOC(r, an + 1);

	cy = mpn_add_1(rp, a->_mp_d, an, b);
	rp[an] = cy;
	an += cy;

	return an;
}

/* Subtract from the absolute value. Returns new size, (or -1 on underflow),
   but doesn't store it. */
static mp_size_t mpz_abs_sub_ui(mpz_t r, const mpz_t a, unsigned long b)
{
	mp_size_t an = GMP_ABS(a->_mp_size);
	mp_ptr rp = MPZ_REALLOC(r, an);

	if (an == 0) {
		rp[0] = b;
		return -(b > 0);
	} else if (an == 1 && a->_mp_d[0] < b) {
		rp[0] = b - a->_mp_d[0];
		return -1;
	} else {
		gmp_assert_nocarry(mpn_sub_1(rp, a->_mp_d, an, b));
		return mpn_normalized_size(rp, an);
	}
}

void mpz_add_ui(mpz_t r, const mpz_t a, unsigned long b)
{
	if (a->_mp_size >= 0)
		r->_mp_size = mpz_abs_add_ui(r, a, b);
	else
		r->_mp_size = -mpz_abs_sub_ui(r, a, b);
}

void mpz_sub_ui(mpz_t r, const mpz_t a, unsigned long b)
{
	if (a->_mp_size < 0)
		r->_mp_size = -mpz_abs_add_ui(r, a, b);
	else
		r->_mp_size = mpz_abs_sub_ui(r, a, b);
}

void mpz_ui_sub(mpz_t r, unsigned long a, const mpz_t b)
{
	if (b->_mp_size < 0)
		r->_mp_size = mpz_abs_add_ui(r, b, a);
	else
		r->_mp_size = -mpz_abs_sub_ui(r, b, a);
}

static mp_size_t mpz_abs_add(mpz_t r, const mpz_t a, const mpz_t b)
{
	mp_size_t an = GMP_ABS(a->_mp_size);
	mp_size_t bn = GMP_ABS(b->_mp_size);
	mp_ptr rp;
	mp_limb_t cy;

	if (an < bn) {
		MPZ_SRCPTR_SWAP(a, b);
		MP_SIZE_T_SWAP(an, bn);
	}

	rp = MPZ_REALLOC(r, an + 1);
	cy = mpn_add(rp, a->_mp_d, an, b->_mp_d, bn);

	rp[an] = cy;

	return an + cy;
}

static mp_size_t mpz_abs_sub(mpz_t r, const mpz_t a, const mpz_t b)
{
	mp_size_t an = GMP_ABS(a->_mp_size);
	mp_size_t bn = GMP_ABS(b->_mp_size);
	int cmp;
	mp_ptr rp;

	cmp = mpn_cmp4(a->_mp_d, an, b->_mp_d, bn);
	if (cmp > 0) {
		rp = MPZ_REALLOC(r, an);
		gmp_assert_nocarry(mpn_sub(rp, a->_mp_d, an, b->_mp_d, bn));
		return mpn_normalized_size(rp, an);
	} else if (cmp < 0) {
		rp = MPZ_REALLOC(r, bn);
		gmp_assert_nocarry(mpn_sub(rp, b->_mp_d, bn, a->_mp_d, an));
		return -mpn_normalized_size(rp, bn);
	} else
		return 0;
}

void mpz_add(mpz_t r, const mpz_t a, const mpz_t b)
{
	mp_size_t rn;

	if ((a->_mp_size ^ b->_mp_size) >= 0)
		rn = mpz_abs_add(r, a, b);
	else
		rn = mpz_abs_sub(r, a, b);

	r->_mp_size = a->_mp_size >= 0 ? rn : -rn;
}

void mpz_sub(mpz_t r, const mpz_t a, const mpz_t b)
{
	mp_size_t rn;

	if ((a->_mp_size ^ b->_mp_size) >= 0)
		rn = mpz_abs_sub(r, a, b);
	else
		rn = mpz_abs_add(r, a, b);

	r->_mp_size = a->_mp_size >= 0 ? rn : -rn;
}


/* MPZ multiplication */
void mpz_mul_si(mpz_t r, const mpz_t u, long int v)
{
	if (v < 0) {
		mpz_mul_ui(r, u, GMP_NEG_CAST(unsigned long int, v));

		mpz_neg(r, r);
	} else
		mpz_mul_ui(r, u, (unsigned long int)v);
}

void mpz_mul_ui(mpz_t r, const mpz_t u, unsigned long int v)
{
	mp_size_t un, us;
	mp_ptr tp;
	mp_limb_t cy;

	us = u->_mp_size;

	if (us == 0 || v == 0) {
		r->_mp_size = 0;
		return;
	}

	un = GMP_ABS(us);

	tp = MPZ_REALLOC(r, un + 1);
	cy = mpn_mul_1(tp, u->_mp_d, un, v);
	tp[un] = cy;

	un += (cy > 0);
	r->_mp_size = (us < 0) ? -un : un;
}

void mpz_mul(mpz_t r, const mpz_t u, const mpz_t v)
{
	int sign;
	mp_size_t un, vn, rn;
	mpz_t t;
	mp_ptr tp;

	un = u->_mp_size;
	vn = v->_mp_size;

	if (un == 0 || vn == 0) {
		r->_mp_size = 0;
		return;
	}

	sign = (un ^ vn) < 0;

	un = GMP_ABS(un);
	vn = GMP_ABS(vn);

	mpz_init2(t, (un + vn) * GMP_LIMB_BITS);

	tp = t->_mp_d;
	if (un >= vn)
		mpn_mul(tp, u->_mp_d, un, v->_mp_d, vn);
	else
		mpn_mul(tp, v->_mp_d, vn, u->_mp_d, un);

	rn = un + vn;
	rn -= tp[rn - 1] == 0;

	t->_mp_size = sign ? -rn : rn;
	mpz_swap(r, t);
	mpz_clear(t);
}

void mpz_mul_2exp(mpz_t r, const mpz_t u, mp_bitcnt_t bits)
{
	mp_size_t un, rn;
	mp_size_t limbs;
	unsigned shift;
	mp_ptr rp;

	un = GMP_ABS(u->_mp_size);
	if (un == 0) {
		r->_mp_size = 0;
		return;
	}

	limbs = bits / GMP_LIMB_BITS;
	shift = bits % GMP_LIMB_BITS;

	rn = un + limbs + (shift > 0);
	rp = MPZ_REALLOC(r, rn);
	if (shift > 0) {
		mp_limb_t cy = mpn_lshift(rp + limbs, u->_mp_d, un, shift);

		rp[rn - 1] = cy;
		rn -= (cy == 0);
	} else
		mpn_copyd(rp + limbs, u->_mp_d, un);

	mpn_zero(rp, limbs);

	r->_mp_size = (u->_mp_size < 0) ? -rn : rn;
}

void mpz_addmul_ui(mpz_t r, const mpz_t u, unsigned long int v)
{
	mpz_t t;

	mpz_init(t);
	mpz_mul_ui(t, u, v);
	mpz_add(r, r, t);
	mpz_clear(t);
}

void mpz_submul_ui(mpz_t r, const mpz_t u, unsigned long int v)
{
	mpz_t t;

	mpz_init(t);
	mpz_mul_ui(t, u, v);
	mpz_sub(r, r, t);
	mpz_clear(t);
}

void mpz_addmul(mpz_t r, const mpz_t u, const mpz_t v)
{
	mpz_t t;

	mpz_init(t);
	mpz_mul(t, u, v);
	mpz_add(r, r, t);
	mpz_clear(t);
}

void mpz_submul(mpz_t r, const mpz_t u, const mpz_t v)
{
	mpz_t t;

	mpz_init(t);
	mpz_mul(t, u, v);
	mpz_sub(r, r, t);
	mpz_clear(t);
}

/* Allows q or r to be zero. Returns 1 iff remainder is non-zero. */
static int mpz_div_qr(mpz_t q, mpz_t r, const mpz_t n, const mpz_t d, enum mpz_div_round_mode mode)
{
	mp_size_t ns, ds, nn, dn, qs;

	ns = n->_mp_size;
	ds = d->_mp_size;

	if (ds == 0)
		gmp_die("mpz_div_qr: Divide by zero.");

	if (ns == 0) {
		if (q)
			q->_mp_size = 0;
		if (r)
			r->_mp_size = 0;
		return 0;
	}

	nn = GMP_ABS(ns);
	dn = GMP_ABS(ds);

	qs = ds ^ ns;

	if (nn < dn) {
		if (mode == GMP_DIV_CEIL && qs >= 0) {
			/* q = 1, r = n - d */
			if (r)
				mpz_sub(r, n, d);
			if (q)
				mpz_set_ui(q, 1);
		} else if (mode == GMP_DIV_FLOOR && qs < 0) {
			/* q = -1, r = n + d */
			if (r)
				mpz_add(r, n, d);
			if (q)
				mpz_set_si(q, -1);
		} else {
			/* q = 0, r = d */
			if (r)
				mpz_set(r, n);
			if (q)
				q->_mp_size = 0;
		}
		return 1;
	} else {
		mp_ptr np, qp;
		mp_size_t qn, rn;
		mpz_t tq, tr;

		mpz_init_set(tr, n);
		np = tr->_mp_d;

		qn = nn - dn + 1;

		if (q) {
			mpz_init2(tq, qn * GMP_LIMB_BITS);
			qp = tq->_mp_d;
		} else
			qp = NULL;

		mpn_div_qr(qp, np, nn, d->_mp_d, dn);

		if (qp) {
			qn -= (qp[qn - 1] == 0);

			tq->_mp_size = qs < 0 ? -qn : qn;
		}
		rn = mpn_normalized_size(np, dn);
		tr->_mp_size = ns < 0 ? -rn : rn;

		if (mode == GMP_DIV_FLOOR && qs < 0 && rn != 0) {
			if (q)
				mpz_sub_ui(tq, tq, 1);
			if (r)
				mpz_add(tr, tr, d);
		} else if (mode == GMP_DIV_CEIL && qs >= 0 && rn != 0) {
			if (q)
				mpz_add_ui(tq, tq, 1);
			if (r)
				mpz_sub(tr, tr, d);
		}

		if (q) {
			mpz_swap(tq, q);
			mpz_clear(tq);
		}
		if (r)
			mpz_swap(tr, r);

		mpz_clear(tr);

		return rn != 0;
	}
}

void mpz_cdiv_qr(mpz_t q, mpz_t r, const mpz_t n, const mpz_t d)
{
	mpz_div_qr(q, r, n, d, GMP_DIV_CEIL);
}

void mpz_fdiv_qr(mpz_t q, mpz_t r, const mpz_t n, const mpz_t d)
{
	mpz_div_qr(q, r, n, d, GMP_DIV_FLOOR);
}

void mpz_tdiv_qr(mpz_t q, mpz_t r, const mpz_t n, const mpz_t d)
{
	mpz_div_qr(q, r, n, d, GMP_DIV_TRUNC);
}

void mpz_cdiv_q(mpz_t q, const mpz_t n, const mpz_t d)
{
	mpz_div_qr(q, NULL, n, d, GMP_DIV_CEIL);
}

void mpz_fdiv_q(mpz_t q, const mpz_t n, const mpz_t d)
{
	mpz_div_qr(q, NULL, n, d, GMP_DIV_FLOOR);
}

void mpz_tdiv_q(mpz_t q, const mpz_t n, const mpz_t d)
{
	mpz_div_qr(q, NULL, n, d, GMP_DIV_TRUNC);
}

void mpz_cdiv_r(mpz_t r, const mpz_t n, const mpz_t d)
{
	mpz_div_qr(NULL, r, n, d, GMP_DIV_CEIL);
}

void mpz_fdiv_r(mpz_t r, const mpz_t n, const mpz_t d)
{
	mpz_div_qr(NULL, r, n, d, GMP_DIV_FLOOR);
}

void mpz_tdiv_r(mpz_t r, const mpz_t n, const mpz_t d)
{
	mpz_div_qr(NULL, r, n, d, GMP_DIV_TRUNC);
}

void mpz_mod(mpz_t r, const mpz_t n, const mpz_t d)
{
	mpz_div_qr(NULL, r, n, d, d->_mp_size >= 0 ? GMP_DIV_FLOOR : GMP_DIV_CEIL);
}

static void mpz_div_q_2exp(mpz_t q, const mpz_t u, mp_bitcnt_t bit_index, enum mpz_div_round_mode mode)
{
	mp_size_t un, qn;
	mp_size_t limb_cnt;
	mp_ptr qp;
	int adjust;

	un = u->_mp_size;
	if (un == 0) {
		q->_mp_size = 0;
		return;
	}
	limb_cnt = bit_index / GMP_LIMB_BITS;
	qn = GMP_ABS(un) - limb_cnt;
	bit_index %= GMP_LIMB_BITS;

	if (mode == ((un > 0) ? GMP_DIV_CEIL : GMP_DIV_FLOOR))	/* un != 0 here. */
		/* Note: Below, the final indexing at limb_cnt is valid because at
		   that point we have qn > 0. */
		adjust = (qn <= 0 || !mpn_zero_p(u->_mp_d, limb_cnt)
			  || (u->_mp_d[limb_cnt]
			      & (((mp_limb_t) 1 << bit_index) - 1)));
	else
		adjust = 0;

	if (qn <= 0)
		qn = 0;

	else {
		qp = MPZ_REALLOC(q, qn);

		if (bit_index != 0) {
			mpn_rshift(qp, u->_mp_d + limb_cnt, qn, bit_index);
			qn -= qp[qn - 1] == 0;
		} else {
			mpn_copyi(qp, u->_mp_d + limb_cnt, qn);
		}
	}

	q->_mp_size = qn;

	if (adjust)
		mpz_add_ui(q, q, 1);
	if (un < 0)
		mpz_neg(q, q);
}

static void mpz_div_r_2exp(mpz_t r, const mpz_t u, mp_bitcnt_t bit_index, enum mpz_div_round_mode mode)
{
	mp_size_t us, un, rn;
	mp_ptr rp;
	mp_limb_t mask;

	us = u->_mp_size;
	if (us == 0 || bit_index == 0) {
		r->_mp_size = 0;
		return;
	}
	rn = (bit_index + GMP_LIMB_BITS - 1) / GMP_LIMB_BITS;
	assert(rn > 0);

	rp = MPZ_REALLOC(r, rn);
	un = GMP_ABS(us);

	mask = GMP_LIMB_MAX >> (rn * GMP_LIMB_BITS - bit_index);

	if (rn > un) {
		/* Quotient (with truncation) is zero, and remainder is
		   non-zero */
		if (mode == ((us > 0) ? GMP_DIV_CEIL : GMP_DIV_FLOOR)) {	/* us != 0 here. */
			/* Have to negate and sign extend. */
			mp_size_t i;
			mp_limb_t cy;

			for (cy = 1, i = 0; i < un; i++) {
				mp_limb_t s = ~u->_mp_d[i] + cy;

				cy = s < cy;
				rp[i] = s;
			}
			assert(cy == 0);
			for (; i < rn - 1; i++)
				rp[i] = GMP_LIMB_MAX;

			rp[rn - 1] = mask;
			us = -us;
		} else {
			/* Just copy */
			if (r != u)
				mpn_copyi(rp, u->_mp_d, un);

			rn = un;
		}
	} else {
		if (r != u)
			mpn_copyi(rp, u->_mp_d, rn - 1);

		rp[rn - 1] = u->_mp_d[rn - 1] & mask;

		if (mode == ((us > 0) ? GMP_DIV_CEIL : GMP_DIV_FLOOR)) {	/* us != 0 here. */
			/* If r != 0, compute 2^{bit_count} - r. */
			mp_size_t i;

			for (i = 0; i < rn && rp[i] == 0; i++) ;
			if (i < rn) {
				/* r > 0, need to flip sign. */
				rp[i] = ~rp[i] + 1;
				while (++i < rn)
					rp[i] = ~rp[i];

				rp[rn - 1] &= mask;

				/* us is not used for anything else, so we can modify it
				   here to indicate flipped sign. */
				us = -us;
			}
		}
	}
	rn = mpn_normalized_size(rp, rn);
	r->_mp_size = us < 0 ? -rn : rn;
}

void mpz_cdiv_q_2exp(mpz_t r, const mpz_t u, mp_bitcnt_t cnt)
{
	mpz_div_q_2exp(r, u, cnt, GMP_DIV_CEIL);
}

void mpz_fdiv_q_2exp(mpz_t r, const mpz_t u, mp_bitcnt_t cnt)
{
	mpz_div_q_2exp(r, u, cnt, GMP_DIV_FLOOR);
}

void mpz_tdiv_q_2exp(mpz_t r, const mpz_t u, mp_bitcnt_t cnt)
{
	mpz_div_q_2exp(r, u, cnt, GMP_DIV_TRUNC);
}

void mpz_cdiv_r_2exp(mpz_t r, const mpz_t u, mp_bitcnt_t cnt)
{
	mpz_div_r_2exp(r, u, cnt, GMP_DIV_CEIL);
}

void mpz_fdiv_r_2exp(mpz_t r, const mpz_t u, mp_bitcnt_t cnt)
{
	mpz_div_r_2exp(r, u, cnt, GMP_DIV_FLOOR);
}

void mpz_tdiv_r_2exp(mpz_t r, const mpz_t u, mp_bitcnt_t cnt)
{
	mpz_div_r_2exp(r, u, cnt, GMP_DIV_TRUNC);
}

void mpz_divexact(mpz_t q, const mpz_t n, const mpz_t d)
{
	gmp_assert_nocarry(mpz_div_qr(q, NULL, n, d, GMP_DIV_TRUNC));
}

int mpz_divisible_p(const mpz_t n, const mpz_t d)
{
	return mpz_div_qr(NULL, NULL, n, d, GMP_DIV_TRUNC) == 0;
}

int mpz_congruent_p(const mpz_t a, const mpz_t b, const mpz_t m)
{
	mpz_t t;
	int res;

	/* a == b (mod 0) iff a == b */
	if (mpz_sgn(m) == 0)
		return (mpz_cmp(a, b) == 0);

	mpz_init(t);
	mpz_sub(t, a, b);
	res = mpz_divisible_p(t, m);
	mpz_clear(t);

	return res;
}

unsigned long mpz_div_qr_ui(mpz_t q, mpz_t r, const mpz_t n, unsigned long d, enum mpz_div_round_mode mode)
{
	mp_size_t ns, qn;
	mp_ptr qp;
	mp_limb_t rl;
	mp_size_t rs;

	ns = n->_mp_size;
	if (ns == 0) {
		if (q)
			q->_mp_size = 0;
		if (r)
			r->_mp_size = 0;
		return 0;
	}

	qn = GMP_ABS(ns);
	if (q)
		qp = MPZ_REALLOC(q, qn);
	else
		qp = NULL;

	rl = mpn_div_qr_1(qp, n->_mp_d, qn, d);
	assert(rl < d);

	rs = rl > 0;
	rs = (ns < 0) ? -rs : rs;

	if (rl > 0 && ((mode == GMP_DIV_FLOOR && ns < 0)
		       || (mode == GMP_DIV_CEIL && ns >= 0))) {
		if (q)
			gmp_assert_nocarry(mpn_add_1(qp, qp, qn, 1));
		rl = d - rl;
		rs = -rs;
	}

	if (r) {
		r->_mp_d[0] = rl;
		r->_mp_size = rs;
	}
	if (q) {
		qn -= (qp[qn - 1] == 0);
		assert(qn == 0 || qp[qn - 1] > 0);

		q->_mp_size = (ns < 0) ? -qn : qn;
	}

	return rl;
}

unsigned long mpz_cdiv_qr_ui(mpz_t q, mpz_t r, const mpz_t n, unsigned long d)
{
	return mpz_div_qr_ui(q, r, n, d, GMP_DIV_CEIL);
}

unsigned long mpz_fdiv_qr_ui(mpz_t q, mpz_t r, const mpz_t n, unsigned long d)
{
	return mpz_div_qr_ui(q, r, n, d, GMP_DIV_FLOOR);
}

unsigned long mpz_tdiv_qr_ui(mpz_t q, mpz_t r, const mpz_t n, unsigned long d)
{
	return mpz_div_qr_ui(q, r, n, d, GMP_DIV_TRUNC);
}

unsigned long mpz_cdiv_q_ui(mpz_t q, const mpz_t n, unsigned long d)
{
	return mpz_div_qr_ui(q, NULL, n, d, GMP_DIV_CEIL);
}

unsigned long mpz_fdiv_q_ui(mpz_t q, const mpz_t n, unsigned long d)
{
	return mpz_div_qr_ui(q, NULL, n, d, GMP_DIV_FLOOR);
}

unsigned long mpz_tdiv_q_ui(mpz_t q, const mpz_t n, unsigned long d)
{
	return mpz_div_qr_ui(q, NULL, n, d, GMP_DIV_TRUNC);
}

unsigned long mpz_cdiv_r_ui(mpz_t r, const mpz_t n, unsigned long d)
{
	return mpz_div_qr_ui(NULL, r, n, d, GMP_DIV_CEIL);
}

unsigned long mpz_fdiv_r_ui(mpz_t r, const mpz_t n, unsigned long d)
{
	return mpz_div_qr_ui(NULL, r, n, d, GMP_DIV_FLOOR);
}

unsigned long mpz_tdiv_r_ui(mpz_t r, const mpz_t n, unsigned long d)
{
	return mpz_div_qr_ui(NULL, r, n, d, GMP_DIV_TRUNC);
}

unsigned long mpz_cdiv_ui(const mpz_t n, unsigned long d)
{
	return mpz_div_qr_ui(NULL, NULL, n, d, GMP_DIV_CEIL);
}

unsigned long mpz_fdiv_ui(const mpz_t n, unsigned long d)
{
	return mpz_div_qr_ui(NULL, NULL, n, d, GMP_DIV_FLOOR);
}

unsigned long mpz_tdiv_ui(const mpz_t n, unsigned long d)
{
	return mpz_div_qr_ui(NULL, NULL, n, d, GMP_DIV_TRUNC);
}

unsigned long mpz_mod_ui(mpz_t r, const mpz_t n, unsigned long d)
{
	return mpz_div_qr_ui(NULL, r, n, d, GMP_DIV_FLOOR);
}

void mpz_divexact_ui(mpz_t q, const mpz_t n, unsigned long d)
{
	gmp_assert_nocarry(mpz_div_qr_ui(q, NULL, n, d, GMP_DIV_TRUNC));
}

int mpz_divisible_ui_p(const mpz_t n, unsigned long d)
{
	return mpz_div_qr_ui(NULL, NULL, n, d, GMP_DIV_TRUNC) == 0;
}


/* GCD */
static mp_limb_t mpn_gcd_11(mp_limb_t u, mp_limb_t v)
{
	unsigned shift;

	assert((u | v) > 0);

	if (u == 0)
		return v;
	else if (v == 0)
		return u;

	gmp_ctz(shift, u | v);

	u >>= shift;
	v >>= shift;

	if ((u & 1) == 0)
		MP_LIMB_T_SWAP(u, v);

	while ((v & 1) == 0)
		v >>= 1;

	while (u != v) {
		if (u > v) {
			u -= v;
			do
				u >>= 1;
			while ((u & 1) == 0);
		} else {
			v -= u;
			do
				v >>= 1;
			while ((v & 1) == 0);
		}
	}
	return u << shift;
}

unsigned long mpz_gcd_ui(mpz_t g, const mpz_t u, unsigned long v)
{
	mp_size_t un;

	if (v == 0) {
		if (g)
			mpz_abs(g, u);
	} else {
		un = GMP_ABS(u->_mp_size);
		if (un != 0)
			v = mpn_gcd_11(mpn_div_qr_1(NULL, u->_mp_d, un, v), v);

		if (g)
			mpz_set_ui(g, v);
	}

	return v;
}

static mp_bitcnt_t mpz_make_odd(mpz_t r)
{
	mp_bitcnt_t shift;

	assert(r->_mp_size > 0);
	/* Count trailing zeros, equivalent to mpn_scan1, because we know that there is a 1 */
	shift = mpn_common_scan(r->_mp_d[0], 0, r->_mp_d, 0, 0);
	mpz_tdiv_q_2exp(r, r, shift);

	return shift;
}

void mpz_gcdext(mpz_t g, mpz_t s, mpz_t t, const mpz_t u, const mpz_t v)
{
	mpz_t tu, tv, s0, s1, t0, t1;
	mp_bitcnt_t uz, vz, gz;
	mp_bitcnt_t power;

	if (u->_mp_size == 0) {
		/* g = 0 u + sgn(v) v */
		signed long sign = mpz_sgn(v);

		mpz_abs(g, v);
		if (s)
			mpz_set_ui(s, 0);
		if (t)
			mpz_set_si(t, sign);
		return;
	}

	if (v->_mp_size == 0) {
		/* g = sgn(u) u + 0 v */
		signed long sign = mpz_sgn(u);

		mpz_abs(g, u);
		if (s)
			mpz_set_si(s, sign);
		if (t)
			mpz_set_ui(t, 0);
		return;
	}

	mpz_init(tu);
	mpz_init(tv);
	mpz_init(s0);
	mpz_init(s1);
	mpz_init(t0);
	mpz_init(t1);

	mpz_abs(tu, u);
	uz = mpz_make_odd(tu);
	mpz_abs(tv, v);
	vz = mpz_make_odd(tv);
	gz = GMP_MIN(uz, vz);

	uz -= gz;
	vz -= gz;

	/* Cofactors corresponding to odd gcd. gz handled later. */
	if (tu->_mp_size < tv->_mp_size) {
		mpz_swap(tu, tv);
		MPZ_SRCPTR_SWAP(u, v);
		MPZ_PTR_SWAP(s, t);
		MP_BITCNT_T_SWAP(uz, vz);
	}

	/* Maintain
	 *
	 * u = t0 tu + t1 tv
	 * v = s0 tu + s1 tv
	 *
	 * where u and v denote the inputs with common factors of two
	 * eliminated, and det (s0, t0; s1, t1) = 2^p. Then
	 *
	 * 2^p tu =  s1 u - t1 v
	 * 2^p tv = -s0 u + t0 v
	 */

	/* After initial division, tu = q tv + tu', we have
	 *
	 * u = 2^uz (tu' + q tv)
	 * v = 2^vz tv
	 *
	 * or
	 *
	 * t0 = 2^uz, t1 = 2^uz q
	 * s0 = 0,    s1 = 2^vz
	 */

	mpz_setbit(t0, uz);
	mpz_tdiv_qr(t1, tu, tu, tv);
	mpz_mul_2exp(t1, t1, uz);

	mpz_setbit(s1, vz);
	power = uz + vz;

	if (tu->_mp_size > 0) {
		mp_bitcnt_t shift;

		shift = mpz_make_odd(tu);
		mpz_mul_2exp(t0, t0, shift);
		mpz_mul_2exp(s0, s0, shift);
		power += shift;

		for (;;) {
			int c;

			c = mpz_cmp(tu, tv);
			if (c == 0)
				break;

			if (c < 0) {
				/* tv = tv' + tu
				 *
				 * u = t0 tu + t1 (tv' + tu) = (t0 + t1) tu + t1 tv'
				 * v = s0 tu + s1 (tv' + tu) = (s0 + s1) tu + s1 tv' */

				mpz_sub(tv, tv, tu);
				mpz_add(t0, t0, t1);
				mpz_add(s0, s0, s1);

				shift = mpz_make_odd(tv);
				mpz_mul_2exp(t1, t1, shift);
				mpz_mul_2exp(s1, s1, shift);
			} else {
				mpz_sub(tu, tu, tv);
				mpz_add(t1, t0, t1);
				mpz_add(s1, s0, s1);

				shift = mpz_make_odd(tu);
				mpz_mul_2exp(t0, t0, shift);
				mpz_mul_2exp(s0, s0, shift);
			}
			power += shift;
		}
	}

	/* Now tv = odd part of gcd, and -s0 and t0 are corresponding
	   cofactors. */

	mpz_mul_2exp(tv, tv, gz);
	mpz_neg(s0, s0);

	/* 2^p g = s0 u + t0 v. Eliminate one factor of two at a time. To
	   adjust cofactors, we need u / g and v / g */

	mpz_divexact(s1, v, tv);
	mpz_abs(s1, s1);
	mpz_divexact(t1, u, tv);
	mpz_abs(t1, t1);

	while (power-- > 0) {
		/* s0 u + t0 v = (s0 - v/g) u - (t0 + u/g) v */
		if (mpz_odd_p(s0) || mpz_odd_p(t0)) {
			mpz_sub(s0, s0, s1);
			mpz_add(t0, t0, t1);
		}
		mpz_divexact_ui(s0, s0, 2);
		mpz_divexact_ui(t0, t0, 2);
	}

	/* Arrange so that |s| < |u| / 2g */
	mpz_add(s1, s0, s1);
	if (mpz_cmpabs(s0, s1) > 0) {
		mpz_swap(s0, s1);
		mpz_sub(t0, t0, t1);
	}
	if (u->_mp_size < 0)
		mpz_neg(s0, s0);
	if (v->_mp_size < 0)
		mpz_neg(t0, t0);

	mpz_swap(g, tv);
	if (s)
		mpz_swap(s, s0);
	if (t)
		mpz_swap(t, t0);

	mpz_clear(tu);
	mpz_clear(tv);
	mpz_clear(s0);
	mpz_clear(s1);
	mpz_clear(t0);
	mpz_clear(t1);
}

int mpz_invert(mpz_t r, const mpz_t u, const mpz_t m)
{
	mpz_t g, tr;
	int invertible;

	if (u->_mp_size == 0 || mpz_cmpabs_ui(m, 1) <= 0)
		return 0;

	mpz_init(g);
	mpz_init(tr);

	mpz_gcdext(g, tr, NULL, u, m);
	invertible = (mpz_cmp_ui(g, 1) == 0);

	if (invertible) {
		if (tr->_mp_size < 0) {
			if (m->_mp_size >= 0)
				mpz_add(tr, tr, m);
			else
				mpz_sub(tr, tr, m);
		}
		mpz_swap(r, tr);
	}

	mpz_clear(g);
	mpz_clear(tr);
	return invertible;
}

void mpz_powm(mpz_t r, const mpz_t b, const mpz_t e, const mpz_t m)
{
	mpz_t tr;
	mpz_t base;
	mp_size_t en, mn;
	mp_srcptr mp;
	struct gmp_div_inverse minv;
	unsigned shift;
	mp_ptr tp = NULL;

	en = GMP_ABS(e->_mp_size);
	mn = GMP_ABS(m->_mp_size);
	if (mn == 0)
		gmp_die("mpz_powm: Zero modulo.");

	if (en == 0) {
		mpz_set_ui(r, 1);
		return;
	}

	mp = m->_mp_d;
	mpn_div_qr_invert(&minv, mp, mn);
	shift = minv.shift;

	if (shift > 0) {
		/* To avoid shifts, we do all our reductions, except the final
		   one, using a *normalized* m. */
		minv.shift = 0;

		tp = gmp_xalloc_limbs(mn);
		gmp_assert_nocarry(mpn_lshift(tp, mp, mn, shift));
		mp = tp;
	}

	mpz_init(base);

	if (e->_mp_size < 0) {
		if (!mpz_invert(base, b, m))
			gmp_die("mpz_powm: Negative exponent and non-invertible base.");
	} else {
		mp_size_t bn;

		mpz_abs(base, b);

		bn = base->_mp_size;
		if (bn >= mn) {
			mpn_div_qr_preinv(NULL, base->_mp_d, base->_mp_size, mp, mn, &minv);
			bn = mn;
		}

		/* We have reduced the absolute value. Now take care of the
		   sign. Note that we get zero represented non-canonically as
		   m. */
		if (b->_mp_size < 0) {
			mp_ptr bp = MPZ_REALLOC(base, mn);

			gmp_assert_nocarry(mpn_sub(bp, mp, mn, bp, bn));
			bn = mn;
		}
		base->_mp_size = mpn_normalized_size(base->_mp_d, bn);
	}
	mpz_init_set_ui(tr, 1);

	while (--en >= 0) {
		mp_limb_t w = e->_mp_d[en];
		mp_limb_t bit;

		bit = GMP_LIMB_HIGHBIT;
		do {
			mpz_mul(tr, tr, tr);
			if (w & bit)
				mpz_mul(tr, tr, base);
			if (tr->_mp_size > mn) {
				mpn_div_qr_preinv(NULL, tr->_mp_d, tr->_mp_size, mp, mn, &minv);
				tr->_mp_size = mpn_normalized_size(tr->_mp_d, mn);
			}
			bit >>= 1;
		}
		while (bit > 0);
	}

	/* Final reduction */
	if (tr->_mp_size >= mn) {
		minv.shift = shift;
		mpn_div_qr_preinv(NULL, tr->_mp_d, tr->_mp_size, mp, mn, &minv);
		tr->_mp_size = mpn_normalized_size(tr->_mp_d, mn);
	}
	if (tp)
		gmp_free(tp);

	mpz_swap(r, tr);
	mpz_clear(tr);
	mpz_clear(base);
}

void mpz_powm_ui(mpz_t r, const mpz_t b, unsigned long elimb, const mpz_t m)
{
	mpz_t e;

	mpz_powm(r, b, mpz_roinit_n(e, &elimb, 1), m);
}

/* Primality testing */
static int gmp_millerrabin(const mpz_t n, const mpz_t nm1, mpz_t y, const mpz_t q, mp_bitcnt_t k)
{
	assert(k > 0);

	/* Caller must initialize y to the base. */
	mpz_powm(y, y, q, n);

	if (mpz_cmp_ui(y, 1) == 0 || mpz_cmp(y, nm1) == 0)
		return 1;

	while (--k > 0) {
		mpz_powm_ui(y, y, 2, n);
		if (mpz_cmp(y, nm1) == 0)
			return 1;
		/* y == 1 means that the previous y was a non-trivial square root
		   of 1 (mod n). y == 0 means that n is a power of the base.
		   In either case, n is not prime. */
		if (mpz_cmp_ui(y, 1) <= 0)
			return 0;
	}
	return 0;
}

/* This product is 0xc0cfd797, and fits in 32 bits. */
#define GMP_PRIME_PRODUCT \
  (3UL*5UL*7UL*11UL*13UL*17UL*19UL*23UL*29UL)

/* Bit (p+1)/2 is set, for each odd prime <= 61 */
#define GMP_PRIME_MASK 0xc96996dcUL

int mpz_probab_prime_p(const mpz_t n, int reps)
{
	mpz_t nm1;
	mpz_t q;
	mpz_t y;
	mp_bitcnt_t k;
	int is_prime;
	int j;

	/* Note that we use the absolute value of n only, for compatibility
	   with the real GMP. */
	if (mpz_even_p(n))
		return (mpz_cmpabs_ui(n, 2) == 0) ? 2 : 0;

	/* Above test excludes n == 0 */
	assert(n->_mp_size != 0);

	if (mpz_cmpabs_ui(n, 64) < 0)
		return (GMP_PRIME_MASK >> (n->_mp_d[0] >> 1)) & 2;

	if (mpz_gcd_ui(NULL, n, GMP_PRIME_PRODUCT) != 1)
		return 0;

	/* All prime factors are >= 31. */
	if (mpz_cmpabs_ui(n, 31 * 31) < 0)
		return 2;

	/* Use Miller-Rabin, with a deterministic sequence of bases, a[j] =
	   j^2 + j + 41 using Euler's polynomial. We potentially stop early,
	   if a[j] >= n - 1. Since n >= 31*31, this can happen only if reps >
	   30 (a[30] == 971 > 31*31 == 961). */

	mpz_init(nm1);
	mpz_init(q);
	mpz_init(y);

	/* Find q and k, where q is odd and n = 1 + 2**k * q.  */
	nm1->_mp_size = mpz_abs_sub_ui(nm1, n, 1);
	k = mpz_scan1(nm1, 0);
	mpz_tdiv_q_2exp(q, nm1, k);

	for (j = 0, is_prime = 1; is_prime & (j < reps); j++) {
		mpz_set_ui(y, (unsigned long)j * j + j + 41);
		if (mpz_cmp(y, nm1) >= 0) {
			/* Don't try any further bases. This "early" break does not affect
			   the result for any reasonable reps value (<=5000 was tested) */
			assert(j >= 30);
			break;
		}
		is_prime = gmp_millerrabin(n, nm1, y, q, k);
	}
	mpz_clear(nm1);
	mpz_clear(q);
	mpz_clear(y);

	return is_prime;
}


/* Logical operations and bit manipulation. */

/* Numbers are treated as if represented in two's complement (and
   infinitely sign extended). For a negative values we get the two's
   complement from -x = ~x + 1, where ~ is bitwise complement.
   Negation transforms

     xxxx10...0

   into

     yyyy10...0

   where yyyy is the bitwise complement of xxxx. So least significant
   bits, up to and including the first one bit, are unchanged, and
   the more significant bits are all complemented.

   To change a bit from zero to one in a negative number, subtract the
   corresponding power of two from the absolute value. This can never
   underflow. To change a bit from one to zero, add the corresponding
   power of two, and this might overflow. E.g., if x = -001111, the
   two's complement is 110001. Clearing the least significant bit, we
   get two's complement 110000, and -010000. */
int mpz_tstbit(const mpz_t d, mp_bitcnt_t bit_index)
{
	mp_size_t limb_index;
	unsigned shift;
	mp_size_t ds;
	mp_size_t dn;
	mp_limb_t w;
	int bit;

	ds = d->_mp_size;
	dn = GMP_ABS(ds);
	limb_index = bit_index / GMP_LIMB_BITS;
	if (limb_index >= dn)
		return ds < 0;

	shift = bit_index % GMP_LIMB_BITS;
	w = d->_mp_d[limb_index];
	bit = (w >> shift) & 1;

	if (ds < 0) {
		/* d < 0. Check if any of the bits below is set: If so, our bit
		   must be complemented. */
		if (shift > 0 && (w << (GMP_LIMB_BITS - shift)) > 0)
			return bit ^ 1;
		while (--limb_index >= 0)
			if (d->_mp_d[limb_index] > 0)
				return bit ^ 1;
	}
	return bit;
}

static void mpz_abs_add_bit(mpz_t d, mp_bitcnt_t bit_index)
{
	mp_size_t dn, limb_index;
	mp_limb_t bit;
	mp_ptr dp;

	dn = GMP_ABS(d->_mp_size);

	limb_index = bit_index / GMP_LIMB_BITS;
	bit = (mp_limb_t) 1 << (bit_index % GMP_LIMB_BITS);

	if (limb_index >= dn) {
		mp_size_t i;

		/* The bit should be set outside of the end of the number.
		   We have to increase the size of the number. */
		dp = MPZ_REALLOC(d, limb_index + 1);

		dp[limb_index] = bit;
		for (i = dn; i < limb_index; i++)
			dp[i] = 0;
		dn = limb_index + 1;
	} else {
		mp_limb_t cy;

		dp = d->_mp_d;

		cy = mpn_add_1(dp + limb_index, dp + limb_index, dn - limb_index, bit);
		if (cy > 0) {
			dp = MPZ_REALLOC(d, dn + 1);
			dp[dn++] = cy;
		}
	}

	d->_mp_size = (d->_mp_size < 0) ? -dn : dn;
}

static void mpz_abs_sub_bit(mpz_t d, mp_bitcnt_t bit_index)
{
	mp_size_t dn, limb_index;
	mp_ptr dp;
	mp_limb_t bit;

	dn = GMP_ABS(d->_mp_size);
	dp = d->_mp_d;

	limb_index = bit_index / GMP_LIMB_BITS;
	bit = (mp_limb_t) 1 << (bit_index % GMP_LIMB_BITS);

	assert(limb_index < dn);

	gmp_assert_nocarry(mpn_sub_1(dp + limb_index, dp + limb_index, dn - limb_index, bit));
	dn = mpn_normalized_size(dp, dn);
	d->_mp_size = (d->_mp_size < 0) ? -dn : dn;
}

void mpz_setbit(mpz_t d, mp_bitcnt_t bit_index)
{
	if (!mpz_tstbit(d, bit_index)) {
		if (d->_mp_size >= 0)
			mpz_abs_add_bit(d, bit_index);
		else
			mpz_abs_sub_bit(d, bit_index);
	}
}

void mpz_clrbit(mpz_t d, mp_bitcnt_t bit_index)
{
	if (mpz_tstbit(d, bit_index)) {
		if (d->_mp_size >= 0)
			mpz_abs_sub_bit(d, bit_index);
		else
			mpz_abs_add_bit(d, bit_index);
	}
}

void mpz_combit(mpz_t d, mp_bitcnt_t bit_index)
{
	if (mpz_tstbit(d, bit_index) ^ (d->_mp_size < 0))
		mpz_abs_sub_bit(d, bit_index);
	else
		mpz_abs_add_bit(d, bit_index);
}

void mpz_com(mpz_t r, const mpz_t u)
{
	mpz_neg(r, u);
	mpz_sub_ui(r, r, 1);
}

void mpz_and(mpz_t r, const mpz_t u, const mpz_t v)
{
	mp_size_t un, vn, rn, i;
	mp_ptr up, vp, rp;

	mp_limb_t ux, vx, rx;
	mp_limb_t uc, vc, rc;
	mp_limb_t ul, vl, rl;

	un = GMP_ABS(u->_mp_size);
	vn = GMP_ABS(v->_mp_size);
	if (un < vn) {
		MPZ_SRCPTR_SWAP(u, v);
		MP_SIZE_T_SWAP(un, vn);
	}
	if (vn == 0) {
		r->_mp_size = 0;
		return;
	}

	uc = u->_mp_size < 0;
	vc = v->_mp_size < 0;
	rc = uc & vc;

	ux = -uc;
	vx = -vc;
	rx = -rc;

	/* If the smaller input is positive, higher limbs don't matter. */
	rn = vx ? un : vn;

	rp = MPZ_REALLOC(r, (int)(rn + rc));

	up = u->_mp_d;
	vp = v->_mp_d;

	i = 0;
	do {
		ul = (up[i] ^ ux) + uc;
		uc = ul < uc;

		vl = (vp[i] ^ vx) + vc;
		vc = vl < vc;

		rl = ((ul & vl) ^ rx) + rc;
		rc = rl < rc;
		rp[i] = rl;
	}
	while (++i < vn);
	assert(vc == 0);

	for (; i < rn; i++) {
		ul = (up[i] ^ ux) + uc;
		uc = ul < uc;

		rl = ((ul & vx) ^ rx) + rc;
		rc = rl < rc;
		rp[i] = rl;
	}
	if (rc)
		rp[rn++] = rc;
	else
		rn = mpn_normalized_size(rp, rn);

	r->_mp_size = rx ? -rn : rn;
}

void mpz_ior(mpz_t r, const mpz_t u, const mpz_t v)
{
	mp_size_t un, vn, rn, i;
	mp_ptr up, vp, rp;

	mp_limb_t ux, vx, rx;
	mp_limb_t uc, vc, rc;
	mp_limb_t ul, vl, rl;

	un = GMP_ABS(u->_mp_size);
	vn = GMP_ABS(v->_mp_size);
	if (un < vn) {
		MPZ_SRCPTR_SWAP(u, v);
		MP_SIZE_T_SWAP(un, vn);
	}
	if (vn == 0) {
		mpz_set(r, u);
		return;
	}

	uc = u->_mp_size < 0;
	vc = v->_mp_size < 0;
	rc = uc | vc;

	ux = -uc;
	vx = -vc;
	rx = -rc;

	/* If the smaller input is negative, by sign extension higher limbs
	   don't matter. */
	rn = vx ? vn : un;

	rp = MPZ_REALLOC(r, (int)(rn + rc));

	up = u->_mp_d;
	vp = v->_mp_d;

	i = 0;
	do {
		ul = (up[i] ^ ux) + uc;
		uc = ul < uc;

		vl = (vp[i] ^ vx) + vc;
		vc = vl < vc;

		rl = ((ul | vl) ^ rx) + rc;
		rc = rl < rc;
		rp[i] = rl;
	}
	while (++i < vn);
	assert(vc == 0);

	for (; i < rn; i++) {
		ul = (up[i] ^ ux) + uc;
		uc = ul < uc;

		rl = ((ul | vx) ^ rx) + rc;
		rc = rl < rc;
		rp[i] = rl;
	}
	if (rc)
		rp[rn++] = rc;
	else
		rn = mpn_normalized_size(rp, rn);

	r->_mp_size = rx ? -rn : rn;
}

void mpz_xor(mpz_t r, const mpz_t u, const mpz_t v)
{
	mp_size_t un, vn, i;
	mp_ptr up, vp, rp;

	mp_limb_t ux, vx, rx;
	mp_limb_t uc, vc, rc;
	mp_limb_t ul, vl, rl;

	un = GMP_ABS(u->_mp_size);
	vn = GMP_ABS(v->_mp_size);
	if (un < vn) {
		MPZ_SRCPTR_SWAP(u, v);
		MP_SIZE_T_SWAP(un, vn);
	}
	if (vn == 0) {
		mpz_set(r, u);
		return;
	}

	uc = u->_mp_size < 0;
	vc = v->_mp_size < 0;
	rc = uc ^ vc;

	ux = -uc;
	vx = -vc;
	rx = -rc;

	rp = MPZ_REALLOC(r, (int)(un + rc));

	up = u->_mp_d;
	vp = v->_mp_d;

	i = 0;
	do {
		ul = (up[i] ^ ux) + uc;
		uc = ul < uc;

		vl = (vp[i] ^ vx) + vc;
		vc = vl < vc;

		rl = (ul ^ vl ^ rx) + rc;
		rc = rl < rc;
		rp[i] = rl;
	}
	while (++i < vn);
	assert(vc == 0);

	for (; i < un; i++) {
		ul = (up[i] ^ ux) + uc;
		uc = ul < uc;

		rl = (ul ^ ux) + rc;
		rc = rl < rc;
		rp[i] = rl;
	}
	if (rc)
		rp[un++] = rc;
	else
		un = mpn_normalized_size(rp, un);

	r->_mp_size = rx ? -un : un;
}

mp_bitcnt_t mpz_scan1(const mpz_t u, mp_bitcnt_t starting_bit)
{
	mp_ptr up;
	mp_size_t us, un, i;
	mp_limb_t limb, ux;

	us = u->_mp_size;
	un = GMP_ABS(us);
	i = starting_bit / GMP_LIMB_BITS;

	/* Past the end there's no 1 bits for u>=0, or an immediate 1 bit
	   for u<0. Notice this test picks up any u==0 too. */
	if (i >= un)
		return (us >= 0 ? ~(mp_bitcnt_t) 0 : starting_bit);

	up = u->_mp_d;
	ux = 0;
	limb = up[i];

	if (starting_bit != 0) {
		if (us < 0) {
			ux = mpn_zero_p(up, i);
			limb = ~limb + ux;
			ux = -(mp_limb_t) (limb >= ux);
		}

		/* Mask to 0 all bits before starting_bit, thus ignoring them. */
		limb &= (GMP_LIMB_MAX << (starting_bit % GMP_LIMB_BITS));
	}

	return mpn_common_scan(limb, i, up, un, ux);
}

/* MPZ base conversion. */

size_t mpz_sizeinbase(const mpz_t u, int base)
{
	mp_size_t un;
	mp_srcptr up;
	mp_ptr tp;
	mp_bitcnt_t bits;
	struct gmp_div_inverse bi;
	size_t ndigits;

	assert(base >= 2);
	assert(base <= 36);

	un = GMP_ABS(u->_mp_size);
	if (un == 0)
		return 1;

	up = u->_mp_d;

	bits = (un - 1) * GMP_LIMB_BITS + mpn_limb_size_in_base_2(up[un - 1]);
	switch (base) {
	case 2:
		return bits;
	case 4:
		return (bits + 1) / 2;
	case 8:
		return (bits + 2) / 3;
	case 16:
		return (bits + 3) / 4;
	case 32:
		return (bits + 4) / 5;
		/* FIXME: Do something more clever for the common case of base
		   10. */
	}

	tp = gmp_xalloc_limbs(un);
	mpn_copyi(tp, up, un);
	mpn_div_qr_1_invert(&bi, base);

	ndigits = 0;
	do {
		ndigits++;
		mpn_div_qr_1_preinv(tp, tp, un, &bi);
		un -= (tp[un - 1] == 0);
	}
	while (un > 0);

	gmp_free(tp);
	return ndigits;
}
