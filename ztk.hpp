// Copyright (c) 2019, Anders Dalskov
//
// This program is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this program.  If not, see <https://www.gnu.org/licenses/>.

// Header only efficient implementation of the Ring Z/Z_2^K for K <= 128, and
// with limited support for the Galois Ring GR(2^K, 4), i.e., a degree 4
// extension of Z/Z_2^K (for K <= 128).
//
// Reacts to 2 symbol definitions.
//
// If ZTK_NO_GCC_UINT128 is defined, then operations on Z/Z_2^K elements for K >
// 64 relies on inline asm. Decrement (-=) is not currently supported in this
// case and use of this flag is not recommended.
//
// If ZTK_NO_GALOIS_RINGS is defined then no support for GR(2^K, 4) is
// provided. This (I'm guessing) will reduce code size.

#ifndef _ZTK_HPP
#define _ZTK_HPP

#include <cstdint>
#include <cstring>
#include <array>
#include <sstream>
#include <cassert>

namespace ztk {

// A Z/Z_2^K element is a vector of u64s
typedef uint64_t limb_t;

// Macros for converting to and from a length 2 array of u64 ints to a single
// __uint128 value.
#ifndef ZTK_NO_GCC_UINT128
#define LOAD_U128(r, x) do { (r) = ((__uint128_t)((x)[1]) << 64) | (x)[0]; } while (0)
#define STORE_U128(r, x) do { (r)[0] = (limb_t)(x); (r)[1] = (limb_t)((x) >> 64); } while (0)
#endif

// To circumvent C++'s rule on partial specialization, we use these helper
// functions for performing arithmetic, and provide specializations depending on
// the number of limbs that an element uses.
template<bool B> inline void mask_if(limb_t&, const limb_t);
template<size_t N> inline void op_add(limb_t[N], const limb_t[N], const limb_t[N]);
template<size_t N> inline void op_inc(limb_t[N], const limb_t[N]);
template<size_t N> inline void op_sub(limb_t[N], const limb_t[N], const limb_t[N]);
template<size_t N> inline void op_dec(limb_t[N], const limb_t[N]);
template<size_t N> inline void op_mul(limb_t[N], const limb_t[N], const limb_t[N]);

// A Z/Z_2^K element is templated by its bit-length K.
template<size_t K>
class Z2k {
public:

    template<size_t L>
    friend class Z2k;

    // Constants that show up in various places.
    static constexpr Z2k<K> Zero  = Z2k<K>((const limb_t)0, Z2k<K>::mask);
    static constexpr Z2k<K> One   = Z2k<K>(1, Z2k<K>::mask);
    static constexpr Z2k<K> Two   = Z2k<K>(2, Z2k<K>::mask);
    static constexpr Z2k<K> Three = Z2k<K>(3, Z2k<K>::mask);
    static constexpr Z2k<K> Four  = Z2k<K>(4, Z2k<K>::mask);
    static constexpr Z2k<K> Five  = Z2k<K>(5, Z2k<K>::mask);

    // Element size in bits.
    //
    // @return size in bits (i.e., the template parameter).
    static constexpr size_t SizeInBits() {
	return K;
    };

    // Element size in bytes.
    //
    // The return value is the closest multiple of 8 that fits K bits.
    //
    // @return size of element in bytes.
    static constexpr size_t SizeInBytes() {
	return byte_size;
    };

    // Element size in limbs.
    //
    // @return closest multiple of sizeof(limb_t) that fits K bits.
    static constexpr size_t SizeInLimbs() {
	return number_of_limbs;
    };

    // Indicates whether elements need masking.
    //
    // If K does not fit exactly into sizeof(limb_t) bytes, then we need to mask
    // it after e.g., arithmetic operations. Masking corresponds to performing a
    // reduction modulo 2^K.
    //
    // For elements that fit into a limb_t value exactly, no explicit reduction
    // is needed.
    //
    // @return true if an explicit reduction needs to be computed.
    static constexpr bool NeedsMasking() {
	return mask + 1;
    };

    // The mask for a particular K.
    //
    // @return a mask (i.e., a string "000 ... 000111 ... 111").
    static constexpr limb_t GetMask() {
	return mask;
    };

    // Default constructor.
    //
    // Constructs an instance of the additive identity (i.e., 0).
    constexpr Z2k() {};

    // Constructor for 1 limb with explicit mask.
    //
    // Constructs the element `x & mask'.
    //
    // @param x a limb.
    // @param mask mask to be applied on x.
    constexpr Z2k(const limb_t x, const limb_t mask) {
	limbs[0] = x & mask;
    };

    // Constructor for array of limbs with explicit mask.
    //
    // Construct the element
    //
    //   x_1 x_2 x_3 ... x_N & mask
    //
    // where x_i is the i'th entry in `limbs'.
    //
    // @param limbs an array of limbs
    // @param mask mask to be applied to the last element of `limbs'.
    constexpr Z2k(const limb_t limbs[SizeInLimbs()], const limb_t mask) {
	for (size_t i = 0; i < SizeInLimbs(); i++)
	    this->limbs[i] = limbs[i];
	this->limbs[SizeInLimbs()-1] &= mask;
    };

    // Constructor for 1 limb.
    //
    // Will be masked to fit into K bits.
    //
    // @param x a limb
    Z2k(const limb_t x) {
	limbs[0] = x;
	mask_if<NeedsMasking()>(limbs[0], mask);
    };

    // Constructor for array of limbs.
    //
    // Constructs the element
    //
    //  x_1 x_2 x_3 ... x_N & GetMask()
    //
    // where x_i is the i'th element of `limbs'.
    //
    // @param limbs array of limbs
    Z2k(const limb_t limbs[SizeInLimbs()]) {
	for (size_t i = 0; i < SizeInLimbs(); i++)
	    this->limbs[i] = limbs[i];
	mask_if<NeedsMasking()>(this->limbs[SizeInLimbs()-1], mask);
    };

    // Copy constructor
    //
    // Constructs a new element identical to `other'.
    //
    // @param other a Z2k object with K bits.
    Z2k(const Z2k<K> &other)
	: Z2k{other.limbs} {};

    // Copy constructor from different K
    //
    // Constructs a new K-bit ring element from an element `other' of L-bits. If
    // L > K, then `other' will be truncated by L-K bits when masking is
    // applied. Otherwise (L <= K) the new element and `other' will be
    // identical.
    //
    // @param other a Z2k element with L bits.
    template<size_t L>
    Z2k(const Z2k<L> &other) {
	const auto nbytes = L > K ? SizeInBytes() : other.SizeInBytes();
	memcpy(this->limbs, other.limbs, nbytes);
	mask_if<NeedsMasking()>(limbs[SizeInLimbs()-1], mask);
    };

    // Constructor from sequence of bytes.
    //
    // Exactly `Z2k<K>::SizeInBytes()' bytes of buf will be read.
    //
    // @param a pointer to at least `Z2k<K>::SizeInBytes()' bytes.
    Z2k(const unsigned char *buf)
	: Z2k{(const limb_t *)buf} {};

    // TODO(not implemented)
    template<size_t Base>
    Z2k(const std::string &str);

    // Assignment.
    //
    // Assigns `this' element with the value of `x'.
    //
    // @param x the element whose value will be compied into `this'.
    //
    // @return reference to `this' with the value of `x'.
    Z2k<K>& operator=(const Z2k<K> &x) {
	// memcpy(this->limbs, x.limbs, SizeInBytes());
	for (size_t i = 0; i < SizeInLimbs(); i++)
	    this->limbs[i] = x.limbs[i];
	return *this;
    };

    // Addition.
    //
    // @param x, y K-bit Z2k objects.
    //
    // @return x + y mod 2^K
    friend Z2k<K> operator+(const Z2k<K> &x, const Z2k<K> &y) {
	Z2k<K> r;
	op_add<SizeInLimbs()>(r.limbs, x.limbs, y.limbs);
	mask_if<NeedsMasking()>(r.limbs[SizeInLimbs()-1], mask);
	return r;
    };

    // Increment.
    //
    // Computes `this + x mod 2^K' and assigns the result to `this'.
    //
    // @param x a K-bit Z2k objects.
    //
    // @return `this' with the value `this + x mod 2^K'.
    Z2k<K> operator+=(const Z2k<K> &x) {
	op_inc<SizeInLimbs()>(limbs, x.limbs);
	mask_if<NeedsMasking()>(limbs[SizeInLimbs()-1], mask);
	return *this;
    };

    // Subtraction.
    //
    // Computes `x - y mod 2^K'.
    //
    // @param x, y K-bit Z2k objects.
    //
    // @return x - y mod 2^K.
    friend Z2k<K> operator-(const Z2k<K> &x, const Z2k<K> &y) {
    	Z2k<K> r;
	op_sub<SizeInLimbs()>(r.limbs, x.limbs, y.limbs);
	mask_if<NeedsMasking()>(r.limbs[SizeInLimbs()-1], mask);
	return r;
    };

    // Decrement.
    //
    // Compute `this - x mod 2^K' storing the result in `this'
    //
    // @param x a K-bit Z2k objects.
    //
    // @return `this' with the new value `this - x mod 2^k'.
    Z2k<K> operator-=(const Z2k<K> &x) {
	op_dec<SizeInLimbs()>(limbs, x.limbs);
	mask_if<NeedsMasking()>(limbs[SizeInLimbs()-1], mask);
	return *this;
    };

    // Negation.
    //
    // Equivalent to `0 - x mod 2^K'
    //
    // @param x a K-bit Z2k objects.
    //
    // @return -x mod 2^K.
    friend Z2k<K> operator-(const Z2k<K> &x) {
	Z2k<K> r;
	op_sub<SizeInLimbs()>(r.limbs, Z2k<K>::Zero.limbs, x.limbs);
	mask_if<NeedsMasking()>(r.limbs[SizeInLimbs()-1], mask);
	return r;
    };

    // Multiplication.
    //
    // Computes `x * y mod 2^K'.
    //
    // @param x, y K-bit Z2k objects.
    //
    // @return `x * y mod 2^K'.
    friend Z2k<K> operator*(const Z2k<K> &x, const Z2k<K> &y) {
	Z2k<K> r;
	op_mul<SizeInLimbs()>(r.limbs, x.limbs, y.limbs);
	mask_if<NeedsMasking()>(r.limbs[SizeInLimbs()-1], mask);
	return r;
    };

    // Multiply-assign
    //
    // Compute `this * x mod 2^K' storing the result in `this'.
    //
    // @param x a K-bit Z2k objects.
    //
    // @return this with value `this * x mod 2^K'.
    Z2k<K> operator*=(const Z2k<K> &x) {
	op_mul<SizeInLimbs()>(limbs, limbs, x.limbs);
	mask_if<NeedsMasking()>(limbs[SizeInLimbs()-1], mask);
	return *this;
    };

    // Equality.
    //
    // @param x a K-bit Z2k object against which to check equality.
    //
    // @return true if `this' represents the same value as `x'.
    bool operator==(const Z2k<K> &x) const {
	limb_t m = 0;
	for (size_t i = 0; i < SizeInLimbs(); i++) {
	    m |= limbs[i] ^ x.limbs[i];
	}
	return m == 0;
    };

    // Inequality.
    //
    // @param x a K-bit Z2k object against which to check inequality.
    //
    // @return true if `this' and `x' represent different values.
    bool operator!=(const Z2k<K> &x) const {
	return !(*this == x);
    };

    // Parity (odd).
    //
    // @return true if LSB of this element is 1.
    inline bool IsOdd() const {
	return (bool)(limbs[0] & 1);
    };

    // Parity (even).
    //
    // @return true if LSB of this element is 0.
    inline bool IsEven() const {
	return !IsOdd();
    };

    // Check an inverse of `this' exists.
    //
    // Recall that an inverse `x mod L' exists if `gcd(x, L) = 1'. Since L = 2^K
    // in our case, inverses exists for all elements that are odd.
    //
    // @return true if there exists an x s.t. `x * this = 1 mod 2^K'.
    inline bool IsInvertible() const {
	return IsOdd();
    }

    // Compute the inverse of `this'.
    //
    // This function computes an inverse of `this' modulo a power of 2 in
    // constant time using a technique described in
    //
    //  https://marc-b-reynolds.github.io/math/2017/09/18/ModInverse.html
    //
    // Uses at least 7 multiplications and at most 11, depending on K.
    //
    // This function fails with an assertion error if `this' is not invertible.
    //
    // @return x such that `this * x = 1 mod 2^K'.
    Z2k<K> Invert() const {
	assert (this->IsInvertible());

	// See https://marc-b-reynolds.github.io/math/2017/09/18/ModInverse.html
	// z = 3*x ^ 2
	const auto x = *this;
	auto z = x * Z2k<K>::Three;
	z.limbs[0] ^= 2;

	z = z * (Z2k<K>::Two - x * z);  // 10 bits
	z = z * (Z2k<K>::Two - x * z);  // 20 bits
	z = z * (Z2k<K>::Two - x * z);  // 40 bits
	if (K <= 40)
	    return z;
	z = z * (Z2k<K>::Two - x * z);  // 80 bits
	if (K <= 80)
	    return z;
	z = z * (Z2k<K>::Two - x * z);  // 160 bits

	return z;
    };

    // Division.
    //
    // Computes `x / y = x * y^-1 = 1 mod 2^K'. Notice y must be invertible for
    // this operation to be well defined. If `y.IsInvertible()' returns false,
    // this funcation fails with an assertion error.
    //
    // @param x, y K-bit Z2k objects
    //
    // @return x/y mod 2^K.
    friend Z2k<K> operator/(const Z2k<K> &x, const Z2k<K> &y) {
    	return x * y.Invert();
    };

    // Packs `this' into a buf.
    //
    // @param buf pointer to at least Z2k<K>::SizeInBytes() available memory.
    void Pack(unsigned char *buf) const {
	memcpy(buf, this->limbs, SizeInBytes());
    };

    // Apply a function on this element as a raw buffer.
    //
    // @param f the function to apply. This function is called with limbs, the
    //        size of this element in bytes and arg1
    // @param arg1 additional argument that is passed to f. Can be used to
    //        circumvent issues related to captures.
    void Apply(void (*f)(unsigned char*, size_t, void*), void *arg1) {
	f((unsigned char *)limbs, SizeInBytes(), arg1);
    };

    void Apply(void (*f)(unsigned char*, size_t)) {
	f((unsigned char *)limbs, SizeInBytes());
    };

    // Computes a string representation of `this'.
    //
    // An element is represented in string form as
    //
    //  "{x_1, x_2, ..., x_N}"
    //
    // where x_i is the i'th limb.
    //
    // @return a string representation of `this'.
    std::string ToString() const;

    // Support for cout << x syntax when x is a Z2k object.
    //
    // @param os an output stream
    // @param x an L-bit Z2k object.
    template<size_t L>
    friend std::ostream& operator<<(std::ostream &os, const Z2k<L> &x);

#ifdef TESTING

    // For testing, it is helpful to be able to read each individual limb of
    // `this'. Only exposed if compiled with -DTESTING.
    const limb_t *GetLimbs() const {
	return limbs;
    };

#endif

private:

    static constexpr size_t byte_size = (K + 7) / 8;
    static constexpr size_t limb_size = sizeof(limb_t);
    static constexpr size_t limb_size_bits = limb_size * 8;
    static constexpr size_t number_of_limbs = (byte_size + (limb_size - 1)) / limb_size;
    static constexpr limb_t mask =
	uint64_t(-1LL) >> (limb_size_bits - 1 - (K - 1) % limb_size_bits);

    limb_t limbs[number_of_limbs] = {0};
};

template<size_t K> const Z2k<K> Z2k<K>::Zero;
template<size_t K> const Z2k<K> Z2k<K>::One;
template<size_t K> const Z2k<K> Z2k<K>::Two;
template<size_t K> const Z2k<K> Z2k<K>::Three;
template<size_t K> const Z2k<K> Z2k<K>::Four;
template<size_t K> const Z2k<K> Z2k<K>::Five;

template<>
inline void mask_if<true>(limb_t &r, const limb_t mask) {
    r &= mask;
}

template<>
inline void mask_if<false>(limb_t &r, const limb_t mask) {
    (void)r;
    (void)mask;
}

template<>
inline void op_add<1>(limb_t r[1], const limb_t x[1], const limb_t y[1]) {
    r[0] = x[0] + y[0];
}

template<>
inline void op_add<2>(limb_t r[2], const limb_t x[2], const limb_t y[2]) {
#ifdef ZTK_NO_GCC_UINT128
    asm ("movq	%3, %1 \n\t"
    	 "movq	%2, %0 \n\t"
    	 "addq	%5, %1 \n\t"
    	 "adcq	%4, %0"
    	 : "+r" (r[1]), "+r" (r[0])
    	 : "r" (x[1]), "r" (x[0]), "r" (y[1]), "r" (y[0]) : "cc"
    	);
#else
    __uint128_t _x, _r;
    LOAD_U128(_x, x);
    LOAD_U128(_r, y);
    _r += _x;
    STORE_U128(r, _r);
#endif
}

template<>
inline void op_inc<1>(limb_t r[1], const limb_t x[1]) {
    r[0] += x[0];
}

template<>
inline void op_inc<2>(limb_t r[2], const limb_t x[2]) {
#ifdef ZTK_NO_GCC_UINT128
    asm ("addq %3, %1 \n\t"
	 "adcq %2, %0 \n\t"
	 : "+r" (r[1]), "+r" (r[0])
	 : "r" (x[1]), "r" (x[0]) : "cc"
	);
#else
    __uint128_t _x, _r;
    LOAD_U128(_x, x);
    LOAD_U128(_r, r);
    _r += _x;
    STORE_U128(r, _r);
#endif
}

template<>
inline void op_sub<1>(limb_t r[1], const limb_t x[1], const limb_t y[1]) {
    r[0] = x[0] - y[0];
}

template<>
inline void op_sub<2>(limb_t r[2], const limb_t x[2], const limb_t y[2]) {
#ifdef ZTK_NO_GCC_UINT128
    asm ("movq %3, %1 \n\t"
	 "movq %2, %0 \n\t"
	 "subq %5, %1 \n\t"
	 "sbbq %4, %0"
	 : "+r" (r[1]), "+r" (r[0])
	 : "r" (x[1]), "r" (x[0]), "r" (y[1]), "r" (y[0]) : "cc"
	);
#else
    __uint128_t _y, _r;
    LOAD_U128(_y, y);
    LOAD_U128(_r, x);
    _r -= _y;
    STORE_U128(r, _r);
#endif
}

template<>
inline void op_dec<1>(limb_t r[1], const limb_t x[1]) {
    r[0] -= x[0];
}

template<>
inline void op_dec<2>(limb_t r[2], const limb_t x[2]) {
#ifdef ZTK_NO_GCC_UINT128
    // TODO(Implement)
    throw std::runtime_error("sub of 2 limbs not supported without __uint128");
#else
    __uint128_t _x, _r;
    LOAD_U128(_x, x);
    LOAD_U128(_r, r);
    _r -= _x;
    STORE_U128(r, _r);
#endif
}

template<>
inline void op_mul<1>(limb_t r[1], const limb_t x[1], const limb_t y[1]) {
    r[0] = x[0] * y[0];
}

template<>
inline void op_mul<2>(limb_t r[2], const limb_t x[2], const limb_t y[2]) {
#ifdef ZTK_NO_GCC_UINT128
    asm ("mulx  %5, %1, %%r8     \n\t"
         "mulx  %4, %%r9, %%r10  \n\t"
         "addq  %%r9, %%r8       \n\t"
         "adcq  $0, %%r10        \n\t"
         "movq  %2, %%rdx        \n\t"
         "mulx  %5, %5, %%r9     \n\t"
         "adcq  %5, %%r8         \n\t"
         "movq  %%r8, %0         \n\t"
	 : "+r" (r[1]), "+r" (r[0])
	 : "r" (x[1]), "d" (x[0]), "r" (y[1]), "r" (y[0])
	 : "cc", "%r8", "%r9", "%r10"
	);
#else
    __uint128_t _x, _r;
    LOAD_U128(_x, x);
    LOAD_U128(_r, y);
    _r *= _x;
    STORE_U128(r, _r);
#endif
}

template<size_t K>
std::string Z2k<K>::ToString() const {
    std::string s = "{";
    for (size_t i = 0; i < SizeInLimbs()-1; i++)
	s += std::to_string(limbs[i]) + ", ";
    s += std::to_string(limbs[SizeInLimbs()-1]) + "}";
    return s;
}

template<size_t K>
std::ostream& operator<<(std::ostream &os, const Z2k<K> &x) {
    os << x.ToString();
    return os;
}

#ifndef ZTK_NO_GALOIS_RINGS

////////////////////////////////////////////////////////////////////////////////
// Galois rings
////////////////////////////////////////////////////////////////////////////////

// A galois ring is represented as an array of Z2k elements
template<size_t K, size_t D>
using gr_coeff = std::array<Z2k<K>, D>;

template<size_t K, size_t D>
class GR {
public:

    // Since a GR element is essentially a vector of Z_2^K elements, the bit
    // size is K*D.
    static constexpr size_t SizeInBits() {
	return K * D;
    };

    static constexpr size_t Degree() {
	return D;
    };

    static constexpr size_t CoeffSizeInBytes() {
	return Z2k<K>::SizeInBytes();
    };

    static constexpr size_t CoeffSizeInBits() {
	return K;
    };

    // Size of element in bytes
    static constexpr size_t SizeInBytes() {
	return Z2k<K>::SizeInBytes() * D;
    };

    constexpr GR() {};

    GR(const gr_coeff<K, D> &coeff) : coeff{coeff} {};

    GR(const GR<K, D> &x) : coeff{x.coeff} {};

    GR(const Z2k<K> &x) {
	coeff[0] = x;
    };

    Z2k<K> Project() const {
	return this->coeff[0];
    };

    GR<K, D>& operator=(const GR<K, D> &x) {
	for (size_t i = 0; i < D; i++)
	    this->coeff[i] = x.coeff[i];
	return *this;
    };

    friend GR<K, D> operator+(const GR<K, D> &x, const GR<K, D> &y) {
	gr_coeff<K, D> rcoeff;
	for (size_t i = 0; i < D; i++)
	    rcoeff[i] = x.coeff[i] + y.coeff[i];
	return GR<K, D>{rcoeff};
    };

    GR<K, D> operator+=(const GR<K, D> &x) {
	for (size_t i = 0; i < D; i++)
	    this->coeff[i] += x.coeff[i];
	return GR<K, D>{coeff};
    };

    friend GR<K, D> operator-(const GR<K, D> &x, const GR<K, D> &y) {
	gr_coeff<K, D> rcoeff;
	for (size_t i = 0; i < D; i++)
	    rcoeff[i] = x.coeff[i] - y.coeff[i];
	return GR<K, D>{rcoeff};
    };

    GR<K, D> operator-=(const GR<K, D> &x) {
	for (size_t i = 0; i < D; i++)
	    this->coeff[i] -= x.coeff[i];
	return GR<K, D>{coeff};
    };

    bool operator==(const GR<K, D> &x) const {
    	bool b = true;
    	for (size_t i = 0; i < D; i++)
    	    b &= this->coeff[i] == x.coeff[i];
    	return b;
    };

    bool operator !=(const GR<K, D> &x) const {
    	return !(*this == x);
    };

    bool IsInvertible() const {
	bool b = false;
	for (const auto &x: coeff)
	    b = x.IsOdd() or b;
	return b;
    };

    void Pack(unsigned char *buf) const {
	auto p = buf;
	for (size_t i = 0; i < D; i++) {
	    coeff[i].Pack(p);
	    p += Z2k<K>::SizeInBytes();
	}
    };

    // Apply f on each coefficient
    void Apply(void (*f)(unsigned char*, size_t, void*), void *arg1) {
	for (auto &c : coeff)
	    c.Apply(f, arg1);
    };

    void Apply(void (*f)(unsigned char*, size_t)) {
	for (auto &c : coeff)
	    c.Apply(f);
    };

    std::string ToString() const;

    template<size_t L, size_t H>
    friend std::ostream& operator<<(std::ostream &os, const GR<K, D> &x);

#ifdef TESTING

    const gr_coeff<K, D> GetCoeff() const {
	return this->coeff;
    };

#endif

protected:

    gr_coeff<K, D> coeff;

};

template<size_t K>
class GR4 : public GR<K, 4> {
public:

    GR4() {};
    GR4(const gr_coeff<K, 4> &coeff) : GR<K, 4>{coeff} {};
    GR4(const GR<K, 4> &x) : GR<K, 4>{x} {};

    friend GR4<K> operator*(const GR4<K> &x, const GR4<K> &y) {

    	const auto v = x.coeff;
    	const auto u = y.coeff;

    	auto v0 = v[0]; auto v1 = v[1]; auto v2 = v[2]; auto v3 = v[3];
    	auto u0 = u[0]; auto u1 = u[1]; auto u2 = u[2]; auto u3 = u[3];

    	gr_coeff<K, 4> r;

    	r[0] = v0*u0 - v3*u1 - v2*u2 - v1*u3;
    	r[1] = v1*u0 + v0*u1 - v3*(u1 + u2) - v2*u2 - v2*u3 - v1*u3;
    	r[2] = v2*u0 - v2*u3 + v1*u1 + u2*(v0 - v3) - v3*u3;
    	r[3] = v3*u0 + v2*u1 + v1*u2 + v0*u3 - v3*u3;

    	return GR4<K>{r};
    };

    GR4<K> Invert() const {
	assert (this->IsInvertible());
    	gr_coeff<K, 4> r;

	const auto v0 = this->coeff[0]; const auto v1 = this->coeff[1];
	const auto v2 = this->coeff[2]; const auto v3 = this->coeff[3];

	auto d = Denom().Invert();

	r[0] = (v0*v0*v0 - v1*v1*v1 + Z2k<K>::Three*v0*v1*v2 - v1*v1*v2 + v0*v2*v2
		+ v2*v2*v2 - Z2k<K>::Three*v0*v0*v3 + Z2k<K>::Two*v0*v1*v3
		- Z2k<K>::Three*v1*v2*v3 + v2*v2*v3 + Z2k<K>::Three*v0*v3*v3
		- Z2k<K>::Two*v1*v3*v3 + v2*v3*v3 - v3*v3*v3)*d;

	r[1] = (-d)*(v0*v0*v1 - v1*v2*v2 + v1*v1*v3 + v3*v3*v3 + v0*(v2*v2 - v1*v3 + Z2k<K>::Two*v2*v3));

	r[2] = d*(-v0*v0*v2 - v2*v2*v2 + Z2k<K>::Two*v1*v2*v3 + v3*v3*v3 + v0*(v1*v1 + (v2 - v3)*v3));

	r[3] = (-d)*(v1*v1*v1 - Z2k<K>::Two*v0*v1*v2 - v2*v2*v2 - v2*v2*v3
		     + (v0 - v3)*(v0-v3)*v3 + v1*v3*(Z2k<K>::Three*v2 + v3));

	return GR4<K>{r};
    };

    friend GR4<K> operator/(const GR4<K> &x, const GR4<K> &y) {
	auto y_ = y.Invert();
	return x * y_;
    };

private:

    Z2k<K> Denom() const {
	const auto v0 = this->coeff[0]; const auto v1 = this->coeff[1];
	const auto v2 = this->coeff[2]; const auto v3 = this->coeff[3];

	return v0*v0*v0*v0 + v1*v1*v1*v1 - v1*v2*v2*v2 + v2*v2*v2*v2
	    - Z2k<K>::Three*v0*v0*v0*v3 + v1*(Z2k<K>::Three*v1 - Z2k<K>::Four*v2)*v2*v3
	    + Z2k<K>::Two*v1*v1*v3*v3 + (v1 - v2)*v3*v3*v3 + v3*v3*v3*v3
	    + v0*v0*(Z2k<K>::Three*v1*v2 + Z2k<K>::Two*v2*v2 + Z2k<K>::Four*v1*v3 + Z2k<K>::Three*v3*v3)
	    - v0*(v1*v1*v1 + Z2k<K>::Four*v1*v1*v2 - v2*v2*v2 + (Z2k<K>::Three*v1 - v2)*v2*v3
		  + (Z2k<K>::Five*v1 - Z2k<K>::Four*v2)*v3*v3 + v3*v3*v3);
    };
};

template<size_t K, size_t D>
std::string GR<K, D>::ToString() const {
    std::stringstream ss;
    ss << "{";
    for (size_t i = 0; i < D - 1; i++)
	ss << this->coeff[i] << ", ";
    ss << this->coeff[D - 1] << "}";
    return ss.str();
}

template<size_t K, size_t D>
std::ostream& operator<<(std::ostream &os, const GR<K, D> &x) {
    os << x.ToString();
    return os;
}

#endif // ZTK_NO_GALOIS_RINGS

} // ztk

#endif // _ZTK_HPP
