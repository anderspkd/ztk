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

// ztk.hpp provides a header only implemention for the ring Z/Z_2^K for K <=
// 128.

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

namespace {

// To circumvent C++'s rule on partial specialization, we use these helper
// functions for performing arithmetic, and provide specializations depending on
// the number of limbs that an element uses.
template<bool B> inline void mask_if(limb_t&, const limb_t);
template<std::size_t N> inline void op_add(limb_t[N], const limb_t[N], const limb_t[N]);
template<std::size_t N> inline void op_inc(limb_t[N], const limb_t[N]);
template<std::size_t N> inline void op_sub(limb_t[N], const limb_t[N], const limb_t[N]);
template<std::size_t N> inline void op_dec(limb_t[N], const limb_t[N]);
template<std::size_t N> inline void op_mul(limb_t[N], const limb_t[N], const limb_t[N]);

}

// A Z/Z_2^K element is templated by its bit-length K.
template<std::size_t K>
class Z2k {
public:

    template<std::size_t L>
    friend class Z2k;

    // Constants that show up in various places.
    static constexpr Z2k<K> zero  = Z2k<K>((const limb_t)0, Z2k<K>::_mask);
    static constexpr Z2k<K> one   = Z2k<K>(1, Z2k<K>::_mask);
    static constexpr Z2k<K> two   = Z2k<K>(2, Z2k<K>::_mask);
    static constexpr Z2k<K> three = Z2k<K>(3, Z2k<K>::_mask);
    static constexpr Z2k<K> four  = Z2k<K>(4, Z2k<K>::_mask);
    static constexpr Z2k<K> five  = Z2k<K>(5, Z2k<K>::_mask);

    // Element size in bits.
    //
    // @return size in bits (i.e., the template parameter).
    static constexpr std::size_t bit_size() {
	return K;
    };

    // Element size in bytes.
    //
    // The return value is the closest multiple of 8 that fits K bits.
    //
    // @return size of element in bytes.
    static constexpr std::size_t byte_size() {
	return _byte_size;
    };

    // Element size in limbs.
    //
    // @return closest multiple of sizeof(limb_t) that fits K bits.
    static constexpr std::size_t size() {
	return _number_of_limbs;
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
    static constexpr bool needs_masking() {
	return _mask + 1;
    };

    // The mask for a particular K.
    //
    // @return a mask (i.e., a string "000 ... 000111 ... 111").
    static constexpr limb_t mask() {
	return _mask;
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
	_limbs[0] = x & mask;
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
    constexpr Z2k(const limb_t limbs[size()], const limb_t mask) {
	for (std::size_t i = 0; i < size(); i++)
	    this->_limbs[i] = limbs[i];
	this->_limbs[size()-1] &= mask;
    };

    // Constructor for 1 limb.
    //
    // Will be masked to fit into K bits.
    //
    // @param x a limb
    Z2k(const limb_t x) {
	_limbs[0] = x;
	mask_if<needs_masking()>(_limbs[0], _mask);
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
    Z2k(const limb_t limbs[size()]) {
	for (std::size_t i = 0; i < size(); i++)
	    this->_limbs[i] = limbs[i];
	mask_if<needs_masking()>(this->_limbs[size()-1], _mask);
    };

    // Copy constructor
    //
    // Constructs a new element identical to `other'.
    //
    // @param other a Z2k object with K bits.
    Z2k(const Z2k<K> &other)
	: Z2k{other._limbs} {};

    // Copy constructor from different K
    //
    // Constructs a new K-bit ring element from an element `other' of L-bits. If
    // L > K, then `other' will be truncated by L-K bits when masking is
    // applied. Otherwise (L <= K) the new element and `other' will be
    // identical.
    //
    // @param other a Z2k element with L bits.
    template<std::size_t L>
    Z2k(const Z2k<L> &other) {
	const auto nbytes = L > K ? byte_size() : other.byte_size();
	memcpy(this->_limbs, other._limbs, nbytes);
	mask_if<needs_masking()>(this->_limbs[size()-1], _mask);
    };

    // Constructor from sequence of bytes.
    //
    // Exactly `Z2k<K>::SizeInBytes()' bytes of buf will be read.
    //
    // @param a pointer to at least `Z2k<K>::SizeInBytes()' bytes.
    Z2k(const unsigned char *buf)
	: Z2k{(const limb_t *)buf} {};

    // TODO(not implemented)
    template<std::size_t B>
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
	for (std::size_t i = 0; i < size(); i++)
	    this->_limbs[i] = x._limbs[i];
	return *this;
    };

    // Addition.
    //
    // @param x, y K-bit Z2k objects.
    //
    // @return x + y mod 2^K
    friend Z2k<K> operator+(const Z2k<K> &x, const Z2k<K> &y) {
	Z2k<K> r;
	op_add<size()>(r._limbs, x._limbs, y._limbs);
	mask_if<needs_masking()>(r._limbs[size()-1], _mask);
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
	op_inc<size()>(_limbs, x._limbs);
	mask_if<needs_masking()>(_limbs[size()-1], _mask);
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
	op_sub<size()>(r._limbs, x._limbs, y._limbs);
	mask_if<needs_masking()>(r._limbs[size()-1], _mask);
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
	op_dec<size()>(_limbs, x._limbs);
	mask_if<needs_masking()>(_limbs[size()-1], _mask);
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
	op_sub<size()>(r._limbs, Z2k<K>::zero._limbs, x._limbs);
	mask_if<needs_masking()>(r._limbs[size()-1], _mask);
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
	op_mul<size()>(r._limbs, x._limbs, y._limbs);
	mask_if<needs_masking()>(r._limbs[size()-1], _mask);
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
	op_mul<size()>(_limbs, _limbs, x._limbs);
	mask_if<needs_masking()>(_limbs[size()-1], _mask);
	return *this;
    };

    // Equality.
    //
    // @param x a K-bit Z2k object against which to check equality.
    //
    // @return true if `this' represents the same value as `x'.
    bool operator==(const Z2k<K> &x) const {
	limb_t m = 0;
	for (std::size_t i = 0; i < this->size(); i++) {
	    m |= _limbs[i] ^ x._limbs[i];
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
    inline bool is_odd() const {
	return (bool)(_limbs[0] & 1);
    };

    // Parity (even).
    //
    // @return true if LSB of this element is 0.
    inline bool is_even() const {
	return !is_odd();
    };

    // Check an inverse of `this' exists.
    //
    // Recall that an inverse `x mod L' exists if `gcd(x, L) = 1'. Since L = 2^K
    // in our case, inverses exists for all elements that are odd.
    //
    // @return true if there exists an x s.t. `x * this = 1 mod 2^K'.
    inline bool is_invertible() const {
	return is_odd();
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
    Z2k<K> invert() const {
	assert (this->is_invertible());

	// See https://marc-b-reynolds.github.io/math/2017/09/18/ModInverse.html
	// z = 3*x ^ 2
	const auto x = *this;
	auto z = x * Z2k<K>::three;
	z._limbs[0] ^= 2;

	z = z * (Z2k<K>::two - x * z);  // 10 bits
	z = z * (Z2k<K>::two - x * z);  // 20 bits
	z = z * (Z2k<K>::two - x * z);  // 40 bits
	if (K <= 40)
	    return z;
	z = z * (Z2k<K>::two - x * z);  // 80 bits
	if (K <= 80)
	    return z;
	z = z * (Z2k<K>::two - x * z);  // 160 bits

	return z;
    };

    // Division.
    //
    // Computes `x / y = x * y^-1 mod 2^K'. Notice y must be invertible for this
    // operation to be well defined. If `y.is_invertible()' returns false, this
    // funcation fails with an assertion error.
    //
    // @param x, y K-bit Z2k objects
    //
    // @return x/y mod 2^K.
    friend Z2k<K> operator/(const Z2k<K> &x, const Z2k<K> &y) {
    	return x * y.invert();
    };

    // Packs `this' into a buf.
    //
    // @param buf pointer to at least Z2k<K>::SizeInBytes() available memory.
    void pack(unsigned char *buf) const {
	memcpy(buf, this->_limbs, byte_size());
    };

    // Return a pointer to the internal representation of this element.
    limb_t* limbs() {
	return _limbs;
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
    std::string to_string() const;

    // Support for cout << x syntax when x is a Z2k object.
    //
    // @param os an output stream
    // @param x an L-bit Z2k object.
    template<std::size_t L>
    friend std::ostream& operator<<(std::ostream &os, const Z2k<L> &x);

private:

    static constexpr std::size_t _byte_size       = (K + 7) / 8;
    static constexpr std::size_t _limb_size       = sizeof(limb_t);
    static constexpr std::size_t _limb_size_bits  = _limb_size * 8;
    static constexpr std::size_t _number_of_limbs =
        (_byte_size + (_limb_size - 1)) / _limb_size;

    static constexpr limb_t _mask =
	limb_t(-1LL) >> (_limb_size_bits - 1 - (K - 1) % _limb_size_bits);

    limb_t _limbs[_number_of_limbs] = {0};
};

template<std::size_t K> const Z2k<K> Z2k<K>::zero;
template<std::size_t K> const Z2k<K> Z2k<K>::one;
template<std::size_t K> const Z2k<K> Z2k<K>::two;
template<std::size_t K> const Z2k<K> Z2k<K>::three;
template<std::size_t K> const Z2k<K> Z2k<K>::four;
template<std::size_t K> const Z2k<K> Z2k<K>::five;

namespace {

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

}

template<std::size_t K>
std::string Z2k<K>::to_string() const {
    std::string s = "{";
    for (std::size_t i = 0; i < size()-1; i++)
	s += std::to_string(_limbs[i]) + ", ";
    s += std::to_string(_limbs[size()-1]) + "}";
    return s;
}

template<std::size_t K>
std::ostream& operator<<(std::ostream &os, const Z2k<K> &x) {
    os << x.to_string();
    return os;
}

#ifndef ZTK_NO_GALOIS_RINGS

////////////////////////////////////////////////////////////////////////////////
// Galois rings
////////////////////////////////////////////////////////////////////////////////

// A galois ring is represented as an array of Z2k elements
template<std::size_t K, std::size_t D>
using gr_coeff = std::array<Z2k<K>, D>;

namespace {

template<std::size_t K> static inline void gr_deg4_mul(gr_coeff<K, 4>&, const gr_coeff<K, 4>&, const gr_coeff<K, 4>&);
template<std::size_t K> static inline void gr_deg4_inv(gr_coeff<K, 4>&, const gr_coeff<K, 4>&);
template<std::size_t K> static inline Z2k<K> gr_deg4_den(const Z2k<K>&, const Z2k<K>&, const Z2k<K>&, const Z2k<K>&);

}

template<std::size_t K, std::size_t D>
class GR {
public:

    // GR element size in bits.
    //
    // @return element size in bits.
    static constexpr std::size_t bit_size() {
	return K * D;
    };

    // Degree of extension.
    //
    // @return degree of extension (template parameter D).
    static constexpr std::size_t degree() {
	return D;
    };

    // Coefficient size in bytes.
    //
    // @return byte size of a single Z/Z_2^K element.
    static constexpr std::size_t coeff_byte_size() {
	return Z2k<K>::byte_size();
    };

    // Coefficient size in bits.
    //
    // @return size of coefficienct in bits (template parameter K).
    static constexpr std::size_t coeff_bit_size() {
	return K;
    };

    // Size of element in bytes.
    //
    // @return Size of an element in bytes.
    static constexpr std::size_t byte_size() {
	return Z2k<K>::byte_size() * D;
    };

    // Zero element.
    //
    // @return additive identity in GR(2^K, D).
    static const GR<K, D> zero() {
	static GR<K, D> zero {Z2k<K>::zero};
	return zero;
    };

    // One element.
    //
    // @return multiplicative identity in GR(2^K, D).
    static const GR<K, D> one() {
	static GR<K, D> one {Z2k<K>::one};
	return one;
    };

    // Default constructor.
    //
    // Constructs an element with all 0 coefficients.
    constexpr GR() {};

    // Constructor from Z/Z_2^K element.
    //
    // Construct the constant polynomial P(X) = c.
    //
    // @param c a Z2k element of bit-length K.
    constexpr GR(const Z2k<K> &c) {
	_coeff[0] = c;
    };

    // Constructor from a set of coefficients.
    //
    // Constructs the polynomial
    //   P(X) = coeff[0] + coeff[1].X + ... + coeff[D-1].X^{D-1}
    //
    // @param coeff an array of D Z2k<K> elements.
    GR(const gr_coeff<K, D> &coeff)
	: _coeff{coeff}
	{};

    // Copy constructor.
    //
    // @param x a GR element.
    GR(const GR<K, D> &x)
	: _coeff{x._coeff}
	{};

    // Constructor from a buffer.
    //
    // Use data pointed to by buf to construct a new GR element.
    //
    // @param buf, a pointer to byte_size() bytes which will be interpreted as
    //        an array of Z2k<K> elements.
    GR(const unsigned char *buf) {
	std::size_t offset = 0;
	for (std::size_t i = 0; i < D; i++) {
	    _coeff[i] = Z2k<K>{buf + offset};
	    offset += Z2k<K>::byte_size();
	}
    };

    // Project the current GR element to Z/Z_2^K.
    //
    // @return constant term of GR element.
    Z2k<K> project() const {
	return this->_coeff[0];
    };

    // Assignment operator
    //
    // @param x a GR element that this is assigned to.
    //
    // @return reference to this.
    GR<K, D>& operator=(const GR<K, D> &x) {
	for (std::size_t i = 0; i < D; i++)
	    this->_coeff[i] = x._coeff[i];
	return *this;
    };

    // Addition operator.
    //
    // @param x first summand.
    // @param y second summand.
    //
    // @return x + y.
    friend GR<K, D> operator+(const GR<K, D> &x, const GR<K, D> &y) {
	gr_coeff<K, D> rcoeff;
	for (std::size_t i = 0; i < D; i++)
	    rcoeff[i] = x._coeff[i] + y._coeff[i];
	return GR<K, D>{rcoeff};
    };

    // Addition of a constant.
    //
    // Corresponds to adding x with the constant polynomial P(X) = c.
    //
    // @param x GR element
    // @param c a Z/Z_2^K constant
    //
    // @return x + P(X) where P(X) = c. I.e., c is addded to the constant term
    //         of x.
    friend GR<K, D> operator+(const GR<K, D> &x, const Z2k<K> &c) {
	gr_coeff<K, D> rcoeff {x._coeff};
	rcoeff[0] += c;
	return GR<K, D>{rcoeff};
    };

    // Addition of a constant.
    //
    // See above. (Addition commutes in GR(2^K, D).)
    friend GR<K, D> operator+(const Z2k<K> &x, const GR<K, D> &c) {
	return c + x;
    };

    // Increment.
    //
    // Add x coefficient wise to this.
    //
    // @param x a GR element to increment this with.
    //
    // @return this incremented by x.
    GR<K, D> operator+=(const GR<K, D> &x) {
	for (std::size_t i = 0; i < D; i++)
	    this->_coeff[i] += x._coeff[i];
	return *this;
    };

    // Increment by constant.
    //
    // Adds c to the constant term of this.
    //
    // @param a Z/Z_2^K constant.
    //
    // @return this.
    GR<K, D> operator+=(const Z2k<K> &c) {
	this->_coeff[0] += c;
	return *this;
    };

    // Subtraction.
    friend GR<K, D> operator-(const GR<K, D> &x, const GR<K, D> &y) {
	gr_coeff<K, D> rcoeff;
	for (std::size_t i = 0; i < D; i++)
	    rcoeff[i] = x._coeff[i] - y._coeff[i];
	return GR<K, D>{rcoeff};
    };

    friend GR<K, D> operator-(const GR<K, D> &x, const Z2k<K> &c) {
	gr_coeff<K, D> rcoeff {x._coeff};
	rcoeff[0] -= c;
	return GR<K, D>{rcoeff};
    };

    friend GR<K, D> operator-(const Z2k<K> &c, const GR<K, D> &x) {
	return -(x - c);
    };

    friend GR<K, D> operator-(const GR<K, D> &x) {
	gr_coeff<K, D> rcoeff;
	for (std::size_t i = 0; i < D; i++)
	    rcoeff[i] = -x._coeff[i];
	return GR<K, D>{rcoeff};
    };

    GR<K, D> operator-=(const GR<K, D> &x) {
	for (std::size_t i = 0; i < D; i++)
	    this->_coeff[i] -= x._coeff[i];
	return *this;
    };

    GR<K, D> operator-=(const Z2k<K> &c) {
	this->_coeff[0] -= c;
	return *this;
    };

    friend GR<K, D> operator*(const GR<K, D> &x, const GR<K, D> &y) {
	gr_coeff<K, D> rcoeff;
	if (D == 4) {
	    gr_deg4_mul(rcoeff, x._coeff, y._coeff);
	    return GR<K, D>{rcoeff};
	}
	throw std::runtime_error("multiplication of this degree not supported");
    };

    friend GR<K, D> operator*(const GR<K, D> &x, const Z2k<K> &y) {
	gr_coeff<K, D> rcoeff;
	for (std::size_t i = 0; i < D; i++)
	    rcoeff[i] = x._coeff[i] * y;
	return GR<K, D>{rcoeff};
    };

    friend GR<K, D> operator*(const Z2k<K> &x, const GR<K, D> &y) {
	return y * x;
    };

    GR<K, D> invert() const {
	assert (this->is_invertible());
	gr_coeff<K, D> rcoeff;
	if (D == 4) {
	    gr_deg4_inv(rcoeff, this->_coeff);
	    return GR<K, D>{rcoeff};
	}
	throw std::runtime_error("inversion of this degree not supported");
    };

    friend GR<K, D> operator/(const GR<K, D> &x, const GR<K, D> &y) {
	auto y_ = y.invert();
	return x * y_;
    };

    Z2k<K>& operator[](const std::size_t idx) {
	return _coeff[idx];
    };

    const Z2k<K>& operator[](const std::size_t idx) const {
	return _coeff[idx];
    };

    bool operator==(const GR<K, D> &x) const {
    	bool b = true;
    	for (std::size_t i = 0; i < D; i++)
    	    b &= this->_coeff[i] == x._coeff[i];
    	return b;
    };

    bool operator !=(const GR<K, D> &x) const {
    	return !(*this == x);
    };

    bool is_invertible() const {
	bool b = false;
	for (const auto &x: _coeff)
	    b |= x.is_odd();
	return b;
    };

    void pack(unsigned char *buf) const {
	unsigned char *p = buf;
	for (std::size_t i = 0; i < D; i++) {
	    _coeff[i].pack(p);
	    p += Z2k<K>::byte_size();
	}
    };

    // Apply f on each coefficient
    void apply(void (*f)(unsigned char*, std::size_t, void*), void *arg1) {
	for (auto &c : _coeff)
	    c.apply(f, arg1);
    };

    void apply(void (*f)(unsigned char*, std::size_t)) {
	for (auto &c : _coeff)
	    c.apply(f);
    };

    std::string to_string() const;

    template<std::size_t L, std::size_t H>
    friend std::ostream& operator<<(std::ostream &os, const GR<K, D> &x);

private:

    gr_coeff<K, D> _coeff;

};

namespace {

template<std::size_t K>
inline void gr_deg4_mul(gr_coeff<K, 4> &r, const gr_coeff<K, 4> &x, const gr_coeff<K, 4> &y) {

    auto x0 = x[0]; auto x1 = x[1]; auto x2 = x[2]; auto x3 = x[3];
    auto y0 = y[0]; auto y1 = y[1]; auto y2 = y[2]; auto y3 = y[3];

    r[0] = x0*y0 - x3*y1 - x2*y2 - x1*y3;
    r[1] = x1*y0 + x0*y1 - x3*(y1 + y2) - x2*y2 - x2*y3 - x1*y3;
    r[2] = x2*y0 - x2*y3 + x1*y1 + y2*(x0 - x3) - x3*y3;
    r[3] = x3*y0 + x2*y1 + x1*y2 + x0*y3 - x3*y3;
}

template<std::size_t K>
inline void gr_deg4_inv(gr_coeff<K, 4> &r, const gr_coeff<K, 4> &x) {
	const auto v0 = x[0]; const auto v1 = x[1];
	const auto v2 = x[2]; const auto v3 = x[3];

	auto d = gr_deg4_den(v0, v1, v2, v3).invert();

	r[0] = (v0*v0*v0 - v1*v1*v1 + Z2k<K>::three*v0*v1*v2 - v1*v1*v2 + v0*v2*v2
		+ v2*v2*v2 - Z2k<K>::three*v0*v0*v3 + Z2k<K>::two*v0*v1*v3
		- Z2k<K>::three*v1*v2*v3 + v2*v2*v3 + Z2k<K>::three*v0*v3*v3
		- Z2k<K>::two*v1*v3*v3 + v2*v3*v3 - v3*v3*v3)*d;

	r[1] = (-d)*(v0*v0*v1 - v1*v2*v2 + v1*v1*v3 + v3*v3*v3 + v0*(v2*v2 - v1*v3 + Z2k<K>::two*v2*v3));

	r[2] = d*(-v0*v0*v2 - v2*v2*v2 + Z2k<K>::two*v1*v2*v3 + v3*v3*v3 + v0*(v1*v1 + (v2 - v3)*v3));

	r[3] = (-d)*(v1*v1*v1 - Z2k<K>::two*v0*v1*v2 - v2*v2*v2 - v2*v2*v3
		     + (v0 - v3)*(v0-v3)*v3 + v1*v3*(Z2k<K>::three*v2 + v3));
}

template<std::size_t K>
inline Z2k<K> gr_deg4_den(const Z2k<K> &v0, const Z2k<K> &v1, const Z2k<K> &v2, const Z2k<K> &v3) {
    return v0*v0*v0*v0 + v1*v1*v1*v1 - v1*v2*v2*v2 + v2*v2*v2*v2
	- Z2k<K>::three*v0*v0*v0*v3 + v1*(Z2k<K>::three*v1 - Z2k<K>::four*v2)*v2*v3
	+ Z2k<K>::two*v1*v1*v3*v3 + (v1 - v2)*v3*v3*v3 + v3*v3*v3*v3
	+ v0*v0*(Z2k<K>::three*v1*v2 + Z2k<K>::two*v2*v2 + Z2k<K>::four*v1*v3 + Z2k<K>::three*v3*v3)
	- v0*(v1*v1*v1 + Z2k<K>::four*v1*v1*v2 - v2*v2*v2 + (Z2k<K>::three*v1 - v2)*v2*v3
	      + (Z2k<K>::five*v1 - Z2k<K>::four*v2)*v3*v3 + v3*v3*v3);
}

}

template<std::size_t K, std::size_t D>
std::string GR<K, D>::to_string() const {
    std::stringstream ss;
    ss << "{";
    for (std::size_t i = 0; i < D - 1; i++)
	ss << this->_coeff[i] << ", ";
    ss << this->_coeff[D - 1] << "}";
    return ss.str();
}

template<std::size_t K, std::size_t D>
std::ostream& operator<<(std::ostream &os, const GR<K, D> &x) {
    os << x.to_string();
    return os;
}

#endif // ZTK_NO_GALOIS_RINGS

} // ztk

#endif // _ZTK_HPP
