// Copyright Anders Dalskov, 2019
//
// Credit of multiplication asm goes to Diego Aranha.
//
// This header implements the ring Z/Z_2^k for k <= 128. Supports addition,
// subtraction, multiplication and division (provided denominator is odd).
//
// The header reacts to one flag:
//
// If ZTK_GCC_UINT128 is defined then arithmetic operations (except division)
// for elements with k > 64 will use the GCC __uint128 extension.

#ifndef _ZTK_HPP
#define _ZTK_HPP

#include <cstdint>
#include <string>

namespace ztk {

// An element is represented as an array of limb_t types.
typedef uint64_t limb_t;

#ifdef ZTK_GCC_UINT128
#define LOAD_U128(r, x) do { (r) = ((__uint128_t)((x)[1]) << 64) | (x)[0]; } while (0)
#define STORE_U128(r, x) do { (r)[0] = (limb_t)(x); (r)[1] = (limb_t)((x) >> 64); } while (0)
#endif

template<bool B> inline void mask_if(limb_t&, const limb_t);
template<size_t N> inline void op_add(limb_t[N], const limb_t[N], const limb_t[N]);
template<size_t N> inline void op_inc(limb_t[N], const limb_t[N]);
template<size_t N> inline void op_sub(limb_t[N], const limb_t[N], const limb_t[N]);
template<size_t N> inline void op_dec(limb_t[N], const limb_t[N]);
template<size_t N> inline void op_mul(limb_t[N], const limb_t[N], const limb_t[N]);

// Z/Z_2^k
template<size_t K>
class Z2k {
public:

    template<size_t L>
    friend class Z2k;

    // Useful constants
    static constexpr Z2k<K> Zero  = Z2k<K>((const limb_t)0, Z2k<K>::mask);
    static constexpr Z2k<K> One   = Z2k<K>(1, Z2k<K>::mask);
    static constexpr Z2k<K> Two   = Z2k<K>(2, Z2k<K>::mask);
    static constexpr Z2k<K> Three = Z2k<K>(3, Z2k<K>::mask);

    // Size of an element in bits.
    static constexpr size_t SizeInBits() {
	return K;
    };

    // Size of an element in bytes.
    static constexpr size_t SizeInBytes() {
	return byte_size;
    };

    // Size of element in # of limbs
    static constexpr size_t SizeInLimbs() {
	return number_of_limbs;
    };

    // True if the result of an arithmetic operation needs to be masked
    static constexpr bool NeedsMasking() {
	return mask + 1;
    };

    // Constructors.

    // Construct the value 0
    constexpr Z2k() {};

    // Constexpr constructor from constant.
    constexpr Z2k(const limb_t x, const limb_t mask) {
	limbs[0] = x & mask;
    };

    // Constexpr constructor from constant.
    constexpr Z2k(const limb_t x[SizeInLimbs()], const limb_t mask) {
	limbs[SizeInLimbs()-1] = x[SizeInLimbs()-1] & mask;
	for (size_t i = SizeInLimbs()-2; i > 0; i--)
	    limbs[i] = x[i];
    };

    // Construct from a single limb.
    Z2k(const limb_t x) {
	limbs[0] = x;
	mask_if<NeedsMasking()>(limbs[0], mask);
    };

    // Construct frin a set of limbs.
    Z2k(const limb_t limbs[SizeInLimbs()]) {
	memcpy(this->limbs, limbs, SizeInBytes());
	mask_if<NeedsMasking()>(this->limbs[SizeInLimbs()-1], mask);
    };

    // Copy constructor
    Z2k(const Z2k<K> &other)
	: Z2k{other.limbs} {};

    // Copy constructor with a different modulus. When L <= K then the new
    // element is equivalent with the old. In the other situation, the new
    // element corresponds to the old one with the top K - L bits set to 0.
    template<size_t L>
    Z2k(const Z2k<L> &other) {
	const auto nbytes = L > K ? SizeInBytes() : other.SizeInBytes();
	memcpy(this->limbs, other.limbs, nbytes);
	mask_if<NeedsMasking()>(limbs[SizeInLimbs()-1], mask);
    };

    // Construct a new element from a raw data. buf is assumed to point to
    // SizeInBytes bytes.
    Z2k(const unsigned char *buf)
	: Z2k{(const limb_t *)buf} {};

    // Construct from a string representation in a particular base
    // TODO(implement me)
    Z2k(const std::string &str, const size_t base = 10);

    // Operators

    // Assignment.
    Z2k<K>& operator=(const Z2k<K> &x) {
	memcpy(this->limbs, x.limbs, SizeInBytes());
	return *this;
    };

    // Addition
    friend Z2k<K> operator+(const Z2k<K> &x, const Z2k<K> &y) {
	Z2k<K> r;
	op_add<SizeInLimbs()>(r.limbs, x.limbs, y.limbs);
	mask_if<NeedsMasking()>(r.limbs[SizeInLimbs()-1], mask);
	return r;
    };

    // Increment
    Z2k<K> operator+=(const Z2k<K> &x) {
	op_inc<SizeInLimbs()>(limbs, x.limbs);
	mask_if<NeedsMasking()>(limbs[SizeInLimbs()-1], mask);
	return *this;
    };

    // Subtraction
    friend Z2k<K> operator-(const Z2k<K> &x, const Z2k<K> &y) {
    	Z2k<K> r;
	op_sub<SizeInLimbs()>(r.limbs, x.limbs, y.limbs);
	mask_if<NeedsMasking()>(r.limbs[SizeInLimbs()-1], mask);
	return r;
    };

    // Decrement
    Z2k<K> operator-=(const Z2k<K> &x) {
	op_dec<SizeInLimbs()>(limbs, x.limbs);
	mask_if<NeedsMasking()>(limbs[SizeInLimbs()-1], mask);
	return *this;
    }

    // Multiplication
    friend Z2k<K> operator*(const Z2k<K> &x, const Z2k<K> &y) {
	Z2k<K> r;
	op_mul<SizeInLimbs()>(r.limbs, x.limbs, y.limbs);
	mask_if<NeedsMasking()>(r.limbs[SizeInLimbs()-1], mask);
	return r;
    };

    // Equality test
    bool operator==(const Z2k<K> &x) const {
	limb_t m = 0;
	for (size_t i = 0; i < SizeInLimbs(); i++) {
	    m |= limbs[i] ^ x.limbs[i];
	}
	return m == 0;
    };

    // Inequality test
    bool operator!=(const Z2k<K> &x) const {
	return !(*this == x);
    };

    // Parity (odd)
    inline bool IsOdd() const {
	return (bool)(limbs[0] & 1);
    };

    // Parity (even)
    inline bool IsEven() const {
	return !IsOdd();
    };

    // Check if element is invertible
    inline bool IsInvertible() const {
	return IsOdd();
    }

    // Compute inverse of element. Fails if element is not invertible.
    Z2k<K> Invert() const {
	assert (IsInvertible());

	// See https://marc-b-reynolds.github.io/math/2017/09/18/ModInverse.html
	const auto x = *this;
	auto z = x * Z2k<K>::Three;
	z.limbs[0] ^= 2;

	z = z * (Z2k<K>::Two - x * z);  // 5 bits
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

    // Division
    friend Z2k<K> operator/(const Z2k<K> &x, const Z2k<K> &y) {
    	return x * y.Invert();
    };

    // Misc functionality

    // Serialization. buf is assumed to point to at least SizeInBytes() bytes of
    // memory.
    void Pack(unsigned char *buf) const {
	memcpy(buf, this->limbs, SizeInBytes());
    };

    // To string
    std::string ToString() const;

    // Printing
    template<int L>
    friend std::ostream& operator<<(std::ostream &os, const Z2k<L> &x);

#ifdef TESTING

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
#ifdef ZTK_GCC_UINT128
    __uint128_t _x, _r;
    LOAD_U128(_x, x);
    LOAD_U128(_r, y);
    _r += _x;
    STORE_U128(r, _r);
#else
    asm ("movq	%3, %1 \n\t"
    	 "movq	%2, %0 \n\t"
    	 "addq	%5, %1 \n\t"
    	 "adcq	%4, %0"
    	 : "+r" (r[1]), "+r" (r[0])
    	 : "r" (x[1]), "r" (x[0]), "r" (y[1]), "r" (y[0]) : "cc"
    	);
#endif
}

template<>
inline void op_inc<1>(limb_t r[1], const limb_t x[1]) {
    r[0] += x[0];
}

template<>
inline void op_inc<2>(limb_t r[2], const limb_t x[2]) {
#ifdef ZTK_GCC_UINT128
    __uint128_t _x, _r;
    LOAD_U128(_x, x);
    LOAD_U128(_r, r);
    _r += _x;
    STORE_U128(r, _r);
#else
    asm ("addq %3, %1 \n\t"
	 "adcq %2, %0 \n\t"
	 : "+r" (r[1]), "+r" (r[0])
	 : "r" (x[1]), "r" (x[0]) : "cc"
	);
#endif
}

template<>
inline void op_sub<1>(limb_t r[1], const limb_t x[1], const limb_t y[1]) {
    r[0] = x[0] - y[0];
}

template<>
inline void op_sub<2>(limb_t r[2], const limb_t x[2], const limb_t y[2]) {
#ifdef ZTK_GCC_UINT128
    __uint128_t _x, _r;
    LOAD_U128(_x, x);
    LOAD_U128(_r, y);
    _r -= _x;
    STORE_U128(r, _r);
#else
    asm ("movq %3, %1 \n\t"
	 "movq %2, %0 \n\t"
	 "subq %5, %1 \n\t"
	 "sbbq %4, %0"
	 : "+r" (r[1]), "+r" (r[0])
	 : "r" (x[1]), "r" (x[0]), "r" (y[1]), "r" (y[0]) : "cc"
	);
#endif
}

template<>
inline void op_dec<1>(limb_t r[1], const limb_t x[1]) {
    r[0] -= x[0];
}

template<>
inline void op_mul<1>(limb_t r[1], const limb_t x[1], const limb_t y[1]) {
    r[0] = x[0] * y[0];
}

template<>
inline void op_mul<2>(limb_t r[2], const limb_t x[2], const limb_t y[2]) {
#ifdef ZTK_GCC_UINT128
    __uint128_t _x, _r;
    LOAD_U128(_x, x);
    LOAD_U128(_r, y);
    _r *= _x;
    STORE_U128(r, _r);
#else
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
std::ostream &operator<<(std::ostream &os, const Z2k<K> &x) {
    os << x.ToString();
    return os;
}

} // ztk

#endif // _ZTK_HPP
