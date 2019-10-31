#ifndef _ZTK_HPP
#define _ZTK_HPP

#include <cstdint>
#include <string>

namespace ztk {

typedef uint64_t limb_t;

template<bool B>
inline void mask_if(limb_t&, const limb_t);

template<size_t N>
inline void op_add(limb_t[N], const limb_t[N], const limb_t[N]);

template<size_t N>
inline void op_inc(limb_t[N], const limb_t[N]);

template<size_t N>
inline void op_sub(limb_t[N], const limb_t[N], const limb_t[N]);

template<size_t N>
inline void op_dec(limb_t[N], const limb_t[N]);

template<size_t N>
inline void op_mul(limb_t[N], const limb_t[N], const limb_t[N]);

// Z_2^K
template<size_t K>
class Z2k {
public:

    static const Z2k<K> Zero;
    static constexpr Z2k<K> One {1};
    static constexpr Z2k<K> Two {2};
    static constexpr Z2k<K> Three {3};

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
    }

    // Constructor. Default to 0
    constexpr Z2k() {};

    // .. From another Z2k element
    Z2k(const Z2k<K> &other) : Z2k{other.limbs} {};

    // .. From another element with different modulus
    template<size_t L>
    Z2k(const Z2k<L> &other);

    // .. From a buffer
    Z2k(const unsigned char *buf) : Z2k{(const limb_t *)buf} {};

    // .. From a string
    Z2k(const std::string &str, const size_t base = 10);

    // .. From a limb (i.e., uint64_t)
    Z2k(const limb_t x) {
	limbs[0] = x;
	mask_if<NeedsMasking()>(limbs[0], mask);
    };

    // .. From a set of limbs
    Z2k(const limb_t limbs[SizeInLimbs()]) {
	memcpy(this->limbs, limbs, SizeInBytes());
	mask_if<NeedsMasking()>(this->limbs[SizeInLimbs()-1], mask);
    };

    Z2k<K> operator=(const Z2k<K> &x) {
	memcpy(this->limbs, x.limbs, SizeInBytes());
	return *this;
    };

    void Pack(unsigned char *buf) const {
	memcpy(buf, this->limbs, SizeInBytes());
    };

    // arithmetic operators
    friend Z2k<K> operator+(const Z2k<K> &x, const Z2k<K> &y) {
	Z2k<K> r;
	op_add<SizeInLimbs()>(r.limbs, x.limbs, y.limbs);
	mask_if<NeedsMasking()>(r.limbs[SizeInLimbs()-1], mask);
	return r;
    };

    Z2k<K> operator+=(const Z2k<K> &x) {
	op_inc<SizeInLimbs()>(limbs, x.limbs);
	mask_if<NeedsMasking()>(limbs[SizeInLimbs()-1], mask);
	return *this;
    };

    friend Z2k<K> operator-(const Z2k<K> &x, const Z2k<K> &y) {
    	Z2k<K> r;
	op_sub<SizeInLimbs()>(r.limbs, x.limbs, y.limbs);
	mask_if<NeedsMasking()>(r.limbs[SizeInLimbs()-1], mask);
	return r;
    };

    Z2k<K> operator-=(const Z2k<K> &x) {
	op_dec<SizeInLimbs()>(limbs, x.limbs);
	mask_if<NeedsMasking()>(limbs[SizeInLimbs()-1], mask);
	return *this;
    }

    friend Z2k<K> operator*(const Z2k<K> &x, const Z2k<K> &y) {
	Z2k<K> r;
	op_mul<SizeInLimbs()>(r.limbs, x.limbs, y.limbs);
	mask_if<NeedsMasking()>(r.limbs[SizeInLimbs()-1], mask);
	return r;
    };

    bool operator==(const Z2k<K> &x) const {
	bool b = true;
	for (size_t i = 0; i < SizeInLimbs(); i++) {
	    b &= limbs[i] == x.limbs[i];
	}
	return b;
    };

    bool operator!=(const Z2k<K> &x) const {
	return !(*this == x);
    };

    inline bool IsOdd() const {
	return (bool)(limbs[0] & 1);
    };

    inline bool IsEven() const {
	return !IsOdd();
    };

    inline bool IsZero() const {
	return *this == Z2k<K>::Zero;
    };

    inline bool IsInvertible() const {
	// IsOdd() == true implies non-zero.
	return IsOdd();
    }

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

    friend Z2k<K> operator/(const Z2k<K> &x, const Z2k<K> &y) {
    	return x * y.Invert();
    };

    // template<int L>
    // friend std::ostream& operator<<(std::ostream &os, const Z2k<L> &x);

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

// Addition code

template<>
inline void op_add<1>(limb_t r[1], const limb_t x[1], const limb_t y[1]) {
    r[0] = x[0] + y[0];
}

template<>
inline void op_add<2>(limb_t r[2], const limb_t x[2], const limb_t y[2]) {
#ifdef ZTK_GCC_UINT128
    __uint128_t _x = ((__uint128_t)(x[1]) << 64) | x[0];
    __uint128_t _y = ((__uint128_t)(y[1]) << 64) | y[0];
    auto _r = _x + _y;
    r[0] = (limb_t)_r;
    r[1] = (limb_t)(_r >> 64);
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
    asm ("addq %3, %1 \n\t"
	 "adcq %2, %0 \n\t"
	 : "+r" (r[1]), "+r" (r[0])
	 : "r" (x[1]), "r" (x[0]) : "cc"
	);
}

// Subtraction code

template<>
inline void op_sub<1>(limb_t r[1], const limb_t x[1], const limb_t y[1]) {
    r[0] = x[0] - y[0];
}

template<>
inline void op_sub<2>(limb_t r[2], const limb_t x[2], const limb_t y[2]) {
#ifdef ZTK_GCC_UINT128
    __uint128_t _x = ((__uint128_t)(x[1]) << 64) | x[0];
    __uint128_t _y = ((__uint128_t)(y[1]) << 64) | y[0];
    auto _r = _x - _y;
    r[0] = (limb_t)_r;
    r[1] = (limb_t)(_r >> 64);
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
    __uint128_t _x = ((__uint128_t)(x[1]) << 64) | x[0];
    __uint128_t _y = ((__uint128_t)(y[1]) << 64) | y[0];
    auto _r = _x * _y;
    r[0] = (limb_t)_r;
    r[1] = (limb_t)(_r >> 64);
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

} // ztk

#endif // _ZTK_HPP
