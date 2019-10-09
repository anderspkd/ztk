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
inline void op_sub(limb_t[N], const limb_t[N], const limb_t[N]);

template<size_t N>
inline void op_mul(limb_t[N], const limb_t[N], const limb_t[N]);

// Z_2^K
template<size_t K>
class Z2k {
public:

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
    Z2k(const unsigned char *buf);

    // .. From a string
    Z2k(const std::string &str, const size_t base = 10);

    // .. From a limb (i.e., uint64_t)
    Z2k(const limb_t x) {
	limbs[0] = x;
	mask_if<NeedsMasking()>(limbs[0], mask);
    };

    Z2k(const limb_t limbs[SizeInLimbs()]) {
	memcpy(this->limbs, limbs, SizeInBytes());
	mask_if<NeedsMasking()>(this->limbs[SizeInLimbs()-1], mask);
    };

    // arithmetic operators
    friend Z2k<K> operator+(const Z2k<K> &x, const Z2k<K> &y) {
	Z2k<K> r;
	op_add<SizeInLimbs()>(r.limbs, x.limbs, y.limbs);
	mask_if<NeedsMasking()>(r.limbs[SizeInLimbs()-1], mask);
	return r;
    };

    friend Z2k<K> operator-(const Z2k<K> &x, const Z2k<K> &y) {
    	Z2k<K> r;
	op_sub<SizeInLimbs()>(r.limbs, x.limbs, y.limbs);
	mask_if<NeedsMasking()>(r.limbs[SizeInLimbs()-1], mask);
	return r;
    };

    friend Z2k<K> operator*(const Z2k<K> &x, const Z2k<K> &y) {
	Z2k<K> r;
	op_mul<SizeInLimbs()>(r.limbs, x.limbs, y.limbs);
	mask_if<NeedsMasking()>(r.limbs[SizeInLimbs()-1], mask);
	return r;
    };

    bool IsZero() const { return limbs[0] == (limb_t)0; };
    bool IsOne() const  { return limbs[0] == (limb_t)1; };

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

    limb_t limbs[number_of_limbs];
};

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
    asm ("movq	%3, %1 \n\t"						\
    	 "movq	%2, %0 \n\t"						\
	 "addq	%5, %1 \n\t"						\
	 "adcq	%4, %0"							\
	 : "+r" (r[1]), "+r" (r[0])					\
	 : "r" (x[1]), "r" (x[0]), "r" (y[1]), "r" (y[0]) : "cc"	\
	);
}

template<>
inline void op_sub<1>(limb_t r[1], const limb_t x[1], const limb_t y[1]) {
    r[0] = x[0] - y[0];
}

template<>
inline void op_mul<1>(limb_t r[1], const limb_t x[1], const limb_t y[1]) {
    r[0] = x[0] * y[0];
}

} // ztk

#endif // _ZTK_HPP
