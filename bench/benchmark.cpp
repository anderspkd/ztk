#include <sodium.h>
#include <gmp.h>
#include <cstdint>
#include <cassert>
#include <iostream>

#include "bench.h"
#include "../ztk.hpp"

typedef uint64_t u64;
typedef __uint128_t u128;

using ztk::Z2k;

template<typename T>
static void randomize(T &x) {
    randombytes_buf((T*)&x, sizeof(T));
}

template<size_t K>
static inline void add_ztk(volatile Z2k<K> &z, Z2k<K> &x, Z2k<K> &y) {
    z = x + y;
    assert(!z.IsOne());
}

int main(int argc, char **argv) {

    (void)argc;
    (void)argv;

    u64 x[2], y[2], z[4];

    u128 a, b;
    volatile u128 c = 0;

    Z2k<128> x_, y_;
    volatile Z2k<128> z_;

    if (sodium_init() < 0)
	assert(0);

    BENCH_BEGIN("add_128_gmp") {
    	for (size_t i = 0; i < 2; i++) {
    	    randomize<u64>(x[i]);
    	    randomize<u64>(y[i]);
    	}
    	BENCH_ADD(mpn_add_n(z, x, y, 2));
    } BENCH_END;

    BENCH_BEGIN("add_128_gcc") {
    	randomize<u128>(a);
    	randomize<u128>(b);
    	BENCH_ADD(c = a + b);
    } BENCH_END;

    BENCH_BEGIN("add_128_ztk") {
	for (size_t i = 0; i < 2; i++) {
	    randomize<u64>(x[i]);
	    randomize<u64>(y[i]);
	}
	x_ = x;
	y_ = y;
	BENCH_ADD(volatile Z2k<128> z_ {x_ + y_});
    } BENCH_END;

    BENCH_BEGIN("mul_128_gmp") {
    	for (size_t i = 0; i < 2; i++) {
    	    randomize<u64>(x[i]);
    	    randomize<u64>(y[i]);
    	}
    	BENCH_ADD(mpn_mul_n(z, x, y, 2));
    } BENCH_END;

    BENCH_BEGIN("mul_128_gcc") {
    	randomize<u128>(a);
    	randomize<u128>(b);
    	BENCH_ADD(c = a * b);
    } BENCH_END;

    BENCH_BEGIN("mul_128_ztk") {
    	for (size_t i = 0; i < 2; i++) {
    	    randomize<u64>(x[i]);
    	    randomize<u64>(y[i]);
    	}
    	ztk::Z2k<128> x_ {x};
    	ztk::Z2k<128> y_ {y};
    	BENCH_ADD(volatile ztk::Z2k<128> r {x_ * y_});
    } BENCH_END;
}
