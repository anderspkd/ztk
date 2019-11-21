#include "test-main.cpp"

#include "../ztk.hpp"

#include <sodium.h>
#include <gmp.h>

using namespace ztk;

template<typename T>
void randomize(T* x, size_t n) {
    if (sodium_init() < 0)
	assert(0);
    randombytes_buf(x, sizeof(T) * n);
}

template<size_t K>
Z2k<K> random_Z2k() {
    auto size = Z2k<K>::SizeInLimbs();
    limb_t e[size];
    randomize<limb_t>(e, size);
    return Z2k<K>{e};
}

template<size_t K>
Z2k<K> bench_mul(Z2k<K> x, Z2k<K> y) {
    return x * y;
}

template<size_t L, size_t K>
limb_t bench_mul_gmp(limb_t r[2*L], limb_t *x, limb_t *y) {
    mpn_mul_n(r, x, y, L);
    return r[0] + r[1];
}

TEST_CASE("multiplication") {
    auto x = random_Z2k<64>();
    auto y = random_Z2k<64>();

    BENCHMARK("ztk") {
	return bench_mul<64>(x, y);
    };

    auto xlimbs = (limb_t *)x.GetLimbs();
    auto ylimbs = (limb_t *)y.GetLimbs();
    limb_t r[2*Z2k<64>::SizeInLimbs()];

    BENCHMARK("gmp") {
	return bench_mul_gmp<Z2k<64>::SizeInLimbs(), 64>(r, xlimbs, ylimbs);
    };
}
