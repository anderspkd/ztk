#define CATCH_CONFIG_MAIN

#include "catch.hpp"

#include <cstdint>
#include <sodium.h>
#include <gmp.h>

#include <iostream>
#include <sstream>

#include "../ztk.hpp"


using ztk::Z2k;
using ztk::GR;
using ztk::GR4;
using std::string;

template<typename T>
void randomize(T &x) {
    if (sodium_init() < 0)
	assert(0);
    randombytes_buf((T*)&x, sizeof(T));
}

template<size_t K>
Z2k<K> get_rand() {
    auto size = Z2k<K>::SizeInLimbs();
    ztk::limb_t e[size];
    for (size_t i = 0; i < size; i++) {
	randomize<ztk::limb_t>(e[i]);
    }
    return Z2k<K>{e};
}

template<size_t K, size_t D>
GR<K, D> get_rand_gr() {
    std::array<Z2k<K>, D> x;
    for (size_t i = 0; i < D; i++)
	x[i] = get_rand<K>();
    return GR<K, D>{x};
}

TEST_CASE("parameters") {
    Z2k<64> z64;
    Z2k<62> z62;
    Z2k<66> z66;
    Z2k<128> z128;

    // All are 0
    REQUIRE(z64 == Z2k<64>::Zero);
    REQUIRE(z62 == Z2k<62>::Zero);
    REQUIRE(z66 == Z2k<66>::Zero);
    REQUIRE(z128 == Z2k<128>::Zero);

    // check bit/byte/limb sizes
    REQUIRE(z64.SizeInBits() == 64);
    REQUIRE(z62.SizeInBits() == 62);
    REQUIRE(z66.SizeInBits() == 66);
    REQUIRE(z128.SizeInBits() == 128);

    REQUIRE(z64.SizeInBytes() == 8);
    REQUIRE(z62.SizeInBytes() == 8);
    REQUIRE(z66.SizeInBytes() == 9);
    REQUIRE(z128.SizeInBytes() == 16);

    REQUIRE(z64.SizeInLimbs() == 1);
    REQUIRE(z62.SizeInLimbs() == 1);
    REQUIRE(z66.SizeInLimbs() == 2);
    REQUIRE(z128.SizeInLimbs() == 2);

    // check masking
    REQUIRE(z64.NeedsMasking() == false);
    REQUIRE(z62.NeedsMasking() == true);
    REQUIRE(z66.NeedsMasking() == true);
    REQUIRE(z128.NeedsMasking() == false);
}

TEST_CASE("equality") {
    Z2k<64> x0 {12345};
    Z2k<64> x1 {12345};
    Z2k<64> x2 {123};

    REQUIRE(x0 == x1);
    REQUIRE(x0 != x2);
}

TEST_CASE("unary minus") {
    SECTION("64") {
	Z2k<64> x {1234567};
	REQUIRE(-x == -(1234567));
    }

    SECTION("2 limbs") {
	Z2k<96> x = get_rand<96>();
	uint64_t t[2];
	mpn_neg(t, x.GetLimbs(), 2);

	auto r = -x;

	REQUIRE(r.GetLimbs()[0] == t[0]);
	REQUIRE(r.GetLimbs()[1] == (t[1] & Z2k<96>::GetMask()));
    }

    SECTION("2 limbs det") {
	Z2k<96> x (13);
	uint64_t t[2];
	mpn_neg(t, x.GetLimbs(), 2);

	auto r = -x;

	REQUIRE(r.GetLimbs()[0] == t[0]);
	REQUIRE(r.GetLimbs()[1] == (t[1] & Z2k<96>::GetMask()));
    }
}

TEST_CASE("constructors") {
    // 0 constructor tested above

    // test with explicit mask
    SECTION("explicit mask") {
	// Even if the modulus allows for a larger value, the mask argument lets
	// us restrict the element arbitrarily.
	Z2k<64> x {(ztk::limb_t)255, (ztk::limb_t)0x0F};
	Z2k<64> t {(ztk::limb_t)15};
	REQUIRE(x == t);
    }

    SECTION("from limb_t") {
	Z2k<64> x64 {(ztk::limb_t)255};
	REQUIRE(x64.GetLimbs()[0] == 255);
	Z2k<8> x8 {(ztk::limb_t)256};
	REQUIRE(x8.GetLimbs()[0] == 0);
    }

    SECTION("from limb_t[]") {
	ztk::limb_t x[2] = {1234, 555};
	Z2k<128> y {x};
	REQUIRE(y.GetLimbs()[0] == 1234);
	REQUIRE(y.GetLimbs()[1] == 555);
    }

    SECTION("copy") {
	Z2k<64> x {555};
	Z2k<64> y {x};
	REQUIRE(y.GetLimbs()[0] == 555);
	REQUIRE(x == y);

	Z2k<8>  w {y};
	REQUIRE(w.GetLimbs()[0] == (555 & 255));

	// for a larger modulus, the new value is always the same as the old
	// one.
	for (size_t i = 0; i < 100; i++) {
	    Z2k<100> z = get_rand<64>();
	    Z2k<128> zz {z};
	    REQUIRE(z.GetLimbs()[0] == zz.GetLimbs()[0]);
	    REQUIRE(z.GetLimbs()[1] == zz.GetLimbs()[1]);
	}
    }
}

TEST_CASE("printing") {
    Z2k<8> x {123};
    std::string t = "{123}";
    std::stringstream ss;
    ss << x;
    REQUIRE(t == ss.str());
}

TEST_CASE("addition") {
    SECTION("32") {
	Z2k<32> x = get_rand<32>();
	Z2k<32> y = get_rand<32>();
	auto z = x + y;
	uint32_t x_ = (uint32_t)x.GetLimbs()[0];
	uint32_t y_ = (uint32_t)y.GetLimbs()[0];
	auto z_ = x_ + y_;

	REQUIRE(z.GetLimbs()[0] == z_);
    }

    SECTION("64") {
	Z2k<64> x = get_rand<64>();
	Z2k<64> y = get_rand<64>();
	auto z = x + y;
	uint64_t x_ = (uint64_t)x.GetLimbs()[0];
	uint64_t y_ = (uint64_t)y.GetLimbs()[0];
	auto z_ = x_ + y_;

	REQUIRE(z.GetLimbs()[0] == z_);
    }
}

TEST_CASE("increment") {
    SECTION("64") {
	Z2k<64> x = get_rand<64>();
	Z2k<64> y = get_rand<64>();
	uint64_t x_ = (uint64_t)x.GetLimbs()[0];
	uint64_t y_ = (uint64_t)y.GetLimbs()[0];
	x += y;
	x_ += y_;

	REQUIRE(x.GetLimbs()[0] == x_);
	REQUIRE(y.GetLimbs()[0] == y_);
    }

    SECTION("128") {
	Z2k<128> x = get_rand<128>();
	Z2k<128> y = get_rand<128>();

	uint64_t y0 = y.GetLimbs()[0];
	uint64_t y1 = y.GetLimbs()[1];

	uint64_t t[2];
	mpn_add_n(t, x.GetLimbs(), y.GetLimbs(), 2);

	x += y;

	REQUIRE(x.GetLimbs()[0] == t[0]);
	REQUIRE(x.GetLimbs()[1] == t[1]);
	REQUIRE(y0 == y.GetLimbs()[0]);
	REQUIRE(y1 == y.GetLimbs()[1]);
    }
}

TEST_CASE("addition gr") {
    GR<64, 4> x = get_rand_gr<64, 4>();
    auto xcoeff = x.GetCoeff();
    GR<64, 4> y = get_rand_gr<64, 4>();
    auto ycoeff = y.GetCoeff();
    GR<64, 4> r {x + y};
    auto rcoeff = r.GetCoeff();

    REQUIRE(rcoeff[0] == xcoeff[0] + ycoeff[0]);
    REQUIRE(rcoeff[1] == xcoeff[1] + ycoeff[1]);
    REQUIRE(rcoeff[2] == xcoeff[2] + ycoeff[2]);
    REQUIRE(rcoeff[3] == xcoeff[3] + ycoeff[3]);
}

TEST_CASE("assignment gr") {
    GR<123, 5> x = get_rand_gr<123, 5>();
    GR<123, 5> y;

    y = x;

    REQUIRE(x == y);
}

// TEST_CASE("mul gr") {
//     GR4<64> a {{1, 2, 4, 8}};
//     GR4<64> b {{1, 1, 1, 1}};

//     auto v = a * b;

//     // blegh
//     Z2k<64> t0 (13);
//     Z2k<64> t1 (23);
//     GR4<64> c {{-t0, -t1, -t0, 7}};

//     REQUIRE(v == c);
// }
