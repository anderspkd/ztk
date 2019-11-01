#define CATCH_CONFIG_MAIN

#include "catch.hpp"

#include <cstdint>
#include <sodium.h>
#include <gmp.h>

#include <iostream>
#include <sstream>

#include "../ztk.hpp"


using ztk::Z2k;
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
