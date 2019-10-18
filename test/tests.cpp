#define CATCH_CONFIG_MAIN

#include "catch.hpp"

#include <cstdint>
#include <sodium.h>
#include <gmp.h>

#include "../ztk.hpp"


typedef uint32_t u32;
typedef uint64_t u64;

using ztk::Z2k;
using std::string;


template<typename T>
static void randomize(T &x) {
    if (sodium_init() < 0)
	assert(0);
    randombytes_buf((T*)&x, sizeof(T));
}


TEST_CASE("Z2k construct 0", "[constructor]") {
    SECTION("32") {
	Z2k<32> x {(int)0};
	CHECK(x.IsZero());
	CHECK(x.SizeInLimbs() == 1);
    }
    SECTION("35") {
	Z2k<35> x {(int)0};
	CHECK(x.IsZero());
	CHECK(x.SizeInLimbs() == 1);
    }
    SECTION("64") {
	Z2k<64> x {(int)0};
	CHECK(x.IsZero());
	CHECK(x.SizeInLimbs() == 1);
    }
    SECTION("67") {
	Z2k<67> x {(int)0};
	CHECK(x.IsZero());
	CHECK(x.SizeInLimbs() == 2);
    }
    SECTION("100") {
	Z2k<100> x {(int)0};
	CHECK(x.IsZero());
	CHECK(x.SizeInLimbs() == 2);
    }
    SECTION("128") {
	Z2k<128> x {(int)0};
	CHECK(x.IsZero());
	CHECK(x.SizeInLimbs() == 2);
    }
}

TEST_CASE("Z2k addition", "[add]") {
    SECTION("32") {
	u32 x, y, z;
	randomize<u32>(x);
	randomize<u32>(y);
	z = x + y;

	Z2k<32> x_ {x};
	Z2k<32> y_ {y};
	Z2k<32> z_ {x_ + y_};

	CHECK(z_.NeedsMasking() == true);
	CHECK(z_.SizeInLimbs() == 1);
	CHECK((u32)z_.GetLimbs()[0] == z);
    }

    SECTION("35") {
	u64 x, y, z;
	randomize<u64>(x);
	randomize<u64>(y);
	z = x + y;

	size_t mask = 0x7ffffffff;

	Z2k<35> x_ {x & mask};
	Z2k<35> y_ {y & mask};
	Z2k<35> z_ {x_ + y_};

	z &= mask;

	CHECK(z_.NeedsMasking() == true);
	CHECK(z_.SizeInLimbs() == 1);
	CHECK((u64)z_.GetLimbs()[0] == z);
    }

    SECTION("64") {
	u64 x, y, z;
	randomize<u64>(x);
	randomize<u64>(y);
	z = x + y;

	Z2k<64> x_ {x};
	Z2k<64> y_ {y};
	Z2k<64> z_ {x_ + y_};

	CHECK(z_.NeedsMasking() == false);
	CHECK(z_.SizeInLimbs() == 1);
	CHECK((u64)z_.GetLimbs()[0] == z);
    }

    SECTION("128") {
    	u64 x[2], y[2], t[2];
	for (size_t i = 0; i < 2; i++) {
	    randomize<u64>(x[i]);
	    randomize<u64>(y[i]);
	}

	mpn_add_n(t, x, y, 2);

    	Z2k<128> x_ {(ztk::limb_t *)x};
    	Z2k<128> y_ {(ztk::limb_t *)y};

	CHECK(x_.GetLimbs()[0] == x[0]);
	CHECK(x_.GetLimbs()[1] == x[1]);
	CHECK(y_.GetLimbs()[0] == y[0]);
	CHECK(y_.GetLimbs()[1] == y[1]);

    	Z2k<128> z_ {x_ + y_};
	Z2k<128> t_ {t};

    	CHECK(z_.NeedsMasking() == false);
    	CHECK(z_.SizeInLimbs() == 2);
	CHECK(t_.GetLimbs()[0] == t[0]);
	CHECK(t_.GetLimbs()[1] == t[1]);
    	CHECK(z_.GetLimbs()[0] == t[0]);
    	CHECK(z_.GetLimbs()[1] == t[1]);
    }

    SECTION("101") {
    	u64 x[2], y[2], t[2];
	for (size_t i = 0; i < 2; i++) {
	    randomize<u64>(x[i]);
	    randomize<u64>(y[i]);
	}

	size_t mask = 0x1fffffffff;

	mpn_add_n(t, x, y, 2);

    	Z2k<101> x_ {(ztk::limb_t *)x};
    	Z2k<101> y_ {(ztk::limb_t *)y};

	CHECK(x_.GetLimbs()[0] == x[0]);
	CHECK(x_.GetLimbs()[1] == (x[1] & mask));
	CHECK(y_.GetLimbs()[0] == y[0]);
	CHECK(y_.GetLimbs()[1] == (y[1] & mask));

    	Z2k<101> z_ {x_ + y_};

	t[1] &= mask;

	Z2k<101> t_ {t};

    	CHECK(z_.NeedsMasking() == true);
    	CHECK(z_.SizeInLimbs() == 2);
	CHECK(t_.GetLimbs()[0] == t[0]);
	CHECK(t_.GetLimbs()[1] == t[1]);
    	CHECK(z_.GetLimbs()[0] == t[0]);
    	CHECK(z_.GetLimbs()[1] == t[1]);
    }
}

TEST_CASE("Z2k multiplication", "[mul]") {
    SECTION("64") {
	u64 x, y, z;
	randomize<u64>(x);
	randomize<u64>(y);

	z = x * y;

	Z2k<64> x_ {x};
	Z2k<64> y_ {y};
	Z2k<64> z_ {x_ * y_};

	CHECK(z_.GetLimbs()[0] == z);
    }

    // SECTION("128") {
    // 	u64 x[2], y[2], z[4];
    // 	for (size_t i = 0; i < 2; i++) {
    // 	    randomize<u64>(x[i]);
    // 	    randomize<u64>(y[i]);
    // 	}

    // 	mpn_mul_n(z, x, y, 2);

    // 	Z2k<128> x_ {x};
    // 	Z2k<128> y_ {y};
    // 	Z2k<128> z_ {x_ * y_};

    // 	CHECK(z_.GetLimbs()[0] == z[0]);
    // 	CHECK(z_.GetLimbs()[1] == z[1]);
    // }
}
