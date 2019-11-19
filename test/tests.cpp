#include "catch.hpp"

#include "../ztk.hpp"

// helpers
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


TEST_CASE("constructors") {
    SECTION("0") {
	Z2k<64> x64;
	Z2k<32> x32;
	Z2k<123> x123;
	Z2k<128> x128;

	REQUIRE(x64.GetLimbs()[0] == 0);
	REQUIRE(x32.GetLimbs()[0] == 0);
	REQUIRE(x123.GetLimbs()[0] == 0);
	REQUIRE(x123.GetLimbs()[1] == 0);
	REQUIRE(x128.GetLimbs()[0] == 0);
	REQUIRE(x128.GetLimbs()[1] == 0);
    }

    SECTION("from limb_t* and mask") {
	limb_t e[1] = {123456};
	Z2k<64> x {e, Z2k<64>::GetMask()};
	REQUIRE(x.GetLimbs()[0] == 123456);

	limb_t r[2] = {123456, 345671};
	Z2k<65> y {r, Z2k<65>::GetMask()};
	REQUIRE(y.GetLimbs()[0] == 123456);
	REQUIRE(y.GetLimbs()[1] == 1);
    }

    SECTION("from limb_t and mask") {
	limb_t e = 1231231;
	Z2k<64> x {e, Z2k<64>::GetMask()};
	REQUIRE(x.GetLimbs()[0] == e);
    }

    SECTION("from limb_t") {
	limb_t e = 0;
	randomize<limb_t>(&e, 1);
	REQUIRE(e != 0);  // sanity check
	Z2k<64> x {e};
	REQUIRE(x.GetLimbs()[0] == e);

	Z2k<32> y {e};
	// Shouldn't be the same as limb_t is 64 bits and y is only 32
	REQUIRE(y.NeedsMasking());
	REQUIRE(y.GetLimbs()[0] != e);
	REQUIRE(y.GetLimbs()[0] == (e & 0xFFFFFFFF));
	REQUIRE(y.GetLimbs()[0] != (e & 0xFFFFFFFF00000000));
	REQUIRE((y.GetLimbs()[0] & 0xFFFFFFFF00000000) == 0);
    }

    SECTION("from limb_t*") {
	limb_t e[1];
	randomize<limb_t>(e, 1);
	REQUIRE(e[0] != 0);
	Z2k<64> x {e};

	REQUIRE(x.GetLimbs()[0] == e[0]);

	limb_t f[2];
	randomize<limb_t>(f, 2);
	REQUIRE(f[0] != 0);
	REQUIRE(f[1] != 0);

	Z2k<128> y {f};
	REQUIRE(y.GetLimbs()[0] == f[0]);
	REQUIRE(y.GetLimbs()[1] == f[1]);

	Z2k<96> w {f};
	REQUIRE(w.GetLimbs()[0] == f[0]);
	REQUIRE(w.GetLimbs()[1] == (f[1] & 0xFFFFFFFF));
	REQUIRE(w.GetLimbs()[1] != (f[1] & 0xFFFFFFFF00000000));
	REQUIRE((w.GetLimbs()[1] & 0xFFFFFFFF00000000) == 0);

    }

    SECTION("copy") {
	limb_t e[2];
	randomize<limb_t>(e, 2);

	Z2k<128> x {e};
	Z2k<128> y {x};

	REQUIRE(x == y);
	REQUIRE(x.GetLimbs()[0] == y.GetLimbs()[0]);
	REQUIRE(x.GetLimbs()[1] == y.GetLimbs()[1]);
    }

    SECTION("copy with diff modulus") {
	limb_t e[2];
	randomize<limb_t>(e, 2);

	Z2k<128> x {e};
	Z2k<96> y {x};

	REQUIRE(x != y);
	REQUIRE(x.GetLimbs()[0] == y.GetLimbs()[0]);
	REQUIRE(x.GetLimbs()[1] != y.GetLimbs()[1]);
	REQUIRE((x.GetLimbs()[1] & 0xFFFFFFFF) == y.GetLimbs()[1]);
    }

    SECTION("from a raw buffer") {
	unsigned char buf[sizeof(limb_t) * 2];
	randomize<unsigned char>(buf, sizeof(limb_t) * 2);

	Z2k<128> x {buf};

	auto limb_buf = (limb_t *)buf;

	REQUIRE(x.GetLimbs()[0] == limb_buf[0]);
	REQUIRE(x.GetLimbs()[1] == limb_buf[1]);
	REQUIRE(x.GetLimbs()[0] != 0);
	REQUIRE(x.GetLimbs()[1] != 0);
    }

    // construction from a string is not currently supported
}

TEST_CASE("assignment") {
    limb_t e[2];
    randomize<limb_t>(e, 2);
    limb_t f[2];
    randomize<limb_t>(f, 2);

    Z2k<128> x {e};
    Z2k<128> y {f};

    REQUIRE(x != y);

    y = x;

    REQUIRE(x == y);
}
