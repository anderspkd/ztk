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

template<size_t K, size_t D>
GR<K, D> random_gr() {
    gr_coeff<K, D> coeff;
    for (size_t i = 0; i < D; i++)
	coeff[i] = random_Z2k<K>();
    return GR<K, D>{coeff};
}

template<size_t K>
GR4<K> random_gr4() {
    return static_cast<GR4<K>>(random_gr<K, 4>());
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

template<size_t K>
void test_addition() {
    auto a = random_Z2k<K>();
    auto b = random_Z2k<K>();

    auto c = a + b;

    auto ln = Z2k<K>::SizeInLimbs();
    limb_t t[ln];
    mpn_add_n(t, a.GetLimbs(), b.GetLimbs(), ln);

    for (size_t i = 0; i < ln - 1; i++)
	REQUIRE(c.GetLimbs()[i] == t[i]);

    REQUIRE(c.GetLimbs()[ln-1] == (t[ln-1] & Z2k<K>::GetMask()));

    a += b;

    REQUIRE(a == c);
}

TEST_CASE("addition") {
    SECTION("32") {
	test_addition<32>();
    }
    SECTION("64") {
	test_addition<64>();
    }
    SECTION("123") {
	test_addition<123>();
    }
    SECTION("128") {
	test_addition<128>();
    }
    SECTION("42") {
	test_addition<42>();
    }
    SECTION("96") {
	test_addition<96>();
    }
}

template<size_t K>
void test_subtraction() {
    auto a = random_Z2k<K>();
    auto b = random_Z2k<K>();

    auto c0 = a - b;
    auto c1 = b - a;

    auto ln = Z2k<K>::SizeInLimbs();
    limb_t t0[ln];
    limb_t t1[ln];

    mpn_sub_n(t0, a.GetLimbs(), b.GetLimbs(), ln);
    mpn_sub_n(t1, b.GetLimbs(), a.GetLimbs(), ln);

    for (size_t i = 0; i < ln - 1; i++) {
	REQUIRE(c0.GetLimbs()[i] == t0[i]);
	REQUIRE(c1.GetLimbs()[i] == t1[i]);
    }

    REQUIRE(c0.GetLimbs()[ln-1] == (t0[ln-1] & Z2k<K>::GetMask()));
    REQUIRE(c1.GetLimbs()[ln-1] == (t1[ln-1] & Z2k<K>::GetMask()));

    auto d = a;

    a -= b;
    REQUIRE(a == c0);

    b -= d;  // b -= a;
    REQUIRE(b == c1);
}

TEST_CASE("subtraction") {
    SECTION("32") {
	test_subtraction<32>();
    }
    SECTION("64") {
	test_subtraction<64>();
    }
    SECTION("123") {
	test_subtraction<123>();
    }
    SECTION("128") {
	test_subtraction<128>();
    }
    SECTION("42") {
	test_subtraction<42>();
    }
    SECTION("96") {
	test_subtraction<96>();
    }
}

template<size_t K>
void test_multiplication() {
    auto a = random_Z2k<K>();
    auto b = random_Z2k<K>();

    auto c = a * b;

    auto ln = Z2k<K>::SizeInLimbs();
    limb_t t[ln * 2];

    mpn_mul_n(t, a.GetLimbs(), b.GetLimbs(), ln);

    for (size_t i = 0; i < ln - 1; i++)
	REQUIRE(c.GetLimbs()[i] == t[i]);

    REQUIRE(c.GetLimbs()[ln-1] == (t[ln-1] & Z2k<K>::GetMask()));
}

TEST_CASE("multiplication") {
    SECTION("32") {
    	test_multiplication<32>();
    }
    SECTION("64") {
	test_multiplication<64>();
    }
    SECTION("123") {
    	test_multiplication<123>();
    }
    SECTION("128") {
    	test_multiplication<128>();
    }
    SECTION("42") {
    	test_multiplication<42>();
    }
    SECTION("96") {
    	test_multiplication<96>();
    }
}

template<size_t K>
void test_division() {
    auto a = random_Z2k<K>();
    while (!a.IsInvertible())
	a = random_Z2k<K>();

    auto b = a;

    REQUIRE(a.IsInvertible());  // sanitiy check

    for (size_t i = 0; i < Z2k<K>::SizeInLimbs(); i++)
	REQUIRE(b.GetLimbs()[i] != 0);

    auto c = a / b;

    REQUIRE(c == Z2k<K>::One);
}

TEST_CASE("division") {
    SECTION("32") {
    	test_division<32>();
    }
    SECTION("64") {
	test_division<64>();
    }
    SECTION("123") {
    	test_division<123>();
    }
    SECTION("128") {
    	test_division<128>();
    }
    SECTION("42") {
    	test_division<42>();
    }
    SECTION("96") {
    	test_division<96>();
    }
}

TEST_CASE("equality") {
    SECTION("1 limb") {
	auto a = random_Z2k<64>();
	auto b = random_Z2k<64>();
	REQUIRE(a != b);
	auto c = a;
	REQUIRE(a == c);
    }
    SECTION("2 limbs") {
	auto a = random_Z2k<123>();
	auto b = random_Z2k<123>();
	REQUIRE(a != b);
	auto c = a;
	REQUIRE(a == c);
    }

    limb_t x = 43;
    Z2k<32> a{x};
    REQUIRE(a.IsOdd());
    REQUIRE(a.IsInvertible());

    limb_t y = 42;
    Z2k<32> b{y};
    REQUIRE(b.IsEven());
}

TEST_CASE("packing/serialization") {
    auto a = random_Z2k<123>();
    unsigned char buf[a.SizeInBytes()];
    a.Pack(buf);
    Z2k<123> b{buf};
    REQUIRE(a == b);
}

////////////////////////////////////////////////////////////////////////////////
// Galois rings tests
////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Constructors") {
    SECTION("0") {
	GR<32, 4> x;
	REQUIRE(x.GetCoeff()[0] == Z2k<32>::Zero);
	REQUIRE(x.GetCoeff()[1] == Z2k<32>::Zero);
	REQUIRE(x.GetCoeff()[2] == Z2k<32>::Zero);
	REQUIRE(x.GetCoeff()[3] == Z2k<32>::Zero);
    }
    SECTION("from coefficients") {
	gr_coeff<32, 4> a = {Z2k<32>::One, Z2k<32>::Two, Z2k<32>::Three, Z2k<32>::Four};
	GR<32, 4> x {a};
	REQUIRE(x.GetCoeff()[0] == Z2k<32>::One);
	REQUIRE(x.GetCoeff()[1] == Z2k<32>::Two);
	REQUIRE(x.GetCoeff()[2] == Z2k<32>::Three);
	REQUIRE(x.GetCoeff()[3] == Z2k<32>::Four);
    }
    SECTION("from random coefficients") {
	gr_coeff<64, 5> a;
	for (size_t i = 0; i < 5; i++)
	    a[i] = random_Z2k<64>();
	GR<64, 5> x {a};
	for (size_t i = 0; i < 5; i++)
	    REQUIRE(a[i] == x.GetCoeff()[i]);
    }
    SECTION("from other") {
	auto a = random_gr<64, 4>();
	GR<64, 4> b {a};

	for (size_t i = 0; i < 4; i++)
	    REQUIRE(a.GetCoeff()[i] == b.GetCoeff()[i]);
    }
    SECTION("from z2k") {
	const auto a = random_Z2k<64>();
	auto b = random_Z2k<64>();
	GR<64, 4> x (a);
	GR<64, 4> y (b);
	REQUIRE(x.GetCoeff()[0] == a);
	REQUIRE(x.GetCoeff()[1] == Z2k<32>::Zero);
	REQUIRE(x.GetCoeff()[2] == Z2k<32>::Zero);
	REQUIRE(x.GetCoeff()[3] == Z2k<32>::Zero);
	REQUIRE(y.GetCoeff()[0] == b);
	REQUIRE(y.GetCoeff()[1] == Z2k<32>::Zero);
	REQUIRE(y.GetCoeff()[2] == Z2k<32>::Zero);
	REQUIRE(y.GetCoeff()[3] == Z2k<32>::Zero);
    }
}

TEST_CASE("gr assignment") {
    auto a = random_gr<64, 4>();
    auto b = a;

    for (size_t i = 0; i < 4; i++)
	REQUIRE(a.GetCoeff()[i] == b.GetCoeff()[i]);
}

template<size_t K, size_t D>
void test_addition_and_subtraction() {
    // doing two for one in this function

    auto a = random_gr<K, D>();
    auto b = random_gr<K, D>();

    auto c = a + b;
    auto d = a - b;

    for (size_t i = 0; i < D; i++) {
	REQUIRE(c.GetCoeff()[i] == (a.GetCoeff()[i] + b.GetCoeff()[i]));
	REQUIRE(d.GetCoeff()[i] == (a.GetCoeff()[i] - b.GetCoeff()[i]));
    }
}

TEST_CASE("gr addition subtraction") {
    SECTION("32, 4") {
	test_addition_and_subtraction<32, 4>();
    }
    SECTION("42, 4") {
	test_addition_and_subtraction<32, 4>();
    }
    SECTION("64, 4") {
	test_addition_and_subtraction<32, 4>();
    }
    SECTION("96, 4") {
	test_addition_and_subtraction<32, 4>();
    }
    SECTION("123, 4") {
	test_addition_and_subtraction<32, 4>();
    }
    SECTION("128, 4") {
	test_addition_and_subtraction<32, 4>();
    }
}

TEST_CASE("invert") {
    GR4<123> a = random_gr4<123>();
    while (!a.IsInvertible())
	a = random_gr4<123>();
    GR4<123> b = a;

    GR4<123> c = a / b;

    REQUIRE(c.GetCoeff()[0] == Z2k<123>::One);
    REQUIRE(c.GetCoeff()[1] == Z2k<123>::Zero);
    REQUIRE(c.GetCoeff()[2] == Z2k<123>::Zero);
    REQUIRE(c.GetCoeff()[3] == Z2k<123>::Zero);
}
