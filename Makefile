CXX      = g++
CXXFLAGS = -Wall -Wextra -Werror -march=native -fpie -O3

TEST_LDFLAGS = -lsodium -lgmp

ifeq ($(GCC_UINT128), 1)
	CXXFLAGS += -DZTK_GCC_UINT128
endif

tests: ztk.hpp
	$(CXX) $(CXXFLAGS) -DTESTING test-main.o test/tests.cpp -o test_asm $(TEST_LDFLAGS)
	$(CXX) $(CXXFLAGS) -DTESTING -DZTK_GCC_UINT128 test-main.o test/tests.cpp -o test_uint128 $(TEST_LDFLAGS)
	@echo -e "\ntesting with GCC __uint128_t extension"
	@./test_uint128
	@echo "Testing with inline asm"
	@./test_asm

test/test-main.o: test/test-main.cpp
	$(CXX) $(CXXFLAGS) test/test-main.cpp -c

clean:
	rm -f test_uint128 test_asm

.PHONE: tests benchmark
