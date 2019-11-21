CXX      = g++
CXXFLAGS = -Wall -Wextra -Werror -march=native -fpie

TEST_LDFLAGS = -lsodium -lgmp

tests: test/test-main.o ztk.hpp
	$(CXX) $(CXXFLAGS) -DTESTING test/test-main.o test/tests.cpp -o test_uint128 $(TEST_LDFLAGS)
	@echo -e "\ntesting with GCC __uint128_t extension"
	@./test_uint128

benchmark: test/bench-main.o ztk.hpp
	$(CXX) $(CXXFLAGS) -DTESTING -O3 test/bench-main.o test/bench.cpp -o bench $(TEST_LDFLAGS)
	@./bench

test/test-main.o: test/test-main.cpp
	$(CXX) $(CXXFLAGS) test/test-main.cpp -c -o test/test-main.o

test/bench-main.o: test/bench-main.cpp
	$(CXX) $(CXXFLAGS) test/bench-main.cpp -c -o test/bench-main.o

clean:
	rm -f test_uint128 test_asm test/test-main.o

.PHONE: tests benchmark
