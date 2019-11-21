CXX      = g++
CXXFLAGS = -Wall -Wextra -Werror -march=native -fpie -O3

TEST_LDFLAGS = -lsodium -lgmp

tests: test/test-main.o ztk.hpp
	$(CXX) $(CXXFLAGS) -DTESTING test/test-main.o test/tests.cpp -o run_test $(TEST_LDFLAGS)
	@echo -e "\ntesting with GCC __uint128_t extension"
	@./run_test

benchmark: test/test-main.o ztk.hpp
	$(CXX) $(CXXFLAGS) -DTESTING test/test-main.o test/bench.cpp -o run_bench $(TEST_LDFLAGS)
	@./run_bench

test/test-main.o: test/test-main.cpp
	$(CXX) $(CXXFLAGS) test/test-main.cpp -c -o test/test-main.o

# test/bench-main.o: test/bench-main.cpp
# 	$(CXX) $(CXXFLAGS) test/bench-main.cpp -c -o test/bench-main.o

clean:
	rm -f test/test-main.o run_test run_bench

.PHONE: tests benchmark
