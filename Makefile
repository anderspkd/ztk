CXX      = g++
CXXFLAGS = -Wall -Wextra -Werror -march=native -fpie

TEST_LDFLAGS = -lsodium -lgmp

TEST_EXEC  = runtest
BENCH_EXEC = benchmark

tests: ztk.hpp
	$(CXX) $(CXXFLAGS) -O3 -DTESTING test/tests.cpp -o $(TEST_EXEC) $(TEST_LDFLAGS)
	./$(TEST_EXEC)

benchmark: ztk.hpp bench/benchmark.cpp
	$(CXX) -march=native -O3 bench/bench.c bench/benchmark.cpp -o $(BENCH_EXEC) -lgmp -lsodium

clean:
	rm -f $(TEST_EXEC) $(BENCH_EXEC)

.PHONE: tests benchmark
