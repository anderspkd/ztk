CXX      = g++
CXXFLAGS = -Wall -Wextra -Werror -march=native -fpie

TEST_LDFLAGS = -lsodium -lgmp

TEST_EXEC = runtest

tests: ztk.hpp
	$(CXX) $(CXXFLAGS) -O3 -DTESTING test/tests.cpp -o $(TEST_EXEC) $(TEST_LDFLAGS)
	./$(TEST_EXEC)

clean:
	rm -f $(TEST_EXEC)

.PHONE: tests
