CXX=g++

CXXFLAGS=-fopenmp -O3 -Wall -Wextra

# --- Uncomment the following line to enable debugging ---
# CXXFLAGS=-fopenmp -O3 -Wall -Wextra -g

# --- Uncomment the following line to enable sanitizers ---
# CXXFLAGS=-fopenmp -O3 -Wall -Wextra -g -fsanitize=address -fsanitize=undefined -fsanitize=leak -fsanitize=pointer-subtract -fsanitize=pointer-compare -fsanitize=pointer-overflow

all: pulp libpulp

pulp.o:
	$(CXX) $(CXXFLAGS) -c pulp.cpp

libpulp: pulp.o
	ar rvs libpulp.a pulp.o

pulp: pulp.o
	$(CXX) $(CXXFLAGS) -o pulp pulp_main.cpp pulp.o

clean:
	rm -f pulp *.o *.a
