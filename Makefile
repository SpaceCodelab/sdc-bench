# sdc-bench
#
# -ffp-contract=off forbids fused multiply-add, so every product and sum is rounded
# separately; together with the deterministic generators in src/common.hpp this makes the
# records bit-identical on any x86-64 Linux machine with a C++17 compiler.

CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++17 -ffp-contract=off -fno-fast-math -pthread \
            -Wall -Wno-misleading-indentation -Wno-unused-variable

BIN := bin/campaign bin/verify bin/bench_timing

.PHONY: all quick full test clean

all: $(BIN)

bin/campaign: src/campaign.cpp src/common.hpp | bin
	$(CXX) $(CXXFLAGS) -o $@ src/campaign.cpp

bin/verify: src/verify.cpp src/common.hpp | bin
	$(CXX) $(CXXFLAGS) -o $@ src/verify.cpp

bin/bench_timing: src/bench_timing.cpp | bin
	$(CXX) $(CXXFLAGS) -o $@ src/bench_timing.cpp

bin:
	mkdir -p bin

# Reduced campaign, about a minute.
quick: all
	scripts/run_all.sh quick

# Full campaign as reported in the paper.
full: all
	scripts/run_all.sh full

# Determinism regression test: the quick campaign must reproduce the reference checksums.
test: all
	tests/check_quick.sh

clean:
	rm -rf bin
