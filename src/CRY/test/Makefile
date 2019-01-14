#
# make       = run short test
# make all   = run all tests (including root)
# make runRoot  = compile and run specified test
#
# Doug Wright

 include ../Makefile.common  #....define compiler and any other common settings
-include ../Makefile.local   #....local definitions for this host (if any)

EXE = testMain testLong testOut testRoot
CXXFLAGS += -I../src

default: runMain
all: runMain  runOut runRoot runLong

#....rule to make test executables
test%: test%.cc ../lib/libCRY.a
	$(CXX) $(CXXFLAGS) -o $@ $@.cc -L../lib -lCRY
 
#....testRoot is different
testRoot: testRoot.cc ../lib/libCRY.a
	$(CXX) $(CXXFLAGS) -I$(ROOTSYS)/include -o $@ $@.cc -L../lib -lCRY $(shell root-config --libs)

clean:
	rm -f *.o $(EXE)

run%: test%
	@echo
	@echo "Run $+"
	@[ $+ = "testLong" ] && echo "(takes 1 min on aztec)" || true
	@[ $+ = "testRoot" ] && echo "(makes test.root; use root compare_root.cc to see differences)" || true
	@[ $+ = "testOut"  ] && echo "(writes cry particles to shower.out)" || true
	./$+ setup.file > $+.out
	-diff -s $+.out $+.ref
	@[ $+ = "testRoot" ] && (root -l -q -b compare_root.cc; echo check compare_root.pdf) || true

#....make sure these executables are not deleted as intermediate files
.PRECIOUS: $(EXE)