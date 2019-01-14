
 include ../Makefile.common  #....define compiler and any other common settings
-include ../Makefile.local   #....local definitions for this host (if any)

SRC = $(wildcard *.cc)
OBJ = $(SRC:%.cc=%.o)
DEP = $(SRC:%.cc=.%.d)  
LIB = ../lib/libCRY.a

CPPFLAGS += -Wall

default: depend $(LIB)

depend: $(DEP)

#....generate dependencies for each source code file
#    (make dependency files hidden to keep directory neater)
.%.d: %.cc
	$(CXX) $(CPPFLAGS) -MM -MG -MT "$@ $*.o" -MF $@ $<

$(LIB): $(OBJ)
	@if (test -f $(LIB) ) then (rm $(LIB)) fi
	@ar -crs $(LIB) $(OBJ)
	@echo "Done making $(LIB)"

clean:
	rm -f *.o .*.d ../lib/libCRY.a 

.PHONY: depend clean

-include $(DEP)
