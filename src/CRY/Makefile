#
# make = Create CRY library and run quick test
#
# Doug Wright

default: setup lib test

lib:
	$(MAKE) -C src

test:
	$(MAKE) -C test

clean:
	$(MAKE) -C src clean
	$(MAKE) -C test clean

setup:
	./setup.create

.PHONY: lib test clean

MAKEFLAGS += --no-print-directory

include src/Makefile.copyright
