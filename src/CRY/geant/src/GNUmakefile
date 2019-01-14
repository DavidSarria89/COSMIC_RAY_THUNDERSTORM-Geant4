# $Id: GNUmakefile,v 1.1 1999/01/07 16:05:40 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := cosmic
G4TARGET := $(name)
G4EXLIB := true

EXTRALIBS += -L$(CRYHOME)/lib -lCRY
CPPFLAGS  += -I$(CRYHOME)/src

#....make sure geant uses static (not shared) libraries
export G4LIB_BUILD_STATIC=1
export G4LIB_BUILD_SHARED=

.PHONY: all
all: check lib bin

ifdef G4INSTALL
  include $(G4INSTALL)/config/binmake.gmk
endif

#
#....Check that GEANT and CRY are properly installed
#
check:
ifndef G4INSTALL
	@echo G4INSTALL not defined, please setup GEANT4 environment.; false
endif
ifndef CRYHOME
	@echo CRYHOME not defined, please setup CRY using: source ../setup; false
endif
