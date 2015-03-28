# this is <Makefile>
# ----------------------------------------------------------------------------
#   ($Id: Makefile 4001 2011-06-10 12:10:09Z lrehor $)
# 
# 25/10/2000 by Thomas Forbriger (IMGF Frankfurt)
# 
# libseife Makefile
#
# This library contains code which was extracted from the program seife.f
# by Erhard Wielandt. The original version was written in 1984 at ETH Zurich.
# The code can be obtained through 
# http://www.software-for-seismometry.de/
#
# The code from seife.f was extracted to a library by Wolfgang Friederich in
# 1996. The current version of libseife is a an extended code collection based
# on the abovementioned code.
# 
# REVISIONS and CHANGES
#    25/10/2000   V1.0   Thomas Forbriger
#    14/12/2007   V1.1   g77 compilation is the default now
#    17/12/2007   V1.2   moved to gfortran
#    18/09/2010   V1.3   start migration to new SVN scheme
#                        discard f2c option (fallback)
#                        check variables
#    15/11/2010   V1.4   do not use fdep.sh
#    17/11/2010   V1.5   libseifemon.a is out of use
#    17/01/2011   V1.6   distinguish library creation and installation
#    10/06/2011   V1.7   preparing Makefile for export of cseife to  
#                        the DENISE code by Lisa Rehor
# 
# ============================================================================
# 

.PHONY: all
all: install doc

.PHONY: doc
doc: 

LIBRARIES=libcseife.a
.PHONY: install
install: 
	$(MAKE) $(LIBRARIES)
	 
# ============================================================================

CLIBSRC=$(wildcard cseife*.c)
HEADERS=$(shell find . -name \*.h)

.PHONY: clean edit
clean:
	-find . -name \*.d | xargs --no-run-if-empty /bin/rm -v
	-find . -name \*.o | xargs --no-run-if-empty /bin/rm -v
	-find . -name \*.bak | xargs --no-run-if-empty /bin/rm -v
	-/bin/rm -fv flist
	-/bin/rm -fv *.a

edit: flist ; vim flist
flist: $(wildcard Makefile *.c *.h)
	echo $^ | tr ' ' '\n' | sort > $@

AR=ar
AS=as
RANLIB=ranlib
CC=mpicc

FLAGS += $(MYFLAGS)
CFLAGS += -O3 $(FLAGS)
CPPFLAGS=$(addprefix -I,.) \
  $(FLAGS)

.c.o:
	$(CC) $(CFLAGS) $< -c $(CPPFLAGS)


#======================================================================
# reinstall target
# is necessary in case of header file problems (e.g. remcmmnt not installed)
.PHONY: reinstall
reinstall:
	$(MAKE) clean
	$(MAKE) install

#======================================================================
# pattern rules
# -------------

%.d: %.c
	$(SHELL) -ec '$(CC) -M $(CPPFLAGS) $< \
      | sed '\''s,\($(notdir $*)\)\.o[ :]*,$(dir $@)\1.o $@ : ,g'\'' \
        > $@; \
      [ -s $@ ] || rm -f $@'

#======================================================================
# library part
# ------------

include $(patsubst %.c,%.d,$(CLIBSRC))

libcseife.a: $(patsubst %.c,%.o,$(CLIBSRC))

libcseife.a: 
	$(AR) rcv $@ $^
	$(RANLIB) $@

# 
# ----- END OF Makefile ----- 
