include $(LAWA_HOME)/Makefile.common

SOURCEFILES  = $(wildcard *.cc)
TARGETS      = $(patsubst %.cc,%,$(SOURCEFILES))

LDFLAGS += -L$(LAWA_HOME) -llawa -lextensionsflens -lextensionssparsegrid
EX_LDFLAGS += -bind_at_load

ifdef GMP_HOME
  CXXFLAGS += -DMPFR -I$(MPFR_HOME)/include
  LDFLAGS += -L$(MPFR_HOME)/lib -lgmpxx -lgmp
endif

CXXFLAGS += -DFLENS_IO_WITH_RANGES  -std=c++11

all: $(TARGETS) 

clean:
	$(RM) *.o
	$(RMDIR) *.dSYM
	$(RM) $(TARGETS)
    
