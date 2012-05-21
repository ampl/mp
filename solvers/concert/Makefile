SYSTEM    = x86-64_sles10_4.1
LIBFORMAT = static_pic

CPLEXSTUDIODIR = /opt/ibm/ILOG/CPLEX_Studio124
CONCERTDIR     = $(CPLEXSTUDIODIR)/concert
CPLEXDIR       = $(CPLEXSTUDIODIR)/cplex
CPOPTIMIZERDIR = $(CPLEXSTUDIODIR)/cpoptimizer
AMPLDIR        = ..

CONCERT_INC=$(CONCERTDIR)/include
CONCERT_LIB=$(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CPLEX_INC=$(CPLEXDIR)/include
CPLEX_LIB=$(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CPOPTIMIZER_INC=$(CPOPTIMIZERDIR)/include
CPOPTIMIZER_LIB=$(CPOPTIMIZERDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
#AMPL_INC=$(AMPLDIR)/ampl_include
AMPL_INC=$(AMPLDIR)
#AMPL_LIB=$(AMPLDIR)/amplsolver.a $(AMPLDIR)/funcadd0.o
AMPL_LIB=$(AMPLDIR)/amplsolver.a
LIBS=-lnsl -lm -ldl -lpthread

#DEBUG= -w -xO2 -xtarget=ultra -xarch=v8plus -DNDEBUG
DEBUG = -g
MTFLAGS = -DILCUSEMT
LDMTFLAGS = -mt

CXX_FLAGS = $(DEBUG) -I$(CONCERT_INC) -I$(CPLEX_INC) -I$(CPOPTIMIZER_INC) \
            -I$(AMPL_INC) $(MTFLAGS) -DIL_STD
LD_FLAGS  = -L$(CPOPTIMIZER_LIB) -lcp -L$(CPLEX_LIB) -lilocplex -lcplex \
            -L$(CONCERT_LIB) -lconcert $(LDMTFLAGS) $(AMPL_LIB) $(LIBS)

BINARIES = concert.o build_expr.o build_constr.o build_numberof.o main.o
EXTRAS = objconst0.o

cplexcp: $(BINARIES) $(EXTRAS)
	$(CXX) -o $@ $(BINARIES) $(EXTRAS) $(LD_FLAGS)

%.o: %.cpp
	$(CXX) $(CXX_FLAGS) -c $<
%.o: %.c
	$(CXX) $(CXX_FLAGS) -c $<
