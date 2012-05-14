SYSTEM=ultrasparc_5_5.0
LIBFORMAT=static_pic_mt

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
LIBS=-lsocket -lnsl -lm -lsunmath -ldl

#DEBUG= -w -xO2 -xtarget=ultra -xarch=v8plus -DNDEBUG
DEBUG = -g
MTFLAGS = -DILCUSEMT
LDMTFLAGS = -mt

CXX_FLAGS = $(DEBUG) -I$(CONCERT_INC) -I$(CPLEX_INC) -I$(CPOPTIMIZER_INC) \
            -I$(AMPL_INC) $(MTFLAGS) -DIL_STD
LD_FLAGS  = -L$(CPOPTIMIZER_LIB) -lhybrid -lsolver -L$(CPLEX_LIB) -lilocplex \
            -lcplex -L$(CONCERT_LIB) -lconcert $(LDMTFLAGS) $(AMPL_LIB) $(LIBS)

BINARIES = concert.o build_expr.o build_constr.o build_numberof.o
EXTRAS = objconst0.o

concert: $(BINARIES) $(EXTRAS)
	$(CXX) -o concert $(BINARIES) $(EXTRAS) $(LD_FLAGS)

%.o: %.cpp
	$(CXX) $(CXX_FLAGS) -c $<
%.o: %.c
	$(CXX) $(CXX_FLAGS) -c $<
