SYSTEM=ultrasparc_5_5.0
LIBFORMAT=static_pic_mt

CONCERTDIR=../concert10
CPLEXDIR=../cplex70
SOLVERDIR=../solver50
#AMPLDIR=../../../ampl
AMPLDIR=/home/scherzo/dmg/nsolvers

CONCERT_INC=$(CONCERTDIR)/include
CONCERT_LIB=$(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CPLEX_INC=$(CPLEXDIR)/include
CPLEX_LIB=$(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
SOLVER_INC=$(SOLVERDIR)/include
SOLVER_LIB=$(SOLVERDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
#AMPL_INC=$(AMPLDIR)/ampl_include
AMPL_INC=$(AMPLDIR)
#AMPL_LIB=$(AMPLDIR)/amplsolver.a $(AMPLDIR)/funcadd0.o
AMPL_LIB=$(AMPLDIR)/amplsolver.a
LIBS=-lsocket -lnsl -lm -lsunmath -ldl

CXX=CC
#DEBUG= -w -xO2 -xtarget=ultra -xarch=v8plus -DNDEBUG
DEBUG = -g
MTFLAGS = -DILCUSEMT -mt
LDMTFLAGS = -mt

CXX_FLAGS = $(DEBUG) -I$(CONCERT_INC) -I$(CPLEX_INC) -I$(SOLVER_INC) -I$(AMPL_INC) $(MTFLAGS)
LD_FLAGS  = -L$(SOLVER_LIB) -lhybrid -lsolver -L$(CPLEX_LIB) -lilocplex -lcplex -L$(CONCERT_LIB) -lconcert $(LDMTFLAGS) $(AMPL_LIB) $(LIBS)


BINARIES = concert.o build_expr.o build_constr.o build_numberof.o
EXTRAS = objconst0.o

concert: $(BINARIES) $(EXTRAS)
	$(CXX) -o concert $(BINARIES) $(EXTRAS) $(LD_FLAGS)

%.o: %.cpp
	$(CXX) $(CXX_FLAGS) -c $<
%.o: %.c
	$(CXX) $(CXX_FLAGS) -c $<
