################################################################################
# This configuration assumes GCC 6.2 or better.
# If you are using GCC 4.x, then 'std=gnu++0x' flag needs to be added to CCC.
#
# Configuration below is for CPLEX 128.
# Please set CPLEXROOT variable below to the correct path
################################################################################


CPLEXROOT		= /YOUR CPLEX PATH - CHANGE THIS TO YOUR "..."/cplex128 DIRECTORY PATH


SYSTEM			= x86-64_linux
LIBFORMAT		= static_pic
CPLEXDIR		= ${CPLEXROOT}/cplex
CONCERTDIR		= ${CPLEXROOT}/concert
CCC				= g++ -O4
CCOPT			= -m64 -O -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD
CPLEXBINDIR		= $(CPLEXDIR)/bin/$(BINDIST)
CPLEXJARDIR		= $(CPLEXDIR)/lib/cplex.jar
CPLEXLIBDIR		= $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR	= $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CCLNDIRS		= -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)
CLNDIRS			= -L$(CPLEXLIBDIR)
CCLNFLAGS		= -lconcert -lilocplex -lcplex -lm -lpthread -ldl
CONCERTINCDIR	= $(CONCERTDIR)/include
CPLEXINCDIR		= $(CPLEXDIR)/include
CCFLAGS			= $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) 

EXE1 = mcsc
OBJ1 = mcsc.o
SRC1 = mcsc.cpp

EXE2 = mcsi
OBJ2 = mcsi.o
SRC2 = mcsi.cpp

all: $(EXE1) $(EXE2)

only_mcsi: $(EXE2)

clean:
	rm -f *.o $(EXE1) $(EXE2)

$(EXE1): $(OBJ1)
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o $(EXE1) $(OBJ1) $(CCLNFLAGS)
$(OBJ1): $(SRC1)
	$(CCC) -c $(CCFLAGS) $(SRC1) -o $(OBJ1)

$(EXE2): $(OBJ2)
	$(CCC) -o $(EXE2) $(OBJ2)
$(OBJ2): $(SRC2)
	$(CCC) -c $(SRC2) -o $(OBJ2)
