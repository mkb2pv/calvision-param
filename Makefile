ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
ROOTGLIBS  = $(shell root-config --glibs)
CXXFLAGS  += $(ROOTCFLAGS)
LIBS       = $(ROOTLIBS) 
GLIBS      = $(ROOTGLIBS)
GXX	   = g++ -Wall -O3

param: param.cpp
	$(GXX) -o param param.cpp $(ROOTCFLAGS) $(LIBS) $(ROOTGLIBS)

param2: param2.cpp
	$(GXX) -o param2 param2.cpp $(ROOTCFLAGS) $(LIBS) $(ROOTGLIBS)
