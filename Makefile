ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
ROOTGLIBS  = $(shell root-config --glibs)
CXXFLAGS  += $(ROOTCFLAGS)
LIBS       = $(ROOTLIBS) 
GLIBS      = $(ROOTGLIBS)
GXX	   = g++ -Wall -O3

test_conv: test_conv.cpp
	$(GXX) -o test_conv test_conv.cpp $(ROOTCFLAGS) $(LIBS) $(ROOTGLIBS)

param: param.cpp
	$(GXX) -o param param.cpp $(ROOTCFLAGS) $(LIBS) $(ROOTGLIBS)
