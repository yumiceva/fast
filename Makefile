CC = g++
CFLAGS = `root-config --cflags`
INC = -I`root-config --incdir`
LDFLAGS = `root-config --libs`
# setup boost libraries
# check version using scram tool info boost
#INC += -I/cvmfs/cms-ib.cern.ch/nweek-02480/slc6_amd64_gcc530/external/boost/1.63.0-mlhled2/include
INC += -I$(CMSSW_FWLITE_INCLUDE_PATH)
#LDFLAGS += -L$(CMSSW_RELEASE_BASE)/external/slc6_amd64_gcc491/lib/ -lboost_program_options -#lboost_date_time -lboost_iostreams
LDFLAGS += -L$(LD_LIBRARY_PATH) -lTreePlayer -lboost_program_options -lboost_date_time -lboost_iostreams

SOURCES = PUReweight.cc
#OBJECTS = $(SOURCES:.C=.o)
#OBJECTS += $(SOURCES:.cc=.o)
OBJECTS = PUReweight.o
#EXECUTABLE = tanalyzer


all: $(SOURCES) fast

print-%  : ; @echo $* = $($*)

fast: $(OBJECTS)
	$(CC) -o fast $(LDFLAGS) $(INC) $(CFLAGS) $(OBJECTS) fast.cc

#$(EXCUTABLE): $(OBJECTS)
#	$(CC) -o $(EXECUTABLE) $(LDFLAGS) $(INC) $(CFLAGS) $(OBJECTS)

.C.o:
	$(CC) -c $(INC) $(CFLAGS) $< -o $@

.cc.o:
	$(CC) -c $(INC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o

cleanall: clean
	rm -f fast
