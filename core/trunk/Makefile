CC=g++
RM=rm -rf
CFLAGS=-c -Wall -O -g
LDFLAGS=-g
SOURCES=ACEFile.cpp Configuration.cpp Contig.cpp ContigPrinter.cpp CSVWriter.cpp HaploType.cpp Logger.cpp QualitySNPpp.cpp SAMContig.cpp SAMRead.cpp SAMFile.cpp SeqRead.cpp Variation.cpp ContigProvider.cpp ReadGroup.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=QSNPng

all: $(SOURCES) $(EXECUTABLE)
    
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	$(RM) $(OBJECTS) $(EXECUTABLE)


