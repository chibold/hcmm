# MakeFile

CC_FLAGS= -O2 #-Wall 
FOLDER = ./
CC= g++ $(CC_FLAGS)
REMOVE = rm

OBJECTS = communityMovGen-rel0.3.o socialNet.o xmlParser.o
OBJECTS_O = $(FOLDER)communityMovGen-rel0.3.o $(FOLDER)socialNet.o $(FOLDER)xmlParser.o

EXE=hcmm2

all: $(EXE) #clean

$(EXE): $(OBJECTS)
	$(CC) -o $(FOLDER)$@ $(OBJECTS_O)

communityMovGen-rel0.3.o: communityMovGen-rel0.3.cc
	$(CC) -o $(FOLDER)$@ -c communityMovGen-rel0.3.cc

socialNet.o: socialNet.cc
	$(CC) -o $(FOLDER)$@ -c socialNet.cc
	
xmlParser.o: xmlParser.cpp
	$(CC) -o $(FOLDER)$@ -c xmlParser.cpp
	
clean:
	$(REMOVE) $(OBJECTS_O)
