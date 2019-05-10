SharedLibraryPath=../../../ShareLibrary/

BAMTOOLS=../../../ShareLibrary/bamtools/install/include/bamtools
BAMTOOLS_LD=../../../ShareLibrary/bamtools/install/lib64

CC = g++
INC_DIR += $(SharedLibraryPath) #Identify include path
INC_DIR += $(BAMTOOLS) #Identify include path
CFLAGS = -O3 -Wall -std=c++0x -I$(BAMTOOLS) -L$(BAMTOOLS_LD) -Wl,-rpath,$(BAMTOOLS_LD)
LIBS += -lz
LIBS += -lpthread

OBJ = ./obj/

all: EigenDel

clsbasealgorithm.o : $(SharedLibraryPath)clsbasealgorithm.cpp $(SharedLibraryPath)clsbasealgorithm.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clsbasealgorithm.cpp

clsfastareader.o: $(SharedLibraryPath)clsfastareader.cpp $(SharedLibraryPath)clsfastareader.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clsfastareader.cpp

clsfastqreader.o: $(SharedLibraryPath)clsfastqreader.cpp $(SharedLibraryPath)clsfastqreader.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clsfastqreader.cpp


clsbwa.o : $(SharedLibraryPath)clsbwa.cpp $(SharedLibraryPath)clsbwa.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clsbwa.cpp

clsvcf1000genome.o : $(SharedLibraryPath)clsvcf1000genome.cpp $(SharedLibraryPath)clsvcf1000genome.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clsvcf1000genome.cpp

clsreadconfigini.o : $(SharedLibraryPath)clsreadconfigini.cpp $(SharedLibraryPath)clsreadconfigini.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clsreadconfigini.cpp

clsmuscle.o : $(SharedLibraryPath)clsmuscle.cpp $(SharedLibraryPath)clsmuscle.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clsmuscle.cpp

clsvelvet.o : $(SharedLibraryPath)clsvelvet.cpp $(SharedLibraryPath)clsvelvet.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clsvelvet.cpp

clsblast.o : $(SharedLibraryPath)clsblast.cpp $(SharedLibraryPath)clsblast.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clsblast.cpp

clskmeralgorithm.o : $(SharedLibraryPath)clskmeralgorithm.cpp $(SharedLibraryPath)clskmeralgorithm.h
	$(CC) $(CFLAGS) -c $(SharedLibraryPath)clskmeralgorithm.cpp

clsparsebam.o: ../EigenDel/clsparsebam.cpp ../EigenDel/clsparsebam.h
	$(CC) $(CFLAGS) -c ../EigenDel/clsparsebam.cpp

clssvdeldetect.o: ../EigenDel/clssvdeldetect.cpp ../EigenDel/clssvdeldetect.h
	$(CC) $(CFLAGS) -c ../EigenDel/clssvdeldetect.cpp

clslearning.o: ../EigenDel/clslearning.cpp ../EigenDel/clslearning.h
	$(CC) $(CFLAGS) -c ../EigenDel/clslearning.cpp

clscomparison.o: ../EigenDel/clscomparison.cpp ../EigenDel/clscomparison.h
	$(CC) $(CFLAGS) -c ../EigenDel/clscomparison.cpp

clsdebug.o: ../EigenDel/clsdebug.cpp ../EigenDel/clsdebug.h
	$(CC) $(CFLAGS) -c ../EigenDel/clsdebug.cpp

clsconfig.o: ../EigenDel/clsconfig.cpp ../EigenDel/clsconfig.h
	$(CC) $(CFLAGS) -c ../EigenDel/clsconfig.cpp

main.o: ../EigenDel/main.cpp 
	$(CC) $(CFLAGS) -c ../EigenDel/main.cpp

EigenDel: clsbasealgorithm.o clsfastareader.o clsfastqreader.o clsbwa.o clsvcf1000genome.o clsreadconfigini.o clsmuscle.o clsvelvet.o clsblast.o clskmeralgorithm.o clsparsebam.o clssvdeldetect.o clslearning.o clscomparison.o clsdebug.o clsconfig.o main.o 
	$(CC) $(CFLAGS) -o EigenDel clsbasealgorithm.o clsfastareader.o clsfastqreader.o clsbwa.o clsvcf1000genome.o clsreadconfigini.o clsmuscle.o clsvelvet.o clsblast.o clskmeralgorithm.o clsparsebam.o clssvdeldetect.o clslearning.o clscomparison.o clsdebug.o clsconfig.o main.o -lz -lbamtools -lpthread


RM-F := rm -f
.PHONY : all deps objs clean veryclean rebuild info
clean :	
	@$(RM-F) *.o
	@$(RM-F) *.d

veryclean: clean
	@$(RM-F) CheckScafInfo


