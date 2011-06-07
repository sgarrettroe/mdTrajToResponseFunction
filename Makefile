LOCAL_DIR = ~/Projects/mdTrajToResponseFunction/
SRC_DIR = src/

VPATH = $(SRC_DIR)
CC = gcc
CFLAGS_FAST=-Wall -O3 -ffast-math -lm

H_FILES = globalArgs.h mymath.h
O_FILES = globalArgs.o mymath.o

%.fast: %.c
	$(CC) $(CFLAGS_FAST) $< -o $*

%.O3: %.c
	$(CC) -O3 -lm $< -o $*

%.O2: %.c
	$(CC) -O2 -lm $< -o $*

%.wall: %.c
	$(CC) -Wall -lm $< -o $*

%.gdb: %.c
	$(CC) -ggdb -Wall -lm $< -o $*

tools: mdTrajToFreq.fast mdTrajToFreq.fast freqTrajToResponse.fast freqTrajStats.fast

mdTrajToFreq: mdTrajToFreq.o $(H_FILES) $(O_FILES) 
	$(CC) -Wall -lm mdTrajToFreq.o $(O_FILES) -o mdTrajToFreq

freqTrajToR5: freqTrajToR5.o $(H_FILES) $(O_FILES)
	$(CC) -Wall -lm freqTrajToR5.o $(O_FILES) -o freqTrajToR5

clean:
	rm -f *{.log,.aux,.dvi,.bbl,.blg,.lof,.lot,.toc,.o}

purge:
	rm -rf [!Makefile] 
