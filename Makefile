LOCAL_DIR = ~/Projects/mdTrajToResponseFunction/
SRC_DIR = src/

VPATH = $(SRC_DIR)
CC = gcc
CFLAGS = -Wall -Wextra -ggdb
CFLAGS = -Wall -Wextra -O3 -ffast-math
LDFLAGS = -lm
#CFLAGS_FAST=-Wall -Wextra -O3 -ffast-math

H_FILES = globalArgs.h mymath.h
O_FILES = globalArgs.o mymath.o

#tools: mdTrajToFreq.fast mdTrajToFreq.fast freqTrajToResponse.fast freqTrajStats.fast
tools: mdTrajToFreq freqTrajToR5

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

mdTrajToFreq: mdTrajToFreq.o $(H_FILES) $(O_FILES) 
	$(CC) $(CFLAGS) $(LDFLAGS) mdTrajToFreq.o $(O_FILES) -o mdTrajToFreq

freqTrajToR5: freqTrajToR5.o $(H_FILES) $(O_FILES)
	$(CC) $(CFLAGS) $(LDFLAGS) freqTrajToR5.o $(O_FILES) -o freqTrajToR5

clean:
	rm -f *{.log,.aux,.dvi,.bbl,.blg,.lof,.lot,.toc,.o}

purge:
	rm -rf [!Makefile] 
