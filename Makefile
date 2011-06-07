LOCAL_DIR = ~/Projects/mdTrajToResponseFunction/
SRC_DIR = src/

VPATH = $(SRC_DIR)
CC = gcc
CFLAGS_FAST=-Wall -O3 -ffast-math -lm

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

clean:
	rm -f *{.log,.aux,.dvi,.bbl,.blg,.lof,.lot,.toc,.o}

purge:
	rm -rf [!Makefile] 
