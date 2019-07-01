CC = g++
#CC = cc
CFLAGS = -O3 -std=c++11 -lm -fopenmp -Wall -g -DWIN_SIZE=31 -DWIN_SIZE_PLUS=33 #-DCM_DEBUG_PRINT -DCM_DEBUG
#CFLAGS = -std=c++11 -O3 -lm -fopenmp -Wall -g -DWIN_SIZE=43 -DWIN_SIZE_PLUS=45 #-DCM_DEBUG_PRINT #-DCM_DEBUG #-DSD_COLL_CHECK
SRC = $(wildcard *.cpp)
HDR = $(wildcard *.h)
NAME = kmer
INCLUDE = -Iuthash

all: $(NAME)

$(NAME) : $(SRC) $(HDR)
	$(CC) $(CFLAGS) $(INCLUDE) -o $(NAME) $? -fopenmp 
clean:
	rm -f $(NAME)

distclean: clean
