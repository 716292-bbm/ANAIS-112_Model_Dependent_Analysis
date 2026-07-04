USR := $(shell whoami)
CC=gcc 
GSLPATH=/usr/local/include
LD_LIBRARY_PATH=/usr/local/lib
export LD_LIBRARY_PATH

ifneq (,$(findstring andrewcheek, $(USR)))
	LPATH=/sw/lib
	CFLAGS= -Wall -I$(GSLPATH) -std=gnu99 -m64 -O3

else 
	LPATH=/usr/local/lib
	CFLAGS= -Wall -I$(GSLPATH) -std=gnu99 -m64 -O4 -fpic

endif
#LPATH=/usr/local/lib
LDFLAGS=-L$(LPATH) -lgsl -lgslcblas -DHCUBATURE -lm 
SOURCES=main.c hcubature.c source/halo.c source/effFormFact.c source/difCrosSec.c source/difRateGen.c source/WriteData.c source/coeffs_eft.c source/phys_consts.c source/binning_general.c source/MenenFF.c

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=main

all: $(SOURCES) $(EXECUTABLE) 
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(CFLAGS) $(OBJECTS) $(LDFLAGS) -o $@
	
.c:
	$(CC) $(CFLAGS) $< -c $@
clean: 
	rm *.o $(EXECUTABLE) source/*.o
