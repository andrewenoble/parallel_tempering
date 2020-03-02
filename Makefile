# check the hostname (hashed lines are comments)
# HN := $(shell /bin/hostname -s)

# set executable name
EXECNAME = parallel_tempering

# default compiler settings
CC = gcc
OPT = -O2 -fopenmp
DEBUG = -g -Wall -ansi
LDFLAGS = -lm -lgsl -lmpc
INCLUDE = -I/usr/include -I/usr/local/include -I../include

# ifeq ($(HN), moo)
#	CC = icc
#	OPT = -O2 -static
# endif

SRC = *.c
OBJS = $*(SRC).o

# generic compilation for data production
$(EXECNAME):
	$(CC) $(OPT) $(SRC) $(INCLUDE) -o $(EXECNAME) $(LDFLAGS)
	/bin/rm -rf *.o

# debugging compilation
db:
	$(CC) $(DEBUG) $(SRC) $(INCLUDE) -o $(EXECNAME).db $(LDFLAGS)
	/bin/rm -rf *.o

# clean up
clean:
	rm -rf *.o core *~ $(EXECNAME) $(EXECNAME).db $(EXECNAME).db.dSYM

