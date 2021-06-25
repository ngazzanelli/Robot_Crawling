#---------------------------------------------
# Target file to be compiled by default
#---------------------------------------------
MAIN = crawler
#---------------------------------------------
# CC is the compiler to be used
#---------------------------------------------
CC = gcc
#---------------------------------------------
# CFLAGS are the options passed to the compiler
#---------------------------------------------
CFLAGS = -Wall 
#---------------------------------------------
# OBJS are the object files to be linked
#---------------------------------------------
OBJ1 = crawler
OBJ2 = qlearn
OBJ3 = ptask
OBJ4 = dynamics
OBJ5 = interface
OBJS = $(OBJ1).o $(OBJ2).o $(OBJ3).o $(OBJ4).o $(OBJ5).o 
#---------------------------------------------
# LIBS are the external libraries to be used
#--------------------------------------------
LIBS = -lrt -lm `allegro-config --libs`
#---------------------------------------------
# Dependencies
#--------------------------------------------
$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) -o $(MAIN) $(OBJS) $(LIBS)

$(OBJ1).o: $(OBJ1).c qlearn.h ptask.h
	$(CC) $(CFLAGS) -c $(OBJ1).c

$(OBJ2).o: $(OBJ2).c qlearn.h
	$(CC) $(CFLAGS) -c $(OBJ2).c
	
$(OBJ3).o: $(OBJ3).c ptask.h
	$(CC) $(CFLAGS) -c $(OBJ3).c

$(OBJ4).o: $(OBJ4).c ptask.h
	$(CC) $(CFLAGS) -c $(OBJ4).c

$(OBJ5).o: $(OBJ5).c ptask.h 
	$(CC) $(CFLAGS) -c $(OBJ5).c
#-------------------------------------------
# Command that can be specified inline: make clean
#-------------------------------------------
clean: 
	rm -rf *.o $(MAIN)
