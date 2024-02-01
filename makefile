CC = g++
OBJ = *.o
EXE = heat
FLAGS = -lm

all:${EXE}

heat: heat_propagation.cpp
	$(CC) -o $@ $^ $(FLAGS) 
clean:
	rm -f $(OBJ) $(EXE)