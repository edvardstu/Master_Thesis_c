CC=gcc
FLAGS := -O3
LINKING := -lblas -fopenmp -lgsl -lgslcblas -lm -lhdf5 -I/usr/include/hdf5/serial/ -L/usr/lib/x86_64-linux-gnu/hdf5/serial
DEPS = utilities.h interactions.h str_builder.h systems.h hdf5ReaderWriter.h
OBJ = main.o utilities.o interactions.o str_builder.o systems.o hdf5ReaderWriter.o



.PHONY: run clean

%.o: %.c $(DEPS)
	$(CC) $(FLAGS) -c $< -o $@ $(LINKING)

main: $(OBJ)
	$(CC) $(FLAGS) -o $@ $^ $(LINKING)

run: main
	./main

clean:
	rm -f main *.o
