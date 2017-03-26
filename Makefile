assn3: nBody.c submit.c
	mpicc -g -o assn3 nBody.c submit.c -lm

clean:
	rm *.o
	rm assn3
