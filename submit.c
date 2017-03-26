/*
Assignment 3
Team Member 1 : Simon Huynh
Team Member 2 : Marco Simone
*/

#include "nBody.h"
#include <math.h>
#include <time.h>

void readnbody(double** s, double** v, double* m, int n) {
	int myrank;
	int nprocs;
	int i;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	// Open the input file
	FILE * fp;
	fp = fopen ("./input.txt", "r");
	if (fp == NULL) {
	  fprintf(stderr, "error, cannot find input file");
	}

	// This is an example of reading the body parameters from the input file.
	//if (myrank == 0) {
		for (i = 0; i < n; i++) {
			double x, y, z, vx, vy, vz, mass;

			int result = fscanf(fp, INPUT_BODY, &x, &y, &z, &vx, &vy, &vz, &mass);
			if (result != 7) {
				fprintf(stderr, "error reading body %d. Check if the number of bodies is correct.\n", i);
				exit(0);
			}
			if (i % nprocs == myrank) {
				s[i/nprocs][0] = x;
				s[i/nprocs][1] = y;
				s[i/nprocs][2] = z;
				v[i/nprocs][0] = vx;
				v[i/nprocs][1] = vy;
				v[i/nprocs][2] = vz;
				m[i/nprocs] = mass;
			}
		}
	//}
	fclose(fp);
}

void gennbody(double** s, double** v, double* m, int n) {
	//printf("Generate nBody initial condition here.\n");
	int myrank;
	int nprocs;
	int i;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	double x, y, z, mass, distance, theta;
	srand(time(NULL)-myrank);  //hallelujah
	for (i = 0; i < n/nprocs; i++) {
		mass = 1E30 * ( (double)rand() / (double)RAND_MAX );
		distance = 0.5E13 * ( (double)rand() / (double)RAND_MAX );
		theta = 2 * M_PI * ( (double)rand() / (double)RAND_MAX );
		//printf("mass: %f, distance: %f, theta: %f\n", mass, distance, theta);
		x = distance * cos(theta);
		y = distance * sin(theta);
		z = 1E11 * ( (double)rand() / (double)RAND_MAX ) - 0.5;
		//printf("x: %f, y: %f, z: %f\n", x, y, z);

		s[i][0] = x;
		s[i][1] = y;
		s[i][2] = z;
		m[i] = mass;
	}
}

void nbody(double** s, double** v, double* m, int n, int iter, int timestep) {
	int myrank;
	int nprocs;
	int i, j, k, mhs, w;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	double G = 6.674E-11;
	double dist, forcescalar;
	double tempvec[3] = {0, 0, 0};
	double force[n/nprocs][3];

	int counter;
	int counter2;
	int copyrino;

	double recvbuf[4*n/nprocs];
	double sendbuf[4*n/nprocs];
	double* temp;
	double* sendptr;
	double* recvptr;
	sendptr = sendbuf;
	recvptr = recvbuf;
	for (i = 0; i < iter; i++) {
		for (counter = 0; counter < n/nprocs; counter++) {
			for (counter2 = 0; counter2 < 3; counter2++) {
				force[counter][counter2]=0;
			}
		}
		//calculations with processor's local bodies
		for (k = 0; k < n/nprocs; k++) {
			for (mhs = 0; mhs < n/nprocs; mhs++) {
				if (k != mhs) {
					tempvec[0] = s[mhs][0] - s[k][0];
					tempvec[1] = s[mhs][1] - s[k][1];
					tempvec[2] = s[mhs][2] - s[k][2];
					dist = sqrt( tempvec[0]*tempvec[0] + tempvec[1]*tempvec[1] + tempvec[2]*tempvec[2] );
					if(isfinite(dist)){
						forcescalar = (G * m[k] * m[mhs]) / pow(dist,2);
						force[k][0] = force[k][0] + forcescalar*(tempvec[0]/dist);
						force[k][1] = force[k][1] + forcescalar*(tempvec[1]/dist);
						force[k][2] = force[k][2] + forcescalar*(tempvec[2]/dist);
					}
				}
			}
		}
		int barf;
		for (barf = 0; barf < n/nprocs; barf++) {
			sendbuf[4*barf]   = s[barf][0];
			sendbuf[4*barf+1] = s[barf][1];
			sendbuf[4*barf+2] = s[barf][2];
			sendbuf[4*barf+3] = m[barf];
		}
		//calculations + communication with foreign bodies
		for (j = 0; j < nprocs-1; j++) {
			if(myrank%2==0){
				if(myrank==0){
						MPI_Send(&sendbuf, 4*n/nprocs, MPI_DOUBLE, nprocs-1, 0, MPI_COMM_WORLD);
						MPI_Recv(&recvbuf, 4*n/nprocs, MPI_DOUBLE, myrank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}else{
						MPI_Send(&sendbuf, 4*n/nprocs, MPI_DOUBLE, myrank-1, 0, MPI_COMM_WORLD);
						MPI_Recv(&recvbuf, 4*n/nprocs, MPI_DOUBLE, myrank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
			}else{
				if(myrank==nprocs-1){
						MPI_Recv(&recvbuf, 4*n/nprocs, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						MPI_Send(&sendbuf, 4*n/nprocs, MPI_DOUBLE, myrank-1, 0, MPI_COMM_WORLD);
				}else{
						MPI_Recv(&recvbuf, 4*n/nprocs, MPI_DOUBLE, myrank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						MPI_Send(&sendbuf, 4*n/nprocs, MPI_DOUBLE, myrank-1, 0, MPI_COMM_WORLD);
				}
			}

			for (k = 0; k < n/nprocs; k++) {
				for (mhs = 0; mhs < n/nprocs; mhs++) {
					tempvec[0] = recvbuf[mhs*4] - s[k][0];
					tempvec[1] = recvbuf[mhs*4+1] - s[k][1];
					tempvec[2] = recvbuf[mhs*4+2] - s[k][2];
					dist = sqrt( tempvec[0]*tempvec[0] + tempvec[1]*tempvec[1] + tempvec[2]*tempvec[2] );
					if(isfinite(dist)){
						forcescalar = (G * m[k] * recvbuf[mhs*4+3]) / pow(dist,2);
						force[k][0] = force[k][0] + forcescalar*(tempvec[0]/dist);
						force[k][1] = force[k][1] + forcescalar*(tempvec[1]/dist);
						force[k][2] = force[k][2] + forcescalar*(tempvec[2]/dist);
					}
				}
			}
			//temp = sendptr;
			//sendptr = recvptr;
			//recvptr = temp;
			for (copyrino = 0; copyrino < 4*n/nprocs; copyrino++) {
				sendbuf[copyrino] = recvbuf[copyrino];
			}
		}

		for (w = 0; w < n/nprocs; w++) {
			//tempvec holds acceleration
			tempvec[0] = force[w][0]/m[w];
		 	tempvec[1] = force[w][1]/m[w];
	 		tempvec[2] = force[w][2]/m[w];
			//update velocities
			v[w][0] = v[w][0] + timestep * tempvec[0];
			v[w][1] = v[w][1] + timestep * tempvec[1];
			v[w][2] = v[w][2] + timestep * tempvec[2];
 			//update positions
			s[w][0] = s[w][0] + timestep * v[w][0];
			s[w][1] = s[w][1] + timestep * v[w][1];
			s[w][2] = s[w][2] + timestep * v[w][2];
		}

	}

	// This is an example of printing the body parameters to the stderr. Your code should print out the final body parameters
	// in the exact order as the input file. Since we are writing to the stderr in this case, rather than the stdout, make
	// sure you dont add extra debugging statements in stderr.

	//if (myrank == 0) {
	for (i = 0; i < n / nprocs; i++) {
		fprintf(stderr, OUTPUT_BODY, s[i][0], s[i][1], s[i][2], v[i][0], v[i][1], v[i][2], m[i], myrank);
	}
	//}
}
