/*
Example of coupling C program and FreemFem++ script with mmap and semaphore

The c code is FFMaster.c
The FreeFem++ code is FFSlave.edp
Here FreeFem++ is a slave process

The compilation steps are
cc -c libff-mmap-semaphore.c
cc FFMaster.c -o FFMaster libff-mmap-semaphore.o -g -pthread
ff-c++ -auto ff-mmap-semaphore.cpp
./FFMaster

F. Hecht, Feb. 2018 - frederic.hecht@upmc.fr
*/

#include "libff-mmap-semaphore.h"
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
ff_Psem sem_ff, sem_c; //the semaphore for mutex

int main(int argc, const char ** argv)
{
	int debug = 0;
	ff_Pmmap shd;
	double cff, rff;
	long status;
	int i;
	if (argc > 1) debug = atoi(argv[1]);
	ff_mmap_sem_verb = debug;

	sem_ff = ffsem_malloc();
	sem_c = ffsem_malloc();
	shd = ffmmap_malloc();

	ffsem_init(sem_ff, "ff-slave1", 1);
	ffsem_init(sem_c, "ff-master1", 1);
	ffmmap_init(shd, "shared-data", 1024);

	status = 1;
	ffmmap_write(shd, &status, sizeof(status), 8);
	ffmmap_msync(shd, 0, 32);

	char ff[1024];
	sprintf(ff, "FreeFem++ FFSlave.edp -nw -ns -v %d&", debug);
	system(ff); //lauch FF++ in batch no graphics
	if(debug) printf("cc: before wait\n");

	if(debug) printf("cc: before wait 0 ff\n");
	ffsem_wait(sem_ff);

	for (i = 0; i < 10; ++i){
		printf(" iter : %d \n", i);
		cff = 10+i;
		ffmmap_write(shd, &cff, sizeof(cff), 0);
		ffsem_post(sem_c);

		if(debug) printf(" cc: before wait 2\n");
		ffsem_wait(sem_ff);
		ffmmap_read(shd, &rff, sizeof(rff), 16);
		printf(" iter = %d rff= %f\n", i, rff);
	}

	status = 0; //end
	ffmmap_write(shd, &status, sizeof(status), 8);
	ffsem_post(sem_c);
	printf("End Master \n");
	ffsem_wait(sem_ff);
	ffsem_del(sem_ff);
	ffsem_del(sem_c);
	ffmmap_del(shd);
	return 0;
}
