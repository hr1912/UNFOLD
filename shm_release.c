#include<stdio.h>
#include<stdlib.h>
#include<sys/ipc.h>
#include<sys/shm.h>
#include<sys/types.h>

#define MY_ID 13579

int main(int argc,char *argv[])
{
	int share_memory_id;
	if((share_memory_id=shmget(MY_ID,0,0))==-1)
	{
		printf("Share memory don't exist.\n");
		exit(0);
	}
	shmctl(share_memory_id,IPC_RMID,NULL);
	printf("Remove share memory.\n");
	return 1;
}
