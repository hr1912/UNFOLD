#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include<sys/ipc.h>
#include<sys/shm.h>
#include<sys/types.h>

#define BUFFER_LENGTH 1024
#define MAP_LENGTH 0x100000000
#define CHROM_NUMBER 24
#define CHROM_LENGTH 250000000

#define MY_ID 13579

int *kmer_position;
char *kmer_chrom_index;
char *chrom[CHROM_NUMBER+1];

int read_reference(char *reference);

void useage()
{
	
	printf("\n \033[33;1m%s\033[0;m\tversion: \033[36;1m%.1f\033[0;m\n", "shm_initialize", 2.0);	
	printf(" Last modified on : 2013.11.4\n");
	printf("\n Author: \033[34;1m%s\033[0;m\n","Pan Kai <pankai@novogene.cn>");
	printf(" Useage:\n");
	printf("     -r  hg19 reference file (All bases should be all upcased)\n");
	printf("     -h  help\n");
}

int main(int argc,char *argv[])
{
	int option,para_bits=0,share_memory_id,flag,i;
	char reference[BUFFER_LENGTH];

	if(argc==1)
	{
		//printf("The number of parameters is not right!\n");
		useage();
		exit(0);
	}

	while((option=getopt(argc,argv,"r:h"))!=-1) 
	{
		switch(option)
		{
		case 'r':
			sscanf(optarg,"%s",reference);
			para_bits=para_bits|0x1;
			break;
		case 'h':
			useage();
			break;
		default:
			//printf("The parameter is not right!\n");
			useage();
		}
	}

	if((para_bits&0x1)!=1)
	{
		//printf("The number of parameters is not right!\n");
		useage();
		exit(0);
	}

	if((share_memory_id=shmget(MY_ID,MAP_LENGTH*sizeof(int)+MAP_LENGTH*sizeof(char)+6000000000+1024,IPC_CREAT|IPC_EXCL|0666))==-1)
	{
		printf("Create share memory failed!\n");
		exit(0);
	}

	printf("Share memory ID is %d\n",share_memory_id);

	kmer_position=(int*)shmat(share_memory_id,NULL,0);
	memset(kmer_position,0xff,MAP_LENGTH*sizeof(int));

	kmer_chrom_index=(char*)(kmer_position)+MAP_LENGTH*sizeof(int);
	memset(kmer_chrom_index,0xff,MAP_LENGTH*sizeof(char));

	chrom[1]=kmer_chrom_index+MAP_LENGTH*sizeof(char);
	for(i=2;i<=CHROM_NUMBER;i++)
	{
		chrom[i]=chrom[i-1]+CHROM_LENGTH;
	}

	flag=read_reference(reference);

	if(flag==0)
	{
		shmdt(kmer_position);
		shmctl(share_memory_id,IPC_RMID,NULL);
		return 0;
	}

	shmdt(kmer_position);

	return 1;
}

int read_reference(char *reference)
{
	FILE *fp;
	char buffer[BUFFER_LENGTH],command[BUFFER_LENGTH];
	char *file_flag;
	char *p;
	int temp_length,chrom_length,i,j,k,flag;
	unsigned int number;
	long l,sum;

	temp_length=strlen(reference);
	if((reference[temp_length-3]=='.')&&(reference[temp_length-2]=='g')&&(reference[temp_length-1]=='z'))
		flag=1;
	else flag=0;
	if(flag==1)
	{
		sprintf(command,"gzip -dc %s",reference);
		if((fp=popen(command,"r"))==NULL)
		{
			printf("Can't open %s\n",reference);
			return 0;
		}
	}
	else
	{
		if((fp=fopen(reference,"r"))==NULL)
		{
			printf("Can't open %s\n",reference);
			return 0;
		}
	}

	file_flag=fgets(buffer,BUFFER_LENGTH,fp);
	if(file_flag==NULL)
	{
		printf("file is empty\n");
		if(flag==1)
			pclose(fp);
		else fclose(fp);
		return 0;
	}
	while(buffer[0]!='>')
	{
		file_flag=fgets(buffer,BUFFER_LENGTH,fp);
		if(file_flag==NULL)
		{
			if(flag==1)
				pclose(fp);
			else fclose(fp);
			return 0;
		}
	}
	
	for(i=1;i<=CHROM_NUMBER;i++)
	{
		p=chrom[i];
		chrom_length=0;
		file_flag=fgets(buffer,BUFFER_LENGTH,fp);
		while(buffer[0]!='>')
		{
			temp_length=strlen(buffer);
			if(buffer[temp_length-1]=='\n')
				temp_length--;
			memcpy(p,buffer,temp_length);
			p+=temp_length;
			chrom_length+=temp_length;
			file_flag=fgets(buffer,BUFFER_LENGTH,fp);
		}
		printf("Chrom%d read finished %d\n",i,chrom_length);
		j=0;
		p=chrom[i];
		while(1)
		{
LOOP:		while((j<=chrom_length-16)&&(p[j]=='N'))
			{
				j++;
			}
						
			if(j>chrom_length-16)
				goto end;
			for(k=0;k<16;k++)
				if(p[j+k]=='N')
					break;
			if(k<16)
			{
				j+=k+1;
				goto LOOP;
			}
			
			number=(((unsigned int)(p[j])&0x6)>>1)|(((unsigned int)(p[j+1])&0x6)<<1)|(((unsigned int)(p[j+2])&0x6)<<3)|(((unsigned int)(p[j+3])&0x6)<<5)
				  |(((unsigned int)(p[j+4])&0x6)<<7)|(((unsigned int)(p[j+5])&0x6)<<9)|(((unsigned int)(p[j+6])&0x6)<<11)|(((unsigned int)(p[j+7])&0x6)<<13)
				  |(((unsigned int)(p[j+8])&0x6)<<15)|(((unsigned int)(p[j+9])&0x6)<<17)|(((unsigned int)(p[j+10])&0x6)<<19)|(((unsigned int)(p[j+11])&0x6)<<21)
				  |(((unsigned int)(p[j+12])&0x6)<<23)|(((unsigned int)(p[j+13])&0x6)<<25)|(((unsigned int)(p[j+14])&0x6)<<27)|(((unsigned int)(p[j+15])&0x6)<<29);
			while(1)
			{				
				if(kmer_chrom_index[number]<0)
				{
					kmer_position[number]=j;
					kmer_chrom_index[number]=i;
				}
				else if(kmer_chrom_index[number]>0)
				{
					kmer_position[number]=-1;
					kmer_chrom_index[number]=0;
				}
				else
				{
				}
				j++;
				if(j>chrom_length-16)
					goto end;
				if(p[j+15]!='N')
				{
					number=(number>>2)|(((unsigned int)(p[j+15])&0x6)<<29);
				}
				else
				{
					j+=16;
					goto LOOP;
				}
			}
		}
end:
		printf("Chrom%d process finished.\n",i);
	}

	if(flag==1)
		pclose(fp);
	else fclose(fp);

	sum=0;
	for(l=0;l<0x100000000;l++)
	{
		sum+=(kmer_chrom_index[l]>0);
	}
	printf("unique %ld\n",sum);

	sum=0;
	for(l=0;l<0x100000000;l++)
	{
		sum+=(kmer_chrom_index[l]==0);
	}
	printf("repeat %ld\n",sum);

	return 1;
}
