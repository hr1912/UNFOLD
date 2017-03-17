#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include<sys/ipc.h>
#include<sys/shm.h>
#include<sys/types.h>
#include<time.h>

#define BUFFER_LENGTH 1024
#define MAP_LENGTH 0x100000000
#define CHROM_NUMBER 24
#define CHROM_LENGTH 250000000
#define MAX_THREAD 14
#define READS_NUMBER 1024
#define SEQ_LENGTH 128
#define BIN_NUMBER 5000
#define MAX_SAMPLE_NUMBER 100
#define DUPLICATION_LENGTH 750000000
#define NAME_LENGTH 256

#define MY_ID 13579

char *sample_name;
int *kmer_position;
char *kmer_chrom_index;
char *chrom[CHROM_NUMBER+1];
char *total_thread_buffer;
char *total_thread_barcode;
char *total_thread_quality;
long *total_thread_result;
long *total_thread_sample_sum;
long *map_sum_all;
long *sample_sum_all;
long *result_all;
unsigned char *duplication_table;
long *duplication_all;
long *total_thread_gc_num;
long *gc_sum_all;
long *total_thread_q20;
long *q20_sum_all;
long *total_thread_q30;
long *q30_sum_all;
long *total_thread_ch_gc_num;
long *ch_gc_sum_all;


typedef struct _PARAM
{
	char filename[BUFFER_LENGTH];
	char barcodelist[BUFFER_LENGTH];
	char out_prefix[BUFFER_LENGTH];
	char exp_file[BUFFER_LENGTH];
	char bin_file[BUFFER_LENGTH];
	int thread_num;
	int sample_num;
	int para_bits;
	int start;
}PARAM;
PARAM pk;

int read_barcodelist(char *barcodelist);
int read_reads_1(char *filename);
int read_reads_2(char *filename);
int filetest(char *filename);

void useage()
{
	printf("\n \033[33;1m%s\033[0;m\tversion: \033[36;1m%.1f\033[0;m\n", "uacd_map", 2.0);
	printf(" Last modified on : 2014.01.08\n");
	printf("\n Author: \033[34;1m%s\033[0;m\n","Pan Kai <pankai@novogene.cn> Ruan Hang <rogerholmes1036@gmail.com>");
	printf("\n Useage:\n");
	//printf("	-s  reads file\n");
	printf("	-l  reads file list\n");
	//printf("	-b  barcode list file\n");
	printf("	-o  output dir\n");
	//printf("	-B  bin stats file\n");
	//printf("	-e  experience value file\n");
	printf("	-t(optional,default 14)  thread num\n");
	printf("	-r(optional,default 0)  start position(start from 0)\n");
	printf("	-d(optional,default OFF)  double file mode\n");
	printf("	-h  help\n");
}

int main(int argc,char *argv[])
{
	int option,share_memory_id,i;
	long start,end;

	if(argc==1)
	{
		//printf("The number of parameters is not right!\n");
		useage();
		exit(0);
	}

	pk.para_bits=0;
	pk.thread_num=MAX_THREAD;
	pk.start=1;

	while((option=getopt(argc,argv,"s:l:b:t:o:e:B:r:dh"))!=-1) 
	{
		switch(option)
		{
/*		case 's':
			sscanf(optarg,"%s",pk.filename);
			pk.para_bits=pk.para_bits|0x1;
			break;*/
		case 'l':
			sscanf(optarg,"%s",pk.filename);
			pk.para_bits=pk.para_bits|0x2;
			break;
/*		case 'b':
			sscanf(optarg,"%s",pk.barcodelist);
			pk.para_bits=pk.para_bits|0x4;
			break;*/
		case 't':
			sscanf(optarg,"%d",&pk.thread_num);
			pk.para_bits=pk.para_bits|0x8;
			if(pk.thread_num>MAX_THREAD)
				pk.thread_num=MAX_THREAD;
			break;
		case 'o':
			sscanf(optarg,"%s",pk.out_prefix);
			pk.para_bits=pk.para_bits|0x10;
			break;
/*		case 'e':
			sscanf(optarg,"%s",pk.exp_file);
			pk.para_bits=pk.para_bits|0x20;
			break;
		case 'B':
			sscanf(optarg,"%s",pk.bin_file);
			pk.para_bits=pk.para_bits|0x40;
			break;*/
		case 'r':
			sscanf(optarg,"%d",&pk.start);
			pk.para_bits=pk.para_bits|0x80;
			break;
		case 'd':
			pk.para_bits=pk.para_bits|0x100;
			break;
		case 'h':
			useage();
			break;
		default:
			//printf("The parameter is not right!\n");
			useage();
		}
	}

	//if(((pk.para_bits&0x52)!=0x52)&&((pk.para_bits&0x55)!=0x55))
	if(((pk.para_bits&0x12)!=0x12))
	{
		printf("The number of parameters is not right: 0x%x\n", pk.para_bits);
		useage();
		exit(0);
	}

	char cmd[BUFFER_LENGTH];
	sprintf(cmd, "mkdir %s", pk.out_prefix);
	system(cmd);

	if((share_memory_id=shmget(MY_ID,0,0))==-1)
	{
		printf("Share memory don't exist.\n");
		exit(0);
	}

	printf("Share memory ID is %d\n",share_memory_id);

	kmer_position=(int*)shmat(share_memory_id,NULL,SHM_RDONLY);
	kmer_chrom_index=(char*)(kmer_position)+MAP_LENGTH*sizeof(int);
	chrom[1]=kmer_chrom_index+MAP_LENGTH*sizeof(char);
	for(i=2;i<=CHROM_NUMBER;i++)
	{
		chrom[i]=chrom[i-1]+CHROM_LENGTH;
	}

	start=time(0);
	if((pk.para_bits&0x5)==5)
	{
		if(read_barcodelist(pk.barcodelist)==0)
			goto end;

		total_thread_buffer=(char*)malloc(pk.thread_num*READS_NUMBER*SEQ_LENGTH*sizeof(char));
		total_thread_barcode=(char*)malloc(pk.thread_num*READS_NUMBER*SEQ_LENGTH*sizeof(char));
		total_thread_result=(long*)malloc(pk.thread_num*pk.sample_num*CHROM_NUMBER*BIN_NUMBER*sizeof(long));
		total_thread_sample_sum=(long*)malloc(pk.thread_num*pk.sample_num*sizeof(long));
		map_sum_all=(long*)malloc(pk.sample_num*sizeof(long));
		sample_sum_all=(long*)malloc(pk.sample_num*sizeof(long));
		
		read_reads_1(pk.filename);

		//stats();

		free(total_thread_buffer);
		free(total_thread_barcode);
		free(total_thread_result);
		free(total_thread_sample_sum);
		free(map_sum_all);
		free(sample_sum_all);
	}
	else if((pk.para_bits&0x2)==2)
	{
		if(filetest(pk.filename)==0)
			goto end;
		total_thread_buffer=(char*)malloc(pk.thread_num*READS_NUMBER*SEQ_LENGTH*sizeof(char));
		total_thread_quality=(char*)malloc(pk.thread_num*READS_NUMBER*SEQ_LENGTH*sizeof(char));
		total_thread_result=(long*)malloc(pk.thread_num*CHROM_NUMBER*BIN_NUMBER*sizeof(long));
		total_thread_sample_sum=(long*)malloc(pk.thread_num*sizeof(long));
		map_sum_all=(long*)malloc(pk.sample_num*sizeof(long));
		sample_sum_all=(long*)malloc(pk.sample_num*sizeof(long));
		result_all=(long*)malloc(pk.sample_num*CHROM_NUMBER*BIN_NUMBER*sizeof(long));
		duplication_table=(unsigned char*)malloc(DUPLICATION_LENGTH*sizeof(unsigned char));
		duplication_all=(long*)malloc(pk.sample_num*sizeof(long));
		total_thread_gc_num=(long*)malloc(pk.thread_num*sizeof(long));
		gc_sum_all=(long*)malloc(pk.sample_num*sizeof(long));
		total_thread_q20=(long*)malloc(pk.thread_num*sizeof(long));
		q20_sum_all=(long*)malloc(pk.sample_num*sizeof(long));
		total_thread_q30=(long*)malloc(pk.thread_num*sizeof(long));
		q30_sum_all=(long*)malloc(pk.sample_num*sizeof(long));
		total_thread_ch_gc_num=(long*)malloc(pk.thread_num*CHROM_NUMBER*sizeof(long));
		ch_gc_sum_all=(long*)malloc(pk.sample_num*CHROM_NUMBER*sizeof(long));

		read_reads_2(pk.filename);

		//stats();

		free(total_thread_buffer);
		free(total_thread_quality);
		free(total_thread_result);
		free(total_thread_sample_sum);
		free(map_sum_all);
		free(sample_sum_all);
		free(result_all);
		free(duplication_table);
		free(duplication_all);
		free(total_thread_gc_num);
		free(gc_sum_all);
		free(total_thread_q20);
		free(q20_sum_all);
		free(total_thread_q30);
		free(q30_sum_all);
		free(total_thread_ch_gc_num);
		free(ch_gc_sum_all);
		free(sample_name);
	}
	end=time(0);
	printf("The run time is %lds\n",end-start);

end:
	shmdt(kmer_position);

	return 1;
}

int filetest(char *filename)
{
	int temp_length,flag,i,j,k,sample_index,name_flag=0;
	char command[BUFFER_LENGTH],temp_buffer[BUFFER_LENGTH],name1[NAME_LENGTH],name2[NAME_LENGTH];
	FILE *fp;

	temp_length=strlen(filename);
	if((filename[temp_length-3]=='.')&&(filename[temp_length-2]=='g')&&(filename[temp_length-1]=='z'))
		flag=1;
	else flag=0;
	if(flag==1)
	{
		sprintf(command,"gzip -dc %s",filename);
		if((fp=popen(command,"r"))==NULL)
		{
			printf("Can't open %s\n",filename);
			return 0;
		}
	}
	else
	{
		if((fp=fopen(filename,"r"))==NULL)
		{
			printf("Can't open %s\n",filename);
			return 0;
		}
	}

	sample_index=0;
	while(fgets(temp_buffer,BUFFER_LENGTH,fp)!=NULL)
	{
		sample_index++;
	}

	if(sample_index==0)
	{
		if(flag==1)
			pclose(fp);
		else
			fclose(fp);
		printf("The list file is empty\n");
		return 0;
	}

	if((pk.para_bits&0x100)==0x100)
	{
		if((sample_index&0x1)==0x1)
		{
			if(flag==1)
				pclose(fp);
			else
				fclose(fp);
			printf("The file num is not right\n");
			return 0;
		}
	}

	if((pk.para_bits&0x100)==0x100)
	{
		pk.sample_num=(sample_index>>1);
	}
	else
	{
		pk.sample_num=sample_index;
	}

	sample_name=(char*)malloc(pk.sample_num*NAME_LENGTH*sizeof(char));

	fseek(fp,0,0);

	if((pk.para_bits&0x100)==0x100)
	{
		for(i=0;i<pk.sample_num;i++)
		{
			fgets(temp_buffer,BUFFER_LENGTH,fp);
			temp_length=strlen(temp_buffer);
			for(j=temp_length-1;j>=0;j--)
			{
				if(temp_buffer[j]=='/')
					break;
			}
			for(j++,k=0;j<temp_length;j++,k++)
			{
				if((temp_buffer[j]=='_')||(temp_buffer[j]=='.'))
					break;
				name1[k]=temp_buffer[j];
			}
			name1[k]='\0';

			fgets(temp_buffer,BUFFER_LENGTH,fp);
			temp_length=strlen(temp_buffer);
			for(j=temp_length-1;j>=0;j--)
			{
				if(temp_buffer[j]=='/')
					break;
			}
			for(j++,k=0;j<temp_length;j++,k++)
			{
				if((temp_buffer[j]=='_')||(temp_buffer[j]=='.'))
					break;
				name2[k]=temp_buffer[j];
			}
			name2[k]='\0';

			if(strcmp(name1,name2)==0)
			{
				strcpy(&sample_name[i*NAME_LENGTH],name1);
			}
			else
			{
				printf("line %d and line %d are not right.\n",2*i+1,2*i+2);
				name_flag=1;
			}
		}
		if(name_flag==1)
		{
			free(sample_name);
			return 0;
		}
	}
	else
	{
		for(i=0;i<pk.sample_num;i++)
		{
			fgets(temp_buffer,BUFFER_LENGTH,fp);
			temp_length=strlen(temp_buffer);
			for(j=temp_length-1;j>=0;j--)
			{
				if(temp_buffer[j]=='/')
					break;
			}
			for(j++,k=0;j<temp_length;j++,k++)
			{
				if((temp_buffer[j]=='_')||(temp_buffer[j]=='.'))
					break;
				name1[k]=temp_buffer[j];
			}
			name1[k]='\0';

			strcpy(&sample_name[i*NAME_LENGTH],name1);
		}
	}

	if(flag==1)
		pclose(fp);
	else
		fclose(fp);

	return 1;
}
