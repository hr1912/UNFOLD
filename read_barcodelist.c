#include<stdio.h>
#include<string.h>

#define BUFFER_LENGTH 1024
#define MAX_BARCODE_NUMBER 100
#define MAX_BARCODE_LENGTH 10

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
extern PARAM pk;

char barcode_map[MAX_BARCODE_NUMBER*MAX_BARCODE_LENGTH];
int barcode_table[484375];
int barcode_shift[9]={0,0,0,0,0,0,0,15625,93750};
int atonum[256];

int read_barcodelist(char *barcodelist)
{
	char command[BUFFER_LENGTH],buffer[BUFFER_LENGTH],barcode_name[MAX_BARCODE_LENGTH],barcode_temp[MAX_BARCODE_LENGTH];
	int temp_length,temp_num,barcode_num,flag,i,j;
	FILE *fb;

	temp_length=strlen(barcodelist);
	if((barcodelist[temp_length-3]=='.')&&(barcodelist[temp_length-2]=='g')&&(barcodelist[temp_length-1]=='z'))
		flag=1;
	else flag=0;

	if(flag==1)
	{
		sprintf(command,"gzip -dc %s",barcodelist);
		if((fb=popen(command,"r"))==NULL)
		{
			printf("Can't open %s\n",barcodelist);
			return 0;
		}
	}
	else
	{
		if((fb=fopen(barcodelist,"r"))==NULL)
		{
			printf("Can't open %s\n",barcodelist);
			return 0;
		}
	}

	memset(barcode_table,0xff,484375*sizeof(int));
	barcode_num=0;
	memset(atonum,0,256*sizeof(int));
	atonum[65]=0;
	atonum[67]=1;
	atonum[84]=2;
	atonum[71]=3;
	atonum[78]=4;

	while(fgets(buffer,BUFFER_LENGTH,fb)!=NULL)
	{
		temp_length=strlen(buffer);
		if(buffer[temp_length-1]=='\n')
			temp_length--;
		memcpy(barcode_name,buffer,temp_length);
		barcode_name[temp_length]='\0';
		memcpy(&barcode_map[barcode_num*MAX_BARCODE_LENGTH],barcode_name,temp_length+1);
		memcpy(barcode_temp,barcode_name,temp_length+1);
		temp_num=atonum[(int)(barcode_temp[temp_length-1])];
		for(j=temp_length-2;j>=0;j--)
			temp_num=temp_num*5+atonum[(int)(barcode_temp[j])];
		temp_num+=barcode_shift[temp_length];
		barcode_table[temp_num]=barcode_num;
		for(i=0;i<temp_length;i++)
		{
			if(barcode_name[i]!='A')
			{
				barcode_temp[i]='A';
				temp_num=atonum[(int)(barcode_temp[temp_length-1])];
				for(j=temp_length-2;j>=0;j--)
					temp_num=temp_num*5+atonum[(int)(barcode_temp[j])];
				temp_num+=barcode_shift[temp_length];
				barcode_table[temp_num]=barcode_num;
			}
			if(barcode_name[i]!='T')
			{
				barcode_temp[i]='T';
				temp_num=atonum[(int)(barcode_temp[temp_length-1])];
				for(j=temp_length-2;j>=0;j--)
					temp_num=temp_num*5+atonum[(int)(barcode_temp[j])];
				temp_num+=barcode_shift[temp_length];
				barcode_table[temp_num]=barcode_num;
			}
			if(barcode_name[i]!='C')
			{
				barcode_temp[i]='C';
				temp_num=atonum[(int)(barcode_temp[temp_length-1])];
				for(j=temp_length-2;j>=0;j--)
					temp_num=temp_num*5+atonum[(int)(barcode_temp[j])];
				temp_num+=barcode_shift[temp_length];
				barcode_table[temp_num]=barcode_num;
			}
			if(barcode_name[i]!='G')
			{
				barcode_temp[i]='G';
				temp_num=atonum[(int)(barcode_temp[temp_length-1])];
				for(j=temp_length-2;j>=0;j--)
					temp_num=temp_num*5+atonum[(int)(barcode_temp[j])];
				temp_num+=barcode_shift[temp_length];
				barcode_table[temp_num]=barcode_num;
			}
			if(barcode_name[i]!='N')
			{
				barcode_temp[i]='N';
				temp_num=atonum[(int)(barcode_temp[temp_length-1])];
				for(j=temp_length-2;j>=0;j--)
					temp_num=temp_num*5+atonum[(int)(barcode_temp[j])];
				temp_num+=barcode_shift[temp_length];
				barcode_table[temp_num]=barcode_num;
			}
			barcode_temp[i]=barcode_name[i];
		}
		barcode_num++;
	}

	pk.sample_num=barcode_num;

	if(flag==1)
		pclose(fb);
	else fclose(fb);

	if(pk.sample_num>MAX_BARCODE_NUMBER)
	{
		printf("barcode num is more than max barcode num.\n");
		return 0;
	}

	return 1;
}