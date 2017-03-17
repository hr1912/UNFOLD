#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<pthread.h>
#include<time.h>

#define BUFFER_LENGTH 1024
#define CHROM_NUMBER 24
#define MAX_THREAD 14
#define READS_NUMBER 1024
#define SEQ_LENGTH 128
#define BIN_LENGTH 50000
#define BIN_NUMBER 5000
#define BLOCK_NUMBER 33
#define MAX_BARCODE_LENGTH 10
#define MAX_BARCODE_NUMBER 100
#define MAX_SAMPLE_NUMBER 100
#define DUPLICATION_LENGTH 750000000
#define CHROM_LENGTH 250000000
#define NAME_LENGTH 256

extern char *sample_name;
extern int *kmer_position;
extern char *kmer_chrom_index;
extern char *chrom[CHROM_NUMBER+1];
extern char *total_thread_buffer;
extern char *total_thread_barcode;
extern char *total_thread_quality;
extern long *total_thread_result;
extern long *total_thread_sample_sum;
extern long *map_sum_all;
extern long *sample_sum_all;
extern long *result_all;
extern unsigned char *duplication_table;
extern long *duplication_all;
extern long *total_thread_gc_num;
extern long *gc_sum_all;
extern long *total_thread_q20;
extern long *q20_sum_all;
extern long *total_thread_q30;
extern long *q30_sum_all;
extern long *total_thread_ch_gc_num;
extern long *ch_gc_sum_all;
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

extern int barcode_table[484375];
extern int barcode_shift[9];
extern int atonum[256];

unsigned int reverse_4[256]=
{
     170,     234,      42,     106,     186,     250,      58,     122,
     138,     202,      10,      74,     154,     218,      26,      90,
     174,     238,      46,     110,     190,     254,      62,     126,
     142,     206,      14,      78,     158,     222,      30,      94,
     162,     226,      34,      98,     178,     242,      50,     114,
     130,     194,       2,      66,     146,     210,      18,      82,
     166,     230,      38,     102,     182,     246,      54,     118,
     134,     198,       6,      70,     150,     214,      22,      86,
     171,     235,      43,     107,     187,     251,      59,     123,
     139,     203,      11,      75,     155,     219,      27,      91,
     175,     239,      47,     111,     191,     255,      63,     127,
     143,     207,      15,      79,     159,     223,      31,      95,
     163,     227,      35,      99,     179,     243,      51,     115,
     131,     195,       3,      67,     147,     211,      19,      83,
     167,     231,      39,     103,     183,     247,      55,     119,
     135,     199,       7,      71,     151,     215,      23,      87,
     168,     232,      40,     104,     184,     248,      56,     120,
     136,     200,       8,      72,     152,     216,      24,      88,
     172,     236,      44,     108,     188,     252,      60,     124,
     140,     204,      12,      76,     156,     220,      28,      92,
     160,     224,      32,      96,     176,     240,      48,     112,
     128,     192,       0,      64,     144,     208,      16,      80,
     164,     228,      36,     100,     180,     244,      52,     116,
     132,     196,       4,      68,     148,     212,      20,      84,
     169,     233,      41,     105,     185,     249,      57,     121,
     137,     201,       9,      73,     153,     217,      25,      89,
     173,     237,      45,     109,     189,     253,      61,     125,
     141,     205,      13,      77,     157,     221,      29,      93,
     161,     225,      33,      97,     177,     241,      49,     113,
     129,     193,       1,      65,     145,     209,      17,      81,
     165,     229,      37,     101,     181,     245,      53,     117,
     133,     197,       5,      69,     149,     213,      21,      85
};

unsigned int reverse_1[4]={2,3,0,1};

unsigned int table_4[256]=
{
	 1094795585,     1094795587,     1094795604,     1094795591,     1094796097,     1094796099,     1094796116,     1094796103,
     1094800449,     1094800451,     1094800468,     1094800455,     1094797121,     1094797123,     1094797140,     1094797127,
     1094926657,     1094926659,     1094926676,     1094926663,     1094927169,     1094927171,     1094927188,     1094927175,
     1094931521,     1094931523,     1094931540,     1094931527,     1094928193,     1094928195,     1094928212,     1094928199,
     1096040769,     1096040771,     1096040788,     1096040775,     1096041281,     1096041283,     1096041300,     1096041287,
     1096045633,     1096045635,     1096045652,     1096045639,     1096042305,     1096042307,     1096042324,     1096042311,
     1095188801,     1095188803,     1095188820,     1095188807,     1095189313,     1095189315,     1095189332,     1095189319,
     1095193665,     1095193667,     1095193684,     1095193671,     1095190337,     1095190339,     1095190356,     1095190343,
     1128350017,     1128350019,     1128350036,     1128350023,     1128350529,     1128350531,     1128350548,     1128350535,
     1128354881,     1128354883,     1128354900,     1128354887,     1128351553,     1128351555,     1128351572,     1128351559,
     1128481089,     1128481091,     1128481108,     1128481095,     1128481601,     1128481603,     1128481620,     1128481607,
     1128485953,     1128485955,     1128485972,     1128485959,     1128482625,     1128482627,     1128482644,     1128482631,
     1129595201,     1129595203,     1129595220,     1129595207,     1129595713,     1129595715,     1129595732,     1129595719,
     1129600065,     1129600067,     1129600084,     1129600071,     1129596737,     1129596739,     1129596756,     1129596743,
     1128743233,     1128743235,     1128743252,     1128743239,     1128743745,     1128743747,     1128743764,     1128743751,
     1128748097,     1128748099,     1128748116,     1128748103,     1128744769,     1128744771,     1128744788,     1128744775,
     1413562689,     1413562691,     1413562708,     1413562695,     1413563201,     1413563203,     1413563220,     1413563207,
     1413567553,     1413567555,     1413567572,     1413567559,     1413564225,     1413564227,     1413564244,     1413564231,
     1413693761,     1413693763,     1413693780,     1413693767,     1413694273,     1413694275,     1413694292,     1413694279,
     1413698625,     1413698627,     1413698644,     1413698631,     1413695297,     1413695299,     1413695316,     1413695303,
     1414807873,     1414807875,     1414807892,     1414807879,     1414808385,     1414808387,     1414808404,     1414808391,
     1414812737,     1414812739,     1414812756,     1414812743,     1414809409,     1414809411,     1414809428,     1414809415,
     1413955905,     1413955907,     1413955924,     1413955911,     1413956417,     1413956419,     1413956436,     1413956423,
     1413960769,     1413960771,     1413960788,     1413960775,     1413957441,     1413957443,     1413957460,     1413957447,
     1195458881,     1195458883,     1195458900,     1195458887,     1195459393,     1195459395,     1195459412,     1195459399,
     1195463745,     1195463747,     1195463764,     1195463751,     1195460417,     1195460419,     1195460436,     1195460423,
     1195589953,     1195589955,     1195589972,     1195589959,     1195590465,     1195590467,     1195590484,     1195590471,
     1195594817,     1195594819,     1195594836,     1195594823,     1195591489,     1195591491,     1195591508,     1195591495,
     1196704065,     1196704067,     1196704084,     1196704071,     1196704577,     1196704579,     1196704596,     1196704583,
     1196708929,     1196708931,     1196708948,     1196708935,     1196705601,     1196705603,     1196705620,     1196705607,
     1195852097,     1195852099,     1195852116,     1195852103,     1195852609,     1195852611,     1195852628,     1195852615,
     1195856961,     1195856963,     1195856980,     1195856967,     1195853633,     1195853635,     1195853652,     1195853639
};

int total_thread_reads_num[MAX_THREAD];
volatile int read_flag[MAX_THREAD];
long (*barcode_ch_bin_result)[CHROM_NUMBER][BIN_NUMBER];
long duplication_num;
pthread_spinlock_t lock;

void* thread_match_1(void *param)
{
	char *thread_buffer,*thread_barcode,*buffer,*buffer_1,*buffer_2;
	int temp_index[BLOCK_NUMBER],reverse_temp_index[BLOCK_NUMBER],temp_pos[BLOCK_NUMBER],reverse_temp_pos[BLOCK_NUMBER],temp[BLOCK_NUMBER];
	int thread_index,barcode_index,chrom_index,bin_index;
	int i,j,num,flag,num1,num2,temp_length,temp_num;
	unsigned int number[BLOCK_NUMBER],reverse_number[BLOCK_NUMBER];
	long *thread_result,*thread_sample_sum;
	unsigned long *p1,*p2;
	unsigned int *p3;
	unsigned char *p4;

	thread_index=*(int*)(param);
	thread_buffer=&total_thread_buffer[thread_index*READS_NUMBER*SEQ_LENGTH];
	thread_barcode=&total_thread_barcode[thread_index*READS_NUMBER*SEQ_LENGTH];
	thread_result=&total_thread_result[thread_index*pk.sample_num*CHROM_NUMBER*BIN_NUMBER];
	thread_sample_sum=&total_thread_sample_sum[thread_index*pk.sample_num];
	while(1)
	{
		while(read_flag[thread_index])
		{
		}
		num=total_thread_reads_num[thread_index];
		buffer=thread_buffer+pk.start;
		buffer_1=thread_barcode;
		j=0;
		while(j<num)
		{
			temp_length=strlen(buffer_1);
			buffer_2=buffer_1+temp_length-2;
			temp_num=atonum[(int)(*buffer_2--)];
			for(i=1;i<MAX_BARCODE_LENGTH;i++)
			{
				if((*buffer_2)==':')
					break;
				temp_num=(temp_num<<2)+temp_num+atonum[(int)(*buffer_2--)];
			}
			temp_num+=barcode_shift[i];
			barcode_index=barcode_table[temp_num];
			if(barcode_index==-1)
			{
				goto NEXT;
			}
			thread_sample_sum[barcode_index]++;
			i=0;
			i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');
			i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');
			i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');
			i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');

			i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');
			i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');
			i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');
			i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');

			i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');
			i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');
			i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');
			i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer)=='N');

			buffer-=47;

			if(i>0)
				goto NEXT;
			number[0]=(((unsigned int)(buffer[0])&0x6)>>1)|(((unsigned int)(buffer[1])&0x6)<<1)|(((unsigned int)(buffer[2])&0x6)<<3)|(((unsigned int)(buffer[3])&0x6)<<5)
					 |(((unsigned int)(buffer[4])&0x6)<<7)|(((unsigned int)(buffer[5])&0x6)<<9)|(((unsigned int)(buffer[6])&0x6)<<11)|(((unsigned int)(buffer[7])&0x6)<<13)
					 |(((unsigned int)(buffer[8])&0x6)<<15)|(((unsigned int)(buffer[9])&0x6)<<17)|(((unsigned int)(buffer[10])&0x6)<<19)|(((unsigned int)(buffer[11])&0x6)<<21)
					 |(((unsigned int)(buffer[12])&0x6)<<23)|(((unsigned int)(buffer[13])&0x6)<<25)|(((unsigned int)(buffer[14])&0x6)<<27)|(((unsigned int)(buffer[15])&0x6)<<29);
			reverse_number[32]=(reverse_4[number[0]&0xff]<<24)|(reverse_4[(number[0]&0xff00)>>8]<<16)|(reverse_4[(number[0]&0xff0000)>>16]<<8)|(reverse_4[(number[0]&0xff000000)>>24]);

			buffer+=16;
			number[1]=(number[0]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[31]=(reverse_number[32]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[2]=(number[1]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[30]=(reverse_number[31]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[3]=(number[2]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[29]=(reverse_number[30]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[4]=(number[3]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[28]=(reverse_number[29]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[5]=(number[4]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[27]=(reverse_number[28]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[6]=(number[5]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[26]=(reverse_number[27]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[7]=(number[6]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[25]=(reverse_number[26]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[8]=(number[7]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[24]=(reverse_number[25]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[9]=(number[8]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[23]=(reverse_number[24]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[10]=(number[9]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[22]=(reverse_number[23]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[11]=(number[10]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[21]=(reverse_number[22]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[12]=(number[11]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[20]=(reverse_number[21]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[13]=(number[12]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[19]=(reverse_number[20]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[14]=(number[13]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[18]=(reverse_number[19]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[15]=(number[14]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[17]=(reverse_number[18]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[16]=(number[15]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[16]=(reverse_number[17]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[17]=(number[16]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[15]=(reverse_number[16]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[18]=(number[17]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[14]=(reverse_number[15]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[19]=(number[18]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[13]=(reverse_number[14]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[20]=(number[19]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[12]=(reverse_number[13]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[21]=(number[20]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[11]=(reverse_number[12]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[22]=(number[21]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[10]=(reverse_number[11]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[23]=(number[22]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[9]=(reverse_number[10]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[24]=(number[23]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[8]=(reverse_number[9]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[25]=(number[24]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[7]=(reverse_number[8]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[26]=(number[25]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[6]=(reverse_number[7]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[27]=(number[26]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[5]=(reverse_number[6]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[28]=(number[27]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[4]=(reverse_number[5]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[29]=(number[28]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[3]=(reverse_number[4]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[30]=(number[29]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[2]=(reverse_number[3]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[31]=(number[30]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[1]=(reverse_number[2]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[32]=(number[31]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[0]=(reverse_number[1]<<2)|(reverse_1[(((unsigned int)(*buffer)&0x6)>>1)]);
			
			buffer-=47;

			temp_pos[0]=kmer_position[number[0]];
			temp_index[0]=kmer_chrom_index[number[0]];
			reverse_temp_pos[0]=kmer_position[reverse_number[0]];
			reverse_temp_index[0]=kmer_chrom_index[reverse_number[0]];
			temp_pos[1]=kmer_position[number[1]];
			temp_index[1]=kmer_chrom_index[number[1]];
			reverse_temp_pos[1]=kmer_position[reverse_number[1]];
			reverse_temp_index[1]=kmer_chrom_index[reverse_number[1]];
			temp_pos[2]=kmer_position[number[2]];
			temp_index[2]=kmer_chrom_index[number[2]];
			reverse_temp_pos[2]=kmer_position[reverse_number[2]];
			reverse_temp_index[2]=kmer_chrom_index[reverse_number[2]];
			temp_pos[3]=kmer_position[number[3]];
			temp_index[3]=kmer_chrom_index[number[3]];
			reverse_temp_pos[3]=kmer_position[reverse_number[3]];
			reverse_temp_index[3]=kmer_chrom_index[reverse_number[3]];
			temp_pos[4]=kmer_position[number[4]];
			temp_index[4]=kmer_chrom_index[number[4]];
			reverse_temp_pos[4]=kmer_position[reverse_number[4]];
			reverse_temp_index[4]=kmer_chrom_index[reverse_number[4]];
			temp_pos[5]=kmer_position[number[5]];
			temp_index[5]=kmer_chrom_index[number[5]];
			reverse_temp_pos[5]=kmer_position[reverse_number[5]];
			reverse_temp_index[5]=kmer_chrom_index[reverse_number[5]];
			temp_pos[6]=kmer_position[number[6]];
			temp_index[6]=kmer_chrom_index[number[6]];
			reverse_temp_pos[6]=kmer_position[reverse_number[6]];
			reverse_temp_index[6]=kmer_chrom_index[reverse_number[6]];
			temp_pos[7]=kmer_position[number[7]];
			temp_index[7]=kmer_chrom_index[number[7]];
			reverse_temp_pos[7]=kmer_position[reverse_number[7]];
			reverse_temp_index[7]=kmer_chrom_index[reverse_number[7]];
			temp_pos[8]=kmer_position[number[8]];
			temp_index[8]=kmer_chrom_index[number[8]];
			reverse_temp_pos[8]=kmer_position[reverse_number[8]];
			reverse_temp_index[8]=kmer_chrom_index[reverse_number[8]];
			temp_pos[9]=kmer_position[number[9]];
			temp_index[9]=kmer_chrom_index[number[9]];
			reverse_temp_pos[9]=kmer_position[reverse_number[9]];
			reverse_temp_index[9]=kmer_chrom_index[reverse_number[9]];
			temp_pos[10]=kmer_position[number[10]];
			temp_index[10]=kmer_chrom_index[number[10]];
			reverse_temp_pos[10]=kmer_position[reverse_number[10]];
			reverse_temp_index[10]=kmer_chrom_index[reverse_number[10]];
			temp_pos[11]=kmer_position[number[11]];
			temp_index[11]=kmer_chrom_index[number[11]];
			reverse_temp_pos[11]=kmer_position[reverse_number[11]];
			reverse_temp_index[11]=kmer_chrom_index[reverse_number[11]];
			temp_pos[12]=kmer_position[number[12]];
			temp_index[12]=kmer_chrom_index[number[12]];
			reverse_temp_pos[12]=kmer_position[reverse_number[12]];
			reverse_temp_index[12]=kmer_chrom_index[reverse_number[12]];
			temp_pos[13]=kmer_position[number[13]];
			temp_index[13]=kmer_chrom_index[number[13]];
			reverse_temp_pos[13]=kmer_position[reverse_number[13]];
			reverse_temp_index[13]=kmer_chrom_index[reverse_number[13]];
			temp_pos[14]=kmer_position[number[14]];
			temp_index[14]=kmer_chrom_index[number[14]];
			reverse_temp_pos[14]=kmer_position[reverse_number[14]];
			reverse_temp_index[14]=kmer_chrom_index[reverse_number[14]];
			temp_pos[15]=kmer_position[number[15]];
			temp_index[15]=kmer_chrom_index[number[15]];
			reverse_temp_pos[15]=kmer_position[reverse_number[15]];
			reverse_temp_index[15]=kmer_chrom_index[reverse_number[15]];
			temp_pos[16]=kmer_position[number[16]];
			temp_index[16]=kmer_chrom_index[number[16]];
			reverse_temp_pos[16]=kmer_position[reverse_number[16]];
			reverse_temp_index[16]=kmer_chrom_index[reverse_number[16]];
			temp_pos[17]=kmer_position[number[17]];
			temp_index[17]=kmer_chrom_index[number[17]];
			reverse_temp_pos[17]=kmer_position[reverse_number[17]];
			reverse_temp_index[17]=kmer_chrom_index[reverse_number[17]];
			temp_pos[18]=kmer_position[number[18]];
			temp_index[18]=kmer_chrom_index[number[18]];
			reverse_temp_pos[18]=kmer_position[reverse_number[18]];
			reverse_temp_index[18]=kmer_chrom_index[reverse_number[18]];
			temp_pos[19]=kmer_position[number[19]];
			temp_index[19]=kmer_chrom_index[number[19]];
			reverse_temp_pos[19]=kmer_position[reverse_number[19]];
			reverse_temp_index[19]=kmer_chrom_index[reverse_number[19]];
			temp_pos[20]=kmer_position[number[20]];
			temp_index[20]=kmer_chrom_index[number[20]];
			reverse_temp_pos[20]=kmer_position[reverse_number[20]];
			reverse_temp_index[20]=kmer_chrom_index[reverse_number[20]];
			temp_pos[21]=kmer_position[number[21]];
			temp_index[21]=kmer_chrom_index[number[21]];
			reverse_temp_pos[21]=kmer_position[reverse_number[21]];
			reverse_temp_index[21]=kmer_chrom_index[reverse_number[21]];
			temp_pos[22]=kmer_position[number[22]];
			temp_index[22]=kmer_chrom_index[number[22]];
			reverse_temp_pos[22]=kmer_position[reverse_number[22]];
			reverse_temp_index[22]=kmer_chrom_index[reverse_number[22]];
			temp_pos[23]=kmer_position[number[23]];
			temp_index[23]=kmer_chrom_index[number[23]];
			reverse_temp_pos[23]=kmer_position[reverse_number[23]];
			reverse_temp_index[23]=kmer_chrom_index[reverse_number[23]];
			temp_pos[24]=kmer_position[number[24]];
			temp_index[24]=kmer_chrom_index[number[24]];
			reverse_temp_pos[24]=kmer_position[reverse_number[24]];
			reverse_temp_index[24]=kmer_chrom_index[reverse_number[24]];
			temp_pos[25]=kmer_position[number[25]];
			temp_index[25]=kmer_chrom_index[number[25]];
			reverse_temp_pos[25]=kmer_position[reverse_number[25]];
			reverse_temp_index[25]=kmer_chrom_index[reverse_number[25]];
			temp_pos[26]=kmer_position[number[26]];
			temp_index[26]=kmer_chrom_index[number[26]];
			reverse_temp_pos[26]=kmer_position[reverse_number[26]];
			reverse_temp_index[26]=kmer_chrom_index[reverse_number[26]];
			temp_pos[27]=kmer_position[number[27]];
			temp_index[27]=kmer_chrom_index[number[27]];
			reverse_temp_pos[27]=kmer_position[reverse_number[27]];
			reverse_temp_index[27]=kmer_chrom_index[reverse_number[27]];
			temp_pos[28]=kmer_position[number[28]];
			temp_index[28]=kmer_chrom_index[number[28]];
			reverse_temp_pos[28]=kmer_position[reverse_number[28]];
			reverse_temp_index[28]=kmer_chrom_index[reverse_number[28]];
			temp_pos[29]=kmer_position[number[29]];
			temp_index[29]=kmer_chrom_index[number[29]];
			reverse_temp_pos[29]=kmer_position[reverse_number[29]];
			reverse_temp_index[29]=kmer_chrom_index[reverse_number[29]];
			temp_pos[30]=kmer_position[number[30]];
			temp_index[30]=kmer_chrom_index[number[30]];
			reverse_temp_pos[30]=kmer_position[reverse_number[30]];
			reverse_temp_index[30]=kmer_chrom_index[reverse_number[30]];
			temp_pos[31]=kmer_position[number[31]];
			temp_index[31]=kmer_chrom_index[number[31]];
			reverse_temp_pos[31]=kmer_position[reverse_number[31]];
			reverse_temp_index[31]=kmer_chrom_index[reverse_number[31]];
			temp_pos[32]=kmer_position[number[32]];
			temp_index[32]=kmer_chrom_index[number[32]];
			reverse_temp_pos[32]=kmer_position[reverse_number[32]];
			reverse_temp_index[32]=kmer_chrom_index[reverse_number[32]];

			flag=0;
			num1=num2=0;
			for(i=0;i<BLOCK_NUMBER;i++)
			{
				if(temp_index[i]==0)
					num1++;
				else if(temp_index[i]>0)
				{
					temp[num2]=i;
					num2++;
				}
				else
				{
				}
			}
			if((num1+num2==BLOCK_NUMBER)&&(num2>=1))
			{
				p1=(unsigned long*)(chrom[temp_index[temp[0]]]+temp_pos[temp[0]]-temp[0]);
				p2=(unsigned long*)(buffer);
				i=0;
				i+=((*p1++)==(*p2++));
				i+=((*p1++)==(*p2++));
				i+=((*p1++)==(*p2++));
				i+=((*p1++)==(*p2++));
				i+=((*p1++)==(*p2++));
				i+=((*p1)==(*p2));
				if(i==6)
				{
					flag++;
					chrom_index=temp_index[temp[0]];
					bin_index=(temp_pos[temp[0]]-temp[0]-pk.start)/BIN_LENGTH;
				}
			}

			num1=num2=0;
			for(i=0;i<BLOCK_NUMBER;i++)
			{
				if(reverse_temp_index[i]==0)
					num1++;
				else if(reverse_temp_index[i]>0)
				{
					temp[num2]=i;
					num2++;
				}
				else
				{
				}
			}
			if((num1+num2==BLOCK_NUMBER)&&(num2>=1))
			{
				p3=(unsigned int*)(chrom[reverse_temp_index[temp[0]]]+reverse_temp_pos[temp[0]]-temp[0]);
				p4=(unsigned char*)(reverse_number);
				i=0;
				i+=((*p3++)==table_4[(*p4++)]);
				i+=((*p3++)==table_4[(*p4++)]);
				i+=((*p3++)==table_4[(*p4++)]);
				i+=((*p3++)==table_4[(*p4)]);
				p4=(unsigned char*)(reverse_number+16);
				i+=((*p3++)==table_4[(*p4++)]);
				i+=((*p3++)==table_4[(*p4++)]);
				i+=((*p3++)==table_4[(*p4++)]);
				i+=((*p3++)==table_4[(*p4)]);
				p4=(unsigned char*)(reverse_number+32);
				i+=((*p3++)==table_4[(*p4++)]);
				i+=((*p3++)==table_4[(*p4++)]);
				i+=((*p3++)==table_4[(*p4++)]);
				i+=((*p3)==table_4[(*p4)]);
				if(i==12)
				{
					flag++;
					chrom_index=reverse_temp_index[temp[0]];
					bin_index=(reverse_temp_pos[temp[0]]-temp[0]-pk.start)/BIN_LENGTH;
				}
			}

			if(flag==1)
			{
				thread_result[barcode_index*120000+(chrom_index-1)*5000+bin_index]++;
			}
NEXT:
			buffer_1+=SEQ_LENGTH;
			buffer+=SEQ_LENGTH;
			j++;
		}
		if(num<READS_NUMBER)
			break;
		read_flag[thread_index]=1;
	}
}

void* thread_match_2(void *param)
{
	char *thread_buffer,*thread_quality,*buffer,*gc_buffer,*q_buffer;
	int temp_index[BLOCK_NUMBER],reverse_temp_index[BLOCK_NUMBER],temp_pos[BLOCK_NUMBER],reverse_temp_pos[BLOCK_NUMBER],temp[BLOCK_NUMBER];
	int thread_index,chrom_index,bin_index;
	int i,j,num,flag,num1,num2;
	long t1,t2,gc_num,q20_num,q30_num;
	unsigned int number[BLOCK_NUMBER],reverse_number[BLOCK_NUMBER];
	long *thread_result,*thread_sample_sum,*thread_ch_gc_num;
	unsigned long *p1,*p2;
	unsigned int *p3;
	unsigned char *p4;

	thread_index=*(int*)(param);
	thread_buffer=&total_thread_buffer[thread_index*READS_NUMBER*SEQ_LENGTH];
	thread_quality=&total_thread_quality[thread_index*READS_NUMBER*SEQ_LENGTH];
	thread_result=&total_thread_result[thread_index*CHROM_NUMBER*BIN_NUMBER];
	thread_sample_sum=&total_thread_sample_sum[thread_index];
	thread_ch_gc_num=&total_thread_ch_gc_num[thread_index*CHROM_NUMBER];

	q20_num=0;
	q30_num=0;
	while(1)
	{
		while(read_flag[thread_index])
		{
		}
		num=total_thread_reads_num[thread_index];
		(*thread_sample_sum)+=num;
		buffer=thread_buffer+pk.start;
		gc_buffer=thread_buffer;
		q_buffer=thread_quality;
		j=0;
		while(j<num)
		{
			q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;
			q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;
			q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;
			q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;
			q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;
			q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;
			q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;
			q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;
			q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;
			q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;
			q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;
			q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;
			q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;
			q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;
			q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;
			q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;
			q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;
			q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;
			q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;
			q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;
			q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;
			q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;
			q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;
			q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;
			q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);q_buffer++;q20_num+=(*q_buffer>=53);q30_num+=(*q_buffer>=63);
			q_buffer-=49;
			gc_num=0;
			gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);
			gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);
			gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);
			gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);
			gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);
			gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);
			gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);
			gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);
			gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);
			gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);
			gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);
			gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer++)&0x3)==3);
			gc_num+=(((*gc_buffer++)&0x3)==3);gc_num+=(((*gc_buffer)&0x3)==3);
			gc_buffer-=49;
			total_thread_gc_num[thread_index]+=gc_num;

			i=0;
			i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');
			i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');
			i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');
			i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');

			i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');
			i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');
			i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');
			i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');

			i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');
			i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');
			i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');
			i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer++)=='N');i+=((*buffer)=='N');

			buffer-=47;

			if(i>0)
				goto NEXT;
			number[0]=(((unsigned int)(buffer[0])&0x6)>>1)|(((unsigned int)(buffer[1])&0x6)<<1)|(((unsigned int)(buffer[2])&0x6)<<3)|(((unsigned int)(buffer[3])&0x6)<<5)
					 |(((unsigned int)(buffer[4])&0x6)<<7)|(((unsigned int)(buffer[5])&0x6)<<9)|(((unsigned int)(buffer[6])&0x6)<<11)|(((unsigned int)(buffer[7])&0x6)<<13)
					 |(((unsigned int)(buffer[8])&0x6)<<15)|(((unsigned int)(buffer[9])&0x6)<<17)|(((unsigned int)(buffer[10])&0x6)<<19)|(((unsigned int)(buffer[11])&0x6)<<21)
					 |(((unsigned int)(buffer[12])&0x6)<<23)|(((unsigned int)(buffer[13])&0x6)<<25)|(((unsigned int)(buffer[14])&0x6)<<27)|(((unsigned int)(buffer[15])&0x6)<<29);
			reverse_number[32]=(reverse_4[number[0]&0xff]<<24)|(reverse_4[(number[0]&0xff00)>>8]<<16)|(reverse_4[(number[0]&0xff0000)>>16]<<8)|(reverse_4[(number[0]&0xff000000)>>24]);

			buffer+=16;
			number[1]=(number[0]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[31]=(reverse_number[32]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[2]=(number[1]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[30]=(reverse_number[31]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[3]=(number[2]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[29]=(reverse_number[30]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[4]=(number[3]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[28]=(reverse_number[29]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[5]=(number[4]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[27]=(reverse_number[28]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[6]=(number[5]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[26]=(reverse_number[27]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[7]=(number[6]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[25]=(reverse_number[26]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[8]=(number[7]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[24]=(reverse_number[25]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[9]=(number[8]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[23]=(reverse_number[24]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[10]=(number[9]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[22]=(reverse_number[23]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[11]=(number[10]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[21]=(reverse_number[22]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[12]=(number[11]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[20]=(reverse_number[21]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[13]=(number[12]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[19]=(reverse_number[20]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[14]=(number[13]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[18]=(reverse_number[19]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[15]=(number[14]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[17]=(reverse_number[18]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[16]=(number[15]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[16]=(reverse_number[17]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[17]=(number[16]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[15]=(reverse_number[16]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[18]=(number[17]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[14]=(reverse_number[15]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[19]=(number[18]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[13]=(reverse_number[14]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[20]=(number[19]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[12]=(reverse_number[13]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[21]=(number[20]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[11]=(reverse_number[12]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[22]=(number[21]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[10]=(reverse_number[11]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[23]=(number[22]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[9]=(reverse_number[10]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[24]=(number[23]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[8]=(reverse_number[9]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[25]=(number[24]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[7]=(reverse_number[8]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[26]=(number[25]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[6]=(reverse_number[7]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[27]=(number[26]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[5]=(reverse_number[6]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[28]=(number[27]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[4]=(reverse_number[5]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[29]=(number[28]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[3]=(reverse_number[4]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[30]=(number[29]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[2]=(reverse_number[3]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[31]=(number[30]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[1]=(reverse_number[2]<<2)|(reverse_1[(((unsigned int)(*buffer++)&0x6)>>1)]);
			number[32]=(number[31]>>2)|(((unsigned int)(*buffer)&0x6)<<29);
			reverse_number[0]=(reverse_number[1]<<2)|(reverse_1[(((unsigned int)(*buffer)&0x6)>>1)]);

			buffer-=47;

			temp_pos[0]=kmer_position[number[0]];
			temp_index[0]=kmer_chrom_index[number[0]];
			reverse_temp_pos[0]=kmer_position[reverse_number[0]];
			reverse_temp_index[0]=kmer_chrom_index[reverse_number[0]];
			temp_pos[1]=kmer_position[number[1]];
			temp_index[1]=kmer_chrom_index[number[1]];
			reverse_temp_pos[1]=kmer_position[reverse_number[1]];
			reverse_temp_index[1]=kmer_chrom_index[reverse_number[1]];
			temp_pos[2]=kmer_position[number[2]];
			temp_index[2]=kmer_chrom_index[number[2]];
			reverse_temp_pos[2]=kmer_position[reverse_number[2]];
			reverse_temp_index[2]=kmer_chrom_index[reverse_number[2]];
			temp_pos[3]=kmer_position[number[3]];
			temp_index[3]=kmer_chrom_index[number[3]];
			reverse_temp_pos[3]=kmer_position[reverse_number[3]];
			reverse_temp_index[3]=kmer_chrom_index[reverse_number[3]];
			temp_pos[4]=kmer_position[number[4]];
			temp_index[4]=kmer_chrom_index[number[4]];
			reverse_temp_pos[4]=kmer_position[reverse_number[4]];
			reverse_temp_index[4]=kmer_chrom_index[reverse_number[4]];
			temp_pos[5]=kmer_position[number[5]];
			temp_index[5]=kmer_chrom_index[number[5]];
			reverse_temp_pos[5]=kmer_position[reverse_number[5]];
			reverse_temp_index[5]=kmer_chrom_index[reverse_number[5]];
			temp_pos[6]=kmer_position[number[6]];
			temp_index[6]=kmer_chrom_index[number[6]];
			reverse_temp_pos[6]=kmer_position[reverse_number[6]];
			reverse_temp_index[6]=kmer_chrom_index[reverse_number[6]];
			temp_pos[7]=kmer_position[number[7]];
			temp_index[7]=kmer_chrom_index[number[7]];
			reverse_temp_pos[7]=kmer_position[reverse_number[7]];
			reverse_temp_index[7]=kmer_chrom_index[reverse_number[7]];
			temp_pos[8]=kmer_position[number[8]];
			temp_index[8]=kmer_chrom_index[number[8]];
			reverse_temp_pos[8]=kmer_position[reverse_number[8]];
			reverse_temp_index[8]=kmer_chrom_index[reverse_number[8]];
			temp_pos[9]=kmer_position[number[9]];
			temp_index[9]=kmer_chrom_index[number[9]];
			reverse_temp_pos[9]=kmer_position[reverse_number[9]];
			reverse_temp_index[9]=kmer_chrom_index[reverse_number[9]];
			temp_pos[10]=kmer_position[number[10]];
			temp_index[10]=kmer_chrom_index[number[10]];
			reverse_temp_pos[10]=kmer_position[reverse_number[10]];
			reverse_temp_index[10]=kmer_chrom_index[reverse_number[10]];
			temp_pos[11]=kmer_position[number[11]];
			temp_index[11]=kmer_chrom_index[number[11]];
			reverse_temp_pos[11]=kmer_position[reverse_number[11]];
			reverse_temp_index[11]=kmer_chrom_index[reverse_number[11]];
			temp_pos[12]=kmer_position[number[12]];
			temp_index[12]=kmer_chrom_index[number[12]];
			reverse_temp_pos[12]=kmer_position[reverse_number[12]];
			reverse_temp_index[12]=kmer_chrom_index[reverse_number[12]];
			temp_pos[13]=kmer_position[number[13]];
			temp_index[13]=kmer_chrom_index[number[13]];
			reverse_temp_pos[13]=kmer_position[reverse_number[13]];
			reverse_temp_index[13]=kmer_chrom_index[reverse_number[13]];
			temp_pos[14]=kmer_position[number[14]];
			temp_index[14]=kmer_chrom_index[number[14]];
			reverse_temp_pos[14]=kmer_position[reverse_number[14]];
			reverse_temp_index[14]=kmer_chrom_index[reverse_number[14]];
			temp_pos[15]=kmer_position[number[15]];
			temp_index[15]=kmer_chrom_index[number[15]];
			reverse_temp_pos[15]=kmer_position[reverse_number[15]];
			reverse_temp_index[15]=kmer_chrom_index[reverse_number[15]];
			temp_pos[16]=kmer_position[number[16]];
			temp_index[16]=kmer_chrom_index[number[16]];
			reverse_temp_pos[16]=kmer_position[reverse_number[16]];
			reverse_temp_index[16]=kmer_chrom_index[reverse_number[16]];
			temp_pos[17]=kmer_position[number[17]];
			temp_index[17]=kmer_chrom_index[number[17]];
			reverse_temp_pos[17]=kmer_position[reverse_number[17]];
			reverse_temp_index[17]=kmer_chrom_index[reverse_number[17]];
			temp_pos[18]=kmer_position[number[18]];
			temp_index[18]=kmer_chrom_index[number[18]];
			reverse_temp_pos[18]=kmer_position[reverse_number[18]];
			reverse_temp_index[18]=kmer_chrom_index[reverse_number[18]];
			temp_pos[19]=kmer_position[number[19]];
			temp_index[19]=kmer_chrom_index[number[19]];
			reverse_temp_pos[19]=kmer_position[reverse_number[19]];
			reverse_temp_index[19]=kmer_chrom_index[reverse_number[19]];
			temp_pos[20]=kmer_position[number[20]];
			temp_index[20]=kmer_chrom_index[number[20]];
			reverse_temp_pos[20]=kmer_position[reverse_number[20]];
			reverse_temp_index[20]=kmer_chrom_index[reverse_number[20]];
			temp_pos[21]=kmer_position[number[21]];
			temp_index[21]=kmer_chrom_index[number[21]];
			reverse_temp_pos[21]=kmer_position[reverse_number[21]];
			reverse_temp_index[21]=kmer_chrom_index[reverse_number[21]];
			temp_pos[22]=kmer_position[number[22]];
			temp_index[22]=kmer_chrom_index[number[22]];
			reverse_temp_pos[22]=kmer_position[reverse_number[22]];
			reverse_temp_index[22]=kmer_chrom_index[reverse_number[22]];
			temp_pos[23]=kmer_position[number[23]];
			temp_index[23]=kmer_chrom_index[number[23]];
			reverse_temp_pos[23]=kmer_position[reverse_number[23]];
			reverse_temp_index[23]=kmer_chrom_index[reverse_number[23]];
			temp_pos[24]=kmer_position[number[24]];
			temp_index[24]=kmer_chrom_index[number[24]];
			reverse_temp_pos[24]=kmer_position[reverse_number[24]];
			reverse_temp_index[24]=kmer_chrom_index[reverse_number[24]];
			temp_pos[25]=kmer_position[number[25]];
			temp_index[25]=kmer_chrom_index[number[25]];
			reverse_temp_pos[25]=kmer_position[reverse_number[25]];
			reverse_temp_index[25]=kmer_chrom_index[reverse_number[25]];
			temp_pos[26]=kmer_position[number[26]];
			temp_index[26]=kmer_chrom_index[number[26]];
			reverse_temp_pos[26]=kmer_position[reverse_number[26]];
			reverse_temp_index[26]=kmer_chrom_index[reverse_number[26]];
			temp_pos[27]=kmer_position[number[27]];
			temp_index[27]=kmer_chrom_index[number[27]];
			reverse_temp_pos[27]=kmer_position[reverse_number[27]];
			reverse_temp_index[27]=kmer_chrom_index[reverse_number[27]];
			temp_pos[28]=kmer_position[number[28]];
			temp_index[28]=kmer_chrom_index[number[28]];
			reverse_temp_pos[28]=kmer_position[reverse_number[28]];
			reverse_temp_index[28]=kmer_chrom_index[reverse_number[28]];
			temp_pos[29]=kmer_position[number[29]];
			temp_index[29]=kmer_chrom_index[number[29]];
			reverse_temp_pos[29]=kmer_position[reverse_number[29]];
			reverse_temp_index[29]=kmer_chrom_index[reverse_number[29]];
			temp_pos[30]=kmer_position[number[30]];
			temp_index[30]=kmer_chrom_index[number[30]];
			reverse_temp_pos[30]=kmer_position[reverse_number[30]];
			reverse_temp_index[30]=kmer_chrom_index[reverse_number[30]];
			temp_pos[31]=kmer_position[number[31]];
			temp_index[31]=kmer_chrom_index[number[31]];
			reverse_temp_pos[31]=kmer_position[reverse_number[31]];
			reverse_temp_index[31]=kmer_chrom_index[reverse_number[31]];
			temp_pos[32]=kmer_position[number[32]];
			temp_index[32]=kmer_chrom_index[number[32]];
			reverse_temp_pos[32]=kmer_position[reverse_number[32]];
			reverse_temp_index[32]=kmer_chrom_index[reverse_number[32]];

			flag=0;
			num1=num2=0;
			for(i=0;i<BLOCK_NUMBER;i++)
			{
				if(temp_index[i]==0)
					num1++;
				else if(temp_index[i]>0)
				{
					temp[num2]=i;
					num2++;
				}
				else
				{
				}
			}
			if((num1+num2==BLOCK_NUMBER)&&(num2>=1))
			{
				p1=(unsigned long*)(chrom[temp_index[temp[0]]]+temp_pos[temp[0]]-temp[0]);
				p2=(unsigned long*)(buffer);
				i=0;
				i+=((*p1++)==(*p2++));
				i+=((*p1++)==(*p2++));
				i+=((*p1++)==(*p2++));
				i+=((*p1++)==(*p2++));
				i+=((*p1++)==(*p2++));
				i+=((*p1)==(*p2));
				if(i==6)
				{
					flag++;
					chrom_index=temp_index[temp[0]];
					bin_index=(temp_pos[temp[0]]-temp[0]-pk.start)/BIN_LENGTH;
					t1=(((long)(chrom_index)-1)*CHROM_LENGTH+((long)(temp_pos[temp[0]])-temp[0]-pk.start))>>3;
					t2=(temp_pos[temp[0]]-temp[0]-pk.start)&0x7;
					pthread_spin_lock(&lock);
					if((duplication_table[t1]&(1<<t2))==0)
						duplication_table[t1]=duplication_table[t1]|(1<<t2);
					else
					{
						duplication_num++;
						pthread_spin_unlock(&lock);
						goto NEXT;
					}
					pthread_spin_unlock(&lock);
				}
			}

			num1=num2=0;
			for(i=0;i<BLOCK_NUMBER;i++)
			{
				if(reverse_temp_index[i]==0)
					num1++;
				else if(reverse_temp_index[i]>0)
				{
					temp[num2]=i;
					num2++;
				}
				else
				{
				}
			}
			if((num1+num2==BLOCK_NUMBER)&&(num2>=1))
			{
				p3=(unsigned int*)(chrom[reverse_temp_index[temp[0]]]+reverse_temp_pos[temp[0]]-temp[0]);
				p4=(unsigned char*)(reverse_number);
				i=0;
				i+=((*p3++)==table_4[(*p4++)]);
				i+=((*p3++)==table_4[(*p4++)]);
				i+=((*p3++)==table_4[(*p4++)]);
				i+=((*p3++)==table_4[(*p4)]);
				p4=(unsigned char*)(reverse_number+16);
				i+=((*p3++)==table_4[(*p4++)]);
				i+=((*p3++)==table_4[(*p4++)]);
				i+=((*p3++)==table_4[(*p4++)]);
				i+=((*p3++)==table_4[(*p4)]);
				p4=(unsigned char*)(reverse_number+32);
				i+=((*p3++)==table_4[(*p4++)]);
				i+=((*p3++)==table_4[(*p4++)]);
				i+=((*p3++)==table_4[(*p4++)]);
				i+=((*p3)==table_4[(*p4)]);
				if(i==12)
				{
					flag++;
					chrom_index=reverse_temp_index[temp[0]];
					bin_index=(reverse_temp_pos[temp[0]]-temp[0]-pk.start)/BIN_LENGTH;
				}
			}

			if(flag==1)
			{
				thread_result[(chrom_index-1)*BIN_NUMBER+bin_index]++;
				thread_ch_gc_num[chrom_index-1]+=gc_num;
			}
NEXT:
			buffer+=SEQ_LENGTH;
			gc_buffer+=SEQ_LENGTH;
			q_buffer+=SEQ_LENGTH;
			j++;
		}
		if(num<READS_NUMBER)
			break;
		read_flag[thread_index]=1;
	}

	total_thread_q20[thread_index]+=q20_num;
	total_thread_q30[thread_index]+=q30_num;
}

int process_1(char *reads)
{
	char *thread_buffer,*thread_barcode,command[BUFFER_LENGTH],temp_buffer[BUFFER_LENGTH];
	int temp_length,flag,error,num,i,j,temp[MAX_THREAD];
	pthread_t thread_id[MAX_THREAD];
	pthread_attr_t attr;
	FILE *fp;

	temp_length=strlen(reads);
	if((reads[temp_length-3]=='.')&&(reads[temp_length-2]=='g')&&(reads[temp_length-1]=='z'))
		flag=1;
	else flag=0;
	if(flag==1)
	{
		sprintf(command,"gzip -dc %s",reads);
		if((fp=popen(command,"r"))==NULL)
		{
			printf("Can't open %s\n",reads);
			return 0;
		}
	}
	else
	{
		if((fp=fopen(reads,"r"))==NULL)
		{
			printf("Can't open %s\n",reads);
			return 0;
		}
	}

	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);

	for(i=0;i<pk.thread_num;i++)
	{
		temp[i]=i;
		read_flag[i]=1;
		error=pthread_create(&thread_id[i],&attr,thread_match_1,(void*)(&temp[i]));
		if(error)
		{
			printf("Create thread failed!\n");
			exit(0);
		}
	}

	while(1)
	{
		for(i=0;i<pk.thread_num;i++)
		{
			if(read_flag[i])
			{
				thread_buffer=&total_thread_buffer[i*READS_NUMBER*SEQ_LENGTH];
				thread_barcode=&total_thread_barcode[i*READS_NUMBER*SEQ_LENGTH];
				for(num=0;num<READS_NUMBER;num++)
				{
					if(fgets(&thread_barcode[num*SEQ_LENGTH],SEQ_LENGTH,fp)==NULL)
						goto end;
					fgets(&thread_buffer[num*SEQ_LENGTH],SEQ_LENGTH,fp);
					fgets(temp_buffer,BUFFER_LENGTH,fp);
					fgets(temp_buffer,BUFFER_LENGTH,fp);
				}
				total_thread_reads_num[i]=num;
				read_flag[i]=0;
			}
		}
	}

end:
	total_thread_reads_num[i]=num;
	read_flag[i]=0;

	for(j=0;j<i;j++)
	{
		while(!read_flag[j])
		{
		}
		total_thread_reads_num[j]=0;
		read_flag[j]=0;
	}
	for(j=i+1;j<pk.thread_num;j++)
	{
		while(!read_flag[j])
		{
		}
		total_thread_reads_num[j]=0;
		read_flag[j]=0;
	}

	for(i=0;i<pk.thread_num;i++)
	{
		pthread_join(thread_id[i],NULL);
	}

	pthread_attr_destroy(&attr);

	if(flag==1)
		pclose(fp);
	else fclose(fp);

	return 1;
}

int process_2(char *reads)
{
	char *thread_buffer,*thread_quality,command[BUFFER_LENGTH],temp_buffer[BUFFER_LENGTH];
	int temp_length,flag,error,num,i,j,temp[MAX_THREAD];
	pthread_t thread_id[MAX_THREAD];
	pthread_attr_t attr;
	FILE *fp;

	temp_length=strlen(reads);
	if((reads[temp_length-3]=='.')&&(reads[temp_length-2]=='g')&&(reads[temp_length-1]=='z'))
		flag=1;
	else flag=0;
	if(flag==1)
	{
		sprintf(command,"gzip -dc %s",reads);
		if((fp=popen(command,"r"))==NULL)
		{
			printf("Can't open %s\n",reads);
			return 0;
		}
	}
	else
	{
		if((fp=fopen(reads,"r"))==NULL)
		{
			printf("Can't open %s\n",reads);
			return 0;
		}
	}

	pthread_spin_init(&lock,0);
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);

	for(i=0;i<pk.thread_num;i++)
	{
		temp[i]=i;
		read_flag[i]=1;
		error=pthread_create(&thread_id[i],&attr,thread_match_2,(void*)(&temp[i]));
		if(error)
		{
			printf("Create thread failed!\n");
			exit(0);
		}
	}

	while(1)
	{
		for(i=0;i<pk.thread_num;i++)
		{
			if(read_flag[i])
			{
				thread_buffer=&total_thread_buffer[i*READS_NUMBER*SEQ_LENGTH];
				thread_quality=&total_thread_quality[i*READS_NUMBER*SEQ_LENGTH];
				for(num=0;num<READS_NUMBER;num++)
				{
					if(fgets(temp_buffer,BUFFER_LENGTH,fp)==NULL)
						goto end;
					fgets(&thread_buffer[num*SEQ_LENGTH],SEQ_LENGTH,fp);
					fgets(temp_buffer,BUFFER_LENGTH,fp);
					fgets(&thread_quality[num*SEQ_LENGTH],SEQ_LENGTH,fp);
				}
				total_thread_reads_num[i]=num;
				read_flag[i]=0;
			}
		}
	}

end:
	total_thread_reads_num[i]=num;
	read_flag[i]=0;

	for(j=0;j<i;j++)
	{
		while(!read_flag[j])
		{
		}
		total_thread_reads_num[j]=0;
		read_flag[j]=0;
	}
	for(j=i+1;j<pk.thread_num;j++)
	{
		while(!read_flag[j])
		{
		}
		total_thread_reads_num[j]=0;
		read_flag[j]=0;
	}

	for(i=0;i<pk.thread_num;i++)
	{
		pthread_join(thread_id[i],NULL);
	}

	pthread_spin_destroy(&lock);
	pthread_attr_destroy(&attr);

	if(flag==1)
		pclose(fp);
	else fclose(fp);

	return 1;
}

int read_reads_1(char *filename)
{
	int error,i,j,k,temp,temp_length,flag;
	char reads_file[BUFFER_LENGTH];
	long start,end;
	FILE *fp,*fq;

	temp_length=strlen(filename);
	if((filename[temp_length-3]=='.')&&(filename[temp_length-2]=='g')&&(filename[temp_length-1]=='z'))
		flag=1;
	else flag=0;

	memset(total_thread_result,0,pk.thread_num*pk.sample_num*CHROM_NUMBER*BIN_NUMBER*sizeof(long));
	memset(total_thread_sample_sum,0,pk.thread_num*pk.sample_num*sizeof(long));
	if(flag==1)
	{
		error=process_1(filename);
		if(error==0)
			return 0;
	}
	else
	{
		if((fq=fopen(filename,"r"))==NULL)
		{
			printf("Can't open %s\n",filename);
			return 0;
		}
		while(fgets(reads_file,BUFFER_LENGTH,fq)!=NULL)
		{
			start=time(0);
			temp_length=strlen(reads_file);
			if(reads_file[temp_length-1]=='\n')
				reads_file[temp_length-1]='\0';
			error=process_1(reads_file);
			if(error==0)
			{
				printf("%s process failed\n",reads_file);
			}
			else
			{
				end=time(0);
				printf("%s process finished %lds\n",reads_file,end-start);
			}
		}
		fclose(fq);
	}

	temp=pk.sample_num*CHROM_NUMBER*BIN_NUMBER;
	for(i=0;i<temp;i++)
		for(j=1;j<pk.thread_num;j++)
			total_thread_result[i]+=total_thread_result[j*temp+i];

	barcode_ch_bin_result=(long(*)[CHROM_NUMBER][BIN_NUMBER])total_thread_result;
	for(i=0;i<pk.sample_num;i++)
	{
		map_sum_all[i]=0;
		for(j=0;j<CHROM_NUMBER;j++)
			for(k=0;k<BIN_NUMBER;k++)
				map_sum_all[i]+=barcode_ch_bin_result[i][j][k];
	}

	for(i=0;i<pk.sample_num;i++)
	{
		sample_sum_all[i]=0;
		for(j=0;j<pk.thread_num;j++)
			sample_sum_all[i]+=total_thread_sample_sum[j*pk.sample_num+i];
	}

	fp=fopen("result","w");
	for(i=0;i<pk.sample_num;i++)
	{
		fprintf(fp,">%d_%ld_%ld\n",i,sample_sum_all[i],map_sum_all[i]);
		for(j=0;j<CHROM_NUMBER;j++)
			for(k=0;k<BIN_NUMBER;k++)
				fprintf(fp,"%ld\n",barcode_ch_bin_result[i][j][k]);
	}
	fclose(fp);

	printf("%s process finished\n",filename);

	return 1;
}

int read_reads_2(char *filename)
{
	int temp_length,flag,i,j,k,temp,error,error1,error2,sample_index;
	char command[BUFFER_LENGTH],reads_file[BUFFER_LENGTH];
	long start,end,chrom_base;
	FILE *fp,*fq;

	temp_length=strlen(filename);
	if((filename[temp_length-3]=='.')&&(filename[temp_length-2]=='g')&&(filename[temp_length-1]=='z'))
		flag=1;
	else flag=0;
	if(flag==1)
	{
		sprintf(command,"gzip -dc %s",filename);
		if((fq=popen(command,"r"))==NULL)
		{
			printf("Can't open %s\n",filename);
			return 0;
		}
	}
	else
	{
		if((fq=fopen(filename,"r"))==NULL)
		{
			printf("Can't open %s\n",filename);
			return 0;
		}
	}

	if((pk.para_bits&0x100)==0x100)
	{
		for(sample_index=0;sample_index<pk.sample_num;sample_index++)
		{
			start=time(0);
			memset(total_thread_result,0,pk.thread_num*CHROM_NUMBER*BIN_NUMBER*sizeof(long));
			memset(total_thread_sample_sum,0,pk.thread_num*sizeof(long));
			memset(duplication_table,0,DUPLICATION_LENGTH*sizeof(unsigned char));
			memset(total_thread_gc_num,0,pk.thread_num*sizeof(long));
			memset(total_thread_q20,0,pk.thread_num*sizeof(long));
			memset(total_thread_q30,0,pk.thread_num*sizeof(long));
			memset(total_thread_ch_gc_num,0,pk.thread_num*CHROM_NUMBER*sizeof(long));
			duplication_num=0;
			fgets(reads_file,BUFFER_LENGTH,fq);
			temp_length=strlen(reads_file);
			if(reads_file[temp_length-1]=='\n')
				reads_file[temp_length-1]='\0';
			error1=process_2(reads_file);
			if(error1==0)
			{
				printf("%s process failed\n",reads_file);
			}
			else
			{
				printf("%s process finished\n",reads_file);
			}
			fgets(reads_file,BUFFER_LENGTH,fq);
			temp_length=strlen(reads_file);
			if(reads_file[temp_length-1]=='\n')
				reads_file[temp_length-1]='\0';
			error2=process_2(reads_file);
			if(error2==0)
			{
				printf("%s process failed\n",reads_file);
			}
			else
			{
				printf("%s process finished\n",reads_file);
			}
			if((error1==1)&&(error2==1))
			{
				temp=CHROM_NUMBER*BIN_NUMBER;
				for(i=0;i<temp;i++)
					for(j=1;j<pk.thread_num;j++)
						total_thread_result[i]+=total_thread_result[j*temp+i];

				memcpy(&result_all[sample_index*CHROM_NUMBER*BIN_NUMBER],total_thread_result,CHROM_NUMBER*BIN_NUMBER*sizeof(long));

				map_sum_all[sample_index]=0;
				for(i=0;i<temp;i++)
					map_sum_all[sample_index]+=total_thread_result[i];

				sample_sum_all[sample_index]=0;
				for(i=0;i<pk.thread_num;i++)
					sample_sum_all[sample_index]+=total_thread_sample_sum[i];

				gc_sum_all[sample_index]=0;
				for(i=0;i<pk.thread_num;i++)
					gc_sum_all[sample_index]+=total_thread_gc_num[i];

				q20_sum_all[sample_index]=0;
				for(i=0;i<pk.thread_num;i++)
					q20_sum_all[sample_index]+=total_thread_q20[i];

				q30_sum_all[sample_index]=0;
				for(i=0;i<pk.thread_num;i++)
					q30_sum_all[sample_index]+=total_thread_q30[i];

				for(i=0;i<CHROM_NUMBER;i++)
				{
					ch_gc_sum_all[sample_index*CHROM_NUMBER+i]=0;
					for(j=0;j<pk.thread_num;j++)
						ch_gc_sum_all[sample_index*CHROM_NUMBER+i]+=total_thread_ch_gc_num[j*CHROM_NUMBER+i];
				}

				duplication_all[sample_index]=duplication_num;

				end=time(0);

				printf("%s process finished %lds\n",&sample_name[sample_index*NAME_LENGTH],end-start);
			}
		}
	}
	else
	{
		for(sample_index=0;sample_index<pk.sample_num;sample_index++)
		{
			start=time(0);
			memset(total_thread_result,0,pk.thread_num*CHROM_NUMBER*BIN_NUMBER*sizeof(long));
			memset(total_thread_sample_sum,0,pk.thread_num*sizeof(long));
			memset(duplication_table,0,DUPLICATION_LENGTH*sizeof(unsigned char));
			memset(total_thread_gc_num,0,pk.thread_num*sizeof(long));
			memset(total_thread_q20,0,pk.thread_num*sizeof(long));
			memset(total_thread_q30,0,pk.thread_num*sizeof(long));
			memset(total_thread_ch_gc_num,0,pk.thread_num*CHROM_NUMBER*sizeof(long));
			duplication_num=0;
			fgets(reads_file,BUFFER_LENGTH,fq);
			temp_length=strlen(reads_file);
			if(reads_file[temp_length-1]=='\n')
				reads_file[temp_length-1]='\0';
			error=process_2(reads_file);
			if(error==0)
			{
				printf("%s process failed\n",reads_file);
			}
			else
			{
				temp=CHROM_NUMBER*BIN_NUMBER;
				for(i=0;i<temp;i++)
					for(j=1;j<pk.thread_num;j++)
						total_thread_result[i]+=total_thread_result[j*temp+i];

				memcpy(&result_all[sample_index*CHROM_NUMBER*BIN_NUMBER],total_thread_result,CHROM_NUMBER*BIN_NUMBER*sizeof(long));

				map_sum_all[sample_index]=0;
				for(i=0;i<temp;i++)
					map_sum_all[sample_index]+=total_thread_result[i];

				sample_sum_all[sample_index]=0;
				for(i=0;i<pk.thread_num;i++)
					sample_sum_all[sample_index]+=total_thread_sample_sum[i];

				gc_sum_all[sample_index]=0;
				for(i=0;i<pk.thread_num;i++)
					gc_sum_all[sample_index]+=total_thread_gc_num[i];

				q20_sum_all[sample_index]=0;
				for(i=0;i<pk.thread_num;i++)
					q20_sum_all[sample_index]+=total_thread_q20[i];

				q30_sum_all[sample_index]=0;
				for(i=0;i<pk.thread_num;i++)
					q30_sum_all[sample_index]+=total_thread_q30[i];

				for(i=0;i<CHROM_NUMBER;i++)
				{
					ch_gc_sum_all[sample_index*CHROM_NUMBER+i]=0;
					for(j=0;j<pk.thread_num;j++)
						ch_gc_sum_all[sample_index*CHROM_NUMBER+i]+=total_thread_ch_gc_num[j*CHROM_NUMBER+i];
				}

				duplication_all[sample_index]=duplication_num;

				end=time(0);

				printf("%s process finished\n",reads_file);
				printf("%s process finished %lds\n",&sample_name[sample_index*NAME_LENGTH],end-start);
			}
		}
	}

	if(flag==1)
		pclose(fq);
	else
		fclose(fq);

	barcode_ch_bin_result=(long(*)[CHROM_NUMBER][BIN_NUMBER])result_all;

	char name[BUFFER_LENGTH];
	sprintf(name, "%s/bin.stats",pk.out_prefix);

	if ((fp=fopen(name,"w"))==NULL) {
		printf("Cannot write %s. Now exit to system...\n",name);
		exit(-1);
	}

	for(i=0;i<sample_index;i++)
	{
		fprintf(fp,">%d_%s_%ld_%ld_%ld_%ld_%ld_%ld\n",i,&sample_name[i*NAME_LENGTH],sample_sum_all[i],map_sum_all[i],duplication_all[i],gc_sum_all[i],q20_sum_all[i],q30_sum_all[i]);
		for(j=0;j<CHROM_NUMBER;j++)
			for(k=0;k<BIN_NUMBER;k++)
				fprintf(fp,"%ld\n",barcode_ch_bin_result[i][j][k]);
	}

	fclose(fp);

	sprintf(name, "%s/gc.stats",pk.out_prefix);

	if ((fp=fopen(name,"w"))==NULL) {
		printf("Cannot write %s. Now exit to system...\n",name);
		exit(-1);
	}

	fprintf(fp,"sample_id sample_name");
	for(j=0;j<CHROM_NUMBER;j++)
		fprintf(fp," chr%d_bases gc_bases gc_percent",j+1);
	fprintf(fp,"\n");

	for(i=0;i<sample_index;i++)
	{
		fprintf(fp,"%d %s",i,&sample_name[i*NAME_LENGTH]);
		for(j=0;j<CHROM_NUMBER;j++)
		{
			chrom_base=0;
			for(k=0;k<BIN_NUMBER;k++)
				chrom_base+=barcode_ch_bin_result[i][j][k];
			fprintf(fp," chr%d_%ld %ld %.2lf",j+1,chrom_base*50,ch_gc_sum_all[i*CHROM_NUMBER+j],ch_gc_sum_all[i*CHROM_NUMBER+j]*2.0/chrom_base);
		}
		fprintf(fp,"\n");
	}

	fclose(fp);

	return 1;
}
