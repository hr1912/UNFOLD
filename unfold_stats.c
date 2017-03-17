#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#define FILE_NAME_LENGTH 1024
#define BUFFER_LENGTH 1024
#define NAME_LENGTH 256
#define CHROM_NUMBER 24
#define AUTOSOME_NUMBER 22
#define MAX_SAMPLE 600
#define MAX_BARCODE_LENGTH 10
#define BIN_NUMBER 5000

#define MIN_GC_PER 300
#define MAX_GC_PER 650
#define MED_GC_PER 411

#define NORMAL_VAL 6000000
#define MAX_LOOP 5


int barcode_ch_bin_result[MAX_SAMPLE][CHROM_NUMBER][BIN_NUMBER];
static double barcode_ch_bin_result_gc_corr[MAX_SAMPLE][CHROM_NUMBER][BIN_NUMBER];
int total_barcode_result[MAX_SAMPLE]={0};
int sum_reads[MAX_SAMPLE]={0};
int dup_reads[MAX_SAMPLE]={0};
int gc_bases[MAX_SAMPLE]={0};
int q20_bases[MAX_SAMPLE]={0};
int q30_bases[MAX_SAMPLE]={0};

static double cff_per[MAX_SAMPLE];
static double cff_per_gc[MAX_SAMPLE];

static double cff_per_byX[MAX_SAMPLE];
static double cff_per_gc_byX[MAX_SAMPLE];

static int sex[MAX_SAMPLE];

char *lib_name;

int total_barcode_index; 

char out_prefix[FILE_NAME_LENGTH];

extern int
lowess(double *x, double *y, size_t n,
       double f, size_t nsteps,
       double delta, double *ys, double *rw, double *res);

static double barcode_ch_bin_normal[MAX_SAMPLE][CHROM_NUMBER][BIN_NUMBER];
static double barcode_ch_bin_normal_gc_corr[MAX_SAMPLE][CHROM_NUMBER][BIN_NUMBER];

static double ch_bin_gc[CHROM_NUMBER][BIN_NUMBER];
static double ch_bin_n[CHROM_NUMBER][BIN_NUMBER];

static int bin_valid_flag[CHROM_NUMBER][BIN_NUMBER];

static double count_gc_bin[MAX_SAMPLE][MAX_GC_PER];
static double count_gc_bin_num[MAX_SAMPLE][MAX_GC_PER];
static double corr_f[MAX_SAMPLE][MAX_GC_PER];
 
//static int sample_valid_flag[MAX_SAMPLE];

static int barcode_conclusion[MAX_SAMPLE][10];

char exp_file[FILE_NAME_LENGTH];
char corr_file[FILE_NAME_LENGTH];
char flag_file[FILE_NAME_LENGTH];
int gc_verbose;

int loop;

typedef struct num_chrom_bin
{
	double num;
	int chrom;
	int bin;
} num_chrom_bin;

static void load_bin_gc(void) {

	char name[BUFFER_LENGTH],line[BUFFER_LENGTH];
	int chr_id;
	int bin_id ;
	double n_per,gc_per;	

	FILE *fp;
	
	sprintf(name, "/PROJ/SD/PrenatalDiagnosis/GC_bias_test/files/bin.stat");
	
	if ((fp = fopen(name, "r")) == NULL) {
		printf("Cannot open %s. Now exit to system...\n",name);
		exit(-1);
	}

	while(fgets(line,sizeof(line),fp)!=NULL) {
		if(line[0]=='c') {
			int match =sscanf(line,"chr%d.%d %*d %lf %lf", &chr_id, &bin_id, &n_per, &gc_per); 
			if (match != 4) {
				sscanf(line,"chr%c.%d %*d %lf %lf", &chr_id, &bin_id, &n_per, &gc_per);
				if (chr_id == 'X') 
					chr_id = 23;
				else if (chr_id == 'Y')
					chr_id = 24;
			}

			ch_bin_gc[chr_id-1][bin_id-1]=gc_per;
			ch_bin_n[chr_id-1][bin_id-1]=n_per;
		}
	}
	
	fclose(fp);
}

static int cmp_bin(const void *a, const void *b)
{
	num_chrom_bin *A, *B;
	A = (num_chrom_bin *)a;	
	B = (num_chrom_bin *)b;

	if (A->num > B->num)
		return -1;
	else if (A->num < B->num)
		return 1;
	else if (A->chrom > B->chrom)
		return -1;
	else if (A->chrom < B->chrom)
		return 1;
	else if (A->bin > B->bin)
		return -1;
	else if (A->bin < B->bin)
		return 1;
	else 
		return 0;
}


static void filt_abnormal_bins(void) {
	
	int i,j,k;

	double sample_num = total_barcode_index + 1;
	
	double sum[CHROM_NUMBER][BIN_NUMBER],mean[CHROM_NUMBER][BIN_NUMBER],var[CHROM_NUMBER][BIN_NUMBER],sd[CHROM_NUMBER][BIN_NUMBER];
	double cv[CHROM_NUMBER][BIN_NUMBER];

	long autosome_reads[MAX_SAMPLE];
	
	FILE *fp;
	char name[BUFFER_LENGTH];

	sprintf(name, "%s/bin.cv", out_prefix);

	if ((fp = fopen(name, "w")) == NULL) {
		printf("Cannot write to %s. Now exit to system...\n",name);
		exit(-1);
	}

	for(i=0;i<=total_barcode_index;i++) {
		autosome_reads[i] = 0;
		for(j=0;j<AUTOSOME_NUMBER;j++) 
			for(k=0;k<BIN_NUMBER;k++) 
				autosome_reads[i] += barcode_ch_bin_result[i][j][k];
	}

	for(i=0;i<=total_barcode_index;i++)
		for(j=0;j<AUTOSOME_NUMBER;j++)
			for(k=0;k<BIN_NUMBER;k++) 
				barcode_ch_bin_normal[i][j][k] = (double)barcode_ch_bin_result[i][j][k] / autosome_reads[i] * NORMAL_VAL;

	for (j=0;j<AUTOSOME_NUMBER;j++) 
		for (k=0;k<BIN_NUMBER;k++) {
			if (ch_bin_n[j][k] || !ch_bin_gc[j][k]) 
				continue;
			for(i=0;i<=total_barcode_index;i++) {
				sum[j][k] += barcode_ch_bin_normal[i][j][k];
				var[j][k] += barcode_ch_bin_normal[i][j][k] * barcode_ch_bin_normal[i][j][k];
			}
		}

	for (j=0;j<AUTOSOME_NUMBER;j++)
		for (k=0;k<BIN_NUMBER;k++) {
			if (ch_bin_n[j][k] || !ch_bin_gc[j][k] || sum[j][k] <= 0 )
				continue;
			mean[j][k] = sum[j][k] / sample_num;
			var[j][k] = var[j][k] / sample_num - mean[j][k] * mean[j][k];
			sd[j][k] = sqrt(var[j][k]);
			cv[j][k] = sd[j][k] / mean[j][k];
			if (cv[j][k] < 0.22) 
				bin_valid_flag[j][k] = 1;
			fprintf(fp, "%d\t%d\t%.4lf\t%.4lf\n", j, k, mean[j][k], cv[j][k] * 100);
		}

	fclose(fp);


	//output	

	for(j=AUTOSOME_NUMBER;j<CHROM_NUMBER;j++) 
		for(k=0;k<BIN_NUMBER;k++) 
			bin_valid_flag[j][k] = 1;

	sprintf(name, "%s/bin.flag", out_prefix);

	if ((fp = fopen(name, "w")) == NULL) {
		printf("Cannot write to %s. Now exit to system...\n",name);
		exit(-1);
	}
			
	for(j=0;j<CHROM_NUMBER;j++)
		for(k=0;k<BIN_NUMBER;k++)
			fprintf(fp,"%d %d %d\n", j, k, bin_valid_flag[j][k]);

	fclose(fp);

	if (flag_file[0]) {

	int flag;
	char line[BUFFER_LENGTH];

	if ((fp = fopen(flag_file, "r")) == NULL) {
		printf("Cannot open %s. Now exit to system...\n",flag_file);
		exit(-1);
	}

	while(fgets(line,sizeof(line),fp)!=NULL) {
		sscanf(line,"%d %d %d\n", &j, &k, &flag);
		bin_valid_flag[j][k] = flag;
	}

	}
		
}


static void count_gc_reads(void) {


	int i,j,k;
	long autosome_reads[MAX_SAMPLE];
	
	for(i=0;i<=total_barcode_index;i++) {
		autosome_reads[i] = 0;
		for(j=0;j<AUTOSOME_NUMBER;j++)
			for(k=0;k<BIN_NUMBER;k++) {
				if (!bin_valid_flag[j][k])
					continue; 
				autosome_reads[i] += barcode_ch_bin_result[i][j][k];
			}
	}	
	
	memset((double *)barcode_ch_bin_normal,0,(MAX_SAMPLE*CHROM_NUMBER*BIN_NUMBER)*sizeof(double));

	for(i=0;i<=total_barcode_index;i++) 
		for(j=0;j<CHROM_NUMBER;j++)
			for(k=0;k<BIN_NUMBER;k++) {
				if (!bin_valid_flag[j][k])
					continue;
				double gc_per = ch_bin_gc[j][k];
				barcode_ch_bin_normal[i][j][k] = (double)barcode_ch_bin_result[i][j][k] * NORMAL_VAL / autosome_reads[i];
				count_gc_bin[i][(int)(gc_per * 10)] += barcode_ch_bin_normal[i][j][k];
				count_gc_bin_num[i][(int)(gc_per * 10)] ++;
			}	
			

}
	
static void gc_bias_lowess(int s_id) {
	
	size_t i,j,k;
	size_t n = 0;
	
	double x[MAX_GC_PER];
	double y[MAX_GC_PER];
	double yy[MAX_GC_PER];

	char name[BUFFER_LENGTH];
	FILE *fp = NULL; //FILE *fp,*fp1;

	// filt bins with abnormal mapping per of same gc rate

	struct num_chrom_bin pk[CHROM_NUMBER*BIN_NUMBER];
	double ch_bin_result[CHROM_NUMBER][BIN_NUMBER];
	int pk_num=0;	


	for(j=0;j<CHROM_NUMBER;j++)
		for(k=0;k<BIN_NUMBER;k++) {
			ch_bin_result[j][k]=barcode_ch_bin_normal[s_id][j][k];
		}
		
	for(j=0;j<CHROM_NUMBER;j++)
		for(k=0;k<BIN_NUMBER;k++) {
			if (ch_bin_result[j][k]>0) {
				pk[pk_num].num=ch_bin_result[j][k];
				pk[pk_num].chrom=j;
				pk[pk_num].bin=k;
				pk_num++;
			}
		}

	qsort(pk,pk_num,sizeof(num_chrom_bin),cmp_bin);
	
	int median_chrom=pk[pk_num>>1].chrom;
	int median_bin=pk[pk_num>>1].bin;

	double median_gc = ch_bin_gc[median_chrom][median_bin];	

	printf("median GC %.1f\tmedian for %d bins of sample %d\n", median_gc, pk_num, s_id);	
				
	for(i=MIN_GC_PER; i<=MAX_GC_PER; i++) {
		if (!count_gc_bin_num[s_id][i]) 
			continue;
		x[n] = i;
		y[n] = count_gc_bin[s_id][i] / count_gc_bin_num[s_id][i];
		yy[n] = count_gc_bin[s_id][i];
		n++;
	}
	
	double f = 0.6;
	size_t nsteps = 3;
	
	double delta = 0.3;

	double ys[MAX_GC_PER];
	double rw[MAX_GC_PER];
	double res[MAX_GC_PER];

	lowess(x,y,n,
		f,nsteps,
		delta,ys,rw,res);

	int gc_index=0;

	for (i=0; i<=n; i++) 
		//if ((int)x[i] == (int)(median_gc*10)) {
		if ((int)x[i] == MED_GC_PER){
			gc_index = i;
			break;
		}		

	double median_num;
	median_num = ys[gc_index];

	printf("sample %d median_gc %d median_num %lf\n",s_id, (int)x[gc_index], median_num);

if (gc_verbose) {
	
	sprintf (name,"%s/gc/gc.%d.csv", out_prefix, s_id); 	

	if ((fp = fopen(name, "w")) == NULL) {
		printf("Cannot write to %s. Now exit to system...\n",name);
		exit(-1);
	}

	fprintf(fp,"gc,per,per_s,num\n");
}

/* mute
	sprintf (name,"%s/gc/corr_f.%d.csv", out_prefix, s_id); 	
	if ((fp1 = fopen(name, "w")) == NULL ) {
		printf("Cannot write to %s. Now exit to system...\n",name);
		exit(-1);
	}
*/
	for(i=0; i<=n; i++) {
		//fprintf(stderr,"%lf\t%lf\t%lf\n",x[i],y[i],ys[i]);
		if (ys[i] <=0 || x[i] < MIN_GC_PER)
			continue;
		if (gc_verbose) 
			fprintf(fp,"%lf,%.8lf,%.8lf,%.8lf\n",x[i],y[i],ys[i],yy[i]);
		corr_f[s_id][(int)x[i]] = median_num / ys[i] ;

//		fprintf(fp1, "%lf,%lf\n", x[i], corr_f[s_id][(int)x[i]]);	
	}

	if (gc_verbose)
		fclose(fp);
//	fclose(fp1);

// load experience correction factor
// OBSELETE

/*
	if (corr_file[0]) {
	char line[BUFFER_LENGTH];
	double gc, corr_factor; 
	
	sprintf(name,"%s",corr_file);
	
	if ((fp = fopen(name, "r")) == NULL) {
		printf("Cannot read %s. Now exit to system...\n",name);
		exit(-1);
	} 

	while(fgets(line,sizeof(line),fp)!=NULL) {
		sscanf(line,"%lf,%lf",&gc, &corr_factor);
		corr_f[(int)gc] = corr_factor;
	}	

	fclose(fp);

	}
*/

if (gc_verbose) 
{
	sprintf(name,"%s/gc/plot_gc.%d.R",out_prefix,s_id);

	if ((fp = fopen(name, "w")) == NULL) {
		printf("Cannot write to %s. Now exit to system...\n",name);
		exit(-1);
	}

	fprintf(fp,
"#! /usr/bin/env Rscript\n\
library(ggplot2)\n\
a<-read.csv(\"%s/gc/gc.%d.csv\")\n\n",out_prefix,s_id);

	fprintf(fp, 
"png(type=\"cairo\",\"%s/gc/gc.per.%d.png\")\n\
p<-ggplot(a,aes(x=gc,y=per))+geom_point(aes(),shape=1,size=2)+geom_line(aes(x=gc,y=per_s,color=\"smoothed\"))\n\
p\n\
dev.off()\n\n",out_prefix,s_id);

/*
	fprintf(fp,
"png(type=\"cairo\",\"%s/gc/gc.num.%d.png\")\n\
p<-ggplot(a,aes(x=gc,y=num))+geom_point(aes(),shape=1,size=2)\n\
p\n\
dev.off()\n\n",out_prefix,s_id);
*/
	fclose(fp);
	
	char cmd[1024];
	sprintf(cmd,"chmod 755 %s &&Rscript %s", name,name);
	
	system(cmd);	
}
	
}

static void gc_corr(void) {

	int i,j,k;
	double autosome_reads[MAX_SAMPLE];

	for(i=0;i<=total_barcode_index;i++) {
		for(j=0;j<CHROM_NUMBER;j++) 
			for(k=0;k<BIN_NUMBER;k++) {
				double gc_per = ch_bin_gc[j][k];
				if (corr_f[i][(int)(gc_per*10)])
					barcode_ch_bin_result_gc_corr[i][j][k] = barcode_ch_bin_result[i][j][k] * corr_f[i][(int)(gc_per*10)];
				else {
					fprintf(stderr, "sample %d chrom %d bin %d gc value missing %.1f\n", i, j, k, gc_per);
					barcode_ch_bin_result_gc_corr[i][j][k] = barcode_ch_bin_result[i][j][k];
				}
			}
	}	

	for(i=0;i<=total_barcode_index;i++) {
		autosome_reads[i] = 0;
		for(j=0;j<AUTOSOME_NUMBER;j++)
			for(k=0;k<BIN_NUMBER;k++) {
				if (!bin_valid_flag[j][k]) 
					continue; 
				autosome_reads[i] += barcode_ch_bin_result_gc_corr[i][j][k];
			}
	}		

	for (i=0;i<=total_barcode_index;i++) {
		double tmp=0;
		for (j=0;j<CHROM_NUMBER;j++)
			for (k=0;k<BIN_NUMBER;k++) {
				if (!bin_valid_flag[j][k])
					continue;
				tmp += barcode_ch_bin_result_gc_corr[i][j][k] / autosome_reads[i];
			}
		printf("sample %d tmp %.6f\n", i, tmp);
	}
	
	for(i=0;i<=total_barcode_index;i++) 
		for(j=0;j<CHROM_NUMBER;j++)
			for(k=0;k<BIN_NUMBER;k++)  
				barcode_ch_bin_normal_gc_corr[i][j][k] = barcode_ch_bin_result_gc_corr[i][j][k] *NORMAL_VAL /autosome_reads[i];

}

void stats(void) {

	int i,j,k;

	load_bin_gc();

	filt_abnormal_bins();

	count_gc_reads();

	for (i=0;i<=total_barcode_index;i++)	
		gc_bias_lowess(i);		
	
	gc_corr();

	double barcode_ch_normal[MAX_SAMPLE][CHROM_NUMBER];
	double barcode_ch_normal_gc_corr[MAX_SAMPLE][CHROM_NUMBER];
	
	double sum[CHROM_NUMBER], sum_ir[CHROM_NUMBER], sum_gc[CHROM_NUMBER];
	double var[CHROM_NUMBER], var_ir[CHROM_NUMBER], var_gc[CHROM_NUMBER];
	double sd[CHROM_NUMBER], sd_ir[CHROM_NUMBER], sd_gc[CHROM_NUMBER]; 
	double mean[CHROM_NUMBER], mean_ir[CHROM_NUMBER], mean_gc[CHROM_NUMBER];

	// sex chromsome related variables
 
	int male_num=0,  female_num=0;
	
	double sum_f[CHROM_NUMBER], sum_f_ir[CHROM_NUMBER], sum_f_gc[CHROM_NUMBER];
	double var_f[CHROM_NUMBER], var_f_ir[CHROM_NUMBER], var_f_gc[CHROM_NUMBER];
	double sd_f[CHROM_NUMBER], sd_f_ir[CHROM_NUMBER], sd_f_gc[CHROM_NUMBER]; 
	double mean_f[CHROM_NUMBER], mean_f_ir[CHROM_NUMBER], mean_f_gc[CHROM_NUMBER];

	double sum_m[CHROM_NUMBER], sum_m_ir[CHROM_NUMBER], sum_m_gc[CHROM_NUMBER];
	double var_m[CHROM_NUMBER], var_m_ir[CHROM_NUMBER], var_m_gc[CHROM_NUMBER];
	double sd_m[CHROM_NUMBER], sd_m_ir[CHROM_NUMBER], sd_m_gc[CHROM_NUMBER]; 
	double mean_m[CHROM_NUMBER], mean_m_ir[CHROM_NUMBER], mean_m_gc[CHROM_NUMBER];


	double z_score[MAX_SAMPLE][CHROM_NUMBER];
	double z_score_ir[MAX_SAMPLE][CHROM_NUMBER];
	double z_score_gc[MAX_SAMPLE][CHROM_NUMBER];

	int internal_ref[CHROM_NUMBER] = 
	{
			0,	0,	0, 	0,	0,	0,	0,
		0,	0,	0,	0,	0,	3,	0,	0,
		0,	0,	7,	0,	0,	13,	0,	6,
		6	
	};	

//	some normal value
	double z_score_cut = 3.0;
	double chr_y_per_cut = 0.0197 ;

	double cff_hundred = 0.309144;
	double cff_hundred_gc = 0.320532;
	double cff_zero = 0.017513;
	double cff_zero_gc = 0.017704;

	for (i=0;i<=total_barcode_index;i++) 
		for (j=0;j<CHROM_NUMBER;j++) {
			barcode_ch_normal[i][j] = barcode_ch_normal_gc_corr[i][j] = 0;
			for (k=0;k<BIN_NUMBER;k++) {
				if (!bin_valid_flag[j][k])
					continue;
				barcode_ch_normal[i][j] += barcode_ch_bin_normal[i][j][k];
				barcode_ch_normal_gc_corr[i][j] += barcode_ch_bin_normal_gc_corr[i][j][k];
			}
			if (j == 23 && (barcode_ch_normal_gc_corr[i][j] > (chr_y_per_cut * NORMAL_VAL / 100))) 
				sex[i] = 1;

			fprintf(stdout, "%d\t%d\t%lf\t%lf\n",i,j,barcode_ch_normal[i][j],barcode_ch_normal_gc_corr[i][j]);
		}

	for (i=0;i<total_barcode_index;i++)
		if (sex[i])
			male_num ++;
		else 
			female_num ++;
/*
	for (i=0;i<=total_barcode_index;i++) 
		for (j=AUTOSOME_NUMBER;j<CHROM_NUMBER;j++) {			
			barcode_ch_normal[i][j] = barcode_ch_normal_gc_corr[i][j] = 0;
			for (k=0;k<BIN_NUMBER;k++) {
				if (!bin_valid_flag[j][k])
					continue;
				barcode_ch_normal[i][j] += barcode_ch_bin_normal[i][j][k];
				arcode_ch_normal_gc_corr[i][j] += barcode_ch_bin_normal_gc_corr[i][j][k];
			}
		}
*/
	if (1) { // without experience file
		for (j=0;j<AUTOSOME_NUMBER;j++) {
			sum[j] = var[j] = sum_gc[j] = var_gc[j] = sum_ir[j] = var_ir[j] = 0;

			for (i=0;i<=total_barcode_index;i++) {
				sum[j] += barcode_ch_normal[i][j];
				var[j] += barcode_ch_normal[i][j] * barcode_ch_normal[i][j];
				sum_gc[j] += barcode_ch_normal_gc_corr[i][j];
				var_gc[j] += barcode_ch_normal_gc_corr[i][j] * barcode_ch_normal_gc_corr[i][j];
				if (internal_ref[j]) {
					sum_ir[j] += barcode_ch_normal[i][j] / barcode_ch_normal[i][internal_ref[j]];
					var_ir[j] += barcode_ch_normal[i][j] * barcode_ch_normal[i][j] 
						/ barcode_ch_normal[i][internal_ref[j]] / barcode_ch_normal[i][internal_ref[j]];
				}
			}
		
			mean[j] = sum[j] / (total_barcode_index+1);
			mean_gc[j] = sum_gc[j] / (total_barcode_index+1);
			mean_ir[j] = sum_ir[j] / (total_barcode_index+1);
			sd[j] = sqrt(var[j] / (total_barcode_index+1) - mean[j] * mean[j]);
			sd_gc[j] = sqrt(var_gc[j] / (total_barcode_index+1) - mean_gc[j] * mean_gc[j]);
			sd_ir[j] = sqrt(var_ir[j] / (total_barcode_index+1) - mean_ir[j] * mean_ir[j]);
		}

		for (j=AUTOSOME_NUMBER;j<CHROM_NUMBER;j++) {
			
			sum_f[j] = var_f[j] = sum_f_gc[j] = var_f_gc[j] = sum_f_ir[j] = var_f_ir[j] = 0;
			sum_m[j] = var_m[j] = sum_m_gc[j] = var_m_gc[j] = sum_m_ir[j] = var_m_ir[j] = 0;

			for (i=0;i<=total_barcode_index;i++) {
				if (sex[i] ) {
					sum_m[j] += barcode_ch_normal[i][j];
					var_m[j] += barcode_ch_normal[i][j] * barcode_ch_normal[i][j];
					sum_m_gc[j] += barcode_ch_normal_gc_corr[i][j];
					var_m_gc[j] += barcode_ch_normal_gc_corr[i][j] * barcode_ch_normal_gc_corr[i][j];
					sum_m_ir[j] += barcode_ch_normal[i][j] / barcode_ch_normal[i][internal_ref[j]];
					var_m_ir[j] += barcode_ch_normal[i][j] * barcode_ch_normal[i][j] 
						/ barcode_ch_normal[i][internal_ref[j]] / barcode_ch_normal[i][internal_ref[j]];
				} else {
					sum_f[j] += barcode_ch_normal[i][j];
					var_f[j] += barcode_ch_normal[i][j] * barcode_ch_normal[i][j];
					sum_f_gc[j] += barcode_ch_normal_gc_corr[i][j];
					var_f_gc[j] += barcode_ch_normal_gc_corr[i][j] * barcode_ch_normal_gc_corr[i][j];
					sum_f_ir[j] += barcode_ch_normal[i][j] / barcode_ch_normal[i][internal_ref[j]];
					var_f_ir[j] += barcode_ch_normal[i][j] * barcode_ch_normal[i][j] 
						/ barcode_ch_normal[i][internal_ref[j]] / barcode_ch_normal[i][internal_ref[j]];
				}
			}

			mean_f[j] = sum_f[j] / female_num;
			mean_f_gc[j] = sum_f_gc[j] / female_num;
			mean_f_ir[j] = sum_f_ir[j] / female_num;
			sd_f[j] = sqrt(var_f[j] / female_num - mean_f[j] * mean_f[j]);
			sd_f_gc[j] = sqrt(var_f_gc[j] / female_num - mean_f_gc[j] * mean_f_gc[j]);
			sd_f_ir[j] = sqrt(var_f_ir[j] / female_num - mean_f_ir[j] * mean_f_ir[j]);

			mean_m[j] = sum_m[j] / male_num;
			mean_m_gc[j] = sum_m_gc[j] / male_num;
			mean_m_ir[j] = sum_m_ir[j] / male_num;
			sd_m[j] = sqrt(var_m[j] / male_num - mean_m[j] * mean_m[j]);
			sd_m_gc[j] = sqrt(var_m_gc[j] / male_num - mean_m_gc[j] * mean_m_gc[j]);
			sd_m_ir[j] = sqrt(var_m_ir[j] / male_num - mean_m_ir[j] * mean_m_ir[j]);
		}

		// filt abnormal samples	  

		int sample_num_valid;

		for (k=0;k<loop;k++) 
		{
	
		fprintf(stderr,"\n");		

		for (j=0;j<AUTOSOME_NUMBER;j++) {
			sum[j] = var[j] = sum_gc[j] = var_gc[j] = sum_ir[j] = var_ir[j] = sample_num_valid = 0;
			for (i=0;i<=total_barcode_index;i++) {
				if (barcode_ch_normal[i][j] < mean[j] - 3 * sd[j] || barcode_ch_normal[i][j] > mean[j] + 3 * sd[j]) {
					fprintf(stderr,"abnormal sample %d :%d chr %d value %lf mean %lf sd %lf\n",k, i, j+1, barcode_ch_normal[i][j], mean[j], sd[j]);
					continue;
				}
				sample_num_valid++;
				sum[j] += barcode_ch_normal[i][j];
				var[j] += barcode_ch_normal[i][j] * barcode_ch_normal[i][j];
				sum_gc[j] += barcode_ch_normal_gc_corr[i][j];
				var_gc[j] += barcode_ch_normal_gc_corr[i][j] * barcode_ch_normal_gc_corr[i][j];
				if (internal_ref[j]) {
					sum_ir[j] += barcode_ch_normal[i][j] / barcode_ch_normal[i][internal_ref[j]];
					var_ir[j] += barcode_ch_normal[i][j] * barcode_ch_normal[i][j] 
						/ barcode_ch_normal[i][internal_ref[j]] / barcode_ch_normal[i][internal_ref[j]];
				}
			}
			mean[j] = sum[j] / sample_num_valid;
			mean_gc[j] = sum_gc[j] / sample_num_valid;
			mean_ir[j] = sum_ir[j] / sample_num_valid;
			sd[j] = sqrt(var[j] / sample_num_valid - mean[j] * mean[j]);
			sd_gc[j] = sqrt(var_gc[j] / sample_num_valid - mean_gc[j] * mean_gc[j]);
			sd_ir[j] = sqrt(var_ir[j] / sample_num_valid - mean_ir[j] * mean_ir[j]);
		}

		for (j=AUTOSOME_NUMBER;j<CHROM_NUMBER;j++) {
			
			sum_f[j] = var_f[j] = sum_f_gc[j] = var_f_gc[j] = sum_f_ir[j] = var_f_ir[j] = 0;
			sum_m[j] = var_m[j] = sum_m_gc[j] = var_m_gc[j] = sum_m_ir[j] = var_m_ir[j] = 0;
			male_num = female_num = 0;

			for (i=0;i<=total_barcode_index;i++) {
				if (sex[i] ) {
					if (barcode_ch_normal[i][j] < mean_m[j] - 3 * sd_m[j] || barcode_ch_normal[i][j] > mean_m[j] + 3 * sd_m[j]) {
						fprintf(stderr,"abnormal sample %d :%d chr %d value %lf mean %lf sd %lf\n",k, i, j+1, barcode_ch_normal[i][j], mean_m[j], sd_m[j]);
						continue;
					}
					male_num ++;
					sum_m[j] += barcode_ch_normal[i][j];
					var_m[j] += barcode_ch_normal[i][j] * barcode_ch_normal[i][j];
					sum_m_gc[j] += barcode_ch_normal_gc_corr[i][j];
					var_m_gc[j] += barcode_ch_normal_gc_corr[i][j] * barcode_ch_normal_gc_corr[i][j];
					sum_m_ir[j] += barcode_ch_normal[i][j] / barcode_ch_normal[i][internal_ref[j]];
					var_m_ir[j] += barcode_ch_normal[i][j] * barcode_ch_normal[i][j] 
						/ barcode_ch_normal[i][internal_ref[j]] / barcode_ch_normal[i][internal_ref[j]];
				} else {
					if (barcode_ch_normal[i][j] < mean_f[j] - 3 * sd_f[j] || barcode_ch_normal[i][j] > mean_f[j] + 3 * sd_f[j]) {
						fprintf(stderr,"abnormal sample %d :%d chr %d value %lf mean %lf sd %lf\n",k, i, j+1, barcode_ch_normal[i][j], mean_f[j], sd_f[j]);
						continue;
					}
					female_num ++;
					sum_f[j] += barcode_ch_normal[i][j];
					var_f[j] += barcode_ch_normal[i][j] * barcode_ch_normal[i][j];
					sum_f_gc[j] += barcode_ch_normal_gc_corr[i][j];
					var_f_gc[j] += barcode_ch_normal_gc_corr[i][j] * barcode_ch_normal_gc_corr[i][j];
					sum_f_ir[j] += barcode_ch_normal[i][j] / barcode_ch_normal[i][internal_ref[j]];
					var_f_ir[j] += barcode_ch_normal[i][j] * barcode_ch_normal[i][j] 
						/ barcode_ch_normal[i][internal_ref[j]] / barcode_ch_normal[i][internal_ref[j]];
				}
			}

			mean_f[j] = sum_f[j] / female_num;
			mean_f_gc[j] = sum_f_gc[j] / female_num;
			mean_f_ir[j] = sum_f_ir[j] / female_num;
			sd_f[j] = sqrt(var_f[j] / female_num - mean_f[j] * mean_f[j]);
			sd_f_gc[j] = sqrt(var_f_gc[j] / female_num - mean_f_gc[j] * mean_f_gc[j]);
			sd_f_ir[j] = sqrt(var_f_ir[j] / female_num - mean_f_ir[j] * mean_f_ir[j]);

			mean_m[j] = sum_m[j] / male_num;
			mean_m_gc[j] = sum_m_gc[j] / male_num;
			mean_m_ir[j] = sum_m_ir[j] / male_num;
			sd_m[j] = sqrt(var_m[j] / male_num - mean_m[j] * mean_m[j]);
			sd_m_gc[j] = sqrt(var_m_gc[j] / male_num - mean_m_gc[j] * mean_m_gc[j]);
			sd_m_ir[j] = sqrt(var_m_ir[j] / male_num - mean_m_ir[j] * mean_m_ir[j]);
		}

		}

		FILE *fp_curr;
		char name_curr[BUFFER_LENGTH];

		sprintf(name_curr,"%s/current_samples.stats",out_prefix);
		
		if ((fp_curr = fopen(name_curr, "w")) == NULL) {
			printf("Cannot write %s. Now exit to system...\n",name_curr);
			exit(-1);
		}

		fprintf(fp_curr,"chr\tmethod\tmean\tsd\n");
		
		int chr4test[5] = {13,18,21,23,24};
		
		for(i=0;i<3;i++) {	
			fprintf(fp_curr,"chr%d\tor\t%lf\t%lf\n",chr4test[i],mean[chr4test[i]-1],sd[chr4test[i]-1]);
			fprintf(fp_curr,"chr%d\tir\t%lf\t%lf\n",chr4test[i],mean_ir[chr4test[i]-1],sd_ir[chr4test[i]-1]);	
			fprintf(fp_curr,"chr%d\tgc\t%lf\t%lf\n",chr4test[i],mean_gc[chr4test[i]-1],sd_gc[chr4test[i]-1]);
		}

		for(i=3;i<5;i++) {
			fprintf(fp_curr,"chr%d\tmor\t%lf\t%lf\n",chr4test[i],mean_m[chr4test[i]-1],sd_m[chr4test[i]-1]);
			fprintf(fp_curr,"chr%d\tmir\t%lf\t%lf\n",chr4test[i],mean_m_ir[chr4test[i]-1],sd_m_ir[chr4test[i]-1]);	
			fprintf(fp_curr,"chr%d\tmgc\t%lf\t%lf\n",chr4test[i],mean_m_gc[chr4test[i]-1],sd_m_gc[chr4test[i]-1]);
			fprintf(fp_curr,"chr%d\tfor\t%lf\t%lf\n",chr4test[i],mean_f[chr4test[i]-1],sd_f[chr4test[i]-1]);
			fprintf(fp_curr,"chr%d\tfir\t%lf\t%lf\n",chr4test[i],mean_f_ir[chr4test[i]-1],sd_f_ir[chr4test[i]-1]);	
			fprintf(fp_curr,"chr%d\tfgc\t%lf\t%lf\n",chr4test[i],mean_f_gc[chr4test[i]-1],sd_f_gc[chr4test[i]-1]);
		}
		fclose(fp_curr);		
	} else {		  	

	}
	// load experience value

	if (exp_file[0]) {

	FILE *fp_exp;
	char line[BUFFER_LENGTH];

	if ((fp_exp = fopen(exp_file, "r")) == NULL) {
		printf("Cannot read %s. Now exit to system...\n",exp_file);
		exit(-1);
	}

	int chr_exp;	
	char flag[2]; 
	double mean_exp, sd_exp;

	while(fgets(line,sizeof(line),fp_exp)!=NULL) {
	 	sscanf(line,"chr%d %s %lf %lf", &chr_exp, flag, &mean_exp, &sd_exp);
		if (flag[0] == 'o') {
			mean[chr_exp-1] = mean_exp;
			sd[chr_exp-1] = sd_exp;
		} else if (flag[0] == 'i') {
			mean_ir[chr_exp-1] = mean_exp;
			sd_ir[chr_exp-1] = sd_exp;
		} else if (flag[0] == 'g') {
			mean_gc[chr_exp-1] = mean_exp;
			sd_gc[chr_exp-1] = sd_exp;
		} else if (flag[0] == 'm') {
			if (flag[1] == 'o') {
				mean_m[chr_exp-1] = mean_exp;
				sd_m[chr_exp-1] = sd_exp;
			} else if (flag[1] == 'i') {
				mean_m_ir[chr_exp-1] = mean_exp;
				sd_m_ir[chr_exp-1] = sd_exp;
			} else if (flag[1] == 'g') {
				mean_m_gc[chr_exp-1] = mean_exp;
				sd_m_gc[chr_exp-1] = sd_exp;
			}
		} else if (flag[0] == 'f') {
			if (flag[1] == 'o') {
				mean_f[chr_exp-1] = mean_exp;
				sd_f[chr_exp-1] = sd_exp;
			} else if (flag[1] == 'i') {
				mean_f_ir[chr_exp-1] = mean_exp;
				sd_f_ir[chr_exp-1] = sd_exp;
			} else if (flag[1] == 'g') {
				mean_f_gc[chr_exp-1] = mean_exp;
				sd_f_gc[chr_exp-1] = sd_exp;
			}
		}
	}

	fclose(fp_exp);
  	}	
		
	// calcuate z-score	
	
	for (i=0;i<=total_barcode_index;i++) {
		for (j=0;j<AUTOSOME_NUMBER;j++) {
			z_score[i][j] = (barcode_ch_normal[i][j] - mean[j]) / sd[j];
			z_score_gc[i][j] = (barcode_ch_normal_gc_corr[i][j] - mean_gc[j]) / sd_gc[j];
			if (internal_ref[j]) 
			z_score_ir[i][j] = (barcode_ch_normal[i][j] / barcode_ch_normal[i][internal_ref[j]] - mean_ir[j]) / sd_ir[j];
		}

		for (j=AUTOSOME_NUMBER;j<CHROM_NUMBER;j++) {
			if (sex[i]) { 
				z_score[i][j] = (barcode_ch_normal[i][j] - mean_m[j]) / sd_m[j];
				z_score_gc[i][j] = (barcode_ch_normal_gc_corr[i][j] - mean_m_gc[j]) / sd_m_gc[j];
				z_score_ir[i][j] = (barcode_ch_normal[i][j] / barcode_ch_normal[i][internal_ref[j]] - mean_m_ir[j]) / sd_m_ir[j];
			} else { 
				z_score[i][j] = (barcode_ch_normal[i][j] - mean_f[j]) / sd_f[j];
				z_score_gc[i][j] = (barcode_ch_normal_gc_corr[i][j] - mean_f_gc[j]) / sd_f_gc[j];
				z_score_ir[i][j] = (barcode_ch_normal[i][j] / barcode_ch_normal[i][internal_ref[j]] - mean_f_ir[j]) / sd_f_ir[j];
			}
		}
	}

	char name[BUFFER_LENGTH];
	FILE *fp;

	sprintf(name,"%s/allchrom.csv",out_prefix);

	if ((fp = fopen(name, "w")) == NULL) {
		printf("Cannot write to %s. Now exit to system...\n", name);
		exit(-1);
	}
	
	fprintf(fp, "sample,lib.name");
	for (j=0;j<CHROM_NUMBER;j++) 
		fprintf(fp, ",chr%d_per,z_score,z_score_gc",j+1);
	fprintf(fp, "\n");

	for (i=0;i<=total_barcode_index;i++) {
		fprintf(fp, "%d,%s", i, &lib_name[i*NAME_LENGTH]); 
		for (j=0;j<CHROM_NUMBER;j++)  
			fprintf(fp, ",%.6f,%.4f,%.4f", barcode_ch_normal[i][j] / NORMAL_VAL * 100, z_score[i][j], z_score_gc[i][j]);
		fprintf(fp, "\n");
	}	

	fclose(fp);


	for (i=0;i<=total_barcode_index;i++) {
		barcode_conclusion[i][0] = total_barcode_result[i];
		barcode_conclusion[i][1] = sum_reads[i];
		barcode_conclusion[i][2] = (barcode_ch_normal[i][23] / NORMAL_VAL * 100 > chr_y_per_cut ) ? 1 : 0;
		if (z_score[i][20]>z_score_cut&&z_score_gc[i][20]>z_score_cut) {
			barcode_conclusion[i][3] = 1;
		} else if (z_score[i][20]>z_score_cut||z_score_gc[i][20]>z_score_cut) {
			barcode_conclusion[i][3] = 2;
		} else {
			barcode_conclusion[i][3] = 0;
		}
		if (z_score_gc[i][12]>z_score_cut&&z_score_ir[i][12]>z_score_cut) {
			barcode_conclusion[i][4] = 1;
		} else if (z_score_gc[i][12]>z_score_cut||z_score_ir[i][12]>z_score_cut) {
			barcode_conclusion[i][4] = 2;
		} else {
			barcode_conclusion[i][4] = 0;
		}
		if (z_score_gc[i][17]>z_score_cut&&z_score_ir[i][17]>z_score_cut) {
			barcode_conclusion[i][5] = 1;
		} else if (z_score_gc[i][17]>z_score_cut||z_score_ir[i][17]>z_score_cut) {
			barcode_conclusion[i][5] = 2;
		} else {
			barcode_conclusion[i][5] = 0;
		}
		barcode_conclusion[i][6] = dup_reads[i];
		barcode_conclusion[i][7] = gc_bases[i];
		barcode_conclusion[i][8] = q20_bases[i];
		barcode_conclusion[i][9] = q30_bases[i];
		
		cff_per_byX[i] = (barcode_ch_normal[i][22] / NORMAL_VAL * 100 - cff_zero) * 100 / (cff_hundred - cff_zero);
		cff_per_gc_byX[i] = (barcode_ch_normal_gc_corr[i][22] / NORMAL_VAL * 100 - cff_zero_gc) * 100 / (cff_hundred_gc - cff_zero_gc);

		if (barcode_conclusion[i][2]) {
			cff_per[i] = (barcode_ch_normal[i][23] / NORMAL_VAL * 100 - cff_zero) * 100 / (cff_hundred - cff_zero);
			cff_per_gc[i] = (barcode_ch_normal_gc_corr[i][23] / NORMAL_VAL * 100 - cff_zero_gc) * 100 / (cff_hundred_gc - cff_zero_gc);
		}
	}

	// output_csv
	

	sprintf(name,"%s/results.csv",out_prefix);
	
	if ((fp = fopen(name, "w")) == NULL) {
		printf("Cannot write to %s. Now exit to system...\n", name);
		exit(-1);
	}
	
	fprintf(fp,"Sample,lib.name,duplication.rate,gc.rate,q20.rate,q30.rate,identical.rate,identical.num,chr21.risk,chr21.z.score,chr21.z.score.ir,chr21.z.score.gc,chr13.risk,chr13.z.score,chr13.z.score.ir,chr13.z.score.gc,chr18.risk,chr18.z.score,chr18.z.score.ir,chr18.z.score.gc,chrX.percent,cff.per.byX,chrX.percent.gc,cff.per.gc.byX,chrY.percent,cff.per,chrY.percent.gc,cff.per.gc,sex\n");

	for (i=0;i<=total_barcode_index;i++) {
		fprintf(fp,"%d,%s,%.2f,%.2f,%.2f,%.2f,%.2f,%d,%s,%.2f,%.2f,%.2f,%s,%.2f,%.2f,%.2f,%s,%.2f,%.2f,%.2f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%s\n",
				i,&lib_name[i*NAME_LENGTH],barcode_conclusion[i][6]*100.0/barcode_conclusion[i][0],barcode_conclusion[i][7]*2.0/barcode_conclusion[i][0],
				barcode_conclusion[i][8]*2.0/barcode_conclusion[i][0], barcode_conclusion[i][9]*2.0/barcode_conclusion[i][0],
				barcode_conclusion[i][1]*100.0/barcode_conclusion[i][0],barcode_conclusion[i][1],
				barcode_conclusion[i][3]==1?"high":(barcode_conclusion[i][3]==2?"questionable":"low"),z_score[i][20],z_score_ir[i][20],z_score_gc[i][20],
				barcode_conclusion[i][4]==1?"high":(barcode_conclusion[i][4]==2?"questionable":"low"),z_score[i][12],z_score_ir[i][12],z_score_gc[i][12],
				barcode_conclusion[i][5]==1?"high":(barcode_conclusion[i][5]==2?"questionable":"low"),z_score[i][17],z_score_ir[i][17],z_score_gc[i][17],
				barcode_ch_normal[i][22] / NORMAL_VAL * 100, cff_per_byX[i],
				barcode_ch_normal_gc_corr[i][22] / NORMAL_VAL * 100, cff_per_gc_byX[i],
				barcode_ch_normal[i][23] / NORMAL_VAL * 100, cff_per[i],
				barcode_ch_normal_gc_corr[i][23] / NORMAL_VAL * 100 ,
				cff_per_gc[i], (barcode_conclusion[i][2]?"male":"female")); 	
	}

	fclose(fp);

	sprintf(name,"%s/pre_report.csv",out_prefix);
	
	if ((fp = fopen(name, "w")) == NULL) {
		printf("Cannot write to %s. Now exit to system...\n", name);
		exit(-1);
	}
	
	for (i=0;i<=total_barcode_index;i++) {

	fprintf(fp,"%s,%.4f,%s,%.4f,%s,%.4f,%s\n",
		&lib_name[i*NAME_LENGTH], (z_score_ir[i][12]+z_score_gc[i][12])/2, (z_score_ir[i][12]+z_score_gc[i][12])/2>3?"+":"-",
		(z_score_ir[i][17]+z_score_gc[i][17])/2, (z_score_ir[i][17]+z_score_gc[i][17])/2>3?"+":"-",
		(z_score[i][20]+z_score_gc[i][20])/2,((z_score[i][20]+z_score_gc[i][20])/2>3)?"+":"-");

	}

	fclose(fp);

	// plotting
	
	sprintf(name,"%s/plot_stats.R",out_prefix);

	if ((fp = fopen(name, "w")) == NULL) {
		printf("Cannot write to %s. Now exit to system...\n", name);
		exit(-1);
	}
		
	int chrom_p[3] = {21,13,18};

	fprintf(fp,
"#! /usr/bin/env Rscript\n\
library(ggplot2)\n\
a<-read.csv(\"%s/results.csv\")\n\n",out_prefix); 

	for (i=0;i<3;i++) {
		fprintf(fp,
"png(type=\"cairo\",\"%s/chr%d.z_score.png\")\n\
p<-ggplot(a,aes(x=\"Original\",y=chr%d.z.score))+geom_point(aes(color=chr%d.risk),shape=2,size=2)+geom_hline(color=\"darkblue\",linetype=3,yintercept=%.2f)+labs(list(title=\"chr%d z-score distribution\",x=\"Samples\",y=\"Z-score\"))\n\
p\n\
dev.off()\n\n",out_prefix,chrom_p[i],chrom_p[i],chrom_p[i],z_score_cut,chrom_p[i]);
		fprintf(fp,
"png(type=\"cairo\",\"%s/chr%d.z_score_ir.png\")\n\
p<-ggplot(a,aes(x=\"Internal Reference corrected\",y=chr%d.z.score.ir))+geom_point(aes(color=chr%d.risk),shape=2,size=2)+geom_hline(color=\"darkblue\",linetype=3,yintercept=%.2f)+labs(list(title=\"chr%d z-score distribution\",x=\"Samples\",y=\"Z-score-ir\"))\n\
p\n\
dev.off()\n\n",out_prefix,chrom_p[i],chrom_p[i],chrom_p[i],z_score_cut,chrom_p[i]);
		fprintf(fp,
"png(type=\"cairo\",\"%s/chr%d.z_score_gc.png\")\n\
p<-ggplot(a,aes(x=\"GC corrected\",y=chr%d.z.score.gc))+geom_point(aes(color=chr%d.risk),shape=2,size=2)+geom_hline(color=\"darkblue\",linetype=3,yintercept=%.2f)+labs(list(title=\"chr%d z-score distribution\",x=\"Samples\",y=\"Z-score-gc\"))\n\
p\n\
dev.off()\n\n",out_prefix,chrom_p[i],chrom_p[i],chrom_p[i],z_score_cut,chrom_p[i]);
		fprintf(fp,
"png(type=\"cairo\",\"%s/chr%d.z_score_all.png\")\n\
p<-ggplot(a,aes(x=\"Original\",y=chr%d.z.score))+geom_point(aes(color=chr%d.risk),shape=2,size=2)+geom_point(aes(x=\"Interal Reference corrected\",y=chr%d.z.score.ir,color=chr%d.risk),shape=2,size=2)+geom_point(aes(x=\"GC corrected\",y=chr%d.z.score.gc,color=chr%d.risk),shape=2,size=2)+geom_hline(color=\"darkblue\",linetype=3,yintercept=%.2f)+labs(list(title=\"chr%d z-score distribution\",x=\"Samples\",y=\"Z-score\"))\n\
p\n\
dev.off()\n\n",out_prefix,chrom_p[i],chrom_p[i],chrom_p[i],chrom_p[i],chrom_p[i],chrom_p[i],chrom_p[i],z_score_cut,chrom_p[i]);

	}

	fprintf(fp,
"png(type=\"cairo\",\"%s/chrY.per.png\")\n\
p<-ggplot(a,aes(x=\"Original\",y=chrY.percent))+geom_point(aes(color=sex),shape=2,size=2)+geom_hline(color=\"darkblue\",linetype=3,yintercept=%.4f)+labs(list(title=\"chrY mapping percent distribution\",x=\"Samples\",y=\"chrY %%\"))\n\
p\n\
dev.off()\n\n",out_prefix,chr_y_per_cut);

	fprintf(fp,
"png(type=\"cairo\",\"%s/chrY.per.gc.png\")\n\
p<-ggplot(a,aes(x=\"Original\",y=chrY.percent.gc))+geom_point(aes(color=sex),shape=2,size=2)+geom_hline(color=\"darkblue\",linetype=3,yintercept=%.4f)+labs(list(title=\"chrY mapping percent distribution\",x=\"Samples\",y=\"chrY.gc %%\"))\n\
p\n\
dev.off()\n\n",out_prefix,chr_y_per_cut);
	fclose(fp);
	
	char cmd[BUFFER_LENGTH];
	sprintf(cmd,"chmod 755 %s &&Rscript %s", name, name);
	
	system(cmd);

}

static void display_usage()
{
	printf("\n \033[33;1m%s\033[0;m\tversion: \033[36;1m%.1f\033[0;m\n", "uacd_stats", 3.0);
	printf(" Last modified on : 2014.01.08\n");
	printf("\n Author: \033[34;1m%s\033[0;m\n","Ruan Hang <rogerholmes1036@gmail.com>");
	printf("\n Usage:\n");
	printf("	-i input_statistics\n");
	printf("	-o output_prefix\n");
	printf("	[-e exp_file]\n");
	printf("	[-f bin_valid_flag_file]\n");
	printf("	[-l loop_value | default:5]\n");
	printf("	[-v output gc_graph | default:no ]\n");
	printf("\n");
}

int main(int argc, char *argv[]) {
	int paras;
	extern char *optarg;
	char in_stats[FILE_NAME_LENGTH];
	char temp[FILE_NAME_LENGTH];	

	if (argc==1) {
		display_usage();
		return 0;
	}

	gc_verbose = 0;
	loop = 5;

	while((paras=getopt(argc,argv,"i:o:e:c:f:l:v"))!=EOF) {
		switch(paras) {
		case 'i':
			sscanf(optarg,"%s",in_stats);
			continue;
		case 'o':
			sscanf(optarg,"%s",out_prefix);
			continue;
		case 'e':
			sscanf(optarg,"%s",exp_file);
			continue;
		case 'c':
			sscanf(optarg,"%s",corr_file);
			continue;
		case 'f':
			sscanf(optarg,"%s",flag_file);
			continue;
		case 'v':
			gc_verbose = 1;
			continue;
		case 'l':
			sscanf(optarg,"%s",temp);
			loop = atoi(temp);	
			continue;
		default :
			display_usage();
			exit(-1);
		}
	}

	FILE *fp;
	char line[BUFFER_LENGTH];
	
	if ((fp=fopen(in_stats,"r"))==NULL) {
		printf("Can not open %s\n",in_stats);
		return -1;
	}

	total_barcode_index=-1;	
	int sample_id,total_reads,mapped_reads,duplication_reads,bin_id,bin_reads;	
	int gc,q20,q30;
	char name[NAME_LENGTH];
	
	lib_name = (char *)calloc(1,MAX_SAMPLE*NAME_LENGTH*sizeof(char));
	bin_id = 0;

	while(fgets(line,sizeof(line),fp)!=NULL) {
		if (line[0] == '>') {
			total_barcode_index++;
			sscanf(line,">%d_%[^_]_%d_%d_%d_%d_%d_%d",&sample_id,name,&total_reads,&mapped_reads,&duplication_reads,&gc,&q20,&q30);
			total_barcode_result[total_barcode_index]=total_reads;
			sum_reads[total_barcode_index]=mapped_reads; 	
			dup_reads[total_barcode_index]=duplication_reads;
			gc_bases[total_barcode_index]=gc;
			q20_bases[total_barcode_index]=q20;
			q30_bases[total_barcode_index]=q30;
			strcpy(&lib_name[total_barcode_index*NAME_LENGTH],name);
			bin_id=0;
		} else {
			sscanf(line,"%d",&bin_reads);
			int current_chrom = bin_id / BIN_NUMBER ;
			int current_bin = bin_id  % BIN_NUMBER ;
			barcode_ch_bin_result[total_barcode_index][current_chrom][current_bin] = bin_reads;
			bin_id ++;
		}
	}

	fclose(fp);

	char cmd[BUFFER_LENGTH];
	sprintf(cmd, "mkdir %s", out_prefix);	
	system(cmd);
	
	if (gc_verbose) {
		sprintf(cmd, "mkdir %s/gc", out_prefix);
		system(cmd);
	}

	stats();

	free((void*) lib_name);

	return 1;
}
				
