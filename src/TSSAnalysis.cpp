/* 

Program TSSAnalysis.c to:

	1) input gene, location/incidence data for a set of genes denoted in column 1, location (nucleotide site) in second column, and read number for the npop isolates in the remaining columns  ;

	2) calculate various statistics for each gene -- moments within genes; anova between samples.
	
*****NOTE:***** 
	1) Columns are tab delimited.
	2) Default name for input file is "datain.txt".

*/


#define npop		4					/* number of populations for which there are data in the input file */



#include	<stdio.h>
#include 	<math.h>
#include 	<sys/types.h>
#include	<stdlib.h>
#include	<time.h>


/* Point to the input and output files. */

FILE *instream;
FILE *outstream;


void main()

{

int gene;														/* annotated gene number */
int oldgene;

int site;														/* nucleotide site preceding the first nucleotide of annotated translation-start codon, e.g., the immediately prior site is -1; the first site after ATG codon is +3 */

int ig;															/* counter for samples */

double fileline;												/* line number in data input file being read */

int n[100];														/* the read numbers at the site for the different samples */
double nd[100];													/* shift to double precision from integer counts */

double sited;													/* shift to double precision for site number */

double totsum;													/* total number of reads for the gene */
double sum[100];												/* total number of reads for the different samples for the gene */

double totsum1, totsum2, totsum3;								/* sums for the first four moments around zero for the total sample for a gene */
double sum1[100], sum2[100], sum3[100];							/* sums for the first four moments around zero for the different samples for a gene */

double tot1, tot2, tot3;										/* first four moments around zero for the total sample for a gene */
double samp1[100], samp2[100], samp3[100];						/* first four moments around zero for the different samples for a gene */

double var, sd, skw, k3;										/* variance and skewness statistics for total sample of a gene */
double sampvar[100], sampsd[100], sampskw[100], sampk3[100];



/* Open the output and input files. */

fopen_s(&outstream,"dataout.txt", "w");
fopen_s(&instream, "datain.txt", "r");



/* Set the initial conditions and counters. */

totsum = 0.0;								
totsum1 = 0.0;
totsum2 = 0.0;
totsum3 = 0.0;

for (ig = 1; ig <= npop; ++ig) {
	sum[ig] = 0.0;							
	sum1[ig] = 0.0;
	sum2[ig] = 0.0;
	sum3[ig] = 0.0; }

fileline = 0.0;



/* Start inputing and analyzing the data line by line. */

while (  fscanf(instream, "%i \t%i \t%i \t%i \t%i \t%i", &gene, &site, &n[1], &n[2], &n[3], &n[4])  != EOF) { 

	fileline = fileline + 1.0;
	if (fileline == 1.0) {
		oldgene = gene; }


	 for (ig = 1; ig <= npop; ++ig) {								/* Switch to double precision entries. */
		 nd[ig] = ((double) n[ig]);  }
	 sited = ((double) site);



	 if (gene == oldgene) {											/* If this line of data is a continuation of the same gene, add to the summation terms. */

		 for (ig = 1; ig <= npop; ++ig) {

			 totsum = totsum + nd[ig];								/* Total sums of counts and moments over all samples. */
			 totsum1 = totsum1 + (nd[ig] * sited);
			 totsum2 = totsum2 + (nd[ig] * pow(sited, 2.0));
			 totsum3 = totsum3 + (nd[ig] * pow(sited, 3.0));

			 sum[ig] = sum[ig] + nd[ig];							/* Sums of counts and moments for the gene for sample ig. */
			 sum1[ig] = sum1[ig] + (nd[ig] * sited);
			 sum2[ig] = sum2[ig] + (nd[ig] * pow(sited, 2.0));
			 sum3[ig] = sum3[ig] + (nd[ig] * pow(sited, 3.0));
		 }
	 }


	 else {															/* If line of data is for the next gene, summarize the statistics for the prior gene. */

		 tot1 = totsum1 / totsum;									/* Moments over all samples for the gene. */
		 tot2 = totsum2 / totsum;
		 tot3 = totsum3 / totsum;

		 for (ig = 1; ig <= npop; ++ig) {							/* Moments for each sample for the gene. */
			 samp1[ig] = sum1[ig] / sum[ig];
			 samp2[ig] = sum2[ig] / sum[ig];
			 samp3[ig] = sum3[ig] / sum[ig];  }

																	/* Calculate the variance and skewness statistics. */
																	
		 var = (tot2 - pow(tot1, 2.0)) * (totsum / (totsum - 1.0));
		 if (var <= 0.0) {
			 var = 0.000000001;  }

		 sd = pow(var, 0.5);
		 skw = (tot3 - (3.0 * tot2 * tot1) + (2.0 * pow(tot1, 3.0))) * (   pow(totsum, 2.0) / ((totsum - 1.0) * (totsum - 2.0))   );
		 k3 = skw / pow(var, 1.5);

		 for (ig = 1; ig <= npop; ++ig) {							
			 sampvar[ig] = (samp2[ig] - pow(samp1[ig], 2.0)) * (sum[ig] / (sum[ig] - 1.0));
			 if (sampvar[ig] <= 0.0) {
				 sampvar[ig] = 0.000000001; }

			 sampsd[ig] = pow(sampvar[ig], 0.5);
			 sampskw[ig] = (samp3[ig] - (3.0 * samp2[ig] * samp1[ig]) + (2.0 * pow(samp1[ig], 3.0))) * (pow(sum[ig], 2.0) / ((sum[ig] - 1.0) * (sum[ig] - 2.0)));
			 sampk3[ig] = sampskw[ig] / pow(sampvar[ig], 1.5); 	 }


		printf("%7d ,%s, %12.5f , %12.5f , %12.5f , %10.0f\n", oldgene, "   Total", tot1, sd, k3, totsum);
		fprintf(outstream, "%7d, %s, %12.5f , %12.5f , %12.5f , %10.0f\n", oldgene, "  Total", tot1, sd, k3, totsum); 

		for (ig = 1; ig <= npop; ++ig) {
			printf("%7d, %7d , %12.5f , %12.5f , %12.5f , %10.0f\n", oldgene, ig, samp1[ig], sampsd[ig], sampk3[ig], sum[ig]); 
			fprintf(outstream, "%7d, %6d , %12.5f , %12.5f , %12.5f , %10.0f\n", oldgene, ig, samp1[ig], sampsd[ig], sampk3[ig], sum[ig]);
		}
		printf("\n");
		fprintf(outstream, "\n");


		 for (ig = 1; ig <= npop; ++ig) {

			totsum = totsum + nd[ig];								/* Start the counts and computations for the next gene (current line of data). */
			totsum1 = totsum1 + (nd[ig] * sited);
			totsum2 = totsum2 + (nd[ig] * pow(sited, 2.0));
			totsum3 = totsum3 + (nd[ig] * pow(sited, 3.0));

			sum[ig] = sum[ig] + nd[ig];
			sum1[ig] = sum1[ig] + (nd[ig] * sited);
			sum2[ig] = sum2[ig] + (nd[ig] * pow(sited, 2.0));
			sum3[ig] = sum3[ig] + (nd[ig] * pow(sited, 3.0));  }

		 oldgene = gene;

	 }
 
 }					/* Ends the computations when the end of file is reached. */



/* Do the calculations for the final gene. */

tot1 = totsum1 / totsum;									
tot2 = totsum2 / totsum;
tot3 = totsum3 / totsum;

for (ig = 1; ig <= npop; ++ig) {							
	samp1[ig] = sum1[ig] / sum[ig];
	samp2[ig] = sum2[ig] / sum[ig];
	samp3[ig] = sum3[ig] / sum[ig]; }

var = (tot2 - pow(tot1, 2.0)) * (totsum / (totsum - 1.0));
if (var <= 0.0) {
	var = 0.000000001; }

sd = pow(var, 0.5);
skw = (tot3 - (3.0 * tot2 * tot1) + (2.0 * pow(tot1, 3.0))) * (pow(totsum, 2.0) / ((totsum - 1.0) * (totsum - 2.0)));
k3 = skw / pow(var, 1.5);

for (ig = 1; ig <= npop; ++ig) {
	sampvar[ig] = (samp2[ig] - pow(samp1[ig], 2.0)) * (sum[ig] / (sum[ig] - 1.0));
	if (sampvar[ig] <= 0.0) {
		sampvar[ig] = 0.000000001; 	}

	sampsd[ig] = pow(sampvar[ig], 0.5);
	sampskw[ig] = (samp3[ig] - (3.0 * samp2[ig] * samp1[ig]) + (2.0 * pow(samp1[ig], 3.0))) * (pow(sum[ig], 2.0) / ((sum[ig] - 1.0) * (sum[ig] - 2.0)));
	sampk3[ig] = sampskw[ig] / pow(sampvar[ig], 1.5); }


printf("%7d ,%s, %12.5f , %12.5f , %12.5f , %10.0f\n", oldgene, "   Total", tot1, sd, k3, totsum);
fprintf(outstream, "%7d, %s, %12.5f , %12.5f , %12.5f , %10.0f\n", oldgene, "  Total", tot1, sd, k3, totsum);

for (ig = 1; ig <= npop; ++ig) {
	printf("%7d, %7d , %12.5f , %12.5f , %12.5f , %10.0f\n", oldgene, ig, samp1[ig], sampsd[ig], sampk3[ig], sum[ig]);
	fprintf(outstream, "%7d, %6d , %12.5f , %12.5f , %12.5f , %10.0f\n", oldgene, ig, samp1[ig], sampsd[ig], sampk3[ig], sum[ig]);
}


exit(0);

}










/* 

while (fscanf(instream, "%s\t%s\t%i\t%i\t%i\t%i", id1, id2, &n[1], &n[2], &n[3], &n[4]) != EOF) {

	coverage = n[1] + n[2] + n[3] + n[4];



	printf("%s , %s , %c , %c , %6.5f , %6.5f , %6.5f , %5d , %10.5f\n", id1, id2, nuc[major], nuc[minor], pml, (1.0 - pml), eml, coverage, llstat);

	fprintf(outstream, "%s , %s , %c , %c , %6.5f , %6.5f , %6.5f , %5d , %10.5f\n", id1, id2, nuc[major], nuc[minor], pml, (1.0 - pml), eml, coverage, llstat);

	*/


/* printf("\n"); */


/* 


			printf("%5d, %5d, %10.5f, %10.5f, %10.5f, %10.5f\n", gene, oldgene, totsum, totsum1, totsum2, totsum3 ); 
			getchar();

			*/



/* printf("%10.0f, %5d, %5d, %5d, %5d, %5d, %5d\n", fileline, gene, site, n[1], n[2], n[3], n[4]); */

