#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <zlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <assert.h>
#include <math.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int stk_fqchk(char *argv[]) {
	gzFile fp;
	kseq_t *seq;
    //int offset = 33;

	fp = gzopen(argv[1], "r");
	if (fp == 0) {
		fprintf(stderr, "[E::%s] failed to open the input file/stream.\n", __func__);
		return 1;
	}
	seq = kseq_init(fp);
    FILE *fptr = fopen(argv[2],"w");
	fprintf(fptr,"length, quality\n");
    while (kseq_read(seq) >= 0) {
		if (seq->qual.l == 0) continue;
        double avg_qual = 0;
        for (size_t i = 0; i < seq->qual.l; ++i)
            avg_qual += pow(10,(seq->qual.s[i])/-10);
        avg_qual = -10* log10(avg_qual/seq->qual.l);
	    fprintf(fptr,"%zu, %.3f\n", seq->seq.l, avg_qual);
	}
	kseq_destroy(seq);
	gzclose(fp);
    fclose(fptr);
    return 0;
}

int main(int argc, char *argv[]) {
    if(argc<2) {
        printf("Not enough arguments");
        return 0;
    }
    printf("%d",argc);
	return stk_fqchk(argv);
}
