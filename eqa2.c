/* eqaa.c equalise an amino acid sequence .. from the Nuritas test */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>

#define ROWSZ 100
#define GBUF 4
#define n 4
#define FLAG unsigned char /* some time also called a boolean */

typedef struct /*oc_t */
{
	char l; // letter
	int q; //quantity
} oc_t;

typedef struct /*oc2_t */
{
	char l; // letter
	int q; //quantity
    int fap; // first appearance L to R
    int lap; // last appearance L to R
} oc2_t;

typedef struct /*a_t */
{
	int *els; // elements of the array
	int sz; //sze
	int gb; // buf
} a_t;

void usage(char *pname)
{
	printf("Program \"%s\" to find smallest subsequence where amino-acid sequence has equal numbers of all AAs.\n", "eqaa");
	return;
}

int cmpocs2(const void *ocsa, const void *ocsb) /* compare uoa by occurence */
{
    oc_t *oa = (oc_t*)ocsa; /* cast our void! */
    oc_t *ob = (oc_t*)ocsb; /* cast our void! */
    return ob->q  - oa->q; /* compare quanitities */
}

void prtoc(oc_t *ocs, int ocsz)
{
	int i;
	for(i=0;i<ocsz;++i)
		printf("%4c", ocs[i].l);
	printf("\n"); 
	for(i=0;i<ocsz;++i)
		printf("%4i", ocs[i].q);
	printf("\n"); 
	return;
}

oc2_t *convooc(char *aa, size_t aal, oc_t *ocs, int ocsz, int maxq, int *rnsz) // converted ordreed oc
{
	int i, j, nsz;
	for(i=0;i<ocsz;++i)
		if(ocs[i].q <= maxq) {
            nsz=i; // so cut off occurs here, meaning sz is the index (i.e. one lower)
            break;
        }

	oc2_t *oc2=calloc(nsz, sizeof(oc2_t));
	for(i=0;i<nsz;++i) {
		oc2[i].l = ocs[i].l;
		oc2[i].q = ocs[i].q;
		oc2[i].fap = 0;
		oc2[i].lap = 0;
    }

    unsigned char seenc;
	int *runoc=calloc(nsz, sizeof(int)); // running occurences
	for(i=0; i<aal; i++) {
		seenc=0;
		for(j=0;j<nsz;++j)
			if(aa[i]== oc2[j].l) {
                runoc[j]++;
                seenc=1;
                if(runoc[j] ==maxq+1)
			        oc2[j].fap=i;
                else if(runoc[j] == oc2[j].q)
			        oc2[j].lap=i;
                break;
            }
		if(seenc==0)
            continue;
    }
    // edge case, a letter is just one over the max allowed, that means lap=fap;
	for(j=0;j<nsz;++j)
	    if(oc2[j].q ==maxq+1)
            oc2[j].lap = oc2[j].fap;


    printf("Runoc:\n"); 
    for(i=0;i<nsz;++i) 
        printf("%i ", runoc[i]); 
    printf("\n"); 

    free(runoc);
    *rnsz=nsz;
	return oc2;
}

void expa_r(char *aa, size_t aal, oc2_t *oc2, int nsz, int maxq) /* following Gemini learnings */
{
	int i, j, kold=0;
    int minlen=(int)aal;

	int *exoc=calloc(nsz, sizeof(int)); // excess ocs
    int exsum, numz;

    int ri, k=280, oldri=(int)aal;
    while(k<360) {
	    for(i=0;i<nsz;++i)
            exoc[i] = oc2[i].q - maxq;
        numz=nsz;
        for(i=k; i<aal; i++) {
            for(j=0;j<nsz;++j) {
                if(aa[i]== oc2[j].l) {
                    if(exoc[j] >0) {
                        exoc[j]--;
                        if(exoc[j] ==0) {
                            ri = i;
                            numz--;
                            // if(numz==0) {
                            //     printf("@numz=0, letter=%c, seqpos=%i\n", aa[i], i); 
                            // }
                        }
                    }
                    break;
                }
            }
            if(numz==0) { 
                // kold=k;
                break;
            }
        }
        if(oldri!=ri) {
            printf("Explored window at kold=%i and oldri=%i len= %i\n", kold, oldri, oldri-kold+1); 
            if(minlen>(oldri-kold+1)) {
                minlen = oldri-kold+1;
                printf("new minlen=%i\n", minlen); 
            }
        }

        exsum=0;
        for(i=0;i<nsz;++i) 
            exsum += exoc[i];
        if(exsum!=0)
            printf("inval window at k=%i and ri=%i\n", k, ri); 

        kold=k;
        k++;
        oldri=ri;
    }

    printf("Exoc:\n"); 
    for(i=0;i<nsz;++i) 
        printf("%i ", exoc[i]); 
    printf("\n"); 

    free(exoc);
    return;
}

void prtoc2(oc2_t *ocs, int ocsz)
{
    int i;
    for(i=0;i<ocsz;++i)
        printf("%4c", ocs[i].l);
    printf("\n"); 
    for(i=0;i<ocsz;++i)
        printf("%4i", ocs[i].q);
    printf("\n"); 
    for(i=0;i<ocsz;++i)
        printf("%4i", ocs[i].fap);
    printf("\n"); 
    for(i=0;i<ocsz;++i)
        printf("%4i", ocs[i].lap);
    printf("\n"); 
    return;
}

void ocfirst(char *aa, oc_t **ocs, int *gbp, size_t aal, int *ocsz2)
{
    int i, j, ocsz=0;
    int gb=*gbp;
    oc_t *ocs2=*ocs;
    ocsz=*ocsz2;
    unsigned char seenc;
    for(i=0; aa[i]!='\0'; i++) {
        seenc=0;
        for(j=0;j<ocsz;++j) {
            if(aa[i]== ocs2[j].l) {
                seenc=1;
                ocs2[j].q++;
                break;
            }
        }
        if(seenc==0) {
            if(ocsz==gb-1) {
                gb+=GBUF;
                ocs2=realloc(ocs2, gb*sizeof(oc_t));
                for(j=gb-GBUF; j<gb; j++)
                    ocs2[j].q=0;
            }
            ocs2[ocsz].l=aa[i];
            ocs2[ocsz].q++;
            ocsz++;
        }
    }
    ocs2=realloc(ocs2, ocsz*sizeof(oc_t));
    // printf("equal quan should be %zu\n", aal/ocsz); 
    *ocs=ocs2;
    *ocsz2=ocsz;
    return;
}

int main(int argc, char *argv[])
{
    if(argc!=1) {
        usage(argv[0]);
        exit(EXIT_SUCCESS);
    }
    // introduce sequence as string literal
    char *aa="ASQLDRFRVFLGWDNGMMLVYQGNKTYEPWLNCDMASPTLSLVSKKAPKILKAADINTTLQPCLAFFIELLLKGIDNERIPNSGSGGREMGLLAPTYSSEATLVTRENNMMEGVHGFENMQDVEVIKLKLPEGYSDVCHFMFMLAGILYIVYDLQMHMSSERETGKFPNPLSDEEFDNPKVVVTNSFVLLEFTVTGAGARPSEQGQEPHNLGATKGSLAISSKTPEIHKDTNPASAQFEGKHTESDAKGVSNEDVFITKERDGREVEPTIKLSKKSVLNPMNVVYNPMLQISEGALRKHSMNDEITILNNTLINKERSVDLGAVSFVNDLLIDKLGCVSGKLAVQLNQSAPPEILHGKDPLTLFLGNTIALMLSKMQRIRVWEEYIFLNLHLALAWEPLLGNLKTHDSQKKWALCGFLIISRIRNLFESEGPVHGLRFSAMPCNTDTRQIKALERFPYAPEKPQWHGDELESPCRLVVASKLLASHDGVSIGKTIGSWPLPAQRYNAYVAWAANDSSILSARPGFAVKEDRLGHSLAQESGTIVVRNPQYGVRFINYNKDEHREFKREATFYPKTVVTHLGAIEGTLMFEIGDAAFTMLHLEEATDAEVRELYYMDMLDKKSSLGRACERIRRVLAPGDHKANGLESAIVSGQNGYEGRIRGLQTFQSNPLKKGRTHMAFCTTLHPFGGLKLVSSQLLKKELAVGTYGHQRTVLHSAEYSCPTSIPNLEGLMYNLISAQGEVNSDAKCHYAALAYICLQVRSVSMNQTEASDLRNFLETPILANDALASEQLLGSKKAKS";
    // char *aa="ASQLDRFRVFLGWDNGMMLVYQGNKTYEPWLNCDMA";

    size_t aal=strlen(aa);
    printf("aa length=%zu\n", aal); 

    int ocsz=0;

    // first pass
    int gb=GBUF;
    oc_t *ocs=calloc(gb, sizeof(oc_t));
    ocfirst(aa, &ocs, &gb, aal, &ocsz);
    if(aal%ocsz!=0) {
        printf("Error:length not dvisible by number of different aa's\n"); 
        exit(EXIT_FAILURE);
    }
    int maxq=aal/ocsz;

    qsort(ocs, ocsz, sizeof(oc_t), cmpocs2);

    prtoc(ocs, ocsz);

    // void convooc(char *aa, size_t aal, oc_t *ocs, int ocsz, int maxq) // converted ordreed oc
    int nsz;
    oc2_t *oc2 = convooc(aa, aal, ocs, ocsz, maxq, &nsz);

    prtoc2(oc2, nsz);

    // calc
    int i, fch=aal-1, lch=0; // first and last changes initial values
    for(i=0;i<nsz;++i) {
        if(oc2[i].fap<fch)
            fch=oc2[i].fap;
        if(oc2[i].lap>lch)
            lch=oc2[i].lap;
    }
    printf("lch=%i - fch=%i + 1 = %i\n", lch, fch, lch-fch +1); 

    expa_r(aa, aal, oc2, nsz, maxq);

    free(ocs);
    free(oc2);

    return 0;
}
