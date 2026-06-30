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
    int fap; // first appearance L to R
    int lap; // last appearance L to R
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

int cmpocs(const void *ocsa, const void *ocsb) /* compare uoa by occurence */
{
    oc_t *oa = (oc_t*)ocsa; /* cast our void! */
    oc_t *ob = (oc_t*)ocsb; /* cast our void! */
    return oa->l  - ob->l; /* lexicographic character comparison: returns positive if b > a and negative if a > b: i.e. highest values first */
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
    qsort(ocs, ocsz, sizeof(oc_t), cmpocs2);
	for(i=0;i<ocsz;++i)
		printf("%4c", ocs[i].l);
	printf("\n"); 
	for(i=0;i<ocsz;++i)
		printf("%4i", ocs[i].q);
	printf("\n"); 
	return;
}

void convooc(oc_t *ocs, int ocsz, int maxq) // converted ordreed oc
{
	int i=0, nsz;
	for(i=0;i<ocsz;++i)
		if(ocs[i].q <= maxq) {
            nsz=i; // so cut off occurs here, meaning sz is the index (i.e. one lower)
            break;
        }

	oc2_t *oc2=calloc(nsz, sizeof(oc2_t));
	for(j=0;j<nsz;++j) {
		oc2[i].l = ocs[i].l;
		oc2[i].q = ocs[i].q;
		oc2[i].fap = 0;
		oc2[i].lap = 0;
    }

	for(i=0; aa[i]!='\0'; i++) {
		seenc=0;
		for(j=0;j<nsz;++j) {
			if(aa[i]== oc2[j].l) {
            }
		if(seenc==0)
            break;
            }


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
	return;
}

int prtnoverreps(oc_t *ocs, int ocsz, int maxq)
{
	int i, novreps=0 /* number of over-representatives */;
	for(i=0; i<ocsz; i++)
		if(ocs[i].q>maxq) {
			novreps+=ocs[i].q -maxq;
		}
	printf("Number by which overrepresented chars exceed %i = %i\n", maxq, novreps);
    printf("iThis is essentially the number of changes required for an equalised AA sequence.\n"); 
	return novreps;
}

int minlet(oc_t *ocs, int ocsz, int maxq)
{
	int i;
	int mn=maxq;
	int mnidx;
	for(i=0; i<ocsz; i++)
		if(ocs[i].q<mn) {
			mn=ocs[i].q;
			mnidx=i;
		}
	printf("minlet was %c, its quan was %i and it unsorted idx was %i\n", ocs[mnidx].l, ocs[mnidx].q, i); 
	return mnidx;
}

int *minlet2(oc_t *ocs, int ocsz, int novreps, int maxq)
{
	int i, j, k, times;
	int *isched=calloc(novreps, sizeof(int)); // index schedule
	k=0;
	for(i=0; i<ocsz; i++)
		if(ocs[i].q<maxq) {
			times=maxq - ocs[i].q;
			for(j=0;j<times;++j) 
				isched[k++]=i;
		}
	return isched;
}

a_t *catchl(char *aa, int aal) /* ctach Leucine occurrences */
{
	int i;
	a_t *arr=malloc(sizeof(a_t));
	arr->gb=GBUF;
	arr->sz=0;
	arr->els=malloc(arr->gb*sizeof(int));
	for(i=0;i<aal;++i)
		if(aa[i]=='L') {
			if(arr->sz == arr->gb-1) {
				arr->gb+=GBUF;
				arr->els=realloc(arr->els, arr->gb*sizeof(int));
			}
			arr->els[arr->sz]=i;
			arr->sz++;
		}
	arr->els=realloc(arr->els, arr->sz*sizeof(int));

	printf("Here are the positions of the L's\n"); 
	for(i=0;i<arr->sz;++i) 
		printf("%i ", arr->els[i]);
	printf("\n"); 

	return arr;
}

int *cleverroute(char *aa, int aal)
{
	int i, ib, ie, j, kb /* k beginning */, ke /* k end */;
	a_t *arrl=catchl(aa, aal);
	int *ia=malloc(aal*sizeof(int)); // index array
	ib=0;
	ie=aal-1;
	j=0;
	kb=0;
	ke=arrl->sz-1;
	// int bwin=0, ewin=0;
	for(;;) {
		if(j>aal)
			break;
		if(arrl->els[kb]-ib > ie-arrl->els[ke]) {
			for(i=ib;i<=arrl->els[kb];++i) 
				ia[j++]=i;
			ib=1+arrl->els[kb];
			kb++;
		// } else if(arrl->els[kb]-ib <= ie-arrl->els[ke]) {
		} else {
			for(i=ie;i>=arrl->els[ke];--i) 
				ia[j++]=i;
			ie=arrl->els[ke]-1;
			ke--;
		}
	}
	printf("\n"); 
	free(arrl->els);
	free(arrl);
	return ia;
}


int *cleverroute2(char *aa, int aal)
{
	/* OK, this may look like an easy thing on the outside, but it's not really */
	/* Beware indez 0 will get left out,, it indexing-from-beginning never gets a go */
	int i, ib, ie, j, kb /* k beginning */, ke /* k end */;
	a_t *arrl=catchl(aa, aal);
	int *ia=malloc(aal*sizeof(int)); // index array
	ib=0;
	ie=aal-1;
	j=0;
	kb=0;
	ke=arrl->sz-1;
	int fromend=1;
	int endp;
	for(;;) {
		// if(j>aal-2) {
		// printf("%i:%i \n", kb, ke); 
		if(ke == kb-2) { /* pretty incredibly: yes */
// 			printf("%zu : %i ", aal, j); 
			break;
		}
		if(fromend) {
			endp=(ke==kb-1)?arrl->els[ke]+1:arrl->els[ke];
			// printf("From %i to %i:", ie, endp);
			for(i=ie;i>=endp;--i) {
				ia[j++]=i;
				// printf("%i ", ia[j-1]);
			}
			// printf("\n"); 
			ie=arrl->els[ke]-1; // the new i-end for next time
			ke--;
			fromend++;
			if(fromend==2)
				fromend=0;
		} else {
			endp=(ke==kb-1)?arrl->els[kb]-1:arrl->els[kb];
			// printf("From %i to %i:", ib, endp);
			for(i=ib;i<=endp;++i) {
				ia[j++]=i;
				// printf("%i ", ia[j-1]);
			}
			// printf("\n");
			ib=1+arrl->els[kb];
			kb++;
			fromend=1;
		}
	}
	free(arrl->els);
	free(arrl);
	return ia;
}

void printaaoc(char *aa, oc_t *ocs2, int ocsz)
{
	int i;
	printf("%i aa's were seen\n", ocsz); 
	for(i=0;i<ocsz;++i)
		printf("%4c", ocs2[i].l);
	printf("\n"); 
	for(i=0;i<ocsz;++i)
		printf("%4i", ocs2[i].q);
	printf("\n"); 

	return;
}

void printcmp2str(char *aa, char *aa2)
{
	int i, j;
	size_t aal1=strlen(aa);
	size_t aal2=strlen(aa2);
	if(aal1 != aal2) {
		printf("Sorry, this comparison only allow equal sized strings.\n"); 
		exit(EXIT_FAILURE);
	}
	int numchnks=aal1/ROWSZ;
	int partchnk=aal1%ROWSZ;
	for(i=0;i<aal1;++i) 
		if(aa[i]==aa2[i])
			aa2[i]='.';
	if(partchnk)
		numchnks++;
	int *chnka=calloc(numchnks, sizeof(int));
	for(i=0;i<numchnks-1;++i) 
		chnka[i]=ROWSZ;
	chnka[numchnks-1]=(partchnk)?partchnk:ROWSZ;
	for(i=0;i<numchnks;++i) {
		for(j=0;j<chnka[i];++j)
			putchar(aa[i*ROWSZ+j]);
		printf("\n"); 
		for(j=0;j<chnka[i];++j)
			putchar(aa2[i*ROWSZ+j]);
		printf("\n\n"); 
	}

	free(chnka);
	return;
}

void ocgo(char *aa)
{
	int i, j, ocsz=0;
	int gb=GBUF;
	oc_t *ocs=calloc(gb, sizeof(oc_t));
	unsigned char seenc;
	for(i=0; aa[i]!='\0'; i++) {
		seenc=0;
		for(j=0;j<ocsz;++j) {
			if(aa[i]== ocs[j].l) {
				seenc=1;
				ocs[j].q++;
				break;
			}
		}
		if(seenc==0) {
			if(ocsz==gb-1) {
				gb+=GBUF;
				ocs=realloc(ocs, gb*sizeof(oc_t));
				for(j=gb-GBUF; j<gb; j++)
					ocs[j].q=0;
			}
			ocs[ocsz].l=aa[i];
			ocs[ocsz].q++;
			ocsz++;
		}
	}
	ocs=realloc(ocs, ocsz*sizeof(oc_t));
	free(ocs);
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

void ocfirst0(char *aa, oc_t **ocs, int *gbp, size_t aal, int *ocsz2, int maxq)
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
				if(ocs2[j].q == maxq+1)
				    ocs2[j].fap == i;
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

void ocfirst2(char *aa, oc2_t **ocs, int *gbp, size_t aal, int *ocsz2, int maxq)
{
	int i, j, ocsz=0;
	int gb=*gbp;
	oc2_t *ocs2=*ocs;
	ocsz=*ocsz2;
	unsigned char seenc;
	for(i=0; aa[i]!='\0'; i++) {
		seenc=0;
		for(j=0;j<ocsz;++j) {
			if(aa[i]== ocs2[j].l) {
				seenc=1;
				ocs2[j].q++;
                if(ocs2[j].q == maxq+1)
				    ocs2[j].fap = i;
				break;
			}
		}
		if(seenc==0) {
			if(ocsz==gb-1) {
				gb+=GBUF;
				ocs2=realloc(ocs2, gb*sizeof(oc2_t));
				for(j=gb-GBUF; j<gb; j++) {
					ocs2[j].q=0;
					ocs2[j].fap=0;
                }
			}
			ocs2[ocsz].l=aa[i];
			ocs2[ocsz].q++;
			ocsz++;
		}
	}
	ocs2=realloc(ocs2, ocsz*sizeof(oc2_t));
	// printf("equal quan should be %zu\n", aal/ocsz); 
	*ocs=ocs2;
	*ocsz2=ocsz;
	return;
}

char *occheckdynbe(char *aa, int aal, oc_t **ocs, int *ocsz2, int *gbp, int mxq, int *nchanges) // check from both endds .. dynamically.
{
	/* in this version we jump from one end of the string to the other in an attempt to corral in the substring from both sides */
	int i, j, k;
	// a_t *arrl=catchl(aa, aal);
	int *ia=cleverroute2(aa, aal);
#ifdef DBG
	for(i=0;i<aal;++i) 
		printf("%i ", ia[i]); 
	printf("\n"); 
#endif
	char *aa2=malloc((aal+1)*sizeof(char));

	strcpy(aa2, aa);
	int gb=*gbp;
	oc_t *ocs2=*ocs;
	int ocsz=*ocsz2;
	prtoc(ocs2, ocsz);
	unsigned char seenc; /* seen character .. has it already been seen? or is ithis a first time? */
	int minp=aal, maxp=0;
	*nchanges=0;
	for(i=0; i<aal; i++) {
		seenc=0;
		for(j=0;j<ocsz;++j) {
			if(aa2[ia[i]]== ocs2[j].l) {
				seenc=1;
				if(ocs2[j].q > mxq) {
					printf("Max quan found at %i with %c giving %i  ...", ia[i], ocs2[j].l, ocs2[j].q); 
					k=minlet(ocs2, ocsz, mxq);
					printf("minchar %c:%i\n", ocs2[k].l, ocs2[k].q);
					aa2[ia[i]]=ocs2[k].l;
					ocs2[k].q++;
					/* in context of bouncing index, not obvious how to mark extreme lower and extreme upper index .. but here's how: */
					if(ia[i] < minp)
						minp=ia[i];
					if(ia[i] > maxp)
						maxp=ia[i];
					ocgo(aa2);
					(*nchanges)++;
					printf("=> nchanges=%i\n", *nchanges); 
				} else
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
			ocs2[ocsz].l=aa2[ia[i]];
			ocs2[ocsz].q++;
			ocsz++;
		}
	}
	ocs2=realloc(ocs2, ocsz*sizeof(oc_t));
	printf("Summary: Num changes to input string=%i, first change at %i and last change at %i\n", *nchanges, minp, maxp); 
	int szsub=maxp-minp+1;
	printf("Length of smallest substring for changes = %i, %2.2f times over the min num changes\n", szsub, (float)szsub/(*nchanges));
	*ocs=ocs2;
	aa2[aal]='\0';
	free(ia);
	return aa2;
}

char *occheckdynbe2(char *aa, int aal, oc_t **ocs, size_t *ocsz2, int *gbp, int mxq, int novreps, int *nchanges, int *isched) // check from both endds .. dynamically.
{
	/* in this version we jump from one end of the string to the other in an attempt to corral in the substring from both sides */
	int i, ii=0, j, k;
	// a_t *arrl=catchl(aa, aal);
	int *ia=cleverroute2(aa, aal);
#ifdef DBG
	for(i=0;i<aal;++i) 
		printf("%i ", ia[i]); 
	printf("\n"); 
#endif
	char *aa2=malloc((aal+1)*sizeof(char));
	strcpy(aa2, aa); // aa fully copied over to aa2 so 

	int gb=*gbp;
	oc_t *ocs2=*ocs;
	int ocsz=*ocsz2;
	unsigned char seenc; /* seen character .. has it already been seen? or is ithis a first time? */
	int minp=aal, maxp=0;
	*nchanges=0;
	for(i=0; i<aal; i++) {
		seenc=0;
		for(j=0;j<ocsz;++j) {
			if(aa2[ia[i]]== ocs2[j].l) {
				seenc=1;
				if(ocs2[j].q > mxq) {
					// printf("%c over maxq found at idx %i with %c giving %i  ...", ia[i], ocs2[j].l, ocs2[j].q); 
					k=isched[ii++];
					printf("minchar %c:%i\n", ocs2[k].l, ocs2[k].q);
					aa2[ia[i]]=ocs2[k].l;
					ocs2[k].q++;
					/* in context of bouncing index, not obvious how to mark extreme lower and extreme upper index .. but here's how: */
					if(ia[i] < minp)
						minp=ia[i];
					if(ia[i] > maxp)
						maxp=ia[i];
					ocgo(aa2);
					prtoc(ocs2, ocsz);
					(*nchanges)++;
					printf("=> nchanges=%i\n", *nchanges); 
				} else
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
			ocs2[ocsz].l=aa2[ia[i]];
			ocs2[ocsz].q++;
			ocsz++;
		}
	}
	ocs2=realloc(ocs2, ocsz*sizeof(oc_t));
	printf("Summary: Num changes to input string=%i, first change at %i and last change at %i\n", *nchanges, minp, maxp); 
	int szsub=maxp-minp+1;
	printf("Length of smallest substring for changes = %i, %2.2f times over the min num changes\n", szsub, (float)szsub/(*nchanges));
	*ocs=ocs2;
	aa2[aal]='\0';
	free(ia);
	return aa2;
}

char *occheckbe(char *aa, int aal, oc_t **ocs, int *gbp, int mxq, int *firstch, int *lastch, int *nchanges) // check from both ends
{
	/* in this version we jump from one end of the string to the other in an attempt to corral in the substring from both sides */
	int i, j, k, ocsz=0;
	char *aa2=malloc((aal+1)*sizeof(char));
	int *ia=malloc(aal*sizeof(int)); // index array
	j=0;
	int halfway=aal/2;
	for(i=aal-1;i>=halfway;i--) {
		ia[j]=i;
		j+=2;
	}
	j=1;
	for(i=0;i<halfway;i++) {
		ia[j]=i;
		j+=2;
	}
	FLAG firstchangeeven, nolastch;
	strcpy(aa2, aa);
	int gb=*gbp;
	oc_t *ocs2=*ocs;
	unsigned char seenc;
	*nchanges=0;
	for(i=0; i<aal; i++) {
		seenc=0;
		for(j=0;j<ocsz;++j) {
			if(aa2[ia[i]]== ocs2[j].l) {
				seenc=1;
				if(ocs2[j].q == mxq) {
#ifdef DBG
					printf("Max quan found at %i with %c\n", ia[i], ocs2[j].l); 
#endif
					k=minlet(ocs2, ocsz, mxq);
#ifdef DBG
					printf("minchar %c:%i\n", ocs2[k].l, ocs2[k].q);
#endif
					aa2[ia[i]]=ocs2[k].l;
					ocs2[k].q++;
					if(*nchanges==0) {
						*firstch=ia[i];
						firstchangeeven=(i%2)?0:1; // note how I work here with i, not ia[i] because i have in mind how to deal with evens and unevens
						nolastch=1;
					} else if( (!firstchangeeven) & nolastch & !(i%2)) {
						*lastch=ia[i]; // last change in this case need much cleverer calculation it the next opposit change!
						nolastch=0;
					} else if( firstchangeeven & nolastch & (i%2)) {
						*lastch=ia[i]; // last change in this case need much cleverer calculation it the next opposit change!
						nolastch=0;
					}
					(*nchanges)++;
				} else
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
			ocs2[ocsz].l=aa2[ia[i]];
			ocs2[ocsz].q++;
			ocsz++;
		}
	}
	ocs2=realloc(ocs2, ocsz*sizeof(oc_t));
	for(i=0;i<ocsz;++i)
		printf("%4c", ocs2[i].l);
	printf("\n"); 
	for(i=0;i<ocsz;++i)
		printf("%4i", ocs2[i].q);
	printf("\n"); 
	*ocs=ocs2;
	aa2[aal]='\0';
	free(ia);
	return aa2;
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

	int i, ocsz=0;

	// first pass
	int gb=GBUF;
	oc_t *ocs=calloc(gb, sizeof(oc_t));
	ocfirst(aa, &ocs, &gb, aal, &ocsz);
	if(aal%ocsz!=0) {
		printf("Error:length not dvisible by number of different aa's\n"); 
		exit(EXIT_FAILURE);
	}
	int maxq=aal/ocsz;
	int novreps=prtnoverreps(ocs, ocsz, maxq);
	prtoc(ocs, ocsz);

	int gb2=GBUF, oc2sz=0;
	oc2_t *oc2=calloc(gb2, sizeof(oc2_t));
	ocfirst2(aa, &oc2, &gb2, aal, &oc2sz, maxq);
	prtoc2(oc2, oc2sz);

	// second pass
	gb=GBUF;
	int *isched=minlet2(ocs, ocsz, novreps, maxq);
	for(i=0;i<novreps;++i) 
		printf("%i ", isched[i]); 
	printf("\n"); 

    return 0;

	int nchanges=0;
	size_t ocsz2=0;
	oc_t *ocs2=calloc(gb, sizeof(oc_t));
	char *aa2=occheckdynbe2(aa, aal, &ocs2, &ocsz2, &gb, maxq, novreps, &nchanges, isched);

	// printaaoc(aa2, ocs2, ocsz);
	// printcmp2str(aa, aa2);


	//printf("Summary: Num changes to input string=%i, first change at %i and last change at %i\n", nchanges, firstch, lastch); 
	//int szsub=firstch-lastch+1; // measured this and saw that if you subtract indices or positions you need to add 1.
	//printf("Length of smallest substring for changes = %i, %2.2f times over the min num changes\n", szsub, (float)szsub/nchanges);

	free(aa2);
	free(isched);
	free(ocs2);

	return 0;
}
