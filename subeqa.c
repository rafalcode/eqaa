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

typedef struct /*cx_t */
{
	char l; // letter
	int x; // excess over 40
} cx_t;

void usage(char *pname)
{
	printf("this prog extracts a subsequence of famous Nurtias string, given begin and end index integes (must be less than 800)\n");
	return;
}

int main(int argc, char *argv[])
{
    if(argc!=3) {
        usage(argv[0]);
        exit(EXIT_SUCCESS);
    }
    // introduce sequence as string literal
    char *aa="ASQLDRFRVFLGWDNGMMLVYQGNKTYEPWLNCDMASPTLSLVSKKAPKILKAADINTTLQPCLAFFIELLLKGIDNERIPNSGSGGREMGLLAPTYSSEATLVTRENNMMEGVHGFENMQDVEVIKLKLPEGYSDVCHFMFMLAGILYIVYDLQMHMSSERETGKFPNPLSDEEFDNPKVVVTNSFVLLEFTVTGAGARPSEQGQEPHNLGATKGSLAISSKTPEIHKDTNPASAQFEGKHTESDAKGVSNEDVFITKERDGREVEPTIKLSKKSVLNPMNVVYNPMLQISEGALRKHSMNDEITILNNTLINKERSVDLGAVSFVNDLLIDKLGCVSGKLAVQLNQSAPPEILHGKDPLTLFLGNTIALMLSKMQRIRVWEEYIFLNLHLALAWEPLLGNLKTHDSQKKWALCGFLIISRIRNLFESEGPVHGLRFSAMPCNTDTRQIKALERFPYAPEKPQWHGDELESPCRLVVASKLLASHDGVSIGKTIGSWPLPAQRYNAYVAWAANDSSILSARPGFAVKEDRLGHSLAQESGTIVVRNPQYGVRFINYNKDEHREFKREATFYPKTVVTHLGAIEGTLMFEIGDAAFTMLHLEEATDAEVRELYYMDMLDKKSSLGRACERIRRVLAPGDHKANGLESAIVSGQNGYEGRIRGLQTFQSNPLKKGRTHMAFCTTLHPFGGLKLVSSQLLKKELAVGTYGHQRTVLHSAEYSCPTSIPNLEGLMYNLISAQGEVNSDAKCHYAALAYICLQVRSVSMNQTEASDLRNFLETPILANDALASEQLLGSKKAKS";
    // char *aa="ASQLDRFRVFLGWDNGMMLVYQGNKTYEPWLNCDMA";
    int b=atoi(argv[1]);
    int e=atoi(argv[2]);
    int i, j;

    // L   S   A   E   G   K   V   N   T
    cx_t cx[9] = { {'L', 52}, {'S', 21}, {'A', 20}, {'E', 18}, {'G', 17}, {'K', 8}, {'V',6}, {'N',3}, {'T',1}};
    for(i=0;i<9;++i) 
        printf("{%c,%i}, ", cx[i].l, cx[i].x); 
    printf("\n"); 

    size_t aal=strlen(aa);
    printf("aa length=%zu\n", aal); 


    unsigned char seen;
    for(i=b;i<e+1;++i) {
        seen=0;
        for(j=0;j<9;++j)
            if(aa[i]==cx[j].l) {
                seen=1;
                if(cx[j].x>0) {
                    cx[j].x--;
                    if(cx[j].x==0)
                        printf("AA %c excess hits zero first at aaidx %i, 3gram=%c%c%c\n", cx[j].l, i, cx[j-1].l, cx[j].l, cx[j+1].l); 
                    break;
                }
            }
        if(seen==0)
            continue;
    }

    for(i=0;i<9;++i) 
        printf("{%c,%i}, ", cx[i].l, cx[i].x); 
    printf("\n"); 

    return 0;
}
