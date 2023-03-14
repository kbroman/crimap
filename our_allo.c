#if vms
#include stdio
#else
#include <stdio.h>
#endif

#include "defs.h"

typedef double ALIGN;

union header {
	struct {
		union header *ptr;
		ALLOC size;
	} s;
	ALIGN  x;
};

typedef union header HEADER;

static HEADER base;
static HEADER *allocp = NULL;	

char *our_alloc(nbytes)
        ALLOC nbytes;
{
	HEADER *morecore();
	HEADER *p;
	HEADER *q;
	ALLOC nunits;

	nunits = 1 +(nbytes + sizeof(HEADER) - 1)/ sizeof(HEADER);
	if ((q = allocp) == NULL) {
		base.s.ptr = allocp = q = &base;
		base.s.size = 0;
	}
	for (p = q->s.ptr;;q = p, p = p->s.ptr) {
		if (p->s.size >= nunits){
			if(p->s.size == nunits) q->s.ptr = p->s.ptr;
			else {
				p->s.size -= nunits;
				p += p->s.size;
				p->s.size = nunits;
			}
			allocp = q;
			return ((char *)(p + 1));
		}
		if (p == allocp)
			if ((p = morecore(nunits)) == NULL) return (NULL);
	}
}

find_frags()
{
	SHORT num_frags;
	HEADER *p;

	for (p = allocp->s.ptr, num_frags = 0; p != allocp;
            p = p->s.ptr, num_frags++);
        printf("\n number of allocation fragments = %d \n", num_frags);
}

#define  NALLOC  63000

static HEADER *morecore(nu)
        ALLOC nu;
{
	char *malloc();
	char *cp;
	HEADER *up;
	ALLOC rnu;
        SHORT our_free();

	rnu = NALLOC * ((nu + NALLOC - 1) / NALLOC);
        printf("\n%ld bytes allocated in morecore\n", rnu * sizeof(HEADER));
	cp = malloc(rnu * sizeof(HEADER));
	if (!cp) {
           printf("\n ERROR: ALLOCATION FAILED IN MORECORE\n");
           exit();  
        }
	up = (HEADER *)cp;
	up->s.size = rnu;
	our_free((char *)(up + 1));
	return(allocp);
}

SHORT our_free(ap)
char *ap;
{
	HEADER *p;
	HEADER *q;
	
	p = (HEADER *)ap - 1;
	for (q = allocp; !(p > q && p < q->s.ptr); q = q->s.ptr)
		if (q>= q->s.ptr && (p > q || p < q->s.ptr)) break;

	if ( p + p->s.size == q->s.ptr){
		p->s.size += q->s.ptr->s.size;
		p->s.ptr = q->s.ptr->s.ptr;
	} else p->s.ptr = q->s.ptr;
	if (q + q->s.size == p){
		q->s.size += p->s.size;
		q->s.ptr = p->s.ptr;
	} else q->s.ptr = p;
	allocp = q;
        return(0);
}
 
