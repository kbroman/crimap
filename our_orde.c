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

static HEADER orders_base;
static HEADER *orders_allocp = NULL;

#define  NALLOC  6300
static HEADER *orders_morecore(nu)
        ALLOC nu;
{
    char *malloc();
    char *cp;
    HEADER *up;
    ALLOC rnu;
        SHORT our_orders_free();

    rnu = NALLOC * ((nu + NALLOC - 1) / NALLOC);
        printf("\n%ld bytes allocated in orders_morecore\n",
                                         rnu * sizeof(HEADER));
    cp = malloc(rnu * sizeof(HEADER));
    if (!cp) {
         printf("\nALLOC FAILED IN OUR_ORDERS_ALLOC\n");
         exit (1);
        }
    up = (HEADER *)cp;
    up->s.size = rnu;
    our_orders_free((char *)(up + 1));
    return(orders_allocp);
}


char *our_orders_alloc(nbytes)
        ALLOC nbytes;
{
    HEADER *orders_morecore();
    HEADER *p, *q;
    ALLOC nunits;

    nunits = 1 +(nbytes + sizeof(HEADER) - 1)/ sizeof(HEADER);
    if ((q = orders_allocp) == NULL) {
        orders_base.s.ptr = orders_allocp = q = &orders_base;
        orders_base.s.size = 0;
    }
    for (p = q->s.ptr;;q = p, p = p->s.ptr) {
        if (p->s.size >= nunits){
            if(p->s.size == nunits) q->s.ptr = p->s.ptr;
            else {
                p->s.size -= nunits;
                p += p->s.size;
                p->s.size = nunits;
            }
            orders_allocp = q;
            return ((char *)(p + 1));
        }
        if (p == orders_allocp)
            if ((p = orders_morecore(nunits)) == NULL) return (NULL);
    }
}

find_orders_frags()
{
    SHORT num_frags;
    HEADER *p;

    for (p = orders_allocp->s.ptr, num_frags = 0; p != orders_allocp;
            p = p->s.ptr, num_frags++);
        printf("\n number of allocation orders_fragments = %d \n", num_frags);
}


SHORT our_orders_free(ap)
char *ap;
{
    HEADER *p, *q;

    p = (HEADER *)ap - 1;
    for (q = orders_allocp; !(p > q && p < q->s.ptr); q = q->s.ptr)
        if (q>= q->s.ptr && (p > q || p < q->s.ptr)) break;

    if ( p + p->s.size == q->s.ptr){
        p->s.size += q->s.ptr->s.size;
        p->s.ptr = q->s.ptr->s.ptr;
    } else p->s.ptr = q->s.ptr;
    if (q + q->s.size == p){
        q->s.size += p->s.size;
        q->s.ptr = p->s.ptr;
    } else q->s.ptr = p;
    orders_allocp = q;
        return(0);
}
