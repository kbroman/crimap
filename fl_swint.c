#if vms
#include stdio
#else
#include <stdio.h>
#endif

#include "defs.h"
#include "var1.h"
#include "var2.h"

static SHORT int_index = 0;
static SHORT tsw_index = 0;
static SHORT int_ptr_index = 0;
/*
static SHORT fl_tsw_index = 0;
*/
struct int_block{
     struct intervals *int_ptr;
     struct int_block *next_block;
};

struct int_ptr_block{
     struct interval_ptrs *int_ptr;
     struct int_ptr_block *next_block;
};

struct tsw_block{
     struct tswitchs *tsw_ptr;
     struct tsw_block *next_block;
};
/*
struct fl_tsw_block{
     struct flank_tswitchs *fl_tsw_ptr;
     struct fl_tsw_block *next_block;
};
*/
static struct int_block *orig_int_block;
static struct int_block *curr_int_block;
static struct int_ptr_block *orig_int_ptr_block;
static struct int_ptr_block *curr_int_ptr_block;
static struct tsw_block *orig_tsw_block;
static struct tsw_block *curr_tsw_block;
/*
static struct fl_tsw_block *orig_fl_tsw_block;
static struct fl_tsw_block *curr_fl_tsw_block;
*/
reset_int_block()
{
     int_index = 0;
     curr_int_block = orig_int_block;
}

struct int_block *alloc_int_block()
{
         char *our_alloc();
         struct int_block *int_block;

         int_block = (struct int_block *)our_alloc((ALLOC)sizeof(struct int_block));
         int_block->int_ptr = (struct intervals *)our_alloc
             (512L * sizeof(struct intervals)); 
         int_block->next_block = 0; 
         return (int_block);
}

struct intervals *append_interval(first)
/* allocates a new empty interval, and appends to a interval_list */
    struct interval_list *first;
{
    struct intervals *interval;
    struct int_block *alloc_int_block();

    if (int_index > 511){
       if (!curr_int_block->next_block) 
            curr_int_block->next_block = alloc_int_block();
       curr_int_block = curr_int_block->next_block;
       int_index = 0;
    }
    interval = (curr_int_block->int_ptr) + (int_index++);
    interval->next_interval = first->first_interval;
    interval->prev_interval = 0;
    if (interval->next_interval) interval->next_interval->prev_interval = interval;
    first->first_interval = interval;
    first->num_intervals++;
    return(interval);
}

reset_tsw_block()
{
     tsw_index = 0;
     curr_tsw_block = orig_tsw_block;
}

struct tsw_block *alloc_tsw_block()
{
         char *our_alloc();
         struct tsw_block *tsw_block;

         tsw_block = (struct tsw_block *)our_alloc((ALLOC)sizeof(struct tsw_block));
         tsw_block->tsw_ptr = (struct tswitchs *)our_alloc
             (512L * sizeof(struct tswitchs)); 
         tsw_block->next_block = 0; 
         return (tsw_block);
}

struct tswitchs *append_tswitch(first)
/* allocates a new empty tswitch, and appends to a tswitch_list */
    struct tswitch_list *first;
{
    struct tswitchs *tswitch;
    struct tsw_block *alloc_tsw_block();

    if (tsw_index > 511){
       if (!curr_tsw_block->next_block) 
            curr_tsw_block->next_block = alloc_tsw_block();
       curr_tsw_block = curr_tsw_block->next_block;
       tsw_index = 0;
    }

    tswitch = (curr_tsw_block->tsw_ptr) + (tsw_index++);
    tswitch->num_intervals = 0;
    tswitch->first_ptr = 0;
    tswitch->next_tswitch = first->first_tswitch;
    tswitch->prev_tswitch = 0;
    if (tswitch->next_tswitch) tswitch->next_tswitch->prev_tswitch = tswitch;
    first->first_tswitch = tswitch;
    first->num_tswitchs++;
    return(tswitch);
}

reset_int_ptr_block()
{
     int_ptr_index = 0;
     curr_int_ptr_block = orig_int_ptr_block;
}

struct int_ptr_block *alloc_int_ptr_block()
{
         char *our_alloc();
         struct int_ptr_block *int_ptr_block;

         int_ptr_block = (struct int_ptr_block *)our_alloc
             ((ALLOC)sizeof(struct int_ptr_block));
         int_ptr_block->int_ptr = (struct interval_ptrs *)our_alloc
             (512L * sizeof(struct interval_ptrs)); 
         int_ptr_block->next_block = 0; 
         return (int_ptr_block);
}

struct interval_ptrs *append_interval_ptr(first)
/* allocates a new empty interval_ptr, and appends to a tswitch interval_ptr_list */
    struct tswitchs *first;
{
    struct interval_ptrs *interval_ptr;
    struct int_ptr_block *alloc_int_ptr_block();

    if (int_ptr_index > 511){
       if (!curr_int_ptr_block->next_block) 
            curr_int_ptr_block->next_block = alloc_int_ptr_block();
       curr_int_ptr_block = curr_int_ptr_block->next_block;
       int_ptr_index = 0;
    }

    interval_ptr = (curr_int_ptr_block->int_ptr) + (int_ptr_index++);
    interval_ptr->next_ptr = first->first_ptr;
    interval_ptr->prev_ptr = 0;
    if (interval_ptr->next_ptr) interval_ptr->next_ptr->prev_ptr = interval_ptr;
    first->first_ptr = interval_ptr;
    first->num_intervals++;
    return(interval_ptr);
}
/*
reset_fl_tsw_block()
{
     fl_tsw_index = 0;
     curr_fl_tsw_block = orig_fl_tsw_block;
}

struct fl_tsw_block *alloc_fl_tsw_block()
{
         char *our_alloc();
         struct fl_tsw_block *fl_tsw_block;

         fl_tsw_block = (struct fl_tsw_block *)our_alloc((ALLOC)sizeof(struct fl_tsw_block));
         fl_tsw_block->fl_tsw_ptr = (struct flank_tswitchs *)our_alloc
             (512L * sizeof(struct flank_tswitchs)); 
         fl_tsw_block->next_block = 0; 
         return (fl_tsw_block);
}

struct flank_tswitchs *append_flank_tswitch(first, i_fl)
    struct flank_list *first;
    SHORT i_fl; 
{
    struct flank_tswitchs *fl_sw;
    struct fl_tsw_block *alloc_fl_tsw_block();

    if (fl_tsw_index > 511){
       if (!curr_fl_tsw_block->next_block) 
            curr_fl_tsw_block->next_block = alloc_fl_tsw_block();
       curr_fl_tsw_block = curr_fl_tsw_block->next_block;
       fl_tsw_index = 0;
    }

    fl_sw = (curr_fl_tsw_block->fl_tsw_ptr) + (fl_tsw_index++);
    fl_sw->next_tswitch = first->first_tswitch[i_fl];
    fl_sw->in_l_list = 1;
    fl_sw->in_r_list = 1;
    fl_sw->m_in_l_list = 0;
    fl_sw->m_in_r_list = 0;
    fl_sw->affects_l_interval = 0;
    fl_sw->affects_r_interval = 0;
    first->first_tswitch[i_fl] = fl_sw;
    first->num_tswitchs[i_fl] += 1;
    return(fl_sw);
}
*/

make_int_switch_vec()
{
     char *our_alloc();
     struct int_block *alloc_int_block();
     struct tsw_block *alloc_tsw_block();
     struct int_ptr_block *alloc_int_ptr_block();
/*
     struct fl_tsw_block *alloc_fl_tsw_block();
*/
     orig_int_block = alloc_int_block();
     orig_tsw_block = alloc_tsw_block();
     orig_int_ptr_block = alloc_int_ptr_block();
/*
     orig_fl_tsw_block = alloc_fl_tsw_block();
*/
}


elim_tswitch(tswitch,first) 
    struct tswitchs *tswitch;
    struct tswitch_list *first;
{
    if (tswitch->prev_tswitch) tswitch->prev_tswitch->next_tswitch = tswitch->next_tswitch;
    else first->first_tswitch = tswitch->next_tswitch;
    if (tswitch->next_tswitch) tswitch->next_tswitch->prev_tswitch = tswitch->prev_tswitch;
    first->num_tswitchs--;
}

elim_interval(interval,first)
    struct intervals *interval;
    struct interval_list *first;
{
    if (interval->prev_interval) interval->prev_interval->next_interval = interval->next_interval;
    else first->first_interval = interval->next_interval;
    if (interval->next_interval) interval->next_interval->prev_interval = interval->prev_interval;
    first->num_intervals--;
}

elim_interval_ptr(interval_ptr,first)
    struct interval_ptrs *interval_ptr;
    struct tswitchs *first;
{
    if (interval_ptr->prev_ptr) interval_ptr->prev_ptr->next_ptr = interval_ptr->next_ptr;
    else first->first_ptr = interval_ptr->next_ptr;
    if (interval_ptr->next_ptr) interval_ptr->next_ptr->prev_ptr = interval_ptr->prev_ptr;
    first->num_intervals--;
}
