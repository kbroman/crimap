#define MAX_LEN_LOC_NAME 40
#define MAX_LEN_PKU_NAME 80
#define MAX_LEN_ORDS_NAME 80
#define NUM_PRE_COMP     20
#define MAX_DECREASE .000000001

/* MAX_DECREASE is the maximum decrease allowed in an EM step (necessary to allow
for occasional decreases due to rounding error or convergence to boundary) */

typedef long ALLOC;

/* if change type SHORT to anything larger than short need to read values
   in scanf & fscanf into a dummy variable of type short, then assign the
   dummy variable to the true variable.
*/

typedef short SHORT;
typedef long ID;
typedef long LINDEX;
typedef int INT;

struct loci_orders{
     SHORT num_loci;
     LINDEX num_orders;
     SHORT **orders;
     struct loci_orders *next_orders;	
};


struct loci_data{
     SHORT n, num_types;
     double ***data;
};


struct chrom_data{
     SHORT num_loci, num_fams;
     SHORT *num_chroms;  /* different number of chroms for each fam */
     struct phase **phase_choices;
     char ***chrom_array;
     char **fam_nums_array;
     char **locus_names;
};

struct phase{
     SHORT num_switches;
     SHORT *locus_nums;
     char **array;
};


struct int_array{
       SHORT n;
       struct intervals **list;
};

struct intervals {
       SHORT chrom_num, i, j, k, index, sort_index;
       unsigned rank : 4;
       unsigned r : 1;
  /* i = sex; j, k are starting, ending loci; r = reference recomb. count; index, sort_index give position in list */
       struct intervals *next_interval, *prev_interval;
};

struct interval_list{
       SHORT num_intervals;
       struct intervals *first_interval;
};

struct interval_ptrs{
       struct intervals *interval;
       struct interval_ptrs *next_ptr, *prev_ptr;
};

struct tswitchs{
       SHORT chrom_tswitch_index, num_intervals;
       struct interval_ptrs *first_ptr;
       struct tswitchs *next_tswitch, *prev_tswitch;
};

struct tswitch_list {
       SHORT num_tswitchs;
       struct tswitchs *first_tswitch;
};       

struct flank_tswitchs {
/*   struct tswitchs *tswitch;   ELIMINATE -- AND MAKE CORRESPONDING CHANGES IN PROGRAM? */
   LINDEX m_in_l_list, m_in_r_list;
   unsigned affects_l_interval : 1; 
   unsigned affects_r_interval : 1; 
   unsigned in_l_list : 1; 
   unsigned in_r_list : 1; 
   struct flank_tswitchs *next_tswitch;
};

struct flank_list {
   SHORT num_flanks;
   SHORT *num_tswitchs;
   LINDEX *m_left_off, *m_right_off;
   LINDEX *a_l_int, *a_r_int, *l_id_mask, *l_sh_mask;
   SHORT *n_in_l_list, *n_in_r_list;
   char *l_sh_dir;
   struct flank_tswitchs **first_tswitch;
};


struct individual{
     ID id, moth_id, fath_id, moth, fath;
     char sex;
     SHORT **a;
};

struct data{
     SHORT num_fams, num_loci;
     SHORT *num_mems;
     char **fam_id;
     char **locus_names;
     struct individual ***ind;
};


