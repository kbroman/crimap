#if vms
#include stdio
#include ctype
#else
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#endif

#include "defs.h"

#define MAX_LEN_FILE_NAME 800
#define MAX_NUM_LOCI 20000
#define  MAX_STR_LEN 3000
#define  MAX_HAP_SIZE 4000

struct loci_data *recs, *nrecs, *theta, *theta_t, *recs_temp,
   *nrecs_temp,  *num_mei, *num_mei_split, *pk_recs, *pk_nrecs, *theta_1_t;


double TOL;
SHORT PUK_NUM_ORDERS_TOL;
SHORT PK_NUM_ORDERS_TOL;

double PUK_LIKE_TOL;
double PK_LIKE_TOL;

SHORT SEX_EQ;
SHORT **FIXED_INTERVALS;
double *fam_likes;

/*following are globals for this file only -- except gen_file used in
chrompic.c */

ALLOC nb_our_alloc;
SHORT use_ord_file, write_ord_file, use_haps;
SHORT num_inserted, num_ordered;
char ord_file[MAX_LEN_FILE_NAME],par_file[MAX_LEN_FILE_NAME],
  dat_file[MAX_LEN_FILE_NAME],gen_file[MAX_LEN_FILE_NAME];
SHORT inserted_loci[MAX_NUM_LOCI];
SHORT *ordered_loci;

/*default values for global parameters*/
#define D_nb_our_alloc  300000000
#define D_SEX_EQ 1
#define D_use_ord_file 0
#define D_write_ord_file 1
#define D_use_haps 1
#define D_TOL 0.01
#define D_PUK_NUM_ORDERS_TOL 6
#define D_PK_NUM_ORDERS_TOL  8
#define D_PUK_LIKE_TOL 3.0
#define D_PK_LIKE_TOL 3.0

#define PREPARE      0
#define BUILD        1
#define INSTANT      2
#define QUICK        3
#define FIXED        4
#define FLIPS        5
#define ALL          6
#define TWOPOINT     7
#define CHROMPIC     8

main(argc,argv)
SHORT argc;
char *argv[];
{
     struct loci_orders *orders_orig, *orders, *orders_list;
     SHORT *sorted_list, *new_order;
     FILE *fp_data, *fp_param, *fp_orders;
     SHORT i, j, num_to_flip, num_loci, num_its;
     double *likelihoods;
     SHORT **orders_array;
     char *ctemp, *cp;
     struct chrom_data *pk_chrom_data, *puk_chrom_data;
     SHORT choice;
     SHORT dum_index;
     char curr_string[MAX_STR_LEN], var_name[MAX_STR_LEN];
     double *get_likelihoods();
     double mle(), chrompics_mle();
     SHORT num_tswitch_elim;
     struct loci_orders *get_all_orders(), *comp_ords_build_map(), *make_loci_orders();
     char *our_orders_alloc(), *our_alloc();
     SHORT our_free(), our_orders_free(), free_orders(), delete_hap();
     FILE *fopen();
     double llike;
     SHORT *insert_hap();
     char read_variables();

     /*Tests for compiler compatibility */
     if(NULL){
       printf("\n\nNULL is non-zero for your compiler;program will not run properly\n");
       exit(1);
     }
     if(sizeof(INT) != sizeof(ALLOC)){
       printf("\n\nYour compiler uses a different size for integers; see documentation\n");
       printf("for changes that will have to be made in the source code\n\n");
       exit(1);
     }

     /*Determine which option was selected*/
     if (!strcmp(*(argv+2),"prepare")) choice = PREPARE;
     else if (!strcmp(*(argv+2),"build")) choice = BUILD;
     else if (!strcmp(*(argv+2),"instant")) choice = INSTANT;
     else if (!strcmp(*(argv+2),"quick")) choice = QUICK;
     else if (!strcmp(*(argv+2),"fixed")) choice = FIXED;
     else if (!strcmp(*(argv+2),"all")) choice = ALL;
     else if (!strcmp(*(argv+2),"twopoint")) choice = TWOPOINT;
     else if (!strcmp(*(argv+2),"chrompic")) choice = CHROMPIC;
     else if (!strcmp(*(argv+2),"merge")) {
       merge_m();
       return;
     }
     else {
       strcpy(var_name, *(argv+2));
       var_name[5] = 0;
       if (!strcmp(var_name,"flips")) {
    /*for option flipsn, parse to get no. to flip*/
     choice = FLIPS;
     if (5 == strlen(*(argv+2))) num_to_flip = 2;
     else sscanf(*(argv+2),"flips%hd", &num_to_flip);
       }
       else {
     printf("\nERROR: %s is not an option\n",*(argv+2));
     exit(1);
       }
     }

     /*set default values for parameters*/
     cp = *(argv+1);
     set_defaults(cp);

     /*read in .par file*/
     if(!(fp_param = fopen(par_file,"r")) && choice != PREPARE){
       printf("\nERROR: Can't open .par file:  %s ; Run prepare first to create it\n", par_file);
       exit(1);
     }
     num_ordered = num_inserted = 0;
     ordered_loci = (SHORT *)our_orders_alloc((ALLOC)MAX_NUM_LOCI * sizeof(SHORT));
     set_null_hap();
     set_null_fixed();
     if (fp_param) {
       strcpy(var_name, "null");
       fscanf(fp_param, "%s", curr_string);
       for (; strcmp(curr_string, "END"); fscanf(fp_param, "%s", curr_string))
     if (read_variables(curr_string, var_name)) {
       if (choice == PREPARE) break;
       else {
         printf("\nERROR: Run prepare to create a new .par file\n");
         exit(1);
       }
     }
       fclose(fp_param);
     }

     if (choice == PREPARE) {
       prepare();
       return;
     }
     /*modify values in .par file: haplotypes used only if use_haps = 1;
       fixed_distances only apply to FIXED and CHROMPIC; use_ord_file value
       only applies to FLIPS and ALL; write_ord_file value only applies to BUILD */
     if (!use_haps) set_null_hap();
     if (choice != FIXED && choice != CHROMPIC) set_null_fixed();
     if (choice == BUILD || choice == QUICK || choice == INSTANT)
       use_ord_file = 1;
     else if (choice == FIXED || choice == CHROMPIC || choice == TWOPOINT)
       use_ord_file = 0;
     write_ord_file = write_ord_file && choice == BUILD;

     ctemp = (char *)our_alloc((ALLOC)nb_our_alloc);/*allocate initial memory block*/
     our_free(ctemp);

     /*read in data from .dat file*/
     if(!(fp_data = fopen(dat_file,"r"))){
       printf("\nERROR: Can't open .dat file: %s ; run prepare to convert .gen file to .dat file\n",dat_file);
       exit(1);
     }
     pk_chrom_data = (struct chrom_data *)our_alloc((ALLOC)sizeof(struct chrom_data));
     read_file(fp_data,pk_chrom_data);
     puk_chrom_data = (struct chrom_data *)our_alloc((ALLOC)sizeof(struct chrom_data));
     read_file(fp_data,puk_chrom_data);
     fclose(fp_data);

     /*print out values to be used in the analysis (except for loci)*/
     printf("\nOption chosen: %s\n",*(argv+2));
     print_variables();
     print_haps(puk_chrom_data->locus_names);
     print_fixed();

     /*delete secondary haplotyped loci from locus lists*/
     num_ordered = delete_hap(ordered_loci, num_ordered);
     num_inserted = delete_hap(inserted_loci, num_inserted);
     num_loci = num_ordered + num_inserted;

     if (use_ord_file) { /*read .ord file into a list of orders objects*/
       if(!(fp_orders = fopen(ord_file, "r"))){
     printf("\nERROR: can't open %s\n",ord_file);
     exit(1);
       }
       read_orders_file(fp_orders, &orders_list);
       fclose(fp_orders);
     }
     if (write_ord_file) { /*open .ord file for writing*/
        if(!(fp_orders = fopen(ord_file, "w"))){
          printf("\nERROR: can't write to the .ord file %s; run prepare to create it\n",ord_file);
          exit(1);
        }
    fclose(fp_orders);
      }

     /*initialize the interval and switch objects, allocate global parameters*/
     make_int_switch_vec();
     malloc_global(2 - SEX_EQ, 1);
     fam_likes = (double *)our_alloc((ALLOC)(puk_chrom_data->num_fams) *
                     sizeof(double));
     if (choice == BUILD || choice == ALL || choice == INSTANT || choice == QUICK) {
       if (num_ordered < 2){
     printf("\n\nERROR: fewer than 2 ordered loci after deleting 2dary haplotyped loci");
     exit(1);
       }
       /*set up initial orders object (from the ordered loci) */
       orders_array = (SHORT **)our_orders_alloc((ALLOC)sizeof(SHORT *));
       *orders_array = ordered_loci;
       orders = make_loci_orders(num_ordered, (LINDEX)1, orders_array);
       orders_orig = (struct loci_orders *)our_orders_alloc((ALLOC)sizeof(struct loci_orders));
       copy_orders(orders,orders_orig);

       /*print out a sorted list of the loci in the analysis, followed by the
     indices of the ordered and inserted loci*/
       sorted_list = (SHORT *)our_orders_alloc((ALLOC)num_loci*sizeof(SHORT));
       for(i = 0; i < num_ordered; i++) sorted_list[i] = orders_array[0][i];
       for(j = 0; j < num_inserted; j++, i++) sorted_list[i] = inserted_loci[j];
       short_sort(sorted_list, num_loci);
       print_names(puk_chrom_data->locus_names,sorted_list,num_loci,1,1);
       printf("\nordered loci:\n");
       print_array(ordered_loci, num_ordered);
       printf("\ninserted loci:\n");
       print_array(inserted_loci,num_inserted);
       printf("\n\n");
     }
     /*print names of loci to be used in the analysis*/
     else if (choice == FLIPS || choice == CHROMPIC)
       print_names(puk_chrom_data->locus_names,ordered_loci,num_ordered,1,1);
     else if (choice == TWOPOINT) {
       if (num_ordered)
     print_names(puk_chrom_data->locus_names,ordered_loci,num_ordered,5,0);
       if (num_ordered && num_inserted)
     printf("\n\nAGAINST:\n\n");
       if (num_inserted)
     print_names(puk_chrom_data->locus_names,inserted_loci,num_inserted,5,0);
     }

     if (choice == BUILD) {
       comp_ords_build_map(write_ord_file ? ord_file: "", orders, pk_chrom_data,
                puk_chrom_data,	inserted_loci,num_inserted,&orders_list);
       printf("\n\n\n");
       print_results1(orders_orig, puk_chrom_data, inserted_loci, num_inserted,
            orders_list,1,sorted_list,(char)SEX_EQ);
      }
      else if (choice == INSTANT || choice == QUICK) {
         free_orders(orders);
         print_results1(orders_orig, puk_chrom_data, inserted_loci, num_inserted,
                orders_list, (choice == INSTANT),sorted_list,(char)SEX_EQ);
      }
      else if (choice == ALL) {
         orders = get_all_orders(orders,inserted_loci,num_inserted );
     if (use_ord_file) test_and_compress(orders, orders_list);
         likelihoods = get_likelihoods(mle,orders,puk_chrom_data);
         print_best_orders(orders,likelihoods, PUK_LIKE_TOL);
       }
      else if (choice == FIXED || choice == CHROMPIC) {
    new_order = insert_hap(ordered_loci, num_ordered, &num_ordered, &dum_index);
    if (choice == FIXED) llike = mle(&num_its, new_order, num_ordered,
                     puk_chrom_data, &num_tswitch_elim);
    else llike = chrompics_mle(&num_its, new_order, num_ordered,
                   puk_chrom_data, &num_tswitch_elim);
    print_map(puk_chrom_data->locus_names,new_order, (choice == CHROMPIC));
         printf("\n\nlog10_like = %.3f\n\n",llike);
      }
      else if (choice == FLIPS) {
    printf("\nnumber of loci to flip = %d\n", num_to_flip);
    flipsn(ordered_loci,num_ordered,puk_chrom_data,num_to_flip,
           use_ord_file, orders_list);
       }
      else if (choice == TWOPOINT) {
    twopoint(ordered_loci, num_ordered, inserted_loci, num_inserted, puk_chrom_data);
       }
     return;
}

  /* set default parameter values here */
set_defaults(cp)
     char *cp;
{
  SHORT i;

  /*default file names*/
  if (isdigit(*cp)) {
    printf("\n\nchromosome %s\n", cp);
    sprintf(par_file,"chr%s.par",cp);
    sprintf(ord_file,"chr%s.ord",cp);
    sprintf(dat_file,"chr%s.dat",cp);
    sprintf(gen_file,"chr%s.gen",cp);
  }
  else {
    sprintf(par_file,"%s",cp);
    for (i = 0; cp[i] != '.'; i++);
    cp[i + 1] = '\0';
    sprintf(ord_file,"%sord",cp);
    sprintf(dat_file,"%sdat",cp);
    sprintf(gen_file,"%sgen",cp);
  }

  nb_our_alloc = D_nb_our_alloc;
  use_ord_file = D_use_ord_file;
  write_ord_file = D_write_ord_file;
  use_haps = D_use_haps;
  SEX_EQ = D_SEX_EQ;
  TOL = D_TOL;
  PUK_NUM_ORDERS_TOL = D_PUK_NUM_ORDERS_TOL;
  PK_NUM_ORDERS_TOL = D_PK_NUM_ORDERS_TOL;
  PUK_LIKE_TOL = D_PUK_LIKE_TOL;
  PK_LIKE_TOL = D_PK_LIKE_TOL;
}

print_variables()
{
  printf("\nCurrent values for parameters:\n");
  printf("\npar_file = %s",par_file);
  printf("\ndat_file = %s",dat_file);
  printf("\ngen_file = %s",gen_file);
  printf("\nord_file = %s",ord_file);
  printf("\nnb_our_alloc = %ld    [# bytes reserved for our_alloc]",nb_our_alloc);
  printf("\nSEX_EQ = %d   [0 = sex specific analysis, 1 = sex equal]", SEX_EQ);
  printf("\nTOL = %f", TOL);
  printf("\nPUK_NUM_ORDERS_TOL = %d", PUK_NUM_ORDERS_TOL);
  printf("\nPK_NUM_ORDERS_TOL = %d", PK_NUM_ORDERS_TOL);
  printf("\nPUK_LIKE_TOL = %.3f", PUK_LIKE_TOL);
  printf("\nPK_LIKE_TOL = %.3f", PK_LIKE_TOL);
  printf("\nuse_ord_file = %d", use_ord_file);
  printf("\nwrite_ord_file = %d", write_ord_file);
  printf("\nuse_haps = %d", use_haps);
  printf("\n\n");
}

write_variables(fp)
     FILE *fp;
{
  SHORT i;

  fprintf(fp,"dat_file  %s *",dat_file);
  fprintf(fp,"\ngen_file  %s *",gen_file);
  fprintf(fp,"\nord_file  %s *",ord_file);
  fprintf(fp,"\nnb_our_alloc  %ld *",nb_our_alloc);
  fprintf(fp,"\nSEX_EQ  %d *", SEX_EQ);
  fprintf(fp,"\nTOL  %f *", TOL);
  fprintf(fp,"\nPUK_NUM_ORDERS_TOL  %d *", PUK_NUM_ORDERS_TOL);
  fprintf(fp,"\nPK_NUM_ORDERS_TOL  %d *", PK_NUM_ORDERS_TOL);
  fprintf(fp,"\nPUK_LIKE_TOL  %.3f *", PUK_LIKE_TOL);
  fprintf(fp,"\nPK_LIKE_TOL  %.3f *", PK_LIKE_TOL);
  fprintf(fp,"\nuse_ord_file  %d *", use_ord_file);
  fprintf(fp,"\nwrite_ord_file  %d *", write_ord_file);
  fprintf(fp,"\nuse_haps  %d *", use_haps);

  if (num_ordered){
    fprintf(fp,"\nordered_loci ");
    for (i = 0; i < num_ordered; i++) fprintf(fp,"%d ",ordered_loci[i]);
    fprintf(fp," *");
  }
  if (num_inserted){
    fprintf(fp,"\ninserted_loci ");
    for (i = 0; i < num_inserted; i++) fprintf(fp,"%d ",inserted_loci[i]);
    fprintf(fp," *");
  }
}

/*following function is used to interpret strings in .par file or from keyboard
  (in prepare)*/
char read_variables(curr_string, var_name)
     char *curr_string, *var_name;
{
  static int num_in_hap, num_fixed;
  static double fixed_dist;
  static SHORT hap_list[MAX_HAP_SIZE], fixed_list[3];
  static char fix_zero;

  if (!strcmp(var_name, "null")) {
    strcpy(var_name, curr_string);
    num_in_hap = num_fixed = fix_zero = 0;
    if (!strcmp(var_name, "hap_sys0")) {
      strcpy(var_name,"hap_sys");
      fix_zero = 1;
    }
  }
  else if (!strcmp(curr_string, "*")) {
    if (!strcmp(var_name, "hap_sys")) read_haps(hap_list, num_in_hap, fix_zero);
    else if (!strcmp(var_name, "fixed_dist")) read_fixed(fixed_list, fixed_dist, num_fixed);
    strcpy(var_name, "null");
  }
  else if (!strcmp(curr_string, "="));
  else if (!strcmp(var_name, "nb_our_alloc")) sscanf(curr_string,"%ld",&nb_our_alloc);
  else if (!strcmp(var_name, "TOL")) sscanf(curr_string, "%lf", &TOL);
  else if (!strcmp(var_name, "PUK_NUM_ORDERS_TOL")) sscanf(curr_string, "%hd", &PUK_NUM_ORDERS_TOL);
  else if (!strcmp(var_name, "PK_NUM_ORDERS_TOL")) sscanf(curr_string, "%hd", &PK_NUM_ORDERS_TOL);
  else if (!strcmp(var_name, "PUK_LIKE_TOL")) sscanf(curr_string, "%lf", &PUK_LIKE_TOL);
  else if (!strcmp(var_name, "PK_LIKE_TOL")) sscanf(curr_string, "%lf", &PK_LIKE_TOL);
  else if (!strcmp(var_name, "dat_file")) sscanf(curr_string,"%s",dat_file);
  else if (!strcmp(var_name, "gen_file")) sscanf(curr_string,"%s",gen_file);
  else if (!strcmp(var_name, "ord_file")) sscanf(curr_string,"%s",ord_file);
  else if (!strcmp(var_name, "SEX_EQ")) sscanf(curr_string,"%hd",&SEX_EQ);
  else if (!strcmp(var_name, "use_ord_file")) sscanf(curr_string,"%hd",&use_ord_file);
  else if (!strcmp(var_name, "write_ord_file")) sscanf(curr_string,"%hd",&write_ord_file);
  else if (!strcmp(var_name, "use_haps")) sscanf(curr_string,"%hd",&use_haps);
  else if (!strcmp(var_name, "ordered_loci")) {
    sscanf(curr_string,"%hd", ordered_loci + num_ordered);
    num_ordered++;
  }
  else if (!strcmp(var_name, "inserted_loci")) {
    sscanf(curr_string,"%hd", inserted_loci + num_inserted);
    num_inserted++;
  }
  else if (!strcmp(var_name, "hap_sys")) {
    sscanf(curr_string,"%hd", hap_list + num_in_hap);
    num_in_hap++;
  }
  else if (!strcmp(var_name, "fixed_dist")) {
    if (!num_fixed) sscanf(curr_string,"%lf", &fixed_dist);
    else sscanf(curr_string,"%hd", fixed_list + num_fixed - 1);
    num_fixed++;
  }
  else {
    printf("\n\n ERROR: parameter name %s unknown\n\n", var_name);
    return (1);
  }
  return 0;
}

short_sort(v, n)
     SHORT *v;
     SHORT n;
/* Shell sort from p.58 of Kernighan-Ritchie */
{
  SHORT gap,i,j;
  SHORT temp;

  for(gap= n/2; gap > 0; gap /=2)
    for(i=gap; i < n; i++)
      for(j= i-gap; j >=0 && v[j]>v[j+gap]; j -=gap){
    temp = v[j];
    v[j] = v[j+gap];
    v[j+gap] = temp;
      }
}

print_names(names,indices,num_names,num_per_line, flag)
      char **names;
      SHORT *indices;
      SHORT num_names, num_per_line, flag;
/* prints list of names, with num_per_line on each line; if flag = 0,
      no indices are displayed, otherwise they are */
{
      SHORT i;

      printf("\n\n");
      for (i = 0; i < num_names; i++){
          if (flag) printf("%3d   %-20s",indices[i], names[indices[i]]);
          else printf("%-20s",names[indices[i]]);
          if (0 == (i+1)%num_per_line) printf("\n");
      }
      printf("\n\n");
}

print_map(names,indices,flag)
      char **names;
      SHORT *indices;
      SHORT flag;
/* displays map of loci; if flag = 0, original indices are displayed, otherwise
     origin 1 numbering starting with first locus in list
*/
{
      SHORT i, j, k;
      double cum[2];
      double kosambi();

      cum[0] = cum[1] = 0.0;
      if (theta->num_types == 2) printf("\n\nSex-specific ");
      else printf("\n\nSex_averaged ");
      printf("map (recomb. frac., Kosambi cM");
      if (theta->num_types == 2) printf(" -- female, male ");
      printf("):\n\n");
      for (i = 0; i < theta->n; i++){
          if (!flag)printf("%3d   %-20s",indices[i], names[indices[i]]);
          else
              printf("%3d   %-20s",i + 1, names[indices[i]]);

          for(k=0; k< theta->num_types; k++){
             printf("       %7.2f                    ",cum[k]);
             cum[k]+= kosambi(theta->data[k][i][0]);
          }
          printf("\n");
          for (j=0; j < theta->num_types; j++){
             printf("                    %5.3f%c%7.2f",theta->data[j][i][0],
            FIXED_INTERVALS[j][i]?'*':' ',kosambi(theta->data[j][i][0]) );
          }
          printf("\n");
      }
      if (!flag)printf("%3d   %-20s",indices[i], names[indices[i]]);
      else printf("%3d   %-20s",i + 1, names[indices[i]]);

      for (j=0; j < theta->num_types; j++)
         printf("       %7.2f                    ",cum[j]);

      printf("\n\n* denotes recomb. frac. held fixed in this analysis\n");
}

prepare()
{
  FILE *gp, *fp, *fopen();
  SHORT i, flag;
  SHORT ordered_flag, inserted_flag;
  struct chrom_data *chrom_data, *pk_chrom_data;
  struct data *data;
  SHORT *locus_nums, *pk_info, *puk_info;
  char ans[10];
  char loc_file[MAX_LEN_FILE_NAME];
  int choice;
  SHORT our_free();
  char *strcat();
/*  char *strcpy(); this needed with some compilers? */
  char *our_alloc(), *our_orders_alloc();
  SHORT *sort_by_info();
  SHORT *count_meioses();
  char curr_string[MAX_STR_LEN], var_name[MAX_STR_LEN];
  char read_variables();

  /*create new .dat file from .gen file, if none exists or if user desires it*/
  flag = 0;
  if (!(fp = fopen(dat_file, "r"))) printf("\nNo .dat file named %s\n",dat_file);
  else fclose(fp);
  if (!(gp = fopen(gen_file, "r"))) printf("\nNo .gen file named %s\n",gen_file);
  if(!fp && !gp) {
    printf("\nERROR: No data files were found ; You need at least a .gen file\n");
    exit(1);
  }
  if(fp && gp){
    printf("\n\nUse existing .dat file %s? (y/n) ",dat_file);
    scanf("%s",ans);
    flag = (ans[0] != 'y' && ans[0] != 'Y');
    if(flag){
      printf("\nAre you sure you want to delete file %s? (y/n) ",dat_file);
      scanf("%s",ans);
      flag = (ans[0] == 'y' || ans[0] == 'Y');
    }
    if (!flag) fclose(gp);
  }
  else if (!fp && gp) flag = 1;

  chrom_data = (struct chrom_data *)our_alloc((ALLOC)sizeof(struct chrom_data));
  pk_chrom_data = (struct chrom_data *)our_alloc((ALLOC)sizeof(struct chrom_data));
  data = (struct data *)our_alloc((ALLOC)sizeof(struct data));

  if (flag) {
    printf("\n\nCreating .dat file %s from .gen file %s\n",dat_file,gen_file);
    read_gen_file(gp,data,(char)1);
    alloc_c_dat(data, chrom_data, pk_chrom_data);
    fclose(gp);
    derive_gen(data,chrom_data,pk_chrom_data);
    printf("\n\nWriting file %s \n",dat_file);
    fp = fopen(dat_file, "w");
    write_file(fp,data,chrom_data,pk_chrom_data);
    fclose(fp);
    pk_info = count_meioses(pk_chrom_data);
    puk_info = count_meioses(chrom_data);
    printf("\n\nFinished writing %s\n",dat_file);
    strcpy(loc_file,gen_file);
    i = strlen(loc_file);
    loc_file[--i] = 'c';
    loc_file[--i] = 'o';
    loc_file[--i] = 'l';
    printf("\nWriting locus names to %s\n\n",loc_file);
    print_loci(data,loc_file,puk_info,pk_info);
    our_free(pk_info);
    our_free(puk_info);
  }

  printf("\n\n");
  /*display default (or previous .par file) values for parameters, and change them
    if user desires*/
  print_variables();
  printf("Do you wish to change any of these values? (y/n) ");
  scanf("%s",ans);
  if(ans[0] == 'y' || ans[0] == 'Y'){
    printf("\nTo change a value, enter the parameter name, the new value, and an asterisk,");
    printf("\nall separated by spaces; for example:\n\nTOL .001 * [Return]\n\n");
    printf("Type  done  when you are finished.\n");

    strcpy(var_name, "null");
    scanf("%s", curr_string);
    for (; strcmp(curr_string, "done"); scanf("%s", curr_string))
      read_variables(curr_string, var_name);

    print_variables();
  }

  if (flag) chrom_data = pk_chrom_data;
  else {
    if(!(fp = fopen(dat_file, "r"))){
      printf("\nCan't open %s\n",dat_file);
      exit(1);
    }
    read_file(fp, chrom_data);
    fclose(fp);
  }

  /*display locus names*/
  locus_nums= (SHORT *)our_alloc((ALLOC)chrom_data->num_loci * sizeof(SHORT));
  for(i = 0; i < chrom_data->num_loci; i++) locus_nums[i] = i;
  printf("\nThe loci and their indices are:");
  print_names(chrom_data->locus_names,locus_nums,chrom_data->num_loci,3,1);

  /*display haplotyped systems, and add new ones if user desires*/
  print_haps(chrom_data->locus_names);
  printf("\n\nDo you wish to enter any new haplotyped systems? (y/n) ");
  scanf("%s",ans);
  if(ans[0] == 'y' || ans[0] == 'Y'){
    printf("\nFor each new haplotyped system which you wish to enter, type either");
    printf("\nhap_sys0 (if distances between the loci are to be forced to equal 0)");
    printf("\nor hap_sys (if they aren't), followed by the indices of the loci to be haplotyped");
    printf("\n(separated by spaces), followed by * and a carriage return.  Example:\n\nhap_sys 2 0 5 * [Return]\n");
    printf("\nWhen you are done, type\n\ndone [Return]\n\nTo modify or delete a previously");
    printf("\nentered system, you will need to edit the .par file later with a text editor.\nReady: ");
    strcpy(var_name, "null");
    scanf("%s", curr_string);
    for (; strcmp(curr_string, "done"); scanf("%s", curr_string))
      read_variables(curr_string, var_name);
    print_haps(chrom_data->locus_names);
  }

  /*display fixed distances, and add new ones if user desires*/
  print_fixed();
  printf("\nDo you wish to hold any additional recombination fractions fixed (NB these will");
  printf("\nonly be used with the options FIXED and CHROMPIC, and only when the loci in question");
  printf("\n are adjacent)? (y/n) ");
  scanf("%s",ans);
  if(ans[0] == 'y' || ans[0] == 'Y'){
    printf("\nFor each new fixed rec. frac. which you wish to enter, type fixed_dist, followed");
    printf("\nby the recombination fraction, followed by the indices of the loci");
    printf("\nfollowed by the sex (0 for females, 1 for males) if the forced distance is to apply to one");
    printf("\nsex only, followed by * and a carriage return.  Example:\n\nfixed_dist 0.5 0 5 * [Return]\n");
    printf("\nWhen you are done, type\n\ndone [Return]\n\nTo modify or delete a previously");
    printf("\nentered system, you will need to edit the .par file later with a text editor.\nReady: ");
    strcpy(var_name, "null");
    scanf("%s", curr_string);
    for (; strcmp(curr_string, "done"); scanf("%s", curr_string))
      read_variables(curr_string, var_name);
    print_fixed();
  }

  /*determine which option will be run next*/
  printf("\nThe crimap options are:\n\n[1] build  [2] instant  [3] quick  [4] fixed \n");
  printf("\n[5] flips  [6] all  [7] twopoint  [8] chrompic\n");

  printf("\n\nEnter the number of the option you will be running next: ");
  scanf("%d",&choice);
  printf("\n\n");

  /*determine ordered_loci and inserted_loci*/
  ordered_flag = inserted_flag = 0;
  if (choice == BUILD || choice == INSTANT || choice == QUICK) {
    printf("\nDo you wish to build the map incorporating ALL loci in decreasing order of ");
    printf("\ntheir informativeness (the usual procedure)? (y/n) ");
    scanf("%s",ans);
    if(ans[0] == 'y' || ans[0] == 'Y'){
      our_free(locus_nums);
      locus_nums = sort_by_info(chrom_data);
      num_ordered = 2;
      num_inserted = chrom_data->num_loci -2;
      for(i = 0; i<num_ordered; i++) ordered_loci[i] = locus_nums[i];
      for(i = 0; i<num_inserted; i++) inserted_loci[i] = locus_nums[i+2];
    }
    else ordered_flag = inserted_flag = 1;
  }

  printf("\nThe loci and their indices are:");
  print_names(chrom_data->locus_names,locus_nums,chrom_data->num_loci,3,1);
  if (choice == ALL) ordered_flag = inserted_flag = 1;
  else if (choice == FLIPS || choice == FIXED || choice == CHROMPIC) {
    if (num_ordered > 2) {
      printf("\n The current ordered loci are ");
      print_names(chrom_data->locus_names, ordered_loci, num_ordered,3,1);
      printf("\n\nDo you wish to use these? (y/n)");
      scanf("%s",ans);
      if(ans[0] != 'y' && ans[0] != 'Y') ordered_flag = 1;
    }
    else ordered_flag = 1;
  }
  else if (choice == TWOPOINT){
    printf("\nDo you wish to compute LOD tables for ALL pairs of loci? (y/n)");
    scanf("%s",ans);
    if(ans[0] == 'y' || ans[0] == 'Y'){
      num_inserted = 0;
      ordered_loci = locus_nums;
      num_ordered = chrom_data->num_loci;
    }
    else {
      printf("\nYou may specify two separate groups of loci, ordered and inserted. If");
      printf("\nboth groups are non-empty, twopoint will only compute LOD tables for pairs");
      printf("\nconsisting of one locus from each group.  If one group is empty, LOD tables");
      printf("\nfor all pairs from the other group will be computed");
      ordered_flag = inserted_flag = 1;
    }
  }
  if (ordered_flag) {
    printf("\nType the indices of the ordered loci (separated by spaces), followed by a * :\n ");
    num_ordered = 0;
    strcpy(var_name, "ordered_loci");
    scanf("%s", curr_string);
    for (; strcmp(curr_string, "*"); scanf("%s", curr_string))
      read_variables(curr_string, var_name);
    printf("\n\nOrdered loci");
    print_names(chrom_data->locus_names, ordered_loci, num_ordered,3,1);
  }
  if (inserted_flag) {
    printf("\nType indices of loci to insert, followed by a *\n ");
    num_inserted = 0;
    strcpy(var_name, "inserted_loci");
    scanf("%s", curr_string);
    for (; strcmp(curr_string, "*"); scanf("%s", curr_string))
      read_variables(curr_string, var_name);
    printf("\nInserted loci");
    print_names(chrom_data->locus_names,inserted_loci,num_inserted,3,1);
  }

/*
  printf("\n Now specify the loci to be used in the analysis.  Each kind of analysis
 includes ordered_loci and/or inserted_loci; you specify them by typing (e.g.)
ordered_loci, then the locus indices, then
an asterisk.  The ordered_loci should be in their (assumed) chromosomal order. Example:
\n\ninserted_loci 2 5 7 1 0 * [Return]\n\n
For the options fixed, flips, and chrompic, all loci are ordered_loci. For twopoint,
if the lists of inserted and ordered loci are both non-empty, lods are only computed for pairs of
loci with one locus from each list; if only one list is non-empty, lods are computed for all pairs within
that list.  As a quick way to force the inserted_loci to contain all loci (in numerical order),
enter all_inserted_loci.  For the options build, instant, and quick, one normally sorts the loci by decreasing
informativeness and takes the first two as the ordered_loci (to start the map), and the remaining ones as
the inserted loci.  To force this, enter sort_loci.");

  printf("\n\n Now specify loci as above (type done when you are finished)");
  strcpy(var_name, "null");
  scanf("%s", curr_string);

  for (; strcmp(curr_string, "done"); scanf("%s", curr_string))
    read_variables(curr_string, var_name);
*/

  if (choice == BUILD) {
    /* now set up orders_file */
    if(!(fp = fopen(ord_file,"r"))) flag = 1;
    else {
      flag = 0;
      fclose(fp);
      printf("\n\nUse existing orders file? (y/n) ");
      scanf("%s",ans);
      if(ans[0] != 'y' && ans[0] != 'Y'){
    printf("\nAre you sure you want to delete the existing file? (y/n) ");
    scanf("%s",ans);
    if(ans[0] == 'y' || ans[0] == 'Y') flag = 1;
      }
    }
    if (flag && num_ordered > 2) {
      printf("\nDo you really want to assume this order for the %d ordered loci? (y/n) ",
         num_ordered);
      scanf("%s",ans);
      if(ans[0] != 'y' && ans[0] != 'Y') {
    printf("\nNo new orders file will being created now. Use prepare to respecify");
    printf("\nthe ordered loci");
    flag = 0;
      }
    }
    if (flag) {
      printf("\n\nCreating orders file %s\n\n",ord_file);
      fp = fopen(ord_file,"w");
      fprintf(fp,"1\n\n1  %d\n",num_ordered);
      for(i=0;i<num_ordered;i++) fprintf(fp,"%3d  ",ordered_loci[i]);
      fclose(fp);
      printf("\nDone.\n");
    }
  }

  /* Write new .par file */
  printf("\n\n\nOK to set up new parameter file? (y/n) ");
  scanf("%s",ans);
  if(ans[0] == 'y' || ans[0] == 'Y'){
    if (!(fp = fopen(par_file,"w"))) {
      printf("\n\nERROR: Can't create file %s\n\n", par_file);
      exit(1);
    }
    write_variables(fp);
    write_haps(fp);
    write_fixed(fp);
    fprintf(fp,"\nEND\n");
    fclose(fp);
    printf("\n%s has been created; use text editor for further modifications, if needed\n",par_file);
  }
  else printf("\n\naborted\n");

  return;
}

print_loci(data,out_file,puk_info,pk_info)
     struct data *data;
     char *out_file;
     SHORT *puk_info, *pk_info;
{
  SHORT i;
  FILE *fp, *fopen();

  if(!(fp = fopen(out_file,"w")))
    printf("\nUnable to create file with locus names\n");
  else {
    fprintf(fp,"Genotype file %s\n\n",gen_file);
    fprintf(fp,"Number of loci:  %d\n\n",data->num_loci);
    fprintf(fp,"                       #inf. mei.    #inf. mei.(phase known)\n");
    for(i = 0; i < data->num_loci; i++)
      fprintf(fp,"%3d  %-20s %4d          %4d\n", i, data->locus_names[i], puk_info[i],
          pk_info[i]);
    fclose(fp);
  }
  return;
}
