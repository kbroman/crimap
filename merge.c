#if vms
#include stdio
#include ctype
#else
#include <stdio.h>
#include <ctype.h>
#endif

#include "defs.h"

merge_m()
{
  struct data *data, *data1, *data2;
  struct data *merge();
  char file_name[20];
  FILE *fopen(), *fp;
  char *our_alloc();
  SHORT our_free();

  printf("\n\n\nfirst input file = ");
  scanf("%s",file_name);
  printf("\n\ninput file = %s\n",file_name);
  data1 = (struct data *)our_alloc((ALLOC)sizeof(struct data));
  if(!(fp = fopen(file_name, "r"))){
    printf("\nUnable to open %s\n",file_name);
    exit(1);
  }
  read_gen_file(fp, data1, 1);
  fclose(fp);

  printf("\n\n\nsecond input file = ");
  scanf("%s",file_name);
  printf("\n\ninput file = %s\n",file_name);
  data2 = (struct data *)our_alloc((ALLOC)sizeof(struct data));
  if(!(fp = fopen(file_name, "r"))){
    printf("\nUnable to open %s\n",file_name);
    exit(1);
  }
  read_gen_file(fp, data2, 1);
  fclose(fp);

  data = merge(data1,data2);

  printf("\n\n\noutput file = ");
  scanf("%s",file_name);
  printf("\n\noutput file = %s\n",file_name);

  if(!(fp = fopen(file_name, "w"))){
    printf("\nUnable to write to output file %s\n",file_name);
    exit(1);
  }

  gen_write(fp, data);
  fclose(fp);

  return;
}


gen_write(fp, data)
     FILE *fp;
     struct data *data;
{
  SHORT i,j,k;

  printf("\n\nwriting file ...\n");
  fprintf(fp,"%d  ", data->num_fams);
  fprintf(fp,"%d\n", data->num_loci);

  for(k = 0; k < data->num_loci; k++){
    fprintf(fp,"%s\n", data->locus_names[k]);
  }

  for(i=0; i<data->num_fams; i++){
    fprintf(fp,"%s\n", data->fam_id[i]);
    fprintf(fp,"%d\n", data->num_mems[i]);
    for(j = 0; j<data->num_mems[i]; j++){
      fprintf(fp,"%ld ", data->ind[i][j]->id);
      fprintf(fp,"%ld ", data->ind[i][j]->moth_id);
      fprintf(fp,"%ld ", data->ind[i][j]->fath_id);
      fprintf(fp,"%d\n", data->ind[i][j]->sex);

      for(k = 0; k<data->num_loci; k++){
    fprintf(fp,"%d ", data->ind[i][j]->a[k][0]);
    fprintf(fp,"%d ", data->ind[i][j]->a[k][1]);
      }
      fprintf(fp,"\n");
    }  /* close j */

  } /* close i */
  return;
}

print_gen( data)
     struct data *data;
{
  SHORT i,j,k;

  printf("\n\nPrint_gen\n\n");
  printf("%d\n", data->num_fams);
  printf("%d\n", data->num_loci);

  for(k = 0; k < data->num_loci; k++){
    printf("%s\n", data->locus_names[k]);
  }

  for(i=0; i<data->num_fams; i++){
    printf("%s\n", data->fam_id[i]);
    printf("%d\n", data->num_mems[i]);
    for(j = 0; j<data->num_mems[i]; j++){
      printf("%ld ", data->ind[i][j]->id);
      printf("%ld ", data->ind[i][j]->moth_id);
      printf("%ld ", data->ind[i][j]->fath_id);
      printf("%d\n", data->ind[i][j]->sex);

      for(k = 0; k<data->num_loci; k++){
    printf("%d ", data->ind[i][j]->a[k][0]);
    printf("%d ", data->ind[i][j]->a[k][1]);
      }
      printf("\n");
    }  /* close j */

  } /* close i */
  return;
}


struct data *merge(data1,data2)
     struct data *data1, *data2;
{
  SHORT nf, nl, i, j, k, m, flag, num_fams, num_loci, i2, ifam2, jind2, j1, i1;
  SHORT *loc2_pos;
  struct data *data;
  struct individual *dat_ind;
  char *our_alloc();

  printf("\n\nmerge\n");
  data = (struct data *)our_alloc((ALLOC)sizeof(struct data));
  nf = (num_fams = data1->num_fams) + data2->num_fams;
  nl = (num_loci = data1->num_loci) + data2->num_loci;
  data->locus_names = (char **)our_alloc((ALLOC)nl * sizeof(char *));
  loc2_pos = (SHORT *)our_alloc((ALLOC)data2->num_loci * sizeof(SHORT));

  data->num_mems = (SHORT *)our_alloc((ALLOC)nf * sizeof(SHORT));
  data->fam_id = (char **)our_alloc((ALLOC)nf * sizeof(char *));
  data->ind = (struct individual ***)our_alloc((ALLOC)nf * sizeof(struct individual **));

  for (i = 0; i < data1->num_loci; i++)
    data->locus_names[i] = data1->locus_names[i];
  for (i = 0; i < data2->num_loci; i++){
    flag = 0;
    for (j = 0; j < data1->num_loci; j++){
      if (!strcmp(data2->locus_names[i], data1->locus_names[j])) {
    flag = 1;
    loc2_pos[i] = j;
    break;
      }
    }
    if (!flag) {
      loc2_pos[i] = num_loci;
      data->locus_names[num_loci++] = data2->locus_names[i];
    }
  }
  data->num_loci = num_loci;

  ifam2 = -1;
  for (i = 0; i < nf; i++){
    if (i < data1->num_fams) {
      data->fam_id[i] = data1->fam_id[i];
      data->num_mems[i] = data1->num_mems[i];

      for (i2 = 0; i2 < data2->num_fams; i2++)
    if (!strcmp(data2->fam_id[i2],data1->fam_id[i])) {
      nf--;
      data->num_mems[i] += data2->num_mems[i2];
      for (k = 0; k < data2->num_mems[i2]; k++)
        for (m = 0; m < data1->num_mems[i]; m++)
          if (data2->ind[i2][k]->id == data1->ind[i][m]->id) {
        data->num_mems[i] -= 1;
        if (data2->ind[i2][k]->moth_id != data1->ind[i][m]->moth_id
            ||data2->ind[i2][k]->fath_id != data1->ind[i][m]->fath_id
            ||data2->ind[i2][k]->sex != data1->ind[i][m]->sex)
          printf("\n\n MISMATCH IN PEDIGREE %ld, individual %ld ",
             data->fam_id[i],data2->ind[i2][k]->id);
        break;

          }

      break;
    }
    }
    else {
      ifam2++;
      do {
    flag = 0;
    for (i1 = 0; i1 < data1->num_fams; i1++)
      if (!strcmp(data2->fam_id[ifam2], data1->fam_id[i1])){
        ifam2++;
        flag = 1;
        break;
      }

      } while (flag);
      data->fam_id[i] = data2->fam_id[ifam2];
      data->num_mems[i] = data2->num_mems[ifam2];
    }

    data->ind[i] = (struct individual **)our_alloc((ALLOC)data->num_mems[i] * sizeof(struct individual *));

    jind2 = -1;
    for (j = 0; j < data->num_mems[i]; j++){

      data->ind[i][j] = (struct individual *)our_alloc((ALLOC)sizeof(struct individual));
      data->ind[i][j]->a = (SHORT **)our_alloc((ALLOC)num_loci * sizeof(SHORT *));

      for (k = 0; k < num_loci; k++){
    data->ind[i][j]->a[k] = (SHORT *)our_alloc((ALLOC)2 * sizeof(SHORT));
    data->ind[i][j]->a[k][0] = data->ind[i][j]->a[k][1] = 0;
      }

      flag = 0;
      if (i < data1->num_fams && j < data1->num_mems[i]) {
    dat_ind = data1->ind[i][j];
    for (k = 0; k < data1->num_loci; k++){
      data->ind[i][j]->a[k][0] = dat_ind->a[k][0];
      data->ind[i][j]->a[k][1] = dat_ind->a[k][1];
    }
    if (i2 < data2->num_fams){
      for (m = 0; m < data2->num_mems[i2]; m++)
        if (data2->ind[i2][m]->id == data1->ind[i][j]->id) break;

      if (m < data2->num_mems[i2]){
        dat_ind = data2->ind[i2][m];
        flag = 1;
      }
    }
      }
      else if (i >= data1->num_fams) {
    dat_ind = data2->ind[ifam2][j];
    flag = 1;
      }
      else {
    jind2++;
    do {
      flag = 0;
      for (j1 = 0; j1 < data1->num_mems[i]; j1++)
        if (data2->ind[i2][jind2]->id == data1->ind[i][j1]->id){
          jind2++;
          flag = 1;
          break;
        }
    } while (flag);
    dat_ind = data2->ind[i2][jind2];
    flag = 1;
      }
      data->ind[i][j]->id = dat_ind->id;
      data->ind[i][j]->moth_id = dat_ind->moth_id;
      data->ind[i][j]->fath_id = dat_ind->fath_id ;
      data->ind[i][j]->sex = dat_ind->sex;
      if (flag)
    for (k = 0; k < data2->num_loci; k++){
      data->ind[i][j]->a[loc2_pos[k]][0] = dat_ind->a[k][0];
      data->ind[i][j]->a[loc2_pos[k]][1] = dat_ind->a[k][1];
    }

    }
  }

  data->num_fams = nf;
  return (data);
}
