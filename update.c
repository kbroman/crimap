#if vms
#include stdio
#else
#include <stdio.h>
#endif

#include "defs.h"
#include "var1.h"

/* this routine changes theta given nrecs and recs   */

double change_theta(ad_num_its)
     SHORT *ad_num_its;                          /* num_its */
{
  double log_like, prev_log_like;
  double  get_log_like();
  SHORT num_its;

  num_its = 1;
  log_like = get_log_like(recs, nrecs);

  do {
    copy(recs, recs_temp);
    update();

    /***  find log_like every 10 iterations ***/
    if( !(num_its % 10)){
      prev_log_like = log_like;
      log_like = get_log_like(recs, nrecs);
    }
    num_its ++;
  } while( num_its < 11 || (log_like - prev_log_like) > TOL );

  if(log_like < prev_log_like - MAX_DECREASE){
    printf("\n\nLOG_LIKE DECREASED IN CHANGE_THETA  DIFF = %f\n\n",
       log_like - prev_log_like);
  }

  *ad_num_its = num_its;
  return(log_like);
}

make_num_mei_split()
{
     SHORT n, num_types, i, j, k;
     double temp;

     n = theta->n;
     num_types = theta->num_types;

     copy(num_mei, num_mei_split);

     /* now begin splitting num_mei to build num_mei_split */
     for(k = 0; k < num_types; k++)
        for(i = 0; i < n; i++)
            for(j = 1; j < n-i; j++){
                temp = num_mei_split->data[k][i][j];
                num_mei_split->data[k][i][0] += temp;
                num_mei_split->data[k][i+1][j-1] += temp;
            }
}

make_theta_1_t()
{
     SHORT n, num_types, i, j, k;
     double temp;

     n = theta->n;
     num_types = theta->num_types;
     for(k = 0; k < num_types; k++)
        for(i = 0; i < n; i++)
            for(j = 0; j < n-i; j++){
                temp = theta->data[k][i][j];
                theta_1_t->data[k][i][j] = temp / (1 - temp);
            }
}


update()
{

    SHORT i, j, k, n;
    double nms1, nms3, tm, prx1, pnrx1, r1;
    double num_recs, num_nrecs, t1, t2, t3;
    double **temp_step, **rec_t, **thet_t;
    double *nms;

    n = theta->n;
    for(k = 0; k < theta->num_types; k++) {
      rec_t = recs_temp->data[k];
      thet_t = theta->data[k];
      for(i = 0; i < n; i++) {
        nms = num_mei_split->data[k][i];
        nms1 = nms[0];
        if (nms1) {
          t1 = thet_t[i][0];
          if (t1) {
             r1 = rec_t[i][0];
             for (j = n - i - 1; j > 0; j--) {
              nms3 = nms[j];
              if (nms3) {
                t3 = thet_t[i][j];
                    t2 = thet_t[i + 1][j - 1];
                num_recs = rec_t[i][j];
                prx1 = num_recs * t1 * (1 - t2) / t3;
                pnrx1 = (nms3 - num_recs) * t1 * t2 / (1 - t3);
                r1 += prx1 + pnrx1;
                rec_t[i + 1][j - 1] += num_recs + pnrx1 - prx1;
              }
             }
         if(!FIXED_INTERVALS[k][i]){
           tm = r1 / nms1;
           if (tm < 0.0001) tm = 0.0;
           else if (tm > 0.5) tm = 0.5;
           thet_t[i][0] = tm;
         }
         tm = thet_t[i][0];
         t2 = 1-tm-tm;
         for(j = 0; j < i; j++)
           thet_t[j][i-j] = tm+t2*thet_t[j][i-j-1];
          }
          else {
             for (j = n - i - 1; j > 0; j--) rec_t[i + 1][j - 1] += rec_t[i][j];
         for(j = 0; j < i; j++)
           thet_t[j][i-j] = thet_t[j][i-j-1];
          }
        }
        else {
      if(!FIXED_INTERVALS[k][i]) {
        thet_t[i][0] = 0.0;
         for(j = 0; j < i; j++)
           thet_t[j][i-j] = thet_t[j][i-j-1];
      }
      else {
        tm = thet_t[i][0];
        t2 = 1-tm-tm;
        for(j = 0; j < i; j++)
          thet_t[j][i-j] = tm+t2*thet_t[j][i-j-1];
      }
    }
     }
   }
/*
   fill_theta();
*/
}

/* this routine computes theta(i,j)  (the recombination fraction
   between locus i and locus i+j+1) given the recombination fractions
   theta(i,0)                           */

fill_theta()
{
  SHORT i, j, k, n;
  double **data;
  double t1,t2,t3;

  n = theta->n;
  for(k = 0; k < theta->num_types; k++){
    data = theta->data[k];
    for(i = 0; i < n; i++){
      t1 = data[i][0];
      for(j = 1; j < n-i; j++){
    t2 = data[i+j][0];
    t3 = t1*t2;
    data[i][j] = t1 = t1+t2-t3-t3;
      }
    }
  }
}
