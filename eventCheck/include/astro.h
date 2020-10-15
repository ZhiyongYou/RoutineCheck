/***************************************************
 *    This routine is suitable for time conversion *
 *    among rabbittime, local time and MJD.        *
 *    This program is only for analysis of WFCTA.  *
 *    If you find any bug please send email to     *
 *    yinlq@ihep.ac.cn and youzhiyong@ihep.ac.cn.  *
 *    				--yinlq 20200313   *
 **************************************************/


#include <stdint.h>
#ifndef M_2PI
#define M_2PI (M_PI + M_PI)
#endif

#ifndef NULL
#define NULL ((void *)0)
#endif

/**************************************************************/
/* Convert date/time to/from Julian day                       */
/**************************************************************/
int astro_local(int local);

double lt2mjd(int year, int month, int day, int hour, int min, int sec);

void mjd2lt(double jd,
  int *year, int *month, int *day, int *hour, int *min, int *sec, char *zone);

/**************************************************************/
/* Convert rabbit time (TAI) to/from date                     */
/**************************************************************/
void rbtime2lt(int64_t rabbitTime, double rabbittime, 
  int *year, int *month, int *day, int *hour, int *min, int *sec );

int64_t lt2rbTime ( int iy, int im, int id, int ih, int in, int is );

/**************************************************************/
/* Convert rabbit time (TAI) to/from Modified Julian day      */
/**************************************************************/
void mjd2rbtime(double mjd, int64_t *rbtime, int64_t *rabbitTime, double *rabbittime);

double rbtime2mjd(int64_t rabbitTime, double rabbittime);

