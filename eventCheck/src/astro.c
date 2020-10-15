/***************************************************
 *    This routine is suitable for time conversion *
 *    among rabbittime, local time and MJD.        *
 *    This program is only for analysis of WFCTA.  *
 *    If you find any bug please send email to     *
 *    yinlq@ihep.ac.cn and youzhiyong@ihep.ac.cn.  *
 *                              --yinlq 20200313   *
 **************************************************/


#include <string.h>
#include <math.h>
#include <time.h>
#include "astro.h"

/* All routines in this library use MODIFIED julian day */

/* modified julian day for 00:00:00 Jan 1, 1970 */
#define JD0 2440587.5
/* convert modified julian day to time_t */
#define timejd(jd) ((time_t)((jd - JD0) * 86.4e3))

#define MJD19700101 40587
#define TAI2UTC 37

static int local_time = 0;

/**************************************************************/
/* Set time conversion mode (0 = UTC, non-zero = LOCAL)       */
/**************************************************************/
int astroLocal(int local)
{
  int old = local_time;
  local_time = local;
  return old;
}

int astro_local_(int *local)
{
  return astroLocal(*local);
}

/**************************************************************/
/* Convert date/time to Julian day                              */
/**************************************************************/
double lt2mjd(int year, int month, int day, int hour, int min, int sec)
{
  double jd; int a, b;
  if (month <= 2) year -= 1, month += 12;
  a = year / 100;
  b = 2 - a + (a / 4);
  jd = (double)((int)(365.25 * year) + (int)(30.6001 * (month + 1)) + day + b)
    + 1720994.5 + (((sec / 60.0) + min) / 60.0 + hour) / 24.0;
  if (local_time)
  { /* correct julian day if local time */
    time_t s = timejd(jd);
    struct tm *t = gmtime(&s);
    min = t->tm_hour * 60 + t->tm_min;
    t = localtime(&s);
    min -= t->tm_hour * 60 + t->tm_min;
    if (min > 720) min -= 720;
    else if (min <= -720) min += 720;
    jd += min / 1440.0;
  }
  return jd;
}

/**************************************************************/
/* Convert Julian day to date/time                            */
/**************************************************************/
void mjd2lt(double jd,
  int *year, int *month, int *day, int *hour, int *min, int *sec, char *zone)
{
  time_t s = timejd(jd);
  struct tm *t = (local_time ? localtime : gmtime)(&s);
  if (year != NULL) *year = t->tm_year + 1900;
  if (month != NULL) *month = t->tm_mon + 1;
  if (day != NULL) *day = t->tm_mday;
  if (hour != NULL) *hour = t->tm_hour;
  if (min != NULL) *min = t->tm_min;
  if (sec != NULL) *sec = t->tm_sec;
  if (zone != NULL) strncpy(zone, local_time ? tzname[(t->tm_isdst > 0)] : "UTC", 4);
}

/**************************************************************/
/* Convert rabbit time (TAI) to date/time                     */
/**************************************************************/
void rbtime2lt(int64_t rabbitTime, double rabbittime,
  int *year, int *month, int *day, int *hour, int *min, int *sec)
{
  time_t rawtime;
  time_t utc = rabbitTime-TAI2UTC;
  struct tm *info;
  info = localtime(&utc);
  if (year != NULL) *year = 1900 + info->tm_year;
  if (month != NULL) *month = info->tm_mon+1;
  if (day != NULL) *day = info->tm_mday;
  if (hour != NULL) *hour = info->tm_hour;
  if (min != NULL) *min = info->tm_min;
  if (sec != NULL) *sec = info->tm_sec;
}

/**************************************************************/
/* Convert date/time to rabbit time (TAI)                     */
/**************************************************************/
int64_t lt2rbTime ( int iy, int im, int id, int ih, int in, int is )
{
   long iyL, imL;
   double djm;
   int j;

   /* Month lengths in days */
   static int mtab[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
   /* Validate year */
   if ( iy < -4699 ) { j = 1; return 1; }
   /* Validate month */
   if ( ( im < 1 ) || ( im > 12 ) ) { j = 2; return 1; }
   /* Allow for leap year */
   mtab[1] = ( ( ( iy % 4 ) == 0 ) &&
                   ( ( ( iy % 100 ) != 0 ) || ( ( iy % 400 ) == 0 ) ) ) ?
           29 : 28;

   /* Validate day */
   j = ( id < 1 || id > mtab[im-1] ) ? 3 : 0;
   /* Lengthen year and month numbers to avoid overflow */
   iyL = (long) iy;
   imL = (long) im;
   /* Perform the conversion */
   djm = (double)
           ( ( 1461L * ( iyL - ( 12L - imL ) / 10L + 4712L ) ) / 4L
                 + ( 306L * ( ( imL + 9L ) % 12L ) + 5L ) / 10L
                 - ( 3L * ( ( iyL - ( 12L - imL ) / 10L + 4900L ) / 100L ) ) / 4L
                 + (long) id - 2399904L );

   //long double rb_time;
   //(mjd - MJD19700101)*86400 + TAI2UTC = rabbitTime + rabbittime*20/1000000000.;
   int64_t rabbitTime = (int64_t)((djm - MJD19700101)*86400 + TAI2UTC + ih*3600 + in*60 + is - 8*3600);

   return rabbitTime;
}

/**************************************************************/
/* Convert rabbit time (TAI) from Modified Julian day         */
/**************************************************************/
void mjd2rbtime(double mjd, int64_t *rbtime, int64_t *rabbitTime, double *rabbittime)
{
  int64_t temp1 = (int64_t)((( mjd - MJD19700101 ) * 86400 + TAI2UTC) * 1000000000. );

  int64_t temp2 = (int64_t)(temp1 / 1000000000);
  double temp3 = (temp1 / 1000000000. - temp2 )  * 1000000000. / 20;

  *rbtime = temp1;
  *rabbitTime = temp2;
  *rabbittime = temp3;
}

/**************************************************************/
/* Convert rabbit time (TAI) to Modified Julian day           */
/**************************************************************/
double rbtime2mjd(int64_t rabbitTime, double rabbittime)
{
  double mjd;
  mjd = MJD19700101 + (rabbitTime + rabbittime*20/1000000000. - TAI2UTC)/86400;
  return mjd;
}

