/* =================================================================
 * SUBJECT: Prints a statusbar in the R console
 * AUTHOR: Tobias Schoch, September 12, 2011
 * LICENSE: GPL > 2
 * ================================================================= */

#include <R.h>

void statusbar( int *pos, int *total )
{
   /*some definitions*/
   unsigned short int diff, i, j, first;
   unsigned short int thispos;
   unsigned short int lastpos;
   const double positions = 50.0;
   double u, v;
   /*definitions and castings */
   double dpos;
   dpos = *pos * 1.0;
   double dtotal;
   dtotal = *total * 1.0;
   /* */
   /*print the status bar only at the first iteration*/
   if ( *pos == 1 ){
      Rprintf("0   10   20   30   40   50   60   70   80   90   100\n");
      Rprintf("|----|----|----|----|----|----|----|----|----|----|\n");
      v = positions / dtotal;
      first = (int)v;
      for (j = 1; j <= first; j++){
	 Rprintf("*");
      }
   }
   else
   {
      /*compute relative position in status bar*/
      u = ( dpos / dtotal ) * positions;
      thispos = (int)u;
      u = ( (dpos - 1) / dtotal ) * positions;
      lastpos = (int)u;
      diff = thispos - lastpos;
      for ( i = 1; i <= diff; i++ ){
	 Rprintf("*");
      }
      if ( *pos == *total)
      {
	 Rprintf("*\nDONE\n");
      }
   }
}
