/* Simple function to calculate soil temperature from air temperature */
/* This function will only run after a lag of 48 hours */

#include <math.h>

double stemp(int hr, double *atemp, int idx, double acoef){

	int lag = 48;
	int i = idx - lag;
	int j = 0;
	double stmpm = 0.0, stmpmax = -60.0, stmpmin = 60.0, amp = 0.0;
	double ans;
	double hour = hr;
	/* double acoef = 8.0; */
	double offset = 0.0;
	double bcoef = 3.66;

	for(i; j<lag; i++){
		stmpm += atemp[i] / 48.0;
		if(atemp[i] > stmpmax) stmpmax = atemp[i];
		if(atemp[i] < stmpmin) stmpmin = atemp[i];
		j += 1;
	}
	amp = (stmpmax - stmpmin) / acoef;
	ans = stmpm + amp * -1 * sin((hour - offset)/bcoef);

	return(ans);

}
 
