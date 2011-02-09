#include <cpgplot.h>

void plot(float results[], int nchan, int nsamp, int nbits, int xStart, int xEnd, char *plot_device, char *filename)
{
	/*int z;
	for (z = 0; z < 256; z++) {
		//printf("%f\n", results[512 + z]);
		printf("%f\n", results[z]);
	}*/


	int i=0;
	float minx,maxx,miny,maxy;
	float tr[6];

	minx = 0;
	maxx = nchan;
	miny = 0;
	maxy = nsamp;

	tr[0] = minx;
	tr[1] = 0.0;
	tr[2] = (maxx-minx)/nchan;
	tr[3] = miny;
	tr[4] = (maxy-miny)/nsamp;
	tr[5] = 0.0;

    if (cpgopen(plot_device) < 0) {
        printf("Error: Could not open plot device\n");
		cpgbeg(0,"?",1,1);
        //exit(1);
	}

		cpgask(0);

	int k = 20;
	float x = 1.0 / (float) k;
	float r = 0.0;
	float g = 0.0;
	float b  = 0.0;
	int j = 21;

	for( i=0; i<k; i++ ) {
		cpgscr( j, r, g, b );
		r = r + x;
		j++;
	}

	r = 1.0;
	for( i=0; i<k; i++ ) {
		cpgscr( j, r, g, b );
		g = g + x;
		j++;
	}

	g = 1.0;
	for( i=0; i<k; i++ ) {
		cpgscr( j, r, g, b );
		b = b + x;
		j++;
	}

	cpgenv(xStart, xEnd, 1, nchan, 0, 1);

	//cpgenv(550,850,1,nchan,0,1);
	//cpgenv( 1550.0, 1850.0, 1, nchan, 0, 1 );  /* For s080311_044427.sf */
	//cpgenv( 1550.0, 1850.0, 0.0, 511.0, 0, 0 );  /* For s080311_044427.sf */
	//cpgenv( 1700.0, 2000.0, 0.0, 511.0, 0, 0 );  /* For s080311_044903.sf */
	//cpgenv( 400.0, 1000.0, 0.0, 511.0, 0, 0 );  /* For s080313_041333.sf */

	cpgscir(20, 80);

	if (nbits == 8) {
		cpgwedg( "RI", 2.0, 2.0, -8, 18, " " );
		cpgimag(&results[0], nchan, nsamp, 1, nchan, 1, nsamp, -8, 18, tr);
	} else if (nbits == 4) {
		cpgwedg( "RI", 2.0, 2.0, -8, 8, " " );
		cpgimag(&results[0], nchan, nsamp, 1, nchan, 1, nsamp, -8, 8, tr);
	} else if (nbits == 2) {
		cpgwedg( "RI", 2.0, 2.0, -3, 3, " " );
		cpgimag(&results[0], nchan, nsamp, 1, nchan, 1, nsamp, -3, 3, tr);
	} else {
		cpgwedg( "RI", 2.0, 2.0, -1, 1, " " );
		cpgimag(&results[0], nchan, nsamp, 1, nchan, 1, nsamp, -1, 1, tr);
	}

	cpglab("Time  (200 usec samples)", "Frequency (Channels)", filename);
	cpgend();
}

void drawHistogram(float values[], int nbin, float rms, float mean, int scaled, char* filename, char* plot_device)
{
	/*if (nbin == 20) {
		nbin = 24;
	}*/

	float xVals[nbin];
	float xVal;

	if (scaled) {
	  xVal = -(nbin/2 - 0.5);
	} else {
		//xVal = -(nbin/2 + 0.5);
		xVal = -(nbin/2);
	}

	int i;

	if (nbin == 4 && scaled) {
		xVals[0] = -2.0;
		xVals[1] = -0.5;
		xVals[2] = 0.5;
		xVals[3] = 2.0;

	} else {
		if (nbin == 96) {
			xVal = -32;
		}

		for (i = 0; i < nbin; i++) {
			xVals[i] = xVal++;
		}
	}

	int maxY = values[0];
	for (i = 1; i < nbin; i++) {
		if (maxY < values[i]) {
			maxY = values[i];
		}
	}

	if (nbin == 4 && scaled) {
		for (i = 0; i < nbin; i++) {
			printf("Samples at %.1f: %.0f\n", xVals[i], values[i]);
		}
	/*} else if (nbin == 24) {
		float z = -8.0;
		for (i = 0; i < nbin; i++) {
			printf("Samples at %.1f: %.0f\n", z++, values[i]);
		}
	} else {*/
	} else if (nbin == 96) {
		float z = -31.5;
		for (i = 0; i < nbin; i++) {
			printf("Samples at %.1f: %.0f\n", z++, values[i]);
		}
	} else {
		if (scaled) {
			float z = -(nbin/2) + 0.5;
			for (i = 0; i < nbin; i++) {
				printf("Samples at %.1f: %.0f\n", z++, values[i]);
			}
		} else {
			float z = -(nbin/2);
			for (i = 0; i < nbin; i++) {
				printf("Samples at %.1f: %.0f\n", z++, values[i]);
			}
		}
	}


    if (cpgopen(plot_device) < 0) {
        printf("Error: Could not open specified plot device\n");
		cpgbeg(0,"?",1,1);
	}

	cpgask(0);
	cpgenv(-(nbin/2 + 2), nbin/2 + 2, 0, (maxY + 0.1*maxY), 0, 0);

	float tempXvals[130];
	float tempYvals[130];

	for (i = 0; i < 130; i++) {
		tempXvals[i] = 0.0;
		tempYvals[i] = 0.0;
	}

	// set up the first x value

	tempXvals[0] = xVals[0] - 1;

	// set up the last x value
	tempXvals[nbin+1] = xVals[nbin-1] + 1;

	for (i = 1; i <= nbin; i++) {
		tempXvals[i] = xVals[i-1];
		tempYvals[i] = values[i-1];
	}

	/*for (i = 0; i < nbin + 2; i++) {
		printf("%d x: %f y: %f x: %f y: %f\n", i, tempXvals[i], tempYvals[i], xVals[i], values[i]);
	}*/

	cpgbin(nbin + 2, tempXvals, tempYvals, 1);
	//cpgbin(nbin, xVals, values, 1);

	//void cpghist(int n, const float *data, float datmin, float datmax, \
	// int nbin, int pgflag);
	//cpghist(4096, values, -8.0, 7.0, nbin, 1);

	char title[100];
	if (scaled) {
		sprintf(title, "File: %s (scaled)", filename);
	} else {
		sprintf(title, "File: %s (raw)", filename);
	}

	cpglab("Samples", "Number in Each Sample", title);

	float xPoints[2] = {mean, mean};
	float yPoints[2];

	yPoints[0] = 0.0;
	yPoints[1] = maxY + 0.1*maxY;

	cpgsci(2);
	cpgline(2, xPoints, yPoints);

	char meanText[100];
	sprintf(meanText, "mean: %.3f", mean);
	cpgtext(xPoints[1] + 0.2, yPoints[1] - yPoints[1]/20.0, meanText);

	printf("\nmean: %f\n", mean);
	printf("rms: %f\n\n", rms);

	float xCurve[100];
	float yCurve[100];
	float xCurveVal = tempXvals[0];

	//printf("tempXvals[0]: %f\n", tempXvals[0]);

	for (i = 0; i < 100; i++) {
		xCurve[i] = xCurveVal;

		/*if (scaled) {
			if (nbin == 2) {
				yCurve[i] = maxY * exp(-pow((xCurveVal-(mean)), 2) / (2*rms*rms));
			} else {
				yCurve[i] = maxY * exp(-pow((xCurveVal-(mean)), 2) / (2*rms*rms));
			}
		} else {*/
			yCurve[i] = maxY * exp(-pow((xCurveVal-(mean)), 2) / (2*rms*rms));
		//}

		/*if (nbin == 24) {
			xCurveVal += 64.0/100.0;
		} else {*/
			//xCurveVal += ((fabs(tempXvals[0]) * 2.0) + 1) / 100.0;
			xCurveVal += 72.0 / 100.0;
		//}
	}

	xPoints[0] = mean - rms;
	xPoints[1] = xPoints[0];
	cpgsci(7);
	cpgline(2, xPoints, yPoints);

	char rmsText[100];
	sprintf(rmsText, "rms: %.3f", rms);
	cpgtext(mean + rms + 0.2, yPoints[1]/2, rmsText);

	xPoints[0] = mean + rms;
	xPoints[1] = xPoints[0];
	cpgline(2, xPoints, yPoints);

	cpgsci(3);
	//cpgline(100, xCurve, yCurve);
	cpgend();
}
