// gcc -g -o delay_template4 delay_template4.c -lcfitsio -lcpgplot -lpgplot -L/usr/X11R6/lib  -lX11 -lpng -lg2c

/*
 *  Create PAT templates for pulsar backend delay measurements
 *  RNM 30 June 2008
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cpgplot.h"
#include "fitsio.h"

#define MAX 64
#define MAX_INT16 16383.0

typedef int16_t int16;

void usage(int status){
  printf("\nComputes PAT timing template for CalDelay measurements\n");
  printf("Usage:\n");
  printf("delay_template [-h] [-a config_file] [-b impulse_file] \n"); 
  printf("   [-g plotdev] [-i] [-j] [-k] [-z phs1 phs2] [-v] backend\n");
  printf("-a: Over-ride default config list file\n"); 
  printf("-b: Over-ride default impulse response file\n");
  printf("-t: Trigger repetition rate (Hz) (def=1000)\n"); 
  printf("-c: Compression factor for impulse function (def=2048)\n"); 
  printf("-g: Plot device (def=/xs)\n"); 
  printf("-i: Plot impulse response\n");
  printf("-j: Plot template profile\n");
  printf("-k: Plot convolved template profile\n");
  printf("-f: Plot final profile\n");
  printf("-z: Zoom in pulse phase\n");
  printf("-v: Verbose output\n");
  printf("-r: RMS of noise to be added\n");
  printf("Backends: WBC, PDFB1, PDFB2, PDFB3, CPSR2\n");
  printf("\n");
  exit(status);
}

void modify_tdim( fitsfile *cfp, int colnum, int naxes, long *axes, int *status );
void remove_unused_tables(fitsfile *fp);
float get_random_number();

int main(int argc, char *argv[]){
  int verbose,i,j,jn,jc,k,kn,kk,knp,knc,knrc,nbn,bw,nch,cnbn[MAX],cbw[MAX],cnch[MAX];
  int icmp,oversamp,nbrst[MAX],nburst,ntrig,kev,plimp,plprf,plprc,plpff;
  int knr,km,ii,kz,kz1,kz2;
  int addrand;
  float impv,imp[65536],imax,cfrq[MAX],cfreq,trig,ttrig,pspn[MAX],psepn;
  float pdcy[MAX],pdcyc,minprd,mnprd[MAX],prf[65537],pp;
  //float tev[512],tton,ttoff,tp,dtimp1,dtimp,tbin,tt,xp[65537],impr[1024];
  float tev[512*512],tton,ttoff,tp,dtimp1,dtimp,tbin,tt,xp[65537],impr[1024];
  float xmin,xmax,ymin,ymax,kstep,prc[65537],pr[65537],impn,zph1,zph2,zph;
  float rand;
  char *s,cfgfile[32],ncfgfile[32],backend[8],nimpfn[32],impfn[32],plotdev[16];
  char key[8],cfg[32],ccfg[MAX][32],tplbl[64];
  FILE *fptr;

  verbose=0;
  strcpy(ncfgfile,"default");
  strcpy(nimpfn,"default");
  strcpy(plotdev,"/xs");
  icmp=4;
  oversamp=32;
  trig=1000.00;   /* Cal trigger repetition rate (Hz) */
  plimp=0;
  plprf=0;
  plprc=0;
  plpff=0;
  zph1=0.0;
  zph2=0.0;
  addrand=0;

  argc--; argv++;
  if(argc <= 1)usage(1);

  while(argc > 1){                /* Get command line inputs */
    if((*argv)[0] == '-'){
      s=argv[0]+1;
      argc--; argv++;
      switch(*s){
      case 'h':
      case '?':
    usage(0);
      case 'a':
    if(sscanf(*argv,"%s",ncfgfile) != 1)usage(1);
    argc--; argv++;
    break;
      case 'b':
    if(sscanf(*argv,"%s",nimpfn) != 1)usage(1);
    argc--; argv++;
    break;
      case 't':
    if(sscanf(*argv,"%f",&trig) != 1)usage(1);
    argc--; argv++;
    break;
      case 'g':
    if(sscanf(*argv,"%s",plotdev) != 1)usage(1);
    argc--; argv++;
    break;
      case 'i':
    plimp=1;
    break;
      case 'j':
    plprf=1;
    break;
      case 'k':
    plprc=1;
    break;
      case 'f':
    plpff=1;
    break;
      case 'r':
    if(sscanf(*argv,"%f",&rand) != 1)usage(1);
    argc--; argv++;
    addrand = 1;
    break;
      case 'z':
    if(sscanf(*argv,"%f",&zph1) != 1)usage(1);
    argc--; argv++;
    if(sscanf(*argv,"%f",&zph2) != 1)usage(1);
    argc--; argv++;
    break;
      case 'v':
    verbose=1;
    break;
      default:
    usage(1);
      }
    }
  } 

  if(sscanf(*argv,"%s",backend) != 1)usage(1);

  /* Set up plot */

  cpgbeg(0,plotdev,1,1);
  cpgsvp(0.15,0.9,0.15,0.9);
  cpgslw(2);
  cpgscf(2);
  cpgsch(1.3);
  if(strstr(plotdev,"xs") != NULL)cpgask(1); 
  cpgpage();

  /* Open config list file */
  if(strstr(ncfgfile,"default") == NULL){
    strcpy(cfgfile,ncfgfile);
  }
  else{
    strcpy(cfgfile,"DelayList.txt");
  }

  if((fptr = fopen(cfgfile,"r")) == NULL){    
    printf("Failed to open %s\n",cfgfile);
    exit(1);
  }

  j=0;
  while(fscanf(fptr,"%s %s %f %d %d %d %f %f %d %f",
           key,cfg,&minprd,&nbn,&bw,&nch,&cfreq,&psepn,&nburst,&pdcyc) != EOF){
    if(strstr(key,backend) != NULL){
      strcpy(ccfg[j],cfg);
      mnprd[j]=minprd;
      cnbn[j]=nbn;
      cbw[j]=bw;
      cnch[j]=nch;
      cfrq[j]=cfreq;
      pspn[j]=psepn;
      nbrst[j]=nburst;
      pdcy[j]=pdcyc;
      if(verbose==1)printf("%3d %-20s %8.4f %6d %6d %6d %10.4f %8.3f %4d %6.3f\n",
        j,ccfg[j],mnprd[j],cnbn[j],cbw[j],cnch[j],cfrq[j],pspn[j],nbrst[j],pdcy[j]);
      j++;
    }
  }
  fclose(fptr);
  jn=j;

  /* Get power impulse response function */

  if(strstr(nimpfn,"default") == NULL){
    strcpy(impfn,nimpfn);
  }
  else{
    strcpy(impfn,"Impulse.txt");
  }
  if((fptr = fopen(impfn,"r")) == NULL){    
    printf("Failed to open %s\n",impfn);
    exit(1);
  }
  printf("\n%s opened\n",impfn);

  /* Impulse time step (microsec) for 1 MHz chan bw */
  dtimp1=(float)icmp/8192.0;

  j=0;
  k=0;
  while(fscanf(fptr,"%f",&impv) != EOF){
    if(j%icmp == 0){
      imp[k]=impv*impv;
      k++;
    }
    j++;
  }
  kn=k;
  fclose(fptr);
  /* Scale to max 1.0 */
  imax=0.0;
  for(k=0;k<kn;k++){
    if(imp[k]>imax)imax=imp[k];
  }
  for(k=0;k<kn;k++){
    imp[k]/=imax;
    xp[k]=((float)k-0.5*kn)*dtimp1;
  }

  /* Plot impulse response */
  if(plimp==1){
    xmin=-0.52*kn*dtimp1;
    xmax=-xmin;
    ymin=-0.05;
    ymax=1.05;
    cpgswin(xmin,xmax,ymin,ymax);
    cpgbox("BCNST",1.0,5,"BCNST",0.5,5);
    cpglab("Time (microsec)","Normalised Response",
                "PDFB3 Impulse Response: ChBW = 1 MHz");
    cpgline(kn,xp,imp);
  }
 
  /* Begin config loop */

  for(jc=0;jc<jn;jc++){

    cfreq=cfrq[jc];     /* Cal frequency [Hz] */
    psepn=pspn[jc];     /* Pulse separation (microsec) */
    nburst=nbrst[jc];   /* Nr of pulses per trigger */
    pdcyc=pdcy[jc];     /* Pulse duty cycle */

    printf("Config: %d %-20s %8.3f %10.4f %8.3f %4d %8.3f\n",
       jc,ccfg[jc],mnprd[jc],cfreq,psepn,nburst,pdcyc);

    /* Form cal profile */
    
    /* Assume leading edge of first pulse is coincident with leading
       edge of first bin. In fact, it is delayed by 400 +/- 5 ns wrt the 
       trigger pulse */
    
    /* Set event list (in microsec from bin 0) */
    k=0;
    ntrig=trig/cfreq; 
    for(i=0; i<ntrig; i++){
      ttrig=1.e6*i/trig;
      for(j=0; j<nburst; j++){
    tton=j*psepn;
    ttoff=tton + pdcyc*psepn;
    tev[k]=ttrig+tton;
    if(verbose==1)printf("tev: %d %f\n",k,tev[k]);
    k++;
    tev[k]=ttrig+ttoff;
    if(verbose==1)printf("tev: %d %f\n",k,tev[k]);
    k++;
      }
    }

    tev[k]=1.e6/cfreq;
    if(verbose==1)printf("tev: %d %f\n",k,tev[k]);
    kev=k;
    k++;
    tev[k]=tev[k-1]+psepn;
    if(verbose==1)printf("tev: %d %f\n",k,tev[k]);
  
    nbn=oversamp*cnbn[jc];
    tp=0.;
    k=1;
    pp=1.0;
    tbin=1.e6/cfreq/nbn;

    if(verbose==1)printf("%d %f %f\n",nbn,cfreq,tbin);

    /* Fill template profile */

    for(i=0;i<nbn;i++){
      tp+=tbin;
      if(tp <= tev[k]){
	prf[i]=pp;
      }
      else{
	tt=tp-tev[k];
	if(pp > 0.5){
	  prf[i]=tt/tbin;
	  pp=0.0;
	}
	else{
	  if(k<kev){
	    prf[i]=1.0-tt/tbin;
	    pp=1.0;
	  }
	}
	if(verbose==1)printf("%d %d %f %f %f %f %f\n",i,k,tev[k],tp,tt,pp,prf[i]);
	k++;
      }
    }

    /* Plot profile */

    if(plprf == 1){
      ymin=-0.05;
      ymax=1.05;
      if(zph1==0.0 && zph2==0.0){
    for(k=0;k<=nbn;k++){
      xp[k]=(float)k/nbn;
      pr[k]=prf[k];
    }
    pr[nbn]=pr[0];
    knp=nbn+1;
    xmin=-0.05;
    xmax=1.05;
      }
      else{
    kz1=fabs(zph1)*nbn+0.5;
    if(zph1<0.0)kz1=-kz1;
    zph1=(float)kz1/nbn;
    kz2=zph2*nbn+0.5;
    zph2=(float)kz2/nbn;
    zph=zph2-zph1;
    for(kz=kz1;kz<=kz2;kz++){
      kk=kz-kz1;
      k=kz;
      if(k < 0)k+=nbn;
      pr[kk]=prf[k];
      xp[kk]=zph1+(float)kk/nbn;
      if(verbose==1 && (xp[kk]>-0.005 && xp[kk]<0.005))
        printf("%d %d %d %f %f\n",kz,k,kk,xp[kk],pr[kk]);
    }
    knp=kz2-kz1+1;
    xmin=zph1-0.05*zph;
    xmax=zph2+0.05*zph;
      }

      cpgswin(xmin,xmax,ymin,ymax);
      cpgbox("BCNST",0.,0,"BCNST",0.5,5);
      sprintf(tplbl,"Delay Template: %s",ccfg[jc]);
      cpglab("Pulse Phase","Intensity",tplbl);
      cpgline(knp,xp,pr);
      cpgpage();
    }
    
  /* Resample impulse fn at tbin (profile step) */

    dtimp=dtimp1/cbw[jc]*cnch[jc];
    kstep=tbin/dtimp;
    knr=kn/kstep + 0.5;
    knr=knr/2*2 +1;
    knc=kn/2;
    knrc=knr/2;
    impr[knrc]=imp[knc];
    for(k=1;k<=knr/2;k++){
      kk=k*kstep + 0.5;
      impr[knrc-k]=imp[knc-kk];
      impr[knrc+k]=imp[knc+kk];
    }    

    if(verbose==1){
      printf("%d %f %f %d %d %f\n",nbn,tbin,dtimp1,kn,knr,kstep);
      for(k=0;k<knr;k++){
    printf("%d %f\n",k,impr[k]);
      }
    }
    
    /* Convolve template with impulse fn */
    km=knr/2;
    for(i=0;i<nbn;i++){
      prc[i]=0.0;
    }
    for(i=0;i<nbn;i++){
      for(k=0;k<knr;k++){


    ii=i+k-km;
    if(ii < 0)ii+=nbn;
    if(ii >= nbn)ii-=nbn;


    prc[ii]+=prf[i]*impr[k];
      }
    }

    /* Renormalise */
        
    impn=0.;
    for(k=0;k<knr;k++){
      impn+=impr[k];
    }
    for(i=0;i<nbn;i++){
      prc[i] /= impn;
    }

    /* Plot convolved profile */

    if(plprc == 1){
      ymin=-0.05;
      ymax=1.05;
      if(zph1==0.0 && zph2==0.0){
    for(k=0;k<=nbn;k++){
      xp[k]=(float)k/nbn;
      pr[k]=prc[k];
      //printf("prc[%d]: %f, pr[]: %f\n", k, prc[k], pr[k]);
    }
    pr[nbn]=pr[0];
    knp=nbn+1;
    xmin=-0.05;
    xmax=1.05;
      }
      else{
    kz1=fabs(zph1)*nbn+0.5;
    if(zph1<0.0)kz1=-kz1;
    zph1=(float)kz1/nbn;
    kz2=zph2*nbn+0.5;
    zph2=(float)kz2/nbn;
    zph=zph2-zph1;
    for(kz=kz1;kz<=kz2;kz++){
      kk=kz-kz1;
      k=kz;
      if(k < 0)k+=nbn;
      pr[kk]=prc[k];
      xp[kk]=zph1+(float)kk/nbn;
      if(verbose==1 && (xp[kk]>-0.005 && xp[kk]<0.005))
        printf("%d %d %d %f %f\n",kz,k,kk,xp[kk],pr[kk]);
    }
    knp=kz2-kz1+1;
    xmin=zph1-0.05*zph;
    xmax=zph2+0.05*zph;
      }

      cpgswin(xmin,xmax,ymin,ymax);
      cpgbox("BCNST",0.,0,"BCNST",0.5,5);
      sprintf(tplbl,"Convolved Delay Template: %s",ccfg[jc]);
      cpglab("Pulse Phase","Intensity",tplbl);
      cpgline(knp,xp,pr);
      cpgpage();
    }

    /* Convert to original number of bins */
    nbn=cnbn[jc];
    for(i=0;i<nbn;i++){
      prf[i]=0.0;
      for(k=0;k<oversamp;k++){
    prf[i]+=prc[i*oversamp+k];
      }
      prf[i]/=oversamp;
    }

    if(verbose==2){
      for(i=0;i<nbn;i++){
    printf("%d %f\n",i,prf[i]);
      }
    }

    if(plpff == 1){
      ymin=-0.05;
      ymax=1.05;
      if(zph1==0.0 && zph2==0.0){
    for(k=0;k<=nbn;k++){
      xp[k]=(float)k/nbn;
      pr[k]=prf[k];
    }
    pr[nbn]=pr[0];
    knp=nbn+1;
    xmin=-0.05;
    xmax=1.05;
      }
      else{
    kz1=fabs(zph1)*nbn+0.5;
    if(zph1<0.0)kz1=-kz1;
    zph1=(float)kz1/nbn;
    kz2=zph2*nbn+0.5;
    zph2=(float)kz2/nbn;
    zph=zph2-zph1;
    for(kz=kz1;kz<=kz2;kz++){
      kk=kz-kz1;
      k=kz;
      if(k < 0)k+=nbn;
      pr[kk]=prf[k];
      xp[kk]=zph1+(float)kk/nbn;
    }
    knp=kz2-kz1+1;
    xmin=zph1-0.05*zph;
    xmax=zph2+0.05*zph;
      }

      cpgswin(xmin,xmax,ymin,ymax);
      cpgbox("BCNST",0.,0,"BCNST",0.5,5);
      sprintf(tplbl,"Final Delay Template: %s",ccfg[jc]);
      cpglab("Pulse Phase","Intensity",tplbl);
      cpgpt(knp,xp,pr,17);
      cpgline(knp,xp,pr);
      cpgpage();
    }

    // make psrfits archive here

    int status = 0;
    char filename[256];

    sprintf(filename, "!%s.tml(/pulsar/psr/linux/share/psrheader.fits)", ccfg[jc]);

    fitsfile *fp;
    fits_create_file(&fp, filename, &status);
    fits_report_error(stderr, status);

    fits_movabs_hdu(fp, 1, NULL, &status );
    fits_report_error(stderr, status);
    fits_update_key(fp, TSTRING, "SRC_NAME", "CalDelay", NULL, &status);
    fits_report_error(stderr, status);

    fits_update_key(fp, TSTRING, "OBS_MODE", "CAL", NULL, &status);
    fits_report_error(stderr, status);

    double obsfreq = 1400.0;
    fits_update_key(fp, TDOUBLE, "OBSFREQ", &obsfreq, NULL, &status);
    fits_report_error(stderr, status);

    fits_update_key(fp, TSTRING, "CAL_MODE", "FREQ", NULL, &status);
    fits_report_error(stderr, status);

    double cal_freq = 200;
    fits_update_key(fp, TDOUBLE, "CAL_FREQ", &cal_freq, NULL, &status);
    fits_report_error(stderr, status);

    double cal_dcyc = 0.5;
    fits_update_key(fp, TDOUBLE, "CAL_DCYC", &cal_dcyc, NULL, &status);
    fits_report_error(stderr, status);

    double cal_phs = 0.25;
    fits_update_key(fp, TDOUBLE, "CAL_PHS", &cal_phs, NULL, &status);
    fits_report_error(stderr, status);

    int itemp = 53819;
    fits_update_key(fp, TINT, "STT_IMJD", &itemp, NULL, &status );
    fits_report_error(stderr, status);

    itemp = 19440;
    fits_update_key(fp, TINT, "STT_SMJD", &itemp, NULL, &status );
    fits_report_error(stderr, status);

    double offs = 0.00123;

    fits_update_key(fp, TDOUBLE, "STT_OFFS", &offs, NULL, &status );
    fits_report_error(stderr, status);

    fits_update_key(fp, TSTRING, "FD_POLN", "LIN", NULL, &status);
    fits_report_error(stderr, status);

    fits_update_key(fp, TSTRING, "COORD_MD", "J2000", NULL, &status);
    fits_report_error(stderr, status);

    fits_movnam_hdu (fp, BINARY_TBL, "FEEDPAR", 0, &status);
    fits_report_error(stderr, status);

    fits_update_key(fp, TSTRING, "CAL_MTHD", "unknown", NULL, &status);
    fits_report_error(stderr, status);

    fits_movnam_hdu (fp, BINARY_TBL, "SUBINT", 0, &status);
    fits_report_error(stderr, status);

    double period = 2.3;
    fits_update_key(fp, TDOUBLE, "PERIOD", &period, NULL, &status);
    fits_report_error(stderr, status);

    int npol = 1;
    fits_update_key(fp, TINT, "NPOL", &npol, NULL, &status);
    fits_report_error(stderr, status);

    int nbin = nbn;
    fits_update_key(fp, TINT, "NBIN", &nbin, NULL, &status);

    int nchan = 1;
    fits_update_key(fp, TINT, "NCH_FILE", &nchan, NULL, &status);

    long long ltemp = nchan;
    int colnum = 0;

    if (!fits_get_colnum(fp, CASEINSEN, "DAT_FREQ", &colnum, &status))
    	fits_modify_vector_len(fp, colnum, ltemp, &status );

    if (!fits_get_colnum(fp, CASEINSEN, "DAT_WTS", &colnum, &status))
    	fits_modify_vector_len(fp, colnum, ltemp, &status );

    if (!fits_get_colnum(fp, CASEINSEN, "DAT_OFFS", &colnum, &status))
    	fits_modify_vector_len(fp, colnum, ltemp, &status );

    if (!fits_get_colnum(fp, CASEINSEN, "DAT_SCL", &colnum, &status))
    	fits_modify_vector_len(fp, colnum, ltemp, &status );

    ltemp = nbin * nchan;

    if (!fits_get_colnum(fp, CASEINSEN, "DATA", &colnum, &status)) {
    	fits_modify_vector_len(fp, colnum, ltemp, &status );
    	long naxes[3];
    	naxes[0] = nbin;
    	naxes[1] = nchan;
    	naxes[2] = npol;

    	modify_tdim(fp, colnum, 3, naxes, &status);
    }

    float max_short = pow(2.0,15.0) - 1.0;
    float min = prf[0];
    float max = prf[0];

    for (i = 0; i < nbin; ++i) {
    	float random = addrand ? get_random_number() * rand : 0.0;
    	prf[i] += random;
    	/*if (prf[i] < min)
    		min = prf[i];
    	else if (prf[i] > max)
    		max = prf[i];*/

      //printf("prf: %f max: %f\n", prf[i], max);
    }

    //printf("--- max: %f min: %f\n", max, min);

    max = 1.0;
    min = 0.0;

    // hard-coded stuff here
    //float offset = 0.5 * (max + min);
    float offset = 0.5;

    //float scale_factor = (max - min) / max_short;
    float scale_factor = max_short;

    short int binned_data[npol][nchan][nbin];

    for (i = 0; i < npol; ++i) {
    	for (j = 0; j < nchan; ++j) {
    		for (k = 0; k < nbin; ++k) {
    			//binned_data[i][j][k] = (short int)((prf[k] - offset) / scale_factor);
    			binned_data[i][j][k] = (short int)((prf[k] - offset) * scale_factor);
                /*if (prf[k] < 0.5)
                    binned_data[i][j][k] = 0;
                else
                    binned_data[i][j][k] = 1;*/

                //printf("%d\n", binned_data[i][j][k]);
    		}
    	}
    }

    if (!fits_get_colnum(fp, CASEINSEN, "DATA", &colnum, &status))
    	fits_write_col(fp, TSHORT, colnum, 1, 1, nbin, binned_data, &status );

    float mean = (max - min)/2.0;
    float dat_scl = mean/MAX_INT16;
    if (!fits_get_colnum(fp, CASEINSEN, "DAT_SCL", &colnum, &status))
    	fits_write_col(fp, TFLOAT, colnum, 1, 1, 1, &dat_scl, &status);

    float dat_offs = mean;
    if (!fits_get_colnum(fp, CASEINSEN, "DAT_OFFS", &colnum, &status))
    	fits_write_col(fp, TFLOAT, colnum, 1, 1, 1, &dat_offs, &status );

    remove_unused_tables(fp);

    fits_close_file(fp, &status);
    fits_report_error(stderr, status);
  }
  cpgend();
  exit(0);
}

void remove_unused_tables(fitsfile *fp)
{
    int status = 0;
    int num_hdus = -1;
    fits_get_num_hdus(fp, &num_hdus, &status);
    int finish = num_hdus;
    int i;

    for (i = 2; i <= finish; ++i) {
    	fits_movabs_hdu(fp, i, NULL, &status);
    	long nrows = -1;
    	fits_get_num_rows(fp, &nrows, &status);

    	if (!nrows) {
    		fits_delete_hdu(fp, &i, &status);
    		fits_report_error(stderr, status);
    		--i;
    		--finish;
    	}
    }
}

void modify_tdim(fitsfile *cfp, int colnum, int naxes, long *axes, int *status)
{
    int i;
    char astr[48], bstr[16];

    strcpy( astr, "(" );
    for( i=0; i<naxes; i++ ) {
    	sprintf( bstr, "%ld,", axes[i] );
    	strcat( astr, bstr );
    }
    /* Strip off final ',' */
    astr[ strlen( astr ) - 1 ] = '\0';
    strcat( astr, ")" );

    /* Construct key */
    sprintf( bstr, "TDIM%d", colnum );

    fits_update_key( cfp, TSTRING, bstr, astr, NULL, status );
}

float get_random_number()
{
    return drand48();
}

