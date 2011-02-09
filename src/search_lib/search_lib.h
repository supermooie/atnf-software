#include "fitsio.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

#define MAX_CHANS 1024
//#define MAX_CHANS 512
#define MAX_BINS 1024
#define MAX_INT16 16383.0
#define MAX_SUBINTS 1024
#define MAX_POLNS 4
#define MAX_CHARS 180

void eightBits(int eight_bit_number, float results[], int *index);
void fourBits(int eight_bit_number, float results[], int *index);
void twoBits(int eight_bit_number, float results[], int *index);
void fourBits(int eight_bit_number, float results[], int *index);
void oneBit(int eight_bit_number, float results[], int *index);

float fourBitNumber(int num);
float eightBitNumber(int num);
float twoBitNumber(int num);
float oneBitNumber(int num);

float raw_oneBitNumber(int num);
float raw_eightBitNumber(int num);
float raw_fourBitNumber(int num);
float raw_twoBitNumber(int num);

void reverse_channels(unsigned char values[], unsigned char reversed_values[],
        int nbits, int nchan);

void minmax(float profile[], int nbin, float *min, float *max);

struct arch_params {
    fitsfile *fp;
    char stt_date[MAX_CHARS];
    char stt_time[MAX_CHARS];
    int imjd;
    int smjd;
    int nsamp;
    int nchan;
    int nbits;
    float tbin;
    float equinox;
    float cal_freq;
    float cal_dcyc;
    float cal_phs;
    int npol;
    int nsub;
    int nrcvr;
    int be_phase;
    int be_dcc;
    int be_delay;
    int tcycle;
    int obsnchan;
    int scanlen;
    float fa_req;
    double offs;
    double lst;
    double ant_x;
    double ant_y;
    double ant_z;
    double bmaj;
    double bmin;
    double bpa;
    char date[MAX_CHARS];
    char cal_mode[MAX_CHARS];
    char fd_mode[MAX_CHARS];
    char trk_mode[MAX_CHARS];
    char stt_crd1[MAX_CHARS];
    char stt_crd2[MAX_CHARS];
    char stp_crd1[MAX_CHARS];
    char stp_crd2[MAX_CHARS];
    char ra[MAX_CHARS];
    char dec[MAX_CHARS];
    char src_name[MAX_CHARS];
    char obs_mode[MAX_CHARS];
    char beconfig[MAX_CHARS];
    char date_obs[MAX_CHARS];
    char observer[MAX_CHARS];
    char proj_id[MAX_CHARS];
    char telescope[MAX_CHARS];
    char frontend[MAX_CHARS];
    char backend[MAX_CHARS];
    char fd_poln[MAX_CHARS];
    char fd_hand[MAX_CHARS];
    char fd_sang[MAX_CHARS];
    char fd_xyph[MAX_CHARS];
    char coord_md[MAX_CHARS];
    char scale[MAX_CHARS];
    double chan_bw;
    float freq0;
    float obsfreq;
    float obsbw;
    double lst_subs[MAX_SUBINTS];
};

void get_params(struct arch_params *ap);
void modify_tdim(fitsfile *cfp, int colnum, int naxes, long *axes, int *status);
fitsfile *create_archive(struct arch_params ap, int nbin, char *eph_file,
        char *filename);

void make_subint(fitsfile* afp, struct arch_params ap, int subint_cnt, int nbin,
        //float profile[MAX_CHANS][MAX_BINS], double offs_sub, double lst_sub);
        float profile[MAX_POLNS][MAX_CHANS][MAX_BINS], double offs_sub, double lst_sub);

void readHeader(fitsfile *fp, int *nsblk, int *nchan, int *nbits, float *tsamp,
        int *npol, int *nsubint)
{
    int status = 0;
    fits_movnam_hdu(fp, BINARY_TBL, "SUBINT", 0, &status);
    if (status) {
        printf("Unable to move to subint table in FITS file\n");
        exit(1);
    }

    long numrows;
    char chout[100];
    char comment[1000];
    fits_get_num_rows (fp, &numrows, &status);
    printf("Number of subints (NSUB)                = %d\n",(int)numrows);
    *nsubint = (int)numrows;

    fits_read_keyword(fp, "NSBLK", chout, comment, &status); 
    sscanf(chout,"%d",nsblk);
    printf("Number of samples per subint (NSBLK)    = %d\n", *nsblk);
    fits_read_keyword(fp, "NCHAN", chout, comment, &status); 
    sscanf(chout,"%d",nchan);
    printf("Number of channels per sample (NCHAN)   = %d\n", *nchan);
    fits_read_keyword(fp, "NBITS", chout, comment, &status); 
    sscanf(chout, "%d", nbits);
    printf("Number of bits per datum                = %d\n", *nbits);
    fits_read_keyword(fp, "NPOL", chout, comment, &status); 
    sscanf(chout, "%d", npol);
    printf("Number of poln (NPOL)                   = %d\n", *npol);
    fits_read_keyword(fp, "TBIN", chout, comment, &status); 
    sscanf(chout,"%f",tsamp);
    printf("Sample time (TBIN)                      = %g sec\n",*tsamp);
    printf("Length of each subint (TBIN*NSBLK)      = %g sec\n",*tsamp**nsblk);
    printf("Length of observation (TBIN*NSBLK*NSUB) = %g sec\n",
            *tsamp**nsblk*numrows);

    printf("\n");
}

void displayData(fitsfile *fp, int subint, int pol, int colnum, int nchan,
        int samplesperbyte, float results[])
{
    int status = 0;
    int initflag = 0;
    unsigned char nval='*',val2;
    int s = 1, p = 1, i;
    int index = 0;

    fits_get_colnum (fp, CASEINSEN, "DATA", &colnum, &status);
    for (i = 1; i <= 10; i++) {
        fits_read_col_byt(fp, colnum, subint, (s-1)*2*nchan/samplesperbyte + 
                (p-1)* nchan/samplesperbyte + i, 1, nval, &val2,  &initflag, 
                &status);

        fourBits(val2, results, &index);
    }
}

float getChan(fitsfile *fp, int sample,int subint, int npol, int pol,
        int samplesperbyte, int nchan, int chan)
{
    int status = 0;
    int initflag = 0;
    unsigned char val2, nval = '0';

    int colnum;
    fits_get_colnum (fp, CASEINSEN, "DATA", &colnum, &status);
    fits_read_col_byt(fp, colnum, subint, 1, 1, nval, &val2,  &initflag, 
            &status);

    int mod = (samplesperbyte - 1) - (chan % samplesperbyte);
    int shiftedNumber = val2 >> (mod * 8/samplesperbyte);

    switch (samplesperbyte) {
        case 1:
            return eightBitNumber(shiftedNumber);

        case 2:
            return fourBitNumber(shiftedNumber);

        case 4:
            return twoBitNumber(shiftedNumber);

        case 8:
            return oneBitNumber(shiftedNumber);
        default:
            return 0.0;
    }
}

void bytesToValues(int samplesperbyte, int nchan, unsigned char val2[], 
        float values[])
{
    int i;
    int index = 0;
    for (i = 0; i < nchan/samplesperbyte; i++) {
        void (*pointer)(int, float[], int*);
        switch (samplesperbyte)
        {
            case 1:
                pointer = eightBits;
                break;

            case 2:
                pointer = fourBits;
                break;

            case 4:
                pointer = twoBits;
                break;

            case 8:
                pointer = oneBit;
                break;
        }
        pointer(val2[i], values, &index);
    }
}

void getChans(fitsfile *fp, struct arch_params ap, int sample, int subint, 
        int pol, float values[], int colnum)
{
    int samplesperbyte = 8/ap.nbits;
    int status = 0;
    int initflag = 0;
    unsigned char val2[ap.nchan/samplesperbyte];
    unsigned char nval = '0';

    fits_read_col_byt(fp, colnum, subint, sample*ap.npol*ap.nchan/samplesperbyte +
            pol*ap.nchan/samplesperbyte + 1, ap.nchan/samplesperbyte, nval, val2,
            &initflag, &status);

    bytesToValues(samplesperbyte, ap.nchan, val2, values);
    //values[0] = 10;
}

void printChan(fitsfile *fp, int subint, int pol, int samplesperbyte, int nchan,
        int chan_begin, int chan_end, int nsamp, int print, float values[],
        float *rms, float *mean, int scaled)
{
    int counter = 0;

    int status = 0;
    int initflag = 0;
    unsigned char nval = '*', val2;
    int sample;
    int npol = 1;
    int chan;

    int colnum;
    fits_get_colnum (fp, CASEINSEN, "DATA", &colnum, &status);

    for (sample= 1; sample <= nsamp; sample++) {
        for (chan = chan_begin; chan <= chan_end; chan++) {

            fits_read_col_byt(fp, colnum, subint, (sample-1)*npol*nchan/
                    samplesperbyte + (pol-1)*nchan/samplesperbyte + chan/
                    samplesperbyte, 1, nval, &val2, &initflag, &status);

            int mod = chan % samplesperbyte;
            int shiftedNumber = val2 >> (mod * 8/samplesperbyte);

            switch (samplesperbyte) {
                case 1: 
                {
                    float result;
                    if (scaled) {
                        result = eightBitNumber(shiftedNumber);
                        values[(int)(result + 31.5)]++;
                    } else {
                        result = raw_eightBitNumber(shiftedNumber);
                        values[(int)(result + 7.0)]++;
                    }
                    *mean += result;
                    *rms += pow(result, 2);
                } break;
                case 2:
                {
                    float result;
                    if (scaled) {
                        result = fourBitNumber(shiftedNumber);
                        values[(int)(result + 7.5)]++;
                    } else {
                        result = raw_fourBitNumber(shiftedNumber);
                        values[(int)(result + 8.0)]++;
                    }
                    *mean += result;
                    *rms += pow(result, 2);

                } break;
                case 4: {
                    float result;
                    if (scaled) {
                        result = twoBitNumber(shiftedNumber);
                        if (result == -2.0)
                            values[0]++;
                        else if (result == -0.5)
                            values[1]++;
                        else if (result == 0.5)
                            values[2]++;
                        else if (result == 2.0)
                            values[3]++;
                    } else {
                        result = raw_twoBitNumber(shiftedNumber);
                        values[(int)(result + 2.0)]++;
                    }

                    *mean += result;
                    *rms += pow(result, 2);
                } break;
                case 8: {
                    float result;
                    if (scaled) {
                        result = oneBitNumber(shiftedNumber);
                        values[(int)(result + 0.5)]++;
                    } else {
                        result = raw_oneBitNumber(shiftedNumber);
                        values[(int)(result + 1.0)]++;
                    }
                    *mean += result;
                    *rms += pow(result, 2);
                } break;
            }
            counter++;
        }
    }

    *mean = *mean / (float)(counter);
    *rms = *rms / (float)(counter);
    *rms = sqrt(*rms);
}

// raw numbers
float raw_oneBitNumber(int num)
{
    float result = num & 1;
    if (result)
        return -1.0;
    else
        return 0.0;
}

float raw_eightBitNumber(int num)
{
    return num - 136.0;
}

float raw_fourBitNumber(int num)
{
    float result = num & 15;
    if (result >= 8)
        return result - 16;
    else
        return result;
}

float raw_twoBitNumber(int num)
{
    return num & 3;
}


// scaled
float oneBitNumber(int num)
{
    return (num & 1) - 0.5;
}

float eightBitNumber(int num)
{
    return num - 31.5;
}

float fourBitNumber(int num)
{
    float result = num & 15;
    return result - 7.5;
}

float twoBitNumber(int num)
{
    num = num & 3;
    switch (num) {
        case 0:
            return -2.0;

        case 1:
            return -0.5;

        case 2:
            return 0.5;

        case 3:
            return 2.0;

        default:
            return 0.0;
    }
}

void getSample(fitsfile *fp, int npol, int subint, int pol, int sample,
        int samplesperbyte, int nchan, float results[])
{
    int status = 0;
    int initflag = 0;
    unsigned char nval = '*', val2;
    int index = 0;
    int c;
    int colnum;

    fits_get_colnum (fp, CASEINSEN, "DATA", &colnum, &status);

    for (c = 1; c <= nchan/samplesperbyte; c++) {
        fits_read_col_byt(fp, colnum, subint, (sample-1)*npol*nchan/
                samplesperbyte + (pol-1)*nchan/samplesperbyte + c, 1, nval, 
                &val2,  &initflag, &status);

        switch (samplesperbyte) {
            case 1:
                eightBits(val2, results, &index);
                break;
            case 2:
                fourBits(val2, results, &index);
                break;
            case 4:
                twoBits(val2, results, &index);
                break;
            case 8:
                oneBit(val2, results, &index);
                break;
        }
    }
}
void eightBits(int eight_bit_number, float results[], int *index)
{
    // 135 subtracted from the original number to make the range from -8 to something
    // 0.5 is then added to remove bias

    //results[(*index)++] = eight_bit_number - 135.5;
    results[(*index)++] = eight_bit_number -31.5;
}

void fourBits(int eight_bit_number, float results[], int *index)
{
    // anding the least significant 4 bits with 0000 1111 will give the (signed) 4-bit number
    // shifting right 4 bits will produce the other (first) 4-bit numbers
    // 0.5 is added to each number to compensate for the bias, making the range -7.5 -> +7.5

    float tempResults[2];
    int i;
    for (i = 0; i < 2; i++) {
        int andedNumber = eight_bit_number & 15;
        tempResults[i] = andedNumber - 7.5;
        eight_bit_number = eight_bit_number >> 4;
    }

    for (i = 1; i >= 0; i--)
        results[(*index)++] = tempResults[i];
}

void twoBits(int eight_bit_number, float results[], int *index)
{
    // anding the least significant 2 bits with 0000 0011 will give the (signed 2-bit number
    // shifting right 2 bits will produce the next (previous) 2-bit number
    // the numbers are adjusted to compensate for bias and to give a rms of ~0.9
    //
    // 1 -> 2.0
    // 0 -> +0.5
    // -1 -> -0.5
    // -2 -> -2.0

    float tempResults[4];
    int i;

    for (i = 0; i < 4; i++) {
        int andedNumber = eight_bit_number & 3;
        switch (andedNumber) {
            case 0:
                tempResults[i] = -2.5;
                break;
            case 1:
                tempResults[i] = -0.5;
                break;
            case 2:
                tempResults[i] = 0.5;
                break;
            case 3:
                tempResults[i] = 2.5;
                break;
        }
        eight_bit_number = eight_bit_number >> 2;
    }

    for (i = 3; i >= 0; i--)
        results[(*index)++] = tempResults[i];
}

void oneBit(int eight_bit_number, float results[], int *index)
{
    // anding the least significant bit with 0000 0001 will give the (signed 1-bit number
    // shifting right 1 bit will produce the next (previous) 1-bit number
    //
    // 0.5 is added to each number to compensate for bias
    //

    float tempResults[8];
    int i;

    for (i = 0; i < 8; i++) {
        int andedNumber = eight_bit_number & 1;
        tempResults[i] = andedNumber ? 0.5 : -0.5;
        eight_bit_number = eight_bit_number >> 1;
    }

    for (i = 7; i >= 0; i--)
        results[(*index)++] = tempResults[i];
}

fitsfile *create_archive(struct arch_params ap, int nbin, char *eph_file,
        char *filename)
{
#define GSTR_LEN 160
#define MAX_BLKS 20
    char gstr[GSTR_LEN];
    char *pch[MAX_BLKS];
    char fname[100];
    int status = 0;
    fitsfile *afp;
    int j, k;
    double dx;
    int colnum;
    FILE *fp;

    //sprintf(fname, "!%s.ar(psrheader_gh.fits)", strtok(filename, "."));
    sprintf(fname, "!%s.ar(/pulsar/psr/cvshome/kho018/psrchive/Base/Formats/PSRFITS/psrheader.fits)", strtok(filename, "."));

    fits_create_file(&afp, fname, &status);
    if (fits_movnam_hdu(afp, BINARY_TBL, "COHDDISP", 0, &status) == 0)
        fits_delete_hdu(afp, &k, &status);
    else status = 0;  /* In case the table was not found */

    if (fits_movnam_hdu(afp, BINARY_TBL, "FLUX_CAL", 0, &status) == 0)
        fits_delete_hdu(afp, &k, &status);
    else status = 0;  /* In case the table was not found */

    if (fits_movnam_hdu(afp, BINARY_TBL, "CAL_POLN", 0, &status) == 0)
        fits_delete_hdu(afp, &k, &status);
    else status = 0;  /* In case the table was not found */

    if (fits_movnam_hdu(afp, BINARY_TBL, "FEEDPAR", 0, &status) == 0)
        fits_delete_hdu(afp, &k, &status);
    else status = 0;  /* In case the table was not found */

    //float obsfreq = ap.freq0 + ap.nchan/2*ap.chan_bw;

    fits_movabs_hdu(afp, 1, NULL, &status);
    fits_update_key(afp, TSTRING, "OBS_MODE", "PSR", NULL, &status);
    fits_update_key(afp, TSTRING, "STT_DATE", ap.stt_date, NULL, &status);
    fits_update_key(afp, TSTRING, "STT_TIME", ap.stt_time, NULL, &status);
    fits_update_key(afp, TSTRING, "SRC_NAME", ap.src_name, NULL, &status);
    fits_update_key(afp, TSTRING, "BECONFIG", ap.beconfig, NULL, &status);
    fits_update_key(afp, TSTRING, "RA", ap.ra, NULL, &status);
    fits_update_key(afp, TSTRING, "DEC", ap.dec, NULL, &status);
    fits_update_key(afp, TSTRING, "STT_CRD1", ap.stt_crd1, NULL, &status);
    fits_update_key(afp, TSTRING, "STT_CRD2", ap.stt_crd2, NULL, &status);
    fits_update_key(afp, TSTRING, "STP_CRD1", ap.stp_crd1, NULL, &status);
    fits_update_key(afp, TSTRING, "STP_CRD2", ap.stp_crd2, NULL, &status);
    fits_update_key(afp, TSTRING, "TRK_MODE", ap.trk_mode, NULL, &status);
    fits_update_key(afp, TSTRING, "FD_MODE", ap.fd_mode, NULL, &status);
    fits_update_key(afp, TSTRING, "CAL_MODE", ap.cal_mode, NULL, &status);
    fits_update_key(afp, TFLOAT, "FA_REQ", &(ap.fa_req), NULL, &status);
    fits_update_key(afp, TFLOAT, "CAL_FREQ", &(ap.cal_freq), NULL, &status);
    fits_update_key(afp, TFLOAT, "CAL_DCYC", &(ap.cal_dcyc), NULL, &status);
    fits_update_key(afp, TFLOAT, "CAL_PHS", &(ap.cal_phs), NULL, &status);
    fits_update_key(afp, TINT, "SCANLEN", &(ap.scanlen), NULL, &status);
    fits_update_key(afp, TINT, "STT_IMJD", &(ap.imjd), NULL, &status);
    fits_update_key(afp, TINT, "STT_SMJD", &(ap.smjd), NULL, &status);
    fits_update_key(afp, TINT, "NRCVR", &(ap.nrcvr), NULL, &status);
    fits_update_key(afp, TINT, "BE_PHASE", &(ap.be_phase), NULL, &status);
    fits_update_key(afp, TINT, "BE_DCC", &(ap.be_dcc), NULL, &status);
    fits_update_key(afp, TINT, "BE_DELAY", &(ap.be_delay), NULL, &status);
    fits_update_key(afp, TINT, "TCYCLE", &(ap.tcycle), NULL, &status);
    fits_update_key(afp, TINT, "OBSNCHAN", &(ap.nchan), NULL, &status);
    fits_update_key(afp, TDOUBLE, "STT_OFFS", &(ap.offs), NULL, &status);
    fits_update_key(afp, TDOUBLE, "STT_LST", &(ap.lst), NULL, &status);
    fits_update_key(afp, TDOUBLE, "ANT_X", &(ap.ant_x), NULL, &status);
    fits_update_key(afp, TDOUBLE, "ANT_Y", &(ap.ant_y), NULL, &status);
    fits_update_key(afp, TDOUBLE, "ANT_Z", &(ap.ant_z), NULL, &status);
    fits_update_key(afp, TDOUBLE, "BMAJ", &(ap.bmaj), NULL, &status);
    fits_update_key(afp, TDOUBLE, "BMIN", &(ap.bmin), NULL, &status);
    fits_update_key(afp, TDOUBLE, "BPA", &(ap.bpa), NULL, &status);
    fits_update_key(afp, TSTRING, "DATE", ap.date, NULL, &status);
    //fits_update_key(afp, TSTRING, "OBS_MODE", ap.obs_mode, NULL, &status);
    fits_update_key(afp, TSTRING, "DATE-OBS", ap.date_obs, NULL, &status);
    fits_update_key(afp, TSTRING, "OBSERVER", ap.observer, NULL, &status);
    fits_update_key(afp, TSTRING, "PROJID", ap.proj_id, NULL, &status);
    fits_update_key(afp, TSTRING, "TELESCOP", ap.telescope, NULL, &status);
    fits_update_key(afp, TSTRING, "FRONTEND", ap.frontend, NULL, &status);
    fits_update_key(afp, TSTRING, "BACKEND", ap.backend, NULL, &status);
    fits_update_key(afp, TSTRING, "FD_POLN", ap.fd_poln, NULL, &status);
    fits_update_key(afp, TSTRING, "FD_HAND", ap.fd_hand, NULL, &status);
    fits_update_key(afp, TSTRING, "FD_SANG", ap.fd_sang, NULL, &status);
    fits_update_key(afp, TSTRING, "FD_XYPH", ap.fd_xyph, NULL, &status);
    fits_update_key(afp, TSTRING, "COORD_MD", ap.coord_md, NULL, &status);
    fits_update_key(afp, TSTRING, "SCALE", ap.scale, NULL, &status);
    fits_update_key(afp, TFLOAT, "OBSFREQ", &(ap.obsfreq), NULL, &status);
    fits_update_key(afp, TFLOAT, "EQUINOX", &(ap.equinox), NULL, &status);
    fits_update_key(afp, TFLOAT, "OBSBW", &(ap.obsbw), NULL, &status);
    fits_movnam_hdu(afp, BINARY_TBL, "SUBINT", 0, &status);
    fits_update_key(afp, TINT, "NPOL", &(ap.npol), NULL, &status);
    fits_update_key(afp, TINT, "NBIN", &nbin, NULL, &status);
    fits_update_key(afp, TDOUBLE, "TBIN", &(ap.tbin), NULL, &status);
    fits_update_key(afp, TINT, "NCHAN", &(ap.nchan), NULL, &status);
    fits_update_key(afp, TDOUBLE, "CHAN_BW", &(ap.chan_bw), NULL, &status);

    if (fits_movnam_hdu(afp, BINARY_TBL, "POLYCO", 0, &status) == 0)
        fits_delete_hdu(afp, &k, &status);
    else status=0;
    fits_movnam_hdu(afp, BINARY_TBL, "T2PREDICT", 0, &status);

    if (!fits_get_colnum(afp, CASEINSEN, "PREDICT", &colnum, &status)) {
        k = 1;
        if (!(fp = fopen("t2pred.dat","r"))) {
            printf("Error opening t2pred.dat\n");
            exit(1);
        }

        while (fgets(gstr, GSTR_LEN, fp) != NULL) {
            j = strlen(gstr);
            if (gstr[j-1] == '\n') gstr[j-1] = '\0';

            pch[0] = gstr;
            fits_write_col(afp, TSTRING, colnum, k, 1, 1, pch, &status);
            k++;
        }
    }
    fits_report_error(stderr,status);
    fflush(stderr);

    /* Write to Pulsar Ephemeris BINTABLE */
    /* Move to required HDU */
    int add_fract_part, f0_column;
    char msg[100];
    char *qch, *rch, *sch, *tch;
    char astr[80];
    int coltype;
    long repeat, width;
    long lk;
    void *pv;
    short int sj;

    fits_movnam_hdu(afp, BINARY_TBL, "PSREPHEM", 0, &status);
    if ((fp = fopen(eph_file, "r")) == NULL) {
        printf("ERROR opening file '%s'\n", eph_file);
        exit(1);
    }

    while( fgets( gstr, GSTR_LEN, fp ) != NULL ) {

        if( ( qch = strtok( gstr, " \t" ) ) != NULL ) {

          add_fract_part = 0;
            f0_column = 0;

            /* Following gets a name change - From PSRJ to PSR_NAME */
            if( strcmp( qch, "PSRJ" ) == 0 ) strcpy( astr, "PSR_NAME" );
            /* Following two get split into integer and fractional parts */
            else if( ( strcmp( qch, "F" ) == 0 ) || ( strcmp( qch, "F0" ) == 0 ) ) {
                strcpy( astr, "IF0" );
                add_fract_part = 1;
                f0_column = 1;
            }
            else if( strcmp( qch, "TZRMJD" ) == 0 ) {
                strcpy( astr, "TZRIMJD" );
                add_fract_part = 1;
            }
            else {
                strcpy( astr, qch );
            }

            if( fits_get_colnum( afp, CASEINSEN, astr, &colnum, &status ) == 0 ) {

                fits_get_coltype( afp, colnum, &coltype, &repeat, &width, &status );

                /* Get next parameter */
                if( ( rch = strtok( NULL, " \t" ) ) == NULL ) {
                    //sprintf( msg, "ERROR - No value for PSREPHEM entry '%s'", astr );
                    printf("ERROR - No value for PSREPHEM entry '%s'", astr );
                    exit(1);
                }

                if( coltype == TDOUBLE || add_fract_part ) {
                    /* Find exponent and change any D,d,e exponent to E */
                    if( ( sch = strrchr( rch, 'D' ) ) != NULL ) *sch = 'E';
                    else if( ( sch = strrchr( rch, 'd' ) ) != NULL ) *sch = 'E';
                    else if( ( sch = strrchr( rch, 'e' ) ) != NULL ) *sch = 'E';
                    else sch = strrchr( rch, 'E' );

                    tch = sch;

                    if( add_fract_part ) {
                        /* Convert to integer + fractional part AND
                           if F0 convert from Hz to mHz i.e. multiply by 1000 */
                        /* Get exponent */
                        if( sch != NULL ) {
                            if( sscanf( ++sch, "%d", &j ) != 1 ) {
                                printf("exit\n");
                                exit(1);
                            }
                        }
                        else j = 0;

                        /* get position of decimal point */
                        if( ( sch = strrchr( rch, '.' ) ) == NULL ) {
                            /* Integer part only */
                            if( sscanf( rch, "%ld", &lk ) != 1 )  {
                                printf("exit\n");
                                exit(1);
                            }
                        }
                        else j = 0;

                        /* get position of decimal point */
                        if( ( sch = strrchr( rch, '.' ) ) == NULL ) {
                            /* Integer part only */
                            if( sscanf( rch, "%ld", &lk ) != 1 ) {
                                printf("exit\n");
                                exit(1);
                            }
                            dx = 0.0;
                        }
                        else {
                            /* Decimal point present */
                            if( f0_column ) j += 3;  /* multiply by 1000 - sec. to ms. */

                            if( j >= 0 ) {
                                k = (int) ( sch - rch );
                                if( k > 0 ) {
                                    strncpy( astr, rch, k );
                                    astr[k] = '\0';
                                }
                                else strcpy( astr, "0" );

                                if( j > 0 ) { /* multiply by 10**j, i.e. move decimal point */
                                    strncat( astr, (sch+1), j );
                                    *(sch += j) = '.';
                                    *--sch = '0';
                                }
                                /* Now astr contains integer part, sch points to fractional part */
                                if( sscanf( astr, "%ld", &lk ) != 1 )  {
                                    printf("exit\n");
                                    exit(1);
                                }
                                /* Don't want exponent - reduced to fraction */
                                if( tch != NULL ) *tch = '\0';
                                if( sscanf( sch, "%lf", &dx ) != 1 ) {
                                    printf("exit\n");
                                    exit(1);
                                }
                            }
                            else {
                                /* Fractional part only */
                                lk = 0;
                                if( sscanf( sch, "%lf", &dx ) != 1 ) {
                                    printf("exit\n");
                                    exit(1);
                                }
                                if( f0_column ) dx *= 1000.0;
                            }

                        }

                        pv = &lk;
                        j = TLONG;

                    }
                    else {

                        if( sscanf( rch, "%lf", &dx ) != 1 ) {
                            printf("exit\n");
                            exit(1);
                        }

                        j = TDOUBLE;
                        pv = &dx;
                    }
                }
                else if( coltype == TSTRING ) {
                    j = TSTRING;
                    /* Strip off any trailing newline */
                    if( ( sch = strrchr( rch, '\n' ) ) != NULL ) *sch = '\0';
                    pch[0] = rch;
                    pv = pch;
                }
                else if( coltype == TSHORT ) {
                    j = TSHORT;
                    if( sscanf( rch, "%d", &k ) != 1 ) {
                        printf("exit\n");
                        exit(1);
                    }
                    sj = k;
                    pv = &sj;
                }
                else if( coltype == TINT ) {
                    j = TINT;
                    if( sscanf( rch, "%d", &k ) != 1 ) {
                        printf("exit\n");
                        exit(1);
                    }
                    pv = &k;
                }
                else {
                    printf("exit\n");
                    exit(1);
                }

                fits_write_col( afp, j, colnum, 1, 1, 1, pv, &status );

                /* Write fractional part */
                if( add_fract_part ) {
                    j = TDOUBLE;
                    pv = &dx;
                    fits_write_col( afp, j, colnum+1, 1, 1, 1, pv, &status );

                    /* printf( "\nCFITS: In Ephemeris BINTABLE: %s = %ld %e", qch, lk, dx ); */
                }

            }
            else {
                //sprintf( msg, "ERROR - PSREPHEM entry '%s' NOT in table", astr );
                printf("ERROR - PSREPHEM entry '%s' NOT in table", astr );
                exit(1);
            }
        }
    }

    fits_report_error(stderr,status);
    fclose(fp);
    return afp;
}

void make_subint(fitsfile* afp, struct arch_params ap, int subint_cnt, int nbin,
        float profile[MAX_POLNS][MAX_CHANS][MAX_BINS], double offs_sub, double lst_sub)
{
#define GSTR_LEN 160
#define MAX_BLKS 20
    int status = 0;
    long long lj;
    int colnum;
    int nscl = ap.nchan*ap.npol;
    int ndata = nbin*nscl;
    int j,k;
    long naxes[3];
    double dx;
    int i;

    fits_movnam_hdu(afp, BINARY_TBL, "SUBINT", 0, &status);
    fits_report_error(stderr, status);

    lj = ap.nchan;  // Number of channels
    if(fits_get_colnum(afp, CASEINSEN, "DAT_FREQ", &colnum, &status) == 0)
        fits_modify_vector_len(afp, colnum, lj, &status);
    if(fits_get_colnum(afp, CASEINSEN, "DAT_WTS", &colnum, &status) == 0)
        fits_modify_vector_len(afp, colnum, lj, &status);

    lj = nscl;
    if(fits_get_colnum(afp, CASEINSEN, "DAT_OFFS", &colnum, &status) == 0)
        fits_modify_vector_len(afp, colnum, lj, &status);
    if(fits_get_colnum(afp, CASEINSEN, "DAT_SCL", &colnum, &status) == 0)
        fits_modify_vector_len(afp, colnum, lj, &status);
    fits_update_key(afp, TINT, "NBIN", &nbin, NULL, &status);
    fits_update_key(afp, TINT, "NCHAN", &(ap.nchan), NULL, &status);

    lj = ndata;
    if(fits_get_colnum(afp, CASEINSEN, "DATA", &colnum, &status) == 0) {
        fits_modify_vector_len(afp, colnum, lj, &status);
        naxes[0] = nbin;
        naxes[1] = ap.nchan;
        naxes[2] = ap.npol;

        modify_tdim(afp, colnum, 3, naxes, &status);
    }

    fits_report_error(stderr, status);

    /* Add keywords */
    fits_update_key(afp, TSTRING, "INT_TYPE", "TIME", NULL, &status);
    fits_update_key(afp, TSTRING, "INT_UNIT", "SEC", NULL, &status);
    fits_update_key(afp, TINT, "NCH_FILE", &(ap.nchan), NULL, &status);
    j = 0;
    fits_update_key(afp, TINT, "NCH_STRT", &j, NULL, &status);
    fits_update_key(afp, TINT, "NPOL", &(ap.npol), NULL, &status);
    fits_update_key(afp, TINT, "NSBLK", &(ap.nsamp), NULL, &status);

    /* [s] Length of subint */
    if(fits_get_colnum(afp, CASEINSEN, "TSUBINT", &colnum, &status) == 0) {
        dx = (double)ap.tbin * (double)ap.nsamp;
        fits_write_col(afp, TDOUBLE, colnum, subint_cnt, 1, 1, &dx, &status);
    }

    fits_get_colnum(afp, CASEINSEN, "OFFS_SUB", &colnum, &status);
    fits_write_col(afp, TDOUBLE, colnum, subint_cnt, 1, 1, &offs_sub, &status);

    fits_get_colnum(afp, CASEINSEN, "LST_SUB", &colnum, &status);
    fits_write_col(afp, TDOUBLE, colnum, subint_cnt, 1, 1, &lst_sub, &status);
    // SHOULD SET HEAPS OF PARAMETERS LIKE TEL_ZEN ETC.

    /* Centre freq. for each channel - NCHAN floats */
    if(fits_get_colnum(afp, CASEINSEN, "DAT_FREQ", &colnum, &status) == 0)
    {
        float binned_freq[ap.nchan];
        float binned_weight[ap.nchan];

        for (i=0;i<ap.nchan;i++) {
            binned_freq[i] = ap.freq0+ap.chan_bw*i;
            binned_weight[i] = 1.0;
        }

        fits_write_col(afp, TFLOAT, colnum, subint_cnt, 1, ap.nchan, binned_freq, &status);
        if(fits_get_colnum(afp, CASEINSEN, "DAT_WTS", &colnum, &status) == 0)
            fits_write_col(afp, TFLOAT, colnum, subint_cnt, 1, ap.nchan, binned_weight, &status);

    }

    float max_short = pow(2.0,15.0) - 1.0;
    float dat_offss[ap.npol][ap.nchan];
    float dat_scls[ap.npol][ap.nchan];

    short int binned_data[ap.npol][ap.nchan][nbin];
    for (i = 0; i < ap.npol; ++i) {
        for (j = 0; j < ap.nchan; ++j) {
            float min;
            float max;

            minmax(profile[i][j], nbin, &min, &max);

            const float offset = 0.5 * (max + min);
            const float scale_factor = (max - min) / max_short;
            const float dat_offs = (max - min)/2.0;

            dat_offss[i][j] = (max - min)/2.0;
            dat_scls[i][j] = dat_offs / ((max_short/2.0) - 1);

            for (k = 0; k < nbin; ++k) {
                binned_data[i][j][k] = (short int)((profile[i][j][k] - offset) /
                        scale_factor);
            }
        }
    }

    if(fits_get_colnum(afp, CASEINSEN, "DATA", &colnum, &status) == 0)
        fits_write_col(afp, TSHORT, colnum, subint_cnt, 1, ndata, binned_data, &status);

    if(fits_get_colnum(afp, CASEINSEN, "DAT_OFFS", &colnum, &status) == 0)
        fits_write_col(afp, TFLOAT, colnum, subint_cnt, 1, ap.nchan * ap.npol, dat_offss, &status);

    if(fits_get_colnum(afp, CASEINSEN, "DAT_SCL", &colnum, &status) == 0)
        fits_write_col(afp, TFLOAT, colnum, subint_cnt, 1, ap.nchan * ap.npol, dat_scls, &status);

    fits_report_error(stderr, status);

    /* Now FLUSH any internal buffers to the file */
    fits_flush_file(afp, &status);
}

void modify_tdim(fitsfile *cfp, int colnum, int naxes, long *axes, int *status)
{
    int i;
    char astr[48], bstr[16];

    strcpy(astr, "(");
    for(i=0; i<naxes; i++) {
        sprintf(bstr, "%ld,", axes[i]);
        strcat(astr, bstr);
    }
    /* Strip off final ',' */
    astr[ strlen(astr) - 1 ] = '\0';
    strcat(astr, ")");

    /* Construct key */
    sprintf(bstr, "TDIM%d", colnum);

    fits_update_key(cfp, TSTRING, bstr, astr, NULL, status);
}

void copy_history_params(fitsfile *arch_fp, fitsfile *srch_fp)
{
    int status = 0;
    fits_movnam_hdu(arch_fp, BINARY_TBL, "HISTORY", 0, &status);
    fits_movnam_hdu(srch_fp, BINARY_TBL, "HISTORY", 0, &status);
    if (status) {
        printf("Unable to move to history table in FITS file\n");
        exit(1);
    }

    int arch_colnum;
    int srch_colnum;
    int anul = 0;
    char* nul = strdup(" ");
    char strval[1][MAX_CHARS];

    int ival;

    fits_get_colnum(arch_fp, CASEINSEN, "DATE_PRO", &arch_colnum, &status);
    fits_get_colnum(srch_fp, CASEINSEN, "DATE_PRO", &srch_colnum, &status);
    fits_read_col(srch_fp, TSTRING, srch_colnum, 1, 1, 1, &nul, strval, &anul, &status);
    fits_write_col(arch_fp, TSTRING, arch_colnum, 1, 1, 1, strval, &status);

    fits_get_colnum(arch_fp, CASEINSEN, "PROC_CMD", &arch_colnum, &status);
    fits_get_colnum(srch_fp, CASEINSEN, "PROC_CMD", &srch_colnum, &status);
    fits_read_col(srch_fp, TSTRING, srch_colnum, 1, 1, 1, &nul, strval, &anul, &status);
    fits_write_col(arch_fp, TSTRING, arch_colnum, 1, 1, 1, strval, &status);

    fits_get_colnum(arch_fp, CASEINSEN, "SCALE", &arch_colnum, &status);
    fits_get_colnum(srch_fp, CASEINSEN, "SCALE", &srch_colnum, &status);
    fits_read_col(srch_fp, TSTRING, srch_colnum, 1, 1, 1, &nul, strval, &anul, &status);
    fits_write_col(arch_fp, TSTRING, arch_colnum, 1, 1, 1, strval, &status);

    fits_get_colnum(arch_fp, CASEINSEN, "POL_TYPE", &arch_colnum, &status);
    fits_get_colnum(srch_fp, CASEINSEN, "POL_TYPE", &srch_colnum, &status);
    fits_read_col(srch_fp, TSTRING, srch_colnum, 1, 1, 1, &nul, strval, &anul, &status);
    fits_write_col(arch_fp, TSTRING, arch_colnum, 1, 1, 1, strval, &status);
 
    fits_get_colnum(arch_fp, CASEINSEN, "NSUB", &arch_colnum, &status);
    fits_get_colnum(srch_fp, CASEINSEN, "NSUB", &srch_colnum, &status);
    fits_read_col(srch_fp, TINT, srch_colnum, 1, 1, 1, NULL, &ival, NULL, &status);
    fits_write_col(arch_fp, TINT, arch_colnum, 1, 1, 1, &ival, &status);


    /*fits_get_colnum(arch_fp, CASEINSEN, "NPOL", &arch_colnum, &status);
    fits_get_colnum(srch_fp, CASEINSEN, "NPOL", &srch_colnum, &status);
    fits_read_col(srch_fp, TINT, srch_colnum, 1, 1, 1, NULL, &ival, NULL, &status);
    fits_write_col(arch_fp, TINT, arch_colnum, 1, 1, 1, &ival, &status);*/

    fits_get_colnum(arch_fp, CASEINSEN, "NBIN", &arch_colnum, &status);
    fits_get_colnum(srch_fp, CASEINSEN, "NBIN", &srch_colnum, &status);
    fits_read_col(srch_fp, TINT, srch_colnum, 1, 1, 1, NULL, &ival, NULL, &status);
    fits_write_col(arch_fp, TINT, arch_colnum, 1, 1, 1, &ival, &status);

    fits_get_colnum(arch_fp, CASEINSEN, "NBIN_PRD", &arch_colnum, &status);
    fits_get_colnum(srch_fp, CASEINSEN, "NBIN_PRD", &srch_colnum, &status);
    fits_read_col(srch_fp, TINT, srch_colnum, 1, 1, 1, NULL, &ival, NULL, &status);
    fits_write_col(arch_fp, TINT, arch_colnum, 1, 1, 1, &ival, &status);


    double dval;
    fits_get_colnum(arch_fp, CASEINSEN, "TBIN", &arch_colnum, &status);
    fits_get_colnum(srch_fp, CASEINSEN, "TBIN", &srch_colnum, &status);
    fits_read_col(srch_fp, TDOUBLE, srch_colnum, 1, 1, 1, NULL, &dval, NULL, &status);
    fits_write_col(arch_fp, TDOUBLE, arch_colnum, 1, 1, 1, &dval, &status);

    fits_get_colnum(arch_fp, CASEINSEN, "CTR_FREQ", &arch_colnum, &status);
    fits_get_colnum(srch_fp, CASEINSEN, "CTR_FREQ", &srch_colnum, &status);
    fits_read_col(srch_fp, TDOUBLE, srch_colnum, 1, 1, 1, NULL, &dval, NULL, &status);
    fits_write_col(arch_fp, TDOUBLE, arch_colnum, 1, 1, 1, &dval, &status);

    fits_get_colnum(arch_fp, CASEINSEN, "NCHAN", &arch_colnum, &status);
    fits_get_colnum(srch_fp, CASEINSEN, "NCHAN", &srch_colnum, &status);
    fits_read_col(srch_fp, TINT, srch_colnum, 1, 1, 1, NULL, &ival, NULL, &status);
    fits_write_col(arch_fp, TINT, arch_colnum, 1, 1, 1, &ival, &status);

    fits_get_colnum(arch_fp, CASEINSEN, "PR_CORR", &arch_colnum, &status);
    fits_get_colnum(srch_fp, CASEINSEN, "PR_CORR", &srch_colnum, &status);
    fits_read_col(srch_fp, TINT, srch_colnum, 1, 1, 1, NULL, &ival, NULL, &status);
    fits_write_col(arch_fp, TINT, arch_colnum, 1, 1, 1, &ival, &status);

    fits_get_colnum(arch_fp, CASEINSEN, "FD_CORR", &arch_colnum, &status);
    fits_get_colnum(srch_fp, CASEINSEN, "FD_CORR", &srch_colnum, &status);
    fits_read_col(srch_fp, TINT, srch_colnum, 1, 1, 1, NULL, &ival, NULL, &status);
    fits_write_col(arch_fp, TINT, arch_colnum, 1, 1, 1, &ival, &status);

    fits_get_colnum(arch_fp, CASEINSEN, "BE_CORR", &arch_colnum, &status);
    fits_get_colnum(srch_fp, CASEINSEN, "BE_CORR", &srch_colnum, &status);
    fits_read_col(srch_fp, TINT, srch_colnum, 1, 1, 1, NULL, &ival, NULL, &status);
    fits_write_col(arch_fp, TINT, arch_colnum, 1, 1, 1, &ival, &status);

    fits_get_colnum(arch_fp, CASEINSEN, "RM_CORR", &arch_colnum, &status);
    fits_get_colnum(srch_fp, CASEINSEN, "RM_CORR", &srch_colnum, &status);
    fits_read_col(srch_fp, TINT, srch_colnum, 1, 1, 1, NULL, &ival, NULL, &status);
    fits_write_col(arch_fp, TINT, arch_colnum, 1, 1, 1, &ival, &status);

    fits_get_colnum(arch_fp, CASEINSEN, "DEDISP", &arch_colnum, &status);
    fits_get_colnum(srch_fp, CASEINSEN, "DEDISP", &srch_colnum, &status);
    fits_read_col(srch_fp, TINT, srch_colnum, 1, 1, 1, NULL, &ival, NULL, &status);
    fits_write_col(arch_fp, TINT, arch_colnum, 1, 1, 1, &ival, &status);

    fits_get_colnum(arch_fp, CASEINSEN, "DDS_MTHD", &arch_colnum, &status);
    fits_get_colnum(srch_fp, CASEINSEN, "DDS_MTHD", &srch_colnum, &status);
    fits_read_col(srch_fp, TSTRING, srch_colnum, 1, 1, 1, &nul, strval, &anul, &status);
    fits_write_col(arch_fp, TSTRING, arch_colnum, 1, 1, 1, strval, &status);

    fits_get_colnum(arch_fp, CASEINSEN, "SC_MTHD", &arch_colnum, &status);
    fits_get_colnum(srch_fp, CASEINSEN, "SC_MTHD", &srch_colnum, &status);
    fits_read_col(srch_fp, TSTRING, srch_colnum, 1, 1, 1, &nul, strval, &anul, &status);
    fits_write_col(arch_fp, TSTRING, arch_colnum, 1, 1, 1, strval, &status);

    fits_get_colnum(arch_fp, CASEINSEN, "CAL_MTHD", &arch_colnum, &status);
    fits_get_colnum(srch_fp, CASEINSEN, "CAL_MTHD", &srch_colnum, &status);
    fits_read_col(srch_fp, TSTRING, srch_colnum, 1, 1, 1, &nul, strval, &anul, &status);
    fits_write_col(arch_fp, TSTRING, arch_colnum, 1, 1, 1, strval, &status);

    fits_get_colnum(arch_fp, CASEINSEN, "CAL_FILE", &arch_colnum, &status);
    fits_get_colnum(srch_fp, CASEINSEN, "CAL_FILE", &srch_colnum, &status);
    fits_read_col(srch_fp, TSTRING, srch_colnum, 1, 1, 1, &nul, strval, &anul, &status);
    fits_write_col(arch_fp, TSTRING, arch_colnum, 1, 1, 1, strval, &status);

    fits_get_colnum(arch_fp, CASEINSEN, "RFI_MTHD", &arch_colnum, &status);
    fits_get_colnum(srch_fp, CASEINSEN, "RFI_MTHD", &srch_colnum, &status);
    fits_read_col(srch_fp, TSTRING, srch_colnum, 1, 1, 1, &nul, strval, &anul, &status);
    fits_write_col(arch_fp, TSTRING, arch_colnum, 1, 1, 1, strval, &status);

    fits_get_colnum(arch_fp, CASEINSEN, "IFR_MTHD", &arch_colnum, &status);
    fits_get_colnum(srch_fp, CASEINSEN, "IFR_MTHD", &srch_colnum, &status);
    fits_read_col(srch_fp, TSTRING, srch_colnum, 1, 1, 1, &nul, strval, &anul, &status);
    fits_write_col(arch_fp, TSTRING, arch_colnum, 1, 1, 1, strval, &status);
    fits_report_error(stderr, status);
}

void get_params(struct arch_params *ap)
{
    int status = 0;
    fits_read_key(ap->fp, TSTRING, "STT_DATE", ap->stt_date, NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TSTRING, "STT_TIME", ap->stt_time, NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TSTRING, "DATE", ap->date, NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TSTRING, "SRC_NAME", ap->src_name, NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TSTRING, "RA", ap->ra, NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TSTRING, "DEC", ap->dec, NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TSTRING, "STT_CRD1", ap->stt_crd1, NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TSTRING, "STT_CRD2", ap->stt_crd2, NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TSTRING, "STP_CRD1", ap->stp_crd1, NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TSTRING, "STP_CRD2", ap->stp_crd2, NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TSTRING, "TRK_MODE", ap->trk_mode, NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TSTRING, "FD_MODE", ap->fd_mode, NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TSTRING, "CAL_MODE", ap->cal_mode, NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TFLOAT, "FA_REQ", &(ap->fa_req), NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TINT, "SCANLEN", &(ap->scanlen), NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TINT, "STT_IMJD", &(ap->imjd), NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TINT, "STT_SMJD", &(ap->smjd), NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TINT, "NRCVR", &(ap->nrcvr), NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TINT, "BE_PHASE", &(ap->be_phase), NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TINT, "BE_DCC", &(ap->be_dcc), NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TINT, "BE_DELAY", &(ap->be_delay), NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TINT, "TCYCLE", &(ap->tcycle), NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TINT, "OBSNCHAN", &(ap->obsnchan), NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TDOUBLE, "STT_OFFS", &(ap->offs), NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TDOUBLE, "STT_LST", &(ap->lst), NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TDOUBLE, "BMAJ", &(ap->bmaj), NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TDOUBLE, "BMIN", &(ap->bmin), NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TDOUBLE, "BPA", &(ap->bpa), NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TDOUBLE, "ANT_X", &(ap->ant_x), NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TDOUBLE, "ANT_Y", &(ap->ant_y), NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TDOUBLE, "ANT_Z", &(ap->ant_z), NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TFLOAT, "OBSFREQ", &(ap->obsfreq), NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TFLOAT, "EQUINOX", &(ap->equinox), NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TFLOAT, "OBSBW", &(ap->obsbw), NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TFLOAT, "CAL_FREQ", &(ap->cal_freq), NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TFLOAT, "CAL_DCYC", &(ap->cal_dcyc), NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TFLOAT, "CAL_PHS", &(ap->cal_phs), NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TSTRING, "BECONFIG", ap->beconfig, NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TSTRING, "OBS_MODE", ap->obs_mode, NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TSTRING, "DATE-OBS", ap->date_obs, NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TSTRING, "OBSERVER", ap->observer, NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TSTRING, "PROJID", ap->proj_id, NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TSTRING, "TELESCOP", ap->telescope, NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TSTRING, "FRONTEND", ap->frontend, NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TSTRING, "BACKEND", ap->backend, NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TSTRING, "FD_POLN", ap->fd_poln, NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TSTRING, "FD_HAND", ap->fd_hand, NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TSTRING, "FD_SANG", ap->fd_sang, NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TSTRING, "FD_XYPH", ap->fd_xyph, NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TSTRING, "COORD_MD", ap->coord_md, NULL, &status);
    if (status)
        status = 0;
    fits_read_key(ap->fp, TSTRING, "SCALE", ap->scale, NULL, &status);
    if (status)
        status = 0;

    fits_movnam_hdu(ap->fp, BINARY_TBL, "SUBINT", 1, &status);
    if (status) {
        printf("Unable to move to subint table in FITS file\n");
        exit(1);
    }

    long numrows;
    char chout[100];
    char comment[1000];
    fits_get_num_rows(ap->fp, &numrows, &status);
    ap->nsub = (int)numrows;

    fits_read_keyword(ap->fp, "NSBLK", chout, comment, &status);
    sscanf(chout,"%d",&(ap->nsamp));

    fits_read_keyword(ap->fp, "NCHAN", chout, comment, &status);
    sscanf(chout,"%d",&(ap->nchan));

    fits_read_keyword(ap->fp, "NBITS", chout, comment, &status);
    sscanf(chout, "%d", &(ap->nbits));

    fits_read_keyword(ap->fp, "NPOL", chout, comment, &status);
    sscanf(chout, "%d", &(ap->npol));

    fits_read_keyword(ap->fp, "TBIN", chout, comment, &status);
    sscanf(chout,"%f", &(ap->tbin));

    int colnum;

    fits_read_key(ap->fp, TDOUBLE, "CHAN_BW", &(ap->chan_bw), NULL, &status);
    fits_get_colnum(ap->fp, CASEINSEN, "DAT_FREQ", &colnum, &status);
    fits_read_col(ap->fp, TFLOAT, colnum, 1, 1, 1, NULL, &(ap->freq0),
            NULL, &status);

    fits_get_colnum(ap->fp, CASEINSEN, "LST_SUB", &colnum, &status);
    fits_read_col(ap->fp, TDOUBLE, colnum, 1, 1, ap->nsub, NULL, ap->lst_subs,
            NULL, &status);
}

void reverse_channels(unsigned char values[], unsigned char reversed_values[],
        int nbits, int nchan)
{
    int i, k = 0;
    for (i = nchan/8/nbits - 1; i >= 0; i--) {
        reversed_values[k] = values[i];
        k++;
    }
}

void minmax(float profile[], int nbin, float *min, float *max)
{
    *min = profile[0];
    *max = profile[0];

    int i;
    for (i = 0; i < nbin; ++i) {
        if (profile[i] < *min)
            *min = profile[i];
        else if (profile[i] > *max)
            *max = profile[i];
    }
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

//void resetProfile(float profile[MAX_CHANS][MAX_BINS], const uint nchan,
void resetProfile(float profile[MAX_POLNS][MAX_CHANS][MAX_BINS], const uint npol, const uint nchan, const uint nbin)
{
    uint ichan;
    uint ibin;
    uint ipol;

    for (ipol = 0; ipol < npol; ++ipol) {
        for (ichan = 0; ichan < nchan; ++ichan) {
            for (ibin = 0; ibin < nbin; ++ibin) {
                //profile[ichan][ibin] = 0.0;
                profile[ipol][ichan][ibin] = 0.0;
            }
        }
    }
}

void resetValues(float values[MAX_CHANS], const uint nchan)
{
    uint ichan;
    for (ichan = 0; ichan < nchan; ++ichan) {
        values[ichan] = 0.0;
    }
}



