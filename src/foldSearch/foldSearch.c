// gcc -g -o foldSearch foldSearch.c -lm -L$TEMPO2/lib -lcfitsio -I. -I$TEMPO2/include -lg2c -ltempo2pred

#include <stdio.h>
#include "search_lib.h"
#include "search_lib_t2.h"
#include <unistd.h>
#include <string.h>
#include "tempo2pred.h"

#define MAXROWS 1024
#define FILELENGTH 512

void usage()
{
    printf("foldSearch [options] filename\n");
    printf("Options:\n");
    printf("-b                      bins\n");
    printf("-e                      .eph file\n");
    printf("-f                      T2 predictor file\n");
    printf("-h                      help\n");
    printf("-p                      period (if no predictor file specified\n\n");
    exit(1);
}

float profile[MAX_POLNS][MAX_CHANS][MAX_BINS];

int main(int argc, char *argv[])
{
    float values[MAX_CHANS];
    struct arch_params old_ap;
    fitsfile *fp;
    int status = 0;
    int pol = 0;
    int nbin = 0;
    char *pred_file = malloc(FILELENGTH);
    char *eph_file = malloc(FILELENGTH);
    long double period = 0.0;

    int gotc = 0;
    while ((gotc = getopt(argc, argv, "b:e:f:hp:s:")) != -1) {
        switch (gotc) {
            case 'b':
                nbin = atoi(optarg);
                break;

            case 'e':
                eph_file = optarg;
                break;

            case 'f':
                pred_file = optarg;
                break;

            case 'h':
                usage();
                break;

            case 'p':
                period = atof(optarg);
                break;
        }
    }

    fits_open_file(&fp, argv[argc - 1], READONLY, &status);
    if (status) {
        printf("Unable to open file %s\n",argv[argc - 1]);
        exit(1);
    }

    old_ap.fp = fp;
    get_params(&old_ap);

    int col_num = 0;
    long numrows;
    float offs_sub[MAXROWS];
    float lst_sub[MAXROWS];

    fits_movnam_hdu(fp, BINARY_TBL, "SUBINT", 0, &status);
    fits_get_colnum(fp, CASEINSEN, "OFFS_SUB", &col_num, &status);
    fits_get_num_rows(fp, &numrows, &status);

    fits_read_col(fp, TFLOAT, col_num, 1, 1, numrows, NULL, offs_sub, NULL, &status);
    fits_get_colnum(fp, CASEINSEN, "LST_SUB", &col_num, &status);
    fits_get_num_rows(fp, &numrows, &status);
    fits_read_col(fp, TFLOAT, col_num, 1, 1, numrows, NULL, lst_sub, NULL, &status);

    fitsfile *afp = create_archive(old_ap, nbin, eph_file, argv[argc-1]);

    if (!nbin) {
        printf("Enter number of bins: ");
        scanf("%d", &nbin);
    }

    double turn = 0.0;
    int chan;
    int colnum;
    fits_get_colnum (fp, CASEINSEN, "DATA", &colnum, &status);
    int isub;
    int ibin;
    int isamp;
    int ipol;

    T2Predictor t2p = create_T2Predicator(pred_file);

    long double reference_phase = get_phase(old_ap, &t2p, 0, 0);

    for (isub = 1; isub <= old_ap.nsub; ++isub) {
        resetProfile(profile, old_ap.npol, old_ap.nchan, nbin);

        period = get_period(old_ap, &t2p, isub-1);
        const long double temp_phase = get_phase(old_ap, &t2p, isub-1, (double)offs_sub[isub-1]);

        printf("Subint %d/%d\n",isub, old_ap.nsub);
        for (isamp = 0; isamp < old_ap.nsamp; ++isamp) {
            const double phase = turn-floor(turn);
            for (ipol = 0; ipol < old_ap.npol; ++ipol) {
                getChans(fp, old_ap, isamp, isub, ipol+1, values, colnum);
                for (chan = 0; chan < old_ap.nchan; ++chan) {
                    ibin = (int)((double)nbin*phase);
                    profile[ipol][chan][ibin] += values[chan];
                }
            }
            turn += (long double)old_ap.tbin/period;
        }
        float archive_offs_sub = compute_offs_sub(reference_phase, temp_phase,
                period, offs_sub[isub-1]);

        make_subint(afp, old_ap, isub, nbin, profile, archive_offs_sub,
                (double)lst_sub[isub-1]);
    }
    //copy_history_params(afp, old_ap.fp);
    //remove_unused_tables(afp);
    fits_report_error(stderr, status);
    fits_close_file(afp, &status);
    return 0;
}

