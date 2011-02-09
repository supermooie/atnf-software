// g++ -o fix fix.C -lm -L/usr/local/src/cfitsio -lcfitsio -I.                  32 bit
// g++ -o fix fix.C -lm -L/pulsar/psr/linux_64/share/ -lcfitsio -I.             64 bit

/* DTV rfi*/

#include <iostream>
#include <stdlib.h>
#include "fitsio.h"
#include <vector>
#include <cmath>

#define MAXCHARS 512 
#define MAXCHANS 2048

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

struct tv_freq_range {
	float lower_freq;
	float upper_freq;
	double start;
};

const struct tv_freq_range ABC = {652.0, 659.0, 52852};
const struct tv_freq_range Prime = {659.0, 666.0, 54313};
const struct tv_freq_range SBS = {666.0, 673.0, 53310};
const struct tv_freq_range WIN = {680.0, 687.0, 54282};
const struct tv_freq_range Sthn_Cross = {687.0, 694.0, 54313};

	/*							Commencement
Channel  Station   Freq range    Date     MJD
  46       ABC    652.0-659.0   1/8/03   52852
  47      Prime   659.0-666.0   1/8/07   54313
  48       SBS    666.0-673.0   1/11/04  53310
  49       --
  50       WIN    680.0-687.0   1/7/07   54282
  51  Sthn Cross  687.0-694.0   1/8/07   54313
*/

const double MJD_DTV_1_8_03 = 52852;
const double MJD_DTV_1_8_07 = 54313;
const double MJD_DTV_1_11_04 = 53310;
const double MJD_DTV_1_7_07 = 54282;

const double chan_32 = 0.0000012969;
const double chan_128 = 0.0000035547;
const double chan_512 = 0.0000122813;
const double chan_2048 = 0.0000468829;
const double MJD_DFB1_fix_date = 53954.621527777985;

const double MJD_DFB2_5_6_07 = 54256;
const double MJD_DFB2_21_6_07 = 54272;
const double MJD_DFB2_2_12_07 = 54143;
const double MJD_DFB2_16_3_07 = 54175;

const double d1 = 0.000014;
const double d2_256bw_512chan = 0.00001391;
const double d2_256bw_1024chan = 0.00001720;
const double d2_256bw_2048chan = 0.00004315;

const string DFB1_OFFS_CMD = "fitsfix dfb1 offs";
const string DFB2_OFFS_CMD = "fitsfix dfb2 offs";
const string DFB_BW_CMD = "fitsfix dfb bw";

fitsfile *fp;


// fixing methods -- modify FITSArchive

void dfb1_offs();
void dfb2_offs();
void correctOffset(double stt_offs);
void reverse_freqs();
void dfb1_2_bw();
void be_phase();
void CPSR2_names(const char first_letter);
void feed_pol_params();
void replace_receiver();
void write_history(const string cmd);
void DTV_channels(const char first_letter);
void backend_fix(const char first_letter);
void CPSR2_receiver(const char first_letter);

bool isCorrected(const string cmd);
// get methods -- read various headers from FITSfile 

int get_fdhand();
int get_fdsang();
int get_fdxpyh();
unsigned int get_orig_nchan();
string get_receiver();
string get_backend();
string get_beconfig();
string get_source();
double get_delay();
float get_obsfreq();
float get_chan_bw();
float get_obsbw();
void get_MJD(int *imjd, double *stt_offs);
string get_fdpoln();
unsigned int get_nchan();

void zero_weight(unsigned chan);
bool in_range(float freq, float lower, float upper);

void usage(const string prog_name)
{
	cerr << "\nUsage: " << prog_name << " [filename]\n\n";
	exit(0);
}

int main(int argc, char* argv[])
{
	if (argc != 2)
		usage(argv[0]);

	const string filename = argv[1];

	int status = 0;
	fits_open_file(&fp, filename.c_str(), READWRITE, &status);
	fits_report_error(stderr, status);

	backend_fix(filename.at(0));
	dfb1_2_bw();
	dfb1_offs();
	//dfb2_offs(); // removed - Dick doesn't remember what this was for
	replace_receiver();
	CPSR2_names(filename.at(0));
	//	printf("Calling phase\n");
	//	printf("Calling receiver\n");
	CPSR2_receiver(filename.at(0));
	//	be_phase();
	//DTV_channels(filename.at(0)); // removed - zero-weighting is performed elsewhere in the pipeline

	feed_pol_params();

	fits_close_file(fp, &status);
	fits_report_error(stderr, status);
	return 0;
}

void DTV_channels(const char first_letter)
{
	if (first_letter != 'm' && first_letter != 'n')
		return;

	const float obsbw = get_obsbw();
	const float ctr_freq = get_obsfreq();

    vector<struct tv_freq_range> tv_stations;
    tv_stations.push_back(ABC);
    tv_stations.push_back(Prime);
    tv_stations.push_back(SBS);
    tv_stations.push_back(WIN);
    tv_stations.push_back(Sthn_Cross);

    double stt_offs;
    int imjd;
	get_MJD(&imjd, &stt_offs);
	double start_time = (double)imjd + stt_offs;
	const unsigned int nchan = get_nchan();
	const float first_freq = ctr_freq - obsbw/2.0;

    float chan_freqs[nchan];
    float freq = first_freq;
	const float chan_bw = get_chan_bw();

    for (unsigned i = 0; i < nchan; ++i) {
        chan_freqs[i] = freq;
        freq += chan_bw;
    }

    for (unsigned istation = 0; istation < tv_stations.size(); ++istation) {
        if (start_time >= tv_stations[istation].start) {
            for (unsigned i = 0; i < nchan; ++i) {
                if (in_range(chan_freqs[i], tv_stations[istation].lower_freq, tv_stations[istation].upper_freq)) {
                    zero_weight(i);
                }
            }
        }
    }
}

void zero_weight(unsigned chan)
{
	int status = 0;
	long numrows;
	int colnum = 0;
	const unsigned int nchan = get_nchan();
	float weights[MAXCHANS];

	fits_get_colnum (fp, CASEINSEN, "DAT_WTS", &colnum, &status);
	fits_get_num_rows(fp, &numrows, &status);
	fits_read_col(fp, TFLOAT, colnum, 1, 1, nchan, NULL, weights, NULL, &status);

	weights[chan] = 0;

	for (int irow = 1; irow <= numrows; ++irow)
		fits_write_col(fp, TFLOAT, colnum, irow, 1, nchan, weights, &status);
}

bool in_range(float freq, float lower, float upper)
{
	return freq >= lower && freq <= upper;
}


void CPSR2_names(const char first_letter)
{
	if (first_letter != 'm' && first_letter != 'n')
		return;

	string source(get_source());

	if (source[0] != 'J') {
		int status = 0;
		source = "J" + source;
		char *tempstr = const_cast<char*>(source.c_str());
		fits_update_key(fp, TSTRING, "SRC_NAME", tempstr, NULL, &status);
	}
}

double get_delay()
{
	int status = 0;
	double be_delay;
    fits_movabs_hdu(fp, 1, NULL, &status);
	fits_read_key(fp, TDOUBLE, "BE_DELAY", &be_delay, NULL, &status);

	if (status)
		return 0;

	return be_delay;
}

string get_source()
{
	int status = 0;
	char value[MAXCHARS];

    fits_movabs_hdu(fp, 1, NULL, &status);
	fits_read_key(fp, TSTRING, "SRC_NAME", value, NULL, &status);

	if (status)
		return "";

	return value;
}

void be_phase()
{
	const string backend(get_backend());
	const string beconfig(get_beconfig());
	const string last_two(beconfig.substr(beconfig.length()-2, 2));

	int status = 0;
    fits_movabs_hdu(fp, 1, NULL, &status);

	int be_phase = (backend == "WBCORR" && last_two == "_c") ? 1 : -1;
	fits_update_key(fp, TINT, "BE_PHASE", &be_phase, NULL, &status);
}

string get_beconfig()
{
	int status = 0;
	char value[MAXCHARS];

    fits_movabs_hdu(fp, 1, NULL, &status);
	fits_read_key(fp, TSTRING, "BECONFIG", value, NULL, &status);

	if (status)
		return "";

	return value;
}

void backend_fix(const char first_letter)
{
	if (first_letter != 'a')
		return;

	string backend(get_backend());

	if (backend != "WBCORR")
		return;

	backend = "PDFB1";
	char *tempstr = const_cast<char*>(backend.c_str());

	int status = 0;
    fits_movabs_hdu(fp, 1, NULL, &status);
	fits_update_key(fp, TSTRING, "BACKEND", tempstr, NULL, &status);
}

void dfb1_offs()
{
    double stt_offs;
    int imjd;

	get_MJD(&imjd, &stt_offs);
	string backend(get_backend());
	double start_time = (double)imjd + stt_offs;


	if (backend == "PDFB1" && start_time < MJD_DFB1_fix_date && !isCorrected(DFB1_OFFS_CMD))
		correctOffset(stt_offs);
}

void dfb2_offs()
{
	const string backend(get_backend());

	if (backend != "PDFB2" || isCorrected(DFB2_OFFS_CMD))
		return;

    double stt_offs;
    int imjd;

	get_MJD(&imjd, &stt_offs);
	double start_time = (double)imjd + stt_offs;

	unsigned int nchan = get_orig_nchan();
	float bw = get_obsbw();

	cerr << "----------------------\n";
	cerr << "Start time: " << start_time << endl;
	cerr << "Original nchan: " << nchan << endl;
	cerr << "BW: " << bw << endl;

	if (bw != -256 || --nchan != 1024)
		return;

	cerr << "Aboat to change STT_OFFS" << endl;
	
	if (start_time >= MJD_DFB2_16_3_07 && start_time < MJD_DFB2_21_6_07) {
		stt_offs -= d2_256bw_1024chan;
		cerr << "Changed STT_OFFS by -" << d2_256bw_1024chan << endl;
	}

	/*if (start_time < MJD_DFB2_5_6_07) {
		if (nchan == 32)
			stt_offs += chan_32;
		else if (nchan == 128)
			stt_offs += chan_128;
		else if (nchan == 512)
			stt_offs += chan_512;
		else if (nchan == 2048)
			stt_offs += chan_2048;
		else
			cerr << "WARNING --- nchan: != 32, 128, 512 or 2058" << endl;

		stt_offs -= d1;

		if (fabs(bw) != 256)
			return;

		if (nchan == 512)
			stt_offs -= d2_256bw_512chan;
		else if (nchan == 1024)
			stt_offs -= d2_256bw_1024chan;
		else if (nchan == 2048)
			stt_offs -= d2_256bw_2048chan;
		else
			return;

	} else if (start_time > MJD_DFB2_5_6_07 && start_time < MJD_DFB2_21_6_07) {

		if (nchan == 32)
			stt_offs += chan_32;
		else if (nchan == 128)
			stt_offs += chan_128;
		else if (nchan == 512)
			stt_offs += chan_512;
		else if (nchan == 2048)
			stt_offs += chan_2048;
		else
			cerr << "WARNING --- nchan: != 32, 128, 512 or 2058" << endl;

		stt_offs -= d1;

		if (fabs(bw) != 256)
			return;
		if (nchan == 512)
			stt_offs -= 0.5*d2_256bw_512chan;
		else if (nchan == 1024)
			stt_offs -= 0.5*d2_256bw_1024chan;
		else if (nchan == 2048)
			stt_offs -= 0.5*d2_256bw_2048chan;
		else
			return;
	} else if (start_time > MJD_DFB2_21_6_07 && start_time < MJD_DFB2_2_12_07) {
		if (nchan == 512)
			stt_offs += 0.5*d2_256bw_512chan;
		else if (nchan == 1024)
			stt_offs += 0.5*d2_256bw_1024chan;
		else if (nchan == 2048)
			stt_offs += 0.5*d2_256bw_2048chan;
		else
			return;

	} else if (start_time > MJD_DFB2_2_12_07) {

		double be_delay = get_delay();
		stt_offs += 0.5*be_delay;
	} else
		return;*/

	int status = 0;
    fits_movabs_hdu(fp, 1, NULL, &status);
	fits_update_key(fp, TDOUBLE, "STT_OFFS", &stt_offs, "", &status);
	write_history(DFB2_OFFS_CMD);
}

void dfb1_2_bw()
{
	string backend(get_backend());
	if (backend != "PDFB1" && backend != "PDFB2")
		return;

	string frontend(get_receiver());
	if (frontend != "MULTI")
		return;

	float obsbw = get_obsbw();

	if (obsbw != 256)
		return;

	int status = 0;
	obsbw = -256;
    fits_movabs_hdu(fp, 1, NULL, &status);
	fits_update_key(fp, TINT, "OBSBW", &obsbw, NULL, &status);

	fits_movnam_hdu(fp, BINARY_TBL, "SUBINT", 0, &status);
	float chan_bw = get_chan_bw();

	if (chan_bw > 0) {
		chan_bw *= -1.0;
		fits_update_key(fp, TDOUBLE, "CHAN_BW", &chan_bw, NULL, &status);
	}

	reverse_freqs();
	write_history(DFB_BW_CMD);
}

void feed_pol_params()
{
	string receiver(get_receiver());
	float obsfreq = get_obsfreq();
	string fdpoln(get_fdpoln());
	int fdhand = get_fdhand();
	int fdsang = get_fdsang();
	int fdxpyh = get_fdxpyh();

	cerr << "feed pol params" << endl;

	/*cerr << "fdpoln: " << fdpoln << endl;
	cerr << "fdhand: " << fdhand << endl;
	cerr << "fdsang: " << fdsang << endl;
	cerr << "fdxpyh: " << fdxpyh << endl;*/

	unsigned int frontend;

	cerr << "receiver: " << receiver << endl;

	if (receiver == "H-OH")
		frontend = 1;
	else if (receiver == "MULTI")
		frontend = 2;
	else if (receiver == "1050CM")
		frontend = 3;
	else
		return;

	int status = 0;
    fits_movabs_hdu(fp, 1, NULL, &status);

	if (status)
		return;

	int itemp;
    char *lin = "lin";

	switch (frontend) {
		case 1:
			cerr << "HOH" << endl;
			itemp = 0;
			fits_update_key(fp, TINT, "OBSFREQ", &itemp, NULL, &status);
			status = 0;


			fits_update_key(fp, TSTRING, "FD_POLN", lin, NULL, &status);
			status = 0;

			itemp = 1;
			fits_update_key(fp, TINT, "FD_HAND", &itemp, NULL, &status);
			status = 0;

			itemp = -90;
			fits_update_key(fp, TINT, "FD_SANG", &itemp, NULL, &status);
			status = 0;

			itemp = 180;
			fits_update_key(fp, TINT, "FD_XPYH", &itemp, NULL, &status);
			fits_report_error(stderr, status);
			status = 0;
			break;

		case 2:
			itemp = 0;
			fits_update_key(fp, TINT, "OBSFREQ", &itemp, NULL, &status);
			fits_update_key(fp, TSTRING, "FD_POLN", lin, NULL, &status);
			itemp = 1;
			fits_update_key(fp, TINT, "FD_HAND", &itemp, NULL, &status);
			itemp = -90;
			fits_update_key(fp, TINT, "FD_SANG", &itemp, NULL, &status);
			itemp = 0;
			fits_update_key(fp, TINT, "FD_XPYH", &itemp, NULL, &status);
			break;

		case 3:
			/*if (obsfreq != 0) {
				itemp = 0;
				fits_update_key(fp, TINT, "OBSFREQ", &itemp, NULL, &status);
			}*/

			fits_update_key(fp, TSTRING, "FD_POLN", lin, NULL, &status);

			itemp = -1;
			fits_update_key(fp, TINT, "FD_HAND", &itemp, NULL, &status);

			itemp = 0;
			fits_update_key(fp, TINT, "FD_SANG", &itemp, NULL, &status);

			itemp = 0;
			fits_update_key(fp, TINT, "FD_XPYH", &itemp, NULL, &status);
			break;
	}
}

int get_fdxpyh()
{
	int status = 0;
	int fdxpyh;

    fits_movabs_hdu(fp, 1, NULL, &status);
	fits_read_key(fp, TINT, "FD_XPYH", &fdxpyh, NULL, &status);

	if (status)
		return 0;

	return fdxpyh;
}

int get_fdsang()
{
	int status = 0;
	int fdsang;

    fits_movabs_hdu(fp, 1, NULL, &status);
	fits_read_key(fp, TINT, "FD_SANG", &fdsang, NULL, &status);

	if (status)
		return 0;

	return fdsang;
}

string get_fdpoln()
{
	int status = 0;
	char value[MAXCHARS];

    fits_movabs_hdu(fp, 1, NULL, &status);
	fits_read_key(fp, TSTRING, "FD_POLN", value, NULL, &status);

	if (status)
		return "";

	return value;
}

int get_fdhand()
{
	int status = 0;
	int fdhand;

    fits_movabs_hdu(fp, 1, NULL, &status);
	fits_read_key(fp, TINT, "FD_HAND", &fdhand, NULL, &status);

	if (status)
		return 0;

	return fdhand;
}


float get_obsfreq()
{
	int status = 0;
	float obsfreq;
	int colnum = 0;
	long numrows;

	fits_movnam_hdu(fp, BINARY_TBL, "HISTORY", 0, &status);
	fits_report_error(stderr, status);
	fits_get_colnum(fp, CASEINSEN, "CTR_FREQ", &colnum, &status);
	fits_report_error(stderr, status);
	fits_get_num_rows(fp, &numrows, &status);
	fits_report_error(stderr, status);

//	fits_read_col(fp, TFLOAT, colnum, 1, 1, --numrows, NULL, &obsfreq, NULL, &status);
	fits_read_col(fp, TFLOAT, colnum, 1, 1, 1, NULL, &obsfreq, NULL, &status);
	//fits_report_error(stderr, status);

	/*if (status)
		return 0;*/

	return obsfreq;
}

void reverse_freqs()
{
	int status = 0;
	long numrows;
	int colnum = 0;

	fits_get_colnum (fp, CASEINSEN, "DAT_FREQ", &colnum, &status);
	fits_get_num_rows(fp, &numrows, &status);

	unsigned int nchan = get_nchan();

	float orig_freqs[MAXCHANS];
	float reversed_freqs[MAXCHANS];

	fits_read_col(fp, TFLOAT, colnum, 1, 1, nchan, NULL, orig_freqs, NULL, &status);

	for (unsigned int ichan = 0; ichan < nchan; ++ichan)
		reversed_freqs[ichan] = orig_freqs[nchan-ichan-1];

	for (int irow = 1; irow <= numrows; ++irow)
		fits_write_col(fp, TFLOAT, colnum, irow, 1, nchan, reversed_freqs, &status);
	fits_report_error(stderr, status);
}

float get_chan_bw()
{
	float chan_bw = 0;
	int status = 0;

	fits_movnam_hdu(fp, BINARY_TBL, "SUBINT", 0, &status);
	fits_read_key(fp, TFLOAT, "CHAN_BW", &chan_bw, NULL, &status);

	if (status)
		return 0;

	return chan_bw;
}

unsigned int get_nchan()
{
	int status = 0;
	unsigned int nchan = 0;

	fits_movnam_hdu(fp, BINARY_TBL, "SUBINT", 0, &status);
	fits_report_error(stderr, status);
	fits_read_key(fp, TINT, "NCHAN", &nchan, NULL, &status);

	if (status)
		return 0;

	return nchan;
}

float get_obsbw()
{
	float chan_bw = get_chan_bw();
	unsigned int nchan = get_nchan();

	return chan_bw * (float)nchan;
}

string get_backend()
{
	int status = 0;
	char value[MAXCHARS];

    fits_movabs_hdu(fp, 1, NULL, &status);
	fits_read_key(fp, TSTRING, "BACKEND", value, NULL, &status);
	fits_report_error(stderr, status);

	if (status)
		return "";

	return value;
}


void correctOffset(double stt_offs)
{
	unsigned int nchan = get_orig_nchan();

	if (!nchan)
		return;

	if (--nchan == 32)
		stt_offs -= chan_32;
	else if (nchan == 128)
		stt_offs -= chan_128;
	else if (nchan == 512)
		stt_offs -= chan_512;
	else if (nchan == 2048)
		stt_offs -= chan_2048;
	else
		cerr << "WARNING --- nchan: != 32, 128, 512 or 2058" << endl;

	int status = 0;
    fits_movabs_hdu(fp, 1, NULL, &status);
	fits_update_key(fp, TDOUBLE, "STT_OFFS", &stt_offs, "", &status);

	write_history(DFB1_OFFS_CMD);
}

void replace_receiver()
{
	string receiver(get_receiver());

	if (receiver == "MULT_1") {
		int status = 0;
		char *multi = "MULTI";
		fits_movabs_hdu(fp, 1, NULL, &status);
		fits_update_key(fp, TSTRING, "FRONTEND", multi, NULL, &status);
	}
}

string get_receiver()
{
	int status = 0;
	char value[MAXCHARS];

    fits_movabs_hdu(fp, 1, NULL, &status);
	fits_read_key(fp, TSTRING, "FRONTEND", value, NULL, &status);
	fits_report_error(stderr, status);

	if (status)
		return "";

	return value;
}

void get_MJD(int *imjd, double *stt_offs)
{
	int status = 0;
	fits_movabs_hdu(fp, 1, NULL, &status);
    fits_read_key(fp, TDOUBLE, "STT_OFFS", stt_offs, NULL, &status);
	fits_report_error(stderr, status);
    fits_read_key(fp, TINT, "STT_IMJD", imjd, NULL, &status);
	fits_report_error(stderr, status);
}

unsigned int get_orig_nchan()
{
    int status = 0;
    unsigned nchan = 0;

	fits_movnam_hdu(fp, BINARY_TBL, "BANDPASS", 0, &status);
	fits_read_key(fp, TINT, "NCH_ORIG", &nchan, NULL, &status);

	if (!status) {
		return nchan;
	} else {
		fits_report_error(stderr, status);
		return 0;
	}
}

void write_history(const string cmd)
{
	int status = 0;
	int colnum;
	long numrows = 0;

	char* strval = new char[MAXCHARS];

	fits_movnam_hdu(fp, BINARY_TBL, "HISTORY", 0, &status);
	fits_report_error(stderr, status);
	fits_get_colnum(fp, CASEINSEN, "PROC_CMD", &colnum, &status);
	fits_report_error(stderr, status);
	fits_get_num_rows(fp, &numrows, &status);
	fits_report_error(stderr, status);

	fits_insert_rows(fp, 1, 1, &status);
	fits_report_error(stderr, status);

	strcpy(strval, cmd.c_str());
	fits_write_col(fp, TSTRING, colnum, ++numrows, 1, 1, &strval, &status);
}

bool isCorrected(const string cmd)
{
	int status = 0;
	int colnum;
	long numrows;

	fits_movnam_hdu(fp, BINARY_TBL, "HISTORY", 0, &status);
	fits_report_error(stderr, status);
	fits_get_colnum(fp, CASEINSEN, "PROC_CMD", &colnum, &status);
	fits_report_error(stderr, status);
	fits_get_num_rows(fp, &numrows, &status);
	fits_report_error(stderr, status);

	char *proc_cmd  = new char[MAXCHARS];

	if (status) {
		cerr << "Error reading proc history -- terminating" << endl;
		exit(1);
	}

	for (int i = 1; i <= numrows; ++i) {
		fits_read_col(fp, TSTRING, colnum, i, 1, 1, NULL, &proc_cmd, NULL, &status);

		if (proc_cmd == cmd)
			return true;
	}

	return false;
}

void CPSR2_receiver(const char first_letter)
{
	if (first_letter != 'm' && first_letter != 'n')
		return;

	float ctr_freq = get_obsfreq();
	printf("Changing CPSR2 receiver: freq = %g\n",ctr_freq);
	int status = 0;
    fits_movabs_hdu(fp, 1, NULL, &status);

	if (ctr_freq >= 300 && ctr_freq <= 800) {
		char *rcvr = "50CM";
		fits_update_key(fp, TSTRING, "FRONTEND", rcvr, NULL, &status);
	}

	if (ctr_freq >= 2900 && ctr_freq <= 3700) {
		char *rcvr = "10CM";
		fits_update_key(fp, TSTRING, "FRONTEND", rcvr, NULL, &status);
	}

	if (ctr_freq >= 1300 && ctr_freq <= 1500) {
		double stt_offs;
		int imjd;

		get_MJD(&imjd, &stt_offs);
		double start_time = (double)imjd + stt_offs;

		string receiver;
		printf("Start time = %f\n",start_time);
		if (start_time >= 52980 && start_time <= 53250)
		  {
		    if (ctr_freq > 1300 && ctr_freq < 1360) receiver = "mH-OH"; else receiver = "nH-OH";
		  }
		else if (start_time >= 53334 && start_time <= 53328)
		  {
		    if (ctr_freq > 1300 && ctr_freq < 1360) receiver = "mH-OH"; else receiver = "nH-OH";
		  }
		else if (start_time >= 53835 && start_time <= 53837)
		  {
		    if (ctr_freq > 1300 && ctr_freq < 1360) receiver = "mH-OH"; else receiver = "nH-OH";
		  }
		else if (start_time >= 54065 && start_time <= 54228)
		  {
		    if (ctr_freq > 1300 && ctr_freq < 1360) receiver = "mH-OH"; else receiver = "nH-OH";
		  }
		else
		  {
		    if (ctr_freq > 1300 && ctr_freq < 1360) receiver = "mMULTI"; else receiver = "nMULTI";
		  }
		char *tempstr = const_cast<char*>(receiver.c_str());
		printf("Setting receiver = %s\n",tempstr);
		fits_update_key(fp, TSTRING, "FRONTEND", tempstr, NULL, &status);
		printf("COMPLETE Setting receiver\n");
	}
}
