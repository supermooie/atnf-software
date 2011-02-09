// Compile with:  g++ -o addsearch addsearch.C -lcfitsio -I.
/**
 * addsearch: a c++ program to append 1+ PSRFITS search-mode files to a select file.
 *     It requires two or more filenames provided as command-line arguments.
 *
 * Author: Jonathan Khoo
 * Date: 10.06.10
 */

#include <unistd.h>

#include <iostream>
#include <vector>

#include "fitsio.h"

using std::cerr;
using std::cout;
using std::endl;
using std::pair;
using std::vector;
using std::string;

const double SECONDS_IN_A_DAY = 86400.0;

void usage()
{
  cout <<
    "Usage: addsearch [options] <filename 1> <filename 2> [filenames] \n" <<
    "\n" <<
    "   -E ext     Extension of the output filename (default: <file>.add). \n" <<
    "   -h         This very helpful block of text.\n";
}

void SortFiles(vector<string>& filenames);
bool FileSortPredicate(const string& f1, const string& f2);
double GetMJD(const string& filename);

int main(int argc, char *argv[])
{
  string output_file_extension;
  const char* args = "E:h";

  int gotc = 0;
  while ((gotc = getopt(argc, argv, args)) != -1) {
    switch (gotc) {
      case 'E':
        output_file_extension = optarg;
        break;
      case 'h':
        usage();
        return 0;
        break;
      default:
        cerr << "default" << endl;
        break;
    }
  }

  // Extract each filename from the command line.
  vector<string> filenames;
  for (unsigned i = optind; i < argc; ++i) {
    filenames.push_back(argv[i]);
  }

  if (filenames.size() < 2) {
    usage();
    exit(0);
  }

  // Sort the files by start time before appending in case user is lazy.
  SortFiles(filenames);

  string total_filename;
  string total_extension;

  if (output_file_extension.empty()) {
    output_file_extension = "add"; // Default extension = <filename>.add
  }

  const size_t i = filenames[0].find_last_of(".");
  if (i != string::npos) {
    total_filename = filenames[0].substr(0, i) +
      "." +
      output_file_extension;
  }

  cerr << "Writing to: " << total_filename << endl;

  // Prepend '!' - overwrites any existing file with the same output filename.
  total_filename = "!" + total_filename;

  fitsfile* total_fp;
  fitsfile* fp;

  int status = 0;

  // Create new fitsfile with the new output filename and copy the contents
  // of the first (ordered by time) file.
  fits_create_file(&total_fp, (total_filename.c_str()), &status);
  fits_open_file(&fp, filenames[0].c_str(), READONLY, &status);

  fits_copy_file(fp, total_fp, 1, 1, 1, &status);
  fits_close_file(fp, &status);

  // Move to SUBINT table in joined file
  fits_movnam_hdu(total_fp, BINARY_TBL, "SUBINT", 0, &status);

  // Loop over remaining input files (2nd -> last file) and append
  // their SUBINT rows to the total file.
  for (unsigned i = 1; i < filenames.size(); ++i) {
    cerr << "Opening: " << filenames[i] << endl;
    fits_open_file(&fp, filenames[i].c_str(), READONLY, &status);

    // Move to the subint table.
    fits_movnam_hdu(fp, BINARY_TBL, "SUBINT", 0, &status);

    long nrows_in_total = 0;
    fits_get_num_rows(total_fp, &nrows_in_total, &status);

    long nrows_in_current = 0;
    fits_get_num_rows(fp, &nrows_in_current, &status);

    // Increase the total number of rows by however many rows in the current
    // file.
    fits_insert_rows(total_fp, nrows_in_total, nrows_in_current,
        &status);

    // Width of SUBINT table in bytes.
    int table_width;
    fits_read_key(fp, TINT, "NAXIS1", &table_width, NULL, &status);

    // Bytes to be appended to the total file is the number of bytes
    // in SUBINT table of the currently opened file.
    long bytes_to_be_appended = table_width * nrows_in_current;
    vector<unsigned char> data(bytes_to_be_appended);

    // Read the whole table, then write it to the allocated end portion
    // of the total file's SUBINT table.
    fits_read_tblbytes(fp, 1, 1, bytes_to_be_appended, &data[0], &status);
    fits_write_tblbytes(total_fp, nrows_in_total+1, 1, bytes_to_be_appended,
        &data[0], &status);

    fits_close_file(fp, &status);
  }

  fits_close_file(total_fp, &status);

  if (status) {
    fits_report_error(stderr, status);
  }

  return 0;
}

/**
 * Sort the given filenames into chronological order based on their start
 * time.
 */
void SortFiles(vector<string>& filenames)
{
  std::sort(filenames.begin(), filenames.end(), FileSortPredicate);
}

/**
 * Compare two files by their start MJD.
 */
bool FileSortPredicate(const string& f1, const string& f2)
{
  std::pair<double, double> MJDs(GetMJD(f1), GetMJD(f2));

  return MJDs.first < MJDs.second;
}

/**
 * Performs the usual MJD-related values mashup to determine the start time
 * in seconds.
 */
double GetMJD(const string& filename)
{
  fitsfile* fp;
  int status = 0;

  fits_open_file(&fp, filename.c_str(), READONLY, &status);

  int stt_imjd = 0;
  fits_read_key(fp, TINT, "STT_IMJD", &stt_imjd, NULL, &status);

  int stt_smjd = 0;
  fits_read_key(fp, TINT, "STT_SMJD", &stt_smjd, NULL, &status);

  double stt_offs = 0;
  fits_read_key(fp, TDOUBLE, "STT_OFFS", &stt_offs, NULL, &status);

  fits_close_file(fp, &status);

  if (status) {
    fits_report_error(stderr, status);
  }

  const double mjd = static_cast<double>(stt_imjd) +
    (static_cast<double>(stt_smjd) + stt_offs) / SECONDS_IN_A_DAY;

  return mjd;
}
