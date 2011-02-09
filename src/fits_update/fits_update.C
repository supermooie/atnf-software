#include "Pulsar/psrchive.h"
#include "Pulsar/Archive.h"

using std::cerr;
using std::cout;
using std::endl;
using Pulsar::Archive;

int main(int argc, char* argv[])
{
	if (argc != 2) {
		cerr << argv[0] << " [psrfits file]" << endl;
		return 1;
	}

	try {
		Reference::To<Pulsar::Archive> archive = Archive::load(argv[1]);
		archive->unload();
	} catch (Error e) {
		cerr << e << endl;
	}

	return 0;
}
