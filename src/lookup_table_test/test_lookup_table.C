#include <iostream>

#include <ctime>
#include <vector>

#include "LookUpTable.h"

using std::cerr;
using std::cout;
using std::clock_t;
using std::endl;
using std::vector;

void InitialiseInput(vector<unsigned>& v);

int main(int argc, char *argv[])
{
  LookUpTable table = LookUpTable();

  double diff;

  vector<unsigned> input(1024*1024*512);
  InitialiseInput(input);

  vector<float> values(1024*1024*512);

  // Start timer.
  clock_t start = std::clock();

  for (unsigned i = 0; i < input.size(); ++i) {
    values[i] = table.get_value(input[i]);
  }

  // End timer.
  diff = (std::clock() - start) / (double)CLOCKS_PER_SEC;

  cerr << "Time elapsed: " << diff << endl;

  return 0;
}

void InitialiseInput(vector<unsigned>& v)
{
  for (unsigned i = 0; i < v.size(); ++i) {
    v[i] = rand() % NUMBER_OF_BITS + 1;
  }
}

