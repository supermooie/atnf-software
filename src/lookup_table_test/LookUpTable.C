#include "LookUpTable.h"

LookUpTable::LookUpTable()
{
  init();
}

void LookUpTable::init()
{
  table_values.resize(NUMBER_OF_BITS);

  table_values[0] = -2.5;
  table_values[1] = -0.5;
  table_values[2] = 0.5;
  table_values[3] = 2.5;
}


float LookUpTable::get_value(const unsigned n)
{
  return table_values[n];
}

