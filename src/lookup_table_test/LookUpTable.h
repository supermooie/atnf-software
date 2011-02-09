#define NUMBER_OF_BITS 4

#include <vector>

class LookUpTable
{
  public:
    LookUpTable();
    float get_value(const unsigned n);

  private:
    void init();
    std::vector<float> table_values;
};
