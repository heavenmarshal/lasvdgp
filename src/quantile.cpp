#include<cassert>
#include<algorithm>
#include "quantile.hpp"

double quantile(double* vin, double quan, int len)
{
  assert(len>0);
  assert(quan >= 0.0 && quan < 1.0);
  int pos = (int)(quan * len);
  std::nth_element(vin, vin+pos, vin+len);
  return vin[pos];
}
