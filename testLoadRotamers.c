//cc testLoadRotamers.c rotamers.c -o testrota

#include <stdio.h>
#include "rotamers.h"

int main(void)
{
  intialize_AASCRotTable_from_file();
  int ind = getSideChainTemplateIndexFromName("ILE");
  printf("index %d for ILE", ind);
}
