#include "ini_codes.h"

/* *INDENT-OFF* */
/*INDENT OFF*/

char institute[] = "Max Planck Institut für Meteorologie";

Table Table_dim[] = {
  {  1, "lon",    "longitude",            "degrees_E"},
  {  2, "lat",    "Gaussian latitude",    "degrees_N"},
};

int MaxCode_dim = sizeof (Table_dim) / sizeof (Table);
