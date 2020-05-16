#include "ini_codes.h"

/* *INDENT-OFF* */
/*INDENT OFF*/

Table Table_year[] = {
  {  70, "VGRATCLIM", "vegetation ratio",         "fractional"},
  {  71, "VLTCLIM",   "leaf area index",          " "},
  {  72, "TSLCLIM",   "land surface temperature", "K"},
  {  77, "sst",       "sea surface temperature",  "K"},
  {  78, "sic",       "sea ice concentration",    "%"},
  { 198, "VGRAT",     "vegetation ratio",         "fractional"},
  { 200, "VLT",       "leaf area index",          " "},
  {  92, "aflux",     "flux correction",          "W/m^2"},
  {  99, "year",      "annual data",              " "}
};

int MaxCode_year = sizeof (Table_year) / sizeof (Table);
