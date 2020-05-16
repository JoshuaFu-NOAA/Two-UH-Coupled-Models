#include "ini_codes.h"

/* *INDENT-OFF* */
/*INDENT OFF*/

Table Table_surfini[MAXCODE_PTsurfini] = {
  { 129, "GEOSP",  "surface geopotential (orography)",         "m**2/s**2"},
  { 139, "TS",     "surface temperature",                      "K"},
  { 140, "WS",     "soil wetness",                             "m"},
  { 193, "WL",     "skin reservoir content",                   "m"},
  { 141, "SN",     "snow depth",                               "m"},
  { 172, "SLM",    "land sea mask",                            "fractional"},
  { 173, "AZ0",    "surface roughness length",                 "m"},
  { 174, "ALB",    "surface background albedo",                "fractional"},
  {  89, "VARP",   "directional orographic variance (packed)", "m**2"},
  {  90, "OVEW",   " E-W  orographic variance",                "m**2"},
  {  91, "OVNS",   " N-S  orographic variance",                "m**2"},
  {  92, "OVNW",   "NW-SE orographic variance",                "m**2"},
  {  93, "OVNE",   "NE-SW orographic variance",                "m**2"},
  { 199, "VAROR",  "orographic variance",                      "m**2"},
  { 212, "FOREST", "vegetation type",                          "fractional"},
  { 198, "VGRAT",  "vegetation ratio",                         "fractional"},
  { 200, "VLT",    "leaf area index",                          " "},
  { 229, "WSMX",   "field capacity of soil",                   "m"},
  { 226, "FAO",    "FAO data set (soil data flags) [0...5.]",  " "},
  { 232, "GLAC",   "glacier mask",                             "fractional"},
  {  99, "ALAKE",  "lake mask",                                "fractional"},
  {  51, "OROMEA", "Mean orography",                           "m"},
  {  52, "OROSTD", "Orographic standard deviation",            "m"},
  {  53, "OROSIG", "Orographic slope",                         "degree"},
  {  54, "OROGAM", "Orographic anisotropy",                    "degree"},
  {  55, "OROTHE", "Orographic angle",                         "degree"},
  {  56, "OROPIC", "Orographic peaks elevation",               "m"},
  {  57, "OROVAL", "Orographic valleys elevation",             "m"},
};

int MaxCode_surfini = sizeof (Table_surfini) / sizeof (Table);
 
Table Table_specini[] = {
  { 138, "SVO",    "vorticity",           "1/s"},
  { 155, "SD",     "divergence",          " "},
  { 130, "STP",    "temperature",         "K"},
  { 133, "Q",      "specific humidity",   "kg/kg"}
};
 
int MaxCode_specini = sizeof (Table_specini) / sizeof (Table);

/*
int PT_getidx_code(Table *ptable, int code)
{
    int indx;

    for (indx = 0; indx < MAXCODE_PTsurfini; indx++)
      {
	if (ptable->code[indx] == code) return (indx);
      }

    return (-1);
}
*/
