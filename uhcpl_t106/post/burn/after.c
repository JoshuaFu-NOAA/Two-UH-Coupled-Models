#define VERSION 3.71
#define VERDATE "3-Jan-2002"

/* ============================================================= */
/*                                                               */
/* postprocessing program for ECHAM data and ECMWF analysis data */
/*                                                               */
/* Luis     Kornblueh   - MPI    Hamburg                         */
/* Uwe      Schulzweida - MPI    Hamburg                         */
/* Arno     Hellbach    - DKRZ   Hamburg                         */
/* Edilbert Kirk        - MI Uni Hamburg                         */
/* Michael  Ponater     - DLR    Oberpfaffenhofen                */
/*                                                               */
/* ============================================================= */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* REVISION HISTORY                                              */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 2.11 (Arno - 25-Aug-94)                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* -> included processing of GRIB Version 1                      */
/*    neccessary changes:                                        */
/*    - defined global variable GribEdition                      */
/*    - changed readgrib to read full GRIB message either 1 or 0 */
/*    - changed decogb in order to skip 4 bytes after            */
/*      GRIB block 0 and before GRIB block 2 if GribEdition == 1 */
/*                                                               */
/* -> included timesel to hold up to 4 times                     */
/*    neccessary changes:                                        */
/*    - removed global variable timesel,                         */
/*      timesel as input parameter name is still valid !!        */
/*    - defined global constant MAX_HOURS ( = 4)                 */
/*    - defined global field hours[MAX_HOURS]                    */
/*    - defined global variable nrqh as a counter                */
/*    - added routine scantime (adapted from scanlevel)          */
/*    - changed AnalysisProcess and EchamProcess to enable       */
/*      selection of more than one time                          */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 2.12 (Edilbert - 29-Aug-94)                           */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* -> humidity is now optional in model data sets:               */
/*    MakeGeopotHeight() checks if humidity is present and uses  */
/*    either a wet formula (with virtual temperature) or a dry   */
/*    one. Missing humidity in model data is not a fatal error.  */
/*                                                               */
/* -> performance optimization:                                  */
/*    Functions PutnByte() replaced my macros                    */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 2.2 (Edilbert - 04-Jul-95)                            */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* -> The FORTRAN FFT-subroutines are replaced by C-functions    */
/*    Now this is a pure C-program - the new FFT routines are    */
/*    slightly faster on the CRAY than the old ones.             */
/*                                                               */
/* -> The afterburner can now process grib files from models     */
/*    which simulate the atmosphere of planet Mars               */
/*    namelist parameter: mars=1                                 */
/*                                                               */
/* -> The Free function now returns a NULL which can be used     */
/*    to clear the pointer, that was freed.                      */
/*                                                               */
/* -> The number of latitudes is calculated from truncation      */
/*    with Latitudes = 2 * ((Truncation * 3 + 3) / 4)            */
/*    to allow all resolutions                                   */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 2.21 (Arno - 24-Aug-95)                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> used model ID (byte 5 of block 1) to determine wether NCAR */
/*    based atmospheric model is used (retrieved in precntl)     */
/*    model ID's known at DKRZ (so far):                         */
/*                2 HOPE Ocean Model                             */
/*                3 OPYC Ocean Model                             */
/*                4  LSG Ocean Model                             */
/*               21 NCAR based AGCM                              */
/*               31 ECHAM 1/2/3                                  */
/*               50 ECHAM 4                                      */
/*              121 ECMWF Analysis                               */
/*                                                               */
/* -> provided code to decode complex packed spherical harmonics */
/*    changed decogb (added one if branch to be able to get info */
/*    variable PackComplexStart and to get PackComplexScale,     */
/*    indicator is variable PackComplex = TRUE/FALSE)            */
/*    added functions ScaleComplex and ScatterComplex called by  */
/*    *Control functions.                                        */
/*                                                               */
/* -> added code to calculate layer integrals of cloud cover and */
/*    liquid water content (LayerCloud and LayerWater).          */
/*    New codes:                                                 */
/*               34 integral of low cloud cover 750 - 1013 hPa   */
/*               35 integral of mid cloud cover 460 -  730 hPa   */
/*               36 integral of hih cloud cover  50 -  440 hPa   */
/*               37 integral of low LWC         750 - 1013 hPa   */
/*               38 integral of mid LWC         460 -  730 hPa   */
/*               39 integral of hih LWC          50 -  440 hPa   */
/*               40 integral of all LWC          50 - 1013 hPa   */
/*                                                               */
/* -> added function statistics                                  */
/*    statistics prints time and memory counts performed in main */
/*    by calls to clock() and in the *Alloc functions, by calls  */
/*    to mallinfo() to get the number of used blocks             */
/*                                                               */
/* -> changed data type of level[] to double                     */
/*    neccessary to process level below 1 Pa, i.e. for ECHAM/L43 */
/*    or more.                                                   */
/*                                                               */
/* -> changed all occurences of half_press to HalfPress->hybrid  */
/*    since both were the same in the old versions (minimizes    */
/*    memory allocation)                                         */
/*                                                               */
/* -> performance optimization:                                  */
/*                                                               */
/*    included #pragma _CRI inline statements for ExtraT/ExtraZ  */
/*    to force inlining of those functions in Interpolate_T/Z.   */
/*                                                               */
/*    rewrote encoding according to ECHAM4 procedure (E.Kirk).   */
/*    provided code for machines which do not have _pack routine */
/*                                                               */
/*    rewrote decoding according to encoding (used _unpack)      */
/*                                                               */
/*    changed legini in order to read polinomial coefficients    */
/*    instead of calculating them (done for T21/30/42/63/106).   */
/*    if FILE *legpol == NULL, calculation is done as usual.     */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 2.4 (Edilbert - 29-Nov-95)                            */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* -> on SUN and IBM the non-GRIB output file is written in      */
/*    FORTRAN binary form (with controlwords just before and     */
/*    after the data records. (#define FORTRAN_OUTPUT)           */
/*                                                               */
/* -> After ECHAM-GRIB and ECMWF-GRIB a third data format is     */
/*    recognized in the input file:                              */
/*    PUMA (Portable University Model of the Atmosphere) files   */
/*    are processed like ECHAM files                             */
/*                                                               */
/* -> Introduced command line options                            */
/*    -c for list of codes                                       */
/*    -d for additional debug output (previous Debug parameter)  */
/*    -h help (listing of command line syntax)                   */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 2.5 (Arno - 16-Jul-96)                                */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> Changed algorithm for saturation vapour pressure in sh2rh  */
/*    and rh2sh according to ECMWF IFS code, since the reanalysis*/
/*    project has been performed using this code.                */
/*                                                               */
/* -> The reanalysis data comes with 2 new packing types:        */
/*    - complex packed spherical harmonics                       */
/*    - data on quasi-regular (reduced) grids                    */
/*                                                               */
/*    In order to come up with these function decogb was modified*/
/*    mainly in section 2 and 4. 2 global flags were introduced: */
/*    PackComplex and ReducedGrid. There are 3 new functions     */
/*    named ScatterComplex, ScaleComplex and qu2reg2_, the latter*/
/*    one in a separate source file.                             */
/*                                                               */
/* -> Introduced command line option                             */
/*    -s to reduce memory requirements.                          */
/*       Method: swap out (and Free) 3d grid fields and          */
/*               swap them in (and Alloc).                       */
/*                                                               */
/*    3 new functions are invoked: swapini, swapout and swapin.  */
/*    swapini scans environment for TMPDIR and constructs        */
/*            swapname[codeno]                                   */
/*    swapout opens swapfile[codeno], writes 3d grid field and   */
/*            closes swapfile[codeno]. Always used in conjunction*/
/*            with Free.                                         */
/*    swapin  opens swapfile[codeno], reads  3d grid field and   */
/*            closes swapfile[codeno]. Always used in conjunction*/
/*            with DoubleAlloc.                                  */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 2.6 (Arno - 16-Feb-97)                                */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> Changed GRIB read algorithm in order to minimize the fread */
/*    library calls. The idea was to read blocks as large as the */
/*    environment allows. Therefore the token MAX_RBUF_LENGTH is */
/*    defined. The blocks (sizes of multiples of timesteps) are  */
/*    then scanned for "GRIB" strings and a pointer is returned  */
/*    to the calling routine. This minimizes system CPU time     */
/*    extraordinarily.                                           */
/*                                                               */
/* -> In order to make the afterburner pipe-able the "rewind"    */
/*    had to be erased. Due to this the PUMA processing couldn't */
/*    be hold longer inside the afterburner. I apologize.        */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 2.64 (Arno - 22-Jul-97)                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> Invited "level error-handling"                             */
/*    Some ugly combinations of double typed and non existent    */
/*    levels caused dependent routines like rh2sh() to crash,    */
/*    because 2D-slices were used as initialised (that is, 0).   */
/*    The number and location of "levels" in GRIB-input is       */
/*    determined in function precntl(). In scanlevel(), this     */
/*    list is compared with the "requested levels" and the       */
/*    "requested levels" themselves are checked for redundancy.  */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 2.70 (Arno - 15-Oct-97)                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> saturation water vapour pressure formula adopted from IFS  */
/*    (the Integrated Forecast System, ECMWF) was changed in a   */
/*    way that the ice phase is now neglected and even below 0 C */
/*    the pressure is calculated assuming conditions over water. */
/*                                                               */
/* -> bufferlength RBufLength can now be passed as a line argu-  */
/*    ment with "-b". Per default RBufLength = 13 MB which is    */
/*    app. the length of a single-date T106-output. This can be  */
/*    reduced (and therefore memory) if coarser resolutions are  */
/*    to be processed.                                           */
/*                                                               */
/* -> PUMA is in again and also pipe-able now.                   */
/*                                                               */
/* -> An inventory is made in the pre-processing stage. One time-*/
/*    step is scanned and a list of codes available is made.     */
/*    The list is used to determine if any code "selected" shall */
/*    be calculated as a derived code or just decoded.           */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 2.80 (Edilbert - 21-Oct-97)                           */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> The PUMA file format was changed to IEEE                   */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 3.00 (Uwe - 12-Jan-99)                                */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> all INT to int                                             */
/*    all REAL to double                                         */
/*    lgrib to unsigned int                                      */
/*    FileLength,SaveRBufLength and RBufLength to size_t         */
/*    convert to ANSI C                                          */
/*                                                               */
/* -> runs correct with 32 bit integer now                       */
/*      (fix bug with ibits == 32)                               */
/*                                                               */
/* -> performance optimization:                                  */
/*                                                               */
/*      included #pragma loop scalar for vpp7                    */
/*      included #pragma vdir for nec                            */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 3.01 (Luis - 09-Feb-99)                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> replaced gauaw, which failed to converge for truncations   */
/*    larger than T213, and added float.h.                       */
/*                                                               */
/* -> There seems to be some reorgnising of the code necessary.  */
/*    Up to now the whole bunch of data is loaded into memory.   */
/*    For the newest ECMWF analysis this is getting to big ...   */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 3.02 (Uwe - 12-Feb-99)                                */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> reorganize the Buffer handling                             */
/*    (main, AnalysisControl, EchamControl, precntl)             */
/*                                                               */
/* -> performance optimization for NEC by Eckhard Tschirschnitz  */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 3.03 (Uwe - 18-Feb-99)                                */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> use mallinfo if malloc.h is included                       */
/*    use mallinfo for all memory allocation                     */
/* -> FORTRAN_OUTPUT if CRAY is undefinded                       */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 3.04 (Andreas Rhodin 22-Feb-99)                       */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> introduced parameter LAT                                   */
/*    purpose: specify number of latitudes in gridpoint space    */
/*             independently from spectral truncation.           */
/*    works only if input is spherical harmonics,                */
/*    otherwise aborts with 'Latitude/Longitude Conflict'.       */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 3.05 (Uwe - 23-Feb-99)                                */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> MakeGeopotHeight: split loops for save vectorization       */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 3.06 (Uwe - 03-Mar-99)                                */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> new function gp2sp                                         */
/* -> use gp2sp if spectral Divergence,Vorticity or LnPs missing */
/* -> use MaximalProductLength for RBuf                          */
/* -> scancode: code=-1 means all detected codes selected        */
/* -> scanlevel: level=-1 means all detected levels selected     */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 3.07 (Uwe - 09-Mar-99)                                */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> GribFromBuffer: use ItemsInBufr for RBuf                   */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 3.08 (Uwe - 18-Mar-99)                                */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> use RBufCount in CheckAnalyses and CheckContent            */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 3.09 (Uwe+Luis - 19-Apr-99)                           */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> output GRIB version 1                                      */
/* -> set MIN_RBUF_LENGTH to 5MB                                 */
/* -> add leveltype and layer to struct Control                  */
/*    use layer for GRIB output for leveltype LEV_DOWN and       */
/*    LEV_HEIGHT                                                 */
/* -> use TimeU, Time1, T_R_I and Average for GRIB output        */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 3.10 (Uwe - 16-Sep-99)                                */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> Calculate MaximalProductLength after read all codes        */
/* -> option -a    : forces analysis data process                */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 3.11 (Uwe - 18-Feb-2000)                              */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> change (nh > Dim3GP) to (nh >= Dim3GP)  (Edi Kirk)         */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 3.20 (Uwe - 30-Mar-2000)                              */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> remove correction of incomming spectral fields             */
/* -> add 0.5 to the packed GRIB outputfile to get EMOS          */
/*    compatible results                                         */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 3.21 (Uwe - 22-May-2000)                              */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> use HHMM instead of only HH in time record for service     */
/*        output                                                 */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 3.22 (Uwe - 21-Feb-2001)                              */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> T341 support                                               */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 3.40 (Uwe - 21-April-2001)                            */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> split after.c in:                                          */
/*      afterprog.c alloc.c fft.c     interpolation.c  physc.c   */
/*      convert.c   grib.c  service.c util.c           puma.c    */
/* -> Ansi C conform source files:                               */
/*      remove:   RegularFile, FIFO, SymbolicLink                */
/* -> Swap does not work with Mean:                              */
/*      exit if Swap and Mean is defined                         */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 3.50 (Uwe - 22-April-2001)                            */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> remove PUMA format                                         */
/* -> remove buffering of GRIB reading:                          */
/*    remove: GribFromBuffer, SyncBuffer                         */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 3.51 (Uwe - 02-May-2001)                              */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> read news from /pf/m/m214003/doc/afterburner.doc           */
/* -> write stat to  /mf/m/m214003/usage/after                   */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 3.52 (Uwe - 04-May-2001)                              */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> set buffersize of inputfile to FILEBUFFERSIZE              */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 3.53 (Uwe - 31-May-2001)                              */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> Add new Namelist parameter INTERVAL to change the output   */
/*    interval between monthly and daily                         */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 3.54 (Uwe - 5-Jun-2001)                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> compute GeopotHeight only if it is not deteted             */
/* -> add 1500 byte to igrib buffer in codegb                    */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 3.55 (Uwe - 10-Jul-2001)                              */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> set FILEBUFFERSIZE to  16384 bytes                         */
/* -> set INTERVAL to UNLIM_INTERVAL to process all timesteps    */
/* -> option -w to write warnings instead of program termination */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 3.56 (Uwe - 31-Oct-2001)                              */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> T85 and T31 support                                        */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 3.60 (Uwe - 16-Nov-2001)                              */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> new selection parameter >layer< to select multilevel codes */
/*    which are not on hybrid or pressure level                  */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 3.70 (Uwe - 18-Dec-2001)                              */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> change GRIB IO to gribex function layer                    */
/* -> change file buffer size to filesystem blocksize            */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Version 3.71 (Uwe - 3-Jan-2001)                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*                                                               */
/* -> clean up                                                   */
/* -> fix problem with dependencies for u and v (call to scaluv )*/
/*    for presure level data type 70                             */
/*                                                               */
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef _AFTER_H
#define _AFTER_H

#if defined (__uxp__)
#pragma global noalias
#pragma global novrec
#endif
#if defined (SX)
#pragma odir switch,-dv
#endif

/* =============================================== */
/* These include files should be standard on all   */
/* UNIX systems.                                   */
/* =============================================== */
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <pwd.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#if (! defined sun) && (! defined __uxp__)
#include <malloc.h>
#endif

#define MINFILEBUFFERSIZE 20480
char *filebuffer;

#define FT_GRIB  1
#define FT_SERV  2
#define FT_CDF   3

#ifndef CHAR
#define CHAR unsigned char
#endif

#ifndef TRUE
#define TRUE  1
#endif
#ifndef FALSE
#define FALSE 0
#endif

#ifndef M_PI
#define M_PI        3.14159265358979323846
#endif
#ifndef M_SQRT
#define M_SQRT2     1.41421356237309504880
#endif

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

#define MAX_Longitudes 960
#define MAX_Latitudes  480
#define MAX_DimGP     (MAX_Longitudes * MAX_Latitudes)

#define CODES 276

#define GEOSCODE 129
#define    TCODE 130
#define    UCODE 131
#define    VCODE 132
#define   SHCODE 133
#define   PSCODE 134
#define    WCODE 135
#define  VORCODE 138
#define  STRCODE 148
#define  VELCODE 149
#define  SLPCODE 151
#define LNPSCODE 152
#define  DIVCODE 155
#define    ZCODE 156
#define   RHCODE 157
#define  LBCCODE 190
#define  SBCCODE 189
#define  STCCODE 191
#define  LTCCODE 192
#define  NTCCODE 193
#define   LCCODE  34
#define   MCCODE  35
#define   HCCODE  36
#define   LWCODE  37
#define   MWCODE  38
#define   HWCODE  39
#define   AWCODE  40


struct Date
{
   int ct;
   int yr;
   int mo;
   int dy;
   int hr;
   int mn;
   int YYMMDD;
   int HHMM;
};

/* Memory Allocation */

char     *CharAlloc(size_t bytes, char *array_name);
int       *IntAlloc(int words,    char *array_name);
float   *FloatAlloc(int words,    char *array_name);
double *DoubleAlloc(int words,    char *array_name);
void       *FreeMem(void *ptr);

/* UTIL */

void Stars(int n);
void NewLine(void);
void Abort(char *errtext);
void RealCopy(void *destination, void *source, int words);
void IntZero(int *field, int words);
void RealZero(double *field, int words);

/* FFT */

void fft_set(int ifax[], int n);
int rpassc(double *a, double *b, double *c, double *d,
           int inc1, int inc2, int inc3, int inc4,
           int lot , int n   , int ifac, int la  );
int qpassc(double *a, double *b, double *c, double *d,
           int inc1, int inc2, int inc3, int inc4,
           int lot , int n   , int ifac, int la  );
void fc2gp(double *fc, double *gp, int Lat, int Lon, int Lev, int Fou);
void gp2fc(double *gp, double *fc, int Lat, int Lon, int Lev, int Fou);

/* GRIB */

void GribAbort(int errno);
int Get2Byte(CHAR ptr[]);
void codegb (double *data,     /* Pointer to input data */
             int     datalen,  /* length of data record */
             char   *CGrib,    /* packed GRIB record */
             int    *griblen,  /* length of GRIB record */
             int Code, int Level);

void ScaleComplex (double *fpdata, int Truncation);
void ScatterComplex(double *fpdata, int Truncation, int DimSP);

/* SERVICE */

void wrifor(double *field, int dim, int head[8], FILE* gp);

/* Convert */

/* Convert Spectral Array to new resolution */
void ChangeTruncation(double *SPin, int Resin, double *SPout, int Resout);
void sp2fc(double *sa, double *fa, double *poli, int klev, int nlat, int nfc,int nt);
void fc2sp(double *fa, double *sa, double *poli, int klev, int nlat, int nt);

/* Interpolation */
void
geninz (int *interpolation_index, double *pres_of_height,
        double *full_level_pressure, int dimgp, int nrql, int levels);
void
genind (int *interpolation_index, double lv[],
        double *full_level_pressure, int dimgp, int nrql, int levels);
void
Interpolate_X (double *gt, double *pt, double *pf, int *nx,
               double level[], int nrql, int dimgp, int dim3gp);
void
Interpolate_T (double orography[], double gt[], double pt[], double pf[], double ph[], int *nx,
               double level[], int levels, int nrql, int dimgp);
void
Interpolate_Z (double orography[], double gz[], double pz[], double pf[], double ph[], int *nx, double gt[],
               double level[], int levels, int nrql, int dimgp);
void
Interpolate_X_Z (double *gt, double *pt, double *pf, int *nx, double *poh,
                 int nrql, int dimgp, int dim3gp);
void
Interpolate_T_Z (double orography[], double gt[], double pt[], double pf[], double ph[],
                 int *nx, double *poh, int levels, int nrql, int dimgp);
void
Interpolate_Z_Z (double orography[], double gz[], double pz[], double pf[], double ph[],
                 int *nx, double *poh, double gt[], int levels, int nrql, int dimgp);
void Extrap(double *slp, double *aph, double *apf,
            double *Geopotential, double *t, int nhor);
double ExtraT (double PRES, double APH, double APF, double GEOS, double T);
double ExtraZ (double pres, double aph, double apf, double Geopotential, double t);

/* Physc */

void IniQuaSum(double Destination[], double Source[], int Length);
void AddQuaSum(double Destination[], double Source[], int Length);
void VarQuaSum(double Variance[], double Sum[], int Length, int n);
void AddVector(double Destination[], double Source[], int Length);
void MultVectorScalar(double Destination[], double Source[],
                      double Factor, int Length);
void Add2Vectors(double *Destination, double *SourceA,
                 double *SourceB, int Length);
void Sub2Vectors(double *Destination, double *SourceA,
                 double *SourceB, int Length);
void Speed(double *speed, double *u, double *v, int dim3gp);
void dv2ps(double *div, double *pot, int lev, int trunc);
void dv2uv(double *d, double *o, double *u, double *v, double *f, double *g,
           int nt, int nsp, int nlev);
void scaluv(double *fu, double rclat[], int nlat, int lot);
void uv2dv(double *fu, double *fv, double *sd, double *sv,
           double *pol2, double *pol3, int klev, int nlat, int nt);
void geninx(int nt, double *f, double *g);
void theta(double *PThetaF, double *PThetaH, double *PH, double *PS,
           double *TF, double *TS, int Levels, int DimGP, int Dim3GP);
void presh(double *pf, double *php, double *vct, double *ps,
           int Levels, int DimGP, int Dim3GP);
void gauaw (double pa[], double pw[], int nlat);
void phcs(double *PNM, double *HNM, int Waves, double PMU,
          double *ZTEMP1, double *ZTEMP2);
void sh2rh(double *sphum, double *rhum, double *t, int lev,
           int dimgpout, double *level, double *fullpresshybrid);
void rh2sh(double *sphum, double *rhum, double *t, int lev,
           int dimgpout, double *level);
void MakeGeopotHeight(double *geop, double* gt, double *gq, double *ph, int nhor, int nlev);
void LayerWater (double *ww, double *ll, double pmax, double pmin,
                 int DimGP, int HalfLevels, double *vct);
void LayerCloud (double *cc, double *ll, double pmax, double pmin,
                 int DimGP, int HalfLevels, double *vct);
void Derivate(double field[], double derilam[], int levels,
              int Waves, int Latitudes, double DerivationFactor[]);
void h2p(double *PHeight, double lv[], int DimGP, int nrql);


void gribexdp(int *sec0, int *sec1, int *sec2, double *psec2, int *sec3,
              double *psec3, int *sec4, double *psec4, int klenp, int *kgrib,
              int kleng, int *kword, char *hoper, int *kret);

#endif /*  after.h  */
#ifndef _ERROR_H
#define _ERROR_H

#define  _FATAL     1     /* Error flag: exit on error  */
#define  _VERBOSE   2     /* Error flag: report errors  */
#define  _DEBUG     4     /* Error flag: debug          */

extern int _ExitOnError;  /* If set to 1, exit on error (default 1)       */
extern int _Verbose;      /* If set to 1, errors are reported (default 1) */
extern int _Debug;        /* If set to 1, debuggig (default 0)            */

void SysError(const char *caller, const char *fmt, ...);
void    Error(const char *caller, const char *fmt, ...);
void  Warning(const char *caller, const char *fmt, ...);
void  Message(const char *caller, const char *fmt, ...);

#endif  /* _ERROR_H */
typedef struct
{
  int code;
  int level;
  int date;
  int time;
  int nlon;
  int nlat;
  int dispo1;
  int dispo2;
}
SRVHEAD;
#include <stdio.h>
#include <stdarg.h>

int gribFileSeek (FILE *fp, int *recpos);
int gribReadSize (FILE *fp);
int grb_read_record (FILE *fp, void *record, size_t recordsize);
void gribWarning (const char *func, const char *fmt, ...);
#define LEV_SURFACE      1
#define LEV_99          99
#define LEV_ISOBARIC   100
#define LEV_MEANSEA    102
#define LEV_ALTITUDE   103
#define LEV_HEIGHT     105
#define LEV_SIGMA      107
#define LEV_HYBRID     109
#define LEV_DOWN       111
#define LEV_GROUND     112
#define LEV99_MARGIN  1000

#define REP_REGULAR      0
#define REP_GAUSS        4
#define REP_SPECTRAL    50


/*
 *  Macros for the indicator section ( Section 0 )
 */

#define GribLength(sec0)               (sec0[ 0])  /* Number of octets in the GRIB message        */
#define GribEdition(sec0)              (sec0[ 1])  /* GRIB edition number                         */


/*
 *  Macros for the product definition section ( Section 1 )
 */

#define GribCodeTable(sec1)            (sec1[ 0])  /*  Version number of code table               */
#define GribCenterID(sec1)             (sec1[ 1])  /*  Identification of centre                   */
#define GribModelID(sec1)              (sec1[ 2])  /*  Identification of model                    */
#define GribGridDefinition(sec1)       (sec1[ 3])  /*  Grid definition                            */
#define GribBlock2Included(sec1)       (sec1[ 4])  /*  Block 2 included                           */
#define GribParameter(sec1)            (sec1[ 5])  /*  Parameter indicator                        */
#define GribLevelType(sec1)            (sec1[ 6])  /*  Type of level indicator                    */
#define GribLevel1(sec1)               (sec1[ 7])  /*  Level 1                                    */
#define GribLevel2(sec1)               (sec1[ 8])  /*  Level 2                                    */
#define GribYear(sec1)                 (sec1[ 9])  /*  Year of century (YY)                       */
#define GribMonth(sec1)                (sec1[10])  /*  Month (MM)                                 */
#define GribDay(sec1)                  (sec1[11])  /*  Day (DD)                                   */
#define GribHour(sec1)                 (sec1[12])  /*  Hour (HH)                                  */
#define GribMinute(sec1)               (sec1[13])  /*  Minute (MM)                                */
#define GribTimeUnit(sec1)             (sec1[14])  /*  Time unit indicator                        */
#define GribTimePeriod1(sec1)          (sec1[15])  /*  P1 Time period                             */
#define GribTimePeriod2(sec1)          (sec1[16])  /*  P2 Time period                             */
#define GribTimeRange(sec1)            (sec1[17])  /*  Time range indicator                       */
#define GribNumAvg(sec1)               (sec1[18])  /*  Number of products included in an average  */
#define GribNumMiss(sec1)              (sec1[19])  /*  Number of products missing form an average */
#define GribCentury(sec1)              (sec1[20])  /*  Century                                    */
#define GribSubCenterID(sec1)          (sec1[21])  /*  Subcenter identifier                       */
#define GribDecScaleFactor(sec1)       (sec1[22])  /*  Decimal scale factor                       */

#define GribECMWFLocalExtention(sec1)  (sec1[36])
#define GribECMWFClass(sec1)           (sec1[37])


/*
 *  Macros for the grid definition section ( Section 2 )
 */

#define GribGridType(sec2)             (sec2[ 0])  /* Data representation type */

/* Spherical harmonic coeficients */

#define GribPentaJ(sec2)               (sec2[ 1])  /* J pentagonal resolution parameter             */
#define GribPentaK(sec2)               (sec2[ 2])  /* K pentagonal resolution parameter             */
#define GribPentaM(sec2)               (sec2[ 3])  /* M pentagonal resolution parameter             */
#define GribRepType(sec2)              (sec2[ 4])  /* Representation type                           */
#define GribRepMode(sec2)              (sec2[ 5])  /* Representation mode                           */

/* Gaussian grids */

#define GribNumLon(sec2)               (sec2[ 1])  /* Number of points along a parallel (Ni)        */
#define GribNumLat(sec2)               (sec2[ 2])  /* Number of points along a meridian (Nj)        */
#define GribFirstLat(sec2)             (sec2[ 3])  /* Latitude of the first grid point              */
#define GribFirstLon(sec2)             (sec2[ 4])  /* Longitude of the first grid point             */
#define GribResFlag(sec2)              (sec2[ 5])  /* Resolution flag: 128 regular grid             */
#define GribLastLat(sec2)              (sec2[ 6])  /* Latitude of the last grid point               */
#define GribLastLon(sec2)              (sec2[ 7])  /* Longitude of the last grid point              */
#define GribLonIncr(sec2)              (sec2[ 8])  /* i direction increment                         */
#define GribNumPar(sec2)               (sec2[ 9])  /* Number of parallels between a pole and the E. */

#define GribNumVCP(sec2)               (sec2[11])  /* Number of vertical coordinate parameters      */

/*
 *  Macros for the binary data section ( Section 4 )
 */

#define GribNumValues(sec4)            (sec4[ 0])  /* Number of data values for encode/decode     */
#define GribNumBits(sec4)              (sec4[ 1])  /* Number of bits used for each encoded value  */



void gribex  (int *ksec0, int *ksec1, int *ksec2, void *psec2, int *ksec3,
              void *psec3, int *ksec4, void *psec4, int klenp, int *kgrib,
              int kleng, int *kword, char *hoper, int *kret);

void gribexsp(int *ksec0, int *ksec1, int *ksec2, float *psec2, int *ksec3,
              float *psec3, int *ksec4, float *psec4, int klenp, int *kgrib,
              int kleng, int *kword, char *hoper, int *kret);

void gribexdp(int *ksec0, int *ksec1, int *ksec2, double *psec2, int *ksec3,
              double *psec3, int *ksec4, double *psec4, int klenp, int *kgrib,
              int kleng, int *kword, char *hoper, int *kret);

void grprs0  (int *ksec0);

void grprs1  (int *ksec0, int *ksec1);

void grprs2dp(int *ksec0, int *ksec2, double *psec2_dp);
void grprs2sp(int *ksec0, int *ksec2, float  *psec2_sp);
void grprs2  (int *ksec0, int *ksec2, void   *psec2);

void grprs3dp(int *ksec0, int *ksec3, double *psec3_dp);
void grprs3sp(int *ksec0, int *ksec3, float  *psec3_sp);
void grprs3  (int *ksec0, int *ksec3, void   *psec3);

void grprs4dp(int *ksec0, int *ksec4, double *psec4_dp);
void grprs4sp(int *ksec0, int *ksec4, float  *psec4_sp);
void grprs4  (int *ksec0, int *ksec4, void   *psec4);

void grprs4w (int *ksec4);

int _ExitOnError   = 1;	/* If set to 1, exit on error       */
int _Verbose = 1;	/* If set to 1, errors are reported */
int _Debug   = 0;       /* If set to 1, debugging           */


void
SysError(const char *caller, const char *fmt, ...)
{
  va_list args;
	
  va_start(args, fmt);

   fprintf(stderr, "Error (%s) : ", caller);
  vfprintf(stderr, fmt, args);
   fprintf(stderr, "\n");

  va_end(args);

  perror("system error message ");
	
  exit(1);
}

void
Error(const char *caller, const char *fmt, ...)
{
  va_list args;
	
  va_start(args, fmt);

   fprintf(stderr, "Error (%s) : ", caller);
  vfprintf(stderr, fmt, args);
   fprintf(stderr, "\n");

  va_end(args);

  if ( _ExitOnError ) exit(1);
}

void
Warning(const char *caller, const char *fmt, ...)
{
  va_list args;
	
  va_start(args, fmt);

  if ( _Verbose )
    {
       fprintf(stderr, "Warning (%s) : ", caller);
      vfprintf(stderr, fmt, args);
       fprintf(stderr, "\n");
    }

  va_end(args);
}

void
Message(const char *caller, const char *fmt, ...)
{
  va_list args;
	
  va_start(args, fmt);

   fprintf(stderr, "%-18s : ", caller);
  vfprintf(stderr, fmt, args);
   fprintf(stderr, "\n");

  va_end(args);
}

#ifndef CRAY
#define FORTRAN_OUTPUT
#endif

/* GRIBEX */
int    kleng, klenp, kword, kret;
int    sec0[2], sec1[1024], sec2[1024], sec3[2], sec4[512];
double psec2[512], psec3[2];


void Usage(void);

#ifndef STDOUT
#define STDOUT stderr
#endif

#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC 1000000.
#endif

#define MAX_HOURS        4
#define MAX_LEVELS      99
#define MEGABYTE 1024*1024

#define L_TIMES_RHOH2O (-333700000.0)
#define EARTH_RADIUS 6371000.0
#define MARS_RADIUS  3400000.0
#define EARTH_GRAV 9.80665
#define MARS_GRAV  3.7
#define RG   (1.0 / Grav)
#define MARS_RD     (189.0 )

/* ************************************** */
/* Thermodynamical constants adopted from */
/* ECMWF IFS-Code                         */
/* ************************************** */

#define RKBOL (1.380658e-23)
#define RNAVO (6.0221367e+23)
#define R     (RKBOL * RNAVO)
#define RMD   (28.9644)
#define RMV   (18.0153)
#define EARTH_RD    (1000. * R / RMD)
#define RV    (1000. * R / RMV)

int  fileTypeIn  = 0;
int  fileTypeOut = 0;

char   *filename;
char    ifile[128];
char    ofile[128];
char    namelist[1024];
char    AllocName[256];

int    Zero = 0, One = 1, Two = 2;
/* -------------------------------- */
/* Variables related to GRIB record */
/* -------------------------------- */

/*       int    z             ; *//* Counter of GRIB length for output          */
int    Grib1Offset   ; /* Offset in bytes before and after block 1   */
int    CodeGrib      ; /* grib[12]   Block1[ 8]   1 Byte only        */
int    TypeOfLevel   ; /* grib[13]   Block1[ 9]   1 Byte             */
int    LevelGrib     ; /* grib[14]   Block1[10]   1 or 2 Bytes       */
int    LayerGrib     ; /* grib[14]   Block1[10]   1 or 2 Bytes       */

int    RepGrib       ; /* grib[33]   Block2[ 5]   1 Byte             */

int    AnalysisData; /* 0 = ECHAM Data, 1 = ECMWF Spectral Analyses */
int    DayIn       ; /* day increment of infiles if Multi = TRUE    */
int    Debug;

int    Dim3FC     ;
int    Dim3SP     ;
int    Dim3GP     ;
int    DimFC      ;
int    DimGP      ;
int    DimSP      ;
int    DimSP_half ;
int    Dim3FCOut     ;
int    Dim3SPOut     ;
int    Dim3GPOut     ;
int    DimFCOut      ;
int    DimGPOut      ;
int    DimSPOut      ;
int    DimSPOut_half ;
int    ECMWFMars = FALSE;

#define TIMESTEP_INTERVAL  -1
#define MONTHLY_INTERVAL    0
#define DAILY_INTERVAL      1
#define UNLIM_INTERVAL      2

int    OutputInterval;
int    EndOfInterval ;

int    Fouriers      ;
int    FouriersOut   ;
int    Gaussian = FALSE;
double Grav          = EARTH_GRAV;
int    Grib          ;
int    GribErr       ;
int    HalfLevels    ;
int    LatitudesIn   = 0;
int    Latitudes     = 0;
int    LatitudesOut  ;
int    LatitudeFirst ;
int    LatitudeLast  ;
int    LongitudeFirst;
int    LongitudeLast ;
int    LongitudeIncr ;
int    Layers = 0    ;
int    Levels = 0    ;
int    LevelType     ;
int    AnaLevelFactor = 1;
int    Longitudes    ;
int    LongitudesOut ;
double L_times_rhoH2O = L_TIMES_RHOH2O;
int    mars          ;

int    Mean          ;
int    MeanCount     ;
int    Multi         ;
double PlanetRadius = EARTH_RADIUS;
double RD           = EARTH_RD;
int    Representation;
int    Spectral = FALSE;
int    Swap = FALSE  ;
int    TermCount     ;
int    Truncation = 0;
int    TruncationOut ;
int    TruncationGrib;
int    Waves         ;
int    WavesOut      ;

int    labort = TRUE;

struct Date NewDate;
struct Date OldDate;
struct Date StartDate;

struct Control
{
   int   needed    ;
   int   selected  ;
   int   detected  ;
   int   swapped   ;
   int   sfit      ;
   int   surf      ;
   int   hlev      ;
   int   plev      ;
   int   nlay      ;
   int   layers[MAX_LEVELS] ;
   int   leveltype ;
   int   layer     ;
   double *spectral  ;
   double *fourier   ;
   double *hybrid    ;
   double *height    ;
   double *grid      ;
   double *mean      ;
   double *var       ;
};


char *CodeName[CODES];

struct Control All[CODES+5];

struct Control *Geopotential = &All[129];
struct Control  *Temperature = &All[130];
struct Control       *u_wind = &All[131];
struct Control       *v_wind = &All[132];
struct Control     *Humidity = &All[133];
struct Control           *Ps = &All[134];
struct Control        *Omega = &All[135];
struct Control    *Vorticity = &All[138];
struct Control           *Ts = &All[139];
struct Control       *stream = &All[148];
struct Control      *velopot = &All[149];
struct Control          *SLP = &All[151];
struct Control         *LnPs = &All[152];
struct Control   *Divergence = &All[155];
struct Control *GeopotHeight = &All[156];
struct Control    *Rhumidity = &All[157];
struct Control   *sw_bot_clf = &All[189];
struct Control   *lw_bot_clf = &All[190];
struct Control   *sw_top_clf = &All[191];
struct Control   *lw_top_clf = &All[192];
struct Control  *net_top_clf = &All[193];
struct Control    *low_cloud = &All[ 34];
struct Control    *mid_cloud = &All[ 35];
struct Control    *hih_cloud = &All[ 36];
struct Control    *low_water = &All[ 37];
struct Control    *mid_water = &All[ 38];
struct Control    *hih_water = &All[ 39];
struct Control    *all_water = &All[ 40];
struct Control        *speed = &All[259];
struct Control       *precip = &All[260];
struct Control      *net_top = &All[261];
struct Control      *net_bot = &All[262];
struct Control     *net_heat = &All[263];
struct Control    *net_water = &All[264];
struct Control       *sw_clf = &All[265];
struct Control       *lw_clf = &All[266];
struct Control      *net_clf = &All[267];
struct Control       *sw_atm = &All[268];
struct Control       *lw_atm = &All[269];
struct Control      *net_atm = &All[270];
struct Control  *surf_runoff = &All[271];
struct Control        *dpsdx = &All[273];
struct Control        *dpsdy = &All[274];
struct Control  *fresh_water = &All[275];

struct Control    *HalfPress = &All[277];
struct Control    *FullPress = &All[278];
struct Control       *ThetaH = &All[279];
struct Control       *ThetaF = &All[280];

int    *vert_index;
double *dv2uv_f1;
double *dv2uv_f2;
double *Orography;
double *poli;
double *poliOut;
double *pold;
double *pdev;
double *pol2;
double *pol3;
double *p_of_height;

int nrql;
int nfol;
int nrqh;
int nglay = 0;
int nrqlay;

int nvct = 0;

double *vct = NULL;

int SrvHead[8];

CHAR *grib;

int type;
int unitsel;
int hours[MAX_HOURS+1];
double layer[MAX_LEVELS+1];
double level[MAX_LEVELS+1];
int  LevelFound[MAX_LEVELS+1];
int  LayerFound[MAX_LEVELS+1];

double rclat[MAX_Latitudes];
double coslat[MAX_Latitudes];
double rclatOut[MAX_Latitudes];
double coslatOut[MAX_Latitudes];
double DerivationFactor[MAX_Latitudes];
double field[MAX_DimGP];

FILE  *ifileID;
FILE  *ofileID;
FILE  *legpol;
FILE  *legpolOut;

char  *swapname[CODES];


char *FieldName(int code, char *Text)
{
   if (CodeName[code])
      sprintf(AllocName,"%s[%3d].%s",CodeName[code],code,Text);
   else
      sprintf(AllocName,"[%3d].%s",code,Text);

   return AllocName;
}

void FreeSpectral(void)
{
   int code;

   for (code = CODES-1; code >= 0; --code) {
      if (All[code].spectral) {
          All[code].spectral = (double *) FreeMem(All[code].spectral);
      }
   }
}

void FreeFourier(void)
{
   int code;

   for (code = 0; code < CODES; code++) {
      if (All[code].fourier) {
          All[code].fourier = (double *) FreeMem(All[code].fourier);
      }
   }
}

void FreeHybrid(void)
{
   int code;

   for (code = 0; code < CODES; code++) {
      if (All[code].hybrid) {
          All[code].hybrid = (double *) FreeMem(All[code].hybrid);
      }
   }
}

void FreeGrid(void)
{
   int code;

   for (code = 0; code < CODES; code++) {
      if (All[code].grid) {
          All[code].grid = (double *) FreeMem(All[code].grid);
      }
   }
}

void PutGrib(FILE *fileID, double Field[], int Code, int Level)
{
  size_t w;
  int griblen;
  char *CGrib;

  if ( Debug )
    fprintf(stderr, "Code: %d, Level: %d , Pointer: %p\n", Code, Level, (void *) Field);

  CGrib = CharAlloc(DimGPOut*7,"CGrib");

  codegb(Field, DimGPOut, CGrib, &griblen, Code, Level);

  w = fwrite (CGrib, sizeof(char), griblen, fileID);

  if (Debug)
    fprintf(stderr, "Items packed %5d  Items written %5d\n", griblen, w);
  /*
    printf("putgrib: %d %d %d\n", DimGPOut, DimGPOut*7, griblen);
    */
  CGrib = FreeMem(CGrib);
}

void legini(void)
{
   int jgl, jm, jn, n;
   int jsp, pdim;
   char *dumm;
   double *gmu, *gwt, *hnm, *pnm, *ZTEMP1, *ZTEMP2, gmusq;

   pdim = DimSP_half * Latitudes;
   poli = DoubleAlloc(pdim,"poli.permanent");
   if (!AnalysisData) {
      if (type >= 50)    pold = DoubleAlloc(pdim,"pold.permanent");
      if (dpsdy->needed) pdev = DoubleAlloc(pdim,"pdev.permanent");
   }
   if ((Divergence->needed || Vorticity->needed ||
           velopot->needed ||    stream->needed ) && type > 20) {
      pol2 = DoubleAlloc(pdim,"pol2.permanent");
      pol3 = DoubleAlloc(pdim,"pol3.permanent");
   }
   gmu  = DoubleAlloc(Latitudes,  "legini.gmu");
   gwt  = DoubleAlloc(Latitudes,  "legini.gwt");
   hnm  = DoubleAlloc(DimSP,"legini.hnm");
   pnm  = DoubleAlloc(DimSP,"legini.pnm");
   ZTEMP1 = DoubleAlloc(Waves<<1,"ZTEMP1");
   ZTEMP2 = DoubleAlloc(Waves<<1,"ZTEMP2");

   gauaw (gmu, gwt, Latitudes);

   jgl = 0;
   LatitudeFirst  = (int) (180.*1000.*asin(gmu[jgl])/M_PI);
   jgl = Latitudes - 1;
   LatitudeLast   = (int) (180.*1000.*asin(gmu[jgl])/M_PI);
   LongitudeIncr  = (360 * 1000) / Longitudes;
   LongitudeFirst = 0;
   LongitudeLast  = (Longitudes - 1) * LongitudeIncr;

   if (Debug) {
      fprintf(STDOUT,"+    legini: Latitudes\n");
      fprintf(STDOUT,"+           ");
      for (jgl = 0; jgl < Latitudes; jgl++) {
          fprintf(STDOUT,"%7.2f",180.*asin(gmu[jgl])/M_PI);
          if ((jgl+1)%8 == 0) fprintf(STDOUT,"\n+           ");
      }
      fprintf(STDOUT,"\n");
   }

   if (AnalysisData && (type == 70) && Gaussian && !Spectral) {
      if (poli) poli = FreeMem(poli);
      if (pol2) pol2 = FreeMem(pol2);
      if (pol3) pol3 = FreeMem(pol3);
      return;
   }

   if (legpol) {
      dumm = CharAlloc(pdim,"dumm.scratch");
      if (poli) (void) fread(poli,sizeof(double),pdim,legpol);
      else  for (n = 0; n < sizeof(double); n++)
                (void) fread(dumm,sizeof(CHAR),pdim,legpol);
      if (pold) (void) fread(pold,sizeof(double),pdim,legpol);
      else for (n = 0; n < sizeof(double); n++)
                (void) fread(dumm,sizeof(CHAR),pdim,legpol);
      if (pdev) (void) fread(pdev,sizeof(double),pdim,legpol);
      else for (n = 0; n < sizeof(double); n++)
                (void) fread(dumm,sizeof(CHAR),pdim,legpol);
      if (pol2) (void) fread(pol2,sizeof(double),pdim,legpol);
      else for (n = 0; n < sizeof(double); n++)
                (void) fread(dumm,sizeof(CHAR),pdim,legpol);
      if (pol3) (void) fread(pol3,sizeof(double),pdim,legpol);
      else for (n = 0; n < sizeof(double); n++)
                (void) fread(dumm,sizeof(CHAR),pdim,legpol);
      dumm   = FreeMem(dumm);
      rewind(legpol);

      for (jgl = 0; jgl < Latitudes; jgl++) {
         gmusq = 1.0 - gmu[jgl] * gmu[jgl];
         coslat[jgl] =  sqrt(gmusq);
         rclat[jgl] = 1.0 / coslat[jgl];
         DerivationFactor[jgl] = rclat[jgl] / PlanetRadius;
      }
   }
   else
     for (jgl = 0; jgl < Latitudes; jgl++) {
        gmusq = 1.0 - gmu[jgl] * gmu[jgl];
        coslat[jgl] =  sqrt(gmusq);
        rclat[jgl] = 1.0 / coslat[jgl];
        DerivationFactor[jgl] = rclat[jgl] / PlanetRadius;
        phcs(pnm,hnm,Waves,gmu[jgl],ZTEMP1,ZTEMP2);
        jsp = jgl;
        for (jm = 0; jm < Waves; jm++) {
           for (jn = 0; jn < Waves - jm; jn++) {
                        poli[jsp] = pnm[jm*Waves+jn] * 2.0;
              if (pold) pold[jsp] = pnm[jm*Waves+jn] * gwt[jgl];
              if (pdev) pdev[jsp] = hnm[jm*Waves+jn] * 2.0 /
                                    (-PlanetRadius   * sqrt(gmusq));
              if (pol2) pol2[jsp] = hnm[jm*Waves+jn] * gwt[jgl] /
                                    (PlanetRadius    * gmusq);
              if (pol3) pol3[jsp] = pnm[jm*Waves+jn] * gwt[jgl] * jm /
                                    (PlanetRadius    * gmusq);
              jsp += Latitudes;
           }
        }
     }
   ZTEMP2 = FreeMem(ZTEMP2);
   ZTEMP1 = FreeMem(ZTEMP1);
   pnm    = FreeMem(pnm);
   hnm    = FreeMem(hnm);
   gwt    = FreeMem(gwt);
   gmu    = FreeMem(gmu);
}


void leginiOut(void)
{
   int jgl,jm,jn;
   int jsp,pdim;
   double *gmu,*gwt,*hnm,*pnm,*ZTEMP1,*ZTEMP2,gmusq;

   if (TruncationOut == Truncation) {
      poliOut = poli;
      RealCopy(rclatOut ,rclat ,MAX_Latitudes);
      RealCopy(coslatOut,coslat,MAX_Latitudes);
      return;
   }

   pdim = DimSPOut_half * LatitudesOut;
   poliOut = DoubleAlloc(pdim,"poliOut.permanent");
   gmu = DoubleAlloc(LatitudesOut,"leginiOut.gmu");
   gwt = DoubleAlloc(LatitudesOut,"leginiOut.gwt");
   hnm = DoubleAlloc(DimSPOut,"leginiOut.hnm");
   pnm = DoubleAlloc(DimSPOut,"leginiOut.pnm");
   ZTEMP1 = DoubleAlloc(WavesOut<<1,"ZTEMP1Out");
   ZTEMP2 = DoubleAlloc(WavesOut<<1,"ZTEMP2Out");
   gauaw(gmu,gwt,LatitudesOut);

   jgl = 0;
   LatitudeFirst  = (int) (180.*1000.*asin(gmu[jgl])/M_PI);
   jgl = LatitudesOut - 1;
   LatitudeLast   = (int) (180.*1000.*asin(gmu[jgl])/M_PI);
   LongitudeIncr  = (360 * 1000) / LongitudesOut;
   LongitudeFirst = 0;
   LongitudeLast  = (LongitudesOut - 1) * LongitudeIncr;

   if (legpolOut) {
      (void) fread(poliOut,sizeof(double),pdim,legpolOut);
      for (jgl = 0; jgl < LatitudesOut; jgl++) {
         gmusq = 1.0 - gmu[jgl] * gmu[jgl];
         coslatOut[jgl] =  sqrt(gmusq);
         rclatOut[jgl] = 1.0 / coslatOut[jgl];
      }
   }
   else
      for (jgl = 0; jgl < LatitudesOut; jgl++) {
         gmusq = 1.0 - gmu[jgl] * gmu[jgl];
         coslatOut[jgl] =  sqrt(gmusq);
         rclatOut[jgl] = 1.0 / coslatOut[jgl];
         phcs(pnm,hnm,WavesOut,gmu[jgl],ZTEMP1,ZTEMP2);
         jsp = jgl;
         for (jm = 0; jm < WavesOut; jm++) {
            for (jn = 0; jn < WavesOut - jm; jn++) {
               poliOut[jsp] = pnm[jm*WavesOut+jn] * 2.0;
               jsp += LatitudesOut;
            }
         }
      }

   ZTEMP2 = FreeMem(ZTEMP2);
   ZTEMP1 = FreeMem(ZTEMP1);
   pnm    = FreeMem(pnm);
   hnm    = FreeMem(hnm);
   gwt    = FreeMem(gwt);
   gmu    = FreeMem(gmu);
}

char *GribFromFile(FILE *fp)
{
  static char func[] = "GribFromFile";
  static int recpos, recsize;
  static int buffersize = 0;
  static char *gribbuffer = NULL;
  int    GribEdition   ; /* grib[ 3]                1 Byte             */
  int ierr;

  ierr = gribFileSeek(fp, &recpos);
  if ( ierr > 0 )   gribWarning(func, "GRIB record not found\n");
  if (feof(fp)) GribErr = 1;
  if (GribErr) return (NULL);
  if (ierr == -1)
    GribErr = 1;
  else if (ierr == 1)
    GribErr = 2;

  if (GribErr) return (NULL);

  recsize = gribReadSize(fp);

  if ( recsize == 0 ) GribErr = 2;

  if ( GribErr ) return (NULL);

  if ( buffersize != recsize )
    {
      buffersize = recsize;
      kleng = buffersize;
      gribbuffer = realloc(gribbuffer, buffersize);
    }

  if (Debug)
    fprintf(STDOUT,"FilePointer %10d recpos %9d recsize %9d\n",
	    (int) ftell(fp),     recpos,    recsize);

  gribbuffer[0] = 'G';
  gribbuffer[1] = 'R';
  gribbuffer[2] = 'I';
  gribbuffer[3] = 'B';

  ierr = grb_read_record(fp, &gribbuffer[4], recsize-4);

  if ( ierr != 0 )
    {
      GribErr = 2;
      return (NULL);
    }

   GribEdition = gribbuffer[7];
   Grib1Offset = GribEdition * 4;
   return (gribbuffer);
}



void gp2sp(int ccode)
{
   int FieldLength = 0;

   if (All[ccode].spectral == NULL) {
      if (All[ccode].hybrid == NULL) {
         fprintf(STDOUT,"%d.hybrid not found\n",ccode);
         exit(99);
      }
      if (All[ccode].fourier == NULL) {
         FieldLength = DimFC * All[ccode].hlev;
         All[ccode].fourier = DoubleAlloc(FieldLength,"gp2sp.fourier");
         gp2fc(All[ccode].hybrid,All[ccode].fourier,
              Latitudes,Longitudes,All[ccode].hlev,Fouriers);
      }
      All[ccode].spectral = DoubleAlloc(Dim3SP,"gp2sp.spectral");
      fc2sp(All[ccode].fourier,All[ccode].spectral,
            pold,All[ccode].hlev,Latitudes,Truncation);    
   }
}

void OMEGA(void)
{
   int i,j;
   double DeltaHybrid;
   double Cterm;
   double Pterm;
   double *omega = Omega->hybrid;
   double *diver = Divergence->hybrid;
   double *halfp = HalfPress->hybrid;
   double *fullp = FullPress->hybrid;
   double *uwind = u_wind->hybrid;
   double *vwind = v_wind->hybrid;

/* Compute half level part of vertical velocity */

   for (i = 0; i < DimGP ; i++) omega[i] = 0.0;
   for (j = 0; j < Levels; j++) {
      DeltaHybrid = vct[Levels+j+2] - vct[Levels+j+1];
#if defined (SX)
#pragma vdir nodep
#endif
#if defined (__uxp__)
#pragma loop novrec
#endif
      for (i = 0; i < DimGP; i++) {
        omega[DimGP] = *omega
                      - *diver * (halfp[DimGP] - *halfp) - DeltaHybrid
                      * (*uwind * dpsdx->hybrid[i]
                      +  *vwind * dpsdy->hybrid[i]);
        omega++;
        halfp++;
        diver++;
        uwind++;
        vwind++;
      }
   }

/* interpolate to full levels  */

   omega = Omega->hybrid;
#if defined (SX)
#pragma vdir nodep
#endif
#if defined (__uxp__)
#pragma loop novrec
#endif
   for (i = 0; i < Dim3GP; i++)
      omega[i] = 0.5 * (omega[i] + omega[i+DimGP]);

/* compute full level part of vertical velocity */

   omega = Omega->hybrid;
   halfp = HalfPress->hybrid;
   fullp = FullPress->hybrid;
   uwind = u_wind->hybrid;
   vwind = v_wind->hybrid;

   for (j = 0; j < Levels; j++) {
      DeltaHybrid = vct[Levels+j+2] - vct[Levels+j+1];
      if (DeltaHybrid) {
         Cterm = vct[j+1] * vct[Levels+j+1] - vct[j] * vct[Levels+j+2];
#if defined (__uxp__)
#pragma loop novrec
#endif
         for (i = 0; i < DimGP; i++) {
            if (Cterm != 0.0) Pterm = Cterm /
               (halfp[DimGP] - *halfp) * log(halfp[DimGP] / *halfp);
            else Pterm = 0.0;

           *omega += *fullp *
              (*uwind * dpsdx->hybrid[i] + *vwind * dpsdy->hybrid[i])
              / (halfp[DimGP] - *halfp) * (DeltaHybrid + Pterm);
           omega++;
           halfp++;
           fullp++;
           uwind++;
           vwind++;
         }
      }
      else {
         omega += DimGP;
         halfp += DimGP;
         fullp += DimGP;
         uwind += DimGP;
         vwind += DimGP;
      }
   }
}

void CheckDependencies(void)
{

           u_wind->needed |= (  Divergence->needed &&   !Divergence->detected) ||
                             (   Vorticity->needed &&    !Vorticity->detected) ||
                             (     velopot->needed &&      !velopot->detected) ||
                             (      stream->needed &&       !stream->detected) ||
                             (       Omega->needed &&        !Omega->detected) ||
                             (       speed->needed &&        !speed->detected) ||
                             (      v_wind->needed &&       !v_wind->detected);

           v_wind->needed |= (  Divergence->needed &&   !Divergence->detected) ||
                             (   Vorticity->needed &&    !Vorticity->detected) ||
                             (     velopot->needed &&      !velopot->detected) ||
                             (      stream->needed &&       !stream->detected) ||
                             (       Omega->needed &&        !Omega->detected) ||
                             (       speed->needed &&        !speed->detected) ||
                             (      u_wind->needed &&       !u_wind->detected);

       Divergence->needed |= (      u_wind->needed &&       !u_wind->detected) ||
                             (      v_wind->needed &&       !v_wind->detected) ||
                             (      stream->needed &&       !stream->detected) ||
                             (       Omega->needed &&        !Omega->detected) ||
                             (     velopot->needed &&      !velopot->detected);

        Vorticity->needed |= (      u_wind->needed &&       !u_wind->detected) ||
                             (      v_wind->needed &&       !v_wind->detected) ||
                             (      stream->needed &&       !stream->detected) ||
                             (     velopot->needed &&      !velopot->detected);

  if ( AnalysisData )
    {
      if ( Debug )
	fprintf(STDOUT,"CheckDependencies: AnalysisData = %1d\n",AnalysisData);
      Geopotential->needed = (GeopotHeight->needed && !GeopotHeight->detected) ||
	                      Geopotential->selected;
         Rhumidity->needed = (    Humidity->selected &&   !Humidity->detected) ||
                                 Rhumidity->selected;
          Humidity->needed = (   Rhumidity->selected &&  !Rhumidity->detected) ||
                                  Humidity->selected;
       Temperature->needed = (   Rhumidity->needed &&    !Rhumidity->detected) ||
                             (    Humidity->needed &&     !Humidity->detected) ||
                               Temperature->selected;
    }
  else
    {
         Humidity->needed |= (GeopotHeight->needed && !GeopotHeight->detected) ||
                             (   Rhumidity->needed &&    !Rhumidity->detected);

      Temperature->needed |= (GeopotHeight->needed && !GeopotHeight->detected) ||
	                     (   Rhumidity->needed &&    !Rhumidity->detected) ||
                             (         SLP->needed &&          !SLP->detected) ||
                             (      ThetaF->needed &&       !ThetaF->detected);
    }

          All[222].needed = (   all_water->needed &&    !all_water->detected) ||
                            (   low_water->needed &&    !low_water->detected) ||
                            (   mid_water->needed &&    !mid_water->detected) ||
                            (   hih_water->needed &&    !hih_water->detected) ||
                                  All[222].selected;
          All[223].needed = (   low_cloud->needed &&    !low_cloud->detected) ||
                            (   mid_cloud->needed &&    !mid_cloud->detected) ||
                            (   hih_cloud->needed &&    !hih_cloud->detected) ||
                                  All[223].selected;
          All[176].needed = (    net_heat->needed &&     !net_heat->detected) ||
                            (     net_bot->needed &&      !net_bot->detected) ||
                            (     net_atm->needed &&      !net_atm->detected) ||
                            (      sw_atm->needed &&       !sw_atm->detected) ||
                            (  sw_bot_clf->needed &&   !sw_bot_clf->detected) ||
                                  All[176].selected;
          All[177].needed = (    net_heat->needed &&     !net_heat->detected) ||
                            (     net_bot->needed &&      !net_bot->detected) ||
                            (     net_atm->needed &&      !net_atm->detected) ||
                            (      lw_atm->needed &&       !lw_atm->detected) ||
                            (  lw_bot_clf->needed &&   !lw_bot_clf->detected) ||
                                  All[177].selected;
          All[178].needed = (     net_top->needed &&      !net_top->detected) ||
                            (     net_atm->needed &&      !net_atm->detected) ||
                            (     net_clf->needed &&      !net_clf->detected) ||
                            (      sw_atm->needed &&       !sw_atm->detected) ||
                            (      sw_clf->needed &&       !sw_clf->detected) ||
                            (  sw_top_clf->needed &&   !sw_top_clf->detected) ||
                            ( net_top_clf->needed &&  !net_top_clf->detected) ||
                                  All[178].selected;
          All[179].needed = (     net_top->needed &&      !net_top->detected) ||
                            (     net_atm->needed &&      !net_atm->detected) ||
                            (     net_clf->needed &&      !net_clf->detected) ||
                            (      lw_atm->needed &&       !lw_atm->detected) ||
                            (      lw_clf->needed &&       !lw_clf->detected) ||
                            (  lw_top_clf->needed &&   !lw_top_clf->detected) ||
                            ( net_top_clf->needed &&  !net_top_clf->detected) ||
                                  All[179].selected;
          All[185].needed = (  sw_bot_clf->needed &&   !sw_bot_clf->detected) ||
                                  All[185].selected;
          All[186].needed = (  lw_bot_clf->needed &&   !lw_bot_clf->detected) ||
                                  All[186].selected;
          All[187].needed = ( net_top_clf->needed &&  !net_top_clf->detected) ||
                            (  sw_top_clf->needed &&   !sw_top_clf->detected) ||
                                  All[187].selected;
          All[188].needed = ( net_top_clf->needed &&  !net_top_clf->detected) ||
                            (  lw_top_clf->needed &&   !lw_top_clf->detected) ||
                                  All[188].selected;
}

void CheckContent(void)
{
  static char func[] = "CheckContent";
  int code;

  if ( Debug ) Message(func, "TermCount = %d MeanCount = %d", TermCount, MeanCount);

  for (code = 0; code < 256; code++)
    {
      if ( code == GEOSCODE ) continue;
      if ( code ==  SLPCODE ) continue;
      if ( code ==    ZCODE ) continue;
      if ( code ==  STRCODE ) continue;
      if ( code ==  VELCODE ) continue;
      if ( code ==    UCODE ) continue;
      if ( code ==    VCODE ) continue;
      if ( code ==    WCODE ) continue;
      if ( code ==   RHCODE ) continue;
      if ( code ==   LCCODE ) continue;
      if ( code ==   MCCODE ) continue;
      if ( code ==   HCCODE ) continue;
      if ( code ==   PSCODE ) continue;
      if ( code ==   SHCODE )
	{
	  if ( All[code].needed && !All[code].selected &&
	       All[code].spectral == NULL &&
	       All[code].hybrid   == NULL &&
	       All[code].swapped  == FALSE )
	    {
	      Warning(func, "No humidity in data file, set to zero !");
	      All[code].needed = FALSE;
	    }
	}
      else
	{
	  if ( All[code].needed &&
	       All[code].spectral == NULL &&
	       All[code].hybrid   == NULL &&
	       All[code].swapped  == FALSE )
	    {
	      if ( labort )
		Error(func, "Code  %3d not found", code);
	      else
		Warning(func, "Code  %3d not found", code);
	    }
	}
    }
}

void
CheckAnalyses(void)
{
  int code;

  if ( Debug )
    fprintf(STDOUT, "CheckAnalyses: %d %d\n", TermCount, MeanCount);

  for (code = 0; code < 256; code++)
    if ( All[code].needed && code != DIVCODE && code != VORCODE &&
	 code != STRCODE && code != UCODE   && code != SHCODE  &&
	 code != VELCODE && code != VCODE   && code != RHCODE  &&
	 code != ZCODE   && code != PSCODE  &&
	 All[code].spectral == NULL && All[code].grid == NULL )
      {
	fprintf (STDOUT," CheckAnalyses:\n");
	fprintf (STDOUT,"\n ****** E R R O R ******\n");
	fprintf (STDOUT," * Code  %3d not found *\n",code);
	fprintf (STDOUT," ***********************\n");
	exit (1);
      }
}

void AnalysisDependencies(void)
{
   int code;

   for (code = 0; code < CODES; code++)
     All[code].needed = All[code].selected;

           LnPs->needed |=         Ps->needed;

       Humidity->needed |=    Rhumidity->needed;
    Temperature->needed |=    Rhumidity->needed;
      Rhumidity->needed |=     Humidity->needed;
    Temperature->needed |=     Humidity->needed;

   Geopotential->needed |= GeopotHeight->needed;

   CheckDependencies();
}

void EchamDependencies(void)
{
   int code;

   for (code = 0; code < CODES; code++)
     All[code].needed = All[code].selected;

   if (type >= 30) LnPs->needed = TRUE;

   if (type >= 20) {
             LnPs->needed |=        (Omega->needed && !Omega->detected) ||    ThetaF->needed;
            dpsdx->needed |=        (Omega->needed && !Omega->detected);
            dpsdy->needed |=        (Omega->needed && !Omega->detected);
             LnPs->needed |=        (SLP->needed && !SLP->detected);
             LnPs->needed |= (GeopotHeight->needed  && !GeopotHeight->detected) || 
	                        (Rhumidity->needed  && !   Rhumidity->detected);
   }

   LnPs->needed |= HalfPress->needed || dpsdx->needed ||
                   (Ps->needed && !Ps->detected);

   All[139].needed |=      ThetaF->selected;
   All[142].needed |=      precip->selected ||   net_water->selected ||
                      fresh_water->selected || surf_runoff->selected;
   All[143].needed |=      precip->selected ||   net_water->selected ||
                      fresh_water->selected || surf_runoff->selected;
   All[146].needed |=    net_heat->selected;     /* sensible heat */
   All[147].needed |=    net_heat->selected;     /* latent   heat */
   All[160].needed |=   net_water->selected;     /* Runoff        */
   All[176].needed |=    net_heat->selected ||
                          net_bot->selected ||
                          net_atm->selected ||      sw_atm->selected ||
                       sw_bot_clf->selected ;
   All[177].needed |=    net_heat->selected ||
                          net_bot->selected ||
                          net_atm->selected ||      lw_atm->selected ||
                       lw_bot_clf->selected ;
   All[178].needed |=     net_top->selected ||
                          net_atm->selected ||      sw_atm->selected ||
                          net_clf->selected ||      sw_clf->selected ||
                      net_top_clf->selected ||  sw_top_clf->selected;
   All[179].needed |=     net_top->selected ||
                          net_atm->selected ||      lw_atm->selected ||
                          net_clf->selected ||      lw_clf->selected ||
                      net_top_clf->selected ||  lw_top_clf->selected;
   All[182].needed |=   net_water->selected ||
                      fresh_water->selected || surf_runoff->selected;
   All[185].needed |=  sw_bot_clf->selected;
   All[186].needed |=  lw_bot_clf->selected;
   All[187].needed |= net_top_clf->selected ||  sw_top_clf->selected;
   All[188].needed |= net_top_clf->selected ||  lw_top_clf->selected;
   All[218].needed |=    net_heat->selected;     /* snow melt         */
   All[220].needed |=    net_heat->selected;     /* residual term     */
   All[221].needed |= surf_runoff->selected;     /* snow depth change */
   All[222].needed |=   all_water->selected ||   low_water->selected ||
                        mid_water->selected ||   hih_water->selected;
   All[223].needed |=   low_cloud->selected ||
                        mid_cloud->selected ||   hih_cloud->selected;
   All[224].needed |=     net_clf->selected ||      sw_clf->selected;
   All[225].needed |=     net_clf->selected ||      lw_clf->selected;

   CheckDependencies();
}

/* ========================= */
/* Write hybrid-level fields */
/* ========================= */

void wrihyb (int code, double *field, int dim, int lev)
{
   int k, l, SingleLevel, LayerLevel;

   SingleLevel = lev == 1;
   LayerLevel  = All[code].nlay > 0;
   /*
   fprintf (stdout, "lev = %d  nrql = %d\n", lev, nrql);
   */
   if (field == NULL) {
      fprintf(STDOUT," wrihyb:\n");
      fprintf(STDOUT,"\n *** Code %d was not found ***\n",code);
      exit(1);
   }
   SrvHead[0] = code;
   if (SingleLevel) {
      SrvHead[1] = 0;
      wrifor(field, dim, SrvHead, ofileID);
   }
   else if (LayerLevel) {
     for (l = 0; l < Layers; l++)
       for (k = 0; k < All[code].nlay; k++)
	 if (All[code].layers[k] == (int) layer[l])
	   {
	     SrvHead[1] = All[code].layers[k];
	     wrifor(field + k * dim, dim, SrvHead, ofileID);
	   }
   }
   else {
     for (l = 0; l < nrql; l++) {
       SrvHead[1] = k = (int) level[l];
       wrifor(field + (k-1) * dim, dim, SrvHead, ofileID);
     }
   }
}

/* =============================== */
/* Write hybrid zonal mean section */
/* =============================== */

void wrihzm (int code, double *field, int lev)
{
#ifdef FORTRAN_OUTPUT
   int j,dim;
   int headcontrol, fieldcontrol;
   float *FloatField, *ff;
   double  *rf;
#endif
   int k, l, SingleLevel;

   SingleLevel = lev == 1;

   if (field == NULL) {
      if (TermCount == 1) {
         fprintf(STDOUT," wrihzm:\n");
         fprintf(STDOUT,"\n *** Code %d was not found ***\n",code);
      }
      return;
   }

   SrvHead[0] = code;
   SrvHead[1] = -1;
   SrvHead[4] = Latitudes ;
   if (SingleLevel) SrvHead[5] = 1;
   else             SrvHead[5] = nrql;

#ifdef FORTRAN_OUTPUT
   if (SingleLevel) dim = Latitudes;
   else             dim = Latitudes * nrql;
   ff = FloatField = FloatAlloc(dim,"wrihzm.field");
   headcontrol  = sizeof(SrvHead);
   fieldcontrol = dim * sizeof(float);
   if (SingleLevel)
      for (j = 0; j < Latitudes; j++) ff[j] = field[j];
   else for (l = 0; l < nrql; l++) {
      k = (int) level[l];
      rf = field + DimFC * (k-1);
      for (j = 0; j < Latitudes; j++) ff[j] = rf[j];
      ff += Latitudes;
   }
   (void) fwrite (&headcontrol,  sizeof(int),     1, ofileID);
   (void) fwrite (SrvHead,       sizeof(SrvHead[0]), 8, ofileID);
   (void) fwrite (&headcontrol,  sizeof(int),     1, ofileID);
   (void) fwrite (&fieldcontrol, sizeof(int),     1, ofileID);
   (void) fwrite (FloatField,    sizeof(float), dim, ofileID);
   (void) fwrite (&fieldcontrol, sizeof(int),     1, ofileID);
   FloatField = FreeMem(FloatField);
#else
   (void) fwrite((char *)SrvHead, sizeof(SrvHead[0]), 8, ofileID);
   if (SingleLevel)
      (void) fwrite((char *)field, sizeof(double), Latitudes, ofileID);
   else for (l = 0; l < nrql; l++) {
      k = (int) level[l];
      (void) fwrite((char *)field + (k-1) * DimFC, sizeof(double), Latitudes, ofileID);
   }
#endif
}

/* =========================== */
/* Write pressure level fields */
/* =========================== */

void wripre(int code, double *field, int dim, int lev)
{
   int  l, k;

   if (field == NULL) {
      if (TermCount == 1) {
         fprintf(STDOUT," wripre:\n");
         fprintf(STDOUT,"\n *** Code %d was not found ***\n", SrvHead[0]);
      }
      return;
   }

   if (lev == 1) {
      SrvHead[1] = 0;
      wrifor (field, dim, SrvHead, ofileID);
   }
   else if (All[code].nlay > 0) {
     for (l = 0; l < Layers; l++)
       for (k = 0; k < All[code].nlay; k++)
	 if (All[code].layers[k] == (int) layer[l])
	   {
	     SrvHead[1] = All[code].layers[k];
	     wrifor(field + k * dim, dim, SrvHead, ofileID);
	   }
   }   
   else
      for (l = 0; l < nrql; l++) {
        SrvHead[1] = (int) level[l];
        wrifor (field, dim, SrvHead, ofileID);
        field += dim;
      }
}

/* ========================================= */
/* Write pressure level means and deviations */
/* ========================================= */

void wripmv(int code, double *mfield, double *vfield, int dim, int lev)
{
   int  l, k;

   if (mfield == NULL && vfield == NULL) {
      if (TermCount == 1) {
         fprintf(STDOUT," wripmv:\n");
         fprintf(STDOUT,"\n *** Code %d was not found ***\n", SrvHead[0]);
      } 
      return;
   }
   if (lev == 1) {
      SrvHead[1] = 0;
      if (Mean != 2) wrifor (mfield, dim, SrvHead, ofileID);
      if (Mean >= 2) wrifor (vfield, dim, SrvHead, ofileID);
   }
   else if (All[code].nlay > 0) {
     for (l = 0; l < Layers; l++)
       for (k = 0; k < All[code].nlay; k++)
	 if (All[code].layers[k] == (int) layer[l])
	   {
	     SrvHead[1] = All[code].layers[k];
	     if (Mean != 2) wrifor (mfield + k * dim, dim, SrvHead, ofileID);
	     if (Mean >= 2) wrifor (vfield + k * dim, dim, SrvHead, ofileID);
	   }
   }   
   else for (l = 0; l < nrql; l++) {
      SrvHead[1] = (int) level[l];
      if (Mean != 2) wrifor (mfield, dim, SrvHead, ofileID);
      if (Mean >= 2) wrifor (vfield, dim, SrvHead, ofileID);
      mfield += dim;
      vfield += dim;
   }
}

/* ======================================= */
/* Write pressure level zonal mean section */
/* ======================================= */

void wripzm(double *field, int lev, int lats, int dimfc)
{
#ifdef FORTRAN_OUTPUT
   int j,dim;
   int headcontrol,fieldcontrol;
   float *FloatField,*ff;
#endif
   int i,SingleLevel;

   SingleLevel = lev == 1;

   if (field == NULL) {
      if (TermCount == 1) {
         fprintf(STDOUT," wrihzm:\n");
         fprintf(STDOUT,"\n *** Code %d was not found ***\n", SrvHead[0]);
      }
      return;
   }

   SrvHead[1] = -1;
   SrvHead[4] = lats ;
   if (SingleLevel) SrvHead[5] = 1;
   else             SrvHead[5] = nrql;

#ifdef FORTRAN_OUTPUT
   if (SingleLevel) dim = lats;
   else             dim = lats * nrql;
   ff = FloatField = FloatAlloc(dim,"wripzm.field");
   headcontrol = sizeof(SrvHead);
   fieldcontrol = dim * sizeof(float);
   for (i = 0; i<nrql; i++) {
      for (j = 0; j < lats; j++) ff[j] = field[j];
      ff    += lats;
      field += dimfc;
   }
   (void) fwrite (&headcontrol , sizeof(int)    ,  1, ofileID);
   (void) fwrite (SrvHead      , sizeof(SrvHead[0]),  8, ofileID);
   (void) fwrite (&headcontrol , sizeof(int)    ,  1, ofileID);
   (void) fwrite (&fieldcontrol, sizeof(int)    ,  1, ofileID);
   (void) fwrite (FloatField   , sizeof(float)  ,dim, ofileID);
   (void) fwrite (&fieldcontrol, sizeof(int)    ,  1, ofileID);
   FloatField = FreeMem(FloatField);
#else
   (void) fwrite((char *)SrvHead,sizeof(SrvHead[0]),8, ofileID);
   if (SingleLevel)
      (void) fwrite ((char *)field, sizeof(double), lats, ofileID);
   else for (i = 0; i < nrql; i++) {
      (void) fwrite ((char *)field, sizeof(double), lats, ofileID);
      field += dimfc;
   }
#endif
}

void SPuv2SPdv(void)
{
   int i, FieldLength=0;
   double *Div, *DivOut, *Vor, *VorOut;

   Div = DivOut = Divergence->spectral;
   Vor = VorOut =  Vorticity->spectral;
   FieldLength = DimFC * nrql;

   if (u_wind->fourier == NULL)
       u_wind->fourier = DoubleAlloc (FieldLength, "u_wind->fourier");
   if (v_wind->fourier == NULL)
       v_wind->fourier = DoubleAlloc (FieldLength, "v_wind->fourier");

   sp2fc (u_wind->spectral, u_wind->fourier, poli,
          nrql, Latitudes, Fouriers, Truncation);
   sp2fc (v_wind->spectral, v_wind->fourier, poli,
          nrql, Latitudes, Fouriers, Truncation);
   uv2dv (u_wind->fourier, v_wind->fourier, Div, Vor,
          pol2, pol3, nrql, Latitudes, Truncation);

   u_wind->fourier = FreeMem(u_wind->fourier);
   v_wind->fourier = FreeMem(v_wind->fourier);

   for (i = 0; i < nrql; ++i) {
      ChangeTruncation (Div, Truncation, DivOut, TruncationOut);
      ChangeTruncation (Vor, Truncation, VorOut, TruncationOut);
      Div    += DimSP;
      Vor    += DimSP;
      DivOut += DimSPOut;
      VorOut += DimSPOut;
   }
}

/* HUMTEST */
void FCrh2FCsh(void)
{
   int FieldLength;

   FieldLength = DimGPOut * nrql;

   if (Rhumidity->grid == NULL)
       Rhumidity->grid = DoubleAlloc (FieldLength, "Rhumidity->grid");
   if (Temperature->grid == NULL)
       Temperature->grid = DoubleAlloc (FieldLength, "Temperature->grid");
   if (Humidity->grid == NULL)
       Humidity->grid = DoubleAlloc (FieldLength, "Humidity->grid");

   fc2gp (Rhumidity->fourier, Rhumidity->grid,
          LatitudesOut, LongitudesOut, Rhumidity->plev, FouriersOut);
   fc2gp (Temperature->fourier, Temperature->grid,
          LatitudesOut, LongitudesOut, Temperature->plev, FouriersOut);

   rh2sh (Humidity->grid,Rhumidity->grid,Temperature->grid,nrql,
          DimGPOut, level);

   gp2fc (Humidity->grid,Humidity->fourier,
          LatitudesOut, LongitudesOut, Humidity->plev, FouriersOut);

   Humidity->grid    = FreeMem (Humidity->grid);
   Rhumidity->grid   = FreeMem (Rhumidity->grid);
   Temperature->grid = FreeMem (Temperature->grid);
}

void FCsh2FCrh(void)
{
   int FieldLength;

   FieldLength = DimGPOut * nrql;

   if (Rhumidity->grid == NULL)
       Rhumidity->grid = DoubleAlloc (FieldLength, "Rhumidity->grid");
   if (Temperature->grid == NULL)
       Temperature->grid = DoubleAlloc (FieldLength, "Temperature->grid");
   if (Humidity->grid == NULL)
       Humidity->grid = DoubleAlloc (FieldLength, "Humidity->grid");

   fc2gp (Humidity->fourier, Humidity->grid,
          LatitudesOut, LongitudesOut, Humidity->plev, FouriersOut);
   fc2gp (Temperature->fourier, Temperature->grid,
          LatitudesOut, LongitudesOut, Temperature->plev, FouriersOut);

   sh2rh (Humidity->grid, Rhumidity->grid, Temperature->grid, nrql,
          DimGPOut, level, FullPress->hybrid);

   gp2fc (Rhumidity->grid, Rhumidity->fourier,
          LatitudesOut, LongitudesOut, Rhumidity->plev, FouriersOut);

   Humidity->grid    = FreeMem (Humidity->grid);
   Rhumidity->grid   = FreeMem (Rhumidity->grid);
   Temperature->grid = FreeMem (Temperature->grid);
}
/* ENDE HUMTEST */

/* ************** */
/* Swap functions */
/* ************** */

void swapout(int code, int dim, int type)
{
   size_t nwrite;
   FILE *swapfile;

   swapfile = fopen (swapname[code], "w");

   if (swapfile == NULL) {
      All[code].swapped = FALSE;
   } else {
      if (type == 20)
         nwrite = fwrite(All[code].hybrid,sizeof(double),dim,swapfile);
      else
         nwrite = fwrite(All[code].grid  ,sizeof(double),dim,swapfile);
      fclose(swapfile); swapfile = NULL;
      if (nwrite == dim) {
         All[code].swapped = TRUE;
         if (type == 20)
            All[code].hybrid = FreeMem(All[code].hybrid);
         else
            All[code].grid   = FreeMem(All[code].grid);
      } else {
         All[code].swapped = FALSE;
      }
   }
   if (Debug) {
      fprintf(STDOUT,"+   swapout: swap out %8d of %3d ",dim,code);
      if (All[code].swapped) fprintf(STDOUT," done\n");
      else                   fprintf(STDOUT," failed\n");
   }
}

void swapin (int code, int dim, int type)
{
   size_t nread;
   FILE *swapfile;

   swapfile = fopen (swapname[code], "rb");

   if (swapfile != NULL && All[code].swapped) {
      if (type == 20) {
         All[code].hybrid = DoubleAlloc(dim,FieldName(code,"hybrid"));
         nread = fread(All[code].hybrid,sizeof(double),dim,swapfile);
      } else {
         All[code].grid   = DoubleAlloc(dim,FieldName(code,"grid"));
         nread = fread(All[code].grid  ,sizeof(double),dim,swapfile);
      }
      fclose(swapfile); swapfile = NULL;
      if (nread == dim) {
         All[code].swapped = TRUE;
      } else {
         All[code].swapped = FALSE;
         if (type == 20)
            All[code].hybrid = FreeMem(All[code].hybrid);
         else
            All[code].grid   = FreeMem(All[code].grid);
      }
   }
   if (Debug) {
      fprintf(STDOUT,"+   swapin : swap in  %8d of %3d",dim,code);
      if (All[code].swapped) fprintf(STDOUT," done\n");
      else                   fprintf(STDOUT," failed\n");
   }
}

void AnalysisProcess(void)
{
   int code,k,l;
   int FieldLength = 0;
   extern double Grav;

   MeanCount++;
   TermCount++;

   if (MeanCount == 1) {
      CheckAnalyses();
      StartDate = OldDate;
   }
   if (TermCount  > 120) Debug = 0;

/* ============================== */
/* Computations in spectral space */
/* ============================== */

   if (Temperature->needed) {
      Temperature->hlev =    2;
      Temperature->plev = nrql;
      Temperature->sfit = TRUE;
   }

   if (Omega->needed) {
      Omega->hlev =    2;
      Omega->plev = nrql;
      Omega->sfit = TRUE;
   }

   if (GeopotHeight->needed && !GeopotHeight->detected) {
      GeopotHeight->hlev =    2;
      GeopotHeight->plev = nrql;
      GeopotHeight->sfit = TRUE;
   }
   if (GeopotHeight->needed && !GeopotHeight->detected && Geopotential->detected) {
      if (GeopotHeight->spectral == NULL)
          GeopotHeight->spectral = DoubleAlloc(DimSPOut*nrql,
                                             "GeopotHeight.spectral");
      MultVectorScalar(GeopotHeight->spectral,Geopotential->spectral,
                       RG,DimSPOut*nrql);
      Geopotential->needed = Geopotential->selected;
   }

   if (type == 50 && Humidity->needed && Humidity->spectral == NULL) {
      Humidity->plev = nrql;
      Humidity->sfit = TRUE;
      Humidity->spectral = DoubleAlloc(DimSPOut*nrql,"Humidity.spectral");

/*
      SPrh2SPsh();
*/

      Rhumidity->needed = Rhumidity->selected;
      Temperature->needed = Temperature->selected;

   }

   if (u_wind->spectral && v_wind->spectral &&
       (Divergence->needed || Vorticity->needed)) {
      Divergence->hlev = Vorticity->hlev =    2;
      Divergence->plev = Vorticity->plev = nrql;
      Divergence->sfit = Vorticity->sfit = TRUE;
      if (Divergence->spectral == NULL)
          Divergence->spectral = DoubleAlloc(DimSP*nrql,"Divergence.spectral");
      if ( Vorticity->spectral == NULL)
           Vorticity->spectral = DoubleAlloc(DimSP*nrql,"Vorticity.spectral");
      SPuv2SPdv();
   }

   if ((u_wind->needed && !u_wind->detected) ||
       (v_wind->needed && !v_wind->detected)) {
      u_wind->hlev = v_wind->hlev =    2;
      u_wind->plev = v_wind->plev = nrql;
      u_wind->sfit = v_wind->sfit = TRUE;
      if (u_wind->spectral == NULL)
          u_wind->spectral = DoubleAlloc(DimSPOut*nrql,"u_wind.spectral");
      if (v_wind->spectral == NULL)
          v_wind->spectral = DoubleAlloc(DimSPOut*nrql,"v_wind.spectral");
      dv2uv(Divergence->spectral,Vorticity->spectral,
                u_wind->spectral,   v_wind->spectral,
                        dv2uv_f1,           dv2uv_f2,
                        TruncationOut,DimSPOut,nrql);
   }

   if (velopot->needed && !velopot->detected) {
      velopot->hlev =    2;
      velopot->plev = nrql;
      velopot->sfit = TRUE;
      if (velopot->spectral == NULL)
          velopot->spectral = DoubleAlloc(DimSPOut*nrql,"velopot.spectral");
      dv2ps(Divergence->spectral,velopot->spectral,nrql,TruncationOut);
   }

   if (stream->needed && !stream->detected) {
      stream->hlev =    2;
      stream->plev = nrql;
      stream->sfit = TRUE;
      if (stream->spectral == NULL)
          stream->spectral = DoubleAlloc(DimSPOut*nrql,"stream.spectral");
      dv2ps(Vorticity->spectral,stream->spectral,nrql,TruncationOut);
   }

/* --------------------------------- */
/* Output of spectral fields in GRIB */
/* --------------------------------- */

   Representation = REP_SPECTRAL;
   if (AnaLevelFactor == 100) {
       LevelType  = LEV_ISOBARIC;
   } else
   if (AnaLevelFactor ==   1) {
       LevelType  = LEV_HYBRID;
   }
   if (type == 50 && Grib) {
      for (code = 0; code < CODES; code++)
      if (All[code].selected && All[code].spectral)
      for (l = 0; l < All[code].plev; ++l) {
         k = (int) level[l] / AnaLevelFactor;
         PutGrib(ofileID, All[code].spectral+l*DimSPOut, code, k);
      }
      FreeSpectral();
      return;
   }

/* ------------------------- */
/* Output of spectral fields */
/* ------------------------- */

   if (type == 50 && !Grib) {
      SrvHead[4] = DimSPOut;
      SrvHead[5] =     1;
      for (code = 0; code < CODES; code++)
      if (All[code].selected) {
         SrvHead[0] = code;
         wripre(code, All[code].spectral,DimSPOut,2);
      }
      FreeSpectral();
      return;
   }

/* =============================== */
/* Transformation to fourier space */
/* Computations in fourier space   */
/* =============================== */

   if (type >= 60) {
      for (code = 0; code < CODES; code++)
      if (All[code].needed && All[code].spectral) {
         if (All[code].fourier == NULL) {
            FieldLength = All[code].plev * DimFCOut;
            All[code].fourier = DoubleAlloc(FieldLength,
                                          FieldName(code,"fourier"));
         }
         sp2fc(All[code].spectral,All[code].fourier,poliOut,
          All[code].plev,LatitudesOut,FouriersOut,TruncationOut);
      }
      if ( u_wind->needed && u_wind->fourier )
         scaluv(u_wind->fourier, rclatOut, LatitudesOut, FouriersOut*nrql);
      if ( v_wind->needed && v_wind->fourier )
         scaluv(v_wind->fourier, rclatOut, LatitudesOut, FouriersOut*nrql);

 /* HUMTEST */
      if (type < 70 && Humidity->needed && Humidity->fourier == NULL) {
         Humidity->plev = nrql;
         Humidity->sfit = TRUE;
         Humidity->fourier = DoubleAlloc(DimFCOut*nrql,"Humidity.fourier");

         FCrh2FCsh();

         Rhumidity->needed = Rhumidity->selected;
         Temperature->needed = Temperature->selected;
      }

      if (type < 70 && Rhumidity->needed && Rhumidity->fourier == NULL) {
         Rhumidity->plev = nrql;
         Rhumidity->sfit = TRUE;
         Rhumidity->fourier = DoubleAlloc(DimFCOut*nrql,
                                        "Rhumidity.fourier");

         FCsh2FCrh();

         Humidity->needed = Humidity->selected;
         Temperature->needed = Temperature->selected;
      }
 /* ENDE HUMTEST */
   }

   FreeSpectral();

/* ------------------------ */
/* Output of fourier fields */
/* ------------------------ */

   if (type == 60) {
      SrvHead[4] = LatitudesOut ;
      SrvHead[5] = FouriersOut  ;
      for (code = 0; code < CODES; code++)
      if (All[code].selected) {
         SrvHead[0] = code;
         wripre(code, All[code].fourier,DimFCOut,2);
      }
      FreeFourier();
      return;
   }

/* --------------------- */
/* Output of zonal means */
/* --------------------- */

   if (type == 61) {
      for (code = 0; code < CODES; code++)
      if (All[code].selected) {
         SrvHead[0] = code;
         wripzm(All[code].fourier,All[code].hlev,LatitudesOut,DimFCOut);
      }
      FreeFourier();
      return;
   }

/* ============================ */
/* Transformation to gridpoints */
/* ============================ */

   if (Ps->selected && LnPs->grid) {
      if (Ps->grid == NULL) Ps->grid = DoubleAlloc(DimGP,"Ps");
      for (l = 0; l < DimGP; l++) Ps->grid[l] = exp(LnPs->grid[l]);
   }

   if (type >= 70) {
      for (code = 0; code < CODES; code++)
      if (All[code].needed && All[code].fourier) {
         if (All[code].grid == NULL) {
            FieldLength = All[code].plev * DimGPOut;
            All[code].grid = DoubleAlloc(FieldLength,FieldName(code,"grid"));
         }

         fc2gp(All[code].fourier,All[code].grid,
               LatitudesOut,LongitudesOut,All[code].plev,FouriersOut);
      }
   }

   FreeFourier();

/* HUMTEST */
/* ============================== */
/* Computation in gridpoint space */
/* ============================== */

   if(Rhumidity->needed && Rhumidity->grid == NULL) {
      Rhumidity->plev = nrql;
      Rhumidity->sfit = TRUE;
      Rhumidity->grid = DoubleAlloc(DimGPOut*nrql,"Rhumidity.grid");
      sh2rh(Humidity->grid,Rhumidity->grid,Temperature->grid,nrql,
            DimGPOut, level, FullPress->hybrid);
      Humidity->needed = Humidity->selected;
      Temperature->needed = Temperature->selected;
   }

   if(Humidity->needed && Humidity->grid == NULL) {
      Humidity->plev = nrql;
      Humidity->sfit = TRUE;
      Humidity->grid = DoubleAlloc(DimGPOut*nrql,"Humidity.grid");
      rh2sh(Humidity->grid,Rhumidity->grid,Temperature->grid,nrql,
            DimGPOut, level);
      Rhumidity->needed = Rhumidity->selected;
      Temperature->needed = Temperature->selected;
   }
/* HUMTEST ENDE */

/* =========================== */
/* Computation of Means        */
/* =========================== */

   if (Mean)
   for (code = 0; code < CODES; code++)
   if (All[code].needed && All[code].grid) {
      FieldLength = DimGPOut * All[code].plev;
      if (All[code].mean == NULL) {
         All[code].mean = DoubleAlloc(FieldLength,FieldName(code,"mean"));
      }
      if (MeanCount == 1) RealCopy(All[code].mean,All[code].grid,FieldLength);
      else               AddVector(All[code].mean,All[code].grid,FieldLength);
      if (EndOfInterval) MultVectorScalar(All[code].mean,All[code].mean,
					  1.0/MeanCount,FieldLength);
   }

/* ======================== */
/* Computation of Variances */
/* ======================== */

   if (Mean > 1)
   for (code = 0; code < CODES; code++)
   if (All[code].needed && All[code].mean) {
      FieldLength = DimGPOut * All[code].plev;
      if (All[code].var == NULL)
         All[code].var = DoubleAlloc(FieldLength,FieldName(code,"var"));
      if (MeanCount == 1) IniQuaSum(All[code].var,All[code].grid,FieldLength);
      else                AddQuaSum(All[code].var,All[code].grid,FieldLength);
      if (EndOfInterval) VarQuaSum(All[code].var,All[code].mean,FieldLength,
				   MeanCount);
   }

   if (Mean && !EndOfInterval) {
      FreeGrid();
      return;
   }

/* ---------------------------------------------------- */
/* Output of pressure level means and variances in GRIB */
/* ---------------------------------------------------- */

   Representation = REP_GAUSS;
   if (type == 70 && Grib && Mean && EndOfInterval) {
      for (code = 0; code < CODES; code++) {
         if (All[code].selected && All[code].mean != NULL) {
            for (l = 0; l < All[code].plev; l++) {
               if (All[code].surf) {
                  if (All[code].leveltype == LEV_DOWN || All[code].leveltype == LEV_HEIGHT) {
                     LevelType = All[code].leveltype;
                     k         = All[code].layer;
                  }
                  else {
                     LevelType = LEV_SURFACE;
                     k         = 0;
                  }
               }  else {
                  if (level[l] < LEV99_MARGIN) {
                     k = (int) level[l];
                     LevelType = LEV_99;
                  }  else {
                     k = (int) level[l] / AnaLevelFactor;
                     if (AnaLevelFactor == 100) LevelType = LEV_ISOBARIC;
                     else                       LevelType = LEV_HYBRID;
                  }
               }
               if (Mean != 2)
                  PutGrib(ofileID, All[code].mean+l*DimGPOut, code, k);
               if (Mean >= 2)
                  PutGrib(ofileID, All[code].var +l*DimGPOut, code, k);
            }
            All[code].grid = FreeMem(All[code].grid);
         }
      }
      FreeGrid();
      return;
   }

/* -------------------------------------------- */
/* Output of pressure level means and variances */
/* -------------------------------------------- */

   if (type == 70 && !Grib && Mean && EndOfInterval) {
      SrvHead[4] = LongitudesOut;
      SrvHead[5] = LatitudesOut;
      for (code = 0; code < CODES; code++) {
         if (All[code].selected) {
            SrvHead[0] = code;
            wripmv(code, All[code].mean,All[code].var,DimGPOut,All[code].plev);
         }
         All[code].grid = FreeMem(All[code].grid);
      }
      FreeGrid();
      return;
   }


/* -------------------------------------- */
/* Output of pressure level grids in GRIB */
/* -------------------------------------- */

   if (type == 70 && Grib) {
      for (code = 0; code < CODES; code++) {
         if (All[code].selected && All[code].grid) {
            for (l = 0; l < All[code].plev; l++) {
               if (All[code].surf) {
                  if (All[code].leveltype == LEV_DOWN || All[code].leveltype == LEV_HEIGHT) {
                     LevelType = All[code].leveltype;
                     k         = All[code].layer;
                  }
                  else {
                     LevelType = LEV_SURFACE;
                     k         = 0;
                  }
               }  else {
                  if (level[l] < LEV99_MARGIN) {
                     k = (int) level[l];
                     LevelType = LEV_99;
                  }  else {
                     k = (int) level[l] / AnaLevelFactor;
                     if (AnaLevelFactor == 100) LevelType = LEV_ISOBARIC;
                     else                       LevelType = LEV_HYBRID;
                  }
               }
               PutGrib(ofileID, All[code].grid+l*DimGPOut, code, k);
            }
            All[code].grid = FreeMem(All[code].grid);
         }
      }
      FreeGrid();
      return;
   }

/* ------------------------------ */
/* Output of pressure level grids */
/* ------------------------------ */

   if (type == 70 && !Grib) {
      SrvHead[4] = LongitudesOut;
      SrvHead[5] = LatitudesOut;
      for (code = 0; code < CODES; code++)
      if (All[code].selected && All[code].grid) {
         SrvHead[0] = code;
         wripre(code, All[code].grid,DimGPOut,2);
      }
      FreeGrid();
      return;
   }
}

void EchamProcess(void)
{
   int code,k,l,i;
   int FieldLength = 0;
   int VarLev;
   double *var_spectral, *div_spectral, *vor_spectral;
   extern double Grav;

   MeanCount++;
   TermCount++;

   if (MeanCount == 1) {
      CheckContent();
      StartDate = OldDate;
   }
   if (TermCount  > 120) Debug = 0;

/* ====================================== */
/* Swap out hybrid fields if Swap == TRUE */
/* ====================================== */

   if (Swap)
     for (code = 0; code < CODES; code++) {
       if (All[code].hybrid && All[code].hlev > 1)
	 swapout(code,DimGP*All[code].hlev,20);
     }

/* ============================== */
/* Computations in spectral space */
/* ============================== */

   if ((u_wind->needed   ||  v_wind->needed) &&
      (!u_wind->detected || !v_wind->detected)) {
      u_wind->hlev = v_wind->hlev = Divergence->hlev;
      u_wind->plev = v_wind->plev = Divergence->plev;
      u_wind->sfit = v_wind->sfit = TRUE;
      u_wind->spectral = DoubleAlloc(Dim3SP,"u_wind.spectral");
      v_wind->spectral = DoubleAlloc(Dim3SP,"v_wind.spectral");

      if (Divergence->spectral == NULL) gp2sp((int)DIVCODE);
      if (Vorticity->spectral == NULL)  gp2sp((int)VORCODE);

      dv2uv(Divergence->spectral,Vorticity->spectral,
                u_wind->spectral,   v_wind->spectral,
                        dv2uv_f1,           dv2uv_f2,
                   Truncation,DimSP,Divergence->hlev);
   }

   if (velopot->needed && !velopot->detected && type < 30) {
      velopot->hlev = Divergence->hlev;
      velopot->plev = Divergence->plev;
      velopot->spectral = DoubleAlloc(Dim3SP,"velopot.spectral");

      if (Divergence->spectral == NULL) gp2sp((int)DIVCODE);

      dv2ps(Divergence->spectral,velopot->spectral,
            Divergence->hlev,Truncation);
   }

   if (stream->needed && !stream->detected && type < 30) {
      stream->hlev = Vorticity->hlev;
      stream->plev = Vorticity->plev;
      stream->spectral = DoubleAlloc(Dim3SP,"stream.spectral");

      if (Vorticity->spectral == NULL)  gp2sp((int)VORCODE);

      dv2ps(Vorticity->spectral,stream->spectral,
            Vorticity->hlev,Truncation);
   }

   Vorticity->needed = Vorticity->selected && type <= 20;
   if (Vorticity->spectral && !Vorticity->needed)
       Vorticity->spectral = FreeMem( Vorticity->spectral);

   Divergence->needed = (Divergence->selected && type <= 20) ||
                        (Omega->needed && !Omega->detected);
   if (Divergence->spectral && !Divergence->needed)
       Divergence->spectral = FreeMem(Divergence->spectral);

/* ------------------------- */
/* Output of spectral fields */
/* ------------------------- */

   Representation = REP_SPECTRAL;
   if (type == 0) {
      SrvHead[4] = DimSP;
      SrvHead[5] =     1;
      for (code = 0; code < CODES; code++)
      if (All[code].selected)
         wrihyb(code,All[code].spectral,DimSP,All[code].hlev);
      FreeSpectral();
      return;
   }

/* ============================= */
/* Computations in fourier space */
/* ============================= */

   if (type >= 10) {
      for (code = 0; code < CODES; code++)
      if (All[code].needed && All[code].spectral) {
         FieldLength = All[code].hlev * DimFC;
         All[code].fourier = DoubleAlloc(FieldLength,
                                       FieldName(code,"fourier"));
         sp2fc(All[code].spectral,All[code].fourier,poli,
               All[code].hlev,Latitudes,Fouriers,Truncation);
         if (code != LNPSCODE)
            All[code].spectral = FreeMem(All[code].spectral);
      }

/*    if (type < 60) poli = FreeMem(poli); */
/*    if (type < 50) pol2 = FreeMem(pol2); */
/*    if (type < 50) pol3 = FreeMem(pol3); */


      if ( u_wind->needed && u_wind->fourier )
         scaluv(u_wind->fourier, rclat, Latitudes, Fouriers*Levels);
      if ( v_wind->needed && v_wind->fourier )
         scaluv(v_wind->fourier, rclat, Latitudes, Fouriers*Levels);

      if (dpsdx->needed) {
         dpsdx->hlev = 1;
         dpsdx->plev = 1;
         dpsdx->sfit = FALSE;
         dpsdx->fourier = DoubleAlloc(DimFC,"dpsdx->fourier");
         if (LnPs->fourier == NULL)  gp2sp((int)LNPSCODE);
         Derivate(LnPs->fourier,dpsdx->fourier,1,Waves,Latitudes,DerivationFactor);
      }
      if (dpsdy->needed) {
         dpsdy->hlev = 1;
         dpsdy->plev = 1;
         dpsdy->sfit = FALSE;
         dpsdy->fourier = DoubleAlloc(DimFC,"dpsdy->fourier");
         if (LnPs->spectral == NULL)  gp2sp((int)LNPSCODE);
         sp2fc(LnPs->spectral,dpsdy->fourier,pdev,
               dpsdy->hlev,Latitudes,Fouriers,Truncation);
      }
   }
   FreeSpectral();

/* ------------------------ */
/* Output of fourier fields */
/* ------------------------ */

   if (type == 10) {
      SrvHead[4] = Latitudes ;
      SrvHead[5] = Fouriers  ;
      for (code = 0; code < CODES; code++)
      if (All[code].selected)
         wrihyb(code,All[code].fourier,DimFC,All[code].hlev);
      FreeFourier();
      return;
   }

/* --------------------- */
/* Output of zonal means */
/* --------------------- */

   if (type == 11) {
      for (code = 0; code < CODES; code++)
      if (All[code].selected)
         wrihzm(code,All[code].fourier,All[code].hlev);
      FreeFourier();
      return;
   }

/* ============================ */
/* Transformation to gridpoints */
/* ============================ */

   if (type >= 20) {
      for (code = 0; code < CODES; code++)
      if (All[code].needed && All[code].fourier) {
         FieldLength = DimGP * All[code].hlev;
         All[code].hybrid = DoubleAlloc(FieldLength,
                                      FieldName(code,"hybrid"));
         fc2gp(All[code].fourier,All[code].hybrid,
               Latitudes,Longitudes,All[code].hlev,Fouriers);

         if (All[code].hlev > 1 && Swap)
            swapout(code,FieldLength,20);
         All[code].fourier = FreeMem(All[code].fourier);
      }

      if (LnPs->hybrid && Ps->hybrid == NULL) {
         if (Ps->hybrid == NULL) Ps->hybrid = DoubleAlloc(DimGP,"Ps");
         for (l = 0; l < DimGP; l++) Ps->hybrid[l] = exp(LnPs->hybrid[l]);
      }

/*    if (All[1].hybrid) {
         if (Ps->hybrid == NULL) Ps->hybrid = DoubleAlloc(DimGP,"Ps");
         RealCopy(Ps->hybrid,All[1].hybrid,DimGP);
      } */

      LnPs->needed = LnPs->selected;

      if (Orography == NULL) {
         Orography = DoubleAlloc(DimGP , "Orography");
         if (Geopotential->hybrid)
              RealCopy(Orography,Geopotential->hybrid,DimGP);
         else {
            if (Geopotential->selected || type >= 30) {
               fprintf(STDOUT,"!! Orography not found -"
                      " using zero orography !!\n");
               RealZero(Orography,DimGP);
            }
         }
      }
      Geopotential->needed = Geopotential->selected;

      if ((GeopotHeight->needed  && !GeopotHeight->detected) ||
	  (SLP->needed && !SLP->detected) || ThetaF->needed ||
          HalfPress->needed || Rhumidity->needed || (Omega->needed && !Omega->detected) ||
          type >= 30) {
         if  (FullPress->hybrid == NULL)
              FullPress->hybrid = DoubleAlloc(Dim3GP,"FullPress->hybrid");
         HalfPress->hlev  = Levels+1;
         HalfPress->plev  = nrql;
         HalfPress->sfit  = FALSE;
         if  (HalfPress->hybrid == NULL)
              HalfPress->hybrid = DoubleAlloc(Dim3GP+DimGP,"HalfPress->hybrid");
         presh(FullPress->hybrid,HalfPress->hybrid,vct,Ps->hybrid,Levels,DimGP,Dim3GP);

         if (Swap) {
            swapout(277,Dim3GP+DimGP,20);
            swapout(278,Dim3GP,20);
         }
      }

      if (unitsel > 2) FullPress->hybrid = FreeMem(FullPress->hybrid);

      if (ThetaF->needed) {
         ThetaF->hlev = Levels;
         ThetaF->plev = nrql;
         ThetaF->sfit = TRUE;
         if (ThetaF->hybrid == NULL)
            ThetaF->hybrid = DoubleAlloc(Dim3GP,"ThetaF.hybrid");
         if (ThetaH->hybrid == NULL)
            ThetaH->hybrid = DoubleAlloc(Dim3GP,"ThetaH.hybrid");
         theta(ThetaF->hybrid, ThetaH->hybrid, HalfPress->hybrid, Ps->hybrid,
               Temperature->hybrid, Ts->hybrid, Levels, DimGP, Dim3GP);
      }

      if (GeopotHeight->needed  && !GeopotHeight->detected) {
         GeopotHeight->hlev = Levels+1;
         GeopotHeight->plev = nrql;
         GeopotHeight->sfit = TRUE;
         GeopotHeight->hybrid = DoubleAlloc(Dim3GP+DimGP,"GeopotHeight.hybrid");

         if (Swap) {
            swapin (TCODE,Dim3GP,20);
            swapin (SHCODE,Dim3GP,20);
            swapin (277,Dim3GP+DimGP,20);
         }

         RealCopy(GeopotHeight->hybrid+Dim3GP,Orography,DimGP);
         MakeGeopotHeight(GeopotHeight->hybrid,Temperature->hybrid,
                          Humidity->hybrid,HalfPress->hybrid,DimGP,Levels);

         Humidity->needed = Humidity->selected;

         if (Swap) {
            swapout(ZCODE,Dim3GP+DimGP,20);
              HalfPress->hybrid = FreeMem(  HalfPress->hybrid);
               Humidity->hybrid = FreeMem(   Humidity->hybrid);
            Temperature->hybrid = FreeMem(Temperature->hybrid);
         }
      }

      if (dpsdx->needed || dpsdy->needed)
      for (l = 0; l < DimGP; l++) {
         dpsdx->hybrid[l] *= Ps->hybrid[l];
         dpsdy->hybrid[l] *= Ps->hybrid[l];
      }

      if (Omega->needed && !Omega->detected) {
         Omega->hlev = Levels+1;
         Omega->plev = nrql;
         Omega->sfit = TRUE;
         Omega->hybrid = DoubleAlloc(Dim3GP+DimGP,"Omega.hybrid");
         if (Swap) {
            swapin (UCODE,Dim3GP,20);
            swapin (VCODE,Dim3GP,20);
            swapin (155,Dim3GP,20);
            swapin (277,Dim3GP+DimGP,20);
            swapin (278,Dim3GP,20);
         }
         OMEGA();
         dpsdx->needed = dpsdx->selected;
         dpsdy->needed = dpsdy->selected;

         if (Swap) {
            swapout(135,Dim3GP+DimGP,20);
              FullPress->hybrid = FreeMem(  FullPress->hybrid);
              HalfPress->hybrid = FreeMem(  HalfPress->hybrid);
             Divergence->hybrid = FreeMem( Divergence->hybrid);
                 v_wind->hybrid = FreeMem(     v_wind->hybrid);
                 u_wind->hybrid = FreeMem(     u_wind->hybrid);
         }
      }

      if (speed->needed) {
         speed->hlev = Levels;
         speed->plev = nrql;
         speed->sfit = TRUE;
         speed->hybrid = DoubleAlloc(Dim3GP,"speed.hybrid");

         if (Swap) {
            swapin (UCODE,Dim3GP,20);
            swapin (VCODE,Dim3GP,20);
         }

         Speed(speed->hybrid,u_wind->hybrid,v_wind->hybrid,Dim3GP);

         if (Swap) {
            swapout(259,Dim3GP,20);
            v_wind->hybrid = FreeMem(     v_wind->hybrid);
            u_wind->hybrid = FreeMem(     u_wind->hybrid);
         }
      }

      if (Rhumidity->needed) {
         Rhumidity->hlev = Levels;
         Rhumidity->plev = nrql;
         Rhumidity->sfit = FALSE;
         Rhumidity->hybrid = DoubleAlloc(Dim3GP,"Rhumidity.hybrid");

         if (Swap) {
            swapin (TCODE,Dim3GP,20);
            swapin (SHCODE,Dim3GP,20);
            swapin (278,Dim3GP,20);
         }
         sh2rh(Humidity->hybrid,Rhumidity->hybrid, Temperature->hybrid,Levels,
	       DimGPOut, level, FullPress->hybrid);

         if (Swap) {
            swapout(RHCODE,Dim3GP,20);
              FullPress->hybrid = FreeMem(  FullPress->hybrid);
            Temperature->hybrid = FreeMem(Temperature->hybrid);
               Humidity->hybrid = FreeMem(   Humidity->hybrid);
         }

         Temperature->needed = Temperature->selected;
            Humidity->needed =    Humidity->selected;
      }

      if (SLP->needed && !SLP->detected) {
         SLP->hlev = 1;
         SLP->plev = 1;
         SLP->sfit = TRUE;
         SLP->hybrid = DoubleAlloc(DimGP,"SLP.hybrid");

         if (Swap) {
            swapin (277,   Dim3GP+DimGP, 20);
            swapin (278,   Dim3GP, 20);
            swapin (TCODE, Dim3GP, 20);
         }

         Extrap(SLP->hybrid,HalfPress->hybrid + Dim3GP,
                FullPress->hybrid + Dim3GP - DimGP , Orography,
                Temperature->hybrid + Dim3GP - DimGP , DimGP);
         Temperature->needed = Temperature->selected ||
                               GeopotHeight->selected;

         if (Swap) {
            Temperature->hybrid = FreeMem(Temperature->hybrid);
              FullPress->hybrid = FreeMem(  FullPress->hybrid);
              HalfPress->hybrid = FreeMem(  HalfPress->hybrid);
         }
      }

      if (precip->needed) {
         precip->hlev = precip->plev = 1;
         precip->sfit = FALSE;
         precip->hybrid = DoubleAlloc(DimGP,"precip->hybrid");
         Add2Vectors(precip->hybrid,All[142].hybrid,All[143].hybrid,DimGP);
      }

      if (net_top->needed) {
         net_top->hlev = net_top->plev = 1;
         net_top->sfit = FALSE;
         net_top->hybrid = DoubleAlloc(DimGP,"net_top->hybrid");
         Add2Vectors(net_top->hybrid,All[178].hybrid,All[179].hybrid,DimGP);
      }

      if (net_bot->needed) {
         net_bot->hlev = net_bot->plev = 1;
         net_bot->sfit = FALSE;
         net_bot->hybrid = DoubleAlloc(DimGP,"net_bot->hybrid");
         Add2Vectors(net_bot->hybrid,All[176].hybrid,All[177].hybrid,DimGP);
      }

      if (net_heat->needed) {
         net_heat->hlev = net_heat->plev = 1;
         net_heat->sfit = FALSE;
         net_heat->hybrid = DoubleAlloc(DimGP,"net_heat->hybrid");
         MultVectorScalar(net_heat->hybrid,All[218].hybrid,L_times_rhoH2O,DimGP);
         Add2Vectors(net_heat->hybrid,net_heat->hybrid,All[176].hybrid,DimGP);
         Add2Vectors(net_heat->hybrid,net_heat->hybrid,All[177].hybrid,DimGP);
         Add2Vectors(net_heat->hybrid,net_heat->hybrid,All[146].hybrid,DimGP);
         Add2Vectors(net_heat->hybrid,net_heat->hybrid,All[147].hybrid,DimGP);
         Sub2Vectors(net_heat->hybrid,net_heat->hybrid,All[220].hybrid,DimGP);
      }

      if (net_water->needed) {
         net_water->hlev = net_water->plev = 1;
         net_water->sfit = FALSE;
         net_water->hybrid = DoubleAlloc(DimGP,"net_water->hybrid");
         Sub2Vectors(net_water->hybrid,All[182].hybrid  ,All[160].hybrid,DimGP);
         Add2Vectors(net_water->hybrid,net_water->hybrid,All[142].hybrid,DimGP);
         Add2Vectors(net_water->hybrid,net_water->hybrid,All[143].hybrid,DimGP);
      }

      if (Swap)
         if (low_water->needed || mid_water->needed ||
             hih_water->needed || all_water->needed)
            swapin (222,Dim3GP,20);

      if (!low_water->detected && low_water->needed) {
         low_water->hlev = low_water->plev = 1;
         low_water->sfit = FALSE;
         low_water->hybrid = DoubleAlloc(DimGP,"low_water->hybrid");
         LayerWater(All[222].hybrid,low_water->hybrid, 75000.,101300., DimGP, HalfLevels, vct);
      }

      if (!mid_water->detected && mid_water->needed) {
         mid_water->hlev = mid_water->plev = 1;
         mid_water->sfit = FALSE;
         mid_water->hybrid = DoubleAlloc(DimGP,"mid_water->hybrid");
         LayerWater(All[222].hybrid,mid_water->hybrid, 46000., 73000., DimGP, HalfLevels, vct);
      }

      if (!hih_water->detected && hih_water->needed) {
         hih_water->hlev = hih_water->plev = 1;
         hih_water->sfit = FALSE;
         hih_water->hybrid = DoubleAlloc(DimGP,"hih_water->hybrid");
         LayerWater(All[222].hybrid,hih_water->hybrid,  5000., 44000., DimGP, HalfLevels, vct);
      }

      if (!all_water->detected && all_water->needed) {
         all_water->hlev = all_water->plev = 1;
         all_water->sfit = FALSE;
         all_water->hybrid = DoubleAlloc(DimGP,"all_water->hybrid");
         LayerWater(All[222].hybrid,all_water->hybrid,  5000.,101300., DimGP, HalfLevels, vct);
      }


      if (Swap)
         if (low_water->needed || mid_water->needed ||
             hih_water->needed || all_water->needed)
            All[222].hybrid = FreeMem(All[222].hybrid);

      if (Swap)
         if (low_cloud->needed || mid_cloud->needed ||
             hih_cloud->needed )
            swapin (223,Dim3GP,20);

      if (!low_cloud->detected && low_cloud->needed) {
         low_cloud->hlev = low_cloud->plev = 1;
         low_cloud->sfit = FALSE;
         low_cloud->hybrid = DoubleAlloc(DimGP,"low_cloud->hybrid");
         LayerCloud(All[223].hybrid,low_cloud->hybrid, 75000.,101300., DimGP, HalfLevels, vct);
      }

      if (!mid_cloud->detected && mid_cloud->needed) {
         mid_cloud->hlev = mid_cloud->plev = 1;
         mid_cloud->sfit = FALSE;
         mid_cloud->hybrid = DoubleAlloc(DimGP,"mid_cloud->hybrid");
         LayerCloud(All[223].hybrid,mid_cloud->hybrid, 46000., 73000., DimGP, HalfLevels, vct);
      }

      if (!hih_cloud->detected && hih_cloud->needed) {
         hih_cloud->hlev = hih_cloud->plev = 1;
         hih_cloud->sfit = FALSE;
         hih_cloud->hybrid = DoubleAlloc(DimGP,"hih_cloud->hybrid");
         LayerCloud(All[223].hybrid,hih_cloud->hybrid,  5000., 44000., DimGP, HalfLevels, vct);
      }

      if (Swap)
         if (low_cloud->needed || mid_cloud->needed ||
             hih_cloud->needed )
            All[223].hybrid = FreeMem(All[223].hybrid);

      if (sw_clf->needed) {
         sw_clf->hlev = sw_clf->plev = 1;
         sw_clf->sfit = FALSE;
         sw_clf->hybrid = DoubleAlloc(DimGP,"sw_clf->hybrid");
         Sub2Vectors(sw_clf->hybrid,All[178].hybrid,All[224].hybrid,DimGP);
      }

      if (sw_bot_clf->needed && !sw_bot_clf->detected) {
         sw_bot_clf->hlev = sw_bot_clf->plev = 1;
         sw_bot_clf->sfit = FALSE;
         sw_bot_clf->hybrid = DoubleAlloc(DimGP,"sw_bot_clf->hybrid");
         Sub2Vectors(sw_bot_clf->hybrid,All[176].hybrid,All[185].hybrid,DimGP);
      }

      if (sw_top_clf->needed && !sw_top_clf->detected) {
         sw_top_clf->hlev = sw_top_clf->plev = 1;
         sw_top_clf->sfit = FALSE;
         sw_top_clf->hybrid = DoubleAlloc(DimGP,"sw_top_clf->hybrid");
         Sub2Vectors(sw_top_clf->hybrid,All[178].hybrid,All[187].hybrid,DimGP);
      }

      if (lw_clf->needed) {
         lw_clf->hlev = lw_clf->plev = 1;
         lw_clf->sfit = FALSE;
         lw_clf->hybrid = DoubleAlloc(DimGP,"lw_clf->hybrid");
         Sub2Vectors(lw_clf->hybrid,All[179].hybrid,All[225].hybrid,DimGP);
      }

      if (lw_bot_clf->needed && !lw_bot_clf->detected) {
         lw_bot_clf->hlev = lw_bot_clf->plev = 1;
         lw_bot_clf->sfit = FALSE;
         lw_bot_clf->hybrid = DoubleAlloc(DimGP,"lw_bot_clf->hybrid");
         Sub2Vectors(lw_bot_clf->hybrid,All[177].hybrid,All[186].hybrid,DimGP);
      }

      if (lw_top_clf->needed && !lw_top_clf->detected) {
         lw_top_clf->hlev = lw_top_clf->plev = 1;
         lw_top_clf->sfit = FALSE;
         lw_top_clf->hybrid = DoubleAlloc(DimGP,"lw_top_clf->hybrid");
         Sub2Vectors(lw_top_clf->hybrid,All[179].hybrid,All[188].hybrid,DimGP);
      }

      if (net_clf->needed) {
         net_clf->hlev = net_clf->plev = 1;
         net_clf->sfit = FALSE;
         net_clf->hybrid = DoubleAlloc(DimGP,"net_clf->hybrid");
         Add2Vectors(net_clf->hybrid,All[178].hybrid,All[179].hybrid,DimGP);
         Sub2Vectors(net_clf->hybrid,net_clf->hybrid,All[224].hybrid,DimGP);
         Sub2Vectors(net_clf->hybrid,net_clf->hybrid,All[225].hybrid,DimGP);
      }

      if (sw_atm->needed) {
         sw_atm->hlev = sw_atm->plev = 1;
         sw_atm->sfit = FALSE;
         sw_atm->hybrid = DoubleAlloc(DimGP,"sw_atm->hybrid");
         Sub2Vectors(sw_atm->hybrid,All[178].hybrid,All[176].hybrid,DimGP);
      }

      if (lw_atm->needed) {
         lw_atm->hlev = lw_atm->plev = 1;
         lw_atm->sfit = FALSE;
         lw_atm->hybrid = DoubleAlloc(DimGP,"lw_atm->hybrid");
         Sub2Vectors(lw_atm->hybrid,All[179].hybrid,All[177].hybrid,DimGP);
      }

      if (net_atm->needed) {
         net_atm->hlev = net_atm->plev = 1;
         net_atm->sfit = FALSE;
         net_atm->hybrid = DoubleAlloc(DimGP,"net_atm->hybrid");
         Add2Vectors(net_atm->hybrid,All[178].hybrid,All[179].hybrid,DimGP);
         Sub2Vectors(net_atm->hybrid,net_atm->hybrid,All[176].hybrid,DimGP);
         Sub2Vectors(net_atm->hybrid,net_atm->hybrid,All[177].hybrid,DimGP);
      }

      if (surf_runoff->needed) {
         surf_runoff->hlev = surf_runoff->plev = 1;
         surf_runoff->sfit = FALSE;
         surf_runoff->hybrid = DoubleAlloc(DimGP,"surf_runoff->hybrid");
         Sub2Vectors(surf_runoff->hybrid,All[182].hybrid,All[221].hybrid,DimGP);
         Add2Vectors(surf_runoff->hybrid,surf_runoff->hybrid,All[142].hybrid,DimGP);
         Add2Vectors(surf_runoff->hybrid,surf_runoff->hybrid,All[143].hybrid,DimGP);
      }

      if (fresh_water->needed) {
         fresh_water->hlev = fresh_water->plev = 1;
         fresh_water->sfit = FALSE;
         fresh_water->hybrid = DoubleAlloc(DimGP,"fresh_water->hybrid");
         Add2Vectors(fresh_water->hybrid,All[142].hybrid,All[143].hybrid,DimGP);
         Add2Vectors(fresh_water->hybrid,fresh_water->hybrid,All[182].hybrid,DimGP);
      }
   }

   FreeFourier();

/* ============================= */
/* Means on hybrid grids         */
/* ============================= */

   if (Mean && type == 20) {
      for (code = 0; code < CODES; code++)
      if (All[code].selected &&
         (All[code].hybrid   || All[code].swapped)) {
         if (All[code].mean == NULL)
            All[code].mean = DoubleAlloc(All[code].hlev
                           * DimGP,FieldName(code,"mean"));
         if (Mean > 1 && All[code].var == NULL)
            All[code].var = DoubleAlloc(All[code].hlev
                          * DimGP,FieldName(code,"var"));

         if (Swap) {
            swapin (code,All[code].hlev* DimGP,20);
         }

         if (MeanCount == 1) {
            RealCopy(All[code].mean,All[code].hybrid,
                     All[code].hlev * DimGP);
            if (Mean > 1) IniQuaSum(All[code].var,All[code].hybrid,
                                    All[code].hlev * DimGP);
         }  else {
            AddVector(All[code].mean,All[code].hybrid,
                      All[code].hlev * DimGP);
            if (Mean > 1) AddQuaSum(All[code].var,All[code].hybrid,
                                    All[code].hlev * DimGP);
         }
         if (EndOfInterval) {
            MultVectorScalar(All[code].hybrid,All[code].mean,
                             1.0/MeanCount,All[code].hlev * DimGP);
            if (Mean > 1) VarQuaSum(All[code].var,All[code].hybrid,
                                    All[code].hlev * DimGP,MeanCount);
         }
         if (All[code].swapped)
            All[code].hybrid = FreeMem(All[code].hybrid);
      }
   }


/* ------------------------------------ */
/* Output of hybrid level grids in GRIB */
/* ------------------------------------ */

   Representation = REP_GAUSS;
   if (type == 20 && Grib && (Mean == 0 || EndOfInterval)) {
      for (code = 0; code < CODES; code++)
      if (All[code].selected) {
         if (Swap) {
            swapin (code,All[code].hlev* DimGP,20);
         }
         if (All[code].hlev == 1 || All[code].nlay) {
	   if (All[code].nlay) {
	     LevelType = All[code].leveltype;
	     for (l = 0; l < Layers; l++)
	       for (k = 0; k < All[code].nlay; k++)
		 if (All[code].layers[k] == (int) layer[l])
		   {
		     if (Mean != 2) PutGrib(ofileID, All[code].hybrid+k*DimGP, code, All[code].layers[k]);
		     if (Mean >= 2) PutGrib(ofileID, All[code].var   +k*DimGP, code, All[code].layers[k]);
		   }
	   }
	   else {
	     if (All[code].leveltype == LEV_DOWN || All[code].leveltype == LEV_HEIGHT) {
	       LevelType = All[code].leveltype;
	       k         = All[code].layer;
	     }
	     else {
	       LevelType = LEV_SURFACE;
	       k         = 0;
	     }
	     if (Mean != 2) PutGrib(ofileID, All[code].hybrid, code, k);
	     if (Mean >= 2) PutGrib(ofileID, All[code].var   , code, k);
	   }
         }
         else for (l = 0; l < nrql; l++) {
                LevelType = LEV_HYBRID;
                k = (int) level[l];
                if (Mean != 2) PutGrib(ofileID, All[code].hybrid+(k-1)*DimGP, code, k);
                if (Mean >= 2) PutGrib(ofileID, All[code].var   +(k-1)*DimGP, code, k);
         }
         if (All[code].swapped)
             All[code].hybrid = FreeMem(All[code].hybrid);
      }
   }

/* ---------------------------- */
/* Output of hybrid level grids */
/* ---------------------------- */

   if (type == 20 && !Grib && (Mean == 0 || EndOfInterval)) {
      SrvHead[4] = Longitudes;
      SrvHead[5] = Latitudes;
      for (code = 0; code < CODES; code++)
      if (All[code].selected) {
         if (Swap) {
            swapin (code,All[code].hlev* DimGP,20);
         }
         if (Mean != 2) wrihyb(code,All[code].hybrid,DimGP,All[code].hlev);
         if (Mean >= 2) wrihyb(code,All[code].var   ,DimGP,All[code].hlev);
         if (All[code].swapped)
             All[code].hybrid = FreeMem(All[code].hybrid);
      }
   }


   if (type == 20) {
      FreeHybrid();
      return;
   }

/* ====================================== */
/* Vertical interpolation / extrapolation */
/* ====================================== */


   if (type >= 30) {
      if (Swap) {
         swapin (278,Dim3GP,20);
      }
      if (vert_index == NULL)
          vert_index = IntAlloc(nrql*DimGP,"vert_index");
      if (unitsel) {
         if (p_of_height == NULL)
             p_of_height = DoubleAlloc(DimGP*nrql,"p_of_height");
         h2p(p_of_height,level,DimGP,nrql);
         geninz(vert_index,p_of_height,
                FullPress->hybrid,DimGP,nrql,Levels);
      }
      else
         genind(vert_index,level,
                FullPress->hybrid,DimGP,nrql,Levels);

      for (code = 0; code < CODES; code++)
      if (All[code].needed &&
         (All[code].hybrid || All[code].swapped)) {
         if (Swap) {
            swapin (code,All[code].hlev* DimGP,20);
         }
         if (All[code].hlev == 1 || All[code].nlay) {
            if (All[code].grid) FreeMem(All[code].grid);
            All[code].grid = All[code].hybrid;
            All[code].hybrid = NULL;
         }  else {
            if (All[code].grid == NULL) {
               FieldLength = DimGP * nrql;
               All[code].grid = DoubleAlloc(FieldLength,
                                          FieldName(code,"grid"));
            }
            switch (code) {
              case TCODE: if (Swap) swapin (277,Dim3GP+DimGP,20);
                          if (unitsel)
                             Interpolate_T_Z(Orography,
                                             Temperature->hybrid,
                                             Temperature->grid,
                                             FullPress->hybrid,
                                             HalfPress->hybrid,
                                             vert_index,p_of_height,
                                             Levels, nrql, DimGP);
                          else
                             Interpolate_T(  Orography,
                                             Temperature->hybrid,
                                             Temperature->grid,
                                             FullPress->hybrid,
                                             HalfPress->hybrid,
                                             vert_index,
                                             level, Levels, nrql, DimGP);
                          if (HalfPress->swapped)
                              HalfPress->hybrid = FreeMem(HalfPress->hybrid);
                          break;
              case ZCODE: if (Swap) swapin (277,Dim3GP+DimGP,20);
                          if (unitsel)
                             Interpolate_Z_Z(Orography,
                                             GeopotHeight->hybrid,
                                             GeopotHeight->grid,
                                             FullPress->hybrid,
                                             HalfPress->hybrid,
                                             vert_index,p_of_height,
                                             Temperature->hybrid,
                                             Levels, nrql, DimGP);
                          else
                             Interpolate_Z(  Orography,
                                             GeopotHeight->hybrid,
                                             GeopotHeight->grid,
                                             FullPress->hybrid,
                                             HalfPress->hybrid,
                                             vert_index,
                                             Temperature->hybrid,
                                             level, Levels, nrql, DimGP);
                          if (ECMWFMars)
                             for (i = 0; i < FieldLength; i++)
                                 GeopotHeight->grid[i] *= Grav;

                          if (HalfPress->swapped)
                              HalfPress->hybrid = FreeMem(HalfPress->hybrid);
                          break;
                default:  if (unitsel)
                             Interpolate_X_Z(All[code].hybrid,
                                             All[code].grid,
                                             FullPress->hybrid,
                                             vert_index,p_of_height,
                                             nrql, DimGP, Dim3GP);
                          else
                             Interpolate_X(  All[code].hybrid,
                                             All[code].grid,
                                             FullPress->hybrid,
                                             vert_index,
                                             level, nrql, DimGP, Dim3GP);
            }
            if (code != TCODE)
               All[code].hybrid = FreeMem(All[code].hybrid);
            if (Swap) {
               swapout(code,DimGP*nrql,30);
            }
         }
      }
   }
   Temperature->needed = Temperature->selected;
   FreeHybrid();
   if (HalfPress->hybrid)
       HalfPress->hybrid = FreeMem(HalfPress->hybrid);

/* -------------------------------------- */
/* Output of pressure level grids in GRIB */
/* -------------------------------------- */

   if (type == 30 && Grib && Mean == 0) {
      for (code = 0; code < CODES; code++)
      if (All[code].selected &&
         (All[code].grid     || All[code].swapped)) {
         swapin (code,DimGP*nrql,30);

         if (All[code].hlev == 1 || All[code].nlay) {
	   if (All[code].nlay) {
	     LevelType = All[code].leveltype;
	     for (l = 0; l < Layers; l++)
	       for (k = 0; k < All[code].nlay; k++)
		 if (All[code].layers[k] == (int) layer[l])
		   PutGrib(ofileID, All[code].grid+k*DimGP, code, All[code].layers[k]);
	   }
	   else {
	     if (All[code].leveltype == LEV_DOWN || All[code].leveltype == LEV_HEIGHT) {
               LevelType = All[code].leveltype;
               k         = All[code].layer;
	     }
	     else {
               LevelType = LEV_SURFACE;
               k         = 0;
	     }
	     PutGrib(ofileID, All[code].grid, code, k);
	   }
         }
         else for (l = 0; l < nrql; l++) {
            switch (unitsel) {
            case 0: if (level[l] < LEV99_MARGIN) {
                       LevelType = LEV_99;
                       VarLev    = (int) level[l];
                    } else {
                       LevelType = LEV_ISOBARIC;
                       VarLev    = (int) level[l]/100;
                    }
                    break;
            case 1:
            case 2: LevelType = LEV_HEIGHT;
                    VarLev    = (int) level[l];
                    break;
            default: fprintf(STDOUT,"Non allowed value for unitsel: %d\n",
                             unitsel);
                     exit(1);
            }
            PutGrib(ofileID, All[code].grid+l*DimGP, code, VarLev);
         }
         if (All[code].swapped)
             All[code].grid   = FreeMem(All[code].grid);
      }
      FreeGrid();
      return;
   }

/* ------------------------------ */
/* Output of pressure level grids */
/* ------------------------------ */

   if (type == 30 && !Grib && Mean == 0) {
      SrvHead[4] = Longitudes;
      SrvHead[5] = Latitudes;
      for (code = 0; code < CODES; code++)
      if (All[code].selected) {
         swapin (code,DimGP*All[code].plev,30);

         SrvHead[0] = code;
         wripre(code, All[code].grid,DimGP,All[code].plev);
         if (All[code].swapped)
             All[code].grid = FreeMem(All[code].grid);
      }
      FreeGrid();
      return;
   }

/* =========================== */
/* Computation of Means        */
/* =========================== */

   if (type >= 30 && Mean)
   for (code = 0; code < CODES; code++)
   if (All[code].needed &&
      (All[code].grid   || All[code].swapped)) {
      FieldLength = DimGP * All[code].plev;

      swapin (code,FieldLength,30);

      if (All[code].mean == NULL)
         All[code].mean = DoubleAlloc(FieldLength,FieldName(code,"mean"));
      if (MeanCount == 1) RealCopy(All[code].mean,All[code].grid,
                                   FieldLength);
      else                AddVector(All[code].mean,All[code].grid,
                                   FieldLength);
      if (EndOfInterval) MultVectorScalar(All[code].mean,All[code].mean,
                                       1.0/MeanCount,FieldLength);
      if (All[code].swapped)
          All[code].grid = FreeMem(All[code].grid);
   }

/* ======================== */
/* Computation of Variances */
/* ======================== */

   if (type >= 30 && Mean > 1)
   for (code = 0; code < CODES; code++)
   if (All[code].needed && All[code].mean) {
      FieldLength = DimGP * All[code].plev;

      swapin (code,FieldLength,30);

      if (All[code].var == NULL)
         All[code].var = DoubleAlloc(FieldLength,FieldName(code,"var"));
      if (MeanCount == 1) IniQuaSum(All[code].var,All[code].grid,
                                    FieldLength);
      else                AddQuaSum(All[code].var,All[code].grid,
                                    FieldLength);
      if (EndOfInterval)     VarQuaSum(All[code].var,All[code].mean,
                                    FieldLength,MeanCount);
      if (All[code].swapped)
          All[code].grid = FreeMem(All[code].grid);
   }

   if (Mean && !EndOfInterval) {
      FreeGrid();
      return;
   }

/* ---------------------------------------------------- */
/* Output of pressure level means and variances in GRIB */
/* ---------------------------------------------------- */

   if (type == 30 && Grib && Mean) {
      for (code = 0; code < CODES; code++)
      if (All[code].selected && All[code].mean != NULL) {
         if (All[code].hlev == 1 || All[code].nlay) {
	   if (All[code].nlay) {
	     LevelType = All[code].leveltype;
	     for (l = 0; l < Layers; l++)
	       for (k = 0; k < All[code].nlay; k++)
		 if (All[code].layers[k] == (int) layer[l])
		   {
		     if (Mean != 2) PutGrib(ofileID, All[code].mean+k*DimGP, code, All[code].layers[k]);
		     if (Mean >= 2) PutGrib(ofileID, All[code].var +k*DimGP, code, All[code].layers[k]);
		   }
	   }
	   else {
	     if (All[code].leveltype == LEV_DOWN || All[code].leveltype == LEV_HEIGHT) {
               LevelType = All[code].leveltype;
               k         = All[code].layer;
	     }
	     else {
               LevelType = LEV_SURFACE;
               k         = 0;
	     }
	     if (Mean != 2) PutGrib(ofileID, All[code].mean, code, k);
	     if (Mean >= 2) PutGrib(ofileID, All[code].var , code, k);
	   }
         }
         else for (l = 0; l < nrql; l++) {
            switch (unitsel) {
            case 0: if (level[l] < LEV99_MARGIN) {
                       LevelType = LEV_99;
                       VarLev    = (int) level[l];
                    }  else {
                       LevelType = LEV_ISOBARIC;
                       VarLev    = (int) level[l]/100;
                    }
                    break;
            case 1:
            case 2: LevelType = LEV_HEIGHT;
                    VarLev    = (int) level[l];
                    break;
            default: fprintf(STDOUT,"Non allowed value for unitsel: %d\n",
                             unitsel);
                     exit(1);
            }
            if (Mean != 2) PutGrib(ofileID, All[code].mean+l*DimGP, code, VarLev);
            if (Mean >= 2) PutGrib(ofileID, All[code].var +l*DimGP, code, VarLev);
         }
      }
      FreeGrid();
      return;
   }

/* -------------------------------------------- */
/* Output of pressure level means and variances */
/* -------------------------------------------- */

   if (type == 30 && !Grib && Mean) {
      SrvHead[4] = Longitudes;
      SrvHead[5] = Latitudes;
      for (code = 0; code < CODES; code++)
      if (All[code].selected) {
         SrvHead[0] = code;
         wripmv(code, All[code].mean,All[code].var,DimGP,All[code].hlev);
      }
      FreeGrid();
      return;
   }

/* ================ */
/* Free mean fields */
/* ================ */

   if (type >= 40 && Mean)
   for (code = 0; code < CODES; code++)
   if (All[code].mean) {
      if (All[code].var ) All[code].var = FreeMem(All[code].var );
      if (All[code].grid) All[code].grid = FreeMem(All[code].grid);
      All[code].grid = All[code].mean;
      All[code].mean = NULL;
   }

/* =============================== */
/* Transformation to fourier space */
/* =============================== */

   if (type >= 40) {
      for (code = 0; code < CODES; code++)
      if (All[code].needed && (All[code].grid || All[code].swapped) &&
         (All[code].sfit   || type < 70)) {
         if (All[code].fourier == NULL) {
            FieldLength = DimFC * All[code].plev;
            All[code].fourier = DoubleAlloc(FieldLength,FieldName(code,"fourier"));
         }

         swapin (code,DimGP*All[code].plev,30);

         gp2fc(All[code].grid,All[code].fourier,
               Latitudes,Longitudes,All[code].plev,Fouriers);

         if (All[code].grid &&
            (All[code].sfit || type < 70))
             All[code].grid = FreeMem(All[code].grid);
      }
   }

   for (code = 0; code < CODES; code++)
   if (All[code].grid && (All[code].sfit || type < 70))
       All[code].grid = FreeMem(All[code].grid);

/* ------------------------ */
/* Output of fourier fields */
/* ------------------------ */

   if (type == 40) {
      SrvHead[4] = Latitudes ;
      SrvHead[5] = Fouriers  ;
      for (code = 0; code < CODES; code++)
      if (All[code].selected) {
         SrvHead[0] = code;
         wripre(code, All[code].fourier,DimFC,All[code].hlev);
      }
      FreeFourier();
      return;
   }

/* --------------------- */
/* Output of zonal means */
/* --------------------- */

   if (type == 41) {
      for (code = 0; code < CODES; code++)
      if (All[code].selected) {
         SrvHead[0] = code;
         wripzm(All[code].fourier,All[code].plev,Latitudes,DimFC);
      }
      FreeFourier();
      return;
   }

/* ================================ */
/* Transformation to spectral space */
/* ================================ */

   if (type >= 50) {

      if ( u_wind->needed && u_wind->fourier )
         scaluv(u_wind->fourier, coslat, Latitudes, Fouriers*nrql);
      if ( v_wind->needed && v_wind->fourier )
         scaluv(v_wind->fourier, coslat, Latitudes, Fouriers*nrql);

      for (code = 0; code < CODES; code++)
      if (All[code].needed && All[code].fourier) {
         if (All[code].spectral == NULL) {
            FieldLength = All[code].plev * DimSPOut;
            All[code].spectral = DoubleAlloc(FieldLength,
                                           FieldName(code,"spectral"));
         }
         if (TruncationOut < Truncation) {
            FieldLength = All[code].plev * DimSP;
            var_spectral = DoubleAlloc(FieldLength,"var_spectral");
            fc2sp(All[code].fourier,var_spectral,
                  pold,All[code].plev,Latitudes,Truncation);
            for (l = 0; l < All[code].plev; l++)
               ChangeTruncation(var_spectral+l*DimSP,Truncation,
                                All[code].spectral+l*DimSPOut,TruncationOut);
            var_spectral = FreeMem(var_spectral);
         }  else {
            fc2sp(All[code].fourier,All[code].spectral,
                  pold,All[code].plev,Latitudes,Truncation);
         }
      }

      if (Divergence->needed || Vorticity->needed ||
             velopot->needed ||    stream->needed) {
         if (Divergence->spectral == NULL)
             Divergence->spectral = DoubleAlloc(DimSPOut*nrql,
                                              "Divergence.spectral");
         if ( Vorticity->spectral == NULL)
              Vorticity->spectral = DoubleAlloc(DimSPOut*nrql,
                                              "Vorticity.spectral");
         if (TruncationOut < Truncation) {
            FieldLength = DimSP * nrql;
            div_spectral = DoubleAlloc(FieldLength,"div_spectral");
            vor_spectral = DoubleAlloc(FieldLength,"vor_spectral");
            uv2dv(u_wind->fourier,v_wind->fourier,
                  div_spectral,vor_spectral,
                  pol2,pol3,nrql,Latitudes,Truncation);
            for (l = 0; l < nrql; l++) {
               ChangeTruncation(div_spectral+l*DimSP,Truncation,
                                Divergence->spectral+l*DimSPOut,TruncationOut);
               ChangeTruncation(vor_spectral+l*DimSP,Truncation,
                                 Vorticity->spectral+l*DimSPOut,TruncationOut);
            }
            div_spectral = FreeMem(div_spectral);
            vor_spectral = FreeMem(vor_spectral);
         }  else {
            uv2dv(u_wind->fourier,v_wind->fourier,
                  Divergence->spectral,Vorticity->spectral,
                  pol2,pol3,nrql,Latitudes,Truncation);
         }
      }

      if (velopot->needed) {
         velopot->hlev = Divergence->hlev;
         velopot->plev = Divergence->plev;
         velopot->sfit = TRUE;
         if (velopot->spectral == NULL)
             velopot->spectral = DoubleAlloc(DimSPOut*nrql, "velopot.spectral");
         dv2ps(Divergence->spectral,velopot->spectral,nrql,TruncationOut);
      }

      if (stream->needed) {
         stream->hlev = Vorticity->hlev;
         stream->plev = Vorticity->plev;
         stream->sfit = TRUE;
         if ( stream->spectral == NULL)
              stream->spectral = DoubleAlloc(DimSPOut*nrql,
                                           "stream.spectral");
         dv2ps(Vorticity->spectral,stream->spectral,nrql,TruncationOut);
      }
   }

   for (code = 0; code < CODES; code++)
   if (All[code].fourier && (All[code].sfit || type < 61))
       All[code].fourier = FreeMem(All[code].fourier);

/* --------------------------------- */
/* Output of spectral fields in GRIB */
/* --------------------------------- */

   Representation = REP_SPECTRAL;
   if (type == 50 && Grib) {
      for (code = 0; code < CODES; code++)
      if (All[code].selected && All[code].spectral) {
         if (All[code].hlev == 1 || All[code].nlay) {
	   if (All[code].nlay) {
	     LevelType = All[code].leveltype;
	     for (l = 0; l < Layers; l++)
	       for (k = 0; k < All[code].nlay; k++)
		 if (All[code].layers[k] == (int) layer[l])
		   PutGrib(ofileID, All[code].spectral+k*DimSPOut, code, All[code].layers[k]);
	   }
	   else {
	     if (All[code].leveltype == LEV_DOWN || All[code].leveltype == LEV_HEIGHT) {
               LevelType = All[code].leveltype;
               k         = All[code].layer;
	     }
	     else {
               LevelType = LEV_SURFACE;
               k         = 0;
	     }
	     PutGrib(ofileID, All[code].spectral, code, k);
	   }
         }
         else for (l = 0; l < nrql; l++) {
            switch (unitsel) {
            case 0: if (level[l] < LEV99_MARGIN) {
                       LevelType = LEV_99;
                       VarLev    = (int) level[l];
                    }  else {
                       LevelType = LEV_ISOBARIC;
                       VarLev    = (int) level[l]/100;
                    }
                    break;
            case 1:
            case 2: LevelType = LEV_HEIGHT;
                    VarLev    = (int) level[l];
                    break;
            default: fprintf(STDOUT,"Non allowed value for unitsel: %d\n",
                             unitsel);
                     exit(1);
            }
            PutGrib(ofileID, All[code].spectral+l*DimSPOut, code, VarLev);
         }
      }
      FreeSpectral();
      return;
   }

/* ------------------------- */
/* Output of spectral fields */
/* ------------------------- */

   if (type == 50 && !Grib) {
      SrvHead[4] = DimSPOut;
      SrvHead[5] =     1;
      for (code = 0; code < CODES; code++)
      if (All[code].selected) {
         SrvHead[0] = code;
         wripre(code, All[code].spectral,DimSPOut,All[code].hlev);
      }
      FreeSpectral();
      return;
   }

/* ============================= */
/* Computations in fourier space */
/* ============================= */

   if (type >= 60) {
      for (code = 0; code < CODES; code++)
      if (All[code].needed && All[code].spectral) {
         if (All[code].fourier == NULL) {
            FieldLength = All[code].plev * DimFCOut;
            All[code].fourier = DoubleAlloc(FieldLength,
                                          FieldName(code,"fourier"));
         }
         sp2fc(All[code].spectral,All[code].fourier,poliOut,
               All[code].plev,LatitudesOut,FouriersOut,TruncationOut);
      }
      if ( u_wind->needed && u_wind->fourier )
         scaluv(u_wind->fourier, rclatOut, LatitudesOut, FouriersOut*nrql);
      if ( v_wind->needed && v_wind->fourier )
         scaluv(v_wind->fourier, rclatOut, LatitudesOut, FouriersOut*nrql);
   }

   FreeSpectral();

/* ------------------------ */
/* Output of fourier fields */
/* ------------------------ */

   if (type == 60) {
      SrvHead[4] = LatitudesOut;
      SrvHead[5] = FouriersOut ;
      for (code = 0; code < CODES; code++)
      if (All[code].selected) {
         SrvHead[0] = code;
         wripre(code, All[code].fourier,DimFCOut,All[code].hlev);
      }
      FreeFourier();
      return;
   }

/* --------------------- */
/* Output of zonal means */
/* --------------------- */

   if (type == 61) {
      for (code = 0; code < CODES; code++)
      if (All[code].selected) {
         SrvHead[0] = code;
         wripzm(All[code].fourier,All[code].plev,LatitudesOut,DimFCOut);
      }
      FreeFourier();
      return;
   }

/* ============================ */
/* Transformation to gridpoints */
/* ============================ */

   if (type >= 70) {
      for (code = 0; code < CODES; code++)
      if (All[code].needed && All[code].fourier) {
         FieldLength = All[code].plev * DimGPOut;
         if (All[code].grid == NULL) {
            All[code].grid = DoubleAlloc(FieldLength,
                                       FieldName(code,"grid"));
         }
         fc2gp(All[code].fourier,All[code].grid,
               LatitudesOut,LongitudesOut,All[code].plev,FouriersOut);

         if (Swap)
            swapout(code,FieldLength,30);
      }
   }

   FreeFourier();

/* -------------------------------------- */
/* Output of pressure level grids in GRIB */
/* -------------------------------------- */

   Representation = REP_GAUSS;
   if (type == 70 && Grib) {
      for (code = 0; code < CODES; code++)
      if (All[code].selected &&
         (All[code].swapped  || All[code].grid)) {
         if (All[code].swapped)
            swapin (code,DimGPOut*All[code].plev,30);

         if (All[code].hlev == 1 || All[code].nlay) {
	   if (All[code].nlay) {
	     LevelType = All[code].leveltype;
	     for (l = 0; l < Layers; l++)
	       for (k = 0; k < All[code].nlay; k++)
		 if (All[code].layers[k] == (int) layer[l])
		   PutGrib(ofileID, All[code].grid+k*DimGPOut, code, All[code].layers[k]);
	   }
	   else {
	     if (All[code].leveltype == LEV_DOWN || All[code].leveltype == LEV_HEIGHT) {
               LevelType = All[code].leveltype;
               k         = All[code].layer;
	     }
	     else {
               LevelType = LEV_SURFACE;
               k         = 0;
	     }
	     PutGrib(ofileID, All[code].grid, code, k);
	   }
         }
         else for (l = 0; l < nrql; l++) {
            switch (unitsel) {
            case 0: if (level[l] < LEV99_MARGIN) {
                       LevelType = LEV_99;
                       VarLev    = (int) level[l];
                    } else {
                       LevelType = LEV_ISOBARIC;
                       VarLev    = (int) level[l]/100;
                    }
                    break;
            case 1:
            case 2: LevelType = LEV_HEIGHT;
                    VarLev    = (int) level[l];
                    break;
            default: fprintf(STDOUT,"Non allowed value for unitsel: %d\n",
                             unitsel);
                     exit(1);
            }
            PutGrib(ofileID, All[code].grid+l*DimGPOut, code, VarLev);
         }
         if (All[code].swapped)
             All[code].grid = FreeMem(All[code].grid);
      }
      FreeGrid();
      return;
   }

/* ------------------------------ */
/* Output of pressure level grids */
/* ------------------------------ */

   if (type == 70 && !Grib) {
      SrvHead[4] = LongitudesOut;
      SrvHead[5] = LatitudesOut;
      for (code = 0; code < CODES; code++)
      if (All[code].selected &&
         (All[code].swapped  || All[code].grid)) {
         swapin (code,DimGPOut*All[code].plev,30);
         SrvHead[0] = code;
         wripre(code, All[code].grid,DimGPOut,All[code].plev);
         if (All[code].swapped)
             All[code].grid = FreeMem(All[code].grid);
      }
      FreeGrid();
      return;
   }
}

void PostProcess(void)
{
  if ( EndOfInterval )
    {
      if ( OutputInterval == DAILY_INTERVAL )
	fprintf (stderr, " Processed Day %2d   Month %2d   Year %2d",
		 OldDate.dy, OldDate.mo, OldDate.yr);
      else if ( OutputInterval == MONTHLY_INTERVAL )
	fprintf (stderr, " Processed Month %2d   Year %2d", OldDate.mo, OldDate.yr);
      else if ( OutputInterval == UNLIM_INTERVAL )
	fprintf (stderr, " Processed range %02d/%02d/%02d - %02d/%02d/%02d",
		 StartDate.dy, StartDate.mo, StartDate.yr,
		 OldDate.dy, OldDate.mo, OldDate.yr);

      if ( Mean ) fprintf(stderr, "  (Mean of %3d Terms)\n", MeanCount);
      else        fprintf(stderr, "   Terms %3d\n", MeanCount);

      EndOfInterval = FALSE;
      MeanCount  =     0;
    }
}

/* ================= */
/* switch input file */
/* ================= */

void SwitchFile()
{
  int echam4 = FALSE;
  int i, n;
  char y3, y2, y1, y0;
  char         m1, m0;
  char         d1, d0;
  char *rest;

  fclose(ifileID);
  i = strlen (ifile);
  rest = strrchr (ifile, '.');
  if ( rest )
    {
      if ( strlen(rest) == 3 )
	{
	  echam4 = TRUE;
	  y3 = ifile[i-9]; y2 = ifile[i-8];
	  y1 = ifile[i-7]; y0 = ifile[i-6];
	  m1 = ifile[i-5]; m0 = ifile[i-4];
	  d1 = ifile[i-2]; d0 = ifile[i-1];
	}
      else
	{
	  fprintf(stderr, " Not a valid filename: %s \n", ifile);
	  exit(1);
	}
    }
  else
    {
      y3 = ifile[i-6]; y2 = ifile[i-5];
      y1 = ifile[i-4]; y0 = ifile[i-3];
      m1 = ifile[i-2]; m0 = ifile[i-1];
      d1 = '0';        d0 = '1'   ;
    }

   for (n = 0; n < DayIn; n++) {
      if (d0 =='9') { d0 = '0'; d1++; }
      else d0++;
      if (d1 == '3' && d0 > '0') {
         d1 = '0'; d0 = '1';
         if (m1 == '0') {
            if (m0 == '9') { m0 = '0'; m1 = '1'; }
            else m0++;
         }  else {
            if (m0 < '2') m0++; else {
               m1 = '0';  m0 = '1';  y0++;
               if (y0 > '9') { y0 = '0'; y1++; }
               if (y1 > '9') {
                   y1 = (char) '0';
                   if (isdigit((int)y2)) y2++;
                   else           y2 = '1';
                   if (y2 > '9') {
                       y2 = (char) '0';
                       if (isdigit((int)y3)) y3++;
                       else           y3 = '1';
                   }
               }
            }
         }
      }
   }
   if (echam4) {
      ifile[i-9] = y3; ifile[i-8] = y2;
      ifile[i-7] = y1; ifile[i-6] = y0;
      ifile[i-5] = m1; ifile[i-4] = m0;
      ifile[i-2] = d1; ifile[i-1] = d0;
   }  else {
      ifile[i-6] = y3; ifile[i-5] = y2;
      ifile[i-4] = y1; ifile[i-3] = y0;
      ifile[i-2] = m1; ifile[i-1] = m0;
   }

   Multi--;
   fprintf (stderr, "Continuation File: %s\n", ifile);

   ifileID      = fopen(ifile, "rb");
}

void decodeField(double *field);

/*****************/
/* Echam Control */
/*****************/

void EchamControl(void)
{
  static char func[] = "EchamControl";
  int i;
  int righttime;
  int LevelOffset;
  int FieldLength;
  extern int PackComplex;

  while (1)
    {
      if (fileTypeIn == FT_GRIB)
	grib = (CHAR *)GribFromFile(ifileID);
      else
	Error(func, "unsupported file type : %d", fileTypeIn);

      if (GribErr == 1 && Multi > 0) {
         SwitchFile();
         if ( ifileID == NULL ) GribErr = 1;
         else {
            GribErr     = 0;
	    if ( fileTypeIn == FT_GRIB )
	      grib = (CHAR *)GribFromFile(ifileID);
	    else
	      Error(func, "unsupported file type : %d", fileTypeIn);
         }
      }

      if (GribErr == 1)
	{
	  SrvHead[2] = NewDate.YYMMDD;
	  SrvHead[3] = NewDate.HHMM;
	  if ( Mean && OutputInterval == MONTHLY_INTERVAL )
	   {
	     SrvHead[2] -= NewDate.dy;
	     SrvHead[3]  = -1;
	   }

	  EndOfInterval = TRUE;

	  EchamProcess();
	  PostProcess();
	  EchamDependencies();

	  return;
	}

      if (GribErr) GribAbort(GribErr);

      /* Skip records containing unneeded codes */
      /* (code = byte  9 of block 1 = grib[4+Grib1Offset+8])*/
      if (!All[grib[4+Grib1Offset+8]].needed) continue;

      /* Skip records with unselected hour */
      /* (hour = byte 16 of block 1 = grib[4+Grib1Offset+15])*/
      righttime = FALSE;
      for (i = 0; i < nrqh; i++)
	if (hours[i] < 0 || hours[i] == grib[4+Grib1Offset+15])
	  righttime = TRUE;
      if (!righttime) continue;

      decodeField(field);

      if (RepGrib == REP_SPECTRAL && PackComplex) {
         ScatterComplex(field, Truncation, DimSP);
         if (Debug) fprintf(STDOUT,"Complex packed spherical harmonics\n");
         ScaleComplex(field, Truncation);
      }

      if (Debug) {
         if (RepGrib == REP_SPECTRAL) fprintf(STDOUT,"T%d",Truncation);
         fprintf(STDOUT," Code %3d   Level%6d   %2.2d.%2.2d.%2.2d   %2.2d:%2.2d\n",
         CodeGrib,LevelGrib,NewDate.dy,NewDate.mo,NewDate.yr,
         NewDate.hr,NewDate.mn);
      }


      if ( OldDate.mo > 0 && memcmp(&NewDate,&OldDate,sizeof(struct Date)) != 0 )
	{
	  SrvHead[2] = OldDate.YYMMDD;
	  SrvHead[3] = OldDate.HHMM;
	  if ( Mean && OutputInterval == MONTHLY_INTERVAL )
	    {
	      SrvHead[2] -= OldDate.dy;
	      SrvHead[3]  = -1;
	    }

	  if      ( OutputInterval == DAILY_INTERVAL )
	    EndOfInterval = NewDate.dy != OldDate.dy;
	  else if ( OutputInterval == MONTHLY_INTERVAL )
	    EndOfInterval = NewDate.mo != OldDate.mo;
	  else if ( OutputInterval == UNLIM_INTERVAL )
	    EndOfInterval = 0;
	  else
	    Error(func, "output interval %d not implemented!\n", OutputInterval);

	  EchamProcess();
	  PostProcess();
	  EchamDependencies();
	}

      OldDate = NewDate;

      if (RepGrib == REP_SPECTRAL) {

      /* ---------------------------------------------------------- */
      /* Found spectral field ! If needed, allocate memory and copy */
      /* ---------------------------------------------------------- */

         if (All[CodeGrib].needed) {
            All[CodeGrib].sfit = TRUE;
            if (LevelGrib) {
               FieldLength = Dim3SP;
               LevelOffset = DimSP * (LevelGrib - 1);
               All[CodeGrib].surf = FALSE;
               All[CodeGrib].hlev = Levels;
               All[CodeGrib].plev = nrql;
            }  else {
               FieldLength =  DimSP;
               LevelOffset =      0;
               All[CodeGrib].surf = TRUE;
               All[CodeGrib].hlev = 1;
               All[CodeGrib].plev = 1;
               All[CodeGrib].leveltype = TypeOfLevel;
               All[CodeGrib].layer     = LayerGrib;
            }
            if (All[CodeGrib].spectral == NULL)
               All[CodeGrib].spectral = DoubleAlloc(FieldLength,
                                                  FieldName(CodeGrib,"spectral"));

            SrvHead[4] = DimSP;
            SrvHead[5] =     1;
            if (TruncationGrib != Truncation) {
               fprintf(STDOUT," *** Resolution Error ***\n");
               fprintf(STDOUT," Required %d - Found %d\n",
                      Truncation,TruncationGrib);
               exit(1);
            }
	    /*
            Eps = field[3];
            for (i = 1; i < DimSP; i++) field[i] -= Eps;
	    */
            RealCopy(All[CodeGrib].spectral+LevelOffset, field, DimSP);
         }
      }  else {
         if (All[CodeGrib].needed) {
            if (LevelGrib) {
               FieldLength = Dim3GP;
               LevelOffset = DimGP * (LevelGrib - 1);
               All[CodeGrib].surf = FALSE;
               All[CodeGrib].hlev = Levels;
               All[CodeGrib].plev = nrql;
               All[CodeGrib].sfit = TRUE;
	    }
	    else if ( All[CodeGrib].nlay ) {
	       int l;
               FieldLength = DimGP*All[CodeGrib].nlay;
	       for ( l = 0; l < All[CodeGrib].nlay; l++ )
		 if (All[CodeGrib].layers[l] == LayerGrib) break;
	       if ( l == All[CodeGrib].nlay )
		 Error(func, "Layer %d not found for code %d", LayerGrib, CodeGrib);
               LevelOffset = DimGP * l;
               All[CodeGrib].surf = TRUE;
               All[CodeGrib].hlev = All[CodeGrib].nlay;
               All[CodeGrib].plev = All[CodeGrib].nlay;
	       /*
               All[CodeGrib].plev = 1;
	       */
               All[CodeGrib].sfit = FALSE;
            }  else {
               FieldLength = DimGP;
               LevelOffset =     0;
               All[CodeGrib].surf = TRUE;
               All[CodeGrib].hlev = 1;
               All[CodeGrib].plev = 1;
               All[CodeGrib].sfit = FALSE;
               All[CodeGrib].leveltype = TypeOfLevel;
               All[CodeGrib].layer     = LayerGrib;
            }
            if (All[CodeGrib].hybrid == NULL)
                All[CodeGrib].hybrid = DoubleAlloc(FieldLength,
                                       FieldName(CodeGrib,"hybrid"));
            RealCopy(All[CodeGrib].hybrid+LevelOffset, field, DimGP);
            All[CodeGrib].swapped = FALSE;
         }
      }
   }
}


/* ================ */
/* Analysis Control */
/* ================ */

void AnalysisControl(void)
{
  static char func[] = "AnalysisControl";
   int i;
   int righttime;
   int LevelOffset;
   int FieldLength;
   extern int PackComplex;

   while (1) {
      if ( fileTypeIn == FT_GRIB )
	grib = (CHAR *)GribFromFile(ifileID);
      else
	Error(func, "unsupported file type : %d", fileTypeIn);

      if (GribErr == 1 && Multi > 0) {
         SwitchFile();
         if (ifileID == NULL) GribErr = 1;
         else {
            GribErr     = 0;
	    if ( fileTypeIn == FT_GRIB )
	      grib = (CHAR *)GribFromFile(ifileID);
	    else
	      Error(func, "unsupported file type : %d", fileTypeIn);
         }
      }

      if ( GribErr == 1 )
	{
	  SrvHead[2] = NewDate.YYMMDD;
	  SrvHead[3] = NewDate.HHMM;
	  if ( Mean && OutputInterval == MONTHLY_INTERVAL )
	    {
	      SrvHead[2] -= NewDate.dy;
	      SrvHead[3]  = -1;
	    }

	  EndOfInterval = TRUE;

	  AnalysisProcess();
	  PostProcess();
	  AnalysisDependencies();

	  return;
	}

     if (GribErr) GribAbort(GribErr);


      /* Skip records containing unneeded codes */
      /* (code = byte  9 of block 1 = grib[4+Grib1Offset+8])*/
      if (!All[grib[4+Grib1Offset+8]].needed) continue;

      /* Skip records with unselected hour */
      /* (hour = byte 16 of block 1 = grib[4+Grib1Offset+15])*/
      righttime = FALSE;
      for (i = 0; i < nrqh; i++)
      if (hours[i] < 0 || hours[i] == grib[4+Grib1Offset+15])
          righttime = TRUE;
      if (!righttime) continue;

      /* Skip records with unselected levels */
      /* (level = bytes 11+12 of block 1 = grib[4+Grib1Offset+10/11])*/
      LevelGrib   = Get2Byte(grib+4+Grib1Offset+10);

      LevelOffset = -1;
      if (grib[4+Grib1Offset+9] == LEV_99) {
         for (i = 0; i < nrql; ++i) {
             if ( level[i] == LevelGrib) {
                LevelOffset = i;
                break;
             }
         }
      }  else {
         for (i = 0; i < nrql; ++i) {
             if ( level[i] == LevelGrib * 100) {
                LevelOffset = i;
                break;
             }
         }
      }
      if (grib[4+Grib1Offset+9] == LEV_GROUND ||
          grib[4+Grib1Offset+9] == LEV_DOWN   ||
          grib[4+Grib1Offset+9] ==   0        ||
          grib[4+Grib1Offset+9] == LEV_SURFACE ) {
         LevelGrib   = 0;
         LevelOffset = 0;
      }

      if (LevelOffset < 0) continue;

      decodeField(field);

      if (RepGrib == REP_SPECTRAL && PackComplex) {
         ScatterComplex(field, Truncation, DimSP);
         if (Debug) fprintf(STDOUT,"Complex packed spherical harmonics\n");
         ScaleComplex(field, Truncation);
      }

/*    if (NewDate.YYMMDD < 820420 && CodeGrib == RHCODE) CodeGrib = SHCODE; */

      if (Debug) {
         if (RepGrib == REP_SPECTRAL) fprintf(STDOUT,"T%d",Truncation);
         fprintf(STDOUT," Code %3d   Level%6d   %2.2d.%2.2d.%2.2d   %2.2d:%2.2d\n",
         CodeGrib,LevelGrib,NewDate.dy,NewDate.mo,NewDate.yr,
              NewDate.hr,NewDate.mn);
      }

      if ( OldDate.mo > 0 && memcmp(&NewDate,&OldDate,sizeof(struct Date)) != 0 )
	{
	  SrvHead[2] = OldDate.YYMMDD;
	  SrvHead[3] = OldDate.HHMM;
	  if ( Mean && OutputInterval == MONTHLY_INTERVAL )
	    {
	      SrvHead[2] -= OldDate.dy;
	      SrvHead[3]  = -1;
	    }

	  if      ( OutputInterval == DAILY_INTERVAL )
	    EndOfInterval = NewDate.dy != OldDate.dy;
	  else if ( OutputInterval == MONTHLY_INTERVAL )
	    EndOfInterval = NewDate.mo != OldDate.mo;
	  else if ( OutputInterval == UNLIM_INTERVAL )
	    EndOfInterval = 0;
	  else
	    Error(func, "output interval %d not implemented!\n", OutputInterval);

	  AnalysisProcess();
	  PostProcess();
	  AnalysisDependencies();
	}

      OldDate = NewDate;

      if (All[CodeGrib].needed) {
         All[CodeGrib].sfit = TRUE;
         if (RepGrib == REP_SPECTRAL) {
            All[CodeGrib].hlev = nrql;
            All[CodeGrib].plev = nrql;
            SrvHead[4] = DimSPOut;
            SrvHead[5] =     1;
	    /*
            Eps = field[3];
            for (i = 1; i < DimSP; i++) field[i] -= Eps;
	    */
            if (LevelGrib) {
               All[CodeGrib].surf = FALSE;
               if (CodeGrib != UCODE && CodeGrib != VCODE) {
                  FieldLength = Dim3SPOut;
                  if (All[CodeGrib].spectral == NULL)
                      All[CodeGrib].spectral = DoubleAlloc(FieldLength,
                                               FieldName(CodeGrib,"spectral"));
                  ChangeTruncation(field,TruncationGrib,
                     All[CodeGrib].spectral+LevelOffset*DimSPOut,TruncationOut);
               }  else {
                  FieldLength = Dim3SP;
                  if (All[CodeGrib].spectral == NULL)
                      All[CodeGrib].spectral = DoubleAlloc(FieldLength,
                                               FieldName(CodeGrib,"spectral"));
                  RealCopy(All[CodeGrib].spectral+LevelOffset*DimSP,field,DimSP);
               }
            }
         }  else {
            SrvHead[4] = DimGPOut;
            SrvHead[5] =     1;
            if (LevelGrib) {
               FieldLength = Dim3GPOut;
               All[CodeGrib].hlev = nrql;
               All[CodeGrib].plev = nrql;
               All[CodeGrib].surf = FALSE;
               if (All[CodeGrib].grid == NULL)
                   All[CodeGrib].grid = DoubleAlloc(FieldLength,
                                        FieldName(CodeGrib,"grid"));
               RealCopy(All[CodeGrib].grid+LevelOffset*DimGP,field,DimGP);
            }  else {
               FieldLength = DimGPOut;
               All[CodeGrib].hlev = 1;
               All[CodeGrib].plev = 1;
               All[CodeGrib].surf = TRUE;
               All[CodeGrib].leveltype = TypeOfLevel;
               All[CodeGrib].layer     = LayerGrib;
               if (All[CodeGrib].grid == NULL)
                   All[CodeGrib].grid = DoubleAlloc(FieldLength,
                                        FieldName(CodeGrib,"grid"));
               RealCopy(All[CodeGrib].grid,field,DimGP);
            }
         }
      }
   }
}

char *amatch(char *msr, char *sub)
{
  int i,nm,ns;
  nm = strlen(msr);
  ns = strlen(sub);

  for (i = 0; i < nm-ns; i++)
    if (strncmp (msr+i,sub,ns) == 0) return (msr+i+ns);

  return NULL;
}

int scanpar(char *name, int def)
{
  char *cp;
  int value;

  cp = amatch(namelist,name);

  if ( cp == NULL ) value = def;
  else              value = atoi (cp);

  fprintf(stderr, "* %20.20s = %6d ", name, value);
  if ( value == def ) fprintf (stderr, " (default) *\n");
  else                fprintf (stderr, "           *\n");

  return (value);
}

void scantime(void)
{
  char *cp, *icp;
  int time;

  nrqh = 0;

   cp = amatch (namelist, "timesel");
   if (cp == NULL) {
      hours[nrqh++] = -1;
      fprintf (stderr, "* %20.20s = all               *\n","timesel");
      return;
   }

   time = (int) strtol (cp, &icp, 10);

   while ((char *)icp != (char *)cp && nrqh < MAX_HOURS) {
      hours[nrqh++] = time;
      cp = icp;
      time = (int) strtol (cp, &icp, 10);
   }
   fprintf (stderr, "* %20.20s = ", "timesel");
   for (time = 0 ; time < nrqh ; ++time) fprintf (stderr, " %02d", hours[time]);
   NewLine();
}

void PrintCodes(void)
{
   int code;

   fprintf (stderr, "* -------------------------------------- *\n");
   fprintf (stderr, "* Code Name                              *\n");
   fprintf (stderr, "* -------------------------------------- *\n");
   for (code=0 ; code < CODES ; ++code) if (CodeName[code]) {
      fprintf (stderr, "* %4d %-33.33s *\n", code, CodeName[code]);
   }
   fprintf(stderr, "* -------------------------------------- *\n");
}

void scancode(void)
{
   char *cp,*icp;
   int code;

   fprintf(STDOUT,"* ---------------------------------------- *\n");
   cp = amatch(namelist,"code");
   if (cp == NULL) code = -1;
   else code = (int) strtol(cp,&icp,10);
   if (code == -1) {
      fprintf(STDOUT,"* All detected codes selected              *\n");      
      for (code = 0; code < CODES; code++) {
         if (All[code].detected) {
            fprintf(STDOUT,"* Code %3d = %-29.29s *\n",code,CodeName[code]);
            All[code].selected = 1;
         }
      }
   }
   else {
      while (code > 0 && code < CODES) {
         fprintf(STDOUT,"* Code %3d = %-29.29s *\n",code,CodeName[code]);
         All[code].selected = 1;
         cp = icp;
         code = (int) strtol(cp,&icp,10);
      }
   }
}

void scanlevel(void)
{
  static char func[] = "scanlevel";
  char *cp,*icp;
  double lev;
  int n, k, l, Doppelt, NichtDa;
  int RemoveLevel[MAX_LEVELS];

  fprintf(STDOUT,"* ---------------------------------------- *\n");
  nrql = 0;
  cp = amatch(namelist, "level");
  if (cp == NULL) lev = -1;
  else            lev = strtod(cp,&icp);
  if ( lev == -1 )
    {
      fprintf(STDOUT,"* All detected levels selected             *\n");
      for (lev = 1; lev <= Levels; lev++)
	{
	  level[nrql++] = lev;
	  fprintf(STDOUT,"* Level %2d = %14.4f                *\n", nrql, lev);    
	}
    }
  else
    {
      level[nrql++] = lev;
      fprintf(STDOUT,"* Level %2d = %14.4f                *\n", nrql, lev);
      cp = icp;
      lev = strtod(cp, &icp);
      while ( lev > 0 && nrql < MAX_LEVELS )
	{
	  level[nrql++] = lev;
	  fprintf(STDOUT,"* Level %2d = %14.4f                *\n", nrql, lev);
	  cp = icp;
	  lev = strtod(cp, &icp);
	}
    }
  /*
  Stars(44); NewLine();
  */
  if ( AnalysisData || type < 30 )
    for (k = 0; k < nrql; k++)
      {
	LevelGrib = (int) level[k];
	Doppelt = -1;
	for (l = 0; l < nrql; l++)
	  if (LevelGrib == level[l]) Doppelt++;
	if (Doppelt)
	  {
	    fprintf(STDOUT,
		    "* Level %2d = %14.4f double request *\n",k+1,level[k]);
	  }
	NichtDa = TRUE;
	for (l = 0; l < nfol; l++)
	  if (LevelGrib == LevelFound[l]) NichtDa = FALSE;
	if (NichtDa)
	  {
	    fprintf(STDOUT,
		    "* Level %2d = %14.4f not in input   *\n",k+1,level[k]);
	  }
	if (Doppelt || NichtDa)
	  RemoveLevel[k] = TRUE;
	else
	  RemoveLevel[k] = FALSE;
      }
  else
    for (k = 0; k < nrql; k++)
      RemoveLevel[k] = FALSE;

  n = nrql;
  for (k = 0; k < n; k++)
    {
      LevelGrib = (int) level[k];
      if (RemoveLevel[k] == TRUE)
	{
	  for (l = k; l < nrql; l++)
	    if (LevelGrib == level[l])
	      RemoveLevel[l] = FALSE;
	  nrql--;
	  for (l = k; l < nrql; l++)
	    {
	      level[l] = level[l+1];
	      RemoveLevel[l] = RemoveLevel[l+1];
	    }
	}
    }

  if ( n != nrql )
    {
      Stars(44); NewLine();
      if ( labort )
	Error(func, "Inconsistent or invalid level list !");
      else
	Warning(func, "Inconsistent or invalid level list !");
    }
  /*
  NewLine();
  */
}

void scanlayer(void)
{
  static char func[] = "scanlayer";
  char *cp,*icp;
  double lay;
  int n, k, l, Doppelt, NichtDa;
  int RemoveLayer[MAX_LEVELS];

  fprintf(STDOUT,"* ---------------------------------------- *\n");
  nrqlay = 0;
  cp = amatch(namelist, "layer");
  if (cp == NULL) lay = -1;
  else            lay = strtod(cp,&icp);
  if ( lay == -1 )
    {
      fprintf(STDOUT,"* All detected layers selected             *\n");
      for (l = 0; l < Layers; l++)
	{
	  layer[nrqlay++] = LayerFound[l];
	  fprintf(STDOUT,"* Layer %2d = %14.4f                *\n",
		  nrqlay, (double)LayerFound[l]);    
	}
    }
  else
    {
      layer[nrqlay++] = lay;
      fprintf(STDOUT,"* Layer %2d = %14.4f                *\n", nrqlay, lay);
      cp = icp;
      lay = strtod(cp, &icp);
      while ( lay > 0 && nrqlay < MAX_LEVELS )
	{
	  layer[nrqlay++] = lay;
	  fprintf(STDOUT,"* Layer %2d = %14.4f                *\n", nrqlay, lay);
	  cp = icp;
	  lay = strtod(cp, &icp);
	}
    }
  Stars(44); NewLine();

  if ( AnalysisData || type < 30 )
    for (k = 0; k < nrqlay; k++)
      {
	LevelGrib = (int) layer[k];
	Doppelt = -1;
	for (l = 0; l < nrqlay; l++)
	  if (LevelGrib == layer[l]) Doppelt++;
	if (Doppelt)
	  {
	    fprintf(STDOUT,
		    "* Layer %2d = %14.4f double request *\n",k+1,layer[k]);
	  }
	NichtDa = TRUE;
	for (l = 0; l < nglay; l++)
	  if (LevelGrib == LayerFound[l]) NichtDa = FALSE;
	if (NichtDa)
	  {
	    fprintf(STDOUT,
		    "* Layer %2d = %14.4f not in input   *\n",k+1,layer[k]);
	  }
	if (Doppelt || NichtDa)
	  RemoveLayer[k] = TRUE;
	else
	  RemoveLayer[k] = FALSE;
      }
  else
    for (k = 0; k < nrqlay; k++)
      RemoveLayer[k] = FALSE;

  n = nrqlay;
  for (k = 0; k < n; k++)
    {
      LevelGrib = (int) layer[k];
      if (RemoveLayer[k] == TRUE)
	{
	  for (l = k; l < nrqlay; l++)
	    if (LevelGrib == layer[l])
	      RemoveLayer[l] = FALSE;
	  nrqlay--;
	  for (l = k; l < nrqlay; l++)
	    {
	      layer[l] = layer[l+1];
	      RemoveLayer[l] = RemoveLayer[l+1];
	    }
	}
    }

  if ( n != nrqlay )
    {
      Stars(44); NewLine();
      if ( labort )
	Error(func, "Inconsistent or invalid layer list !");
      else
	Warning(func, "Inconsistent or invalid layer list !");
    }

  NewLine();
}

void nameini(void)
{
   int code;

   for (code = 0; code < CODES; code++) CodeName[code] = "--";

   CodeName[129] = "Geopotential";
   CodeName[130] = "Temperature";
   CodeName[131] = "u-velocity";
   CodeName[132] = "v-velocity";
   CodeName[133] = "Specific Humidity";
   CodeName[134] = "Surface Pressure";
   CodeName[138] = "Vorticity";
   CodeName[139] = "Surface Temperature";
   CodeName[140] = "Soil Wetness";
   CodeName[141] = "Snow Depth";
   CodeName[142] = "Large Scale Precipitation";
   CodeName[143] = "Convective Precipitation";
   CodeName[144] = "Snow Fall";
   CodeName[145] = "Boundary Layer Dissipation";
   CodeName[146] = "Surface Sensible Heat Flux";
   CodeName[147] = "Surface Latent Heat Flux";
   CodeName[148] = "Streamfunction";
   CodeName[149] = "Velocity Potential";
   CodeName[151] = "Mean Sea Level Pressure";
   CodeName[152] = "Log Surface Pressure";
   CodeName[153] = "Liquid Water Content";
   CodeName[155] = "Divergence";
   CodeName[156] = "Geopotential Height";
   CodeName[157] = "Relative Humidity";
   CodeName[158] = "Tendency of Surface Pressure";
   CodeName[159] = "U*^3";
   CodeName[160] = "Surface Runoff";
   CodeName[161] = "Liquid Water Content";
   CodeName[162] = "Cloud Cover";
   CodeName[163] = "Total Cloud Cover";
   CodeName[164] = "Total Cloud Cover (mean)";
   CodeName[165] = "10m u-velocity";
   CodeName[166] = "10m v-velocity";
   CodeName[167] = "2m Temperature";
   CodeName[168] = "2m Dew Point Temperature";
   CodeName[222] = "Liquid Water Content (acc.)";
   CodeName[223] = "Cloud Cover (acc.)";
}

void swapini()
{
   int i;
   char *tmpdir;

   tmpdir = getenv("TMPDIR");
   if (tmpdir == NULL) Swap = FALSE;
   else {              Swap = TRUE ;
      for (i = 0; i < CODES; i++) {
          swapname[i] = CharAlloc(L_tmpnam, "swapname[i]");
/*        sprintf(swapname[i],"%s%s%3.3d",tmpdir,"/swp",i); */
          tmpnam(swapname[i]);
      }
   }
}

void parini(void)
{
   int i, c;
   long length;
   extern double Grav;
   extern double RD;

   fseek(stdin,0L,SEEK_END);
   length = ftell(stdin);
   if (length == 0L) {
      fprintf(STDOUT,"\n stdin not connected\n");
      Usage();
   }
   fseek(stdin,0L,SEEK_SET);
   i = 1;
   namelist[0] = ' ';
   c = getchar();
   while ((c != EOF) && i < 1023) {
           if ((c >= '0' && c <= '9') ||
               (c == '-' || c == '.'))  namelist[i++] = c;
      else if  (c >= 'a' && c <= 'z')   namelist[i++] = c;
      else if  (c >= 'A' && c <= 'Z')   namelist[i++] = tolower(c);
      else c = ' ';

      if (c == ' ' && namelist[i-1] != ' ') namelist[i++] = c;
      c = getchar();
   }
   namelist[i] = 0;
   if (Debug) {
      NewLine();
      Stars(77);
      NewLine();
      fprintf(STDOUT,"* Length of namelist:%4d bytes                                             *\n",
             (int) strlen(namelist));

      for (i = 0; i<strlen(namelist); i+=60)
      fprintf(STDOUT,"* namelist[%02d]=%-60.60s *\n",i,namelist+i);
      Stars(77);
      NewLine();
   }

   NewLine();
   Stars(44);
   NewLine();
   type          = scanpar("type" ,0);
   Multi         = scanpar("multi",0);
   Mean          = scanpar("mean" ,1);
   OutputInterval= scanpar("interval",MONTHLY_INTERVAL);
   Grib          = scanpar("grib" ,0);
   SrvHead[6]    = scanpar("head7",0);
   unitsel       = scanpar("unitsel",0);
   DayIn         = scanpar("dayinc" ,30);
   TruncationOut = scanpar("res",0);
   Latitudes     = scanpar("lat",0);
   mars          = scanpar("mars",0);

   if (Multi) --Multi;

   if (mars) {
      Grav         = MARS_GRAV;
      PlanetRadius = MARS_RADIUS;
      RD           = MARS_RD;
   }

   scantime();
   scancode();
   scanlevel();
   scanlayer();
}

void dimcalc(void)
{
  if ( AnalysisData ) Levels = nrql;
  if ( TruncationOut == 0 ) TruncationOut = Truncation;
  if ( TruncationOut > Truncation )
    fprintf(STDOUT,"*** Warning: Output Resolution > Data Resolution ! ***\n");

  if ( Latitudes != 0 )         /* Latitudes specified by namelist */
    { 
      Longitudes    = Latitudes * 2;
      LatitudesOut  = Latitudes;
      LongitudesOut = LatitudesOut * 2;

    }
  else                         /* Latitudes derived from truncation */
    {
      if ( LatitudesIn )
	{
	  Latitudes = LatitudesIn;
	  LatitudesOut  = Latitudes;
	}
      else
	{
	  Latitudes = 2 * ((Truncation*3 + 3) / 4);

	  if ( Truncation == 30 ) Latitudes = 48;

	  LatitudesOut = 2 * ((TruncationOut*3 + 3) / 4);

	  if ( TruncationOut == 30 ) LatitudesOut = 48;
	}
      Longitudes = Latitudes * 2;
      if ( Truncation == 62 ) Longitudes = 192;
      LongitudesOut = LatitudesOut * 2;
      if ( TruncationOut == 62 ) LongitudesOut = 192;
    }

   Waves         = Truncation + 1;
   Fouriers      = Waves * 2;
   DimSP         = (Truncation + 1) * (Truncation + 2);
   DimFC         = Latitudes * Fouriers;
   DimGP         = Latitudes * Longitudes;
   Dim3GP        = Levels * DimGP;
   Dim3FC        = Levels * DimFC;
   Dim3SP        = Levels * DimSP;
   DimSP_half    = DimSP / 2;
   HalfLevels    = Levels + 1;
   WavesOut      = TruncationOut + 1;
   FouriersOut   = WavesOut * 2;
   DimSPOut      = (TruncationOut + 1) * (TruncationOut + 2);
   DimFCOut      = LatitudesOut * FouriersOut;
   DimGPOut      = LatitudesOut * LongitudesOut;
   Dim3GPOut     = Levels * DimGPOut;
   Dim3FCOut     = Levels * DimFCOut;
   Dim3SPOut     = Levels * DimSPOut;
   DimSPOut_half = DimSPOut / 2;

   fprintf(STDOUT," ****************************\n");
   if (AnalysisData)
   fprintf(STDOUT," * Found Ana or Re-Ana Data *\n");
   fprintf(STDOUT," * Truncation        = %4d *\n", Truncation);
   fprintf(STDOUT," * Levels            = %4d *\n", Levels);
   fprintf(STDOUT," * Layers            = %4d *\n", Layers);
   fprintf(STDOUT," * Latitudes         = %4d *\n", Latitudes);
   fprintf(STDOUT," * Longitudes        = %4d *\n", Longitudes);
   if ( TruncationOut != Truncation )
     {
      fprintf(STDOUT," * Output Truncation = %4d *\n", TruncationOut);
      fprintf(STDOUT," * Output Latitudes  = %4d *\n", LatitudesOut);
      fprintf(STDOUT," * Output Longitudes = %4d *\n", LongitudesOut);
     }
   fprintf(STDOUT," ****************************\n\n");
}

void setVCT(int nvct, double *psec2)
{
  int i;

  if ( vct == NULL ) vct = DoubleAlloc(nvct, "vct");

  for (i = 0; i < nvct; i++) vct[i] = psec2[10+i];

  if ( Debug ) for (i = 0; i < nvct/2; i++)
    fprintf(stderr," vct: %10.4f %10.4f\n", vct[i], vct[i+nvct/2]);
}

/* ----------------------------------------------------------- */
/* Extract basic dimension information from GRIB header fields */
/* ----------------------------------------------------------- */

void precntl(void)
{
  static char func[] = "precntl";
  int lats, SchonDa, l;
  int levelType = 0;
  int TimIndex = 9, b, SumOfTime = 0, OldSumOfTime;
  int Level109 = 0;
  int code = 0, level = 0, gridtype = 0;
  int code0 = -1, level0 = -1;
  int irec;

  while ( Truncation == 0 || Levels == 0 )
    {
      if ( fileTypeIn == FT_GRIB )
	grib = (CHAR *)GribFromFile(ifileID);
      else
	Error(func, "unsupported file type : %d", fileTypeIn);

      if ( GribErr == 1 ) break;
      if ( GribErr ) GribAbort(GribErr);

      gribexdp(sec0, sec1, sec2, psec2, sec3,
	       psec3, sec4, field, klenp, (void *) grib,
	       kleng, &kword, "J", &kret);

      code      = GribParameter(sec1);
      level     = GribLevel1(sec1);
      levelType = GribLevelType(sec1);
      gridtype  = GribGridType(sec2);

      if (Debug) fprintf(STDOUT,"Code %3d  Level = %d LevelType %3d Type %3d",
			 code, level, levelType, gridtype);

      if ( code0 == -1 && level0 == -1 )
	{
	  code0  = code;
	  level0 = level;
	}
      else
	{
	  if ( code0 == code && level0 == level ) break;
	}

      SumOfTime = 0;
      for (b = TimIndex; b < TimIndex+12; b++)
	SumOfTime += sec1[b];

      if ( GribNumVCP(sec2) > 0 && nvct == 0 )
	{
	  nvct = GribNumVCP(sec2);
	  Level109 = (nvct >> 1) - 1;
	  if ( nvct > 512-10 ) Error(func, "array psec2 to small. nvct = %d\n", nvct);
	  setVCT(nvct, psec2);
	}

      RepGrib = GribGridType(sec2);
    
      /*    Then, check for horizontal resolution/truncation */
      if      ( RepGrib == REP_SPECTRAL )
	{
	  Truncation   = GribPentaJ(sec2);
	}
      else if ( RepGrib == REP_GAUSS )
	{
	  lats         = GribNumLat(sec2);
	  LatitudesIn  = lats;
	  switch (lats) {
	  case 512: Truncation = 511; break;
	  case 320: 
	    if ( GribCenterID(sec1) == 98 && GribModelID(sec1) < 183 )
	      Truncation = 213;
	    else         
	      Truncation = 319;
	    break;
	  case 160: Truncation = 106; break;
	  case 128: Truncation =  85; break;
	  case  96: Truncation =  63; break;
	  case  94: Truncation =  62; break;
	  case  64: Truncation =  42; break;
	  case  48:
	    if ( GribCenterID(sec1) == 98 && GribSubCenterID(sec1) == 232 )
	      Truncation =  31;
	    else
	      Truncation =  30;
	    break;
	  case  32: Truncation =  21; break;
	  default :
	    fprintf(STDOUT,"%d Gaussian latitudes not supported.\n", lats);
	    fprintf(STDOUT,"Please contact schulzweida@dkrz.de or");
	    fprintf(STDOUT,"                 kornblueh@dkrz.de\n");
	    continue;
	  }
	}
      else
	Abort("precntl: Unknown Grid Representation Type ***\n");

      levelType = GribLevelType(sec1);
      if      (levelType == LEV_HYBRID)
	{
	  Levels = Level109;
	  AnaLevelFactor = 1;
	}
    else if (levelType == LEV_ISOBARIC)
      {
	Levels = 17;
	AnaLevelFactor = 100;
      }
    else if (levelType == LEV_GROUND ||
	     levelType == LEV_SURFACE)
      {
	if (AnalysisData)
	  {
            Levels = 1;
            AnaLevelFactor = 100;
	  }
	else
	  {
            Levels = 0;
            AnaLevelFactor = 0;
	  }
      }

    /*    Re-Analysis, Local extension ECMWF */
    if (levelType == LEV_SURFACE)
      if (GribECMWFClass(sec1) == 3 && GribECMWFLocalExtention(sec1)) AnalysisData = TRUE;

    /*    Leveltype 112 (below land surface) is */
    /*    unique to Re-Analysis Surface Data    */
    if (levelType == LEV_GROUND )
      if (GribECMWFClass(sec1) == 3 && GribECMWFLocalExtention(sec1)) AnalysisData = TRUE;

    /*    Another indicator for Analysis */
    if (levelType == LEV_ISOBARIC ||
	levelType == LEV_99)                 AnalysisData = TRUE;

    if (GribModelID(sec1) ==  21)  /* NCAR base model */
      {
	if (Debug) fprintf(STDOUT,"  NCAR Model ");
	AnalysisData = FALSE;
	Truncation   = 42;
	Levels       = 18;
      }
    if ( Debug ) fprintf(STDOUT, "   T=%3d   L=%2d\n", Truncation, Levels);
  }

  rewind(ifileID);
  GribErr = 0;

  OldSumOfTime = SumOfTime;

  if (Debug) fprintf(STDOUT," CODE CHECK\n");

  irec = 0;
  while ( TRUE )
    {
      if ( fileTypeIn == FT_GRIB )
	grib = (CHAR *)GribFromFile(ifileID);
      else
	Error(func, "unsupported file type : %d", fileTypeIn);

      if ( GribErr == 1 ) break;

      gribexdp(sec0, sec1, sec2, psec2, sec3,
	       psec3, sec4, field, klenp, (void *) grib,
	       kleng, &kword, "J", &kret);

      RepGrib = GribGridType(sec2);

      levelType = GribLevelType(sec1);

      if ( SumOfTime != OldSumOfTime ) break;
      /*
      if ( irec > 0 )
	{
	  if ( code0 == code && level0 == level ) break;
	}
	*/

      irec++;

      if ( RepGrib == REP_SPECTRAL ) Spectral = TRUE;
      else                           Gaussian = TRUE;

      CodeGrib  = GribParameter(sec1);
      LevelGrib = GribLevel1(sec1);
      if (levelType == LEV_SURFACE ||
	  levelType == LEV_GROUND  ||
	  levelType == 0)            LevelGrib = 0;
      if (levelType == LEV_ISOBARIC) LevelGrib *= 100;

      if ( grib[TimIndex+8] == 2 )
	{
	  /*    Products occupying a range of time are excluded */
	  SumOfTime = OldSumOfTime;
	}
      else
	{
	  SumOfTime = 0;
	  for (b = TimIndex; b < TimIndex+12; b++)
	    SumOfTime += sec1[b];
	}

      if ( levelType != LEV_DOWN )
	{
	  if ( nfol == 0 )
	    {
	      LevelFound[nfol++] = LevelGrib;
	    }
	  else
	    {
	      SchonDa = FALSE;
	      for (l = 0; l < nfol; l++)
		{
		  if ( LevelFound[l] == LevelGrib ) SchonDa = TRUE;
		}
	      if ( !SchonDa ) LevelFound[nfol++] = LevelGrib;
	    }
	}

      All[CodeGrib].detected  = TRUE;
      All[CodeGrib].leveltype = levelType;

      if ( levelType == LEV_DOWN )
	{
	  if ( nglay == 0 )
	    {
	      LayerFound[nglay++] = LevelGrib;
	      All[CodeGrib].layers[All[CodeGrib].nlay++] = LevelGrib;
	    }
	  else
	    {
	      SchonDa = FALSE;
	      for (l = 0; l < nglay; l++)
		{
		  if ( LayerFound[l] == LevelGrib ) SchonDa = TRUE;
		}
	      if ( !SchonDa ) LayerFound[nglay++] = LevelGrib;
	      SchonDa = FALSE;
	      for (l = 0; l < All[CodeGrib].nlay; l++)
		{
		  if (All[CodeGrib].layers[l] == LevelGrib) SchonDa = TRUE;
		}
	      if ( !SchonDa ) All[CodeGrib].layers[All[CodeGrib].nlay++] = LevelGrib;
	    }
	}
    }

  Layers = nglay;

  for (code = 0; code < CODES; code++)
    {
      if ( Debug )
	if ( All[code].detected )
	  fprintf(STDOUT," Detected Code %3d with %3d layer ( leveltype = %d )\n",
		  code, All[code].nlay, All[code].leveltype);
    }
}


/* =============================================== */
/* statistics - appends info about memory usage    */
/*              and time consumption to FILE  *sf  */
/*              and echos info to the calling user */
/* =============================================== */

void statistics(char *procpath, double CPUTime)
{
  FILE *sf;
  double MaxMBytes;
  time_t tp;
  long  yy, mm, dd, hh, mi;
  char mtype[12];
  char *proc;
  char *name;
  char  stat_file[128];
  extern int MaxMem;

  (void) time(&tp);
  yy    = gmtime(&tp)->tm_year + 1900;
  mm    = gmtime(&tp)->tm_mon + 1;
  dd    = gmtime(&tp)->tm_mday   ;
  hh    = gmtime(&tp)->tm_hour   ;
  mi    = gmtime(&tp)->tm_min    ;
  name  = getpwuid(getuid())->pw_name;

  proc = strrchr(procpath,'/');
  if (proc == 0) proc = procpath;
  else           proc++         ;

  strcpy(stat_file, "/mf/m/m214003/usage/after");

  if ( Debug ) fprintf(STDOUT,"%s\n",stat_file);

  /* Experience */
  MaxMBytes = 2.16 + 6.27 * (((double) DimGP - 2048.) / 51200.) +
             (double)MaxMem / 1048576.;

  sf = fopen(stat_file,"a");
  if ( sf )
    {
      char *hostname;

      fprintf(sf, "%.7s %4.4ld.%2.2ld.%2.2ld %2.2ld:%2.2ld %4.2f ",
                  name,   yy,   mm,   dd,   hh,   mi,   VERSION);
      if ( (hostname = getenv("HOST")) == NULL )
	fprintf(sf, "%-9.9s %7.1f %7.1f T%3.3d %s\n",
		proc, MaxMBytes, CPUTime, Truncation, "unknown");
      else
        fprintf(sf, "%-9.9s %7.1f %7.1f T%3.3d %s\n",
		proc, MaxMBytes, CPUTime, Truncation, hostname);
    }

#if defined (CRAY)
#  if defined (_CRAYMPP)
     strcpy(mtype, " CRAYMPP --");
#  elif (_MAXVL == 64)
     strcpy(mtype, " CRAYVL64 -");
#  elif (_MAXVL == 128)
     strcpy(mtype, " CRAYVL128 ");
#  else
     strcpy(mtype, " CRAY -----");
#  endif
#elif defined (SX)
     strcpy(mtype, " NECSX ----");
#elif defined (__uxp__)
     strcpy(mtype, " FUJI -----");
#elif defined (sun)
     strcpy(mtype, " SUN ------");
#elif defined (i386)
     strcpy(mtype, " i386 -----");
#elif defined (sgi)
     strcpy(mtype, " sgi ------");
#else
     strcpy(mtype, "-----------");
#endif

  fprintf(STDOUT, "   NORMAL EXIT\n");
  fprintf(STDOUT, " ------   End    after  -%-11.11s- %7.1f sec", mtype, CPUTime);
  fprintf(STDOUT, " --- %7.1f MB ---\n", MaxMBytes);
}

void Usage(void)
{
  fprintf(STDOUT, "\nafterburner [options] <InputFile> <OutputFile>\n");
  fprintf(STDOUT, "     option -c : print available codes and names\n");
  fprintf(STDOUT, "     option -v : debug mode (verbose output)\n");
  fprintf(STDOUT, "     option -h : help (this output)\n");
  fprintf(STDOUT, "     option -a : forces analysis data process\n");
  fprintf(STDOUT, "     option -s : reduces memory at expense of I/O");
  fprintf(STDOUT, " and therefore CPU time\n");
  fprintf(STDOUT, "   <InputFile> : ECHAM or ECMWF Ana or ReAna file\n");
  fprintf(STDOUT, "  <OutputFile> : GRIB or SERVICE format file\n");
  fprintf(STDOUT, "  namelist is read from <stdin>\n");
  fprintf(STDOUT, "  output is written to <stdout>\n\n");

  fprintf(STDOUT, "  dedault Namelist: \n");
  fprintf(STDOUT, "  &SELECT\n");
  fprintf(STDOUT, "    TYPE = 0, CODE = -1, LEVEL = -1, MULTI = 0, DAYIN = 30,\n");
  fprintf(STDOUT, "    MEAN = 0, GRIB = 0, HEAD7 = 0, TIMESEL = -1, UNITSEL = 0 \n");
  fprintf(STDOUT, "  &END\n");

  exit(1);
}

int
read_block (FILE *fp, char block[])
{
  static char func[] = "read_block";  
  size_t nread;

  nread = fread(block, 1, 4, fp);
  if (nread != 4) return (1);

  if (Debug) Message(func, ">%c%c%c%c<", block[0], block[1], block[2], block[3]);

  return (0);
}

int
getFileType (const char *filename)
{
  static char func[] = "getFileType";
  FILE *fp;
  int ierr;
  int filetype = 0;
  unsigned long blocklen = 0;
  char block[4];

  fp = fopen (filename, "r");

  if ( fp )
    {
      ierr = read_block(fp, block);

      if (ierr)
	Error(func, "could not read 4 bytes from %s", filename);

      if (strncmp (block, "GRIB", 4) == 0)
	{
	  filetype = FT_GRIB;
	  if ( Debug ) Message(func, "found GRIB file = %s", filename);
	}
      else if (strncmp (block, "CDF", 3) == 0)
	{
	  filetype = FT_CDF;
	  if ( Debug ) Message(func, "found CDF file = %s", filename);
	}
      else if (strncmp (block, "PUMA", 4) == 0)
	{
	  if ( Debug ) Message(func, "found PUMA file = %s", filename);
	  Error(func, "PUMA files unsupported");
	}
      else
	{
	  int recpos;

	  rewind (fp);
	  ierr = gribFileSeek(fp, &recpos);

	  if (!ierr)
	    {
	      filetype = FT_GRIB;
	      if (Debug) Message(func, "found seek GRIB file = %s", filename);
	    }

	  if ( !filetype ) {		  
	    Error(func, "unimplemented filetype\n"
		  "                file  : %s\n"
		  "                block : %c%c%c%c length : %d\n",
		  filename, block[0], block[1], block[2], block[3],
		  (int) blocklen);
	  }
	}
    }
  else
    {
      Error(func, "open failed on %s", filename);
    }

  if ( fclose(fp) )
    Error(func, "close failed on %s", filename);

  return (filetype);
}

int main(int argc, char *argv[])
{
  static char func[] = "main";
  int   code, i;
  char *proc = argv[0];
  char  Line[132];
  char  legfile[80];
  char  legfileOut[80];
  double CPUTime;
  extern int DebugMem;
  struct stat filestat;
  unsigned long bufferSize;

  (void) clock();

  /*********************/
  /* print information */
  /*********************/

  fprintf(STDOUT,"******************************************\n");
  fprintf(STDOUT,"* afterburner %4.2f (%s)         *\n", VERSION, VERDATE);
  fprintf(STDOUT,"* ECHAM & analyses postprocessor         *\n");
  fprintf(STDOUT,"* -------------------------------------- *\n");
  fprintf(STDOUT,"* Uwe Schulzweida           MPI  Hamburg *\n");
  fprintf(STDOUT,"* Arno Hellbach            DKRZ  Hamburg *\n");
  fprintf(STDOUT,"* Edilbert Kirk       University Hamburg *\n");
  fprintf(STDOUT,"* Michael Ponater   DLR Oberpfaffenhofen *\n");
  fprintf(STDOUT,"* -------------------------------------- *\n");
  fprintf(STDOUT,"* call with option -c for list of codes  *\n");
  fprintf(STDOUT,"******************************************\n");

  if ( sizeof(double) != 8 || sizeof(int) < 4 )
    {
      fprintf(STDOUT, "byte size of type double %d\n", (int) sizeof(double));
      fprintf(STDOUT, "byte size of type int %d\n", (int) sizeof(int));
      fprintf(STDOUT, "byte size of type size_t %d\n", (int) sizeof(size_t));
      return(1);
    }

  ifileID = fopen("/pf/m/m214003/doc/afterburner.doc","r");
  if (ifileID)
    {
      do
	{
	  fgets(Line, 130, ifileID);
	  fprintf(STDOUT, " %s", Line);
	}
      while (!feof(ifileID) && Line[0] == '#');
      fclose(ifileID);
    }

  nameini();

  /***********************/
  /* options & filenames */
  /***********************/

  for (i = 1 ; i < argc ; ++i)
    {
      if ( argv[i][0] == '-' )
	{
	  if      (argv[i][1] == 'c') PrintCodes();
	  else if (argv[i][1] == 'v') Debug = 1;
	  else if (argv[i][1] == 's') Swap = TRUE;
	  else if (argv[i][1] == 'a') AnalysisData = 1;
	  else if (argv[i][1] == 'g') ECMWFMars = TRUE;
	  else if (argv[i][1] == 'w') labort = FALSE;
	  else if (argv[i][1] == 'b')
	    {
	      ++i;
	      Message(func, "\noption -b not longer needed!\n\n");
	    }
	  else Usage();
	}
      else if (ifile[0] == '\0') strcpy(ifile,argv[i]);
      else if (ofile[0] == '\0') strcpy(ofile,argv[i]);
      else if (strcmp("Debug",argv[i]) == 0) Debug = 1;
      else Usage();
    }

  if ( ifile[0] == '\0' || ofile[0] == '\0' )
    {
      fprintf(STDOUT, "*** Missing filename ***\n");
      Usage();
    }
  fprintf(STDOUT, "*  Input File: %-25s *\n", ifile);
  fprintf(STDOUT, "* Output File: %-25s *\n", ofile);
  if ( Debug )
    {
      DebugMem = 1;
      fprintf(STDOUT,"* Debug on!                              *\n");
    }
  if ( Swap ) swapini();
  if ( Swap ) fprintf(STDOUT, "* Swapping on!                           *\n");
  fprintf(STDOUT, "******************************************\n");

  fileTypeIn = getFileType(ifile);

  /***********************/
  /* open in/output file */
  /* get length of input */
  /***********************/

  if ( stat(ifile, &filestat) != 0 )
    SysError(func, ifile);
  /*
  bufferSize = filestat.st_blksize;

  if ( bufferSize < MINFILEBUFFERSIZE ) bufferSize = MINFILEBUFFERSIZE;
  */
  bufferSize = MINFILEBUFFERSIZE;

  ifileID = fopen(ifile, "rb");
  if ( ifileID == 0 )
    SysError(func, ifile);

  filebuffer = (char *) malloc(bufferSize);
  if ( filebuffer == 0 )
    SysError(func, "allocation of file buffer failed");	

  if ( setvbuf(ifileID, filebuffer, _IOFBF, bufferSize) )
    SysError(func, "setvbuf failed");

  if ( Debug )
    Message(func, "file buffer size %lu", bufferSize);

  ofileID = fopen(ofile, "w");
  if ( ofileID == 0 )
    SysError(func, ofile);

  /******************/
  /* pre-processing */
  /******************/

  rewind(ifileID);
  GribErr = 0;

  precntl();

  rewind(ifileID);
  GribErr = 0;

  /*******************/
  /* initializations */
  /*******************/

  parini();

  if ( Swap && Mean )
    {
      fprintf(STDOUT," ::::::::::::::::::::::::::::::::::::::::::::::\n");
      fprintf(STDOUT," -> Don't use Swap with Mean\n");
      fprintf(STDOUT," -> Premature Exit. Sorry.\n");
      exit(1);
    }

  dimcalc();

  if ( type < 50 && AnalysisData )
    {
      fprintf(STDOUT," ::::::::::::::::::::::::::::::::::::::::::::::\n");
      fprintf(STDOUT," -> Type < 50 is not appropriate for Analysis.\n");
      fprintf(STDOUT," -> Please check wether you can use type >= 50.\n");
      fprintf(STDOUT," -> Premature Exit. Sorry.\n");
      exit(1);
    }

  if ( type < 20 && Grib == 1 )
    {
      fprintf(STDOUT," ::::::::::::::::::::::::::::::::::::::::::::::\n");
      fprintf(STDOUT," -> GRIB output for Types < 20 is redundant !\n");
      fprintf(STDOUT," -> Please restart afterburner with GRIB = 0.\n");
      fprintf(STDOUT," -> Premature Exit. Sorry.\n");
      exit(1);
    }

  if ( TruncationOut < Truncation )
    {
      if ( type < 50 )
	{
	  fprintf(STDOUT," ::::::::::::::::::::::::::::::::::::::::::::::\n");
	  fprintf(STDOUT," -> Parameter res:\n");
	  fprintf(STDOUT," -> Please avoid type < 50 for technical reasons.\n");
	  fprintf(STDOUT," -> Premature Exit. Sorry.\n");
	  exit(1);
	}
      if ( Gaussian )
	{
	  fprintf(STDOUT," ::::::::::::::::::::::::::::::::::::::::::::::\n");
	  fprintf(STDOUT," -> Parameter res:\n");
	  fprintf(STDOUT," -> Please avoid 'downmelting' with gaussian\n");
	  fprintf(STDOUT," -> grids. No source code yet !\n");
	  fprintf(STDOUT," -> Premature Exit. Sorry.\n");
	  exit(1);
	}
    }

  if ( TruncationOut > Truncation )
    {
      fprintf(STDOUT," ::::::::::::::::::::::::::::::::::::::::::::::\n");
      fprintf(STDOUT," -> Parameter res\n");
      fprintf(STDOUT," -> Please avoid OutputTruncation > Truncation!\n");
      fprintf(STDOUT," -> No source code yet !\n");
      fprintf(STDOUT," -> Premature Exit. Sorry.\n");
      exit(1);
    }
  if ( unitsel == 2 )
    for (i = 0; i < nrql; i++) level[i] = level[i] * 1000;

  if ( !AnalysisData )
    for (i = 0; i < nrql; i++)
      {
	if ( (level[i] >= 65535) && unitsel )
	  {
	    fprintf(STDOUT,"\n Level %9.2f out of range !\n",level[i]);
	    exit(1);
	  }
	if ( (level[i] == 0) && !unitsel && type >= 20 && nrql > 1 )
	  {
	    fprintf(STDOUT,"\n Level %9.2f illegal for type %d\n",level[i],type);
	    exit(1);
	  }
      }

  filename = strrchr(ifile,'/');
  if (filename == 0) filename = ifile;
  else               filename++ ;
  SrvHead[7] = atoi(filename);
  if (type >= 30 && type < 50 &&
      (Divergence->selected || velopot->selected ||
        Vorticity->selected ||  stream->selected ||
       AnalysisData))
    {
      if (type == 30) type = 70;
      if (type == 40) type = 60;
      if (type == 41) type = 61;
      fprintf(STDOUT,"\n TYPE changed to %d\n",type);
    }
  if ( AnalysisData ) AnalysisDependencies();
  else
    {
      EchamDependencies();
      Geopotential->needed |= type >= 30
                           || (SLP->needed && !SLP->detected) || 
	             (GeopotHeight->needed && !GeopotHeight->detected);
    }
  if ( type > 0 )
    {
      sprintf(legfile,"/pool/GRIB/legini/legini.t%03d",Truncation);
      legpol = fopen(legfile,"rb");
      legini();
      if ( legpol ) fclose(legpol);
    }
  if ( type > 50 )
    {
      sprintf(legfileOut,"/pool/GRIB/legini/legini.t%03d",TruncationOut);
      legpolOut = fopen(legfileOut,"rb");
      leginiOut();
      if (legpolOut) fclose(legpolOut);
    }
  if ( u_wind->needed || v_wind->needed )
    {
      dv2uv_f1 = DoubleAlloc(DimSPOut_half,"dv2uv_f1");
      dv2uv_f2 = DoubleAlloc(DimSPOut_half,"dv2uv_f2");
      geninx(TruncationOut, dv2uv_f1, dv2uv_f2);
    }
  if ( AnalysisData ) AnalysisControl();
  else                EchamControl();

  if ( ofileID ) fclose(ofileID);

  for (code = 0; code < CODES; code++)
    remove(swapname[code]);

  CPUTime = (double) clock() / (double) CLOCKS_PER_SEC;
  statistics(proc,CPUTime);

  return(0);
}

#ifndef _MALLOC_H
#define NOMALLINFO
#endif

int    DebugMem = 0 ;
int    MaxMem = 0;

/* ============================================== */
/* CharAlloc - Allocate space for character array */
/* ============================================== */

char *CharAlloc(size_t bytes, char *array_name)
{
#ifndef NOMALLINFO
  struct mallinfo info;
#endif
  char *result = NULL;
  if ( bytes > 0 ) {
    result = (char *)calloc(bytes, sizeof(char));
#ifndef NOMALLINFO
    info   = mallinfo();
    if (info.uordblks > MaxMem) MaxMem = info.uordblks;
    if ( DebugMem )
      fprintf(stderr,"<   CharAlloc: at%9p %8lu char   for %-28.28s Mem:%10d >\n",
	      (void *) result, (unsigned long) bytes, array_name, info.uordblks);
#else
    if ( DebugMem )
      fprintf(stderr,"<   CharAlloc: at%9p %8lu char   for %-28.28s Mem: - >\n",
	      (void *) result, (unsigned long) bytes, array_name);
#endif
    if ( result == NULL ) Abort(" *** Out of Memory ***\n");
  }
  return (result);
}

/* ======================================= */
/* IntAlloc - Allocate space for int array */
/* ======================================= */

int *IntAlloc(int words,char *array_name)
{
#ifndef NOMALLINFO
   struct mallinfo info;
#endif
   int *result = NULL;
   if (words > 0) {
      result = (int *) malloc(words * sizeof(int));
#ifndef NOMALLINFO
      info   = mallinfo();
      if (info.uordblks > MaxMem) MaxMem = info.uordblks;
      if (DebugMem)
         fprintf(stderr,"<    IntAlloc: at%9p %8d int    for %-28.28s Mem:%10d >\n",
               (void *) result,words,array_name,info.uordblks);
#else
      if (DebugMem)
         fprintf(stderr,"<    IntAlloc: at%9p %8d int    for %-28.28s Mem: - >\n",
               (void *) result,words,array_name);
#endif
      if (result == NULL) Abort(" *** Out of Memory ***\n");
   }
   return(result);
}

/* =========================================== */
/* FloatAlloc - Allocate space for float array */
/* =========================================== */

float *FloatAlloc(int words, char *array_name)
{
#ifndef NOMALLINFO
   struct mallinfo info;
#endif
   float *result = NULL;
   if (words > 0) {
      result = (float *) malloc(words * sizeof(float));
#ifndef NOMALLINFO
      info   = mallinfo();
      if (info.uordblks > MaxMem) MaxMem = info.uordblks;
      if (DebugMem)
         fprintf(stderr,"<  FloatAlloc: at%9p %8d float  for %-28.28s Mem:%10d >\n",
               (void *) result,words,array_name,info.uordblks);
#else
      if (DebugMem)
         fprintf(stderr,"<  FloatAlloc: at%9p %8d float  for %-28.28s Mem: - >\n",
               (void *) result,words,array_name);
#endif
      if (result == NULL) Abort(" *** Out of Memory ***\n");
   }
   return(result);
}

/* ============================================= */
/* DoubleAlloc - Allocate space for double array */
/* ============================================= */

double *DoubleAlloc(int words, char *array_name)
{
#ifndef NOMALLINFO
   struct mallinfo info;
#endif
   double *result = NULL;
   if (words > 0) {
      result = (double *) malloc(words * sizeof(double));
#ifndef NOMALLINFO
      info   = mallinfo();
      if (info.uordblks > MaxMem) MaxMem = info.uordblks;
      if (DebugMem)
         fprintf(stderr,"< DoubleAlloc: at%9p %8d double for %-28.28s Mem:%10d >\n",
               (void *) result,words,array_name,info.uordblks);
#else
      if (DebugMem)
         fprintf(stderr,"< DoubleAlloc: at%9p %8d double for %-28.28s Mem: - >\n",
               (void *) result,words,array_name);
#endif
      if (result == NULL) Abort(" *** Out of Memory ***\n");
   }
   return(result);
}

/* ================ */
/* Free array space */
/* ================ */

void *FreeMem(void *ptr)
{
#ifndef NOMALLINFO
   struct mallinfo info;
#endif
   free(ptr);
#ifndef NOMALLINFO
   info   = mallinfo();
   if (DebugMem)
      fprintf(stderr,"<        free: at%9p                             "
                     "                     Mem:%10d >\n",
             (void *) ptr,info.uordblks);
#else
   if (DebugMem)
      fprintf(stderr,"<        free: at%9p                             "
                     "                     Mem: - >\n",
             (void *) ptr);
#endif
   return NULL;
}

/* ======================================== */
/* Convert Spectral Array to new resolution */
/* ======================================== */

void ChangeTruncation(double *SPin, int Resin, double *SPout, int Resout)
{
   int n, m;

   if (Resout <= Resin) {
      for (n = 0; n <= Resout; n++) {
         for (m = n; m <= Resout; m++) {
            *SPout++ = *SPin++ ;
            *SPout++ = *SPin++ ;
         }
         SPin += 2 * (Resin-Resout);
      }
   }  else {
      for (n = 0; n <= Resin; n++) {
         for (m = n; m <= Resin; m++) {
            *SPout++ = *SPin++ ;
            *SPout++ = *SPin++ ;
         }
         for (m = Resin+1; m <= Resout; ++m) {
            *SPout++ = 0.0;
            *SPout++ = 0.0;
         }
      }
      for (n = Resin+1; n <= Resout; ++n) {
         for (m = n; m <= Resout; ++m) {
            *SPout++ = 0.0;
            *SPout++ = 0.0;
         }
      }
   }
}

#if defined (SX)
#pragma odir switch,-ev
#endif
void sp2fc(double *sa,double *fa,double *poli,int klev,int nlat,int nfc,int nt)
{
   int lev,jmm,jfc,lat;
   double sar,sai;
   double *Far,*fai,*pol;

   RealZero (fa, klev*nlat*nfc);
   for (lev = 0; lev < klev; lev++) {
      pol = poli;
      for (jmm = 0; jmm <= nt; jmm++) {
         for (jfc = jmm; jfc <= nt; jfc++) {
            sar = *sa++     ;
            sai = *sa++     ;
            Far = fa        ;
            fai = fa + nlat ;
            for (lat = 0; lat < nlat; lat++) {
               *Far += *pol * sar;
               *fai += *pol * sai;
               Far++;
               fai++;
               pol++;
            }
         }
         fa += 2 * nlat;
      }
   }
}

#if defined (SX)
#pragma odir switch,-dv
#endif
void fc2sp(double *fa, double *sa, double *poli, int klev, int nlat, int nt)
{
   int lev,jmm,jfc,lat;
   double sar,sai,*Far,*fai,*pol;

   for (lev = 0; lev < klev; lev++) {
      pol = poli;
      for (jmm = 0; jmm <= nt; jmm++) {
         for (jfc = jmm; jfc <= nt; jfc++) {
            Far = fa        ;
            fai = fa + nlat ;
            sar = 0.0       ;
            sai = 0.0       ;
            for (lat = 0; lat < nlat; lat++) {
               sar += *pol * *Far;
               sai += *pol * *fai;
               Far++;
               fai++;
               pol++;
            }
            *sa++ = sar;
            *sa++ = sai;
         }
         fa += 2 * nlat;
      }
   }
}


#define QUA 0.25
#define QT5 0.559016994374947

#define S36 0.587785252292473
#define S60 0.866025403784437
#define S72 0.951056516295154

#define SQ2 0.707106781186547524401

#define D60 (S60+S60)

#define FORK for(k=la;k<=kstop;k+=la){
#define LOOP for(l=0;l<la;++l){i=ibase;j=jbase;for(ijk=0;ijk<lot;++ijk){
#define ENDL i+=inc3;j+=inc4;}ibase+=inc1;jbase+=inc2;}

void fft_set(int ifax[], int n)
{
   int j,k,nfax;
   nfax = 0;
   for (k = 0; k < 9; ++k) ifax[k] = 0;
   ifax[9] = n;
   if    (n % 8 == 0) {ifax[++nfax] = 8; n /= 8;}
   while (n % 6 == 0) {ifax[++nfax] = 6; n /= 6;}
   while (n % 5 == 0) {ifax[++nfax] = 5; n /= 5;}
   while (n % 4 == 0) {ifax[++nfax] = 4; n /= 4;}
   while (n % 3 == 0) {ifax[++nfax] = 3; n /= 3;}
   if    (n % 2 == 0) {ifax[++nfax] = 2; n /= 2;}
   ifax[0] = nfax;
#if defined (CRAY)
#pragma _CRI novector
#endif
#if defined (SX)
#pragma vdir novector
#endif
#if defined (__uxp__)
#pragma loop scalar
#endif
   for (k = 0; k < nfax/2; ++k) {
      j = ifax[k+1];
      ifax[k+1] = ifax[nfax-k];
      ifax[nfax-k] = j;
   }
}


#if defined (SX)
#pragma odir switch,-ev
#endif
int rpassc(double *a, double *b, double *c, double *d,
           int inc1, int inc2, int inc3, int inc4,
           int lot , int n   , int ifac, int la  )
{
/*
   rpassc' - performs one pass through data as part;
   of multiple real fft (fourier synthesis) routine;

   a is first real input vector
   b is equivalent to a + la * inc1
   c is first real output vector
   d is equivalent to c + ifac * la * inc2
   inc1 is the addressing increment for a;
   inc2 is the addressing increment for c;
   inc3 is the increment between input vectors a;
   inc4 is the increment between output vectors c;
   lot is the number of vectors;
   n is the length of the vectors;
   ifac is the current factor of n;
   la is the product of previous factors;
   ierr is an error indicator:;
   0 - pass completed without error;
   2 - ifac not catered for;
   3 - ifac only catered for if la=n/ifac;
*/

  int i0,i1,i2,i3,i4,i5,i6,i7;
  int j0,j1,j2,j3,j4,j5,j6,j7;
  int ia,ib,ic,id,ie,iF;
  int ja,jb,jc,jd,je,jf;
  int i,j,k,ijk,l,m;
  int ibase,jbase;
  int iink,jink;
  int jump;
  int kstop;

  double c1,c2,c3,c4,c5;
  double s1,s2,s3,s4,s5;
  double kpidn;
  double angle;
  double qqrt5;
  double ssin36;
  double ssin72;
  double pin;

  double a10,a11,a20,a21;
  double b10,b11,b20,b21;

  m     = n  / ifac;
  iink  = la * inc1;
  jink  = la * inc2;
  jump  = (ifac-1) * jink;
  kstop = (n-ifac) / (2*ifac);
  pin   = 2.0 * M_PI / n;
  ibase = 0;
  jbase = 0;

  switch (ifac) {
    case 2: {
      double a0m1,b0p1;

      i0 = j0 = 0;
      i1 = i0 + inc1 * (m+m-la);
      j1 = j0 + jink;
      if (la != m) {
        LOOP
          c[j0+j] = a[i0+i] + a[i1+i];
          c[j1+j] = a[i0+i] - a[i1+i];
        ENDL
        i0    += iink;
        iink  += iink;
        i1    -= iink;
        ibase  = 0;
        jbase += jump;
        jump  += jump + jink;

        if (i0 != i1) {
          FORK
            angle = k * pin;
            c1 = cos(angle);
            s1 = sin(angle);
            ibase = 0;
            LOOP
              a0m1 = a[i0+i] - a[i1+i];
              b0p1 = b[i0+i] + b[i1+i];

              c[j0+j] = a[i0+i] + a[i1+i];
              d[j0+j] = b[i0+i] - b[i1+i];
              c[j1+j] = c1 * a0m1 - s1 * b0p1;
              d[j1+j] = s1 * a0m1 + c1 * b0p1;
            ENDL
            i0    += iink;
            i1    -= iink;
            jbase += jump;
          } /* End FORK */
          if (i0 > i1) return 0;
        } /* End (i0 != i1) */
        ibase = 0;
        LOOP
          c[j0+j] =  a[i0+i];
          c[j1+j] = -b[i0+i];
        ENDL
           }
           else /* (la != m) */ {
        LOOP
          c[j0+j] = 2.0 * (a[i0+i] + a[i1+i]);
          c[j1+j] = 2.0 * (a[i0+i] - a[i1+i]);
        ENDL
      }
      return 0;
    }

    case 3: {
      double afa1,a1p2,a1m2,a0mm,a0mp;
      double bfa1,b1p2,b1m2,b0mm,b0mp;

      i0 = j0 = 0 ;
      i1 = i0 + inc1 * (m+m-la);
      i2 = i1;
      j1 = j0 + jink;
      j2 = j1 + jink;

      if (la != m) {
        LOOP
          afa1 = a[i0+i] - 0.5 * a[i1+i];
          bfa1 =           S60 * b[i1+i];

          c[j0+j] = a[i0+i] + a[i1+i];
          c[j1+j] = afa1 - bfa1;
          c[j2+j] = afa1 + bfa1;
        ENDL
        i0    += iink;
        iink  += iink;
        i1    += iink;
        i2    -= iink;
        jbase += jump;
        jump  += jump + jink;

        if (i0 != i2) {
          FORK
            angle = k * pin;
            c1 = cos(angle);
            s1 = sin(angle);
            angle += angle;
            c2 = cos(angle);
            s2 = sin(angle);
            ibase = 0;
            LOOP
              a1p2 = a[i0+i] - 0.5 * (a[i1+i] + a[i2+i]);
              b1m2 = b[i0+i] - 0.5 * (b[i1+i] - b[i2+i]);
              a1m2 =           S60 * (a[i1+i] - a[i2+i]);
              b1p2 =           S60 * (b[i1+i] + b[i2+i]);

              a0mm = a1p2 - b1p2;
              a0mp = a1p2 + b1p2;
              b0mm = b1m2 - a1m2;
              b0mp = b1m2 + a1m2;

              c[j0+j] = a[i0+i] + a[i1+i] + a[i2+i];
              d[j0+j] = b[i0+i] + b[i1+i] - b[i2+i];
              c[j1+j] = c1 * a0mm - s1 * b0mp;
              d[j1+j] = s1 * a0mm + c1 * b0mp;
              c[j2+j] = c2 * a0mp - s2 * b0mm;
              d[j2+j] = s2 * a0mp + c2 * b0mm;
            ENDL
            i0    += iink;
            i1    += iink;
            i2    -= iink;
            jbase += jump;
          } /* End FORK */
          if (i0 > i2) return 0;
        } /* End (i0 != i2) */
        ibase=0;
        LOOP
          a0mp = 0.5 * a[i0+i];
          b0mp = S60 * b[i0+i];

          c[j0+j] = a[i0+i] + a[i1+i];
          c[j1+j] = a0mp - a[i1+i] - b0mp;
          c[j2+j] = a[i1+i] - a0mp - b0mp;
        ENDL
      }
      else /* (la != m) */ {
        LOOP
          a0mp = 2.0 * a[i0+i] - a[i1+i];
          b0mp = D60 * b[i1+i];

          c[j0+j] = 2.0 * (a[i0+i] + a[i1+i]);
          c[j1+j] = a0mp - b0mp;
          c[j2+j] = a0mp + b0mp;
        ENDL
      }
      return 0;
    }

    case 4: {
      double a0m1,a0p2,a1p3,a0m2,a1m3,a0p2ma1p3,a0m2pb1p3,a0m2mb1p3;
      double b0p1,b0p2,b1p3,b0m2,b1m3,b0p2pa1m3,b0p2ma1m3,b0m2mb1m3;

      i0 = j0 = 0;
      i1 = i3 = i0 + inc1 * (m+m-la);
      i2 = i1 + inc1 * (m+m);
      j1 = j0 + jink;
      j2 = j1 + jink;
      j3 = j2 + jink;

      if (la != m) {
        LOOP
          a0p2 = a[i0+i] + a[i2+i];
          a0m2 = a[i0+i] - a[i2+i];

          c[j0+j] = a0p2 + a[i1+i];
          c[j1+j] = a0m2 - b[i1+i];
          c[j2+j] = a0p2 - a[i1+i];
          c[j3+j] = a0m2 + b[i1+i];
        ENDL
        i0    += iink;
        iink  += iink;
        i1    += iink;
        i2    -= iink;
        i3    -= iink;
        jbase += jump;
        jump  += jump + jink;

        if (i1 != i2) {
          FORK
            angle = kpidn = k * pin;
            c1 = cos(angle);
            s1 = sin(angle);
            angle += kpidn;
            c2 = cos(angle);
            s2 = sin(angle);
            angle += kpidn;
            c3 = cos(angle);
            s3 = sin(angle);
            ibase=0;
            LOOP
              a0p2 = a[i0+i] + a[i2+i];
              a0m2 = a[i0+i] - a[i2+i];
              a1p3 = a[i1+i] + a[i3+i];
              a1m3 = a[i1+i] - a[i3+i];
              b0p2 = b[i0+i] + b[i2+i];
              b0m2 = b[i0+i] - b[i2+i];
              b1p3 = b[i1+i] + b[i3+i];
              b1m3 = b[i1+i] - b[i3+i];

              a0p2ma1p3 = a0p2 - a1p3;
              a0m2pb1p3 = a0m2 + b1p3;
              a0m2mb1p3 = a0m2 - b1p3;
              b0p2pa1m3 = b0p2 + a1m3;
              b0p2ma1m3 = b0p2 - a1m3;
              b0m2mb1m3 = b0m2 - b1m3;

              c[j0+j] = a0p2 + a1p3;
              d[j0+j] = b0m2 + b1m3;
              c[j2+j] = c2 * a0p2ma1p3 - s2 * b0m2mb1m3;
              d[j2+j] = s2 * a0p2ma1p3 + c2 * b0m2mb1m3;
              c[j1+j] = c1 * a0m2mb1p3 - s1 * b0p2pa1m3;
              d[j1+j] = s1 * a0m2mb1p3 + c1 * b0p2pa1m3;
              c[j3+j] = c3 * a0m2pb1p3 - s3 * b0p2ma1m3;
              d[j3+j] = s3 * a0m2pb1p3 + c3 * b0p2ma1m3;
            ENDL
            i0    += iink;
            i1    += iink;
            i2    -= iink;
            i3    -= iink;
            jbase += jump;
          } /* End FORK */
          if (i1 > i2) return 0;
        } /* End (i1 != i2) */
        ibase=0;
        LOOP
          a0m1 = a[i0+i] - a[i1+i];
          b0p1 = b[i0+i] + b[i1+i];

          c[j0+j] =  a[i0+i] + a[i1+i];
          c[j2+j] =  b[i1+i] - b[i0+i];

          c[j1+j] =  SQ2 * (a0m1 - b0p1);
          c[j3+j] = -SQ2 * (a0m1 + b0p1);
        ENDL
      }
      else /* (la != m) */ {
        LOOP
          a0p2 = a[i0+i] + a[i2+i];
          a0m2 = a[i0+i] - a[i2+i];

          c[j0+j] = 2.0 * (a0p2 + a[i1+i]);
          c[j1+j] = 2.0 * (a0m2 - b[i1+i]);
          c[j2+j] = 2.0 * (a0p2 - a[i1+i]);
          c[j3+j] = 2.0 * (a0m2 + b[i1+i]);
        ENDL
      }
      return 0;
    }

    case 5: {
      double a1p2,a1m2,a0mm,a0mp,b136,b172,b236,b272;

      i0 = j0 = 0;
      i1 = i4 = i0 + inc1 * (m+m-la);
      i2 = i3 = i1 + inc1 * (m+m);
      j1 = j0 + jink;
      j2 = j1 + jink;
      j3 = j2 + jink;
      j4 = j3 + jink;

      if (la != m) {
        LOOP
          a1p2 = QUA * (a[i1+i] + a[i2+i]);
          a1m2 = QT5 * (a[i1+i] - a[i2+i]);

          a0mp = a[i0+i] - a1p2 + a1m2;
          a0mm = a[i0+i] - a1p2 - a1m2;

          b136 = b[i1+i] * S36;
          b172 = b[i1+i] * S72;
          b236 = b[i2+i] * S36;
          b272 = b[i2+i] * S72;

          c[j0+j] = a[i0+i] + a[i1+i] + a[i2+i];
          c[j1+j] = a0mp - b172 - b236;
          c[j2+j] = a0mm - b136 + b272;
          c[j3+j] = a0mm + b136 - b272;
          c[j4+j] = a0mp + b172 + b236;
        ENDL
        i0    += iink;
        iink  += iink;
        i1    += iink;
        i2    += iink;
        i3    -= iink;
        i4    -= iink;
        jbase += jump;
        jump  += jump + jink;

        if (i1 != i3) {
          FORK
            angle = kpidn = k * pin;
            c1 = cos(angle);
            s1 = sin(angle);
            angle += kpidn;
            c2 = cos(angle);
            s2 = sin(angle);
            angle += kpidn;
            c3 = cos(angle);
            s3 = sin(angle);
            angle += kpidn;
            c4 = cos(angle);
            s4 = sin(angle);
            ibase=0;
            LOOP
              a10=(a[i0+i]-0.25*((a[i1+i]+a[i4+i])+(a[i2+i]+a[i3+i])))
                 +QT5*((a[i1+i]+a[i4+i])-(a[i2+i]+a[i3+i]));
              a20=(a[i0+i]-0.25*((a[i1+i]+a[i4+i])+(a[i2+i]+a[i3+i])))
                 -QT5*((a[i1+i]+a[i4+i])-(a[i2+i]+a[i3+i]));
              b10=(b[i0+i]-0.25*((b[i1+i]-b[i4+i])+(b[i2+i]-b[i3+i])))
                 +QT5*((b[i1+i]-b[i4+i])-(b[i2+i]-b[i3+i]));
              b20=(b[i0+i]-0.25*((b[i1+i]-b[i4+i])+(b[i2+i]-b[i3+i])))
                 -QT5*((b[i1+i]-b[i4+i])-(b[i2+i]-b[i3+i]));

              a11=S72*(b[i1+i]+b[i4+i])+S36*(b[i2+i]+b[i3+i]);
              a21=S36*(b[i1+i]+b[i4+i])-S72*(b[i2+i]+b[i3+i]);
              b11=S72*(a[i1+i]-a[i4+i])+S36*(a[i2+i]-a[i3+i]);
              b21=S36*(a[i1+i]-a[i4+i])-S72*(a[i2+i]-a[i3+i]);

              c[j0+j]=a[i0+i]+((a[i1+i]+a[i4+i])+(a[i2+i]+a[i3+i]));
              d[j0+j]=b[i0+i]+((b[i1+i]-b[i4+i])+(b[i2+i]-b[i3+i]));
              c[j1+j]=c1*(a10-a11)-s1*(b10+b11);
              d[j1+j]=s1*(a10-a11)+c1*(b10+b11);
              c[j4+j]=c4*(a10+a11)-s4*(b10-b11);
              d[j4+j]=s4*(a10+a11)+c4*(b10-b11);
              c[j2+j]=c2*(a20-a21)-s2*(b20+b21);
              d[j2+j]=s2*(a20-a21)+c2*(b20+b21);
              c[j3+j]=c3*(a20+a21)-s3*(b20-b21);
              d[j3+j]=s3*(a20+a21)+c3*(b20-b21);
            ENDL
            i0    += iink;
            i1    += iink;
            i2    += iink;
            i3    -= iink;
            i4    -= iink;
            jbase += jump;
          } /* End FORK */
          if (i1 > i3) return 0;
        } /* End (i1 != i3) */
        ibase=0;
        LOOP
          c[j0+j] = a[i0+i] + a[i1+i] + a[i2+i];
          c[j1+j] = (QT5  * (a[i0+i]-a[i1+i])
             + (0.25  * (a[i0+i]+a[i1+i])-a[i2+i]))
             - (S36 *  b[i0+i]+S72*b[i1+i]);
          c[j4+j] =-(QT5  * (a[i0+i]-a[i1+i])
             + (0.25  * (a[i0+i]+a[i1+i])-a[i2+i]))
             - (S36 *  b[i0+i]+S72*b[i1+i]);
          c[j2+j] = (QT5  * (a[i0+i]-a[i1+i])
             - (0.25  * (a[i0+i]+a[i1+i])-a[i2+i]))
             - (S72 *  b[i0+i]-S36*b[i1+i]);
          c[j3+j] =-(QT5  * (a[i0+i]-a[i1+i])
             - (0.25  * (a[i0+i]+a[i1+i])-a[i2+i]))
             - (S72 *  b[i0+i]-S36*b[i1+i]);
        ENDL
      }  else {
        qqrt5  = 2.0 * QT5 ;
        ssin36 = 2.0 * S36;
        ssin72 = 2.0 * S72;
        LOOP
          c[j0+j]= 2.0 *(a[i0+i]+a[i1+i]+a[i2+i]);
          c[j1+j]=(2.0 *(a[i0+i]-0.25*(a[i1+i]+a[i2+i]))
            +qqrt5*(a[i1+i]-a[i2+i]))-(ssin72*b[i1+i]+ssin36*b[i2+i]);
          c[j2+j]=(2.0 *(a[i0+i]-0.25*(a[i1+i]+a[i2+i]))
            -qqrt5*(a[i1+i]-a[i2+i]))-(ssin36*b[i1+i]-ssin72*b[i2+i]);
          c[j3+j]=(2.0 *(a[i0+i]-0.25*(a[i1+i]+a[i2+i]))
            -qqrt5*(a[i1+i]-a[i2+i]))+(ssin36*b[i1+i]-ssin72*b[i2+i]);
          c[j4+j]=(2.0 *(a[i0+i]-0.25*(a[i1+i]+a[i2+i]))
            +qqrt5*(a[i1+i]-a[i2+i]))+(ssin72*b[i1+i]+ssin36*b[i2+i]);
        ENDL
      }
      return 0;
    }

    case 6: {
      ia = 0;
      ib = ia+(2*m-la)*inc1;
      ic = ib+2*m*inc1;
      id = ic+2*m*inc1;
      ie = ic;
      iF = ib;
      ja = 0;
      jb = ja+jink;
      jc = jb+jink;
      jd = jc+jink;
      je = jd+jink;
      jf = je+jink;

      if (la != m) /* go to 690 */ {
        LOOP
          c[ja+j] =  (a[ia+i]+a[id+i]) +    (a[ib+i]+a[ic+i]) ;
          c[jd+j] =  (a[ia+i]-a[id+i]) -    (a[ib+i]-a[ic+i]) ;
          c[jb+j] = ((a[ia+i]-a[id+i]) +0.5*(a[ib+i]-a[ic+i]))
                     - S60*(b[ib+i]+b[ic+i]);
          c[jf+j] = ((a[ia+i]-a[id+i]) +0.5*(a[ib+i]-a[ic+i]))
                     + S60*(b[ib+i]+b[ic+i]);
          c[jc+j] = ((a[ia+i]+a[id+i]) -0.5*(a[ib+i]+a[ic+i]))
                     - S60*(b[ib+i]-b[ic+i]);
          c[je+j] = ((a[ia+i]+a[id+i]) -0.5*(a[ib+i]+a[ic+i]))
                     + S60*(b[ib+i]-b[ic+i]);
        ENDL
        ia    += iink;
        iink  += iink;
        ib    += iink;
        ic    += iink;
        id    -= iink;
        ie    -= iink;
        iF    -= iink;
        jbase += jump;
        jump  += jump+jink;

        if (ic != id) /* go to 660 */ {
          FORK
            angle = kpidn = k * pin;
            c1 = cos(angle);
            s1 = sin(angle);
            angle += kpidn;
            c2 = cos(angle);
            s2 = sin(angle);
            angle += kpidn;
            c3 = cos(angle);
            s3 = sin(angle);
            angle += kpidn;
            c4 = cos(angle);
            s4 = sin(angle);
            angle += kpidn;
            c5 = cos(angle);
            s5 = sin(angle);
            ibase=0;
            LOOP
              a11 = a[ie+i] + a[ib+i] + a[ic+i] + a[iF+i];
              a20 = a[ia+i] + a[id+i] - 0.5 * a11;
              a21 = S60*((a[ie+i]+a[ib+i])-(a[ic+i]+a[iF+i]));
              b11 = b[ib+i] - b[ie+i] + b[ic+i] - b[iF+i];
              b20 = b[ia+i] - b[id+i] - 0.5 * b11;
              b21 = S60*((b[ib+i]-b[ie+i])-(b[ic+i]-b[iF+i]));

              c[ja+j] = a[ia+i] + a[id+i] + a11;
              d[ja+j] = b[ia+i] - b[id+i] + b11;
              c[jc+j] = c2*(a20-b21)-s2*(b20+a21);
              d[jc+j] = s2*(a20-b21)+c2*(b20+a21);
              c[je+j] = c4*(a20+b21)-s4*(b20-a21);
              d[je+j] = s4*(a20+b21)+c4*(b20-a21);

              a11 = (a[ie+i]-a[ib+i]) + (a[ic+i]-a[iF+i]);
              b11 = (b[ie+i]+b[ib+i]) - (b[ic+i]+b[iF+i]);
              a20 = (a[ia+i]-a[id+i]) - 0.5 * a11;
              a21 = S60 * ((a[ie+i]-a[ib+i]) - (a[ic+i]-a[iF+i]));
              b20 = (b[ia+i]+b[id+i]) + 0.5 * b11;
              b21 = S60 * ((b[ie+i]+b[ib+i]) + (b[ic+i]+b[iF+i]));

              c[jd+j] = c3*(a[ia+i] - a[id+i] + a11)
                 - s3*(b[ia+i] + b[id+i] - b11);
              d[jd+j] = s3*(a[ia+i] - a[id+i] + a11)
                 + c3*(b[ia+i] + b[id+i] - b11);
              c[jb+j] = c1*(a20-b21)-s1*(b20-a21);
              d[jb+j] = s1*(a20-b21)+c1*(b20-a21);
              c[jf+j] = c5*(a20+b21)-s5*(b20+a21);
              d[jf+j] = s5*(a20+b21)+c5*(b20+a21);
            ENDL
            ia    += iink;
            ib    += iink;
            ic    += iink;
            id    -= iink;
            ie    -= iink;
            iF    -= iink;
            jbase += jump;
          }
          if (ic > id) return 0;
        }
        ibase=0;
        LOOP
          c[ja+j]= a[ib+i] + (a[ia+i] + a[ic+i]);
          c[jd+j]= b[ib+i] - (b[ia+i] + b[ic+i]);
          c[jb+j]= (S60*(a[ia+i]-a[ic+i]))-(0.5*(b[ia+i]+b[ic+i])+b[ib+i]);
          c[jf+j]=-(S60*(a[ia+i]-a[ic+i]))-(0.5*(b[ia+i]+b[ic+i])+b[ib+i]);
          c[jc+j]=  S60*(b[ic+i]-b[ia+i]) +(0.5*(a[ia+i]+a[ic+i])-a[ib+i]);
          c[je+j]=  S60*(b[ic+i]-b[ia+i]) -(0.5*(a[ia+i]+a[ic+i])-a[ib+i]);
        ENDL
      }  else {
        LOOP
          c[ja+j]=(2.0*(a[ia+i]+a[id+i]))+(2.0*(a[ib+i]+a[ic+i]));
          c[jd+j]=(2.0*(a[ia+i]-a[id+i]))-(2.0*(a[ib+i]-a[ic+i]));
          c[jb+j]=(2.0*(a[ia+i]-a[id+i])+(a[ib+i]-a[ic+i]))
            -(D60*(b[ib+i]+b[ic+i]));
          c[jf+j]=(2.0*(a[ia+i]-a[id+i])+(a[ib+i]-a[ic+i]))
            +(D60*(b[ib+i]+b[ic+i]));
          c[jc+j]=(2.0*(a[ia+i]+a[id+i])-(a[ib+i]+a[ic+i]))
            -(D60*(b[ib+i]-b[ic+i]));
          c[je+j]=(2.0*(a[ia+i]+a[id+i])-(a[ib+i]+a[ic+i]))
            +(D60*(b[ib+i]-b[ic+i]));
        ENDL
      }
      return 0;
    }

    case 8: {
      double a0p7,a1p5,a2p6,p073,p074,p152;
      double a0m7,a1m5,a2m6,m073,m074,m152;

      if (la != m) return 3;
      i0  = 0;
      i1  = i0 + iink;
      i2  = i1 + iink;
      i3  = i2 + iink;
      i4  = i3 + iink;
      i5  = i4 + iink;
      i6  = i5 + iink;
      i7  = i6 + iink;
      j0  = 0;
      j1  = j0 + jink;
      j2  = j1 + jink;
      j3  = j2 + jink;
      j4  = j3 + jink;
      j5  = j4 + jink;
      j6  = j5 + jink;
      j7  = j6 + jink;

      LOOP
        a0p7 = a[i0+i] + a[i7+i];
        a0m7 = a[i0+i] - a[i7+i];
        a1p5 = a[i1+i] + a[i5+i];
        a1m5 = a[i1+i] - a[i5+i];
        a2p6 = a[i2+i] + a[i6+i];
        a2m6 = a[i2+i] - a[i6+i];

        p073 = a0p7 + a[i3+i];
        m073 = a0p7 - a[i3+i];

        p074 = 2.0 * (a0m7 + a[i4+i]);
        m074 = 2.0 * (a0m7 - a[i4+i]);

        p152 = M_SQRT2 * (a1m5 + a2p6);
        m152 = M_SQRT2 * (a1m5 - a2p6);

        c[j0+j] = 2.0 * (p073 + a1p5);
        c[j4+j] = 2.0 * (p073 - a1p5);
        c[j2+j] = 2.0 * (m073 - a2m6);
        c[j6+j] = 2.0 * (m073 + a2m6);

        c[j1+j] = m074 + m152;
        c[j5+j] = m074 - m152;
        c[j3+j] = p074 - p152;
        c[j7+j] = p074 + p152;
      ENDL
    }
  }
  return 0;
}

#if defined (SX)
#pragma odir switch,-dv
#endif
int qpassc(double *a, double *b, double *c, double *d,
           int inc1, int inc2, int inc3, int inc4,
           int lot , int n   , int ifac, int la  )
{
/*
     qpassc - performs one pass through data as part;
     of multiple real fft (fourier analysis) routine;

     a is first real input vector;
     b is equivalent to a + ifac * la * inc1
     c is first real output vector;
     d is equivalent to c + la * inc2
     inc1 is the addressing increment for a;
     inc2 is the addressing increment for c;
     inc3 is the increment between input vectors a;
     inc4 is the increment between output vectors c;
     lot is the number of vectors;
     n is the length of the vectors;
     ifac is the current factor of n;
     la is the product of previous factors;
     qpassc returns an error indicator:;
       0 - pass completed without error;
       2 - ifac not catered for;
       3 - ifac only catered for if la=n/ifac;
*/

  int i0,i1,i2,i3,i4,i5,i6,i7;
  int j0,j1,j2,j3,j4,j5,j6,j7;
  int ia,ib,ic;
  int ja,jb,jc;
  int i,j,k;
  int ibase,jbase;
  int iink,jink;
  int ijk;
  int jump;
  int kstop;
  int l;
  int m;

  double a0,a1,a2,a3;
  double b0,b1,b2,b3;
  double c1,c2,c3,c4,c5;
  double s1,s2,s3,s4,s5;
  double w,x,y,z;
  double angle,kpidn,pin;

  m     = n  / ifac;
  iink  = la * inc1;
  jink  = la * inc2;
  jump  = (ifac-1) * iink;
  kstop = (n-ifac) / (2*ifac);
  pin   = 2.0 * M_PI / n;
  ibase = 0;
  jbase = 0;

  switch (ifac) {
    case 2: {
      i0 = j0 = 0;
      i1 = i0 + iink;
      j1 = j0 + inc2 * (m+m-la);
      if (la != m) {
        LOOP
          c[j0+j] = a[i0+i] + a[i1+i];
          c[j1+j] = a[i0+i] - a[i1+i];
        ENDL
        j0    += jink;
        jink  += jink;
        j1    -= jink;
        ibase += jump;
        jump  += jump + iink;

        if (j0 != j1) {
          FORK
            angle = k * pin;
            c1 = cos(angle);
            s1 = sin(angle);
            jbase = 0;
            LOOP
              c[j0+j] = a[i0+i] + c1 * a[i1+i] + s1 * b[i1+i];
              c[j1+j] = a[i0+i] - c1 * a[i1+i] - s1 * b[i1+i];
              d[j0+j] = c1 * b[i1+i] - s1 * a[i1+i] + b[i0+i];
              d[j1+j] = c1 * b[i1+i] - s1 * a[i1+i] - b[i0+i];
            ENDL
            j0    += jink;
            j1    -= jink;
            ibase += jump;
          } /* End FORK */
          if (j0 > j1) return 0;
        } /* End (i0 != i1) */
        jbase = 0;
        LOOP
          c[j0+j] =  a[i0+i];
          d[j1+j] = -a[i1+i];
        ENDL
      }
      else /* (la != m) */ {
        z = 1.0 / n;
        LOOP
          c[j0+j] = z * (a[i0+i] + a[i1+i]);
          c[j1+j] = z * (a[i0+i] - a[i1+i]);
        ENDL
      }
      return 0;
    }

    case 3: {
      ia = 0;
      ib = ia + iink;
      ic = ib + iink;

      ja = 0;
      jb = ja + inc2 * (m+m-la);
      jc = jb;

      if (la != m) /* else 390 */ {
        LOOP
          c[ja+j] = a[ia+i] + a[ib+i] + a[ic+i];
          c[jb+j] = a[ia+i] - 0.5 * (a[ib+i] + a[ic+i]);
          d[jb+j] = S60 * (a[ic+i] - a[ib+i]);
        ENDL
        ja    += jink;
        jink  += jink;
        jb    += jink;
        jc    -= jink;
        ibase += jump;
        jump  += jump + iink;

        if (ja != jc) /* else  360 */ {
          FORK
            angle = k * pin;
            c1 = cos(angle);
            s1 = sin(angle);
            angle += angle;
            c2 = cos(angle);
            s2 = sin(angle);
            jbase = 0;
            LOOP
              a1 = c1 * a[ib+i] + s1 * b[ib+i] + c2 * a[ic+i] + s2 * b[ic+i];
              b1 = c1 * b[ib+i] - s1 * a[ib+i] + c2 * b[ic+i] - s2 * a[ic+i];
              a2 = a[ia+i] - 0.5 * a1;
              b2 = b[ia+i] - 0.5 * b1;
              a3 = S60 * (c1 * a[ib+i] + s1 * b[ib+i] - c2 * a[ic+i] - s2 * b[ic+i]);
              b3 = S60 * (c1 * b[ib+i] - s1 * a[ib+i] - c2 * b[ic+i] + s2 * a[ic+i]);

              c[ja+j] = a[ia+i] + a1;
              d[ja+j] = b[ia+i] + b1;
              c[jb+j] = a2 + b3;
              d[jb+j] = b2 - a3;
              c[jc+j] = a2 - b3;
              d[jc+j] =-b2 - a3;
            ENDL
            ja    += jink;
            jb    += jink;
            jc    -= jink;
            ibase += jump;
          } /* End FORK */
          if (ja > jc) return 0;
        } /* End (ia != ic) */
        jbase = 0;
        LOOP
        /* soweit */
          c[ja+j] = a[ia+i] + 0.5 * (a[ib+i] - a[ic+i]);
          d[ja+j] =-S60 * (a[ib+i] + a[ic+i]);
          c[jb+j] = a[ia+i] - a[ib+i] + a[ic+i];
        ENDL
           }
           else /* (la != m) */ {
        z = 1.0   / n;
        y = S60 / n;
        LOOP
          c[ja+j] = z * (a[ia+i] + a[ib+i] + a[ic+i]);
          c[jb+j] = z * (a[ia+i] - 0.5 * (a[ib+i] + a[ic+i]));
          d[jb+j] = y * (a[ic+i] - a[ib+i]);
        ENDL
      }
      return 0;
    }

    case 4: {
      double a0p2,a1p3;

      i0 = 0;
      i1 = i0 + iink;
      i2 = i1 + iink;
      i3 = i2 + iink;
      j0 = 0;
      j1 = j0 + inc2 * (m+m-la);
      j2 = j1 + inc2 * (m+m   );
      j3 = j1;

      if (la != m) /*else go to 490 */ {
        LOOP
          a0p2 = a[i0+i] + a[i2+i];
          a1p3 = a[i1+i] + a[i3+i];

          c[j0+j] = a0p2 + a1p3;
          c[j2+j] = a0p2 - a1p3;

          c[j1+j] = a[i0+i] - a[i2+i];
          d[j1+j] = a[i3+i] - a[i1+i];
        ENDL
        j0    += jink;
        jink  += jink;
        j1    += jink;
        j2    -= jink;
        j3    -= jink;
        ibase += jump;
        jump  += jump + iink;

        if (j1 != j2) /* else go to 460; */ {
          FORK
            angle = kpidn = k * pin;
            c1 = cos(angle);
            s1 = sin(angle);
            angle += kpidn;
            c2 = cos(angle);
            s2 = sin(angle);
            angle += kpidn;
            c3 = cos(angle);
            s3 = sin(angle);
            jbase=0;
            LOOP
              a0 = a[i0+i] + c2 * a[i2+i] + s2 * b[i2+i];
              a2 = a[i0+i] - c2 * a[i2+i] - s2 * b[i2+i];
              b0 = b[i0+i] + c2 * b[i2+i] - s2 * a[i2+i];
              b2 = b[i0+i] - c2 * b[i2+i] + s2 * a[i2+i];

              a1 = c1*a[i1+i] + s1*b[i1+i] + c3*a[i3+i] + s3*b[i3+i];
              a3 = c1*a[i1+i] + s1*b[i1+i] - c3*a[i3+i] - s3*b[i3+i];
              b1 = c1*b[i1+i] - s1*a[i1+i] + c3*b[i3+i] - s3*a[i3+i];
              b3 = c1*b[i1+i] - s1*a[i1+i] - c3*b[i3+i] + s3*a[i3+i];

              c[j0+j] = a0 + a1;
              c[j2+j] = a0 - a1;
              d[j0+j] = b0 + b1;
              d[j2+j] = b1 - b0;
              c[j1+j] = a2 + b3;
              c[j3+j] = a2 - b3;
              d[j1+j] = b2 - a3;
              d[j3+j] =-b2 - a3;
            ENDL
            j0    += jink;
            j1    += jink;
            j2    -= jink;
            j3    -= jink;
            ibase += jump;
          } /* End FORK */
          if (j1 > j2) return 0;
        } /* End (i1 != i2) */
        jbase=0;
        LOOP
          c[j0+j] = a[i0+i] + SQ2 * (a[i1+i] - a[i3+i]);
          c[j1+j] = a[i0+i] - SQ2 * (a[i1+i] - a[i3+i]);
          d[j0+j] =-a[i2+i] - SQ2 * (a[i1+i] + a[i3+i]);
          d[j1+j] = a[i2+i] - SQ2 * (a[i1+i] + a[i3+i]);
        ENDL
      }
      else /* (la != m) */ {
        z = 1.0 / n;
        LOOP
          a0p2 = a[i0+i] + a[i2+i];
          a1p3 = a[i1+i] + a[i3+i];

          c[j0+j] = z * (a0p2 + a1p3);
          c[j2+j] = z * (a0p2 - a1p3);
          c[j1+j] = z * (a[i0+i] - a[i2+i]);
          d[j1+j] = z * (a[i3+i] - a[i1+i]);
        ENDL
      }
      return 0;
    }

    case 5: {
      double a1p4,a2p3,b1p4,b2p3,a025,b025,asps,bsps,a0pq,b0pq;
      double a1m4,a2m3,b1m4,b2m3,aqrt,bqrt,asms,bsms,a0mq,b0mq;

      i0 = 0;
      i1 = i0 + iink;
      i2 = i1 + iink;
      i3 = i2 + iink;
      i4 = i3 + iink;
      j0 = 0;
      j1 = j0 + inc2 * (m+m-la);
      j2 = j1 + inc2 * (m+m);
      j3 = j2;
      j4 = j1;

      if (la != m) /* else go to 590; */ {
        LOOP
          a1p4 = a[i1+i] + a[i4+i];
          a1m4 = a[i1+i] - a[i4+i];
          a2p3 = a[i2+i] + a[i3+i];
          a2m3 = a[i2+i] - a[i3+i];

          a025 = a[i0+i] - 0.25 * (a1p4 + a2p3);
          aqrt =           QT5 * (a1p4 - a2p3);

          c[j0+j] = a[i0+i] + a1p4 + a2p3;
          c[j1+j] = a025 + aqrt;
          c[j2+j] = a025 - aqrt;
          d[j1+j] = -S72 * a1m4 - S36 * a2m3;
          d[j2+j] = -S36 * a1m4 + S72 * a2m3;
        ENDL
        j0    += jink;
        jink  += jink;
        j1    += jink;
        j2    += jink;
        j3    -= jink;
        j4    -= jink;
        ibase += jump;
        jump  += jump + iink;

        if (j1 != j3) {
          FORK
            angle = kpidn = k * pin;
            c1 = cos(angle);
            s1 = sin(angle);
            angle += kpidn;
            c2 = cos(angle);
            s2 = sin(angle);
            angle += kpidn;
            c3 = cos(angle);
            s3 = sin(angle);
            angle += kpidn;
            c4 = cos(angle);
            s4 = sin(angle);
            jbase=0;
            LOOP
              a1p4 = c1*a[i1+i] + s1*b[i1+i] + c4*a[i4+i] + s4*b[i4+i];
              a1m4 = c1*a[i1+i] + s1*b[i1+i] - c4*a[i4+i] - s4*b[i4+i];
              a2p3 = c2*a[i2+i] + s2*b[i2+i] + c3*a[i3+i] + s3*b[i3+i];
              a2m3 = c2*a[i2+i] + s2*b[i2+i] - c3*a[i3+i] - s3*b[i3+i];
              b1p4 = c1*b[i1+i] - s1*a[i1+i] + c4*b[i4+i] - s4*a[i4+i];
              b1m4 = c1*b[i1+i] - s1*a[i1+i] - c4*b[i4+i] + s4*a[i4+i];
              b2p3 = c2*b[i2+i] - s2*a[i2+i] + c3*b[i3+i] - s3*a[i3+i];
              b2m3 = c2*b[i2+i] - s2*a[i2+i] - c3*b[i3+i] + s3*a[i3+i];

              a025 = a[i0+i] - 0.25 * (a1p4 + a2p3);
              aqrt =           QT5  * (a1p4 - a2p3);
              b025 = b[i0+i] - 0.25 * (b1p4 + b2p3);
              bqrt =           QT5  * (b1p4 - b2p3);

              a0pq = a025 + aqrt;
              a0mq = a025 - aqrt;
              b0pq = b025 + bqrt;
              b0mq = b025 - bqrt;

              asps = S72 * a1m4 + S36 * a2m3;
              asms = S36 * a1m4 - S72 * a2m3;
              bsps = S72 * b1m4 + S36 * b2m3;
              bsms = S36 * b1m4 - S72 * b2m3;

              c[j0+j] = a[i0+i] + a1p4 + a2p3;
              c[j1+j] = a0pq + bsps;
              c[j2+j] = a0mq + bsms;
              c[j3+j] = a0mq - bsms;
              c[j4+j] = a0pq - bsps;
              d[j0+j] = b[i0+i] + b1p4 + b2p3;
              d[j1+j] = b0pq - asps;
              d[j2+j] = b0mq - asms;
              d[j3+j] =-b0mq - asms;
              d[j4+j] =-b0pq - asps;
            ENDL
            j0    += jink;
            j1    += jink;
            j2    += jink;
            j3    -= jink;
            j4    -= jink;
            ibase += jump;
          } /* End FORK */
          if (j1 > j3) return 0;
        } /* End (jb != jd) */
        jbase=0;
        LOOP
          a1p4 = a[i1+i] + a[i4+i];
          a1m4 = a[i1+i] - a[i4+i];
          a2p3 = a[i2+i] + a[i3+i];
          a2m3 = a[i2+i] - a[i3+i];

          a025 = a[i0+i] + 0.25 * (a1m4 - a2m3);
          aqrt =           QT5 * (a1m4 + a2m3);

          c[j0+j] = a025 + aqrt;
          c[j1+j] = a025 - aqrt;
          c[j2+j] = a[i0+i] - a1m4 + a2m3;
          d[j0+j] = -S36 * a1p4 - S72 * a2p3;
          d[j1+j] = -S72 * a1p4 + S36 * a2p3;

        ENDL
      }  else {
        z = 1.0 / n;
        y = QT5 / n;
        x = S36 / n;
        w = S72 / n;

        LOOP
          a1p4 = a[i1+i] + a[i4+i];
          a1m4 = a[i1+i] - a[i4+i];
          a2p3 = a[i2+i] + a[i3+i];
          a2m3 = a[i2+i] - a[i3+i];

          a025 = z * (a[i0+i] - 0.25 * (a1p4 + a2p3));
          aqrt = y * (a1p4 - a2p3);

          c[j0+j] = z * (a[i0+i] + a1p4 + a2p3);
          c[j1+j] = a025 + aqrt;
          c[j2+j] = a025 - aqrt;
          d[j1+j] = -w * a1m4 - x * a2m3;
          d[j2+j] =  w * a2m3 - x * a1m4;
        ENDL
      }
      return 0;
    }

    case 6: {
      double ab1a,ab2a,ab3a,ab4a,ab5a;
      double ab1b,ab2b,ab3b,ab4b,ab5b;
      double a0p3,a1p4,a1p5,a2p4,a2p5;
      double a0m3,a1m4,a1m5,a2m4,a2m5;
      double b1p4,b2p5;
      double b1m4,b2m5;
      double ap05,bp05,ap60,bp60;
      double am05,bm05,am60,bm60;

      i0 = 0;
      i1 = i0 + iink;
      i2 = i1 + iink;
      i3 = i2 + iink;
      i4 = i3 + iink;
      i5 = i4 + iink;
      j0 = 0;
      j1 = j0 + inc2 * (m+m-la);
      j2 = j1 + inc2 * (m+m);
      j3 = j2 + inc2 * (m+m);
      j4 = j2;
      j5 = j1;

      if (la != m) {
        LOOP
          a0p3 = a[i0+i] + a[i3+i];
          a0m3 = a[i0+i] - a[i3+i];
          a1p4 = a[i1+i] + a[i4+i];
          a1m4 = a[i1+i] - a[i4+i];
          a2p5 = a[i2+i] + a[i5+i];
          a2m5 = a[i2+i] - a[i5+i];

          c[j0+j] = a0p3 + a1p4 + a2p5;
          c[j3+j] = a0m3 + a2m5 - a1m4;

          c[j1+j] = a0m3 - 0.5 * (a2m5 - a1m4);
          c[j2+j] = a0p3 - 0.5 * (a1p4 + a2p5);

          d[j1+j] = S60 * (-a2m5 - a1m4);
          d[j2+j] = S60 * ( a2p5 - a1p4);
        ENDL
        j0    += jink;
        jink  += jink;
        j1    += jink;
        j2    += jink;
        j3    -= jink;
        j4    -= jink;
        j5    -= jink;
        ibase += jump;
        jump  += jump+iink;

        if (j2 != j3) {
          FORK
            angle = kpidn = k * pin;
            c1 = cos(angle);
            s1 = sin(angle);
            angle += kpidn;
            c2 = cos(angle);
            s2 = sin(angle);
            angle += kpidn;
            c3 = cos(angle);
            s3 = sin(angle);
            angle += kpidn;
            c4 = cos(angle);
            s4 = sin(angle);
            angle += kpidn;
            c5 = cos(angle);
            s5 = sin(angle);
            jbase = 0;
            LOOP
              ab1a = c1 * a[i1+i] + s1 * b[i1+i];
              ab1b = c1 * b[i1+i] - s1 * a[i1+i];
              ab2a = c2 * a[i2+i] + s2 * b[i2+i];
              ab2b = c2 * b[i2+i] - s2 * a[i2+i];
              ab3a = c3 * a[i3+i] + s3 * b[i3+i];
              ab3b = c3 * b[i3+i] - s3 * a[i3+i];
              ab4a = c4 * a[i4+i] + s4 * b[i4+i];
              ab4b = c4 * b[i4+i] - s4 * a[i4+i];
              ab5a = c5 * a[i5+i] + s5 * b[i5+i];
              ab5b = c5 * b[i5+i] - s5 * a[i5+i];

              a1p4 = ab1a + ab4a;
              a1m4 = ab1a - ab4a;
              a2p5 = ab2a + ab5a;
              a2m5 = ab2a - ab5a;

              b1p4 = ab1b + ab4b;
              b1m4 = ab1b - ab4b;
              b2p5 = ab2b + ab5b;
              b2m5 = ab2b - ab5b;

              ap05 = a[i0+i] + ab3a - 0.5 * (a1p4 + a2p5);
              bp05 = b[i0+i] + ab3b - 0.5 * (b1p4 + b2p5);
              am05 = a[i0+i] - ab3a - 0.5 * (a2m5 - a1m4);
              bm05 =-b[i0+i] + ab3b - 0.5 * (b1m4 - b2m5);

              ap60 = S60 * ( a2p5 - a1p4);
              bp60 = S60 * ( b2p5 - b1p4);
              am60 = S60 * (-a2m5 - a1m4);
              bm60 = S60 * (-b2m5 - b1m4);

              c[j0+j] = a[i0+i] + ab3a + a1p4 + a2p5;
              d[j0+j] = b[i0+i] + ab3b + b1p4 + b2p5;
              c[j1+j] = am05 - bm60;
              d[j1+j] = am60 - bm05;
              c[j2+j] = ap05 - bp60;
              d[j2+j] = ap60 + bp05;
              c[j3+j] = a[i0+i] - ab3a - a1m4 + a2m5;
              d[j3+j] =-b[i0+i] + ab3b + b1m4 - b2m5;
              c[j4+j] = ap05 + bp60;
              d[j4+j] = ap60 - bp05;
              c[j5+j] = am05 + bm60;
              d[j5+j] = am60 + bm05;
            ENDL
            j0    += jink;
            j1    += jink;
            j2    += jink;
            j3    -= jink;
            j4    -= jink;
            j5    -= jink;
            ibase += jump;
          }
          if (j2 > j3) return 0;
        }
        jbase = 0;
        LOOP
          a1p5 = a[i1+i] + a[i5+i];
          a1m5 = a[i1+i] - a[i5+i];
          a2p4 = a[i2+i] + a[i4+i];
          a2m4 = a[i2+i] - a[i4+i];

          c[j0+j] = a[i0+i] + 0.5 * a2m4 + S60 * a1m5;
          d[j0+j] =-a[i3+i] - 0.5 * a1p5 - S60 * a2p4;
          c[j1+j] = a[i0+i] - a2m4;
          d[j1+j] = a[i3+i] - a1p5;
          c[j2+j] = a[i0+i] + 0.5 * a2m4 - S60 * a1m5;
          d[j2+j] =-a[i3+i] - 0.5 * a1p5 + S60 * a2p4;
        ENDL
      }  else {
        z = 1.0 / n;
        y = S60 / n;
        LOOP
          a0p3 = a[i0+i] + a[i3+i];
          a0m3 = a[i0+i] - a[i3+i];
          a1p4 = a[i1+i] + a[i4+i];
          a1m4 = a[i1+i] - a[i4+i];
          a2p5 = a[i2+i] + a[i5+i];
          a2m5 = a[i2+i] - a[i5+i];

          c[j0+j] = z * (a0p3 + a1p4 + a2p5);
          c[j3+j] = z * (a0m3 + a2m5 - a1m4);

          c[j1+j] = z * (a0m3 - 0.5 * (a2m5 - a1m4));
          c[j2+j] = z * (a0p3 - 0.5 * (a1p4 + a2p5));

          d[j1+j] = y * (-a2m5 - a1m4);
          d[j2+j] = y * ( a2p5 - a1p4);
        ENDL
      }
      return 0;
    }

    case 8: {
      double a0p4,a1p5,a2p6,a3p7;
      double a0m4,a1m5,a2m6,a3m7;

      if (la != m) return 3;
      i0 = 0;
      i1 = i0 + iink;
      i2 = i1 + iink;
      i3 = i2 + iink;
      i4 = i3 + iink;
      i5 = i4 + iink;
      i6 = i5 + iink;
      i7 = i6 + iink;
      j0 = 0;
      j1 = j0 + jink;
      j2 = j1 + jink;
      j3 = j2 + jink;
      j4 = j3 + jink;
      j5 = j4 + jink;
      j6 = j5 + jink;
      j7 = j6 + jink;
      z  = 1.0      / n;
      y  = SQ2 / n;

      LOOP
        a0p4 = a[i0+i] + a[i4+i];
        a0m4 = a[i0+i] - a[i4+i];
        a1p5 = a[i1+i] + a[i5+i];
        a1m5 = a[i1+i] - a[i5+i];
        a2p6 = a[i2+i] + a[i6+i];
        a2m6 = a[i2+i] - a[i6+i];
        a3p7 = a[i3+i] + a[i7+i];
        a3m7 = a[i3+i] - a[i7+i];

        c[j0+j] = z * (a0p4 + a1p5 + a2p6 + a3p7);
        c[j7+j] = z * (a0p4 - a1p5 + a2p6 - a3p7);

        c[j3+j] = z * (a0p4 - a2p6);
        c[j4+j] = z * (a3p7 - a1p5);

        c[j1+j] = z * a0m4 + y * (a1m5 - a3m7);
        c[j5+j] = z * a0m4 - y * (a1m5 - a3m7);
        c[j2+j] =-z * a2m6 - y * (a1m5 + a3m7);
        c[j6+j] = z * a2m6 - y * (a1m5 + a3m7);
      ENDL
    }
  }
  return 0;
}

/* ====================== */
/* Fast Fourier Transform */
/* ====================== */

void fc2gp(double *fc, double *gp, int Lat, int Lon, int Lev, int Fou)
{
   int Lot;
   int fou;
   int ia;
   int ifac;
   int j;
   int jump;
   int k;
   int la;
   int lat;
   int lev;
   int lon;
   int nfax;
   int rix;
   int wix;

   double *wfc;
   double *wgp;
   double *wpt;

   static int ifax[10];

/* fc2gp performs fourier to gridpoint transforms using       */
/* multiple fast fourier transform of length Lon              */
/*                                                            */
/* fc  - real array of fourier coefficients fc[Lev][Fou][Lat] */
/* gp  - real array of gridpoints           gp[Lev][Lat][Lon] */
/* Lat - Number of latitudes                                  */
/* Lon - Number of longitudes                                 */
/* Lev - Number of levels                                     */
/* Fou - Number of fourier coefficients on 1 latitude         */

/* x(j) = sum(k=0,...,n-1)(c(k)*exp(2*i*j*k*pi/Lon))          */
/*        where c(k) = a(k) + i*b(k) and c(n-k) = a(k)-i*b(k) */

   jump = (Lon + 2) | 1;
   Lot  = Lev * Lat;
   fft_set(ifax,Lon);
   nfax = ifax[0];

   wfc = DoubleAlloc(Lot * jump,"wfc in fc2gp");
   wgp = DoubleAlloc(Lot * jump,"wgp in fc2gp");

   for (lev = 0; lev < Lev; ++lev) {
      for (lat = 0; lat < Lat; ++lat) {
         wix = jump * (lat + lev * Lat);
         rix = lat  +  lev * Lat * Fou;
         for (fou = 0  ; fou < Fou  ; ++fou)
            wfc[wix + fou] = fc[rix + fou * Lat];
         for (fou = Fou; fou < jump; ++fou) wfc[wix + fou] = 0.0;
         wfc[wix+1] = 0.5 * wfc[wix];
      }
   }

   ia = 1;
   la = 1;

   for (k = 0; k < nfax; ++k) {
      ifac = ifax[k+1];
      if (k&1) rpassc(wgp,wgp+la,wfc+ia,wfc+ia+ifac*la,
                      1,1,jump,jump,Lot,Lon,ifac,la);
      else     rpassc(wfc+ia,wfc+ia+la,wgp,wgp+ifac*la,
                      1,1,jump,jump,Lot,Lon,ifac,la);
      la *= ifac;
      ia  = 0;
   }

   if (nfax & 1) wpt = wgp;
   else          wpt = wfc;

   for (j = 0; j < Lot ; ++j)
      for (lon = 0; lon < Lon; ++lon)
         gp[lon + j * Lon] = wpt[lon + j * jump];

   wgp = FreeMem(wgp);
   wfc = FreeMem(wfc);
}

void gp2fc(double *gp, double *fc, int Lat, int Lon, int Lev, int Fou)
{
   int Lot;
   int fou;
   int ia;
   int ifac;
   int jump;
   int k;
   int la;
   int lat;
   int lev;
   int lon;
   int lot;
   int nfax;
   int rix;
   int wix;

   double *wfc;
   double *wgp;
   double *wpt;

   static int ifax[10];

/* fc2gp performs fourier to gridpoint transforms using       */
/* multiple fast fourier transform of length Lon              */
/*                                                            */
/* fc  - real array of fourier coefficients fc[Lev][Fou][Lat] */
/* gp  - real array of gridpoints           gp[Lev][Lat][Lon] */
/* Lat - Number of latitudes                                  */
/* Lon - Number of longitudes                                 */
/* Lev - Number of levels                                     */
/* Fou - Number of fourier coefficients on 1 latitude         */

/* a(k) =  (1/n) * sum(j=0,...,n-1)(x(j) * cos(2*j*k*pi/n))   */
/* b(k) = -(1/n) * sum(j=0,...,n-1)(x(j) * sin(2*j*k*pi/n))   */

   jump = (Lon + 2) | 1;
   Lot  = Lev * Lat;
   fft_set(ifax,Lon);
   nfax = ifax[0];

   wfc = DoubleAlloc(Lot * jump,"wfc in gp2fc");
   wgp = DoubleAlloc(Lot * jump,"wgp in gp2fc");

   rix = 0;
   wix = 0;
   for (lot = 0; lot < Lot; ++lot) {
      for (lon = 0; lon < Lon; ++lon) wgp[wix+lon] = gp[rix+lon];
      wgp[wix+Lon  ] = 0.0;
      wgp[wix+Lon+1] = 0.0;
      rix += Lon;
      wix += jump;
   }

   ia = 0;
   la = Lon;

   for (k = 0; k < nfax; ++k) {
      ifac = ifax[nfax-k];
      la  /= ifac;
      if (k & 1) qpassc(wfc,wfc+ifac*la,wgp+ia,wgp+ia+la,
                        1,1,jump,jump,Lot,Lon,ifac,la);
      else       qpassc(wgp+ia,wgp+ia+ifac*la,wfc,wfc+la,
                        1,1,jump,jump,Lot,Lon,ifac,la);
      ia = 1;
   }

   if (nfax & 1) wpt = wfc;
   else          wpt = wgp+1;

   for (lev = 0; lev < Lev; ++lev) {
      for (lat = 0; lat < Lat; ++lat) {
         rix = jump * (lat + lev * Lat);
         wix = lat + lev * Lat * Fou;
         fc[wix    ] = wpt[rix];
         fc[wix+Lat] = 0.0;
         for (fou = 2  ; fou < Fou  ; ++fou)
            fc[wix + fou * Lat] = wpt[rix + fou - 1];
      }
   }
   wgp = FreeMem(wgp);
   wfc = FreeMem(wfc);
}

int GRB_Fatal   = 1;	/* If set to 1, exit on fatal error */
int GRB_Verbose = 1;	/* If set to 1, errors are reported */
int GRB_Debug   = 0;    /* If set to 1, debugging           */

void
gribError (const char *func, const char *fmt, ...)
{
  va_list args;
	
  va_start (args, fmt);
  fprintf (stderr, "GRB error (%s) : ", func);
  vfprintf (stderr, fmt, args);
  fprintf (stderr, "\n");
  va_end (args);
  if (GRB_Fatal)
    exit(1);
	
  return;
}

void
gribWarning (const char *func, const char *fmt, ...)
{
  va_list args;
	
  va_start (args, fmt);
  if (GRB_Verbose)
    {
      fprintf (stderr, "GRB warning (%s) : ", func);
      vfprintf (stderr, fmt, args);
      fprintf (stderr, "\n");
    }
  va_end (args);
	
  return;
}

void
gribPrint (const char *func, const char *fmt, ...)
{
  va_list args;
	
  va_start (args, fmt);
  fprintf (stderr, "%-18s : ", func);
  vfprintf (stderr, fmt, args);
  fprintf (stderr, "\n");
  va_end (args);
	
  return;
}

int
gribFileSeek (FILE *fp, int *offset)
{
  /* position file pointer after GRIB */
  static char func[] = "gribFileSeek";
  int buffersize = 1024;
  unsigned char buffer[1024];
  int retry = 1024;
  int i;

  *offset = 0;

  buffer[0] = fgetc(fp);
  buffer[1] = fgetc(fp);
  buffer[2] = fgetc(fp);
  buffer[3] = fgetc(fp);

  while ( retry-- )
    {
      for (i = 0; i < buffersize-4; i++)
	{
	  if (buffer[i  ] == 'G' && 
	      buffer[i+1] == 'R' &&
	      buffer[i+2] == 'I' &&
	      buffer[i+3] == 'B')
	    {
	      if ( GRB_Debug )
		gribPrint(func, "record offset = %d", *offset);
	      return (0);
	    }
	  else
	    {
	      if (feof(fp)) return (-1);
	      buffer[i+4] = fgetc (fp);
	      (*offset)++;
	    }
	}
      buffer[0] = buffer[i  ];
      buffer[1] = buffer[i+1];
      buffer[2] = buffer[i+2];
      buffer[3] = buffer[i+3];
    }

  if ( GRB_Debug )
    gribPrint(func, "record offset = %d", *offset);

  return (1);
}

int
gribReadSize (FILE *fp)
{
  static char func[] = "gribReadSize";
  int gribversion, gribsize;
  fpos_t pos;


  fgetpos(fp, &pos); 

  gribsize = (fgetc (fp) << 16) + (fgetc (fp) << 8) + fgetc (fp);

  gribversion = fgetc (fp);

  if ( gribsize == 24 )
    {
      if ( gribversion != 1 ) gribversion = 0;
    }

  if (GRB_Debug)
    gribPrint(func, "gribversion = %d", gribversion);

  if ( gribversion == 0 )
    {
      int pdssize = 0, gdssize = 0, bmssize = 0, bdssize = 0;
      int issize = 4, essize = 4;
      int flag;

      pdssize = gribsize;
      fseek(fp, (long) 3, SEEK_CUR);
      if (GRB_Debug)
	gribPrint(func, "pdssize     = %d", pdssize);
      flag =  fgetc (fp);
      if (GRB_Debug)
	gribPrint(func, "flag        = %d", flag);
  
      fseek(fp, (long) pdssize-8, SEEK_CUR);

      if ( flag & 128 )
	{
	  gdssize = (fgetc (fp) << 16) + (fgetc (fp) << 8) + fgetc (fp);
	  fseek(fp, (long) gdssize-3, SEEK_CUR);
	}
      if (GRB_Debug)
	gribPrint(func, "gdssize     = %d", gdssize);

      if ( flag &  64 )
	{
	  bmssize = (fgetc (fp) << 16) + (fgetc (fp) << 8) + fgetc (fp);
	  fseek(fp, (long) bmssize-3, SEEK_CUR);
	}
      if (GRB_Debug)
	gribPrint(func, "bmssize     = %d", bmssize);

      bdssize = (fgetc (fp) << 16) + (fgetc (fp) << 8) + fgetc (fp);

      if (GRB_Debug)
	gribPrint(func, "bdssize     = %d", bdssize);

      gribsize = issize + pdssize + gdssize + bmssize + bdssize + essize;
    }
  else if ( gribversion != 1 )
    {
      gribsize = 0;
    }

  if ( feof(fp) ) gribsize = 0;

  if (GRB_Debug)
    gribPrint(func, "gribsize    = %d", gribsize);

  fsetpos(fp, &pos);

  return (gribsize);
}

int
grb_read_record (FILE *fp, void *record, size_t recordsize)
{
  static char func[] = "read_gribrecord";
  size_t nread;

  nread = fread(record, 1, recordsize, fp);

  if ( nread < recordsize )
    {
      gribError(func, "can't read record with %d bytes\n", recordsize);
      return (-1);
    }

  return (0);
}

int    PackComplex;         /* Complex packing of sph. harm.         */
int    PackComplexScale   ; /* Scale factor for complex packed data  */
int    PackComplexStart   ; /* Start of complex packed data subfield */

char *GribAbortMsg[] =
{
   "==== G R I B - E R R O R ====\n",
   "Could not read input file\n",
   "Synchronization failed\n",
   "Uncomplete record\n",
   "Record too long\n",
   "Block read error\n",
   "Too many blocks\n"
};

void GribAbort(int errno)
{
   fprintf (stderr, GribAbortMsg[0]);
   fprintf (stderr, GribAbortMsg[errno]);
   fprintf (stderr, GribAbortMsg[0]);
   exit (1);
}

/* ============ */
/* ScaleComplex */
/* ============ */

void ScaleComplex (double *fpdata, int Truncation)
{
   double power;
   double scale[400];
   int  n, m;
   int  index;

   if (PackComplexScale < -10000 || PackComplexScale > 10000) {
      fprintf(stderr, " ScaleComplex: Invalid power given %6d\n",
                        PackComplexScale);
      return;
   }

/* Setup scaling factors = n(n+1)^^p for n = 1 to Truncation */

   if (PackComplexScale == 0) return;

   power = (double) PackComplexScale / 1000.;
   scale[0] = 1.0;

   for (n = 1; n <= Truncation; n++) {
      if (PackComplexScale != 1000)
         scale[n] = 1.0 / pow((double) (n*(n+1)), power);
      else
         scale[n] = 1.0 /     (double) (n*(n+1));
   }

/* Scale the values */

   index = 0;

   for (m = 0; m < PackComplexStart; m++)
   for (n = m; n <= Truncation;      n++) {
      if (n >= PackComplexStart) {
          fpdata[index  ] *= scale[n];
          fpdata[index+1] *= scale[n];
      }
      index += 2;
   }

   for (m = PackComplexStart; m <= Truncation; m++)
   for (n = m;                n <= Truncation; n++) {
       fpdata[index  ] *= scale[n];
       fpdata[index+1] *= scale[n];
       index += 2;
   }
}

/* ============== */
/* ScatterComplex */
/* ============== */

void ScatterComplex(double *fpdata, int Truncation, int DimSP)
{
   double *fphelp = DoubleAlloc(DimSP,"fphelp.sp");
   int  m, n;
   int  index, inext;

   index = inext = 0;

   for (m = 0; m <= PackComplexStart; m++)
   for (n = m; n <= Truncation;       n++) {
       if (PackComplexStart >= n) {
          fphelp[index  ] = fpdata[inext++];
          fphelp[index+1] = fpdata[inext++];
       }
       index += 2;
   }

   index = 0;
   for (m = 0; m <= Truncation;       m++)
   for (n = m; n <= Truncation;       n++) {
       if (n > PackComplexStart) {
          fphelp[index  ] = fpdata[inext++];
          fphelp[index+1] = fpdata[inext++];
       }
       index += 2;
   }

   for (m = 0; m < DimSP; m++) fpdata[m] = fphelp[m];
   fphelp = FreeMem(fphelp);
}

#define ulong unsigned long
#define get2byte(a,b) ((long) ((ulong) (a<<8) + (ulong) b))


char *MoName[13] = {"   ","Jan","Feb","Mar","Apr","May","Jun",
                          "Jul","Aug","Sep","Oct","Nov","Dec"};
int MoSizeReal[13] = {0,31,28,31,30,31,30,31,31,30,31,30,31};
int MoSizeModl[13] = {0,30,30,30,30,30,30,30,30,30,30,30,30};
int *MoSize;

/* =============================== */
/* Time handling routines for GRIB */
/* =============================== */

int leap_year (int year)  {
   int i,y;

   if(MoSize[1] == 30) return (0);
   y = year;

   i = y / 4;
   i = (i*4) - y;
   if (i != 0) return (0);

   i = y / 100;
   i = (i*100) - y;
   if (i != 0) return (1);

   i = y / 400;
   i = (i*400) - y;
   if (i != 0) return (0);

   return (1);
}

/* Add an offset to a time.  Output to date.                          */

struct Date timadj (int hours, struct Date date) {
   int i;
   int cont;

   /* Add  hours directly.  Then normalize
      to days, then normalize extra days to months/years.             */

   date.hr += hours;

   if (date.mn > 59) {
      i = date.mn / 60;
      date.hr += i;
      date.mn = date.mn - (i*60);
   }
   if (date.hr > 23) {
      i = date.hr / 24;
      date.dy += i;
      date.hr = date.hr - (i*24);
   }

   cont = 1;
   while (date.dy > MoSize[date.mo] && cont) {
      if (date.mo == 2 && leap_year(date.yr)) {
         if (date.dy  ==  29) cont = 0;
         else {
           date.dy -= 29;
           date.mo++;
         }
      } else {
         date.dy -= MoSize[date.mo];
         date.mo++;
      }
      while (date.mo > 12) {date.mo -=12; date.yr++;}
   }
   return date;
}

void decodeField(double *field)
{
  int  TReal = 0, Hour;
  int time1 = 0;
  int klenp;
  int kret, kword;
  extern int TruncationGrib;
  extern CHAR *grib;
  extern int    LevelGrib     ; /* grib[14]   Block1[10]   1 or 2 Bytes       */
  extern int    LayerGrib     ; /* grib[14]   Block1[10]   1 or 2 Bytes       */
  extern int    PackComplex;         /* Complex packing of sph. harm.         */
  extern int    PackComplexScale   ; /* Scale factor for complex packed data  */
  extern int    PackComplexStart   ; /* Start of complex packed data subfield */
  extern int    sec0[2], sec1[1024], sec2[1024], sec3[2], sec4[512];
  extern double psec2[512], psec3[2];
  extern struct Date NewDate;
  extern int CodeGrib;
  extern int RepGrib;
  extern int kleng;
  extern int TypeOfLevel;

  klenp = MAX_DimGP;
  /*      kleng : set in GribFromFile */
  gribexdp(sec0, sec1, sec2, psec2, sec3,
	   psec3, sec4, field, klenp, (void *) grib,
	   kleng, &kword, "D", &kret);

  CodeGrib       = GribParameter(sec1);
  TypeOfLevel    = GribLevelType(sec1);

  if ( sec1[6] == LEV_DOWN  || sec1[6] == LEV_HEIGHT )
    {
      LevelGrib = 0;
      LayerGrib = GribLevel1(sec1);
    }
  else
    {
      LevelGrib = GribLevel1(sec1);
      LayerGrib = 0;
    }

  if ( GribLevelType(sec1) == LEV_GROUND ) LevelGrib = 0;

  if ( GribParameter(sec1) == LNPSCODE ) LevelGrib = 0;

  NewDate.yr  = sec1[9];
  NewDate.mo  = sec1[10];
  NewDate.dy  = sec1[11];
  NewDate.hr  = sec1[12];
  NewDate.mn  = sec1[13];

  if ( GribTimeRange(sec1) == 113 ) time1 = GribTimePeriod1(sec1);
  if ( GribTimeRange(sec1) ==   0 ) time1 = GribTimePeriod1(sec1);
  if ( GribTimeRange(sec1) ==   4 ) time1 = GribTimePeriod2(sec1);
  if ( GribTimeRange(sec1) ==   2 ) time1 = GribTimePeriod2(sec1);
  if ( GribTimeRange(sec1) ==  10 ) time1 = get2byte(GribTimePeriod1(sec1), GribTimePeriod2(sec1));

  NewDate.ct = sec1[20];

  NewDate.yr += (NewDate.ct - 1) * 100;

  Hour = 0;
  if ( time1 )
    {
      switch ( GribTimeUnit(sec1) )
	{
	case  1: Hour = time1;
	         break;
	case  2: Hour = time1 * 24;
                 break;
	case  3: Hour = time1 * 24 * 30;
                 break;
	case  4: Hour = time1 * 24 * 30 * 12;
                 break;
	case  5: Hour = time1 * 24 * 30 * 12 *  10;
                 break;
	case  6: Hour = time1 * 24 * 30 * 12 *  30;
                 break;
	case  7: Hour = time1 * 24 * 30 * 12 * 100;
                 break;
	default: Abort("\n *** Unknown GRIB Time-Unit ***\n");
	}
      /*
      TReal = ((igrib[41] == 3) && (BlockLength > 28));
      */
      if ( TReal ) MoSize = &MoSizeReal[0];
      else         MoSize = &MoSizeModl[0];
      NewDate = timadj(Hour, NewDate);
    }

  NewDate.YYMMDD = 10000*NewDate.yr + 100*NewDate.mo + NewDate.dy;
  NewDate.HHMM   = 100*NewDate.hr + NewDate.mn;

  RepGrib = sec2[0] ;

  if ( RepGrib == REP_SPECTRAL )
    TruncationGrib = sec2[1];

  if ( RepGrib == REP_SPECTRAL )
    {
      PackComplex = 0;
      if ( sec4[3] == 64 ) PackComplex = 1;
      PackComplexScale = sec4[16];
      PackComplexStart = sec4[19];
    }
  
  return;
}

void
codegb(double *field,    /* Pointer to input data */
       int     datalen,  /* length of data record */
       char   *CGrib,    /* packed GRIB record */
       int    *griblen,  /* length of GRIB record */
       int code, int level)
{
  static char func[] = "codegb";
  int kword;
  char *grib;
  int klenp, kleng;
  int kret;
  extern int    sec0[2], sec1[1024], sec2[1024], sec3[2], sec4[512];
  extern double psec2[512], psec3[2];
  extern int DimGPOut;
  extern int LevelType;
  extern struct Date OldDate;
  extern int Mean;
  extern int MeanCount;
  extern int nvct;
  extern double   *vct;
  extern int    Representation;
  extern int    TruncationOut;
  extern int    LatitudesOut  ;
  extern int    LatitudeFirst ;
  extern int    LatitudeLast  ;
  extern int    LongitudesOut ;
  extern int    LongitudeFirst;
  extern int    LongitudeLast ;
  extern int    LongitudeIncr ;
  extern int    DimGPOut      ;
  extern int    DimSPOut      ;

  kleng = DimGPOut*7;
  grib  = CGrib;
  klenp = datalen;
  
  GribParameter(sec1) = code;
  GribLevelType(sec1) = LevelType;
  GribLevel1(sec1)    = level;

  GribYear(sec1)   = OldDate.yr % 100;
  GribMonth(sec1)  = OldDate.mo;
  GribDay(sec1)    = OldDate.dy;
  GribHour(sec1)   = OldDate.hr;
  GribMinute(sec1) = OldDate.mn;

  GribCentury(sec1) = OldDate.ct;

  if ( Mean > 0 )
    GribNumAvg(sec1) = MeanCount;

  GribGridType(sec2) = Representation;

  if ( GribLevelType(sec1) == LEV_HYBRID )
    {
      int i;
      GribNumVCP(sec2) = nvct;
      for (i = 0; i < GribNumVCP(sec2); ++i)
	psec2[10+i] = vct[i];
    }
  else
    {
      GribNumVCP(sec2) = 0;
    }

  if ( GribGridType(sec2) == REP_SPECTRAL )
    {
      GribPentaJ(sec2)  = TruncationOut;
      GribPentaK(sec2)  = TruncationOut;
      GribPentaM(sec2)  = TruncationOut;
      GribRepType(sec2) = 1;
      GribRepMode(sec2) = 1;
    }
  else if ( GribGridType(sec2) == REP_GAUSS )
    {
      GribNumLon(sec2)   = LongitudesOut;
      GribNumLat(sec2)   = LatitudesOut;
      GribFirstLat(sec2) = LatitudeFirst;
      GribFirstLon(sec2) = LongitudeFirst;
      GribLastLat(sec2)  = LatitudeLast;
      GribLastLon(sec2)  = LongitudeLast;
      GribLonIncr(sec2)  = LongitudeIncr;
      GribNumPar(sec2)   = GribNumLat(sec2)/2;
    }
  else
    {
      Error(func, "Unsupported Grid type: %d", GribGridType(sec2));
    }

  if ( GribGridType(sec2) == REP_SPECTRAL )
    {
      GribNumValues(sec4) = DimSPOut;
    }
  else if ( GribGridType(sec2) == REP_GAUSS )
    {
      GribNumValues(sec4) = DimGPOut;
    }
  else
    {
      Error(func, "Unsupported Grid type: %d", GribGridType(sec2));
    }

  gribexdp(sec0, sec1, sec2, psec2, sec3,
	   psec3, sec4, field, klenp, (void *) grib,
	   kleng, &kword, "C", &kret);

  *griblen = kword * sizeof(int);
}

void
geninz (int *interpolation_index, double *pres_of_height,
        double *full_level_pressure, int dimgp, int nrql, int levels)
{
   int h,k,l;
   int   *nx;
   double *pressure,*pf;

   pressure = pres_of_height;
   nx       = interpolation_index;
   IntZero (nx, dimgp * nrql);

   for (k = 0; k < nrql; k++) {
      pf       = full_level_pressure;
      for (l = 0; l < levels; l++)
      for (h = 0; h < dimgp ; h++) {
         if (pressure[h] > *pf) nx[h] = l;
         pf++;
      }
      nx       += dimgp;
      pressure += dimgp;
   }
}

void
genind (int *interpolation_index, double lv[],
        double *full_level_pressure, int dimgp, int nrql, int levels)
{
   int  h, k, l;
   int *nx;
   double pressure, *pf;

   nx = interpolation_index;
   IntZero (nx, dimgp * nrql);

   for (k = 0; k < nrql; k++) {
      pressure = lv[k];
      pf       = full_level_pressure;
      for (l = 0; l < levels; l++)
      for (h = 0; h < dimgp ; h++) {
         if (pressure > *pf) nx[h] = l;
         pf++;
      }
      nx += dimgp;
   }
}

void
Extrap (double *slp, double *aph, double *apf,
        double *Geopotential, double *t, int nhor)
{
   double alpha, tstar, tmsl, zprt, zprtal;
   double zrg;
   double zlapse = 0.0065;
   int j;
   extern double RD;
   extern double Grav;

   zrg    = 1.0 / Grav;

   for (j = 0; j < nhor; ++j) {
      if (Geopotential[j] < 0.0001 && Geopotential[j] > -0.0001) slp[j] = aph[j];
      else {
         alpha = RD * zlapse * zrg;
         tstar = (1.0 + alpha * (aph[j]/apf[j] - 1.0)) * t[j];
         if (tstar < 255.0) tstar = 0.5 * (255.0 + tstar);
         tmsl = tstar + zlapse * zrg * Geopotential[j];
         if (tmsl > 290.5 && tstar > 290.5) {
            tstar = 0.5 * (290.5 + tstar);
            tmsl  = tstar;
         }
         if (tmsl-tstar < 0.000001 && tstar-tmsl < 0.000001) alpha = 0.0;
         else if (Geopotential[j] > 0.0001 || Geopotential[j] < -0.0001)
            alpha = RD * (tmsl-tstar) / Geopotential[j];
         zprt   = Geopotential[j] / (RD * tstar);
         zprtal = zprt * alpha;
         slp[j] = aph[j] * exp(zprt*(1.0-zprtal*(0.5-zprtal/3.0)));
      }
   }
}

double
ExtraT (double PRES, double APH, double APF, double GEOS, double T)
{
   double tstar, ztsz, Z1, ZTMSL, ZALPH, PEVAL, ZHTS, ZALP;
   double zrg;
   double zlapse = 0.0065;
   extern double RD;
   extern double Grav;

   zrg    = 1.0 / Grav;
   tstar  = (1.0 + zlapse * RD * zrg * (APH/APF - 1.0)) * T;
   ztsz   = tstar;
   Z1     = tstar + zlapse * zrg * GEOS;

   if (tstar < 255.0) tstar = 0.5 * (255.0 + tstar);
   ZTMSL = tstar + zlapse * zrg * GEOS;

   if (ZTMSL > 290.5 && tstar > 290.5) {
      tstar = 0.5 * (290.5 + tstar);
      ZTMSL  = tstar;
   }
   if (ZTMSL > 290.5 && tstar <= 290.5) ZTMSL=290.5;
   ZALPH=RD*zlapse*zrg;

   if ( ZTMSL-tstar < 0.000001 && tstar-ZTMSL < 0.000001) ZALPH=0.0;
   if ((ZTMSL-tstar > 0.000001 || tstar-ZTMSL > 0.000001) &&
       (GEOS > 0.0001 || GEOS < -0.0001))
      ZALPH=RD*(ZTMSL-tstar)/GEOS;
   if (PRES <= APH)
      PEVAL = ((APH-PRES)*T+ (PRES-APF)*tstar)/ (APH-APF);
   else {
      ZTMSL = Z1;
      tstar = ztsz;
      ZHTS  = GEOS * zrg;
      if (ZHTS > 2000. && Z1 > 298.) {
         ZTMSL=298.;
         if (ZHTS < 2500.) ZTMSL=0.002*((2500.-ZHTS)*Z1+(ZHTS-2000.)*ZTMSL);
      }
      if ((ZTMSL-tstar) < 0.000001)
         ZALPH = 0.;
      else if (GEOS > 0.0001 || GEOS < -0.0001)
         ZALPH = RD*(ZTMSL-tstar)/GEOS;
      else
         ZALPH = RD*zlapse*zrg;

      ZALP  = ZALPH*log(PRES/APH);
      PEVAL = tstar*(1.0+ZALP*(1.0+ZALP*(0.5+0.16666666667*ZALP)));
   }
   return PEVAL;
}

double
ExtraZ (double pres, double aph, double apf, double Geopotential, double t)
{
   double alpha, tstar, tmsl, zalp, zalpal;
   double zrg;
   double zlapse = 0.0065;
   extern double RD;
   extern double Grav;

   zrg   = 1.0 / Grav;
   alpha = RD * zlapse * zrg;
   tstar = (1.0 + alpha * (aph/apf - 1.0)) * t;
   if (tstar < 255.0) tstar = 0.5 * (255.0 + tstar);
   tmsl = tstar + zlapse * zrg * Geopotential;

   if (tmsl > 290.5 && tstar > 290.5) {
     tstar = 0.5 * (290.5 + tstar);
     tmsl  = tstar;
   }
   if (tmsl > 290.5 && tstar <= 290.5) tmsl = 290.5;

   if(tmsl-tstar < 0.000001 && tstar-tmsl < 0.000001)
      alpha = 0.0;
   else if (Geopotential > 0.0001 || Geopotential < -0.0001)
      alpha = RD * (tmsl-tstar) / Geopotential;

   zalp   = log(pres/aph);
   zalpal = zalp * alpha;

   return ((Geopotential - RD*tstar*zalp*(1.0 + zalpal*(0.5 + zalpal/6.0)))*zrg);
}

void
Interpolate_X (double *gt, double *pt, double *pf, int *nx, double level[], int nrql, int dimgp, int dim3gp)
{
   int lp, i;
   int nl, nh;
   double pres;

   for (lp = 0; lp < nrql; lp++) {
      pres = level[lp];
      for (i = 0; i < dimgp; i++) {
         nl = *nx * dimgp + i;
         nh =  nl + dimgp;
         if (nh >= dim3gp)
            *pt = gt[nl];
         else
            *pt =  gt[nl] + (pres-pf[nl])
                * (gt[nh] - gt[nl]) / (pf[nh] - pf[nl]);
         nx++;
         pt++;
      }
   }
}

#if defined (SX)
#pragma odir switch,-ev
#endif
void
Interpolate_T (double orography[], double gt[], double pt[], double pf[], double ph[], int *nx,
               double level[], int levels, int nrql, int dimgp)
{
   int lp, i;
   int nl, nh;
   double pres;
   extern int mars;

   for (lp = 0; lp < nrql; lp++) {
      pres = level[lp];
#if defined (CRAY)
#pragma _CRI inline ExtraT
#endif
      for (i = 0; i < dimgp; i++) {
         nl = *nx++;
         if (nl < 0) pt[lp*dimgp+i] = gt[i];
         else {
            if (nl > levels-2) {
               if (mars) pt[lp*dimgp+i] = gt[(levels-1)*dimgp+i];
               else      pt[lp*dimgp+i] = ExtraT (pres,
                                                  ph[levels*dimgp+i],
                                                  pf[(levels-1)*dimgp+i], orography[i],
                                                  gt[(levels-1)*dimgp+i]);
            }
            else {
               nh = nl + 1;
               pt[lp*dimgp+i] =  gt[nl*dimgp+i] + (pres-pf[nl*dimgp+i])
                              * (gt[nh*dimgp+i] - gt[nl*dimgp+i])
                              / (pf[nh*dimgp+i] - pf[nl*dimgp+i]);
            }
         }
      }
   }
}

void
Interpolate_Z (double orography[], double gz[], double pz[], double pf[], double ph[], int *nx, double gt[],
               double level[], int levels, int nrql, int dimgp)
{
   int lp, i;
   int nl,nh;
   double pres;
   extern int mars;

   for (lp = 0; lp < nrql; lp++) {
      pres = level[lp];
#if defined (CRAY)
#pragma _CRI inline ExtraZ
#endif
      for (i = 0; i < dimgp; i++) {
         nl = *nx++;
         if (pres > ph[(nl+1)*dimgp+i]) nl++;
         if (nl < 0) pz[lp*dimgp+i] = gz[i];
         else {
            if (nl > levels-1) {
               if (mars) pz[lp*dimgp+i] = gt[(levels-1)*dimgp+i];
               else      pz[lp*dimgp+i] = ExtraZ (pres,ph[levels*dimgp+i],
                                                  pf[(levels-1)*dimgp+i], orography[i],
                                                  gt[(levels-1)*dimgp+i]);
            }
            else {
               nh = nl + 1;
               pz[lp*dimgp+i] =  gz[nl*dimgp+i] + (pres-ph[nl*dimgp+i])
                              * (gz[nh*dimgp+i] - gz[nl*dimgp+i])
                              / (ph[nh*dimgp+i] - ph[nl*dimgp+i]);
            }
         }
      }
   }
}

#if defined (SX)
#pragma odir switch,-dv
#endif
void
Interpolate_X_Z (double *gt, double *pt, double *pf, int *nx, double *poh,
                 int nrql, int dimgp, int dim3gp)
{
   int lp, i;
   int nl, nh;
   double pres;

   for (lp = 0; lp < nrql; lp++) {
      for (i = 0; i < dimgp; i++) {
         nl   = *nx * dimgp + i;
         nh   =  nl + dimgp;
         pres = *poh++;
         if (nh >= dim3gp)
            *pt = gt[nl];
         else
            *pt = gt[nl] + (pres-pf[nl]) * (gt[nh] - gt[nl]) / (pf[nh] - pf[nl]);
         nx++;
         pt++;
      }
   }
}

void
Interpolate_T_Z (double orography[], double gt[], double pt[], double pf[], double ph[],
                 int *nx, double *poh, int levels, int nrql, int dimgp)
{
   int lp, i;
   int nl, nh;
   double pres;

   for (lp = 0; lp < nrql; lp++) {
      for (i = 0; i < dimgp; i++) {
         nl   = *nx++;
         pres = *poh++;
         if (nl < 0) pt[lp*dimgp+i] = gt[i];
         else {
            if (nl > levels-2)
               pt[lp*dimgp+i] = ExtraT (pres,ph[levels*dimgp+i],pf[(levels-1)*dimgp+i],
                                        orography[i],gt[(levels-1)*dimgp+i]);
            else {
               nh = nl + 1;
               pt[lp*dimgp+i] =  gt[nl*dimgp+i] + (pres-pf[nl*dimgp+i])
                              * (gt[nh*dimgp+i] - gt[nl*dimgp+i])
                              / (pf[nh*dimgp+i] - pf[nl*dimgp+i]);
            }
         }
      }
   }
}

void
Interpolate_Z_Z (double orography[], double gz[], double pz[], double pf[], double ph[],
                 int *nx, double *poh, double gt[], int levels, int nrql, int dimgp)
{
   int lp, i;
   int nl, nh;
   double pres;

   for (lp = 0; lp < nrql; lp++) {
      for (i = 0; i < dimgp; i++) {
         nl   = *nx++;
         pres = *poh++;
         if (pres > ph[(nl+1)*dimgp+i]) nl++;
         if (nl < 0) pz[lp*dimgp+i] = gz[i];
         else {
            if (nl > levels-1)
               pz[lp*dimgp+i] = ExtraZ (pres,ph[levels*dimgp+i],pf[(levels-1)*dimgp+i],
                                        orography[i],gt[(levels-1)*dimgp+i]);
            else {
               nh = nl + 1;
               pz[lp*dimgp+i] =  gz[nl*dimgp+i] + (pres-ph[nl*dimgp+i])
                              * (gz[nh*dimgp+i] - gz[nl*dimgp+i])
                              / (ph[nh*dimgp+i] - ph[nl*dimgp+i]);
            }
         }
      }
   }
}

#define SCALEHEIGHT (-7000.)
#define SCALESLP 101325.0
#define RCPD  (3.5 * RD)

extern double PlanetRadius;

double ScaleHeight = SCALEHEIGHT;
double ScaleSLP    = SCALESLP;

#define SQUARE_RADIUS (-PlanetRadius * PlanetRadius)


void IniQuaSum(double Destination[], double Source[], int Length)
{
   int i;
   for (i = 0; i<Length; i++) Destination[i] = Source[i] * Source[i];
}

void AddQuaSum(double Destination[], double Source[], int Length)
{
   int i;
   for (i = 0; i<Length; i++) Destination[i] += Source[i] * Source[i];
}

void VarQuaSum(double Variance[], double Sum[], int Length, int n)
{
   int i;
   double rn1;

   rn1 = 1.0 / (n-1);

   for (i = 0; i<Length; i++)
      Variance[i] = (Variance[i] - Sum[i] * Sum[i] * n) * rn1;
   for (i = 0; i<Length; i++)
      if (Variance[i] > 0.0) Variance[i] = sqrt(Variance[i]);
      else                   Variance[i] = 0.0;
}

void AddVector(double Destination[], double Source[], int Length)
{
   int i;
   for (i = 0; i<Length; i++) Destination[i] += Source[i];
}

void MultVectorScalar(double Destination[], double Source[],
                      double Factor, int Length)
{
   int i;
   for (i = 0; i<Length; i++) Destination[i] = Source[i] * Factor;
}

void Add2Vectors(double *Destination, double *SourceA,
                 double *SourceB, int Length)
{
   int i;
   for (i = 0; i<Length; i++)
       Destination[i] = SourceA[i] + SourceB[i];
}

void Sub2Vectors(double *Destination, double *SourceA,
                 double *SourceB, int Length)
{
   int i;
   for (i = 0; i<Length; i++)
       Destination[i] = SourceA[i] - SourceB[i];
}


void Speed(double *speed, double *u, double *v, int dim3gp)
{
   int i;

   for (i = 0; i < dim3gp; i++)
   speed[i] = sqrt(u[i] * u[i] + v[i] * v[i]);
}

void dv2ps(double *div, double *pot, int lev, int trunc)
{
   int l,m,n;

   for (l = 0; l <  lev  ; l++)
   for (m = 0; m <= trunc; m++)
   for (n = m; n <= trunc; n++) {
      if (n > 0) {
         *pot++ = *div++ * SQUARE_RADIUS / (n * n + n);
         *pot++ = *div++ * SQUARE_RADIUS / (n * n + n);
      }  else {
         *pot++ = 0.0;
         *pot++ = 0.0;
         div += 2;
      }
   }
}

void dv2uv(double *d, double *o, double *u, double *v, double *f, double *g,
           int nt, int nsp, int nlev)
{
   /* d(nsp,nlev),o(nsp,nlev)     ! divergence , vorticity       */
   /* u(nsp,nlev),v(nsp,nlev)     ! zonal wind , meridional wind */
   /* f(nsp/2)   ,g(nsp/2)        ! factor tables                */

   int l,m,n;
   int i;

   for (l = 0; l < nlev; l++) {
      i = 0;
      for (m = 0; m < nt-1; m++) {
         /*********/
         /* n = m */
         /*********/

         if (m == 0) {
            *u++ = -g[i+1] * o[2*(i+1)  ];
            *u++ = -g[i+1] * o[2*(i+1)+1];
            *v++ =  g[i+1] * d[2*(i+1)  ];
            *v++ =  g[i+1] * d[2*(i+1)+1];
         }  else {
            *u++ = -f[i] * d[2*i+1] - g[i+1] * o[2*(i+1)  ];
            *u++ =  f[i] * d[2*i  ] - g[i+1] * o[2*(i+1)+1];
            *v++ = -f[i] * o[2*i+1] + g[i+1] * d[2*(i+1)  ];
            *v++ =  f[i] * o[2*i  ] + g[i+1] * d[2*(i+1)+1];
         }
         ++i;

         /****************/
         /* m < n < nt-1 */
         /****************/

         for (n = m+1; n < nt-1; n++) {
            *u++ =  g[i] * o[2*(i-1)  ] - f[i] * d[2*i+1] - g[i+1] * o[2*(i+1)  ];
            *u++ =  g[i] * o[2*(i-1)+1] + f[i] * d[2*i  ] - g[i+1] * o[2*(i+1)+1];
            *v++ = -g[i] * d[2*(i-1)  ] - f[i] * o[2*i+1] + g[i+1] * d[2*(i+1)  ];
            *v++ = -g[i] * d[2*(i-1)+1] + f[i] * o[2*i  ] + g[i+1] * d[2*(i+1)+1];
            ++i;
         }

         /************/
         /* n = nt-1 */
         /************/

         *u++ =  g[i] * o[2*(i-1)  ] - f[i] * d[2*i+1];
         *u++ =  g[i] * o[2*(i-1)+1] + f[i] * d[2*i  ];
         *v++ = -g[i] * d[2*(i-1)  ] - f[i] * o[2*i+1];
         *v++ = -g[i] * d[2*(i-1)+1] + f[i] * o[2*i  ];
         ++i;

         /**********/
         /* n = nt */
         /**********/

         *u++ =  g[i] * o[2*(i-1)  ];
         *u++ =  g[i] * o[2*(i-1)+1];
         *v++ = -g[i] * d[2*(i-1)  ];
         *v++ = -g[i] * d[2*(i-1)+1];
         ++i;
      }

      /***************************/
      /* m = nt-1  and  n = nt-1 */
      /***************************/

      *u++ = -f[i] * d[2*i+1];
      *u++ =  f[i] * d[2*i  ];
      *v++ = -f[i] * o[2*i+1];
      *v++ =  f[i] * o[2*i  ];
      ++i;

      /*************************/
      /* m = nt-1  and  n = nt */
      /*************************/

      *u++ =  g[i] * o[2*(i-1)  ];
      *u++ =  g[i] * o[2*(i-1)+1];
      *v++ = -g[i] * d[2*(i-1)  ];
      *v++ = -g[i] * d[2*(i-1)+1];
      ++i;

      /***********************/
      /* m = nt  and  n = nt */
      /***********************/

      *u++ = 0.0;
      *u++ = 0.0;
      *v++ = 0.0;
      *v++ = 0.0;

      d += nsp;
      o += nsp;
   }
}

void scaluv(double *fu, double rclat[], int nlat, int lot)
{
   int l,lat;

   for (l = 0; l < lot; l++)
   for (lat = 0; lat < nlat; lat++) {
      *fu *= rclat[lat];
      fu++;
   }
}

void uv2dv(double *fu, double *fv, double *sd, double *sv,
           double *pol2, double *pol3, int klev, int nlat, int nt)
{
   int lev,jmm,jfc,lat;
   double dir,dii,vor,voi;
   double *ufr,*ufi,*vfr,*vfi;
   double *po2,*po3;

   for (lev = 0; lev < klev; lev++) {
      po2 = pol2;
      po3 = pol3;
      for (jmm = 0; jmm <= nt; jmm++) {
         for (jfc = jmm; jfc <= nt; jfc++) {
            ufr = fu        ;
            ufi = fu + nlat ;
            vfr = fv        ;
            vfi = fv + nlat ;
            dir = 0.0       ;
            dii = 0.0       ;
            vor = 0.0       ;
            voi = 0.0       ;
            for (lat = 0; lat < nlat; lat++) {
               dir += *vfr * *po2 - *ufi * *po3;
               dii += *vfi * *po2 + *ufr * *po3;
               vor -= *ufr * *po2 + *vfi * *po3;
               voi -= *ufi * *po2 - *vfr * *po3;
               ufr++;
               ufi++;
               vfr++;
               vfi++;
               po2++;
               po3++;
            }
            *sd++ = dir;
            *sd++ = dii;
            *sv++ = vor;
            *sv++ = voi;
         }
         fu += 2 * nlat;
         fv += 2 * nlat;
      }
   }
}


void theta(double *PThetaF, double *PThetaH, double *PH, double *PS,
           double *TF, double *TS, int Levels, int DimGP, int Dim3GP)
{
   int h,l;
   double  Kappa;
   double *ThetaH = PThetaH;
   double *ThetaF = PThetaF;
   extern double RD;

   Kappa = RD / RCPD;

   for (h = 0; h < DimGP; h++) ThetaH[h] = 0.0;
   ThetaH += DimGP;
   for (l = 0; l < Levels - 1; l++) {
      for (h = 0; h < DimGP; h++) {
         ThetaH[h] = 0.5 * (TF[h] + TF[h+DimGP]) * pow((PS[h]/PH[h]),Kappa);
      }
      PH += DimGP;
      TF += DimGP;
      ThetaH += DimGP;
   }
   RealCopy(ThetaH,TS,DimGP);
   ThetaH = PThetaH;
   for (h = 0; h < Dim3GP; h++) {
      ThetaF[h] = 0.5 * (ThetaH[h] + ThetaH[h+DimGP]);
   }
}

void presh(double *pf, double *php, double *vct, double *ps,
           int Levels, int DimGP, int Dim3GP)
{
   int h,l;
   double zp,ze;
   double *ph = php;

   for (l = 0; l<Levels; l++) {
      zp = vct[l];
      ze = vct[l+Levels+1];
      for (h = 0; h<DimGP; h++) ph[h] = zp + ze * ps[h];
      ph += DimGP;
   }
   RealCopy(ph,ps,DimGP);
   ph = php;
   for (h = 0; h<Dim3GP; h++) pf[h] = 0.5 * (ph[h] + ph[h+DimGP]);
}

void geninx(int nt, double *f, double *g)
{
   int m2,n2;
   int m, n ;

   for (m = 0; m <= nt; m++) {
      m2 = m * m;
      for (n = m; n <= nt; n++) {
         n2 = n * n;
         if (n) {
            *g++ = -PlanetRadius / n * sqrt((double)(n2-m2)/(double)(4*n2-1));
            *f++ = -PlanetRadius * m / (double)(n2+n);
         }  else {
            *g++ = 0.0;
            *f++ = 0.0;
         }
      }
   }
}

/* ================================================ */
/* This new version can handle as well (tested) up  */
/* to truncations of T639.                          */
/* ================================================ */

void gauaw (double pa[], double pw[], int nlat)
{
  /*
   * Compute Gaussian latitudes.  On return pa contains the
   * sine of the latitudes starting closest to the north pole and going
   * toward the south
   *
   */

  const int itemax = 20;

  int isym, iter, ins2, jn, j;
  double za, zw, zan;
  double z, zk, zkm1, zkm2, zx, zxn, zldn, zmod;

  /*
   * Perform the Newton loop
   * Find 0 of Legendre polynomial with Newton loop
   */

  ins2 = nlat/2+nlat%2;

  for (j = 0; j < ins2; j++) {
    z = (double)(4*(j+1)-1)*M_PI/(double)(4*nlat+2);
    pa[j] = cos(z+1.0/(tan(z)*(double)(8*nlat*nlat)));
  }

  for (j = 0; j < ins2; j++) {

    za = pa[j];

    iter = 0;
    do {
      iter++;
      zk = 0.0;

      /* Newton iteration step */

      zkm2 = 1.0;
      zkm1 = za;
      zx = za;
      for (jn = 2; jn <= nlat; jn++) {
	zk = ((double)(2*jn-1)*zx*zkm1-(double)(jn-1)*zkm2)/(double)(jn);
	zkm2 = zkm1;
	zkm1 = zk;
      }
      zkm1 = zkm2;
      zldn = ((double)(nlat)*(zkm1-zx*zk))/(1.-zx*zx);
      zmod = -zk/zldn;
      zxn = zx+zmod;
      zan = zxn;

      /* computes weight */

      zkm2 = 1.0;
      zkm1 = zxn;
      zx = zxn;
      for (jn = 2; jn <= nlat; jn++) {
	zk = ((double)(2*jn-1)*zx*zkm1-(double)(jn-1)*zkm2)/(double)(jn);
	zkm2 = zkm1;
	zkm1 = zk;
      }
      zkm1 = zkm2;
      zw = (1.0-zx*zx)/((double)(nlat*nlat)*zkm1*zkm1);
      za = zan;
    } while (iter <= itemax && fabs(zmod) >= DBL_EPSILON);

    pa[j] = zan;
    pw[j] = 2.0*zw;
  }

  for (j = 0; j < nlat/2; j++) {
    isym = nlat-(j+1);
    pa[isym] = -pa[j];
    pw[isym] =  pw[j];
  }

  return;
}

/* ============================================= */
/* phcs - Compute values of Legendre polynomials */
/*        and their meridional derivatives       */
/* ============================================= */

void phcs(double *PNM, double *HNM, int Waves, double PMU,
          double *ZTEMP1, double *ZTEMP2)
{
   int TwoWaves;

   int JK;
   int JN;
   int JM;

   double JNmJK;
   double ZCOS2;
   double Lat;
   double ZAN;
   double ZSINPAR;
   double ZCOSPAR;
   double ZSQP;
   double ZCOSFAK;
   double ZSINFAK;
   double ZQ;
   double ZWM2;
   double ZW;
   double ZWQ;
   double ZQ2M1;
   double ZWM2Q2;
   double Z2Q2;
   double ZCNM;
   double ZDNM;
   double ZENM;

   TwoWaves  = Waves << 1;
   ZCOS2     = sqrt(1.0 - PMU * PMU);
   Lat       = acos(PMU);
   ZAN       = 1.0;
   ZTEMP1[0] = 0.5;

   for (JN = 1; JN < TwoWaves; JN++) {
      ZSQP    = 1.0 / sqrt((double)(JN+JN*JN));
      ZAN    *= sqrt(1.0 - 1.0 / (4 * JN * JN));

      ZCOSPAR = cos(Lat * JN);
      ZSINPAR = sin(Lat * JN) * JN * ZSQP;
      ZCOSFAK = 1.0;

      for (JK = 2; JK < JN; JK += 2) {
         JNmJK = JN - JK;
         ZCOSFAK *= (JK-1) * (JN+JNmJK+2) / (JK * (JN+JNmJK+1));
         ZSINFAK  = ZCOSFAK * (JNmJK) * ZSQP;
         ZCOSPAR += ZCOSFAK * cos(Lat * JNmJK);
         ZSINPAR += ZSINFAK * sin(Lat * JNmJK);
      }

      /*  Code for JK == JN */

      if ((JN & 1) == 0) {
      ZCOSFAK *= (double)((JN-1) * (JN+2)) / (double)(JN * (JN+1));
      ZCOSPAR += ZCOSFAK * 0.5;
      }
      ZTEMP1[JN  ] = ZAN * ZCOSPAR;
      ZTEMP2[JN-1] = ZAN * ZSINPAR;
   }

   memcpy(PNM,ZTEMP1,Waves * sizeof(double));
   PNM += Waves;
   memcpy(PNM,ZTEMP2,Waves * sizeof(double));
   PNM += Waves;

   HNM[0] = 0.0;
   for (JN = 1; JN < Waves; JN++) HNM[JN] =
      JN * (PMU * ZTEMP1[JN] - sqrt((JN+JN+1.0) / (JN+JN-1.0)) * ZTEMP1[JN-1]);
   HNM += Waves;

   HNM[0] = PMU * ZTEMP2[0];
   for (JN = 1; JN < Waves; JN++)
      HNM[JN] = (JN+1) * PMU * ZTEMP2[JN]
              - sqrt((double)((JN+JN+3) * ((JN+1) * (JN+1) - 1))
              / (double)(JN+JN+1)) * ZTEMP2[JN-1];
   HNM += Waves;

   for (JM = 2; JM < Waves; JM++) {
      PNM[0] = sqrt(1.0 + 1.0 / (JM+JM)) * ZCOS2 * ZTEMP2[0];
      HNM[0] = JM * PMU * PNM[0];
#if defined (CRAY)
#pragma _CRI novector
#endif
#if defined (SX)
#pragma vdir novector
#endif
#if defined (__uxp__)
#pragma loop scalar
#endif
      for (JN = 1; JN < TwoWaves-JM; JN++) {
          ZQ      = JM + JM + JN - 1;
          ZWM2    = ZQ+JN;
          ZW      = ZWM2+2;
          ZWQ     = ZW*ZQ;
          ZQ2M1   = ZQ*ZQ-1.;
          ZWM2Q2  = ZWM2*ZQ2M1;
          Z2Q2    = ZQ2M1*2;
          ZCNM    = sqrt((ZWQ*(ZQ-2.))/(ZWM2Q2-Z2Q2));
          ZDNM    = sqrt((ZWQ*(JN+1.))/ZWM2Q2);
          ZENM    = sqrt(ZW * JN /((ZQ+1.0) * ZWM2));
          PNM[JN] = ZCNM * ZTEMP1[JN] - PMU
                  * (ZDNM * ZTEMP1[JN+1] - ZENM * PNM[JN-1]);
          HNM[JN] = (JM + JN) * PMU * PNM[JN]
                  - sqrt(ZW * JN * (ZQ+1) / ZWM2) * PNM[JN-1];
      }
      memcpy(ZTEMP1,ZTEMP2,TwoWaves * sizeof(double));
      memcpy(ZTEMP2,PNM   ,TwoWaves * sizeof(double));
      PNM += Waves;
      HNM += Waves;
   }
}


/* HUMTEST */

/* ************************************** */
/* Thermodynamical constants adopted from */
/* ECMWF IFS-Code                         */
/* ************************************** */

#define RKBOL (1.380658e-23)
#define RNAVO (6.0221367e+23)
#define R     (RKBOL * RNAVO)
#define RMD   (28.9644)
#define RMV   (18.0153)
#define EARTH_RD    (1000. * R / RMD)
#define RV    (1000. * R / RMV)
#define RCPV  (4.0 * RV)
#define RETV  (RV / RD - 1.)
#define RCW   (4218.)
#define RCS   (2106.)
#define RTT   (273.16)
#define RLVTT (2.5008e+6)
#define RLSTT (2.8345e+6)
#define RESTT (611.14)

void sh2rh(double *sphum, double *rhum, double *t, int lev,
            int dimgpout, double *level, double *fullpresshybrid)
{
   int lp,i;
   int lpi,lfp;
   double  es,qsat;
   double *fullp;
   double  RALPW, RBETW, RGAMW;
   double  RALPS, RBETS, RGAMS;
   double  RALP , RBET , RGAM ;
   extern int    AnalysisData;
   extern double RD;

/* ***************************************************** */
/* Define constants for calculation in presence of water */
/* ***************************************************** */
   RGAMW = (RCW - RCPV) / RV;
   RBETW = RLVTT / RV + RGAMW * RTT;
   RALPW = log(RESTT) + RBETW / RTT + RGAMW * log(RTT);

/* ***************************************************** */
/* Define constants for calculation in presence of  ice  */
/* ***************************************************** */
   RGAMS = (RCS - RCPV) / RV;
   RBETS = RLSTT / RV + RGAMS * RTT;
   RALPS = log(RESTT) + RBETS / RTT + RGAMS * log(RTT);

   if (AnalysisData) fullp = level;
   else              fullp = fullpresshybrid;

/***************************************************/
/* Diagnostics of saturation water vapour pressure */
/* over ice makes no sense, therefore ...          */
/* Hint of Michael Ponater                08.10.97 */
/***************************************************/

   RGAM = RGAMW; RBET = RBETW; RALP = RALPW;
   for (lp = 0; lp < lev; lp++) {
      for (i = 0; i < dimgpout; i++) {
         lpi = lp*dimgpout + i;
         lfp = (1 - AnalysisData) * (lp*dimgpout + i) +
                    AnalysisData  *  lp;
/*       if (t[lpi] < RTT) { */
/*          RGAM = RGAMS; RBET = RBETS; RALP = RALPS; */
/*       }  else { */
/*          RGAM = RGAMW; RBET = RBETW; RALP = RALPW; */
/*       } */
         es = (exp(RALP - RBET / t[lpi] - RGAM * log(t[lpi]))) /
               fullp[lfp];
         qsat = es / (1. + RETV * (1. - es));
         rhum[lpi] = sphum[lpi] * 100. / qsat;
      }
   }
}

void rh2sh(double *sphum, double *rhum, double *t, int lev,
           int dimgpout, double *level)
{
   int lp,i;
   int lpi;
   double  es,qsat;
   double  RALPW, RBETW, RGAMW;
   double  RALPS, RBETS, RGAMS;
   double  RALP , RBET , RGAM ;
   extern double RD;

/* ***************************************************** */
/* Define constants for calculation in presence of water */
/* ***************************************************** */
   RGAMW = (RCW - RCPV) / RV;
   RBETW = RLVTT / RV + RGAMW * RTT;
   RALPW = log(RESTT) + RBETW / RTT + RGAMW * log(RTT);

/* ***************************************************** */
/* Define constants for calculation in presence of  ice  */
/* ***************************************************** */
   RGAMS = (RCS - RCPV) / RV;
   RBETS = RLSTT / RV + RGAMS * RTT;
   RALPS = log(RESTT) + RBETS / RTT + RGAMS * log(RTT);

/***************************************************/
/* Diagnostics of saturation water vapour pressure */
/* over ice makes no sense, therefore ...          */
/* Hint of Michael Ponater                08.10.97 */
/***************************************************/

   RGAM = RGAMW; RBET = RBETW; RALP = RALPW;
   for (lp = 0; lp < lev; lp++) {
      for (i = 0; i < dimgpout; i++) {
         lpi = lp*dimgpout + i;
/*       if (t[lpi] < RTT) { */
/*          RGAM = RGAMS; RBET = RBETS; RALP = RALPS; */
/*       }  else { */
/*          RGAM = RGAMW; RBET = RBETW; RALP = RALPW; */
/*       } */
         es = (exp(RALP - RBET / t[lpi] - RGAM * log(t[lpi]))) /
               level[lp];
         qsat = es / (1. + RETV * (1. - es));
         sphum[lpi] = rhum[lpi] * qsat /  100.;
      }
   }
}

/* HUMTEST ENDE */


void MakeGeopotHeight(double *geop, double* gt, double *gq, double *ph, int nhor, int nlev)
{
   int i, j;
   double VTMP;
   double zrg;
   extern double Grav;
   extern double RD;

   VTMP = (RV / RD) - 1.0;
   zrg  = 1.0 / Grav;

   if (gq) /* Humidity is present */ {
      for ( j = nlev ; j > 1 ; j-- ) 
#if defined (SX)
#pragma vdir nodep
#endif
      for ( i = nhor * (j-1) ; i < nhor * j ; i++ )
         geop[i] = geop[i+nhor] + RD * gt[i] * (1.0 + VTMP * gq[i])
                 * log(ph[i+nhor] / ph[i]);

#if defined (SX)
#pragma vdir nodep
#endif
      for (i = 0; i < nhor; i++)
         geop[i] = geop[i+nhor] + RD * gt[i] * (1.0 + VTMP * gq[i])
                 * 2.0 * log(2.0);
   }
   else    /* No humidity */ {
      for ( j = nlev ; j > 1 ; j-- ) 
#if defined (SX)
#pragma vdir nodep
#endif
      for ( i = nhor * (j-1) ; i < nhor * j ; i++ )
         geop[i] = geop[i+nhor] + RD * gt[i] * log(ph[i+nhor] / ph[i]);

#if defined (SX)
#pragma vdir nodep
#endif
      for (i = 0; i < nhor; i++)
         geop[i] = geop[i+nhor] + RD * gt[i] * 2.0 * log(2.0);
   }

#if defined (SX)
#pragma vdir nodep
#endif
   for (i = 0; i < nhor * (nlev+1); i++) geop[i] *= zrg;
}

/* ======================================== */
/* LayerWater integral liquid water content */
/* ======================================== */

#define MAX_LEVELS      99

void LayerWater (double *ww, double *ll, double pmax, double pmin,
                 int DimGP, int HalfLevels, double *vct)
{
   int  i,k;
   int  MaxLev, MinLev;
   double pph[MAX_LEVELS];
   double pps = SCALESLP;
   extern double Grav;

   for (k = 0; k < HalfLevels; k++)
      pph[k] = vct[k] + vct[k+HalfLevels] * pps;
   for (k = 0; k < HalfLevels; k++)
      if (pph[k] > pmax) break;
   MaxLev = k - 1;
   for (k = HalfLevels - 1; k >= 0; k--)
      if (pph[k] < pmin) break;
   MinLev = k;

   RealZero (ll, DimGP);

   for (k = MaxLev; k <= MinLev; k++) {
     for (i = 0;     i < DimGP;  i++)
        ll[i] += ww[i+k*DimGP] * (pph[k+1] - pph[k]);
   }
   for (i = 0; i < DimGP; i++) ll[i] /= Grav;
}

/* ================================================= */
/* LayerCloud calculates random overlap cloud cover */
/* ================================================= */

void LayerCloud (double *cc, double *ll, double pmax, double pmin,
                 int DimGP, int HalfLevels, double *vct)
{
   int  i, k;
   int  MaxLev, MinLev;
   double pph[MAX_LEVELS];
   double pps = SCALESLP;
   double ZEPSEC = 1.0e-12;

   for (k = 0; k < HalfLevels; k++)
      pph[k] = vct[k] + vct[k+HalfLevels] * pps;
   for (k = 0; k < HalfLevels; k++)
      if (pph[k] > pmax) break;
   MaxLev = k - 1;
   for (k  =  HalfLevels - 1; k >=0; k--)
      if (pph[k] < pmin) break;
   MinLev = k;

   for (i = 0; i < DimGP; i++) ll[i] = 1. - cc[i+MaxLev*DimGP];

   for (k = MaxLev + 1; k <= MinLev; k++) {
     for (i = 0;     i < DimGP;  i++)
         ll[i] *= (1. - MAX(cc[i+(k-1)*DimGP],cc[i+k*DimGP]))
                / (1. - MIN(cc[i+(k-1)*DimGP],1.-ZEPSEC));
   }
   for (i = 0; i < DimGP; i++) ll[i] = 1. - ll[i];
}

void Derivate(double field[], double derilam[], int levels,
              int Waves, int Latitudes, double DerivationFactor[])
{
   int l, n, lev;
   int i;

   i = 0;
   for (lev = 0; lev < levels; lev++)
   for (n = 0; n < Waves    ; n++) {
     for (l = 0; l < Latitudes; l++) {
       derilam[i] = -n * field[i+Latitudes] * DerivationFactor[l];
       i++;
     }
     for (l = 0; l < Latitudes; l++) {
       derilam[i] =  n * field[i-Latitudes] * DerivationFactor[l];
       i++;
     }
   }
}

void h2p(double *PHeight, double lv[], int DimGP, int nrql)
{
   double exp_arg;
   int  h,k;
   double Height;

   for (k = 0; k<nrql; k++) {
      Height  = (double) lv[k];
/*    unitsel == 1 : lv[k] is given in meters
      unitsel == 2 : lv[k] is given in kilometers
      h2p needs meters (MKSC-standard) */

      exp_arg = (double) (Height / ScaleHeight);
      for (h = 0; h<DimGP ; h++) {
         *PHeight = ScaleSLP * exp(exp_arg);
          PHeight++;
      }
   }
}

#ifndef CRAY
#define FORTRAN_OUTPUT
#endif

void wrifor (double *field, int dim, int head[8], FILE* gp)
{
#ifdef FORTRAN_OUTPUT
   int i;
   int headcontrol, fieldcontrol;  /* FORTRAN control words */
   float *RealField;

   RealField = FloatAlloc(dim,"wrifor.RealField");
   /*   headcontrol = sizeof(head);*/
   headcontrol  = 32;
   fieldcontrol = dim * sizeof(float);

   for (i = 0; i < dim; i++) RealField[i] = field[i];

   (void) fwrite (&headcontrol , sizeof(int),     1, gp);
   (void) fwrite (head         , sizeof(head[0]), 8, gp);
   (void) fwrite (&headcontrol , sizeof(int),     1, gp);
   (void) fwrite (&fieldcontrol, sizeof(int),     1, gp);
   (void) fwrite (RealField,     sizeof(float), dim, gp);
   (void) fwrite (&fieldcontrol, sizeof(int),     1, gp);

   RealField = FreeMem(RealField);
#else
   (void) fwrite ((char *)head , sizeof(head[0] ),   8, gp);
   (void) fwrite ((char *)field, sizeof(field[0]), dim, gp);
#endif
}

void Stars(int n)
{
   while (n--) fputc('*',stderr);
}

void NewLine(void)
{
   fputc('\n',stderr);
}

/* ==================================== */
/* Abort - Print error message and exit */
/* ==================================== */

void Abort(char *errtext)
{
   Stars(MAX(80,strlen(errtext))); NewLine();
   fprintf(stderr,errtext);
   Stars(MAX(80,strlen(errtext))); NewLine();
   exit(1);
}

/* ==================================== */
/* RealCopy - Copy array of type double */
/* ==================================== */

void RealCopy(void *destination, void *source, int words)
{
   memcpy(destination,source,words * sizeof(double));
}

/* ======================================== */
/* IntZero -  Set array of type int to zero */
/* ======================================== */

void IntZero(int *field, int words)
{
   memset((char *)field,0,words * sizeof(int));
}

/* ============================================ */
/* RealZero -  Set array of type double to zero */
/* ============================================ */

void RealZero(double *field, int words)
{
   memset((char *)field,0,words * sizeof(double));
}


void Error(const char *caller, const char *fmt, ...);
void Warning(const char *caller, const char *fmt, ...);
void Message(const char *caller, const char *fmt, ...);


#undef  FALSE
#define FALSE 0
#undef  TRUE
#define TRUE  1

#undef  UCHAR
#define UCHAR unsigned char

#undef  NINT
#define NINT(x) ((x) < 0 ? (int)((x)-.5) : (int)((x)+.5))

int  BitsPerInt = sizeof(int) * 8;
int *IGrib;
static int    z             ; /* Counter of GRIB length for output          */

/*
 * A version string.
 */

#define LIBVERSION      0.01
#define XSTRING(x)	#x
#define STRING(x)	XSTRING(x)
static const char grb_libvers[] = STRING(LIBVERSION) " of "__DATE__" "__TIME__;

const char *
grbLibraryVersion(void)
{
  return (grb_libvers);
}

double
decfp2(int kexp, int kmant)
{
  /*

    Purpose:
    --------

    Convert GRIB representation of a floating point
    number to machine representation.

    Input Parameters:
    -----------------

    kexp    - 8 Bit signed exponent.
    kmant   - 24 Bit mantissa.

    Output Parameters:
    ------------------

    Return value   - Floating point number represented
                     by kexp and kmant.

    Method:
    -------

    Floating point number represented as 8 bit exponent
    and 24 bit mantissa in integer values converted to
    machine floating point format.

    Externals:
    ----------

    None.

    Reference:
    ----------

    WMO Manual on Codes re GRIB representation.

    Comments:
    ---------

    Rewritten from DECFP, to conform to programming standards.
    Sign bit on 0 value now ignored, if present.
    If using 32 bit reals, check power of 16 is not so small as to
    cause overflows (underflows!); this causes warning to be given
    on Fujitsus.

    Author:
    -------

    John Hennessy   ECMWF   18.06.91

    Modifications:
    --------------

    Uwe Schulzweida   MPIfM   01/04/2001

    Convert to C from EMOS library version 130
     
  */

  static char func[] = "decfp2";
  double pval;
  static int iexp, isign;
  extern int GRB_Debug;

  /* ----------------------------------------------------------------- */
  /*   Section 1 . Convert value of 0.0. Ignore sign bit.              */
  /* ----------------------------------------------------------------- */

  if ( GRB_Debug ) Message(func, "KEXP = %d  KMANT = %d", kexp, kmant);
  /*
  if ( (kexp == 128 || kexp == 0) && kmant == 0 )
  */
  if ( (kexp == 128) || (kexp == 0) || (kexp == 255) )
    {
      pval = 0.0;
      goto LABEL900;
    }

  /* ----------------------------------------------------------------- */
  /*   Section 2 . Convert other values.                               */
  /* ----------------------------------------------------------------- */

  /*  Sign of value. */

  iexp = kexp;
  isign = 1;

  if ( iexp >= 128 )
    {
      iexp += -128;
      isign = -1;
    }

  /*  Decode value. */

  pval = isign * pow(2.0, -24.0) * kmant * pow(16.0, (double)(iexp - 64));

  /* ----------------------------------------------------------------- */
  /*   Section 9. Return to calling routine.                           */
  /* ----------------------------------------------------------------- */

LABEL900:

  if ( GRB_Debug ) Message(func, "Returned value = %f", pval);

  return (pval);
} /* decfp2 */

void
confp3(double pval, int *kexp, int *kmant, int kbits, int kround)
{
  /*

    Purpose:
    --------

    Convert floating point number from machine
    representation to GRIB representation.

    Input Parameters:
    -----------------

       pval    - Floating point number to be converted.
       kbits   - Number of bits in computer word.
       kround  - Conversion type.
                 0 , Closest number in GRIB format less than
                     original number.
                 1 , Closest number in GRIB format to the
                     original number (equal to, greater than or
                     less than original number).

    Output Parameters:
    ------------------

       kexp    - 8 Bit signed exponent.
       kmant   - 24 Bit mantissa.

    Method:
    -------

    Floating point number represented as 8 bit signed
    exponent and 24 bit mantissa in integer values.

    Externals.
    ----------

    decfp2    - Decode from IBM floating point format.

    Reference:
    ----------

    WMO Manual on Codes re GRIB representation.

    Comments:
    ---------

    Routine aborts if an invalid conversion type parameter
    is used or if a 24 bit mantissa is not produced.

    Author:
    -------
     
    John Hennessy   ECMWF   18.06.91

    Modifications:
    --------------

    Uwe Schulzweida   MPIfM   01/04/2001

    Convert to C from EMOS library version 130

  */

  static char func[] = "confp3";
  double zval, rpowref;
  static double zref, zeps;
  static int iexp, isign;
  static int iround;
  extern int GRB_Debug;

  /* ----------------------------------------------------------------- */
  /*   Section 1 . Initialise                                          */
  /* ----------------------------------------------------------------- */

  /*  Check conversion type parameter. */

  iround = kround;
  if ( iround != 0 && iround != 1 )
    {
      Error(func,"Invalid conversion type = %d", iround);

      /*  If not aborting, arbitrarily set rounding to 'up'. */
     iround = 1;
    }

  /* ----------------------------------------------------------------- */
  /*   Section 2 . Convert value of zero.                              */
  /* ----------------------------------------------------------------- */

  if ( ! fabs(pval) > 0.0 )
    {
      *kexp  = 0;
      *kmant = 0;
      iexp   = 0;
      isign  = 0;
      goto LABEL900;
    }

  /* ----------------------------------------------------------------- */
  /*   Section 3 . Convert other values.                               */
  /* ----------------------------------------------------------------- */

  zeps = 1.0e-12;
  if ( kbits == 32 ) zeps = 1.0e-8;
  zref = pval;

  /*  Sign of value. */

  isign = 0;
  if ( zref < 0.0 )
    {
      isign = 128;
      zref  = - zref;
    }

  /*  Exponent. */

  iexp = (int) (log(zref) * (1.0 / log(16.0)) + 64.0 + 1.0 + zeps);

  if ( iexp < 0   ) iexp = 0;
  if ( iexp > 127 ) iexp = 127;

  rpowref  = zref / pow(16.0, (double)(iexp - 70));

  /*  Mantissa. */

  if ( iround == 0 )
    {

      /*  Closest number in GRIB format less than original number. */
      /*    Truncate for positive numbers. */
      /*    Round up for negative numbers. */

      if ( isign == 0 )
	*kmant = (int) rpowref;
      else
	*kmant = NINT(rpowref + 0.5);
    }
  else
    {

      /*  Closest number in GRIB format to the original number   */
      /*  (equal to, greater than or less than original number). */

      *kmant = NINT(rpowref);
    }

  /*  Check that mantissa value does not exceed 24 bits. */
  /*  If it does, adjust the exponent upwards and recalculate */
  /*  the mantissa. */
  /*  16777215 = 2**24 - 1 */

  if ( *kmant > 16777215 )
    {

    LABEL350:

      ++iexp;

      /*  Check for exponent overflow during adjustment  */

      if ( iexp > 127 )
	{
          Message(func, "Exponent overflow");
          Message(func, "Original number = %30.20f", pval);
          Message(func, "Sign = %3d, Exponent = %3d, Mantissa = %12d",
		  isign, iexp, kmant);

	  Error(func, "Exponent overflow");

	  /*  If not aborting, arbitrarily set value to zero  */

          Message(func, "Value arbitrarily set to zero.");
          kexp  = 0;
          kmant = 0;
          iexp  = 0;
          isign = 0;
          goto LABEL900;
	}

      if ( iround == 0 )
	{

	  /*  Closest number in GRIB format less than original number. */
	  /*  Truncate for positive numbers. */
	  /*  Round up for negative numbers. */

	  if ( isign == 0 )
	    *kmant = (int) rpowref;
	  else
	    *kmant = NINT(rpowref + 0.5);
	}
      else
	{

	  /*  Closest number in GRIB format to the original number */
	  /*  (equal to, greater or less than original number). */

	  *kmant = NINT(rpowref);
	}

      /*  Repeat calculation (with modified exponent) if still have */
      /*  mantissa overflow. */

      if ( *kmant > 16777215 ) goto LABEL350;
    }

  /*  Add sign bit to exponent. */

  *kexp = iexp + isign;

  /* ----------------------------------------------------------------- */
  /*   Section 9. Return                                               */
  /* ----------------------------------------------------------------- */

LABEL900:

  if ( GRB_Debug )
    {
      Message(func, "Conversion type parameter = %4d", kround);
      Message(func, "Original number = %30.20f", pval);

      zval = decfp2(*kexp, *kmant);

      Message(func, "Converted to      %30.20f", zval);
      Message(func, "Sign = %3d, Exponent = %3d, Mantissa = %12d",
	      isign, iexp, *kmant);
    }

  return;
} /* confp3 */

void
ref2ibm(double *pref, int kbits)
{
  /*

    Purpose:
    --------

    Code and check reference value in IBM format

    Input Parameters:
    -----------------

    pref       - Reference value
    kbits      - Number of bits per computer word.

    Output Parameters:
    ------------------

    pref       - Reference value

    Method:
    -------

    Codes in IBM format, then decides to ensure that reference 
    value used for packing is not different from that stored
    because of packing differences.

    Externals.
    ----------

    confp3    - Encode into IBM floating point format.
    decfp2    - Decode from IBM floating point format.

    Reference:
    ----------

    None.

    Comments:
    --------

    None.

    Author:
    -------

    J.D.Chambers     ECMWF      17:05:94

    Modifications:
    --------------

    Uwe Schulzweida   MPIfM   01/04/2001

    Convert to C from EMOS library version 130

  */

  static char func[] = "ref2ibm";
  static int itrnd;
  static int kexp, kmant;
  static double ztemp, zdumm;
  extern int GRB_Debug;

  /* ----------------------------------------------------------------- */
  /*   Section 1. Convert to and from IBM format.                      */
  /* ----------------------------------------------------------------- */

  /*  Convert floating point reference value to IBM representation. */

  itrnd = 1;
  zdumm = ztemp = *pref;
  confp3(zdumm, &kexp, &kmant, kbits, itrnd);

  if ( kexp == 0 && kmant == 0 ) return;

  /*  Set reference value to that actually stored in the GRIB code. */

  *pref = decfp2(kexp, kmant);

  /*  If the nearest number which can be represented in */
  /*  GRIB format is greater than the reference value,  */
  /*  find the nearest number in GRIB format lower      */
  /*  than the reference value.                         */

  if ( ztemp < *pref )
    {
      /*  Convert floating point to GRIB representation */
      /*  using truncation to ensure that the converted */
      /*  number is smaller than the original one.      */

      itrnd = 0;
      zdumm = *pref = ztemp;
      confp3(zdumm, &kexp, &kmant, kbits, itrnd);

      /*  Set reference value to that stored in the GRIB code. */

      *pref = decfp2(kexp, kmant);

      if ( ztemp < *pref )
	{
	  if ( GRB_Debug )
	    {
	      Message(func, "Reference value error.");
	      Message(func, "Notify Met.Applications Section.");
	      Message(func, "ZTEMP = ", ztemp);
	      Message(func, "PREF = ", pref);
	    }
	  *pref = ztemp;
	}
    }

  return;
} /* ref2ibm */


int scm0(double *pdl, double *pdr, double *pfl, double *pfr, int *klg)
{
   /* System generated locals */
   int i_1;
   double r_1;

   /* Local variables */
   static double zfac, zeps, zbeta;
   static int jl;
   static double zalpha;

/* **** SCM0   - Apply SCM0 limiter to derivative estimates. */
/* output: */
/*   pdl   = the limited derivative at the left edge of the interval */
/*   pdr   = the limited derivative at the right edge of the interval */
/* inputs */
/*   pdl   = the original derivative at the left edge */
/*   pdr   = the original derivative at the right edge */
/*   pfl   = function value at the left edge of the interval */
/*   pfr   = function value at the right edge of the interval */
/*   klg  = number of intervals where the derivatives are limited */

/*  define constants */

   zeps = (double)1e-12;
   zfac = ((double)1. - zeps) * (double)3.;

   i_1 = *klg;
   for (jl = 0; jl < i_1; ++jl) {
      if ((r_1 = pfr[jl] - pfl[jl], fabs(r_1)) > zeps) {
         zalpha = pdl[jl] / (pfr[jl] - pfl[jl]);
         zbeta  = pdr[jl] / (pfr[jl] - pfl[jl]);
         if (zalpha <= (double)0.) pdl[jl] = (double)0.;
         if (zbeta  <= (double)0.) pdr[jl] = (double)0.;
         if (zalpha > zfac)      pdl[jl] = zfac * (pfr[jl] - pfl[jl]);
         if (zbeta  > zfac)      pdr[jl] = zfac * (pfr[jl] - pfl[jl]);
      } else {
         pdl[jl] = (double)0.;
         pdr[jl] = (double)0.;
      }
   }
   return (0);
} /* scm0_ */

int rowina2(double *p, int *ko, int *ki, double *pw,
	    int kcode, double *pmsval, int *kret)
{
   /* System generated locals */
   int pw_dim1, pw_offset, i_1;

   /* Local variables */
   static double zwt1, zrdi, zpos;
   static int jl, ip;
   static double zdo, zwt;

   /* Parameter adjustments */
   --p;
   pw_dim1 = *ko + 3;
   pw_offset = pw_dim1;
   pw -= pw_offset;

/* **** ROWINA2 - Interpolation of row of values. */
/*     Input Parameters. */
/*     ----------------- */
/*     P      - Row of values to be interpolated. */
/*              Dimension must be at least KO. */
/*     KO     - Number of values required. */
/*     KI     - Number of values in P on input. */
/*     PW     - Working array. */
/*              Dimension must be at least (0:KO+2,3). */
/*     KCODE  - Interpolation required. */
/*              1 , linear. */
/*              3 , cubic. */
/*     PMSVAL - Value used for missing data indicator. */

/*     Output Parameters. */
/*     ------------------ */
/*     P     - Now contains KO values. */
/*     KRET  - Return code */
/*             0, OK */
/*             Non-zero, error */

/*     Author. */
/*     ------- */
/*     J.D.Chambers    ECMWF     22.07.94 */

/*     ********************************    */
/*     Section 1.  Linear interpolation .. */
/*     ********************************    */

   *kret = 0;

   if (kcode == 1) {
/*    Move input values to work array */
      i_1 = *ki;
      for (jl = 1; jl <= i_1; ++jl)
          pw[jl + pw_dim1] = p[jl];

/*    Arrange wrap-around value in work array */
      pw[*ki + 1 + pw_dim1] = p[1];

/*    Set up constants to be used to figure out weighting for */
/*    values in interpolation. */
      zrdi = (double) (*ki);
      zdo = 1. / (double) (*ko);

/*    Loop through the output points */
      i_1 = *ko;
      for (jl = 1; jl <= i_1; ++jl) {

/*    Calculate weight from the start of row */
      zpos = (jl - 1) * zdo;
      zwt = zpos * zrdi;

/*    Get the current array position(minus 1) from the weight - */
/*    note the implicit truncation. */
      ip = (int) zwt;

/*    If the left value is missing, use the right value */
      if (pw[ip + 1 + pw_dim1] == *pmsval) {
         p[jl] = pw[ip + 2 + pw_dim1];

/*    If the right value is missing, use the left value */
      } else if (pw[ip + 2 + pw_dim1] == *pmsval) {
         p[jl] = pw[ip + 1 + pw_dim1];

/*    If neither missing, interpolate ... */
      } else {

/*       Adjust the weight to range (0.0 to 1.0) */
         zwt -= ip;

/*       Interpolate using the weighted values on either side */
/*       of the output point position */
         p[jl] = ((double)1. - zwt) * pw[ip + 1 + pw_dim1] +
                 zwt * pw[ip + 2 + pw_dim1];
      }
   }

/*     *******************************    */
/*     Section 2.  Cubic interpolation .. */
/*     *******************************    */

   } else if (kcode == 3) {
      i_1 = *ki;
      for (jl = 1; jl <= i_1; ++jl) {
          if (p[jl] == *pmsval) {
             fprintf(stderr," ROWINA2: ");
             fprintf(stderr," Cubic interpolation not supported");
             fprintf(stderr," for fields containing missing data.\n");
             *kret = 1;
             goto L900;
          }
          pw[jl + pw_dim1] = p[jl];
      }
      pw[pw_dim1] = p[*ki];
      pw[*ki + 1 + pw_dim1] = p[1];
      pw[*ki + 2 + pw_dim1] = p[2];
      i_1 = *ki;
      for (jl = 1; jl <= i_1; ++jl) {
          pw[jl + (pw_dim1 << 1)] =
                - pw[jl - 1 + pw_dim1] / (double)3. -
                  pw[jl + pw_dim1] * (double).5 +
                  pw[jl + 1 +pw_dim1] - pw[jl + 2 + pw_dim1] / (double)6.;
          pw[jl + 1 + pw_dim1 * 3] =
                  pw[jl - 1 + pw_dim1] / (double)6. -
                  pw[jl + pw_dim1] +
                  pw[jl + 1 + pw_dim1] * (double).5 +
                  pw[jl + 2 + pw_dim1] / (double)3.;
      }
      scm0(&pw[(pw_dim1 << 1) + 1], &pw[pw_dim1 * 3 + 2],
	   &pw[pw_dim1 + 1], &pw[pw_dim1 + 2], ki);
      zrdi = (double) (*ki);
      zdo = (double)1. / (double) (*ko);
      i_1 = *ko;
      for (jl = 1; jl <= i_1; ++jl) {
          zpos = (jl - 1) * zdo;
          zwt = zpos * zrdi;
          ip = (int) zwt + 1;
          zwt = zwt + (double)1. - ip;
          zwt1 = (double)1. - zwt;
          p[jl] = (((double)3. - zwt1 * (double)2.) * pw[ip + pw_dim1] +
                  zwt * pw[ip + (pw_dim1 << 1)]) * zwt1 * zwt1 +
                  (((double)3. - zwt * (double)2.) * pw[ip + 1 + pw_dim1] -
                  zwt1 * pw[ip + 1 + pw_dim1 * 3]) * zwt * zwt;
      }

   } else {

/*    **************************************    */
/*    Section 3.  Invalid interpolation code .. */
/*    **************************************    */
      fprintf(stderr," ROWINA2:");
      fprintf(stderr," Invalid interpolation code = %2d\n",kcode);
      *kret = 2;
   }

L900:
    return 0;
} /* rowina2_ */

int qu2reg2(double *pfield, int *kpoint, int *klat, int *klon,
	    double *ztemp, double *pmsval, int *kret)
{
   /* System generated locals */
   int i_1, i_2;
   int kcode = 1;

   /* Local variables */
   static int ilii, ilio, icode;
   static double *zline = NULL;
   static double zwork[1929];   /* was [643][3] */
   static int iregno, iquano, j210, j220, j230, j240, j225;


   zline = (double *) malloc(*klon*sizeof(double));

   /* Parameter adjustments */
   --pfield;
   --kpoint;

/* **** QU2REG - Convert quasi-regular grid data to regular. */
/*     Input Parameters. */
/*     ----------------- */
/*     PFIELD     - Array containing quasi-regular grid */
/*                  data. */
/*     KPOINT     - Array containing list of the number of */
/*                  points on each latitude (or longitude) of */
/*                  the quasi-regular grid. */
/*     KLAT       - Number of latitude lines */
/*     KLON       - Number of longitude lines */
/*     KCODE      - Interpolation required. */
/*                  1 , linear - data quasi-regular on */
/*                               latitude lines. */
/*                  3 , cubic -  data quasi-regular on */
/*                               latitude lines. */
/*                  11, linear - data quasi-regular on */
/*                               longitude lines. */
/*                  13, cubic -  data quasi-regular on */
/*                               longitude lines. */
/*     PMSVAL     - Value used for missing data indicator. */
/*     Output Parameters. */
/*     ------------------ */
/*     KRET       - return code */
/*                  0 = OK */
/*                  non-zero indicates fatal error */
/*     PFIELD     - Array containing regular grid data. */
/*     Author. */
/*     ------- */
/*     J.D.Chambers     ECMWF      22.07.94 */
/*     J.D.Chambers     ECMWF      13.09.94 */
/*     Add return code KRET and remove calls to ABORT. */


/* ------------------------------ */
/* Section 1. Set initial values. */
/* ------------------------------ */

   *kret = 0;

/* Check input parameters. */

   if (kcode != 1 && kcode != 3 && kcode != 11 && kcode != 13) {
      fprintf(stderr," QU2REG :");
      fprintf(stderr," Invalid interpolation type code = %2d\n",kcode);
      *kret = 1;
      goto L900;
   }

/* Set array indices to 0. */

   ilii = 0;
   ilio = 0;

/* Establish values of loop parameters. */

   if (kcode > 10) {

/*    Quasi-regular along longitude lines. */

      iquano = *klon;
      iregno = *klat;
      icode = kcode - 10;
   } else {

/*    Quasi-regular along latitude lines. */

      iquano = *klat;
      iregno = *klon;
      icode = kcode;
   }

/*     -------------------------------------------------------- */
/**    Section 2. Interpolate field from quasi to regular grid. */
/*     -------------------------------------------------------- */

   i_1 = iquano;
   for (j230 = 1; j230 <= i_1; ++j230) {

      if (iregno != kpoint[j230]) {

/*       Line contains less values than required,so */
/*       extract quasi-regular grid values for a line */

         i_2 = kpoint[j230];
         for (j210 = 1; j210 <= i_2; ++j210) {
            ++ilii;
            zline[j210 - 1] = pfield[ilii];
         }

/*       and interpolate this line. */

         rowina2(zline, &iregno, (int *)&kpoint[j230], zwork, icode, pmsval, kret);
         if (*kret != 0) goto L900;

/*       Add regular grid values for this line to the
         temporary array. */

         i_2 = iregno;
         for (j220 = 1; j220 <= i_2; ++j220) {
            ++ilio;
            ztemp[ilio - 1] = zline[j220 - 1];
         }

      } else {

/*       Line contains the required number of values, so add */
/*       this line to the temporary array. */

         i_2 = iregno;
         for (j225 = 1; j225 <= i_2; ++j225) {
            ++ilio;
            ++ilii;
            ztemp[ilio - 1] = pfield[ilii];
         }
      }
   }

/* Copy temporary array to user array. */

   i_1 = *klon * *klat;
   for (j240 = 1; j240 <= i_1; ++j240) {
      pfield[j240] = ztemp[j240 - 1];
   }

/* -------------------------------------------------------- */
/* Section 9. Return to calling routine. Format statements. */
/* -------------------------------------------------------- */

L900:

   free(zline);

   return 0;
} /* qu2reg2_ */






void PutnZero(int n)
{
   int i;

   for (i = z; i < z+n; i++) IGrib[i] = 0;
   z += n;
}

#define Put1Byte(Value)  (IGrib[z++] = (Value))
#define Put2Byte(Value) ((IGrib[z++] = (Value) >>  8), \
                         (IGrib[z++] = (Value)))
#define Put3Byte(Value) ((IGrib[z++] = (Value) >> 16), \
                         (IGrib[z++] = (Value) >>  8), \
                         (IGrib[z++] = (Value)))


void Put1Real(double Value)
{
   int Exponent, Mantissa;
   int one = 1;

   confp3(Value, &Exponent, &Mantissa, BitsPerInt, one);
   Put1Byte (Exponent);
   Put3Byte (Mantissa);
}

/* GRIB block 0 - indicator block */

void PutGribBlock0(void)
{
  IGrib[0] = 'G';
  IGrib[1] = 'R';
  IGrib[2] = 'I';
  IGrib[3] = 'B';

  /* 
   * IGrib[4]-IGrib[6] contains full length of grib record. 
   * included before finished CODEGB
   */

  z = 7;   
  Put1Byte(1); 
  z = 8;
}

/* GRIB block 5 - end block */

void PutGribBlock5(void)
{
  IGrib[z++] = '7';
  IGrib[z++] = '7';
  IGrib[z++] = '7';
  IGrib[z++] = '7';

  IGrib[4] = z >> 16; IGrib[5] = z >>  8; IGrib[6] = z;
  while (z & 7) IGrib[z++] = 0;
}

/* GRIB block 1 - product definition block. */

void PutGribBlock1(int *sec1)
{
  Put3Byte(28);                       /*  0 Length of Block 1        */
  Put1Byte(GribCodeTable(sec1));      /*  3 Local table number       */
  Put1Byte(GribCenterID(sec1));       /*  4 Identification of centre */
  Put1Byte(GribModelID(sec1));        /*  5 Identification of model  */
  Put1Byte(GribGridDefinition(sec1)); /*  6 Grid definition          */
  Put1Byte(GribBlock2Included(sec1)); /*  7 Block 2 included         */
  Put1Byte(GribParameter(sec1));      /*  8 Parameter Code           */
  Put1Byte(GribLevelType(sec1));      /*  9 Type of level            */
  if ( (GribLevelType(sec1) !=  20) &&
       (GribLevelType(sec1) != 100) &&
       (GribLevelType(sec1) != 103) &&
       (GribLevelType(sec1) != 105) &&
       (GribLevelType(sec1) != 107) &&
       (GribLevelType(sec1) != 109) &&
       (GribLevelType(sec1) != 111) &&
       (GribLevelType(sec1) != 113) &&
       (GribLevelType(sec1) != 115) &&
       (GribLevelType(sec1) != 117) &&
       (GribLevelType(sec1) != 125) &&
       (GribLevelType(sec1) != 127) &&
       (GribLevelType(sec1) != 160) &&
       (GribLevelType(sec1) != 210) )
    {
      Put1Byte(GribLevel1(sec1));
      Put1Byte(GribLevel2(sec1));
    }
  else
    {
      Put2Byte(GribLevel1(sec1));     /* 10 Level                    */    
    }
  Put1Byte(GribYear(sec1));           /* 12 Year of Century          */
  Put1Byte(GribMonth(sec1));          /* 13 Month                    */
  Put1Byte(GribDay(sec1));            /* 14 Day                      */
  Put1Byte(GribHour(sec1));           /* 15 Hour                     */
  Put1Byte(GribMinute(sec1));         /* 16 Minute                   */
  if ( GribNumAvg(sec1) > 0 )
    {
      Put1Byte(2);                    /* 17 Time unit                */
      Put1Byte(0);                    /* 18 Time 1                   */
      Put1Byte(0);                    /* 19 Time 2                   */
      Put1Byte(3);                    /* 20 Timerange flag           */
      Put2Byte(GribNumAvg(sec1));     /* 21 Average                  */
    }
  else
    {
      Put1Byte(GribTimeUnit(sec1));   /* 17 Time unit                */
      if ( GribTimeRange(sec1) == 10 )
	{
	  Put1Byte(GribTimePeriod1(sec1));
	  Put1Byte(GribTimePeriod2(sec1));
	}
      else if ( GribTimeRange(sec1) == 113 || GribTimeRange(sec1) ==   0 )
	{
	  Put1Byte(GribTimePeriod1(sec1));
          Put1Byte(0);
	}
      else if ( GribTimeRange(sec1) ==   4 || GribTimeRange(sec1) ==   2 )
	{
          Put1Byte(0);
	  Put1Byte(GribTimePeriod2(sec1));
	}
      else
	{
          Put1Byte(0);
          Put1Byte(0); 
	}
      Put1Byte(GribTimeRange(sec1));       /* 20 Timerange flag           */
      Put2Byte(GribNumAvg(sec1));          /* 21 Average                  */
   }
   Put1Byte(0);                            /* 23 Missing from averages    */
   Put1Byte(GribCentury(sec1));            /* 24 Century                  */
   Put1Byte(GribSubCenterID(sec1));        /* 25 Subcenter                */
   Put2Byte(GribDecScaleFactor(sec1));     /* 26 Decimal scale factor     */
}

/* GRIB BLOCK 2 - GRID DESCRIPTION BLOCK */

void PutGribBlock2(int *sec1, int *sec2, double *psec2)
{
  int Exponent, Mantissa;
  int one = 1;
  int i;
  int BlockLength = 32;

  if ( GribLevelType(sec1) == LEV_HYBRID ) BlockLength += GribNumVCP(sec2) * 4;

  Put3Byte(BlockLength);   /*  0- 2 Length of Block 2 Byte 0 */

  if ( GribLevelType(sec1) == LEV_HYBRID )
    {
      Put1Byte(GribNumVCP(sec2)); /*  3    NV */
      Put1Byte(          33);   /*  4    PV */
    }
  else
    {
      Put1Byte(           0);   /*  3    0  */
      Put1Byte(         255);   /*  4   255 */
    }

  Put1Byte(GribGridType(sec2));   /*  5    Gauss=4 ,Spectral = 50   */

  if ( GribGridType(sec2) == REP_SPECTRAL )
    {
      Put2Byte(GribPentaJ(sec2));   /*  6- 7 Pentagonal resolution J  */
      Put2Byte(GribPentaK(sec2));   /*  8- 9 Pentagonal resolution K  */
      Put2Byte(GribPentaM(sec2));   /* 10-11 Pentagonal resolution M  */
      Put1Byte(GribRepType(sec2));  /* 12    Representation type      */
      Put1Byte(GribRepMode(sec2));  /* 13    Representation mode      */
      PutnZero(18);                 /* 14-31 reserved                 */
    }  
  else 
    {
      Put2Byte(GribNumLon(sec2));   /*  6- 7 Longitudes               */
      Put2Byte(GribNumLat(sec2));   /*  8- 9 Latitudes                */
      Put3Byte(GribFirstLat(sec2)); /* 10-12 Latitude  of Origin      */
      Put3Byte(GribFirstLon(sec2)); /* 13-15 Longitude of Origin      */
      Put1Byte(           128);     /* 16    Resolution flag          */
      if (GribLastLat(sec2) < 0) GribLastLat(sec2) = 8388608 - GribLastLat(sec2);
      Put3Byte(GribLastLat(sec2));  /* 17-19 Latitude  of Extreme     */
      Put3Byte(GribLastLon(sec2));  /* 20-22 Longitude of Extreme     */
      Put2Byte(GribLonIncr(sec2));  /* 23-24 i - direction increment  */
      Put2Byte(GribNumPar(sec2));   /* 25-26 Latitudes Pole->Equator  */
      Put1Byte(             0);     /* 27    Scanning mode            */
      PutnZero(             4);     /* 28-31 reserved                 */
    }

  if ( GribLevelType(sec1) == LEV_HYBRID )
     for (i = 0; i < GribNumVCP(sec2); ++i)
       {
	 confp3(psec2[10+i], &Exponent, &Mantissa, BitsPerInt, one);
	 Put1Byte(Exponent);
	 Put3Byte(Mantissa);
       }
}

/* GRIB BLOCK 4 - BINARY DATA BLOCK */

void PutGribBlock4(int *sec2, int *sec4, double Field[])
{
  int  i;
  int  BlockLength,PackStart,Dim,Flag;
  int  Scale,pval;
  double ZScale,Factor,Fmin,Fmax;

  Dim = GribNumValues(sec4);

  if ( GribGridType(sec2) == REP_SPECTRAL )
    {
      PackStart   = 1;
      Flag        = 128 + 8;
      BlockLength = 14 + 2 * Dim;
    }
  else
    {
      PackStart   = 0;
      Flag        = 8;
      BlockLength = 12 + 2 * Dim;
    }

  Fmin = Fmax = Field[PackStart];
  for (i = PackStart+1; i < Dim; ++i)
    {
      if (Fmin > Field[i]) Fmin = Field[i];
      if (Fmax < Field[i]) Fmax = Field[i];
    }

   ZScale = (Fmax - Fmin) / (pow(2.0,(double)(16 + 1)) - 1);
   if (ZScale >= 10e-100) ZScale = log(ZScale) / log(2.0) + 2.0;
   else                   ZScale = 0.0;

   Scale  = (int) floor(ZScale);
   Factor = 1;

   if (Fmax > Fmin) {
      if (Scale < 0) Factor = pow(2.0,(double)(-Scale));
      else           Factor = 1.0 / pow(2.0,(double)(Scale));
      if (Scale < 0) Scale = 32768 - Scale;
   }

   ref2ibm(&Fmin, BitsPerInt);

   Put3Byte(BlockLength);    /*  0-2 Length of Block 4        */
   Put1Byte(Flag       );    /*  3   Flag & Unused bits       */
   Put2Byte(Scale      );    /*  4-5 Scale factor             */
   Put1Real(Fmin       );    /*  6-9 Reference value          */
   Put1Byte(         16);    /* 10   Packing size             */

   if (PackStart) Put1Real(Field[0]);

   for (i = PackStart; i < Dim; ++i) {
      pval = (int) ((Field[i] - Fmin) * Factor + 0.5);
      IGrib[z  ] = pval >>  8;
      IGrib[z+1] = pval;
      z += 2;
   }
   Put1Byte(0);              /*  Fillbyte                     */
}




void grbEncodeDP(int *sec0, int *sec1, int *sec2, double *psec2, int *sec3,
		 double *psec3, int *sec4, double *psec4, int klenp, int *kgrib,
		 int kleng, int *kword, char *hoper, int *kret)
{
  int gribLen;
  size_t p, len;
  char *CGrib;

  CGrib = (char *) kgrib;

  len = (klenp * 2 + 2000); /* only for 16 bit packing */
  IGrib = (int *) malloc(len*sizeof(int));

  PutGribBlock0();
  PutGribBlock1(sec1);
  PutGribBlock2(sec1, sec2, psec2);
  PutGribBlock4(sec2, sec4, psec4);
  PutGribBlock5();

  gribLen = z;

  if ( gribLen > len )
    {
      fprintf(stderr, "codegb: IGrib buffer to small!!!\n");
      fprintf(stderr, "codegb: igriblen = %d  gribLen = %d\n", len, gribLen);
      exit(99);
    }

#ifdef _CRAY
  p = _pack((long *)IGrib,(char *)CGrib,gribLen,-1);
#else
  {
    int i;
    for (i = 0; i < gribLen; ++i) CGrib[i] = (char) IGrib[i];

    p = gribLen;
  }
#endif

  free(IGrib); IGrib = NULL;

  *kword = gribLen / sizeof(int);
  if ( *kword * sizeof(int) != gribLen ) *kword += 1;

  *kret = 0;
}


int Get2Byte(UCHAR ptr[])
{
   return ((ptr[0] << 8) + ptr[1]);
}

int Get3Byte(UCHAR ptr[])
{
   return ((ptr[0] << 16) + (ptr[1] << 8) + ptr[2]);
}

void grbDecodeDP(int *sec0, int *sec1, int *sec2, double *psec2, int *sec3,
		 double *psec3, int *sec4, double *psec4, int klenp, int *kgrib,
		 int kleng, int *kword, char *hoper, int *kret)
{
  static char func[] = "grbDecodeDP";
  UCHAR *igrib;
  unsigned int *lgrib;
  int  ReducedGrid = FALSE, VertCoorTab = FALSE;
  int  nolon, nolat;
  int *rowpt = 0;
  int  cplx, ipower, jup, kup, mup;
  int  locnv = 0, locnl, locnd;
  int  jlenl, jlend, jlenc;
  int  i;
  int  iflag,irep,lnil,jscale,imiss;
  int  Flag, Lats, Lons;
  int isLen = 0, pdsLen = 0, gdsLen = 0, bmsLen = 0, bdsLen = 0, esLen = 0;
  int gribLen = 0, grib1Len = 0;
  double rmiss;
  double fmin = 0., zscale = 0.;
  double *ztemp;
  int iexp, imant;
  UCHAR *grib;
  int    grib1offset   ;
  double *fpdata;
  extern int GRB_Debug;

  fpdata = psec4;

  grib = (UCHAR *) &kgrib[0];

  GribEdition(sec0) = grib[7];
  if ( GribEdition(sec0) == 1 ) grib1Len = Get3Byte(grib+4);
  grib1offset = GribEdition(sec0) * 4;

  igrib = grib + 4 + grib1offset;

  /***************************************/
  /* block 1 - product definition block. */
  /***************************************/

  pdsLen = Get3Byte(igrib);

  GribCodeTable(sec1)      = igrib[ 3];
  GribCenterID(sec1)       = igrib[ 4];
  GribModelID(sec1)        = igrib[ 5];
  GribGridDefinition(sec1) = igrib[ 6];
  GribBlock2Included(sec1) = igrib[ 7];
  GribParameter(sec1)      = igrib[ 8];
  GribLevelType(sec1)      = igrib[ 9];

  if ( (GribLevelType(sec1) !=  20) && 
       (GribLevelType(sec1) !=  99) && 
       (GribLevelType(sec1) != 100) && 
       (GribLevelType(sec1) != 103) && 
       (GribLevelType(sec1) != 105) && 
       (GribLevelType(sec1) != 107) && 
       (GribLevelType(sec1) != 109) && 
       (GribLevelType(sec1) != 111) && 
       (GribLevelType(sec1) != 113) && 
       (GribLevelType(sec1) != 115) && 
       (GribLevelType(sec1) != 117) && 
       (GribLevelType(sec1) != 125) && 
       (GribLevelType(sec1) != 127) && 
       (GribLevelType(sec1) != 160) && 
       (GribLevelType(sec1) != 210) )
    {
      GribLevel1(sec1)         = igrib[10];
      GribLevel2(sec1)         = igrib[11];
    }
  else
    {
      GribLevel1(sec1)         = Get2Byte(igrib+10);
      GribLevel2(sec1)         = 0;
    }

  GribYear(sec1)           = igrib[12];
  GribMonth(sec1)          = igrib[13];
  GribDay(sec1)            = igrib[14];
  GribHour(sec1)           = igrib[15];
  GribMinute(sec1)         = igrib[16];
  GribTimeUnit(sec1)       = igrib[17];
  GribTimePeriod1(sec1)    = igrib[18];
  GribTimePeriod2(sec1)    = igrib[19];
  GribTimeRange(sec1)      = igrib[20];
  GribNumAvg(sec1)         = Get2Byte(igrib+21);
  GribNumMiss(sec1)        = igrib[23];

  if ( GribEdition(sec0) == 1 )
    {
      GribCentury(sec1)        = igrib[24];
      GribSubCenterID(sec1)    = igrib[25];
      GribDecScaleFactor(sec1) = igrib[26];
    }
  else
    {
      GribCentury(sec1)        = 1;
      GribSubCenterID(sec1)    = 0;
      GribDecScaleFactor(sec1) = 0;
    }

  if ( pdsLen > 28 )
    {
      for ( i = 0; i < (pdsLen-28); i++)
	sec1[23+i] = igrib[27+i];
    }

  igrib += pdsLen;

  /*************************************/
  /* block 2 - grid description block. */
  /*************************************/

  Flag = GribBlock2Included(sec1)       ;

  if ( Flag & 128 )
    {
      gdsLen = Get3Byte(igrib);
      if ( igrib[4] != 255 && igrib[4] != 0 )
	{ /* Either vct or redgrd */
	  if ( igrib[3] != 0 )
	    { /* we have vct */
	      VertCoorTab = TRUE;
	    }
	  else
	    {
	      VertCoorTab = FALSE;
	    }
	  ReducedGrid = (gdsLen - 32 - 4 * igrib[3]);
	}

      if ( GribEdition(sec0) == 0 )
	{
           if ((gdsLen - 32) > 0) VertCoorTab = TRUE;
           else                   VertCoorTab = FALSE;
	}

      if ( ReducedGrid )
	{
	  locnl       = (igrib[4]) - 1 + (VertCoorTab * 4 * igrib[3]);
	  jlenl       = (gdsLen - locnl)  >> 1;
	  nolon       = ((igrib[25]<<8)+igrib[26]) * 4;
	  nolat       = ((igrib[25]<<8)+igrib[26]) * 2;
	  rmiss       = -999.9;
	  rowpt       = (int *) malloc(jlenl*sizeof(int));
	  for (i = 0; i < jlenl; i++)
            rowpt[i] = (igrib[locnl+2*i] << 8) +
                        igrib[locnl+2*i+1];
	  if ( GRB_Debug )
            fprintf(stderr,"+    decogb: jlenl %3d  locnl %3d\n",jlenl,locnl);
	}

      GribGridType(sec2) = igrib[5];

      if ( GribGridType(sec2) == REP_REGULAR || GribGridType(sec2) == REP_GAUSS )
	{
         if ( ReducedGrid )
	   {
	     Lats = Get2Byte(igrib+8);
	     Lons = 2 * Lats;
	   }
	 else
	   {
	     Lons = Get2Byte(igrib+6);
	     Lats = Get2Byte(igrib+8);
	   }
	 GribNumLon(sec2) = Lons;
	 GribNumLat(sec2) = Lats;

	 GribFirstLat(sec2) = Get3Byte(igrib+10);
	 GribFirstLon(sec2) = Get3Byte(igrib+13);
         GribResFlag(sec2)  = igrib[16];
	 GribLastLat(sec2)  = Get3Byte(igrib+17);
	 if ( GribLastLat(sec2) > 8388608 ) GribLastLat(sec2) = 8388608 - GribLastLat(sec2);
	 GribLastLon(sec2)  = Get3Byte(igrib+20);
	 GribLonIncr(sec2)  = Get2Byte(igrib+23);
	 GribNumPar(sec2)   = Get2Byte(igrib+25);

	 /*
         if ( Lons != Longitudes || Lats != Latitudes )
	   Abort("Latitude/Longitude Conflict");
	   */
	}
      else if ( GribGridType(sec2) == REP_SPECTRAL )
	{
	  GribPentaJ(sec2) = Get2Byte(igrib+ 6); /* Truncation */
	  GribPentaK(sec2) = Get2Byte(igrib+ 8); /* Truncation */
	  GribPentaM(sec2) = Get2Byte(igrib+10); /* Truncation */
	}
      else
	Error(func, "Representation Conflict");

      /*    vertical coordinate parameters for hybrid levels.     */
      /*    get number of vertical coordinate parameters, if any. */

      if ( VertCoorTab == TRUE )
	{
          if ( GribEdition(sec0) == 0 )
	    {
	      locnv = 32;
	      GribNumVCP(sec2) = (gdsLen - 32) >> 2;
	    }
	  else
	    {
	      locnv = (igrib[4]) - 1;
	      GribNumVCP(sec2) =  igrib[3];
	    }

	  for (i = 0; i < GribNumVCP(sec2); i++)
	    {
	      iexp   = (igrib[locnv+4*i  ]);
	      imant  =((igrib[locnv+4*i+1]) << 16) +
		      ((igrib[locnv+4*i+2]) <<  8) +
	  	       (igrib[locnv+4*i+3]);
	      psec2[10+i] = decfp2(iexp,imant);
	    }
	}
      igrib += gdsLen;
    }

  /********************************************************************/
  /*    block 3 - bit map block.                                      */
  /********************************************************************/

  if ( Flag & 64 )
    {
      fprintf(stderr," bit map block not yet defined\n");
      bmsLen = Get3Byte(igrib);
      igrib += bmsLen;
    }

  /********************************************************************/
  /*    block 4 - binary data block.                                  */
  /********************************************************************/

  /* get length of binary data block. */

  bdsLen = Get3Byte(igrib);

  /* 4 bit flag / 4 bit count of unused bits at end of block octet. */

  iflag = igrib[3];

  /* 0------- grid point           */
  /* 1------- spherical harmonics  */

  irep  = iflag >> 7;

  /* -0------  simple packing */
  /* -1------ complex packing */

  cplx  = (iflag >> 6) & 1;

  sec4[3] = 0;
  if ( cplx > 0 ) sec4[3] = 64;

  /* ----++++ number of unused bits at end of section) */

  lnil  = iflag & 15;
  
  /* scale factor (2 bytes) */;

  jscale = ( igrib[4] <<  8 ) + igrib[5];
  if ( jscale > 32767 ) jscale = 32768 - jscale;

  /* get reference value (fmin) in grib format (iexp,imant) */

  iexp  = (igrib[6]);
  imant = (igrib[7] << 16) + (igrib[8] << 8) + igrib[9];

  /* check for missing data indicators. */

  imiss = (jscale == 65535 && iexp == 255 && imant == 16777215);

  /* convert reference value and scale factor. */

  if ( imiss == 0 )
    {
      fmin = decfp2(iexp, imant);
      zscale = pow(2.0, (double)jscale);
    }

  /* get number of bits in each data value. */

  GribNumBits(sec4) = igrib[10];

  /* octet number of start of packed data */
  /* calculated from start of block 4 - 1 */

  locnd = 11;

  /* if data is in spherical harmonic form, distinguish   */
  /* between simple/complex packing (cplx = 0/1)          */

  if ( irep == 1 && cplx == 0 )
    {
      /*    no unpacked binary data present */

      jup = kup = mup = 0;

      /*    octet number of start of packed data */
      /*    calculated from start of block 4 - 1 */

      locnd = 15;

      /*    get real (0,0) coefficient in grib format and     */
      /*    convert to floating point.                        */

      iexp  = (igrib[11]);
      imant = ((igrib[12]) << 16) +
              ((igrib[13]) <<  8) +
               (igrib[14]);

      if (imiss) *fpdata++ = 0.0;
      else       *fpdata++ = decfp2(iexp,imant);
    }

  if (irep == 1 && cplx == 1)
    {
      /*    scaling factor */

      ipower = (igrib[13] << 8) + igrib[14];
      if (ipower > 32767) ipower = 32768 - ipower;

      sec4[16] = ipower;

      /*    pentagonal resolution parameters of the */
      /*    unpacked section of data field          */

      jup   = igrib[15];
      kup   = igrib[16];
      mup   = igrib[17];

      sec4[17] = jup;
      sec4[18] = kup;
      sec4[19] = mup;

      /*    unpacked binary data */

      locnd = 18;
      for (i = 0; i < ((jup+1)*(jup+2)); i++)
	{
          iexp   = (igrib[locnd+4*i  ]);
          imant  =((igrib[locnd+4*i+1]) << 16) +
                  ((igrib[locnd+4*i+2]) <<  8) +
                   (igrib[locnd+4*i+3]);

          if (imiss) *fpdata++ = 0.0;
          else       *fpdata++ = decfp2(iexp,imant);
	}
      locnd = 18 + 4 * (jup+1) * (jup+2);
    }

  /* decode data values to floating point and store in fpdata. */
  /* first calculate the number of data values.                */
  /* Take into account that spherical harmonics can be packed  */
  /* simple (cplx = 0) or complex (cplx = 1)                   */

  jlend = bdsLen - locnd;
  jlend = (jlend * 8 - lnil) / GribNumBits(sec4) ;

  sec4[0] = jlend;
  if ( *hoper == 'J' ) return;

  /* check length of output array. */
  
  if ( ! (*hoper == 'J') )
    if (jlend+irep > klenp)
      {
	fprintf(stderr," hoper %c\n", *hoper);
	fprintf(stderr," values to be decoded are  %d\n", jlend+irep);
	fprintf(stderr," array size:               %d\n", klenp);
	exit(1);
      }

  if ( imiss ) memset((char *)fpdata, 0, jlend*sizeof(double));
  else
    {
      igrib += locnd;
      jlenc = jlend * GribNumBits(sec4) / 8;
      lgrib = (unsigned int *) malloc(jlenc*sizeof(int));
#ifdef _CRAY
      _unpack((char *)igrib,(long *)lgrib,jlenc,-1);
#else
      for (i = 0; i < jlenc; i++) lgrib[i] = igrib[i];
#endif
      if (GribNumBits(sec4) ==  8)
         for (i = 0; i < jlend; i++) {
             fpdata[i] = fmin + zscale *
                           lgrib[i];
         }
      else if (GribNumBits(sec4) == 16)
         for (i = 0; i < jlend; i++) {
             fpdata[i] = fmin + zscale *
                         ((lgrib[2*i  ] <<  8) +  lgrib[2*i+1]);
         }
      else if (GribNumBits(sec4) == 24)
         for (i = 0; i < jlend; i++) {
             fpdata[i] = fmin + zscale *
                         ((lgrib[3*i  ] << 16) + (lgrib[3*i+1] <<  8) +
                           lgrib[3*i+2]);
         }
      else if (GribNumBits(sec4) == 32)
         for (i = 0; i < jlend; i++) {
             fpdata[i] = fmin + zscale *
                         ((lgrib[4*i  ] << 24) + (lgrib[4*i+1] << 16) +
                          (lgrib[4*i+2] <<  8) +  lgrib[4*i+3]);
      }  else {
         fprintf(stderr," Unimplemented packing factor %d\n", GribNumBits(sec4));
         exit(1);
      }
      free(lgrib); lgrib = NULL;
   }
  if (ReducedGrid) {
    ztemp = (double *) malloc(nolon*nolat*sizeof(double));
    (void)qu2reg2(fpdata,rowpt,&nolat,&nolon,ztemp,&rmiss,kret);
    free(ztemp); ztemp = NULL;
    free(rowpt); rowpt = NULL;
  }

  if ( GribEdition(sec0) == 1 ) isLen = 8;
  esLen = 4;

  gribLen = isLen + pdsLen + gdsLen + bmsLen + bdsLen + esLen;

  if ( GribEdition(sec0) == 1 )
    if ( gribLen != grib1Len )
      Warning(func, "grib1Len = %d gribLen = %d", grib1Len, gribLen);

  GribLength(sec0) = gribLen;

  *kword = gribLen / sizeof(int);
  if ( *kword * sizeof(int) != gribLen ) *kword += 1;

  *kret = 0;
}

void
gribexdp(int *sec0, int *sec1, int *sec2, double *psec2, int *sec3,
	 double *psec3, int *sec4, double *psec4, int klenp, int *kgrib,
	 int kleng, int *kword, char *hoper, int *kret)
{
  static char func[] = "gribexdp";

  if ( *hoper == 'D' || *hoper == 'J' )
    grbDecodeDP(sec0, sec1, sec2, psec2, sec3,
		psec3, sec4, psec4, klenp, kgrib,
		kleng, kword, hoper, kret);
  else if ( *hoper == 'C' )
    grbEncodeDP(sec0, sec1, sec2, psec2, sec3,
		psec3, sec4, psec4, klenp, kgrib,
		kleng, kword, hoper, kret);
  else
    {
      Error(func, "gribex: oper %c unsupported\n", *hoper);
      *kret=-9;
    }
}

void
gribexsp(int *sec0, int *sec1, int *sec2, float *psec2, int *sec3,
	 float *psec3, int *sec4, float *psec4, int klenp, int *kgrib,
	 int kleng, int *kword, char *hoper, int *kret)
{
  printf("gribexsp not implemented\n");
  exit (-1);
}
