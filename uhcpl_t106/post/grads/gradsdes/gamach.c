#include <stdio.h>
#include <math.h>


static union {unsigned long l; unsigned char c[sizeof(long)];} u_byteorder = {1};

#undef  IsBigendian
#define IsBigendian()  ( u_byteorder.c[sizeof(long) - 1] )

/* Machine dependent routines.  These routines depend on machine word
   length any byte ordering. */

#ifdef CRAY

/*mf ---------------------------------------

  64-bit machine (cray) routine of gagby 

---------------------------------------- mf*/

/*  routine to get an integer length one to eight and return it as
    a long int.  VERSION FOR CRAY */

int
gagby (char *ch, int ioff, int ilen)
{

  int ival;
  char *ch1;

  ch1 = (char *) (&ival);
  ival = 0;

  if ( IsBigendian() )
    {
      if (ilen == 1)
	*(ch1 + 7) = *(ch + ioff);
      else if (ilen == 2)
	{
	  *(ch1 + 6) = *(ch + ioff);
	  *(ch1 + 7) = *(ch + ioff + 1);
	}
      else if (ilen == 3)
	{
	  *(ch1 + 5) = *(ch + ioff);
	  *(ch1 + 6) = *(ch + ioff + 1);
	  *(ch1 + 7) = *(ch + ioff + 2);
	}
      else
	{
	  *(ch1 + 4) = *(ch + ioff);
	  *(ch1 + 5) = *(ch + ioff + 1);
	  *(ch1 + 6) = *(ch + ioff + 2);
	  *(ch1 + 7) = *(ch + ioff + 3);
	}
    }
  else
    {
      if (ilen == 1)
	*ch1 = *(ch + ioff);
      else if (ilen == 2)
	{
	  *(ch1 + 7) = *(ch + ioff);
	  *ch1 = *(ch + ioff + 7);
	}
      else if (ilen == 3)
	{
	  *(ch1 + 7) = *(ch + ioff + 5);
	  *(ch1 + 6) = *(ch + ioff + 6);
	  *(ch1 + 5) = *(ch + ioff + 4);
	}
      else
	{
	  *(ch1 + 7) = *(ch + ioff + 4);
	  *(ch1 + 6) = *(ch + ioff + 5);
	  *(ch1 + 5) = *(ch + ioff + 6);
	  *(ch1 + 4) = *(ch + ioff + 3);
	}
    }
  return (ival);
}

#else

/*mf ---------------------------------------

  32-bit machine (cray) routine of gagby 

---------------------------------------- mf*/

int
gagby (char *ch, int ioff, int ilen)
{

  int ival;
  char *ch1;

  ch1 = (char *) (&ival);
  ival = 0;

  if ( IsBigendian() )
    {
      if (ilen == 1)
	*(ch1 + 3) = *(ch + ioff);
      else if (ilen == 2)
	{
	  *(ch1 + 2) = *(ch + ioff);
	  *(ch1 + 3) = *(ch + ioff + 1);
	}
      else if (ilen == 3)
	{
	  *(ch1 + 1) = *(ch + ioff);
	  *(ch1 + 2) = *(ch + ioff + 1);
	  *(ch1 + 3) = *(ch + ioff + 2);
	}
      else
	{
	  *ch1 = *(ch + ioff);
	  *(ch1 + 1) = *(ch + ioff + 1);
	  *(ch1 + 2) = *(ch + ioff + 2);
	  *(ch1 + 3) = *(ch + ioff + 3);
	}
    }
  else
    {
      if (ilen == 1)
	*ch1 = *(ch + ioff);
      else if (ilen == 2)
	{
	  *(ch1 + 1) = *(ch + ioff);
	  *ch1 = *(ch + ioff + 1);
	}
      else if (ilen == 3)
	{
	  *(ch1 + 2) = *(ch + ioff);
	  *(ch1 + 1) = *(ch + ioff + 1);
	  *ch1 = *(ch + ioff + 2);
	}
      else
	{
	  *(ch1 + 3) = *(ch + ioff);
	  *(ch1 + 2) = *(ch + ioff + 1);
	  *(ch1 + 1) = *(ch + ioff + 2);
	  *ch1 = *(ch + ioff + 3);
	}
    }
  return (ival);
}

#endif

/* Get bit string value from the character string starting
   at bit offset ioff and number of bits ilen.  ilen should not
   be greater than 24 bits unless byte aligned. */

static char masks[8] = { 0, 127, 63, 31, 15, 7, 3, 1 };

#ifdef CRAY

/*mf ---------------------------------------

  64-bit machine (cray) routine of gagbb 

---------------------------------------- mf*/

int
gagbb (char *ch, int ioff, int ilen)
{
  int ival;
  char *ch1, cc;
  int istrt, iend, cstrt;
  int i, ispac, ioff2, ileav, numb;

  ch1 = (char *) (&ival);
  ival = 0;

  istrt = ioff / 8;
  ispac = ioff - istrt * 8;

  /* Fast path for byte alignment */

  if (ispac == 0 && (ilen == 8 || ilen == 16 || ilen == 24 || ilen == 32))
    {
      if ( IsBigendian() )
	{
	  if (ilen == 8)
	    *(ch1 + 7) = *(ch + istrt);
	  else if (ilen == 16)
	    {
	      *(ch1 + 6) = *(ch + istrt);
	      *(ch1 + 7) = *(ch + istrt + 1);
	    }
	  else if (ilen == 24)
	    {
	      *(ch1 + 6) = *(ch + istrt);
	      *(ch1 + 7) = *(ch + istrt + 1);
	      *(ch1 + 8) = *(ch + istrt + 2);
	    }
	  else
	    {
	      *(ch1 + 4) = *(ch + istrt);
	      *(ch1 + 5) = *(ch + istrt + 1);
	      *(ch1 + 6) = *(ch + istrt + 2);
	      *(ch1 + 7) = *(ch + istrt + 3);
	    }
	}
      else
	{
	  if (ilen == 8)
	    *ch1 = *(ch + istrt);
	  else if (ilen == 16)
	    {
	      *(ch1 + 1) = *(ch + istrt);
	      *ch1 = *(ch + istrt + 1);
	    }
	  else if (ilen == 24)
	    {
	      *(ch1 + 2) = *(ch + istrt);
	      *(ch1 + 1) = *(ch + istrt + 1);
	      *ch1 = *(ch + istrt + 2);
	    }
	  else
	    {
	      *(ch1 + 3) = *(ch + istrt);
	      *(ch1 + 2) = *(ch + istrt + 1);
	      *(ch1 + 1) = *(ch + istrt + 2);
	      *ch1 = *(ch + istrt + 3);
	    }
	}
      return (ival);
    }

  /* Do it the hard way */

  ioff2 = ioff + ilen - 1;
  iend = ioff2 / 8;
  ileav = (iend + 1) * 8 - ioff2;
  numb = iend - istrt;
  if ( IsBigendian () )
    {
      cstrt = 7 - numb;
      if (ispac > 0)
	*(ch1 + cstrt) = *(ch + istrt) & masks[ispac];
      else
	*(ch1 + cstrt) = *(ch + istrt);
      for (i = 1; i <= numb; i++)
	{
	  *(ch1 + cstrt + i) = *(ch + istrt + i);
	}
    }
  else
    {
      if (ispac > 0)
	*(ch1 + numb) = *(ch + istrt) & masks[ispac];
      else
	*(ch1 + numb) = *(ch + istrt);
      for (i = 0; i < numb; i++)
	{
	  *(ch1 + i) = *(ch + iend - i);
	}
    }
  ival = ival >> (ileav - 1);
  return (ival);
}

#else

/*mf ---------------------------------------

  32-bit machine routine of gagbb 

---------------------------------------- mf*/

int
gagbb (char *ch, int ioff, int ilen)
{
  int ival;
  char *ch1, cc;
  int istrt, iend, cstrt;
  int i, ispac, ioff2, ileav, numb;

  ch1 = (char *) (&ival);
  ival = 0;

  istrt = ioff / 8;
  ispac = ioff - istrt * 8;

  /* Fast path for byte alignment */

  if (ispac == 0 && (ilen == 8 || ilen == 16 || ilen == 24 || ilen == 32))
    {
      if ( IsBigendian() )
	{
	  if (ilen == 8)
	    *(ch1 + 3) = *(ch + istrt);
	  else if (ilen == 16)
	    {
	      *(ch1 + 2) = *(ch + istrt);
	      *(ch1 + 3) = *(ch + istrt + 1);
	    }
	  else if (ilen == 24)
	    {
	      *(ch1 + 1) = *(ch + istrt);
	      *(ch1 + 2) = *(ch + istrt + 1);
	      *(ch1 + 3) = *(ch + istrt + 2);
	    }
	  else
	    {
	      *ch1 = *(ch + istrt);
	      *(ch1 + 1) = *(ch + istrt + 1);
	      *(ch1 + 2) = *(ch + istrt + 2);
	      *(ch1 + 3) = *(ch + istrt + 3);
	    }
	}
      else
	{
	  if (ilen == 8)
	    *ch1 = *(ch + istrt);
	  else if (ilen == 16)
	    {
	      *(ch1 + 1) = *(ch + istrt);
	      *ch1 = *(ch + istrt + 1);
	    }
	  else if (ilen == 24)
	    {
	      *(ch1 + 2) = *(ch + istrt);
	      *(ch1 + 1) = *(ch + istrt + 1);
	      *ch1 = *(ch + istrt + 2);
	    }
	  else
	    {
	      *(ch1 + 3) = *(ch + istrt);
	      *(ch1 + 2) = *(ch + istrt + 1);
	      *(ch1 + 1) = *(ch + istrt + 2);
	      *ch1 = *(ch + istrt + 3);
	    }
	}
      return (ival);
    }

  /* Do it the hard way */

  ioff2 = ioff + ilen - 1;
  iend = ioff2 / 8;
  ileav = (iend + 1) * 8 - ioff2;
  numb = iend - istrt;
  if ( IsBigendian () )
    {
      cstrt = 3 - numb;
      if (ispac > 0)
	*(ch1 + cstrt) = *(ch + istrt) & masks[ispac];
      else
	*(ch1 + cstrt) = *(ch + istrt);
      for (i = 1; i <= numb; i++)
	{
	  *(ch1 + cstrt + i) = *(ch + istrt + i);
	}
    }
  else
    {
      if (ispac > 0)
	*(ch1 + numb) = *(ch + istrt) & masks[ispac];
      else
	*(ch1 + numb) = *(ch + istrt);
      for (i = 0; i < numb; i++)
	{
	  *(ch1 + i) = *(ch + iend - i);
	}
    }
  ival = ival >> (ileav - 1);
  return (ival);
}

#endif
