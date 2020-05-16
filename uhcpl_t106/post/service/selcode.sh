#!/bin/sh
#
set -e
#
AFTER=../burn/after
#
DPATH=`dirname $1`
DFILE=`basename $1`
#
if [ ! -f $DPATH/$DFILE ] ; then
  echo "File $DPATH/$DFILE not found"
  exit
fi
#
$AFTER $DPATH/${DFILE} ${DFILE}.srv  << EOR
 &SELECT
  TYPE=70,
  CODE=129,130,131,132,134,135,148,149,151,169
  LEVEL=3000,5000,7000,10000,25000,50000,70000,85000,95000,100000,
  MEAN=1,
  GRIB=0,
 &END
EOR
#
exit
