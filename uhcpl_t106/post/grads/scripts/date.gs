function date ()
*
* return the current date
*
'!date "+%d/%m/%Y" > grads_date'
date=sublin(read(grads_date),2)
'!rm -f grads_date'
say "date: "%date
return date
