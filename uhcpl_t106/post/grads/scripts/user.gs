function user ()
*
* return the current username
*
*'!finger $USER | head -1 | cut -f 3 -d":" > grads_user'
'!ypcat passwd | grep $USER | cut -f 5 -d":" > grads_user'
user=sublin(read(grads_user),2)
'!rm -f grads_user'
say "user: "user
return user
