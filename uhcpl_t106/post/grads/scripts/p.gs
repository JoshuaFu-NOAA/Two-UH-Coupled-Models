* ****************************************************************************
* ******      palette.gs - A script for  selfmade color palettes      ********
* ******                                                              ********
* ******                          written by                          ********
* ******                                                              ********
* ******   Hans-Stefan Bauer, Max-Planck Institut fuer Meteorologie   ********
* ******                  (stefan.bauer@dkrz.de)                      ********
* ******                                                              ********
* ****** Last Edit: 17.3.1997                                         ********
* ****************************************************************************

 function pal(palette)

  mainbut=7
  butlength=1.2
  buthight=0.4
  x=7.9
  y=8.3
  'draw button 21 7.9 7.9 'butlength' 'buthight' rainbow'
  'draw button 22 7.9 7.5 'butlength' 'buthight' PalGrey'
  'draw button 23 7.9 7.1 'butlength' 'buthight' pallai'
  'draw button 24 7.9 6.7 'butlength' 'buthight' C16a'
  'draw button 25 7.9 6.3 'butlength' 'buthight' C16b'
  'draw button 26 7.9 5.9 'butlength' 'buthight' C16c'
  'draw button 27 7.9 5.5 'butlength' 'buthight' C32a'
  'draw button 28 7.9 5.1 'butlength' 'buthight' C32b'
  'draw button 29 7.9 4.7 'butlength' 'buthight' blue_red'
  'draw button 30 7.9 4.3 'butlength' 'buthight' No_10'
  'draw button 31 7.9 3.9 'butlength' 'buthight' farbglet'
  'q pos'
   if (_trace='on'); say result; endif;
   btn=subwrd(result,7)
   palette=_palette
   if (btn=-999); rc = delmen(x,y,i,mainbut,palette); endif
   if (btn>=1 & btn<=9)
     rc = delmen(x,y,i,mainbut,palette)
     rc = choice(btn)
  endif
  _palette = palette
  if (btn>=94 & btn<=99); rc = mainmenu(); endif
  if (btn=21)
    if (_trace='on'); say 'Palette: rainbow'; endif
    _palette='rainbow'
    'set rbcols auto'
  endif
*
* ************************************ Palette PalGrey ************************
*
 if (btn=22)
   'set rgb    16    8    8    8'
   'set rgb    17   16   16   16'
   'set rgb    18   24   24   24'
   'set rgb    19   32   32   32'
   'set rgb    20   40   40   40'
   'set rgb    21   48   48   48'
   'set rgb    22   56   56   56'
   'set rgb    23   64   64   64'
   'set rgb    24   72   72   72'
   'set rgb    25   80   80   80'
   'set rgb    26   88   88   88'
   'set rgb    27   96   96   96'
   'set rgb    28  104  104  104'
   'set rgb    29  112  112  112'
   'set rgb    30  120  120  120'
   'set rgb    31  128  128  128'
   'set rgb    32  136  136  136'
   'set rgb    33  144  144  144'
   'set rgb    34  152  152  152'
   'set rgb    35  160  160  160'
   'set rgb    36  168  168  168'
   'set rgb    37  176  176  176'
   'set rgb    38  184  184  184'
   'set rgb    39  192  192  192'
   'set rgb    40  200  200  200'
   'set rgb    41  208  208  208'
   'set rgb    42  216  216  216'
   'set rgb    43  224  224  224'
   'set rgb    44  232  232  232'
   'set rgb    45  240  240  240'
   'set rgb    46  248  248  248'
   'set rgb    47  255  255  255'
   'set rbcols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47'
   palette='PalGrey'

endif
*
* ************************************ Palette pallai ********************************
* 
 if (btn=23)
   'set rgb  16  0  75  45'
   'set rgb  17  0  80  60'
   'set rgb  18  0  89  68'
   'set rgb  19  0  105  85'
   'set rgb  20  0  119  91'
   'set rgb  21  0  138  105'
   'set rgb  22  0  160  113'
   'set rgb  23  0  189  110'
   'set rgb  24  0  190  130'
   'set rgb  25  0  200  140'
   'set rgb  26  0  210  150'
   'set rgb  27  0  220  160'
   'set rgb  28  0  240  170'
   'set rgb  29  0  255  180'
   'set rgb  30  0  195  140'
   'set rgb  31  0  195  80'
   'set rgb  32  0  216  46'
   'set rgb  33  100  220  0'
   'set rgb  34  115  210  0'
   'set rgb  35  130  200  30'
   'set rgb  36  145  210  25'
   'set rgb  37  160  220  20'
   'set rgb  38  185  230  10'
   'set rgb  39  200  240  0'
   'set rgb  40  220  245  0'
   'set rgb  41  230  250  0'
   'set rgb  42  240  255  0'
   'set rgb  43  253  245  0'
   'set rgb  44  255  251  175'
   'set rgb  45  255  255  200'
   'set rbcols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45'
   palette='pallai'
 endif
*
* ********************************* Palette C16a *****************************
*
 if (btn=24)
   'set rgb  16  122  59  255'
   'set rgb  17  193  0  227'
   'set rgb  18  227  0  189'
   'set rgb  19  255  70  108'
   'set rgb  20  255  98  89'
   'set rgb  21  255  173  122'
   'set rgb  22  255  194  97'
   'set rgb  23  255  250  110'
   'set rgb  24  168  227  0'
   'set rgb  25  0  216  115'
   'set rgb  26  0  192  160'
   'set rgb  27  0  151  250'
   'set rgb  28  59  118  255'
   'set rgb  29  122  139  255'
   'set rbcols 16 17 18 19 20 21 23 23 24 25 26 27 28 29'
   palette='C16a'
  endif
*
* ********************************* Palette C16b *****************************
*
 if (btn=25)
   'set rgb  16  141  0  147'
   'set rgb  17  255  0  255'
   'set rgb  18  255  0  0'
   'set rgb  19  255  56  43'
   'set rgb  20  255  142  89'
   'set rgb  21  255  255  0'
   'set rgb  22  0  189  71'
   'set rgb  23  0  255  0'
   'set rgb  24  0  138  97'
   'set rgb  25  0  160  147'
   'set rgb  26  0  255  255'
   'set rgb  27  0  61  175'
   'set rgb  28  0  0  255'
   'set rgb  29  22  0  60'
   'set rbcols 16 17 18 19 20 21 23 24 25 26 27 28 29'
   palette='C16b'
  endif
*
* ********************************* Palette C16c *****************************
*
 if (btn=26)
   'set rgb  16  35  0  64'
   'set rgb  17  186  60  0'
   'set rgb  18  238  44  0'
   'set rgb  19  255  168  89'
   'set rgb  20  255  209  116'
   'set rgb  21  255  255  0'
   'set rgb  22  160  219  0'
   'set rgb  23  69  156  0'
   'set rgb  24  0  255  0'
   'set rgb  25  0  160  147'
   'set rgb  26  0  255  255'
   'set rgb  27  0  120  198'
   'set rgb  28  0  89  147'
   'set rgb  29  0  56  133'
   'set rbcols 16 17 18 19 20 21 22 23 24 25 26 27 28 29'
   palette='C16c'
  endif
*
* *************************************** Palette C32a ***********************************
*
  if (btn=27)
    'set rgb  16  92  70  255'
    'set rgb  17  0  0  255'
    'set rgb  18  0  59  252'
    'set rgb  19  0  83  230'
    'set rgb  20  0  97  219'
    'set rgb  21  0  112  225' 
    'set rgb  22  0  151  250'
    'set rgb  23  0  186  186'
    'set rgb  24  0  160  123'
    'set rgb  25  0  138  97'
    'set rgb  26  0  151  102'
    'set rgb  27  0  189  110'
    'set rgb  28  0  202  65'
    'set rgb  29  120  219  0'
    'set rgb  30  168  227  0'
    'set rgb  31  255  249  70'
    'set rgb  32  255  191  85'
    'set rgb  33  255  168  89'
    'set rgb  34  255  142  89'
    'set rgb  35  255  100  65'
    'set rgb  36  255  56  43'
    'set rgb  37  225  0  35'
    'set rgb  38  247  0  83'
    'set rgb  39  233  0  134'
    'set rgb  40  199  0  165'
    'set rgb  41  175  0  175'
    'set rgb  42  153  0  195'
    'set rgb  43  118  0  216'
    'set rgb  44  86  0  230'
    'set rgb  45  25  0  68'
    'set rbcols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45'
    palette='C32a'
  endif
*
* *********************************** Palette C32b **************************************
*
  if (btn=28)
    'set rgb  16  198  192  255'
    'set rgb  17  128  137  255'
    'set rgb  18  122  139  255'
    'set rgb  19  122  149  255'
    'set rgb  20  116  154  255'
    'set rgb  21  110  159  255'
    'set rgb  22  104  173  255'
    'set rgb  23  0  186  186'
    'set rgb  24  0  192  160'
    'set rgb  25  0  243  169'
    'set rgb  26  0  216  115'
    'set rgb  27  0  230  74'
    'set rgb  28  105  238  0'
    'set rgb  29  161  240  0'
    'set rgb  30  197  243  0'
    'set rgb  31  255  250  110'
    'set rgb  32  255  209  116'
    'set rgb  33  255  190  119'
    'set rgb  34  255  173  122'
    'set rgb  35  255  162  122'
    'set rgb  36  255  98  89'
    'set rgb  37  255  122  128'
    'set rgb  38  255  122  147'
    'set rgb  39  255  110  170'
    'set rgb  40  255  80  204'
    'set rgb  41  255  70  225'
    'set rgb  42  245  43  255'
    'set rgb  43  216  89  255'
    'set rgb  44  187  110  255'
    'set rgb  45  52  45  62'
    'set rbcols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45'
    palette='C32b'
  endif
*
* **************************************** Palette blue_red *****************************
*
  if (btn =29)
    'set rgb  16  0   37   89'
    'set rgb  17  0   56   133'
    'set rgb  18  0   77   182'
    'set rgb  19  0   97   218'
    'set rgb  20  59   118   255'
    'set rgb  21  80   130   255'
    'set rgb  22  116   154   255'
    'set rgb  23  130   180   255'
    'set rgb  24  150   200   255'
    'set rgb  25  170   220   255'
    'set rgb  26  190   240   255'
    'set rgb  27  210   255   255'
    'set rgb  28  230   255   255'
    'set rgb  29  255   255   255'
    'set rgb  30  255   240   200'
    'set rgb  31  255   220   180'
    'set rgb  32  255   200   160'
    'set rgb  33  255   180   140'
    'set rgb  34  255   160   120'
    'set rgb  35  255   140   100'
    'set rgb  36  255   120   80'
    'set rgb  37  255   100   60'
    'set rgb  38  255   80   40'
    'set rgb  39  255   60   20'
    'set rgb  40  255   40   0'
    'set rgb  41  255   20   0'
    'set rgb  42  255   0   0'
    'set rgb  43  225   0   35'
    'set rgb  44  168   0   26'
    'set rgb  45  130   0   20'
    'set rgb  46  110   0   10'
    'set rgb  47  90    0   0'
    'set rbcols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47'
    palette='blue_red'
  endif
*
* ********************** Palette No_10 ******************
*
 if (btn=30)
   'set rgb 21  150 250 255'
   'set rgb 22 250 50  200'
   'set rgb 23 255 175  255'
   'set rbcols 23 9 14 4 11 5 21 7 12 8 2 6 22 '
   palette='No_10'
 endif


*
* ************************* Palette farbglet **************************
*
  if (btn =31)
    'set rgb  16  255  255  255'
    'set rgb  17  0  104  68'
    'set rgb  18  0  142  53'
    'set rgb  19  195  136  0'
    'set rgb  20  0  39  168'
    'set rgb  21  0  138  105'
    'set rgb  22  0  121  106'
    'set rgb  23  0  186  186'
    'set rgb  24  0  89  147'
    'set rgb  25  247  80  0'
    'set rgb  26  0  120  198'
    'set rgb  27  225  0  35'
    'set rgb  28  255  168  89'
    'set rgb  29  0  202  65'
    'set rgb  30  139  122  255'
    'set rgb  31  255  249  70'
    'set rgb  32  0  243  169'
    'set rgb  33  177  213  255'
    'set rgb  34  255  251  175'
    'set rgb  35  198  192  255'
    'set rbcols 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35'
    palette='farbglet'
  endif

return palette






