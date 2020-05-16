COMMON graph_common, $
  thick,      $        ; line thickness
  xthick,     $        ; xthick thickness
  ythick,     $        ; ythick thickness
  charthick,  $        ; character thickness
  charsize,   $        ; char size
  xstyle,     $        ; x axis style
  ystyle,     $        ; y axis style
  title,      $       
  xtitle,      $
  ytitle,      $
  version,    $    ; version string 
  ntimes,     $
  nitems,     $    ; number of data items (that span time)
  data_nstep, $    ; timestep variable (in case data not every timestep)
  data_sort,  $    ; data fields names
  data_name,  $    ; data fields names
  data_val,   $    ; data values
  data_av,    $    ; data time average
  colour,     $    ; 1 means plot in color
  maxlevs,lev,index, number_list_max, $
  number_list, char_list, lev_list, title_list, $
;
; STRING WIDGET CODES - FAIRLY INEFFICIENT WAY OF DOING THIS!
  string1,string2,string3, string4, string5, string6, string7,string8, $   
  string_label, $
  list_thick,list_xthick,list_ythick,list_charthick,list_charsize, $
  error0, $     ; error message window number
  graph_xsize, $ ; width of data plotting window
  graph_ysize, $ ; height of data plotting window
  ngraph_cols, $ ; nol of columns in draw window
  base, row_base1,col_left1,col_right1,one_row,two_row,  $  ;widget base numbers
  igraph, $     ; current graph
  ngraphs,    $ ; current number of graphs
  graph_win, $  ; list of widget window asscoiated with graphs
  graph_wcur,$  ; window number correspondign to igraph
  ch_graph, $   ; change graph index - equal to igraph if local, array if global
  plot_type, $  ; type of plot, time/data or data/z
  time_index, $ ; index of time
  timestep, $   ; timestep of model
  av_min, av_max, $ ;min and max limits for time averaging
  filename, $   ; filename for postscript file
  ak,bk,    $   ; A's and B's for calculating pressure levels
  press_index, $; index of pressure in results 
  t_index,     $; index of temperature in results 
  q_index,     $; index of moisture in results 
  press_plot, $ ; 0 = plot against level, 1 = plot against height
  xmax_auto, $    ; index for automatic minmax setting
  ymax_auto, $    ; index for automatic minmax setting
  max_index, $  ; index for which of (x/y)(min/max) string applies to 
  max_val,   $  ; values for x/y min/max
  over_max, $   ; max number of overplots -1 
  over_index,$  ; current overplot index
  over_active,$ ; which overplot are active in each case
  over_label, $ ; labels for overplots
  bgroup_over_index, bgroup_over_active,$
  label, $      ; 1=find mouse to label graphs
  label_x,label_y, $
  win_xsize,win_ysize, $
  plot_win, $
  frame_thick, $
  max,nlev           ; large number - max poss data items