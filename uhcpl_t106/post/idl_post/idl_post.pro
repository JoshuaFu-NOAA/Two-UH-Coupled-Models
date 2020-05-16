;=======================================================================
FUNCTION TS_TO_DAY,ts
  @idl_post_comm
  RETURN,FLOAT(ts+data_nstep(0))*timestep/1440.
END
;=======================================================================
FUNCTION DAY_TO_TS,day
  @idl_post_comm
  RETURN,ROUND(day*1440/timestep)-data_nstep(0) 
END
;=======================================================================
FUNCTION STR, n, FMT=fmt
  IF NOT(keyword_set(FMT)) THEN RETURN, STRTRIM(STRING(n),2) $
                           ELSE RETURN, STRTRIM(STRING(n,FORMAT=fmt),2)
END
;=======================================================================
FUNCTION GET_NUMBER,widget_num,error
  @idl_post_comm

  WIDGET_CONTROL, widget_num, GET_VALUE=string
  error=1
  ON_IOERROR, bad_num0
  num = FLOAT(string)
  error=0
  bad_num0: IF error THEN BEGIN
              WIDGET_CONTROL, error0, SET_VALUE='Error: Not a number!'
              print, STRING(7B) ;ring bell
              num = 0
              WIDGET_CONTROL, widget_num, SET_VALUE=STR(num)
            ENDIF
  RETURN, num
END
;=======================================================================
PRO dump_graph
 @idl_post_comm
  OPENW, 1, 'idl_post.graph'

  PRINTF, 1, title(igraph)
  PRINTF, 1, xtitle(igraph)
  PRINTF, 1, ytitle(igraph)
  PRINTF, 1,' Created by idl_post '+version+' on '+SYSTIME(0)
  num_graphs = 0
  FOR j=0,over_max-1 DO $
    IF (over_active(igraph,j)) THEN BEGIN
      num_graphs = num_graphs+1
      PRINTF, 1, over_label(igraph,j)
    ENDIF

  CASE plot_type OF
  0: BEGIN
     x=FLTARR(N_ELEMENTS(data_nstep),num_graphs)
     y=FLTARR(N_ELEMENTS(data_nstep),num_graphs)
     n_elem=N_ELEMENTS(data_nstep)
     ENDCASE
  1: BEGIN
     x=FLTARR(maxlevs,num_graphs)
     y=FLTARR(maxlevs,num_graphs)
     n_elem=maxlevs
     ENDCASE
  2: BEGIN
     x=FLTARR(maxlevs,num_graphs)
     y=FLTARR(maxlevs,num_graphs)
     n_elem=maxlevs
     ENDCASE
  ENDCASE

  FOR j=0,num_graphs-1 DO BEGIN

    CASE plot_type OF
    0: BEGIN
     x(*,j)=data_nstep*timestep[0]/1440
     y(*,j)=data_val(index(igraph,j),lev(igraph,j),*)
     ENDCASE
    1: BEGIN
     x(*,j)=data_val(index(igraph,j),*,time_index)
     CASE press_plot OF
       0: y(*,j)=findgen(maxlevs)+1
       1: y(*,j)=(ak+bk*EXP(data_val(press_index,0,time_index)))/100.
       ENDCASE
     ENDCASE
    2: BEGIN
     x(*,j)=data_av(index(igraph,j),*)
     CASE press_plot OF
       0: y(*,j)=findgen(maxlevs)+1
       1: y(*,j)=(ak+bk*EXP(data_av(press_index,0)))/100.
       ENDCASE
     ENDCASE
    ENDCASE
  ENDFOR

  vector=fltarr(num_graphs+1)
  FOR i=0,n_elem-1 DO BEGIN
    vector[0]=x(i,0)
    vector[1:num_graphs]=y(i,*)
    PRINTF, 1, vector
  ENDFOR

  CLOSE, 1
  WIDGET_CONTROL, error0, SET_VALUE='graph dumped to file idl_post.graph'

END

;========================================================================
PRO auto_title,t,xt,yt
@idl_post_comm

; make sure titles are correct
  CASE plot_type OF
     0: BEGIN
        IF t THEN title(ch_graph)  = 'LEV '+STR(lev(ch_graph,0)+1)
        IF xt THEN xtitle(ch_graph) = 'Time (days) '
        IF yt THEN ytitle(ch_graph) = STRTRIM(data_name(index(ch_graph,0)),2)
        ENDCASE
     1: BEGIN
        IF t THEN title(ch_graph) = 'TIME '+STR(TS_TO_DAY(time_index))
        IF xt THEN xtitle(ch_graph) = STRTRIM(data_name(index(ch_graph,0)),2)
        CASE press_plot OF 
          0: IF yt THEN ytitle(ch_graph) = 'LEV'
          1: IF yt THEN ytitle(ch_graph) = 'Pressure (hPa)'
          ENDCASE
        ENDCASE
     2: BEGIN
        IF t THEN title(ch_graph) = 'Average Day ' + $
           STRTRIM(STRING(TS_TO_DAY(av_min),FORMAT='(F8.2)'))+' to '+ $
           STRTRIM(STRING(TS_TO_DAY(av_max),FORMAT='(F8.2)'))
        IF xt THEN xtitle(ch_graph) = STRTRIM(data_name(index(ch_graph,0)),2)
        CASE press_plot OF 
          0: IF yt THEN ytitle(ch_graph) = 'LEV'
          1: IF yt THEN ytitle(ch_graph) = 'Pressure (hPa)'
          ENDCASE
        ENDCASE
     ENDCASE
  IF yt THEN over_label(ch_graph,over_index)= $
     STRTRIM(data_name(index(ch_graph,over_index)),2)
  WIDGET_CONTROL, string1, SET_VALUE=title(igraph)
  WIDGET_CONTROL, string2, SET_VALUE=xtitle(igraph)
  WIDGET_CONTROL, string3, SET_VALUE=ytitle(igraph)
  WIDGET_CONTROL, string_label, SET_VALUE=over_label(igraph,over_index)
END
;========================================================================
PRO plot_one_graph, i
@idl_post_comm

  IF xmax_auto THEN BEGIN
    max_val(0,i)=1e30
    max_val(1,i)=-1e30
  ENDIF

  CASE plot_type OF
  0: BEGIN
     x=FLTARR(N_ELEMENTS(data_nstep),over_max)
     y=FLTARR(N_ELEMENTS(data_nstep),over_max)
     ENDCASE
  1: BEGIN
     x=FLTARR(maxlevs,over_max)
     y=FLTARR(maxlevs,over_max)
     ENDCASE
  2: BEGIN
     x=FLTARR(maxlevs,over_max)
     y=FLTARR(maxlevs,over_max)
     ENDCASE
  ENDCASE

  FOR j=0,over_max-1 DO BEGIN
  IF (over_active(i,j)) THEN BEGIN

; get data and max/min for each plot type

  CASE plot_type OF
  0: BEGIN
     x(*,j)=data_nstep*timestep[0]/1440
     y(*,j)=data_val(index(i,j),lev(i,j),*)
     IF xmax_auto THEN BEGIN
       max_val(0,i)=TS_TO_DAY(MIN([data_nstep, max_val(0,i)]))
       max_val(1,i)=TS_TO_DAY(MAX([data_nstep, max_val(1,i)]))
     ENDIF
     IF ymax_auto THEN BEGIN
       IF (j eq 0) THEN BEGIN
         max_val(2,i)=1e30
         max_val(3,i)=-1e30
       ENDIF
       mm=MIN(data_val(index(i,j),lev(i,j),*))
       max_val(2,i)=MIN([mm,max_val(2,i)])
       mm=MAX(data_val(index(i,j),lev(i,j),*))
       max_val(3,i)=MAX([mm,max_val(3,i)])
     ENDIF
     ENDCASE
  1: BEGIN
     x(*,j)=data_val(index(i,j),*,time_index)
     CASE press_plot OF
       0: y(*,j)=findgen(maxlevs)+1
       1: y(*,j)=(ak+bk*EXP(data_val(press_index,0,time_index)))/100.
     ENDCASE
     IF xmax_auto THEN BEGIN
       mm=MIN(data_val(index(i,j),*,time_index))
       max_val(0,i)=MIN([mm,max_val(0,i)])
       mm=MAX(data_val(index(i,j),*,time_index))
       max_val(1,i)=MAX([mm,max_val(1,i)])
     ENDIF
     IF ymax_auto THEN BEGIN
       IF (j eq 0) THEN BEGIN
         max_val(2,i)=-1e30
         max_val(3,i)=1e30
       ENDIF
       max_val(2,i)=MAX([y(*),max_val(2,i)])
       max_val(3,i)=MIN([y(*),max_val(3,i)])
     ENDIF
     ENDCASE
  2: BEGIN
     x(*,j)=data_av(index(i,j),*)
     CASE press_plot OF
       0: y(*,j)=findgen(maxlevs)+1
       1: y(*,j)=(ak+bk*EXP(data_av(press_index,0)))/100.
     ENDCASE
     IF xmax_auto THEN BEGIN
       mm=MIN(data_av(index(i,j),*))
       max_val(0,i)=MIN([mm,max_val(0,i)])
       mm=MAX(data_av(index(i,j),*))
       max_val(1,i)=MAX([mm,max_val(1,i)])
     ENDIF
     IF ymax_auto THEN BEGIN
       IF (j eq 0) THEN BEGIN
         max_val(2,i)=-1e30
         max_val(3,i)=1e30
       ENDIF
;       max_val(0,i)=MIN([data_av(index(i,j),*),max_val(0,i)])
;       max_val(1,i)=MAX([data_av(index(i,j),*),max_val(1,i)])
       max_val(2,i)=MAX([y(*),max_val(2,i)])
       max_val(3,i)=MIN([y(*),max_val(3,i)])
     ENDIF
     ENDCASE
  ENDCASE

  ENDIF   ;overplot active
  ENDFOR  ;overplot loop

; plot the base graph

  PLOT, x(*,0), y(*,0), $
     thick=thick(i), ystyle=ystyle(i),xstyle=xstyle(i), $
     xthick=xthick(i), ythick=ythick(i), $
     xtitle=xtitle(i), ytitle=ytitle(i), title=title(i),$
     xrange=[max_val(0,i),max_val(1,i)], $
     yrange=[max_val(2,i),max_val(3,i)],$
     charthick=charthick(i), charsize=charsize(i)

; over plot if necessary

  dummy = WHERE(over_active(i,*) eq 1, count)

  IF (count gt 1) THEN BEGIN
    FOR J=1,over_max-1 DO IF (over_active(i,j)) THEN $
      OPLOT, x(*,j),y(*,j),thick=thick(i), $
      LINESTYLE=j*(1-colour),COLOR=j*colour+1

    yoff=0.1*(max_val(3,i)-max_val(2,i))
    xlen=0.1*(max_val(1,i)-max_val(0,i))
    xt=(label_x(i)-plot_win[0])/(plot_win[2]-plot_win[0])* $
         (max_val(1,i)-max_val(0,i))+max_val(0,i)
    yt=(label_y(i)-plot_win[1])/(plot_win[3]-plot_win[1])* $
         (max_val(3,i)-max_val(2,i))+max_val(2,i)

    FOR J=0,over_max-1 DO $
    IF (over_active(i,j)) THEN BEGIN
      OPLOT, [xt,xt+xlen],[yt,yt],thick=thick(i),$
        LINESTYLE=j*(1-colour),COLOR=j*colour+1*(j ne 0)
      XYOUTS, xt+xlen*1.5, yt, over_label(i,j), $
        charthick=charthick(i),charsize=charsize(i)*0.65
      yt=yt-yoff
    ENDIF
  ENDIF
END
;=======================================================================
PRO plot_graphs
@idl_post_comm

FOR i=0,ngraphs-1 DO BEGIN
  WIDGET_CONTROL, graph_win(i), GET_VALUE = win
  WSET, win
  plot_one_graph, i
ENDFOR

END
;=======================================================================
PRO make_postscript
  @idl_post_comm

  WIDGET_CONTROL, error0 , SET_VALUE='making postscript file '+filename

  SET_PLOT, 'PS'

  nrows=ngraphs/ngraph_cols
  !p.multi = [0,ngraph_cols,nrows]
  !p.font=0
;  DEVICE, /portrait,/helvetica,/bold,/color, file=filename
  DEVICE, /landscape,/helvetica,/bold,color=colour, file=filename
  red = [0, 1, 1, 0, 0, 1]     ;specify the red component of each color.
  green = [0, 1, 0, 1, 0, 1]  ;specify the green component of each color.
  blue = [0, 1, 0, 0, 1, 0]   ;Specify the blue component of each color.
  TVLCT, 255 * red, 255 * green, 255 * blue;load the first six elements of the color table.

  FOR i=0,ngraphs-1 DO plot_one_graph, i

  DEVICE, /CLOSE
  SET_PLOT,'X'
  !p.multi=0
  !p.font=-1
  red = [0, 1, 1, 0, 0, 1]     ;specify the red component of each color.
  green = [0, 1, 0, 1, 0, 1]  ;specify the green component of each color.
  blue = [0, 1, 0, 0, 1, 0]   ;Specify the blue component of each color.
  TVLCT, 255 * red, 255 * green, 255 * blue;load the first six elements of the color table.

END
;=======================================================================
PRO make_tephi
  @idl_post_comm

  WIDGET_CONTROL,error0,SET_VALUE='Must do MAKE_AVERAGE first - AVERAGE '+$
     'saved to DATA.DAT. Use tephi.pro to plot better tephigram'

  OPENW, 1, 'DATA.DAT'
  PRINTF, 1,' Tephigram file for ECHAM model - average between day'$ 
     +STR(TS_TO_DAY(av_min))+$
     ' and '+STR(TS_TO_DAY(av_max))
  PRINTF, 1,' Created by idl_post '+version+' on '+SYSTIME(0)
  PRINTF, 1,' Use tephi.pro to plot a tephigram using this file'
  PRINTF, 1,' Contact Adrian Tompkins at tompkins@dkrz.de '

  FOR i=0,maxlevs-1 DO PRINTF, 1, format='(4(f16.8,x))', $
   (ak(i)+bk(i)*EXP(data_val(press_index,0,time_index)))/100., 9999.0,$
    data_av[t_index,i], MAX([data_av[q_index,i]*1000.,0.00000001])

  CLOSE, 1
END
;=======================================================================
PRO plot_tephi
  @idl_post_comm
  i=FINDFILE('DATA.DAT',COUNT=c)
  IF (c eq 0) THEN BEGIN
    print, STRING(7B) ;ring bell
    WIDGET_CONTROL, error0, SET_VALUE='Error: make tephigram before trying to plot'
  ENDIF ELSE BEGIN
    tephi_auto
    !p.font=-1
    SPAWN, 'gv TEPHI_AUTO.ps &'
  ENDELSE
END
;=========================================================================
PRO load_default
  @idl_post_comm

  WIDGET_CONTROL, error0, SET_VALUE='loading default settings from idl_post.default'
  OPENR, 1, 'idl_post.default', ERROR=err
  IF err ne 0 THEN BEGIN
    print, STRING(7B) ;ring bell
    WIDGET_CONTROL,error0, SET_VALUE='deault file not present'
    RETURN
  ENDIF

  READF, 1, i
  IF (i ne number_list_max) THEN BEGIN
    print, 'incomptible files, you have changed number_list_max'
    CLOSE,1 
    RETURN
  ENDIF

  READF, 1,   $
  thick,      $
  xthick,     $
  ythick,     $
  charthick,  $
  charsize,   $
  xstyle,     $
  ystyle

  i=WHERE(number_list eq thick(0))
  WIDGET_CONTROL, list_thick, SET_DROPLIST_SELECT=i(0)
  i=WHERE(number_list eq xthick(0))
  WIDGET_CONTROL, list_xthick, SET_DROPLIST_SELECT=i(0)
  i=WHERE(number_list eq ythick(0))
  WIDGET_CONTROL, list_ythick, SET_DROPLIST_SELECT=i(0)
  i=WHERE(char_list eq charthick(0))
  WIDGET_CONTROL, list_charthick, SET_DROPLIST_SELECT=i(0)
  i=WHERE(char_list eq charsize(0))
  WIDGET_CONTROL, list_charsize, SET_DROPLIST_SELECT=i(0)

  CLOSE, 1
END
;=========================================================================
PRO save_default
  @idl_post_comm

  WIDGET_CONTROL, error0, SET_VALUE='saving default settings to idl_post.default'

  OPENW, 1, 'idl_post.default'

  PRINTF, 1, number_list_max

  PRINTF, 1,  $
  thick,      $
  xthick,     $
  ythick,     $
  charthick,  $
  charsize,   $
  xstyle,     $
  ystyle

  CLOSE, 1
END
;=========================================================================
PRO save_plots
  @idl_post_comm

  WIDGET_CONTROL, error0, SET_VALUE='saving graph selections to idl_post.plots'

  OPENW, 1, 'idl_post.plots'
  PRINTF, 1, lev, $
             over_active, $
             colour, $
             index, $
             time_index, $
             plot_type, $
             xstyle, $
             ystyle, $
             xmax_auto, $
             ymax_auto, $
             ngraphs

;  PRINTF, 1, title, $
;             xtitle, $
;             ytitle, $
;             over_label
  CLOSE, 1
END
;=========================================================================
PRO load_plots
  @idl_post_comm

  WIDGET_CONTROL, error0, SET_VALUE='loading graph selections from idl_post.plots'

  OPENR, 1, 'idl_post.plots' 
  READF, 1,  lev, $
             over_active, $
             colour, $
             index, $
             time_index, $
             plot_type, $
             xstyle, $
             ystyle, $
             xmax_auto, $
             ymax_auto, $
             ngraphs

;  READF, 1,  title, $
;             xtitle, $
;             ytitle, $
;             over_label
  CLOSE, 1
END
;==========================================================================
PRO list_event, event
  @idl_post_comm

  WIDGET_CONTROL, event.id, GET_UVALUE = widgetName     

; get graph number clicked on
  IF (STRMID(widgetName,0,5) eq 'graph') THEN BEGIN
    igraph = FIX(STRMID(widgetName,5,1))
    widgetName = 'graph'
  ENDIF

  CASE widgetName OF
     "average": BEGIN
        WIDGET_CONTROL, error0 , SET_VALUE='Averaging between timesteps ' $
          +STR(ROUND(av_min))+' and '+STR(ROUND(av_max))
        IF (av_max eq av_min) THEN $
          data_av=data_val[*,*,av_min] $
          ELSE data_av=TOTAL(data_val[*,*,av_min:av_max],3)/(av_max-av_min+1)
;        auto_title
        ENDCASE
     "thick"    : thick(ch_graph)      = number_list(event.index)
     "xthick"   : xthick(ch_graph)     = number_list(event.index)
     "ythick"   : ythick(ch_graph)     = number_list(event.index)
     "charthick": charthick(ch_graph)  = char_list(event.index)
     "charsize" : charsize(ch_graph)   = char_list(event.index)

     "lev": BEGIN
        lev(ch_graph,*) = event.index 
        auto_title,1,0,0
        ENDCASE

     "over_index": BEGIN
        over_index = event.value
        over_active(ch_graph,over_index)=1
        WIDGET_CONTROL,  bgroup_over_active,  SET_VALUE=over_active(igraph,*)
        ENDCASE

     "over_active": BEGIN
        over_active(ch_graph,event.value) = $
           (over_active(ch_graph,event.value) eq 0)
;        auto_title
        ENDCASE

     "colour": colour = event.value

     "label": BEGIN
        WIDGET_CONTROL, error0, SET_VALUE='click on top left to place labels'
        label=1
        ENDCASE

     "index": BEGIN
        index(igraph,over_index) = data_sort(event.index)
        auto_title,1,1,1
        ENDCASE

     "title":BEGIN
        WIDGET_CONTROL, string1, GET_VALUE=string
        title(ch_graph)=string
        WIDGET_CONTROL, error0, SET_VALUE= $
          'Tip: Many functions change titles automatically: '+ $
          'set just before making ps file'
        ENDCASE

     "xtitle":BEGIN
        WIDGET_CONTROL, string2, GET_VALUE=string
        xtitle(ch_graph)=string
        WIDGET_CONTROL, error0, SET_VALUE= $
          'Tip: Many functions change titles automatically: '+ $
          'set just before making ps file'
        ENDCASE

     "ytitle":BEGIN
        WIDGET_CONTROL, string3, GET_VALUE=string
        ytitle(ch_graph)=string
        WIDGET_CONTROL, error0, SET_VALUE= $
          'Tip: Many functions change titles automatically: '+ $
          'set just before making ps file'
        ENDCASE

     "over_label":BEGIN
        WIDGET_CONTROL, string_label, GET_VALUE=string
        over_label(ch_graph,over_index)=string
        WIDGET_CONTROL, error0, SET_VALUE= $
          'Tip: Many functions change titles automatically: '+ $
          'set just before making ps file'
        ENDCASE

     "time_index":BEGIN
        time_index = GET_NUMBER(string4,error)
        IF NOT error THEN BEGIN
          time_index = DAY_TO_TS(time_index)
          time_index = MIN([time_index,ntimes-1])
          time_index = MAX([time_index,0])
          WIDGET_CONTROL, error0 , SET_VALUE='equivalent to timestep '+$
            STR(time_index+data_nstep(0),FMT='(I6)')
          WIDGET_CONTROL, string4 , SET_VALUE=STR(TS_TO_DAY(time_index))
          auto_title,1,0,0
        ENDIF
        ENDCASE

     "stepforward": BEGIN
        time_index = time_index+1
        time_index = MIN([time_index,ntimes-1])
        WIDGET_CONTROL, error0 , SET_VALUE='timestep '+$
            STR(time_index+data_nstep(0),FMT='(I6)')
        WIDGET_CONTROL, string4 , SET_VALUE=STR(TS_TO_DAY(time_index))
        auto_title,1,0,0
        ENDCASE

     "timestep": timestep = GET_NUMBER(string5,error)

     "av_min":BEGIN
        av_min = GET_NUMBER(string6,error)
        IF NOT error THEN BEGIN
          av_min = DAY_TO_TS(av_min)
          av_min = MIN([av_min,ntimes-1])
          av_min = MAX([av_min,0])
          WIDGET_CONTROL, error0 , SET_VALUE='equivalent to timestep '+$
            STR(av_min+data_nstep(0),FMT='(I4)')
          WIDGET_CONTROL, string6 , SET_VALUE=STR(TS_TO_DAY(av_min))
        ENDIF
        ENDCASE

     "av_max":BEGIN
        av_max = GET_NUMBER(string7,error)
        IF NOT error THEN BEGIN
          av_max = DAY_TO_TS(av_max)
          av_max = MIN([av_max,ntimes-1])
          av_max = MAX([av_max,av_min,0])
          WIDGET_CONTROL, error0 , SET_VALUE='equivalent to timestep '+$
            STR(av_max+data_nstep(0),FMT='(I4)')
          WIDGET_CONTROL, string7 , SET_VALUE=STR(TS_TO_DAY(av_max))
        ENDIF
        ENDCASE

     "reload": read_data

     "dump": dump_graph

     "exit": BEGIN
        WIDGET_CONTROL, event.top, /DESTROY
        STOP
        ENDCASE

     "global": IF (event.value) THEN BEGIN
          ch_graph=indgen(number_list_max)
        ENDIF ELSE ch_graph=igraph

     "plot_type": BEGIN
        plot_type = event.value
        auto_title,1,1,1
        ENDCASE

     "xstyle": xstyle(ch_graph) = event.value
     "ystyle": ystyle(ch_graph) = event.value
     "xmax_auto": xmax_auto = event.value
     "ymax_auto": ymax_auto = event.value
     "max_index": max_index = event.value
     "max_val": BEGIN
        val = GET_NUMBER(string8,error)
        IF NOT error THEN max_val(max_index,ch_graph)=val
        ENDCASE

     "pressure": BEGIN
          IF (press_index ne -1) THEN press_plot = event.index $
          ELSE BEGIN
            print, STRING(7B) ;ring bell
            WIDGET_CONTROL, error0, SET_VALUE='Error: pressure data not available'
          ENDELSE
        auto_title,0,0,1       
        ENDCASE

     "default": BEGIN
        CASE event.index OF
          0: load_default
          1: save_default
          ENDCASE
        ENDCASE
     "saveplots": BEGIN
        CASE event.index OF
          0: BEGIN
             load_plots	
             WIDGET_CONTROL,col_right1 , /DESTROY, UPDATE=0
;             ngraphs = number_list(event.index)
             if (ngraphs gt 5) THEN ngraph_cols=(ngraphs-1)/3+1 ELSE ngraph_cols=1
       
             nrows=ngraphs/ngraph_cols
             win_xsize=graph_xsize/ngraph_cols
             win_ysize=graph_ysize/nrows
       
             col_right1 = WIDGET_BASE(row_base1, row=nrows+1)
             igraph = MIN([igraph,ngraphs-1])       ; reset igraph if out of range
             ; not global => change
             IF (N_ELEMENTS(ch_graph) eq 1) THEN ch_graph=igraph  
             graph_win=intarr(ngraphs)
             FOR i=0,ngraphs-1 DO graph_win(i) = $
               WIDGET_DRAW(col_right1, XSIZE=win_xsize, $
                 YSIZE=win_ysize,RET=2, $
                 FRAME=frame_thick*(i eq igraph), $
                 UVALUE='graph'+STR(i),/BUTTON_EVENTS,COLORS=10)
             error0 = WIDGET_LABEL(col_right1, Value="messages will appear here",$
               xsize=graph_xsize) 
             WIDGET_CONTROL, col_right1, /REALIZE, UPDATE=1
             WIDGET_CONTROL, graph_win(igraph), GET_VALUE = graph_wcur
             WSET, graph_wcur
             ENDCASE
          1: save_plots
          ENDCASE
        ENDCASE
     "postscript": BEGIN
        CASE event.index OF
          0: make_postscript
          1: SPAWN, 'gv -swap '+filename+' &'
          ENDCASE
        ENDCASE
     "tephigram": BEGIN
        IF (t_index ne -1 and q_index ne -1) THEN BEGIN
          CASE event.index OF
            0: make_tephi
            1: plot_tephi
            ENDCASE
        ENDIF ELSE BEGIN
          print, STRING(7B) ;ring bell
          WIDGET_CONTROL, error0, SET_VALUE='Error: T or Q data not available'
        ENDELSE
        ENDCASE
     "ngraphs":BEGIN
        WIDGET_CONTROL,col_right1 , /DESTROY, UPDATE=0
        ngraphs = number_list(event.index)
        if (ngraphs gt 5) THEN ngraph_cols=(ngraphs-1)/3+1 ELSE ngraph_cols=1

        nrows=ngraphs/ngraph_cols
        win_xsize=graph_xsize/ngraph_cols
        win_ysize=graph_ysize/nrows

        col_right1 = WIDGET_BASE(row_base1, row=nrows+1)
        igraph = MIN([igraph,ngraphs-1])       ; reset igraph if out of range
        ; not global => change
        IF (N_ELEMENTS(ch_graph) eq 1) THEN ch_graph=igraph  
        graph_win=intarr(ngraphs)
        FOR i=0,ngraphs-1 DO graph_win(i) = $
          WIDGET_DRAW(col_right1, XSIZE=win_xsize, $
            YSIZE=win_ysize,RET=2, $
            FRAME=frame_thick*(i eq igraph), $
            UVALUE='graph'+STR(i),/BUTTON_EVENTS,COLORS=10)
        error0 = WIDGET_LABEL(col_right1, Value="messages will appear here",$
          xsize=graph_xsize) 
        WIDGET_CONTROL, col_right1, /REALIZE, UPDATE=1
        WIDGET_CONTROL, graph_win(igraph), GET_VALUE = graph_wcur
        WSET, graph_wcur
        auto_title,1,1,1
        ENDCASE

     "graph": BEGIN
        IF (event.press eq 1) THEN BEGIN
          ;not global => change
          IF (N_ELEMENTS(ch_graph) eq 1) THEN ch_graph=igraph  
          IF (label) THEN BEGIN    ; positioning labels
            label_x(ch_graph)=float(event.x)/win_xsize
            label_y(ch_graph)=float(event.y)/win_ysize
            label=0
          ENDIF
          WIDGET_CONTROL, graph_win(igraph), GET_VALUE = graph_wcur
          WSET, graph_wcur
          IF (lev(igraph) eq 0) THEN BEGIN   ;hasn't been used before
            lev(igraph) = lev(0,0)
            index(igraph) = index(0,0)
          ENDIF
          WIDGET_CONTROL,col_right1 , /DESTROY, UPDATE=0

          nrows=ngraphs/ngraph_cols
          win_xsize=graph_xsize/ngraph_cols
          win_ysize=graph_ysize/nrows

          col_right1 = WIDGET_BASE(row_base1, row=nrows+1)
          FOR i=0,ngraphs-1 DO graph_win(i) = $
          WIDGET_DRAW(col_right1, XSIZE=win_xsize, $
            YSIZE=win_ysize,RET=2, $
            FRAME=frame_thick*(i eq igraph), $
            UVALUE='graph'+STR(i),/BUTTON_EVENTS,COLORS=10)
          error0 = WIDGET_LABEL(col_right1, $
            Value="messages will appear here",$
            xsize=graph_xsize)
          WIDGET_CONTROL, base, /REALIZE, UPDATE=1
          WIDGET_CONTROL, bgroup_over_active, SET_VALUE=over_active(igraph,*)
        ENDIF
        ENDCASE
     ELSE: $ ; When an event occurs in a widget that has no user value in this
              ; case statement, an error message is shown
         print, "Code Error: Event User Value Not Found"      
  ENDCASE

  plot_graphs

END

PRO CleanUp, base
print, 'bye!'
END

;-----------------------------------------------------------------------
PRO default_vals

; set up default arrays and values
; user defined vars are also here

  @idl_post_comm

;------------------------------
; CONTROL VARIABLES
;------------------------------
  max=8000
  nlev=50
  ngraphs=1
  number_list_max=6
  over_max=4

;------------------------------
; SET UP ARRAYS
;------------------------------
  graph_win =INTARR(ngraphs)

  xstyle  =INTARR(number_list_max)
  ystyle  =INTARR(number_list_max)
  thick   =INTARR(number_list_max)
  xthick  =INTARR(number_list_max)
  ythick  =INTARR(number_list_max)

  charthick =FLTARR(number_list_max)
  charsize  =FLTARR(number_list_max)

  label_x =FLTARR(number_list_max)
  label_y =FLTARR(number_list_max)

  lev     =INTARR(number_list_max,over_max)
  index   =INTARR(number_list_max,over_max)
  over_active =INTARR(number_list_max,over_max)
  over_label=STRARR(number_list_max,over_max)
  title   =STRARR(number_list_max)
  xtitle  =STRARR(number_list_max)
  ytitle  =STRARR(number_list_max)
  max_val =FLTARR(4,number_list_max) ; min max x/y

  igraph=0
  ch_graph = igraph

;------------------------------
; DEFAULT PLOTTING OPTIONS
;------------------------------
  ystyle(*)=1
  xstyle(*)=0
  thick(*)=1.
  xthick(*)=1.
  ythick(*)=1.
  title(*)='title'
  xtitle(*)='Time (days)'
  ytitle(*)='ytitle'
  charthick(*)=1. 
  charsize(*)=1.

  colour=0
  plot_type=0
  time_index=0
  timestep=20

  filename ='idl_post.ps'

  graph_xsize=750. ;default window size when 
  graph_ysize=750. ;ngraph=1

  plot_win=[0.1,0.2,0.9,0.8] ; the proportion of graphic window
                                ; taken by plot
  frame_thick=7 

  index(*)=0
  lev(*)=1
  over_index=0
  over_active(*,0)=1

  av_min=0
  av_max=0
;------------------------------
; END DEFAULT PLOTTING OPTIONS
;------------------------------
  if (ngraphs gt 5) THEN ngraph_cols=(ngraphs-1)/3+1 ELSE ngraph_cols=1

END

;-----------------------------------------------------------------------
; READ DATA     
;-----------------------------------------------------------------------
PRO read_data

  CLOSE, /ALL
  @idl_post_comm

  label = 0
  over_index = 0
  maxlevs=0
  xmax_auto=1
  ymax_auto=1
  max_index=0
  press_index=-1
  t_index=-1
  q_index=-1
  press_plot = 0
  name='1234567890123456'
  rank=lonarr(1)
  shap=lonarr(3)

  nitems=60
  data_name=STRARR(nitems)
  data_nstep=FLTARR(max)
  data_val=FLTARR(nitems,nlev,max)

  nitems=0

  OPENR, 1, 'result2', ERROR=err, /F77_UNFORMATTED
  WHILE (err NE 0) DO BEGIN
    print, STRING(7B) ;ring bell
    read, 'result2 not present, enter data file name: ',string
    OPENR, 1, string, ERROR=err, /F77_UNFORMATTED
  ENDWHILE
  j=long(0)
  print, 'reading data'
  WHILE (eof(1) ne 1) DO BEGIN
    READU, 1, name, rank, shap
    siz=1
    FOR i=0,rank(0)-1 DO siz=siz*shap(i)
    IF (siz eq 0) THEN siz=1
    maxlevs = MAX([maxlevs,siz])
    z=DBLARR(shap(0),siz)
    READU, 1,  z

; FOUND START OF DATA TIMESTEP
    IF (STRMID(name,0,5) eq 'NSTEP') THEN BEGIN
      data_nstep(j)=z(0)
      j=j+1
      print, j
    ENDIF ELSE BEGIN
      found=0
      FOR i=0,nitems-1 DO $
        IF (data_name(i) eq name) THEN BEGIN
          data_val(i,0:siz-1,j)=z   ;found index - store data
          found=1
        ENDIF
      IF NOT found THEN BEGIN ; first time this variable read
        print, name
        IF (STRMID(name,0,4) eq 'ALPS') THEN press_index=nitems
        IF (STRMID(name,0,2) eq 'T ') THEN t_index=nitems
        IF (STRMID(name,0,2) eq 'Q ') THEN q_index=nitems
        data_name(nitems)=name        ; store name
        data_val(nitems,0:siz-1,j)=z  ; store data
        nitems=nitems+1               ; increase index
      ENDIF
    ENDELSE
  endwhile
  close, 1

;----------------------------------------------------------------------- 
; NOTE: distcard first 2 timesteps since they contains rubbish
;-----------------------------------------------------------------------
  ntimes=j-2  ;distcard first timestep

; CUT ARRAYS OFF TO CORRECT SIZE
  data_nstep=data_nstep(2:ntimes+1)
  data_val=data_val(0:nitems-1,0:maxlevs-1,2:ntimes+1)
  data_name=data_name(0:nitems-1)
;  data_sort=SORT(data_name)
  data_sort=findgen(nitems)

; MAKE ARRAY FOR AVERAGES
  data_av=FLTARR(nitems,maxlevs)

  IF (press_index eq -1) THEN BEGIN
    print, STRING(7B) ;ring bell
    print, 'Pressure not found - Should be called ALPS in results file'
    print, 'plots against pressure *NOT* possible, or the tephigram plots'
  ENDIF
  IF (t_index eq -1 or q_index eq -1) THEN BEGIN
    print, STRING(7B) ;ring bell
    print,'temperature or moisture not found - Should be called T & Q in results file'
    print,'tephigram plots *NOT* possible'
  ENDIF
  
  print, 'finished reading data'
  print, 'now reading forcing file for pressure'

;------------------------------
; READ FORCING
;------------------------------
  OPENr, 1, 'forcing', ERROR=err, /f77_unformatted
  WHILE (err NE 0) DO BEGIN
    print, STRING(7B) ;ring bell
    read, 'forcing file not present, enter data file name: ',string
    OPENR, 1, string, ERROR=err, /F77_UNFORMATTED
  ENDWHILE
  erra=1
  errb=1
  j=long(0)
  print, 'reading data'
  WHILE (eof(1) ne 1) DO BEGIN
    READU, 1, name, rank, shap
    siz=1
    print, name
    FOR i=0,rank(0)-1 DO siz=siz*shap(i)
    IF (siz eq 0) THEN siz=1
    z=DBLARR(shap(0),siz)
    READU, 1,  z

; FOUND START OF DATA TIMESTEP
    IF (erra AND STRMID(name,0,2) eq 'AK') THEN BEGIN
      akf=z
      erra=0
    ENDIF
    IF (errb AND STRMID(name,0,2) eq 'BK') THEN BEGIN
      bkf=z
      errb=0
    ENDIF
    IF (STRMID(name,0,5) eq 'DTIME') THEN BEGIN
      timestep=z/60.
    ENDIF
  ENDWHILE
  CLOSE, 1

  IF (erra OR errb) THEN BEGIN
    press_index = -1
    print, STRING(7B) ;ring bell
    print, 'As and Bs for defining pressure levels absent from forcing file'
    print, 'They should be labelled AK and BK'
    print, 'plots against pressure *NOT* possible, or the tephigram plots'
  ENDIF ELSE BEGIN
    ak=FLTARR(maxlevs)
    bk=FLTARR(maxlevs)
    FOR i=0,maxlevs-1 DO BEGIN
      ak(i)=0.5*(akf(i+1)+akf(i))
      bk(i)=0.5*(bkf(i+1)+bkf(i))
    ENDFOR
  ENDELSE
END

;---------------------------------------------------------------------------
PRO screen_stuff, GROUP=group

  @idl_post_comm

  version='v1.1-18.05.00 '

  base = WIDGET_BASE(TITLE='IDL_POST '+version+' - report bugs or suggestions to ' +$
    'Adrian at tompkins@dkrz.de', XOFFSET=10, YOFFSET=10 )

; SET UP WIDGET BASES

  row_base1 = WIDGET_BASE(base, row=1)
  col_left1 = WIDGET_BASE(row_base1, column=2)
  col_right1 = WIDGET_BASE(row_base1, row=2)

  number_list = indgen(number_list_max)+1
  char_list   = findgen(16)*0.1+0.5
  lev_list    = indgen(maxlevs)+1
  over_list   = indgen(over_max)+1
  title_list  = ['Title','X-title','Y-title']


; DATA CONTROL WIDGETS

;  lab1 = WIDGET_LABEL(FONT='6X10',col_left1, Value="Data:", FRAME=1 ) 
  list7 = WIDGET_LIST (col_left1, VALUE=data_name[data_sort], $
    UVALUE='index',ysize=5)
  list6 = WIDGET_DROPLIST( col_left1, VALUE=STRING(lev_list), $
    UVALUE='lev', TITLE='lev:')
;  list7 = WIDGET_DROPLIST (col_left1, VALUE=data_name(data_sort), $
;    UVALUE='index',TITLE='var:')
  list8 = WIDGET_DROPLIST( col_left1, VALUE=STRING(number_list), $
    UVALUE='ngraphs', TITLE='ngraphs:')
  bgroup_over_index = CW_BGROUP(col_left1, STR(over_list), $
    /ROW, /EXCLUSIVE, $
    LABEL_LEFT='Overplot: ', UVALUE='over_index', $
    /RETURN_INDEX, /NO_RELEASE)
  bgroup_over_active = CW_BGROUP(col_left1, STR(over_list), $
    /ROW,/NONEXCLUSIVE, $
    LABEL_LEFT='Activate: ', UVALUE='over_active', $
    /RETURN_INDEX)
  WIDGET_CONTROL,bgroup_over_index, SET_VALUE=over_index 
  WIDGET_CONTROL,bgroup_over_active, SET_VALUE=over_active(0,*)

; TIME SLICES AND AVERAGING NUMBERS

  bgroup_plot_type = CW_BGROUP(col_left1, ['level','time','average'], $
    /ROW, /EXCLUSIVE, $
    LABEL_TOP='Plot data at: ', UVALUE='plot_type', $
    /RETURN_INDEX, /NO_RELEASE)
  WIDGET_CONTROL,bgroup_plot_type, SET_VALUE=0 
  lab0 = WIDGET_LABEL(FONT='6X10',col_left1, Value=" ") 
  lab4 = WIDGET_LABEL(FONT='6X10',col_left1, Value="Time slice (days):", FRAME=1) 
  string4  = WIDGET_TEXT (col_left1, UVALUE='time_index', $
    VALUE='0',XSIZE=10, YSIZE=1, /EDITABLE)
  button_timestep  = WIDGET_BUTTON (FONT='6X10',col_left1, $
     VALUE='step forward',UVALUE='stepforward')

  lab0 = WIDGET_LABEL(FONT='6X10',col_left1, Value=" ") 
  button1  = WIDGET_BUTTON (FONT='6X10',col_left1, VALUE='make average',UVALUE='average')
  string6  = WIDGET_TEXT (FONT='6X10',col_left1,UVALUE='av_min',$
     VALUE='type minimum in here', $
     XSIZE=10, YSIZE=1, /EDITABLE)
  string7  = WIDGET_TEXT (FONT='6X10',col_left1,UVALUE='av_max', $
     VALUE='type maximum in here', $
     XSIZE=10, YSIZE=1, /EDITABLE)

; FILE WIDGETS

  lab0 = WIDGET_LABEL(FONT='6X10',col_left1, Value=" ") 
  lab3 = WIDGET_LABEL(FONT='6X10',col_left1, Value="Files:", FRAME=1 ) 
  list9 = WIDGET_DROPLIST (FONT='6X10',col_left1, $
     VALUE=['Load default','Save default'], $
     UVALUE='default', TITLE='Save/load:')
  list9b = WIDGET_DROPLIST (FONT='6X10',col_left1, $
     VALUE=['Load plots','Save plots'], $
     UVALUE='saveplots')
  list10 = WIDGET_DROPLIST (FONT='6X10',col_left1, $
     VALUE=['Make postscript','View postscript'], $
     UVALUE='postscript', TITLE='Postscript:')
  list11 = WIDGET_DROPLIST (FONT='6X10',col_left1, $
     VALUE=['Make tephi file','View tephigram'], $
     UVALUE='tephigram', TITLE='Tephigram:')

; axis styles
  bgroup_xstyle = CW_BGROUP(FONT='6X10',col_left1, ['x:approx','x:exact'], $
    /ROW, /EXCLUSIVE, $
    LABEL_TOP='axis styles ', UVALUE='xstyle',/FRAME, $
    /RETURN_INDEX ,/NO_RELEASE)
  bgroup_ystyle = CW_BGROUP(FONT='6X10',col_left1, ['y:approx','y:exact'], $
    /ROW, /EXCLUSIVE, UVALUE='ystyle',/FRAME, $
    /RETURN_INDEX ,/NO_RELEASE)
  WIDGET_CONTROL,bgroup_xstyle,SET_VALUE=xstyle(0) 
  WIDGET_CONTROL,bgroup_ystyle,SET_VALUE=ystyle(0) 

; ATTRIBUTES WIDGETS

  lab0 = WIDGET_LABEL(col_left1, Value=" ") 
  bgroup_color = CW_BGROUP(col_left1, ['linestyle','colour'], $
    /ROW, /EXCLUSIVE, UVALUE='colour', $
    /RETURN_INDEX, /NO_RELEASE)
  bgroup1 = CW_BGROUP(col_left1, ['locally','globally'], /ROW, /EXCLUSIVE, $
    LABEL_TOP='Apply Attributes: ', UVALUE='global',/FRAME, $
    /RETURN_INDEX ,/NO_RELEASE)
  WIDGET_CONTROL,bgroup_color,SET_VALUE=colour
  WIDGET_CONTROL, bgroup1, SET_VALUE=0
  list_thick = WIDGET_DROPLIST( col_left1, VALUE=STRING(number_list), $
    UVALUE='thick', TITLE='thick:')
  list_xthick = WIDGET_DROPLIST( col_left1, VALUE=STRING(number_list), $
    UVALUE='xthick', TITLE='xthick:')
  list_ythick = WIDGET_DROPLIST( col_left1, VALUE=STRING(number_list), $
    UVALUE='ythick', TITLE='ythick:')
  list_charthick = WIDGET_DROPLIST( col_left1, $
    VALUE=STRING(char_list,FORMAT='(F4.1)'), $
    UVALUE='charthick', TITLE='charthick:')
  list_charsize = WIDGET_DROPLIST( col_left1, VALUE=STRING(char_list, $
    FORMAT='(F4.1)'), UVALUE='charsize', TITLE='charsize:')
  list12 = WIDGET_DROPLIST (col_left1, VALUE=['level','pressure'], $
         UVALUE='pressure', TITLE='Plot data vs ')

; MIN MAX AUTO

  bgroup_auto=CW_BGROUP(FONT='6X10',col_left1, ['X-Manual','X-Auto'], /ROW, /EXCLUSIVE,$
    LABEL_TOP='Max/Min Setting: ', UVALUE='xmax_auto',/FRAME, $
    /RETURN_INDEX ,/NO_RELEASE)
  WIDGET_CONTROL, bgroup_auto,SET_VALUE=xmax_auto  
  bgroup_auto=CW_BGROUP(FONT='6X10',col_left1, ['Y-Manual','Y-Auto'], /ROW, /EXCLUSIVE,$
    UVALUE='ymax_auto',/FRAME, $
    /RETURN_INDEX ,/NO_RELEASE)
  WIDGET_CONTROL, bgroup_auto,SET_VALUE=ymax_auto  

  bgroup_xymax=CW_BGROUP(col_left1, ['xmin','xmax','ymin','ymax'], $
    /ROW, /EXCLUSIVE, UVALUE='max_index',/FRAME, $
    /RETURN_INDEX ,/NO_RELEASE, FONT='6X10')
  string8  = WIDGET_TEXT (FONT='6X10',col_left1,UVALUE='max_val', $
       VALUE='type min/max in here', $
       XSIZE=10, YSIZE=1, /EDITABLE)

; TITLE WIDGETS

  lab0 = WIDGET_LABEL(FONT='6X10',col_left1, Value=" ") 
  lab2 = WIDGET_LABEL(FONT='6X10',col_left1, Value="Titles:", FRAME=1 ) 
  string1  = WIDGET_TEXT  (FONT='6X10',col_left1, UVALUE='title', $
     VALUE=title(igraph),XSIZE=10, YSIZE=1, /EDITABLE)
  string2  = WIDGET_TEXT  (FONT='6X10',col_left1, UVALUE='xtitle', $
     VALUE=xtitle(igraph),XSIZE=10, YSIZE=1, /EDITABLE)
  string3  = WIDGET_TEXT  (FONT='6X10',col_left1, UVALUE='ytitle', $
     VALUE=ytitle(igraph),XSIZE=10, YSIZE=1, /EDITABLE)
  button_label = WIDGET_BUTTON (FONT='6X10',col_left1, VALUE='Position Labels', $
     UVALUE='label')
  string_label  = WIDGET_TEXT  (FONT='6X10',col_left1, UVALUE='over_label', $
     VALUE=' ',XSIZE=10, YSIZE=1, /EDITABLE)

; TIMESTEP

  lab0 = WIDGET_LABEL(FONT='6X10',col_left1, Value=" ") 
  lab4 = WIDGET_LABEL(FONT='6X10',col_left1, Value="Model timestep (mins):", FRAME=1) 
  string5  = WIDGET_TEXT (FONT='6X10',col_left1,UVALUE='timestep', $
    VALUE=STR(timestep), $
    XSIZE=10, YSIZE=1, /EDITABLE)

; RELOAD DATA
  lab0 = WIDGET_LABEL(FONT='6X10',col_left1, Value=" ") 
  button_load = WIDGET_BUTTON (FONT='6X10',col_left1, VALUE='reload data',UVALUE='reload')

; DUMP TO FILE
  button_dump = WIDGET_BUTTON (FONT='6X10',col_left1, VALUE='dump graph to file',UVALUE='dump')

; EXIT BUTTON
  lab0 = WIDGET_LABEL(FONT='6X10',col_left1, Value=" ") 
  button0 = WIDGET_BUTTON (col_left1, VALUE='EXIT',UVALUE='exit')

; GRAPHIC WIDGETS ;

  nrows=ngraphs/ngraph_cols
  win_xsize=graph_xsize/ngraph_cols
  win_ysize=graph_ysize/nrows

  FOR i=0,ngraphs-1 DO graph_win(i) = WIDGET_DRAW(col_right1,/BUTTON_EVENTS, $
    XSIZE=win_xsize,YSIZE=win_ysize,RET=2,  $
    UVALUE='graph'+STR(i),COLORS=10)
  error0 = WIDGET_LABEL(col_right1, Value='messages will appear here',$
    xsize=graph_xsize)


; make sure char sizes etc start on number 1.0 as default
  i=WHERE(char_list eq charthick(igraph))
  WIDGET_CONTROL, list_charthick, SET_DROPLIST_SELECT=i(0)

  i=WHERE(char_list eq charsize(igraph))
  WIDGET_CONTROL, list_charsize, SET_DROPLIST_SELECT=i(0)

; REALIZE THE BASE WIDGET

  WIDGET_CONTROL, base, /REALIZE

; PLOT THE GRAPHS

  WIDGET_CONTROL, graph_win(igraph), GET_VALUE = graph_wcur
  WSET, graph_wcur
  plot_graphs

  XMANAGER, 'start',base, GROUP_LEADER = group, $
    EVENT_HANDLER = "list_event", CLEANUP="CleanUp", /NO_BLOCK

END

;-----------------------------------------------------------------------
; MAIN ENTRY POINT
;-----------------------------------------------------------------------
@idl_post_comm

; -------------
; RELEASE NOTES
; -------------
;
; v1.00:  First release 14.04.00
;
;         Multi window capability, with overplotting.
;         tephigram and hardplot facilities
; v1.01:  Bugs corrected to hardplot with > 5 windows
;
;         Auto max/min for average plot corrected
;         New data reload button added
;         correction to activate buttons
;         selected a new dataset only changes current graph even if 'global'
;         new dumpfile facility added for loading into XMGR
;         choice between color or linestyle
; V1.1    Release 2 18-5-00
;
;         bug corrected where relabling was lost when changing panel
;
;         read data from file 'result2' 
;         read timestep (DTIME) from 'forcing' (A.Rhodin 3-8-00)
;
COLORSET,decomposed=0
spawn, '/bin/rm DATA.DAT'
print, 'removed old tephigram data files'
device,decomposed=0

default_vals               ; set default values and define arrays
read_data                  ; read data set
screen_stuff,GROUP=group   ; set up widget base

red = [0, 1, 1, 0, 0, 1]     ;specify the red component of each color.
green = [0, 1, 0, 1, 0, 1]  ;specify the green component of each color.
blue = [0, 1, 0, 0, 1, 0]   ;Specify the blue component of each color.
TVLCT, 255 * red, 255 * green, 255 * blue;load the first six elements of the color table.

end

