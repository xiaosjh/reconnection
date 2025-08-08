;pro six_aia_plot
;  ;这个程序是画六个AIA子图的
;  dir='D:/Learning/PHD1st/magnetic_reconnecion/data/AIA2/131/'
;  file_list=file_search(dir+'*.fits')
;  j=88
;  file=file_list[j]
;  read_sdo,file,hdr,data
;  hmi_prep,hdr,data,hdr0,data0
;  index2map,hdr0,data0,map0
;  map0.data=map0.data/map0.dur
;  ;dmap=drot_map(map0,time='2024-06-18T21:17:00')
;  ;dmap.rtime=map0.time
;  !P.BACKGROUND=255
;  !P.COLOR=0
;  sub_map,map0,smap,xrange=[180,410],yrange=[-290,-520]
;  aia_lct,wavelnth='131',/load
;  window,/free,xsize=600,ysize=600
;  plot_map,smap,/log
;end
;;--------------------------------------------------------------------------------------------------------
pro six_aia_plot
cd,'D:\Learning\PHD1st\magnetic_reconnecion\paper\fig\'
;这个程序是针对12s间隔的AIA数据，将其裁剪、旋转、画图、除曝光时间保存为sav
;--------------------------------------------------------------------------------------------------
;nx=3;-------小图片列数
;ny=2;-------小图片行数
;n=nx*ny;----小图片个数
;
;;--------------------------这一段里的数字都是以一为整数来理解
;xoff0=0.1;----------------图片左边界与EPS右边界的距离
;xoff1=0.05;----------------图片右边界与EPS右边界的距离
;yoff1=0.05;----------------图片上边界与EPS上边界的距离
;yoff0=0.4;----------------图片下边界与EPS上边界的距离
;x_inter=0.01;------------图片之间X方向上的间距
;y_inter=0.01;------------图片之间Y方向上的间距
;
;width_eps=19.0;-----------------------------EPS的实际宽度(cm)
;xs=(1-xoff1-xoff0-(nx-1)*x_inter)/nx;-------每张图的宽度（总长为一）
;ys=(1-yoff1-yoff0-(ny-1)*y_inter)/ny;-------每张图的高度（总长为一）
;yx_ratio=(yuc+ydc)/(xrc+xlc);--------------------------------Y方向与X方向长度的缩放比例
;height_eps=width_eps*xs*yx_ratio/ys;--------EPS的高度!!!注意这里是保证最后图形长宽比例与实际相符的重要条件。
;set_plot,'ps'
;device,filename='fig5.eps',xsize=width_eps,ysize=height_eps,xoff=(21.0-width_eps)/2.0,$
;  yoff=(29.7-height_eps)/2.0,/color,bits=8,/times,/bold,isolatin1=1
;-------------------------------------------------------------------------------------------------------------
name_range=['.94','131','171','193','211','335']
value_range=[94,131,171,193,211,335]
map_list=list()
for i=0,5 do begin
  dir='D:/Learning/PHD1st/magnetic_reconnecion/data/AIA2/'+name_range[i]+'/'
  file_list=FILE_SEARCH(dir+'*.fits')
  j=88
  file=file_list[j]
  read_sdo,file,hdr,data
  hmi_prep,hdr,data,hdr0,data0
  index2map,hdr0,data0,map0
  map0.data=map0.data/map0.dur
  sub_map,map0,smap,xrange=[180,410],yrange=[-290,-520]
  datas=smap.data
  map_list.add,smap
endfor
; 设置文件路径和读取数据

; 设置为矢量图输出（EPS）
set_plot, 'ps'
device, filename='six_aia_plot.eps', /color, xsize=18,ysize=12,/encapsulated

!P.BACKGROUND = 255  ; 白底黑字
!P.COLOR = 0
!P.MULTI=[0,3,2]
for i=0,5 do begin
  aia_lct,wavelnth=value_range[i],/load
  plot_map,map_list[i],/log
endfor

device, /close    ; 关闭设备输出
set_plot, 'win'     ; 还原为默认设备（屏幕）
!P.MULTI=0
end

;;--------------------------------------------------------------------------------------------------------
;pro six_aia_plot
;  ; 这个程序是画六个AIA子图的
;  dir = 'D:/Learning/PHD1st/magnetic_reconnecion/data/AIA2/131/'
;  file_list = file_search(dir + '*.fits')
;  j = 88
;  file = file_list[j]
;  read_sdo, file, hdr, data
;  hmi_prep, hdr, data, hdr0, data0
;  index2map, hdr0, data0, map0
;  map0.data = map0.data / map0.dur
;  !P.BACKGROUND = 255
;  !P.COLOR = 0
;  sub_map, map0, smap, xrange=[180,410], yrange=[-290,-520]
;  aia_lct, wavelnth='131', /load
;
;  ; --- 关键修改1：先显示图形检查是否正确 ---
;  ; 先正常显示检查图形
;  window, /free, xsize=600, ysize=600
;  plot_map, smap, /log
;  ; 等待用户确认图形正确（按任意键继续）
;  print, '检查图形是否正确，按任意键继续生成PDF...'
;  void = dialog_message('检查图形是否正确，点击OK生成PDF', /info)
;
;  ; --- 关键修改2：PDF生成设置 ---
;  output_filename = 'aia_plot.pdf'
;  ; 关闭之前的图形窗口
;  window, /free  ; 清除当前窗口
;  ; 设置PS设备
;  set_plot, 'ps'
;  ; 必须设置字体，否则PDF可能无法打开
;  device, font_size=12, font='Helvetica'
;  ; 输出设置
;  device, filename=output_filename, $
;    /color, $
;    /encapsulated, $
;    xsize=8, ysize=8, /inches, $
;    bits_per_pixel=8  ; 重要参数
;
;  ; --- 重新绘制到PDF ---
;  plot_map, smap, /log
;
;  ; --- 关键修改3：正确关闭设备 ---
;  device, /close_file  ; 确保文件完全写入
;  set_plot, 'x'       ; 必须恢复默认设备
;
;  print, 'PDF已生成: ', output_filename
;end
;;------------------------------------------------------------------------------------------------------------------------------------
;  name_range=['.94','131','171','193','211','304','335']
;  value_range=[94,131,171,193,211,304,335]
;  target_time='
;  for i=0,6 do begin
;    dir='D:/Learning/PHD1st/magnetic_reconnecion/data/AIA2/'+name_range[i]+'/'
;    file_list=FILE_SEARCH(dir+'*.fits')
;    ;matrix_struct=replicate({matrix:ptr_new()},101)
;    map_list=list()
;    for j=133,147 do begin
;      file=file_list[j]
;      read_sdo,file,hdr,data
;      hmi_prep,hdr,data,hdr0,data0
;      index2map,hdr0,data0,map0
;      map0.data=map0.data/map0.dur
;      dmap=drot_map(map0,time='2024-06-18T21:17:00')
;      dmap.rtime=map0.time
;      sub_map,dmap,smap,xrange=[180,410],yrange=[-290,-520]
;      datas=smap.data
;      ;matrix_struct[j-32].matrix=ptr_new(datas)
;      map_list.add,smap
;      aia_lct,wavelnth=value_range[i],/load
;      WINDOW,/free,xsize=600,ysize=600
;      plot_map,smap,/log
;      write_bmp,'D:/Learning/PHD1st/magnetic_reconnecion/data/AIA2/'+name_range[i]+'fig/'+strmid(file_list[j],71,17)+'.bmp',tvrd(/true),/rgb
;      print,j
;      wdelete
;    endfor
;    ;    save,matrix_struct,filename='D:/Learning/PHD1st/magnetic_reconnecion/data/AIA2/pre_'+name_range[i]+'.sav'
;    ;    save,map_list,filename='D:/Learning/PHD1st/magnetic_reconnecion/data/AIA2/pre_map_'+name_range[i]+'.sav'
;    ;因为一些原因，这里加了pre_前缀的文件是0-31，不加的是32-132对应的序号
;  endfor