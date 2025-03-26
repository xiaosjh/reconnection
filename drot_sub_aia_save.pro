pro drot_sub_aia_save
  ;这个程序是针对12s间隔的AIA数据，将其裁剪、旋转、画图、除曝光时间保存为sav
  name_range=['.94','131','171','193','211','304','335']
  value_range=[94,131,171,193,211,304,335]
  target_time='
  for i=0,6 do begin
    dir='D:/Learning/PHD1st/magnetic_reconnecion/data/AIA2/'+name_range[i]+'/'
    file_list=FILE_SEARCH(dir+'*.fits')
    matrix_struct=replicate({matrix:ptr_new()},101)
    map_list=list()
    for j=0,31 do begin
      file=file_list[j]
      read_sdo,file,hdr,data
      hmi_prep,hdr,data,hdr0,data0
      index2map,hdr0,data0,map0
      map0.data=map0.data/map0.dur
      dmap=drot_map(map0,time='2024-06-18T21:17:00')
      dmap.rtime=map0.time
      sub_map,dmap,smap,xrange=[180,410],yrange=[-290,-520]
      datas=smap.data
      matrix_struct[j-32].matrix=ptr_new(datas)
      map_list.add,smap
      aia_lct,wavelnth=value_range[i],/load
      WINDOW,/free,xsize=600,ysize=600
      plot_map,smap,/log
      write_bmp,'D:/Learning/PHD1st/magnetic_reconnecion/data/AIA2/'+name_range[i]+'fig/'+strmid(file_list[j],71,17)+'.bmp',tvrd(/true),/rgb
      print,j
      wdelete
    endfor
    save,matrix_struct,filename='D:/Learning/PHD1st/magnetic_reconnecion/data/AIA2/pre_'+name_range[i]+'.sav'
    save,map_list,filename='D:/Learning/PHD1st/magnetic_reconnecion/data/AIA2/pre_map_'+name_range[i]+'.sav'
    ;因为一些原因，这里加了pre_前缀的文件是0-31，不加的是32-132对应的序号
  endfor
;----------------------------------------------------------------------------------------
;下面是处理1600，1700的
;  dir='D:/Learning/PHD1st/magnetic_reconnecion/data/AIA2/1600/'
;  file_list=file_search(dir+'*.fits')
;  matrix_struct=replicate({matrix:ptr_new()},52)
;  map_list=list()
;  for i=16,67 do begin
;    fits2map,file_list[i],map0
;    map0.data=map0.data/map0.dur
;    dmap=drot_map(map0,time='2024-06-18T21:17:00')
;    dmap.rtime=map0.time
;    sub_map,dmap,smap,xrange=[180,410],yrange=[-290,-520]
;    datas=smap.data
;    matrix_struct[i-16].matrix=ptr_new(datas)
;    map_list.add,smap
;    aia_lct,wavelnth=1600,/load
;    WINDOW,/free,xsize=600,ysize=600
;    plot_map,smap,/log
;    write_bmp,'D:/Learning/PHD1st/magnetic_reconnecion/data/AIA2/1600fig/'+strmid(file_list[i],71,17)+'.bmp',tvrd(/true),/rgb
;    print,i
;    wdelete
;  endfor
;;  save,matrix_struct,filename='D:/Learning/PHD1st/magnetic_reconnecion/data/AIA2/1600.sav'
;;  save,map_list,filename='D:/Learning/PHD1st/magnetic_reconnecion/data/AIA2/map_1600.sav'
;  
;  
;  dir='D:/Learning/PHD1st/magnetic_reconnecion/data/AIA2/1700/'
;  file_list=file_search(dir+'*.fits')
;  matrix_struct=replicate({matrix:ptr_new()},52)
;  map_list=list()
;  for i=16,67 do begin
;    fits2map,file_list[i],map0
;    map0.data=map0.data/map0.dur
;    dmap=drot_map(map0,time='2024-06-18T21:17:00')
;    dmap.rtime=map0.time
;    sub_map,dmap,smap,xrange=[180,410],yrange=[-290,-520]
;    datas=smap.data
;    matrix_struct[i-16].matrix=ptr_new(datas)
;    map_list.add,smap
;    aia_lct,wavelnth=1700,/load
;    WINDOW,/free,xsize=600,ysize=600
;    plot_map,smap,/log
;    write_bmp,'D:/Learning/PHD1st/magnetic_reconnecion/data/AIA2/1700fig/'+strmid(file_list[i],71,17)+'.bmp',tvrd(/true),/rgb
;    wdelete
;    print,i
;  endfor
;  save,matrix_struct,filename='D:/Learning/PHD1st/magnetic_reconnecion/data/AIA2/1700.sav'
;  save,map_list,filename='D:/Learning/PHD1st/magnetic_reconnecion/data/AIA2/map_1700.sav'
end