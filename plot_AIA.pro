;这个程序用于画AIA的图像
PRO plot_AIA
  dir = 'D:/Learning/PHD1st/magnetic_reconnecion/data/AIA/AIA.94/'
  file_list = FILE_SEARCH(dir + '*.fits')
  file=file_list[0]
  fits2map,file,map0
  map0.data=map0.data/map0.dur
  sub_map,map0,smap,xrange=[100,600],yrange=[-200,-700]
  aia_lct,wavelnth=94,/load
  plot_map,smap,/log
  ;write_bmp,'D:/111.bmp',tvrd(/true),/rgb
  write_bmp,'D:/Learning/PHD1st/AIA/AIA.94/'+strmid(file_list[0],56,50)+'.bmp',tvrd(/true),/rgb
  for i=0,(size(file_list))[1]-2 DO BEGIN
    file2=file_list[i+1]
    fits2map,file2,map2
    map2.data=map2.data/map2.dur
    dmap=drot_map(map2,ref_map=smap)
    dmap.RTIME=map2.TIME
    ;sub_map,dmap,smap2,xrange=[100,600],yrange=[-200,-700]

    aia_lct,wavelnth=94,/load
    plot_map,dmap,/log
    write_bmp,'D:/Learning/PHD1st/AIA/AIA.94/'+strmid(file_list[i+1],56,50)+'.bmp',tvrd(/true),/rgb
  endfor
  dir = 'D:/Learning/PHD1st/magnetic_reconnecion/data/AIA/AIA131/'
  file_list = FILE_SEARCH(dir + '*.fits')
  file=file_list[0]
  fits2map,file,map0
  map0.data=map0.data/map0.dur
  sub_map,map0,smap,xrange=[100,600],yrange=[-200,-700]
  aia_lct,wavelnth=131,/load
  plot_map,smap,/log
  ;write_bmp,'D:/111.bmp',tvrd(/true),/rgb
  write_bmp,'D:/Learning/PHD1st/AIA/AIA131/'+strmid(file_list[0],56,50)+'.bmp',tvrd(/true),/rgb
  for i=0,1 DO BEGIN
    file2=file_list[i+1]
    fits2map,file2,map2
    map2.data=map2.data/map2.dur  
    dmap=drot_map(map2,ref_map=smap)
    dmap.RTIME=map2.TIME
    ;sub_map,dmap,smap2,xrange=[100,600],yrange=[-200,-700]
    aia_lct,wavelnth=131,/load
    plot_map,dmap,/log
    write_bmp,'D:/Learning/PHD1st/AIA/AIA131/'+strmid(file_list[i+1],56,50)+'.bmp',tvrd(/true),/rgb
  endfor
  ;
  dir = 'D:/Learning/PHD1st/magnetic_reconnecion/data/AIA/AIA171/'
  file_list = FILE_SEARCH(dir + '*.fits')
  file=file_list[0]
  fits2map,file,map0
  map0.data=map0.data/map0.dur
  sub_map,map0,smap,xrange=[-100,700],yrange=[0,-600]
  aia_lct,wavelnth=171,/load
  plot_map,smap,/log
  stop
  ;write_bmp,'D:/111.bmp',tvrd(/true),/rgb
  write_bmp,'D:/Learning/PHD1st/AIA/AIA171/'+strmid(file_list[0],56,50)+'.bmp',tvrd(/true),/rgb
  for i=0,(size(file_list))[1]-2 DO BEGIN
    file2=file_list[i+1]
    fits2map,file2,map2
    map2.data=map2.data/map2.dur
    dmap=drot_map(map2,ref_map=smap)
    dmap.RTIME=map2.TIME
    ;sub_map,dmap,smap2,xrange=[100,600],yrange=[-200,-700]
  
    aia_lct,wavelnth=171,/load
    plot_map,dmap,/log
    write_bmp,'D:/Learning/PHD1st/AIA/AIA171/'+strmid(file_list[i+1],56,50)+'.bmp',tvrd(/true),/rgb
  endfor
  
  dir = 'D:/Learning/PHD1st/magnetic_reconnecion/data/AIA/AIA193/'
  file_list = FILE_SEARCH(dir + '*.fits')
  file=file_list[0]
  fits2map,file,map0
  map0.data=map0.data/map0.dur
  sub_map,map0,smap,xrange=[100,600],yrange=[-200,-700]
  aia_lct,wavelnth=193,/load
  plot_map,smap,/log
  ;write_bmp,'D:/111.bmp',tvrd(/true),/rgb
  write_bmp,'D:/Learning/PHD1st/AIA/AIA193/'+strmid(file_list[0],56,50)+'.bmp',tvrd(/true),/rgb
  for i=0,(size(file_list))[1]-2 DO BEGIN
    file2=file_list[i+1]
    fits2map,file2,map2
    map2.data=map2.data/map2.dur
    dmap=drot_map(map2,ref_map=smap)
    dmap.RTIME=map2.TIME
    sub_map,dmap,smap2,xrange=[100,600],yrange=[-200,-700]
  
    aia_lct,wavelnth=193,/load
    plot_map,dmap,/log
    write_bmp,'D:/Learning/PHD1st/AIA/AIA193/'+strmid(file_list[i+1],56,50)+'.bmp',tvrd(/true),/rgb
  endfor
  
  dir = 'D:/Learning/PHD1st/magnetic_reconnecion/data/AIA/AIA211/'
  file_list = FILE_SEARCH(dir + '*.fits')
  file=file_list[0]
  fits2map,file,map0
  map0.data=map0.data/map0.dur
  sub_map,map0,smap,xrange=[100,600],yrange=[-200,-700]
  aia_lct,wavelnth=211,/load
  plot_map,smap,/log
  ;write_bmp,'D:/111.bmp',tvrd(/true),/rgb
  write_bmp,'D:/Learning/PHD1st/AIA/AIA211/'+strmid(file_list[0],56,50)+'.bmp',tvrd(/true),/rgb
  for i=0,(size(file_list))[1]-2 DO BEGIN
    file2=file_list[i+1]
    fits2map,file2,map2
    map2.data=map2.data/map2.dur
    dmap=drot_map(map2,ref_map=smap)
    dmap.RTIME=map2.TIME
    ;sub_map,dmap,smap2,xrange=[100,600],yrange=[-200,-700]
  
    aia_lct,wavelnth=211,/load
    plot_map,dmap,/log
    write_bmp,'D:/Learning/PHD1st/AIA/AIA211/'+strmid(file_list[i+1],56,50)+'.bmp',tvrd(/true),/rgb
  endfor
  
  dir = 'D:/Learning/PHD1st/magnetic_reconnecion/data/AIA/AIA304/'
  file_list = FILE_SEARCH(dir + '*.fits')
  file=file_list[0]
  fits2map,file,map0
  map0.data=map0.data/map0.dur
  sub_map,map0,smap,xrange=[100,600],yrange=[-200,-700]
  aia_lct,wavelnth=304,/load
  plot_map,smap,/log
  ;write_bmp,'D:/111.bmp',tvrd(/true),/rgb
  write_bmp,'D:/Learning/PHD1st/AIA/AIA304/'+strmid(file_list[0],56,50)+'.bmp',tvrd(/true),/rgb
  for i=0,(size(file_list))[1]-2 DO BEGIN
    file2=file_list[i+1]
    fits2map,file2,map2
    map2.data=map2.data/map2.dur
    dmap=drot_map(map2,ref_map=smap)
    dmap.RTIME=map2.TIME
    ;sub_map,dmap,smap2,xrange=[100,600],yrange=[-200,-700]
  
    aia_lct,wavelnth=304,/load
    plot_map,dmap,/log
    write_bmp,'D:/Learning/PHD1st/AIA/AIA304/'+strmid(file_list[i+1],56,50)+'.bmp',tvrd(/true),/rgb
  endfor
  
  dir = 'D:/Learning/PHD1st/magnetic_reconnecion/data/AIA/AIA335/'
  file_list = FILE_SEARCH(dir + '*.fits')
  file=file_list[0]
  fits2map,file,map0
  map0.data=map0.data/map0.dur
  sub_map,map0,smap,xrange=[100,600],yrange=[-200,-700]
  aia_lct,wavelnth=335,/load
  plot_map,smap,/log
  ;write_bmp,'D:/111.bmp',tvrd(/true),/rgb
  write_bmp,'D:/Learning/PHD1st/AIA/AIA335/'+strmid(file_list[0],56,50)+'.bmp',tvrd(/true),/rgb
  for i=0,(size(file_list))[1]-2 DO BEGIN
    file2=file_list[i+1]
    fits2map,file2,map2
    map2.data=map2.data/map2.dur
    dmap=drot_map(map2,ref_map=smap)
    dmap.RTIME=map2.TIME
    ;sub_map,dmap,smap2,xrange=[100,600],yrange=[-200,-700]
  
    aia_lct,wavelnth=335,/load
    plot_map,dmap,/log
    write_bmp,'D:/Learning/PHD1st/AIA/AIA335/'+strmid(file_list[i+1],56,50)+'.bmp',tvrd(/true),/rgb
  endfor
end