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
  
end