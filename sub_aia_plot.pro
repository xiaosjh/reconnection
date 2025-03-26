pro sub_aia_plot
;这个文件用来画小一点区域的aia图像
name_range=['AIA.94/','AIA131/','AIA171/','AIA193/','AIA211/','AIA304/','AIA335/']
value_range=[94,131,171,193,211,304,335]
for j=0,6 do begin
  dir = 'D:/Learning/PHD1st/magnetic_reconnecion/data/AIA/'+name_range[j]
  file_list = FILE_SEARCH(dir + '*.fits')
  file=file_list[50]
  fits2map,file,map0
  map0.data=map0.data/map0.dur
  sub_map,map0,smap,xrange=[200,420],yrange=[-300,-520]
  aia_lct,wavelnth=value_range[j],/load
  plot_map,smap,/log
  write_bmp,'D:/Learning/PHD1st/AIA/'+name_range[j]+strmid(file_list[50],56,50)+'.bmp',tvrd(/true),/rgb
  
  for i=50,100 DO BEGIN
    file2=file_list[i+1]
    fits2map,file2,map2
    map2.data=map2.data/map2.dur
    dmap=drot_map(map2,ref_map=smap)
    dmap.RTIME=map2.TIME
    aia_lct,wavelnth=value_range[j],/load
    plot_map,dmap,/log
    write_bmp,'D:/Learning/PHD1st/AIA/'+name_range[j]+strmid(file_list[i+1],56,50)+'.bmp',tvrd(/true),/rgb
  endfor
endfor
end