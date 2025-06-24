pro drot_sub_hmi_save
dir='D:\Learning\PHD1st\magnetic_reconnecion\data\HMI3\'
file_list=FILE_SEARCH(dir+'*.fits')
map_list=list()
;预处理
for i=0,53 do begin
  file=file_list[i]
  read_sdo,file,hdr,data
  hmi_prep,hdr,data,hdr0,data0
  index2map,hdr0,data0,map0

  ;旋转，裁剪
  ;map0.data=map0.data/map0.dur;hmi是不是没有exposure_time
  dmap=drot_map(map0,time='2024-06-18T21:18:00')
  img=readfits(file,head)
  dmap.rtime=map0.time
  sub_map,dmap,ssmap,xrange=[200,420],yrange=[-520,-300]
  map_list.add,ssmap
  WINDOW,/free,xsize=600,ysize=600
  plot_map,ssmap
  write_bmp,'D:\Learning\PHD1st\magnetic_reconnecion\data\HMIfig_idl\'+string(i,format='(I2.2)')+'.bmp',tvrd(/true),/rgb
  wdelete
endfor
save,map_list,filename='D:\Learning\PHD1st\magnetic_reconnecion\data_process\drot_sub_hmi.sav'


end