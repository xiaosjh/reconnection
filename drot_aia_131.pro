pro drot_aia_131 ;这个文件用于把20个时间的131图像统一旋转到 '2024-06-18T21:18:06.622'
dir='D:\Learning\PHD1st\magnetic_reconnecion\data\AIA\AIA131\'
file_list=file_search(dir+'*.fits')
fits2map,file_list[76],map0
for i=65,84 DO BEGIN
  if i ne 76 then begin
    file2=file_list[i]
    read_sdo,file2,hdr,data
    hmi_prep,hdr,data,hdr0,data0
    index2map,hdr0,data0,map2
    map2.data=map2.data/map2.dur
    dmap=drot_map(map2,time=map0.time)
    dmap.RTIME=map2.TIME
    map2fits,dmap,'D:\Learning\PHD1st\magnetic_reconnecion\data\AIA\AIA131rot\'+'rot_'+strmid(file_list[i],strlen(file_list[i])-55,55)
    print,i
  endif
    
endfor

end