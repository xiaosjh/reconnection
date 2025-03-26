pro rot_aia;提前把fits文件旋转一下，然后map保存到sav文件中
;这个是131的，然后我把别的也旋转一下吧
;dir1='D:\Learning\PHD1st\magnetic_reconnecion\data\AIA_DEM_test1\'
;file=findfile(dir1+'aia.lev1_euv_12s.*.131.image_lev1.fits')
;maplist=list()
;maplist.add,map0
;n00=0
;num=20
;for i=0,num-1,1 do begin
;  fits2map,file[n00+i],maps
;  if i ne 11 then begin
;    dmap=drot_map(maps,time='2024-06-18T21:18:06.622')
;    dmap.RTIME=maps.TIME
;    maplist.add,dmap
;    print,i
;  end else begin
;    maplist.add,maps
;    print,i
;  end
;endfor
;save,filename='D:\Learning\PHD1st\magnetic_reconnecion\data\all_131_map.sav',maplist
;end
name_list=['94','131','171','193','211','304','335']
for j=0,6 do begin
  str='aia.lev1_euv_12s.*.'+name_list[j]+'.image_lev1.fits'
  dir1='D:\Learning\PHD1st\magnetic_reconnecion\data\AIA_DEM_test1\'
  file=findfile(dir1+str)
  maplist=list()
  n00=0
  num=20
  print,'j=',j
  for i=0,num-1,1 do begin
    fits2map,file[n00+i],maps
    dmap=drot_map(maps,time='2024-06-18T21:18:06.000')
    dmap.RTIME=maps.TIME
    maplist.add,dmap
    print,i
  endfor
  save,filename='D:\Learning\PHD1st\magnetic_reconnecion\data\rotAIA\all_'+name_list[j]+'_map.sav',maplist
endfor

end