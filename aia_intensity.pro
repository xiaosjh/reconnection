pro aia_intensity;本程序用于画AIA这几个波段强度随时间的变化
;aia_int=fltarr(7,20)
;file_list=file_search('D:\Learning\PHD1st\magnetic_reconnecion\data\rotAIA\'+'*.sav')
;for j=0,6 do begin
;  restore,file_list[j]
;  for i=0,19 do begin
;    map0=maplist[i]
;    map0.data=map0.data/map0.dur
;    aia_int[j,i]=total(map0.data)
;  endfor
;endfor
;save,filename='D:\Learning\PHD1st\magnetic_reconnecion\data_process\aia_int.sav',aia_int
;---------------------------------------------------------------------------------------------
;restore,'D:\Learning\PHD1st\magnetic_reconnecion\data\rotAIA\all_131_map.sav'
;for j=0,19 do begin
;  WINDOW,/FREE,XSIZE=400,YSIZE=400
;  map0=maplist[j]
;  map0.data=map0.data/map0.dur
;  sub_map,map0,smap,xrange=[200,420],yrange=[-300,-520]
;  aia_lct,wavelnth=131,/load
;  plot_map,smap,/log
;  write_bmp,'D:\Learning\PHD1st\AIA\shortAIA131\'+string(j)+'.bmp',tvrd(/true),/rgb
;endfor
;-------------------------------------------------------------------------------------------
;file_list=file_search('D:\Learning\PHD1st\magnetic_reconnecion\data\AIA\AIA21_17\'+'*.fits')
;WINDOW, /FREE, XSIZE = 1000, YSIZE = 600
;!P.MULTI=[0,3,2]
;xrange=[1,4,5,3,0,2]
;colors=[94,131,171,193,211,335]
;for i=0,5 do begin
;  fits2map,file_list[xrange[i]],map0
;  aia_lct,wavelnth=colors[i],/load
;  map0.data=map0.data/map0.dur
;  sub_map,map0,smap,xrange=[200,420],yrange=[-300,-520]
;  plot_map,smap,/log
;endfor
;!P.MULTI=0
;write_bmp,'D:\Learning\PHD1st\magnetic_reconnecion\data\AIA\AIA21_17\all_aia.bmp',tvrd(/true),/rgb
;--------------------------------------------------------------------------------------------------
file1='D:\Learning\PHD1st\magnetic_reconnecion\data\HMIsharp\hmi.sharp_cea_720s_nrt_1379026_1379026\SUM23\D1758488751\S00000\Bp.fits'
fits2map,file1,map0
Bp=congrid(map0.data,192,88)
file1='D:\Learning\PHD1st\magnetic_reconnecion\data\HMIsharp\hmi.sharp_cea_720s_nrt_1379026_1379026\SUM23\D1758488751\S00000\Br.fits'
fits2map,file1,map0
Br=congrid(map0.data,192,88)
file1='D:\Learning\PHD1st\magnetic_reconnecion\data\HMIsharp\hmi.sharp_cea_720s_nrt_1379026_1379026\SUM23\D1758488751\S00000\Bt.fits'
fits2map,file1,map0
Bt=congrid(map0.data,192,88)
save,filename='D:\Learning\PHD1st\magnetic_reconnecion\data\HMIsharp\hmi.sharp_cea_720s_nrt_1379026_1379026\SUM23\D1758488751\S00000\congrid_rpt.sav',Bp,Br,Bt
end