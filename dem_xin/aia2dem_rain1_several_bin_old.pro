;这个是我分成了1MK,1.5MK,2MK,4MK,8MK来保存平均温度和EM
dir1='C:\Learning\PHD1st\magnetic_reconnecion\data\AIA2\'
file94  = findfile(dir1+'.94\'+'*.94.image_lev1.fits')
file131 = findfile(dir1+'131\'+'*.131.image_lev1.fits')
file171 = findfile(dir1+'171\'+'*.171.image_lev1.fits')
file193 = findfile(dir1+'193\'+'*.193.image_lev1.fits')
file211 = findfile(dir1+'211\'+'*.211.image_lev1.fits')
file335 = findfile(dir1+'335\'+'*.335.image_lev1.fits')
ntmp=size(file131)
for k=32,132 do begin
;这个算的是从头到尾一共148个AIA文件
;for k=175,195,4 do begin
;fileset=[file131[k],file94[k],file211[k],file193[k],file171[k-1],file335[k]]
filterset=['A131','A94','A211','A193','A171','A335']

reftime = '2024-06-18T21:17:00'
xrange = [180,410]
yrange = [-290,-520]

;;;read AIA data
read_sdo,file131[k],index,data
hmi_prep,index,data,index1,data1
index2map,index1,data1,map1
map1=drot_map(map1,time=reftime)
sub_map,map1,smap1,xrange=xrange,yrange=yrange

read_sdo,file94[k],index,data
hmi_prep,index,data,index2,data2
index2map,index2,data2,map2
map2=drot_map(map2,time=reftime)
sub_map,map2,smap2,xrange=xrange,yrange=yrange

read_sdo,file211[k],index,data
hmi_prep,index,data,index3,data3
index2map,index3,data3,map3
map3=drot_map(map3,time=reftime)
sub_map,map3,smap3,xrange=xrange,yrange=yrange

read_sdo,file193[k],index,data
hmi_prep,index,data,index4,data4
index2map,index4,data4,map4
map4=drot_map(map4,time=reftime)
sub_map,map4,smap4,xrange=xrange,yrange=yrange

read_sdo,file171[k],index,data
hmi_prep,index,data,index5,data5
index2map,index5,data5,map5
map5=drot_map(map5,time=reftime)
sub_map,map5,smap5,xrange=xrange,yrange=yrange

read_sdo,file335[k],index,data
hmi_prep,index,data,index6,data6
index2map,index6,data6,map6
map6=drot_map(map6,time=reftime)
sub_map,map6,smap6,xrange=xrange,yrange=yrange

;;return the data you want to calculate
nn=size(smap1.data)
xs=fix(nn[1]/2) ;;;rebin the data
ys=fix(nn[2]/2)
smap10=rebin_map(smap1,xs,ys)
smap20=rebin_map(smap2,xs,ys)
smap30=rebin_map(smap3,xs,ys)
smap40=rebin_map(smap4,xs,ys)
smap50=rebin_map(smap5,xs,ys)
smap60=rebin_map(smap6,xs,ys)
   ;window,1,xs=xs,ys=ys
   ;plot_image,alog10(smap[0].data)
   ;stop

;;define the maps of temperature and EM
smap_temp=fltarr(xs,ys)
smap_em=fltarr(xs,ys)
;;;setup DEM
print, 'Setting up AIA DEM response...'
aia_dem_setup_x, tresp

; 立即检查
if ~isa(tresp, 'STRUCT') then begin
  print, 'FATAL ERROR: AIA_DEM_SETUP_X failed to return a structure'
  print, 'Check your SSW installation and AIA response files'
  stop
endif

print, 'DEM response setup successful'
help, tresp  ; 打印详细信息

max_t=7.5
min_t=5.5
dt=0.1
mc_iter=0
result={dem_out:dblarr(xs,ys,(max_t-min_t)/dt+1),$
chisq_out:dblarr(xs,ys),$
smap_temp:fltarr(xs,ys),smap_em:fltarr(xs,ys),$
smap_temp_1:fltarr(xs,ys),smap_em_1:fltarr(xs,ys),$
smap_temp_16:fltarr(xs,ys),smap_em_16:fltarr(xs,ys),$
smap_temp_2:fltarr(xs,ys),smap_em_2:fltarr(xs,ys),$
smap_temp_4:fltarr(xs,ys),smap_em_4:fltarr(xs,ys),$
smap_temp_8:fltarr(xs,ys),smap_em_8:fltarr(xs,ys),temp_all:fltarr(xs,ys,(max_t-min_t)/dt+1)}

for i=0,xs-1 do begin
   for j=0,ys-1 do begin
;;;input
obs_val=fltarr(6)
obs_val[0]=(smap10.data)[i,j]/smap1.dur
obs_val[1]=(smap20.data)[i,j]/smap2.dur
obs_val[2]=(smap30.data)[i,j]/smap3.dur
obs_val[3]=(smap40.data)[i,j]/smap4.dur
obs_val[4]=(smap50.data)[i,j]/smap5.dur
obs_val[5]=(smap60.data)[i,j]/smap6.dur

;obs_err=aia_bp_estimate_error(obs_val,[131,94,211,193,171,335])
obs_err = fltarr(6)
obs_err[0] = aia_bp_estimate_error(obs_val[0], 131, /evenorm)
obs_err[1] = aia_bp_estimate_error(obs_val[1], 94,  /evenorm)
obs_err[2] = aia_bp_estimate_error(obs_val[2], 211, /evenorm)
obs_err[3] = aia_bp_estimate_error(obs_val[3], 193, /evenorm)
obs_err[4] = aia_bp_estimate_error(obs_val[4], 171, /evenorm)
obs_err[5] = aia_bp_estimate_error(obs_val[5], 335, /evenorm)

xrt_dem_iterative2,filterset,obs_val,tresp,logt_out,dem_out,$
base_obs=base_obs,mod_obs=mod_obs,chisq=chisq,obs_err=obs_err,max_t=max_t,min_t=min_t,dt=dt,mc_iter=mc_iter

result.dem_out[i,j,*]=dem_out
result.chisq_out[i,j]=chisq
;;;calculate the temperature and em maps

temp=10^logt_out
result.temp_all[i,j,*]=temp
dem=dem_out[*,0]

min_tt=4
max_tt=6
em_1=int_tabulated(temp[min_tt:max_tt],dem[min_tt:max_tt])
temp_1=int_tabulated(temp[min_tt:max_tt],(temp[min_tt:max_tt]*dem[min_tt:max_tt]))
temp_1=temp_1/em_1
result.smap_temp_1[i,j]=temp_1
result.smap_em_1[i,j]=em_1

min_tt=6
max_tt=8
em_16=int_tabulated(temp[min_tt:max_tt],dem[min_tt:max_tt])
temp_16=int_tabulated(temp[min_tt:max_tt],(temp[min_tt:max_tt]*dem[min_tt:max_tt]))
temp_16=temp_16/em_16
result.smap_temp_16[i,j]=temp_16
result.smap_em_16[i,j]=em_16

min_tt=7
max_tt=9
em_2=int_tabulated(temp[min_tt:max_tt],dem[min_tt:max_tt])
temp_2=int_tabulated(temp[min_tt:max_tt],(temp[min_tt:max_tt]*dem[min_tt:max_tt]))
temp_2=temp_2/em_2
result.smap_temp_2[i,j]=temp_2
result.smap_em_2[i,j]=em_2

min_tt=10
max_tt=12
em_4=int_tabulated(temp[min_tt:max_tt],dem[min_tt:max_tt])
temp_4=int_tabulated(temp[min_tt:max_tt],(temp[min_tt:max_tt]*dem[min_tt:max_tt]))
temp_4=temp_4/em_4
result.smap_temp_4[i,j]=temp_4
result.smap_em_4[i,j]=em_4

min_tt=13
max_tt=15
em_8=int_tabulated(temp[min_tt:max_tt],dem[min_tt:max_tt])
temp_8=int_tabulated(temp[min_tt:max_tt],(temp[min_tt:max_tt]*dem[min_tt:max_tt]))
temp_8=temp_8/em_8
result.smap_temp_8[i,j]=temp_8
result.smap_em_8[i,j]=em_8

min_tt=2
max_tt=18
temp=10^logt_out
dem=dem_out[*,0]
em=int_tabulated(temp[min_tt:max_tt],dem[min_tt:max_tt])
temp=int_tabulated(temp[min_tt:max_tt],(temp[min_tt:max_tt]*dem[min_tt:max_tt]))
temp=temp/em
result.smap_temp[i,j]=temp
result.smap_em[i,j]=em
   endfor
print,i,k
endfor

title=strmid(smap1[0].time,0,11)+'_'+strmid(smap1[0].time,12,2)+'_'+strmid(smap1[0].time,15,2)+'_'+strmid(smap1[0].time,18,2)
dir2='C:\Learning\PHD1st\magnetic_reconnecion\data\AIA2\DEM_bin\'
save,result,smap1,filename=dir2+'DEM'+title+'.sav'

print,string(k)+'has been done!!!!!!'

endfor

end
