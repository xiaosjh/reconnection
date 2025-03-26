;flag=1
;
;if flag eq 1 then begin
; dir1='/Volumes/BACKUP2/scar_2nd/null_point/obs/aia_rain/rain_1/'
; dir2='/Volumes/BACKUP2/scar_2nd/null_point/obs/aia_rain/rain1_dem/'
;endif
;
;if flag eq 2 then begin
; dir1='/Volumes/BACKUP2/scar_2nd/null_point/obs/aia_rain/rain_2/'
; dir2='/Volumes/BACKUP2/scar_2nd/null_point/obs/aia_rain/rain2_dem/'
;endif
dir1='D:\Learning\PHD1st\magnetic_reconnecion\data\AIA2\'
file94  = findfile(dir1+'.94\'+'*.94.image_lev1.fits')
file131 = findfile(dir1+'131\'+'*.131.image_lev1.fits')
file171 = findfile(dir1+'171\'+'*.171.image_lev1.fits')
file193 = findfile(dir1+'193\'+'*.193.image_lev1.fits')
file211 = findfile(dir1+'211\'+'*.211.image_lev1.fits')
file335 = findfile(dir1+'335\'+'*.335.image_lev1.fits')
ntmp=size(file131)
for k=133,147 do begin
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

aia_dem_setup_x,tresp
max_t=7.5
min_t=5.5
dt=0.1
mc_iter=0
result={dem_out:dblarr(xs,ys,(max_t-min_t)/dt+1),$
chisq_out:dblarr(xs,ys),$
smap_temp:fltarr(xs,ys),smap_em:fltarr(xs,ys)}

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
obs_err=aia_bp_estimate_error(obs_val,[131,94,211,193,171,335])
xrt_dem_iterative2,filterset,obs_val,tresp,logt_out,dem_out,$
base_obs=base_obs,mod_obs=mod_obs,chisq=chisq,obs_err=obs_err,max_t=max_t,min_t=min_t,dt=dt,mc_iter=mc_iter

result.dem_out[i,j,*]=dem_out
result.chisq_out[i,j]=chisq
;;;calculate the temperature and em maps
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
dir2='D:\Learning\PHD1st\magnetic_reconnecion\data\AIA2\DEM2\'
save,result,smap1,filename=dir2+'DEM'+title+'.sav'

print,string(k)+'has been done!!!!!!'

endfor

end
