; .r aia2dem_1
pro aia2dem_1

dir1='/media/fy/inpre_paper/2023_12_31/AIA+DEM/'
dir2='~/Desktop/dem/dem/5.5-7.5/'


file94  = findfile(dir1+'*.94.image_lev1.fits')
file131 = findfile(dir1+'*.131.image_lev1.fits')
file171 = findfile(dir1+'*.171.image_lev1.fits')
file193 = findfile(dir1+'*.193.image_lev1.fits')
file211 = findfile(dir1+'*.211.image_lev1.fits')
file335 = findfile(dir1+'*.335.image_lev1.fits')
ntmp=size(file131)

for k=59,84,2 do begin

; for k=71,88,2 do begin
;fileset=[file131[k],file94[k],file211[k],file193[k],file171[k-1],file335[k]]
filterset=['A131','A94','A211','A193','A171','A335']

centerx=0  & centery=0

xls1=1070.0  &  xrs1=-970.0   &  yds1=-35.0  &  yus1=115.0
xss1s0=centerx-xls1 & xss1s1=centerx+xrs1 & yss1s0=centery-yds1 & yss1s1=centery+yus1  


xrange=[xss1s0,xss1s1]  &  yrange=[yss1s0,yss1s1]



faia=[file131[k],file94[k+1],file211[k+1],file193[k],file171[k],file335[k+1]]
;;read AIA data
fits2map,faia[0],map1
sub_map,map1,smap1,xrange=xrange,yrange=yrange

fits2map,faia[1],map2
sub_map,map2,smap2,xrange=xrange,yrange=yrange

fits2map,faia[2],map3
sub_map,map3,smap3,xrange=xrange,yrange=yrange

fits2map,faia[3],map4
sub_map,map4,smap4,xrange=xrange,yrange=yrange

fits2map,faia[4],map5
sub_map,map5,smap5,xrange=xrange,yrange=yrange

fits2map,faia[5],map6
sub_map,map6,smap6,xrange=xrange,yrange=yrange

print,k & print,map1.dur & print,map2.dur & print,map3.dur & print,map4.dur & print,map5.dur & print,map6.dur



;;return the data you want to calculate
nn=size(smap1.data)
; xs=nn[1] ;;;rebin the data
; ys=nn[2]
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
min_tt=0
max_tt=20
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
save,result,smap10,smap20,smap30,smap40,smap50,smap60,filename=dir2+'DEM'+title+'.sav'

print,string(k)+'has been done!!!!!!'

endfor

end
