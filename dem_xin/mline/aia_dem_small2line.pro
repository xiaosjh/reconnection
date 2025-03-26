; .r aia_dem_small2line
; xrt_dem_iterative2
; aia_dem_setup_x

cd,'/media/fy/inpre_paper/dem_xin/mline'

dir1='/media/fy/inpre_paper/2023_12_31/AIA+DEM/'
dir2='/media/fy/inpre_paper/dem_xin/mline/'

file94  = findfile(dir1+'*.94.image_lev1.fits')
file131 = findfile(dir1+'*.131.image_lev1.fits')
file171 = findfile(dir1+'*.171.image_lev1.fits')
file193 = findfile(dir1+'*.193.image_lev1.fits')
file211 = findfile(dir1+'*.211.image_lev1.fits')
file335 = findfile(dir1+'*.335.image_lev1.fits')
ntmp=size(file131)

k=69
; k=73
; k=75
; k=79
; 
;fileset=[file131[k],file94[k],file211[k],file193[k],file171[k-1],file335[k]]
filterset=['A131','A94','A211','A193','A171','A335']

centerx=0  & centery=0

xdems1=[1035.0+3,1035.0-8,1035.0-9,1035.0-4]
xdems2=[-1030.0-3,-1030.0+8,-1030.0+9,-1030.0+4]
ydems1=[-44.0,-44.0,-44.0,-44.0]-4.5
ydems2=[48.0,48.0,48.0,48.0]+4.5
xdem0=centerx-xdems1 & xdem1=centerx+xdems2 & ydem0=centery-ydems1 & ydem1=centery+ydems2  


xdemr=centerx-xdems1 & xdeml=centerx+xdems2 & ydemd=centery-ydems1 & ydemu=centery+ydems2  

xrange=[xdemr[0],xdeml[0]]  &  yrange=[ydemd[0],ydemu[0]]
; xrange=[xdemr[1],xdeml[1]]  &  yrange=[ydemd[1],ydemu[1]]
; xrange=[xdemr[2],xdeml[2]]  &  yrange=[ydemd[2],ydemu[2]]
; xrange=[xdemr[3],xdeml[3]]  &  yrange=[ydemd[3],ydemu[3]]



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
xs=nn[1] ;;;rebin the data
ys=nn[2]
; xs=fix(nn[1]/2) ;;;rebin the data
; ys=fix(nn[2]/2)
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
save,result,smap10,smap20,smap30,smap40,smap50,smap60,filename=dir2+'sDEM'+title+'.sav'

print,string(k)+'has been done!!!!!!'


end
