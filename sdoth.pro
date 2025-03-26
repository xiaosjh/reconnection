; .r sdoth

pro sdoth
cd,'D:\Learning\PHD1st\magnetic_reconnecion\program\'

dir1='D:\Learning\PHD1st\magnetic_reconnecion\data\AIA_DEM_test1\'
file=findfile(dir1+'aia.lev1_euv_12s.*.131.image_lev1.fits')
n=n_elements(file)
file=file(sort(file))
print,n

n00=0
n11=19
n304=n11-n00+1
help,n304
sum304=fltarr(n304)  &  t304=strarr(n304)

x1=200 & x2=420 & y1=-300 & y2=-520
  
ov=1
read_sdo,file[12],hdr0,data0
hmi_prep,hdr0,data0,hdr,data
index2map,hdr,data,map0
map0.data=map0.data/map0.dur
rmap0=get_sub_map(map0,xrange=[x1,x2],yrange=[y1,y2])
scale=size(rmap0.data,/dim) & print,scale
xs=scale(0) & ys=scale(1)
aia_lct,wave=131,/load
wdef,0,xs*ov,ys*ov
tv,bytscl(alog(rmap0.data),0.2,6.7)


sp=intarr(2,300)
        i=-1
        repeat begin
            cursor,x0,y0,/down,/device
            if !err eq 1 then print,x0,y0
            i=i+1 & sp[*,i]=[x0,y0]
            if (i gt 0) and (!err eq 1) then plots,[sp[0,i-1],sp[0,i]],[sp[1,i-1],sp[1,i]],color=0,/dev
        endrep until (!err eq 4)
print,'本次共取点：',i,'个。'
sp=sp(*,0:i) & help,sp
dpix=1 & loop_interpol,dpix,sp[0,*],sp[1,*],xx,yy 
help,xx,yy & size=size(xx,/dim) & nm=size[0]


num=n304
slice304=fltarr(num,nm) & time=strarr(num)
restore,'D:\Learning\PHD1st\magnetic_reconnecion\data\all_131_map.sav'
for i=1,num,1 do begin
;for i=0,num-1,1 do begin
    ;fits2map,file[n00+i],map
    ;map=drot_map(map,time=map0.time)
    maps=maplist[n00+i]
    rmap=get_sub_map(maps,xrange=[x1,x2],yrange=[y1,y2])
    rmap.data=rmap.data/maps.dur
    im=congrid(rmap.data,xs*ov,ys*ov,/interp)

    
    wdef,1,xs*ov,ys*ov
    tv,bytscl(alog(im),0.2,6.7)
    time1=rmap.time & time2=strmid(time1,0,11)+' '+strmid(time1,11,9)+' '+'UT'
    xyouts,0.5,0.02,'171 '+time2,charsize=1.5,charthick=1.5,/norm,font=0,alignment=0.5
    
    tvlct,0,255,0,252;green
    plots,xx,yy,color=252,thick=1.5,linestyle=0,/dev
    plots,xx-5,yy-5,color=252,thick=1.5,linestyle=0,/dev
    plots,xx+5,yy+5,color=252,thick=1.5,linestyle=0,/dev
    ;name='D:/Learning/PHD1st/AIA/AIA171_slice/'+'171'+time2+'.bmp'
    name='D:\Learning\PHD1st\AIA\AIA131slice\'+strtrim(i)+'.bmp'
    write_bmp,name,tvrd(/true),/rgb
    ;我这样子修改一下，可以让他变成求一块区域的均值
    ;slice304[i,*]=im[xx,yy]
;    slice_slice=fltarr(5,nm)
;    for i=-2,2,1 do begin
;      slice_slice[i+2,*]=im[xx+i,yy+i]
;    endfor
;    slice304[i,*]=total(slice_slice,1)/5
    ;slice304[i-1,*]=im[xx,yy]
    slice304[i-1,*]=(im[xx-5,yy-5]+im[xx,yy]+im[xx+5,yy+5])/3.; & help,slice193[i,*]
    
    

endfor    

slice304=slice304(0:num-1,*)
scale=size(slice304[0:num-1,*])
mul=10
wdef,2,scale(1)*mul,scale(2)
print,scale(1)*mul

im1=congrid(slice304[0:num-1,*],scale(1)*mul,scale(2),/interp)
help,im1 & pmm,alog(im1)
tv,bytscl(alog(im1),0.2,6.7)


save,filename='D:\Learning\PHD1st\magnetic_reconnecion\data_process\nfb.sav',slice304,sp,ov,mul,xx,yy


end


