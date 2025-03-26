; .pro .r fig4

cd,'D:\Learning\PHD1st\magnetic_reconnecion\program\'
;
;
;centerx=0  & centery=0
;
;xldem=1050.0  &  xrdem=-1000.0   &  yddem=-35.0  &  yudem=80.0
;xsdems0=centerx-xldem & xsdems1=centerx+xrdem & ysdems0=centery-yddem & ysdems1=centery+yudem 
;
;
;file211=findfile('newaia.*.211.image_lev1.fits')
;file211=file211[sort(file211)]

filedem=findfile('DEM*.sav')
filedem=filedem[sort(filedem)]

;fsave=findfile('*.save')
;fsave=fsave[sort(fsave)]


  charsize1=1.0
  charsize2=1.2
  


nx=3;-------小图片列数
ny=2;-------小图片行数
n=nx*ny;----小图片个数

;--------------------------这一段里的数字都是以一为整数来理解
xoff0=0.1;----------------图片左边界与EPS右边界的距离 
xoff1=0.1;----------------图片右边界与EPS右边界的距离 
yoff1=0.1;----------------图片上边界与EPS上边界的距离
yoff0=0.1;----------------图片下边界与EPS上边界的距离
x_inter=0.02;------------图片之间X方向上的间距
y_inter=0.05;------------图片之间Y方向上的间距

width_eps=19.0;-----------------------------EPS的实际宽度(cm)
xs=(1-xoff1-xoff0-(nx-1)*x_inter)/nx;-------每张图的宽度（总长为一）
ys=(1-yoff1-yoff0-(ny-1)*y_inter)/ny;-------每张图的高度（总长为一）
yx_ratio=1;(yudem+yddem)/(xrdem+xldem);--------------------------------Y方向与X方向长度的缩放比例
height_eps=width_eps*xs*yx_ratio/ys;--------EPS的高度!!!注意这里是保证最后图形长宽比例与实际相符的重要条件。


set_plot,'ps'
device,filename='fig4.eps',xsize=width_eps,ysize=height_eps,xoff=(21.0-width_eps)/2.0,$
yoff=(29.7-height_eps)/2.0,/color,bits=8,/times,/bold,isolatin1=1

!p.font = 0


xyouts,0,0,'!3'+' ',charsize=charsize2,color=255,charthick=1,/normal



for i=0,nx*ny-1,1 do begin
  x0=xoff0+(i mod nx)*(xs+x_inter);-----------------think about it,it's easy.
  y0=1-yoff1+y_inter-(byte(i/nx)+1)*(ys+y_inter);---same as above.
  x1=x0+xs
  y1=y0+ys  

;
;;figs211--------------------------------------------------
;  if i eq 0 then begin
;    fits2map,file211[0],map
;    rmap=get_sub_map(map,xrange=[xsdems0,xsdems1],yrange=[ysdems0,ysdems1])
;    print,map.dur
;    rmap=rmap.data/map.dur
;    aia_lct,wave=211,/load
;    tv,bytscl(alog(rmap),5.5,7.0),x0,y0,/normal,xs=xs,ys=ys   ;1-7-1 
;    load,0
;    plots,[x0,x0,x1,x1,x0],[y1,y0,y0,y1,y1],linestyle=0,thick=1.5,color=0,/norm
;    pmm,alog(rmap)
;    ;----------------------------------------------------
;    time0=map.time
;    xyouts,x0+xs*0.75,y1-ys*0.1,strmid(time0,12,8)+ ' UT',charsize=charsize1,color=255,charthick=1.5,/normal,alignment=0.5  
;    ;----------------------------------------------------
;    load,0
;    xyouts,x0+xs*0.5,y1+ys*0.03,'AIA 211 '+string(byte(197)),charsize=charsize1,color=0,charthick=1.5,/normal,font=0,alignment=0.5
;  endif
;  if i eq 3 then begin
;    fits2map,file211[1],map
;    rmap=get_sub_map(map,xrange=[xsdems0,xsdems1],yrange=[ysdems0,ysdems1])
;    print,map.dur
;    rmap=rmap.data/map.dur
;    aia_lct,wave=211,/load
;    tv,bytscl(alog(rmap),5.5,7.0),x0,y0,/normal,xs=xs,ys=ys   ;1-7-1 
;    load,0
;    plots,[x0,x0,x1,x1,x0],[y1,y0,y0,y1,y1],linestyle=0,thick=1.5,color=0,/norm
;    pmm,alog(rmap)
;    ;----------------------------------------------------
;    time0=map.time
;    xyouts,x0+xs*0.75,y1-ys*0.1,strmid(time0,12,8)+ ' UT',charsize=charsize1,color=255,charthick=1.5,/normal,alignment=0.5  
;    ;----------------------------------------------------
;    load,0
;    xyouts,x0+xs*0.5,y1+ys*0.03,'AIA 211 '+string(byte(197)),charsize=charsize1,color=0,charthick=1.5,/normal,font=0,alignment=0.5
;  endif

;tem--------------------------------------------------

if i eq 0 then begin
  restore,filedem[0]
  aia_lct,wave=131,/load
  tv,bytscl(alog10(smap1.data),0.5,2.9),x0,y0,/normal,xs=xs,ys=ys
  load,0
  plots,[x0,x0,x1,x1,x0],[y1,y0,y0,y1,y1],linestyle=0,thick=1.5,color=0,/norm
endif
  if i eq 1 then begin
    restore,filedem[0]
    load,4  
    ; tv,bytscl(result.smap_temp,3.5e+06,9.5e+06),x0,y0,/normal,xs=xs,ys=ys
    tv,bytscl(alog10(result.smap_temp),5.5,6.9),x0,y0,/normal,xs=xs,ys=ys
    load,0
    plots,[x0,x0,x1,x1,x0],[y1,y0,y0,y1,y1],linestyle=0,thick=1.5,color=0,/norm
    ;----------------------------------------------------
    load,0
    xyouts,x0+xs*0.75,y1-ys*0.1,'21:07:06'+ ' UT',charsize=charsize1,color=255,charthick=1.5,/normal,alignment=0.5  
    xyouts,x0+xs*0.55,y1+ys*0.03,'Map of lg Temp. (K)',charsize=charsize1,color=0,charthick=1.5,/normal,font=0,alignment=0.5
  endif
;  if i eq 4 then begin
;    restore,filedem[1]
;    load,4  
;    ; tv,bytscl(result.smap_temp,3.5e+06,9.5e+06),x0,y0,/normal,xs=xs,ys=ys
;    tv,bytscl(alog10(result.smap_temp),6.5,6.9),x0,y0,/normal,xs=xs,ys=ys
;    load,0
;    plots,[x0,x0,x1,x1,x0],[y1,y0,y0,y1,y1],linestyle=0,thick=1.5,color=0,/norm
;    ;----------------------------------------------------
;    load,0
;    xyouts,x0+xs*0.75,y1-ys*0.1,'21:42:42'+ ' UT',charsize=charsize1,color=255,charthick=1.5,/normal,alignment=0.5  
;    xyouts,x0+xs*0.55,y1+ys*0.03,'Map of lg Temp. (K)',charsize=charsize1,color=0,charthick=1.5,/normal,font=0,alignment=0.5
;  endif
; 
;  if i eq 4 then begin
;    load,4
;    colorbar,ncolors=256,position=[x0+xs*0.05,y0-ys*0.04,x1-xs*0.05,y0-ys*0.01],/horizontal,$;range=[0,42000] ,/vertical,/right$
;    ; colorbar,ncolors=256,position=[x1+xs*0.02,y0+ys*0.00,x1+xs*0.035,y1-ys*0.00],/vertical,/right,$;range=[0,42000],/horizontal $
;    divisions=2,ticknames=['6.5','6.7','6.9'],title=' ',color=0,charsize=charsize1-0.4,charthick=1.5,FONT=0,/norm   
;  endif 


;em--------------------------------------------------
  if i eq 2 then begin
    restore,filedem[0]
    load,22  
    ; tv,bytscl(result.smap_em,2.6e+27,2.8e+28),x0,y0,/normal,xs=xs,ys=ys
    tv,bytscl(alog10(result.smap_em),25.3,29.3),x0,y0,/normal,xs=xs,ys=ys
    load,0
    plots,[x0,x0,x1,x1,x0],[y1,y0,y0,y1,y1],linestyle=0,thick=1.5,color=0,/norm
    ;----------------------------------------------------
    load,0
    xyouts,x0+xs*0.75,y1-ys*0.1,'21:07:06'+ ' UT',charsize=charsize1,color=0,charthick=1.5,/normal,alignment=0.5  
    xyouts,x0+xs*0.55,y1+ys*0.03,'Map of lg EM (cm!E-5!X!N)',charsize=charsize1,color=0,charthick=1.5,/normal,font=0,alignment=0.5
  endif
;  if i eq 5 then begin
;    restore,filedem[1]
;    load,22  
;    ; tv,bytscl(result.smap_em,2.6e+27,2.8e+28),x0,y0,/normal,xs=xs,ys=ys
;    tv,bytscl(alog10(result.smap_em),27.3,28.3),x0,y0,/normal,xs=xs,ys=ys
;    load,0
;    plots,[x0,x0,x1,x1,x0],[y1,y0,y0,y1,y1],linestyle=0,thick=1.5,color=0,/norm
;    ;----------------------------------------------------
;    load,0
;    xyouts,x0+xs*0.75,y1-ys*0.1,'21:43:06'+ ' UT',charsize=charsize1,color=0,charthick=1.5,/normal,alignment=0.5  
;    xyouts,x0+xs*0.55,y1+ys*0.03,'Map of lg EM (cm!E-5!X!N)',charsize=charsize1,color=0,charthick=1.5,/normal,font=0,alignment=0.5
;  endif
;  if i eq 5 then begin
;    load,22
;    colorbar,ncolors=256,position=[x0+xs*0.05,y0-ys*0.04,x1-xs*0.05,y0-ys*0.01],/horizontal,$;range=[0,42000] ,/vertical,/right$
;    ; colorbar,ncolors=256,position=[x1+xs*0.02,y0+ys*0.00,x1+xs*0.035,y1-ys*0.00],/vertical,/right,$;range=[0,42000],/horizontal $
;    divisions=2,ticknames=['27.3','27.8','28.3'],title=' ',color=255,charsize=charsize1-0.4,charthick=1.5,FONT=0,/norm   
;  endif 



;  if i le 2 then begin
;    fits2map,file211[0],map0
;    rmap=get_sub_map(map0,xrange=[xsdems0,xsdems1],yrange=[ysdems0,ysdems1])
;    rmap=rmap.data/map0.dur
;    load,0
;    cols=[0,255,0]
;    rmap[*,50:*]=0
;    contour,smooth(alog(rmap),0),xsty=5,ysty=5,pos=[x0,y0,x0+xs,y0+ys],levels=[6.5],color=cols[i],c_thick=[2],/noerase,/norm
;  endif
;
;  if i gt 2 then begin
;    fits2map,file211[1],map0 
;    rmap=get_sub_map(map0,xrange=[xsdems0,xsdems1],yrange=[ysdems0,ysdems1])
;    rmap=rmap.data/map0.dur
;    load,0
;    cols=[0,255,0]
;    rmap[*,50:*]=0
;    ; rmap[*,0:2]=0
;    contour,smooth(alog(rmap),0),xsty=5,ysty=5,pos=[x0,y0,x0+xs,y0+ys],levels=[6.7],color=cols[i-3],c_thick=[2],/noerase,/norm
;  endif
;
;

    load,0
    lab=['!4a1!X','!4a2!X','!4a3!X',$
         '!4b1!X','!4b2!X','!4b3!X']
    xyouts,x0+xs*0.05,y1+ys*0.03,lab[i],charsize=charsize2,color=0,charthick=1.5,/normal,font=0,alignment=0.5

;  if i eq 3 then begin  
;    load,0    
;    plot,[1],[1],pos=[x0,y0,x1,y1],/normal,/nodata,$
;    xtitle='Solar X (arcsec)',xthick=1.5,xrange=[xsdems0,xsdems1],xstyle=1+8,xtickinterval=20,xminor=4,$        
;    ytitle='Solar Y (arcsec)',ythick=1.5,yrange=[ysdems0,ysdems1],ystyle=1+8,ytickinterval=20,yminor=4,$        
;    charsize=charsize1-0.4,charthick=1.5,font=0,background=255,color=0,/noerase,ticklen=-0.03
;  endif   
;  
  
endfor


device,/close
set_plot,'win'



end

