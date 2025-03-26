
; pro .r fig5

; cd,'G:\2023.04.21_eruption\AAA-revised-figs\fig5_down&dis

;cd,'~/Desktop/AAA-revised-figs/fig5_down&dis'
cd,'D:\Learning\PHD1st\magnetic_reconnecion\program'
centerx=0  & centery=-200

xlc=300.0  &  xrc=20.0   &  ydc=250.0  &  yuc=85.0
xss0c=centerx-xlc & xss1c=centerx+xrc & yss0c=centery-ydc & yss1c=centery+yuc


xlb=225.0  &  xrb=-160.0   &  ydb=129.0  &  yub=-86.0
xss0b=centerx-xlb & xss1b=centerx+xrb & yss0b=centery-ydb & yss1b=centery+yub   

  charsize1=0.8
  charsize2=1.2





nx=3;-------小图片列数
ny=2;-------小图片行数
n=nx*ny;----小图片个数

;--------------------------这一段里的数字都是以一为整数来理解
xoff0=0.1;----------------图片左边界与EPS右边界的距离 
xoff1=0.05;----------------图片右边界与EPS右边界的距离 
yoff1=0.05;----------------图片上边界与EPS上边界的距离
yoff0=0.4;----------------图片下边界与EPS上边界的距离
x_inter=0.01;------------图片之间X方向上的间距
y_inter=0.01;------------图片之间Y方向上的间距

width_eps=19.0;-----------------------------EPS的实际宽度(cm)
xs=(1-xoff1-xoff0-(nx-1)*x_inter)/nx;-------每张图的宽度（总长为一）
ys=(1-yoff1-yoff0-(ny-1)*y_inter)/ny;-------每张图的高度（总长为一）
yx_ratio=(yuc+ydc)/(xrc+xlc);--------------------------------Y方向与X方向长度的缩放比例
height_eps=width_eps*xs*yx_ratio/ys;--------EPS的高度!!!注意这里是保证最后图形长宽比例与实际相符的重要条件。

set_plot,'ps'
device,filename='fig5.eps',xsize=width_eps,ysize=height_eps,xoff=(21.0-width_eps)/2.0,$
yoff=(29.7-height_eps)/2.0,/color,bits=8,/times,/bold,isolatin1=1


!p.font = 0


xyouts,0,0,'!3'+' ',charsize=charsize2,color=255,charthick=1,/normal


  
nx=2;-------小图片列数
ny=1;-------小图片行数
n=nx*ny;----小图片个数

;--------------------------这一段里的数字都是以一为整数来理解
xoff0=0.1;----------------图片左边界与EPS右边界的距离 
xoff1=0.05;----------------图片右边界与EPS右边界的距离 
yoff1=0.68;----------------图片上边界与EPS上边界的距离
yoff0=0.08;----------------图片下边界与EPS上边界的距离
x_inter=0.1;------------图片之间X方向上的间距
y_inter=0.00;------------图片之间Y方向上的间距

width_eps=19.0;-----------------------------EPS的实际宽度(cm)
xs=(1-xoff1-xoff0-(nx-1)*x_inter)/nx;-------每张图的宽度（总长为一）
ys=(1-yoff1-yoff0-(ny-1)*y_inter)/ny;-------每张图的高度（总长为一）
yx_ratio=(yuc+ydc)/(xrc+xlc);--------------------------------Y方向与X方向长度的缩放比例
height_eps=width_eps*xs*yx_ratio/ys;--------EPS的高度!!!注意这里是保证最后图形长宽比例与实际相符的重要条件。


for i=0,nx*ny-1,1 do begin
  x0=xoff0+(i mod nx)*(xs+x_inter);-----------------think about it,it's easy.
  y0=1-yoff1+y_inter-(byte(i/nx)+1)*(ys+y_inter);---same as above.
  x1=x0+xs
  y1=y0+ys
  
  
  if i eq 0 then begin
    restore,'D:\Learning\PHD1st\magnetic_reconnecion\program\nfb.sav' ;& iris_lct,'SJI_1400'
;    load,0    
    aia_lct,wave=171,/load
    tv,bytscl(alog(slice304[*,*]),1.0,8.0),x0,y0,/normal,xs=xs,ys=ys 
    load,0
    plots,[x0,x0,x1,x1,x0],[y1,y0,y0,y1,y1],linestyle=0,thick=3.5,color=0,/norm

  endif 
  


  
endfor 


device,/close
 set_plot,'win'



end