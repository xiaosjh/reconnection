
; .r bmp

cd,'~/Desktop/dem/dem/5.5-7.5/line-'


dir1='~/Desktop/dem/dem/5.5-7.5/line-/'
dir2='~/Desktop/dem/dem/5.5-7.5/line-/'



filedem=findfile(dir1+'s*.sav')
filedem=filedem(sort(filedem))

n=n_elements(filedem)

ov=20
; for iii=0,0,1 do begin
; 
for iii=0,n-1,1 do begin
  restore,filedem[iii]
  ss=size(result.dem_out)
  wdef,0,ss(1)*ov*3,ss(2)*ov*2
  print,ss(1)*ov

;1.9e+06 1.7e+07----1.5e06 2.5e06------------------------------------------------------------------------------------
  aia_lct,wave=94,/load
  tv,bytscl(rebin(alog(smap20.data),ss(1)*ov,ss(2)*ov),3.0,5.0),0,0,/normal,xs=ss(1)*ov,ys=ss(2)*ov  
  load,4  
  ; tv,bytscl(rebin(alog10(result.smap_temp),ss(1)*ov,ss(2)*ov),6.8,7.2),1.0/3.0,0,/normal,xs=ss(1)*ov,ys=ss(2)*ov  
  tv,bytscl(rebin(alog10(result.smap_temp),ss(1)*ov,ss(2)*ov),6.7,7.2),1.0/3.0,0,/normal,xs=ss(1)*ov,ys=ss(2)*ov  
  load,22  
  ; tv,bytscl(rebin(result.smap_em,ss(1)*ov,ss(2)*ov),2.5e+27,40e+27),2.0/3.0,0,/normal,xs=ss(1)*ov,ys=ss(2)*ov  
  tv,bytscl(rebin(alog10(result.smap_em),ss(1)*ov,ss(2)*ov),27.3,28.7),2.0/3.0,0,/normal,xs=ss(1)*ov,ys=ss(2)*ov  

; ;--2.3e+27 7.5e+27----------------------------------------------------------------------
; filterset=['A131','A94','A211','A193','A171','A335']
  aia_lct,wave=171,/load
  tv,bytscl(rebin(alog(smap50.data),ss(1)*ov,ss(2)*ov),6.0,7.5),0,1.0/2.0,/normal,xs=ss(1)*ov,ys=ss(2)*ov  
  aia_lct,wave=131,/load
  tv,bytscl(rebin(alog(smap10.data),ss(1)*ov,ss(2)*ov),4.0,7.0),1.0/3,1.0/2.0,/normal,xs=ss(1)*ov,ys=ss(2)*ov  
  aia_lct,wave=335,/load
  tv,bytscl(rebin(alog(smap60.data),ss(1)*ov,ss(2)*ov),3.0,5.0),2.0/3,1.0/2.0,/normal,xs=ss(1)*ov,ys=ss(2)*ov  



TVLCT,255,255,255,0
contour,smooth(rebin(alog10(result.smap_em),ss(1)*ov,ss(2)*ov),0),xsty=5,ysty=5,pos=[0.0/3,1.0/2,1.0/3,2.0/2],levels=[28.3],c_colors=[0],c_thick=[2],/noerase,/norm
contour,smooth(rebin(alog10(result.smap_em),ss(1)*ov,ss(2)*ov),0),xsty=5,ysty=5,pos=[1.0/3,1.0/2,2.0/3,2.0/2],levels=[28.3],c_colors=[0],c_thick=[2],/noerase,/norm
contour,smooth(rebin(alog10(result.smap_em),ss(1)*ov,ss(2)*ov),0),xsty=5,ysty=5,pos=[2.0/3,1.0/2,3.0/3,2.0/2],levels=[28.3],c_colors=[0],c_thick=[2],/noerase,/norm

contour,smooth(rebin(alog10(result.smap_em),ss(1)*ov,ss(2)*ov),0),xsty=5,ysty=5,pos=[0.0/3,0.0/2,1.0/3,1.0/2],levels=[28.3],c_colors=[0],c_thick=[2],/noerase,/norm
contour,smooth(rebin(alog10(result.smap_em),ss(1)*ov,ss(2)*ov),0),xsty=5,ysty=5,pos=[1.0/3,0.0/2,2.0/3,1.0/2],levels=[28.3],c_colors=[0],c_thick=[2],/noerase,/norm
contour,smooth(rebin(alog10(result.smap_em),ss(1)*ov,ss(2)*ov),0),xsty=5,ysty=5,pos=[2.0/3,0.0/2,3.0/3,1.0/2],levels=[28.3],c_colors=[0],c_thick=[2],/noerase,/norm


  load,0
  xyouts,0.45,0.02,strmid(filedem[iii],54,8),charsize=3,color=255,charthick=2,alignment=0 ,/normal ;,font=0
  filename=strmid(filedem[iii],51,17)+'.bmp'
  write_bmp,dir2+filename,tvrd(/true),/rgb  
print,mean(result.smap_em)

endfor

end