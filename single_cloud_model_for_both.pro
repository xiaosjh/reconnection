function one_cloud_p,x,par
  common temp,ha0,light_speed,sp_back
  x0 = ha0
  cloud = par[0]*(1.-exp(-par[1]*exp(-(x-x0-x0/light_speed*par[2])^2/par[3]^2)))+par[4]
  return,cloud
end

function one_cloud_f,x,par
  common temp,ha0,light_speed,sp_back
  cloud = sp_back*par[4]*exp(-par[1]*exp(-(x-ha0-ha0/light_speed*par[2])^2/par[3]^2))+par[0]*(1.-exp(-par[1]*exp(-(x-ha0-ha0/light_speed*par[2])^2/par[3]^2)));+par[4]
  return,cloud
end

pro single_cloud_model_for_both
;;;;obtain the velocity of both filament and prominence based on the single-cloud model 
;;;;Written by Qiu Y.(NJU)
;;;Please cite Qiu Y. et al.(2024) https://ui.adsabs.harvard.edu/abs/2024ApJ...961L..30Q

;;;for prominence
;;;single cloud model
;;;par1[0,1,2,3] = [Source function,\tau_0,doppler_vel,width]

!p.background=255
loadct,0

;;;input;;;;;
dw        = 3           ;in pixel, half window size for calculating the sp_back
thres_cen = 0.1   ;threshold for centroid
thres_rib = 2.5   ;threshold for ribbon
global    = 0
xq        = [320,390]
yq        = [200,250]    ;in pixel,background region
vmin      = -50
vmax      = 50   ; for color
wind      = indgen(81)-81/2+67  ;wavelength windowns for centroid

path_read = 'User_path'
path_read2= 'User_path'
path_sav0  = 'User_path'
file_mkdir,path_sav0
fmap      = file_search(path_read,'*map.sav')

common temp,ha0,light_speed,sp_back
light_speed = 3e5

fcube=file_search(path_read,'data*.sav')
freg = file_search(path_read2,'*.sav',count=num)

;for nn=8,num-1 do begin
for nn=0,num-1 do begin
  path_sav=path_sav0+strtrim(string(nn,format='(i02)'),2)+'/'
  file_mkdir,path_sav
  
  data = !null
  reg_cube = !null  
  ;;obtain the filament region
  restore,freg[nn]   ;reg_cube 
  reg_cube= loc_f 
  sz = size(reg_cube)
  
  ;;;;obtain the datacube
  ;fcube=file_search(path_read+strtrim(string(nn,format='(i02)'),2),'data*.sav')
  ;restore,fcube[0]
  restore,fcube[nn]
  err   = replicate(22.42,header[0].naxis3)    ;standard derivation of dark field
  wave0 = header[0].CDELT3*findgen(header[0].naxis3)+header[0].CRVAL3
  ;;;;;;;get the sp_back and Ha0 at the first series;;;;;
  data = reform(data_cube[*,*,*,0])
  if global then begin
    counts = 0
    sp_back = fltarr(header[0].naxis3)
    zoom   = 1
    while counts eq 0 do begin
      for j = max([0,jj-dw*zoom]),min([jj+dw*zoom,sz[2]-1]) do begin
        for i = max([0,ii-dw*zoom]),min([ii+dw*zoom,sz[1]-1]) do begin
          if reg_cube[i,j,kk] eq 0 && mean(data[i,j,*]) ge 200 then begin
            counts=counts+1
            sp_back=sp_back+reform(data[i,j,*])
          endif
        endfor
      endfor
      zoom = zoom+1
    endwhile
    sp_back = sp_back/counts
    pc0     = wave_centroid2(wind,sp_back[wind],thres_cen)
    ha0     = wave_centroid2(wave0[wind],sp_back[wind],thres_cen)
  endif else begin
    if n_elements(xq) eq 0 then begin
      wdef,1,900,900
      plot_image,reform(data[*,*,67])
      xq = intarr(2)
      yq = intarr(2)
      for i=0,1 do begin
        cursor,a,b
        xq[i] = round(a)
        yq[i] = round(b)
      endfor
      print,'xq =',xq
      print,'yq =',yq
      sp_back = mean(data[min(xq):max(xq),min(yq):max(yq),*],dimension=1)
      sp_back = mean(sp_back,dimension=1)
      pc0     = wave_centroid2(wind,sp_back[wind],thres_cen)
      ha0     = wave_centroid2(wave0[wind],sp_back[wind],thres_cen)
    endif else begin
      sp_back = mean(data[min(xq):max(xq),min(yq):max(yq),*],dimension=1)
      sp_back = mean(sp_back,dimension=1)
      pc0     = wave_centroid2(wind,sp_back[wind],thres_cen)
      ha0     = wave_centroid2(wave0[wind],sp_back[wind],thres_cen)
    endelse
  endelse
 
  eis_colors,/velocity,/dark
  !p.background=127
  for kk=0,sz[3]-1 do begin
    vel  = fltarr(sz[1],sz[2])
    verr = fltarr(sz[1],sz[2])
    soufunc = fltarr(sz[1],sz[2])   ;source function
    chisq_2d = fltarr(sz[1],sz[2])
    
    ;restore,fcube[kk]
    data = reform(data_cube[*,*,*,kk])
    for jj =0,sz[2]-1 do begin
      for ii=0,sz[1]-1 do begin
        
        ;;;;;;;get the sp_back and Ha0;;;;;
        if global then begin
          counts = 0
          sp_back = fltarr(header[0].naxis3)
          zoom   = 1
          while counts eq 0 do begin
            for j = max([0,jj-dw*zoom]),min([jj+dw*zoom,sz[2]-1]) do begin
              for i = max([0,ii-dw*zoom]),min([ii+dw*zoom,sz[1]-1]) do begin
                if reg_cube[i,j,kk] eq 0 && mean(data[i,j,*]) ge 200 then begin
                  counts=counts+1
                  sp_back=sp_back+reform(data[i,j,*])
                endif
              endfor
            endfor
            zoom = zoom+1
          endwhile
          sp_back = sp_back/counts
          pc      = wave_centroid2(wind,sp_back[wind],thres_cen)
          wave    = wave0 - header[0].CDELT3*(pc-pc0)
          ha0     = wave_centroid2(wave[wind],sp_back[wind],thres_cen)
        endif else begin
          sp_back = mean(data[min(xq):max(xq),min(yq):max(yq),*],dimension=1)
          sp_back = mean(sp_back,dimension=1)
          pc     = wave_centroid2(wind,sp_back[wind],thres_cen)
          wave    = wave0 - header[0].CDELT3*(pc-pc0)
          ha0     = wave_centroid2(wave[wind],sp_back[wind],thres_cen)
        endelse
                 
        if reg_cube[ii,jj,kk] ne 0 then begin
          sp = reform(data[ii,jj,*])
          ;;;;for prominence;;;;
          if reg_cube[ii,jj,kk] eq 2 then begin 
            linecore = wave_centroid2(wave[wind],sp[wind],0.5)
            width = spec_width(wave[wind],sp[wind],linecore,0.5)
            vtest = (linecore-ha0)/ha0*light_speed
            start = [max(sp),1,vtest,width,30]
            pi    = replicate({fixed:0,limited:[0,0],limits:[0.D,0.d]},(size(start))[1])
            pi[-1].limited[0]=1
            pi[-1].limits[0]=10
            result = mpfitfun('one_cloud_p',wave,sp,err,start,yfit=yfit,perror=perr,parinfo=pi,bestnorm=chisq,/quiet)
            if chisq ge 20 then begin
              if mean(sp[0:67]) gt mean(sp[67:*]) then begin
                start = [max(sp),1,-50,0.3,100]
              endif else begin
                start = [max(sp),1,50,0.3,100]
              endelse
              result = mpfitfun('one_cloud_p',wave,sp,err,start,yfit=yfit,perror=perr,parinfo=pi,bestnorm=chisq,/quiet)
            endif          
            vel[ii,jj] = result[2]
            verr[ii,jj] = perr[2]
            soufunc[ii,jj] = result[0]
            chisq_2d[ii,jj] = chisq          
          endif
          
          ;;;;for filament;;;;;;;;
          sp_back0=sp_back
          if reg_cube[ii,jj,kk] eq 1 then begin
            vtest = (wave_centroid2(wave[wind],sp[wind],0.2)-ha0)/ha0*light_speed
            start = [min(sp),1,vtest,0.3,1]
            pi    = replicate({fixed:0,limited:[0,0],limits:[0.D,0.d]},(size(start))[1])
            pi[0].limited=[1,1]
            pi[0].limits = [0,max(sp)]
            ;      pi[2].limited=[1,1]
            ;      pi[2].limits = [0,vtest*3]
            pi[3].limited[0] =1
            pi[3].limits[0]=0
            if abs(vtest) ge 10 then begin
              result = mpfitfun('one_cloud_f',wave,sp,err,start,yfit=yfit,perror=perr,parinfo=pi,bestnorm=chisq,/quiet,dof=dof)
            endif else begin
              sp_back = sp_back[wind]
              result = mpfitfun('one_cloud_f',wave[wind],sp[wind],err[wind],start,yfit=yfit,perror=perr,parinfo=pi,bestnorm=chisq,/quiet,dof=dof)
            endelse
            vel[ii,jj] = result[2]
            verr[ii,jj] = perr[2]
            chisq_2d[ii,jj] = chisq
            if result[0] eq 0 then begin
              vel[ii,jj] = vtest
              verr[ii,jj] = 0
              chisq_2d[ii,jj] =0
            endif
          endif
        endif else begin
          sp = reform(data[ii,jj,*])
          if mean(sp) ge 200 then begin
            vel[ii,jj]   = (wave_centroid2(wave[wind],sp[wind],0.2)-ha0)/ha0*light_speed
          endif
        endelse
        
        ;;;;for flare ribbon;;;;;;;;;
        img = reform(data[*,*,67])
        bg_int =    mean(img[min(xq):max(xq),min(yq):max(yq)])  
        if img[ii,jj] ge bg_int*thres_rib then begin
          if n_elements(sp_back) eq n_elements(wind) then begin
            vel[ii,jj]   = (wave_centroid2(wave[wind],sp[wind]-sp_back,0.2,/contrastp)-ha0)/ha0*light_speed
          endif else begin
            vel[ii,jj]   = (wave_centroid2(wave[wind],sp[wind]-sp_back[wind],0.2,/contrastp)-ha0)/ha0*light_speed
          endelse
          verr[ii,jj] = 0
          chisq_2d[ii,jj] =0 
        endif     
      endfor
      print, string(float(jj)/float(sz[2]-1)*100)+'% fitting of one scanning sequence is done.'
    endfor
    save,vel,verr,soufunc,chisq_2d,xran,yran,date,header,sp_back0,filename=path_sav+'velocity'+string(kk,format='(i02)')+'.sav'
    
    ;;;plot_image;;;;;;;
    wdef,1,800,500
    eis_colors,/velocity,/dark
    plot_image,vel,xstyle=5,ystyle=5,position=[.1,.1,.95,.95],min=vmin,max=vmax
    plot,xran,yran,xstyle=1,ystyle=1,position=[.1,.1,.95,.95],charsize=2,xtitle='X (arcsec)',ytitle='Y (arsce)',$
      /noerase,color=0,/nodata,charthick=2
    xyouts,min(xran)+20,min(yran)+20,date[kk],charsize=2.5,color=0,charthick=2   
    write_png,path_sav+'cloud'+string(nn,format='(i02)')+'_'+string(kk,format='(i02)')+'.png',tvrd(true=1)
  endfor
endfor

end
