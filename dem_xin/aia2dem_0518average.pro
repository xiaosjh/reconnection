files = findfile('/Volumes/BACKUP2/scar_2nd/null_point/obs/aia_rain/rain1_dem/*.sav')

nn=size(files)

for kk=0,nn[1]-1 do begin

restore,files[kk]

n2 = size(result.smap_em)

min_tt=5
max_tt=15
logt_out = indgen(11)*0.1+6.0
temp=10.^logt_out

for ii=0,n2[1]-1 do begin
 for jj=0,n2[2]-1 do begin

  dem=result.dem_out[ii,jj,*]
  em=int_tabulated(temp[min_tt-min_tt:max_tt-min_tt],dem[min_tt:max_tt])
  temp1=int_tabulated(temp[min_tt-min_tt:max_tt-min_tt],(temp[min_tt-min_tt:max_tt-min_tt]*dem[min_tt:max_tt]))
  temp1=temp1/em
  result.smap_temp[ii,jj]=temp1
  result.smap_em[ii,jj]=em

 endfor
endfor

title=strmid(smap1[0].time,0,11)+'_'+strmid(smap1[0].time,12,2)+'_'+strmid(smap1[0].time,15,2)+'_'+strmid(smap1[0].time,18,2)
dir1 = '/Volumes/BACKUP2/scar_2nd/null_point/obs/aia_rain/rain1_dem_0518/'
save,result,smap1,filename=dir1+'DEM'+title+'.sav'

print,string(kk)+'has been done!!!!!!'

endfor







files = findfile('/Volumes/BACKUP2/scar_2nd/null_point/obs/aia_rain/rain2_dem/*.sav')

nn=size(files)

for kk=0,nn[1]-1 do begin

restore,files[kk]

n2 = size(result.smap_em)

min_tt=5
max_tt=15
logt_out = indgen(11)*0.1+6.0
temp=10.^logt_out

for ii=0,n2[1]-1 do begin
 for jj=0,n2[2]-1 do begin

  dem=result.dem_out[ii,jj,*]
  em=int_tabulated(temp[min_tt-min_tt:max_tt-min_tt],dem[min_tt:max_tt])
  temp1=int_tabulated(temp[min_tt-min_tt:max_tt-min_tt],(temp[min_tt-min_tt:max_tt-min_tt]*dem[min_tt:max_tt]))
  temp1=temp1/em
  result.smap_temp[ii,jj]=temp1
  result.smap_em[ii,jj]=em

 endfor
endfor

title=strmid(smap1[0].time,0,11)+'_'+strmid(smap1[0].time,12,2)+'_'+strmid(smap1[0].time,15,2)+'_'+strmid(smap1[0].time,18,2)
dir1 = '/Volumes/BACKUP2/scar_2nd/null_point/obs/aia_rain/rain2_dem_0518/'
save,result,smap1,filename=dir1+'DEM'+title+'.sav'

print,string(kk)+'has been done!!!!!!'

endfor

end
