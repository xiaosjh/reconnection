
; .r line

cd,'~/Desktop/dem/dem/5.5-7.5/line-'


dir1='~/Desktop/dem/dem/5.5-7.5/line-/'
dir2='~/Desktop/dem/dem/5.5-7.5/line-/'

filedem=findfile(dir1+'s*.sav')
filedem=filedem(sort(filedem))

n=n_elements(filedem)

ov=3

for iii=0,n-1,1 do begin
; iii=0
  restore,filedem[iii]


xran=5.5+findgen(21)*0.1
mdem=fltarr(21)
  for j=0,21-1,1 do begin
    mdem[j]=total(result.dem_out[*,*,j])*10e-23*0.1
  endfor
  print,mdem
wdef,iii
plot,mdem,psym=10
print,mean(mdem)
; filename=strmid(filedem[iii],51,17)
; save,mdem,filename=dir2+'line'+filename+'.sav'
; wdef,3

endfor




end