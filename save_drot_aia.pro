pro save_drot_aia ;非常无语无法正确保存fits，所以我直接在IDL里把这一小块区域保存为sav，然后再去读
; 创建一个动态大小的结构，用于存储不同大小的矩阵
matrix_struct = REPLICATE({matrix: PTR_NEW()}, 20)  ; 创建包含20个指针的结构数组

dir='D:\Learning\PHD1st\magnetic_reconnecion\data\AIA\AIA131rott\'
file_list=file_search(dir+'*fits')
for i=0,19 do begin
  fits2map,file_list[i],map0
  sub_map,map0,smap,xrange=[200,420],yrange=[-300,-520]
  datas=smap.data
  matrix_struct[i].matrix=ptr_new(datas)
endfor
save,matrix_struct,filename='D:\Learning\PHD1st\magnetic_reconnecion\data\AIA\rot_sub_aia_131.sav'
end