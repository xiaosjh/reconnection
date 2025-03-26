PRO iris_test
;这个文档是学习IRIS的tutorial
;filename='D:\Learning\PHD1st\magnetic_reconnecion\data\test\iris_l2_20140708_114109_3824262996_raster_t000_r00000.fits'
;d=iris_obj(filename)
;d.show_lines
;wave=d.getlam(7)
;data=d.getvar(7,/load)
;mspec=total(total(data,2),2)
;plot,wave,mspec,xrange=[2794,2799],/xst
;;pih,transpose(reform(data[350,*,*])),min=0,max=200,scale=[0.35,0.1667]
;wavecorr = iris_prep_wavecorr_l2(filename)
;for i=0, 399 do for j=0, 1091 do new_data[*, j, i] = interpol(data[*, j, i], wave + wavecorr.corr_nuv[i], wave) endfor endfor
filename='D:\Learning\PHD1st\magnetic_reconnecion\data\IRIS\iris_l2_20240618_163141_3602506433_SJI_1400_t000.fits'
d=iris_obj(filename)
matrix=d.getvar()
help,matrix
;save,matrix,filename='test.sav'
;restore,'test.sav'
;help,matrix
;filename = 'matrix_3d.dat'
;
; 打开文件以二进制写入
OPENW, unit, filename, /BINARY, /GET_LUN  ; 打开文件以二进制写入

; 将矩阵数据写入文件
; 注意：这里我们直接将整个矩阵作为一个大的数据块写入文件
; 在Python中读取这个文件时，你需要知道数据的维度和类型来正确地解析它
WRITEU, unit, matrix

 关闭文件
FREE_LUN, unit
end