pro process_iris
;---------------------------------------------------------------------------------------------------------
;matrix_struct = REPLICATE({matrix: PTR_NEW()}, 990)
;restore,'D:\Learning\PHD1st\magnetic_reconnecion\data_process\iris_1400_map.sav'
;;values=[399*2+1,403*2+1,404*2+1,423*2+1]
;for i=0,989 do begin
;  sub_map,map0[i],smap,xrange=[320,390],yrange=[-352,-462]
;  smap.data=smap.data/smap.dur
;  matrix_struct[i].matrix=ptr_new(smap.data)
;  print,i
;endfor
;save,matrix_struct,filename='D:\Learning\PHD1st\magnetic_reconnecion\data_process\all_sji_image.sav'
;---------------------------------------------------------------------------------------------------------
;;这一部分是用来得到全部时刻狭缝的solar_x,solar_y
;dir='D:\Learning\PHD1st\magnetic_reconnecion\data\IRIS\iris_l2_20240618_163141_3602506433_raster\'
;f=iris_files(path=dir)
;solar_x=fltarr(495,8);58是时间序列上有58幅图像，8是每一幅图有8个solar_x的坐标
;solar_y=fltarr(495,778)
;;取三个时刻，分别是799（399）对应的21：09：27，
;;809（404）对应的21：12：56，847（423）对应的21：26：08
;;values=[399,403,404,423]
;for i=0,494 do begin
;  d=iris_load(f[i])
;  solar_x[i,*]=d.getxpos()
;  solar_y[i,*]=d.getypos()
;  print,i
;endfor
;save,filename='D:\Learning\PHD1st\magnetic_reconnecion\data_process\all_solar_xy.sav',solar_x,solar_y
;----------------------------------------------------------------------------------------------------------
restore,'D:\Learning\PHD1st\magnetic_reconnecion\data_process\iris_1400_map.sav'
 ;假设 data1 已经是一个 3D 数组
dir='D:\Learning\PHD1st\magnetic_reconnecion\data\IRIS\iris_l2_20240618_163141_3602506433_raster\'
f=iris_files(path=dir)
d=iris_load(f[399])
d.show_lines
data1=d.getvar(3,/load)
raster1 = data1[*, 590, 4]  ; 提取一维数组，与 Python 的 data1[:, 590, 4] 对应
raster2 = data1[*,591,4]
; 创建 x 数组
x = FLTARR(165)
FOR i = 0, 164 DO x[i] = 1398.2875 + 0.05058 * i

; 定义 x1 和 x2
x1 = 75
x2 = 95

; 绘制曲线
PLOT, x[x1:x2], (raster1[x1:x2]+raster2[x1:x2])/2, /YSTYLE, TITLE='T21:26:08', XTITLE='X-axis', YTITLE='Y-axis'
print,mean(raster1)

; 添加竖直虚线
xline = 1402.77
PLOTS, [xline, xline], [MIN(raster1[x1:x2]), MAX(raster1[x1:x2])], LINESTYLE=2

end