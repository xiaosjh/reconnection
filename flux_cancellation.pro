pro flux_cancellation
restore,'D:\Learning\PHD1st\magnetic_reconnecion\data_process\drot_sub_hmi.sav'
array_pos=fltarr(54)
array_neg=fltarr(54)
l = 151963542123.29977 * 100
r = l * sin(!pi / (3600 * 180.0))
cm2 = (r * 0.504)^2

for i=0,53 do begin
  mapp=map_list[i]
  sub_map,mapp,ssmap,xrange=[340,371],yrange=[-397,-377]
  ;sub_map,mapp,ssmap,xrange=[338,372],yrange=[-400,-377]
  data=ssmap.data
  array_pos[i]=total(data* (data gt 0))*cm2
  array_neg[i]=total(data* (data lt 0))*cm2
endfor
x=findgen(54)

; 起始时间 -> Julian Date
start_jd = JULDAY(2024, 6, 18, 20, 59, 53)

; 每隔 45 秒
interval_sec = 45.0 / 86400.0  ; 秒转天

; 生成时间序列 (长度 54)
time_jd = start_jd + FINDGEN(54) * interval_sec

; 转字符串 (想在横坐标上写字符串用这个)
; 定义长度54的字符串数组
time_str = STRARR(54)

; 循环把每个 JD 转成字符串
FOR i=0,53 DO BEGIN
  CALDAT, time_jd[i], month, day, year, hour, min, sec
  time_str[i] = STRING(hour, FORMAT='(I2.2)') + ':' + $
                STRING(min, FORMAT='(I2.2)') + ':' + $
                STRING(FIX(sec), FORMAT='(I2.2)')
ENDFOR

; 生成刻度位置
tick_pos = [0,8,17,26,35,44,53]  ; 54个点，均匀选5个
tick_name = time_str[tick_pos]

; 画图
PLOT, FINDGEN(54), array_neg, XTICKS=6, XTICKV=tick_pos, XTICKNAME=tick_name, $
  XTITLE='Time', YTITLE='Value'
save,array_pos,array_neg,filename='D:\Learning\PHD1st\magnetic_reconnecion\data_process\flux_cancellation.sav'
end