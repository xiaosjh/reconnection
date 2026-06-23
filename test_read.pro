pro test_read

  file = 'C:\Learning\PHD1st\magnetic_reconnecion\data\IRIS\iris_l2_20240618_163141_3602506433_raster\iris_l2_20240618_163141_3602506433_raster_t000_r00392.fits'

  data = mrdfits(file, 2, header, /silent, /fscale)
  spectrum = float(data[*, *, 4])
  
  kernel = [[1,1,1], [1,0,1], [1,1,1]]
  
  asm = convol(spectrum, kernel, /norm, /edge_truncate)
diff = spectrum - asm

; 防止对负数开平方
denom = sqrt(asm > 0.0) > 1.0
sigma_array = f_div(diff, denom)

repl = where((sigma_array GT 8.0) AND (spectrum GE 6.0), nreplace)

print, 'First-pass replacements = ', nreplace
writefits, 'C:\Learning\sigma_idl.fits', sigma_array

end