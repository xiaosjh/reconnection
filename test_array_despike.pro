;pro test_array_despike
;
;  file = 'C:\Learning\PHD1st\magnetic_reconnecion\data\IRIS\iris_l2_20240618_163141_3602506433_raster\iris_l2_20240618_163141_3602506433_raster_t000_r00392.fits'
;
;  data = mrdfits(file, 2, header, /silent)
;  spectrum = float(data[*, *, 4])
;
;  ; 清理 NaN、Inf、-200 和负强度
;  bad = where((finite(spectrum) EQ 0) OR (spectrum LT 0), nbad)
;  if nbad GT 0 then spectrum[bad] = 0.0
;
;  print, 'Bad or negative pixels = ', nbad
;
;  spectrum_original = spectrum
;
;  spectrum_clean = array_despike( $
;    spectrum, $
;    sigma=8.0, $
;    threshold=6.0, $
;    itmax=20, $
;    /NOLOW, $
;    /VERBOSE $
;    )
;
;  changed = where(spectrum_clean NE spectrum_original, nchanged)
;
;  ntotal = n_elements(spectrum)
;
;  print, 'Total pixels   = ', ntotal
;  print, 'Changed pixels = ', nchanged
;  print, 'Changed ratio  = ', 100.0d * double(nchanged) / double(ntotal), ' %'
;
;end
pro test_array_despike

  file = 'C:\Learning\PHD1st\magnetic_reconnecion\data\IRIS\iris_l2_20240618_163141_3602506433_raster\iris_l2_20240618_163141_3602506433_raster_t000_r00392.fits'

  ; /FSCALE 会应用 FITS 头中的 BSCALE 和 BZERO
  data = mrdfits(file, 2, header, /silent, /fscale)
  spectrum = data[*, *, 4]
  

  print, min(spectrum, /nan)
  print, max(spectrum, /nan)
  print, mean(spectrum, /nan)

  spectrum_original = spectrum

  spectrum_clean = array_despike( $
      spectrum, $
      sigma=8.0, $
      threshold=6.0, $
      itmax=20, $
      /NOLOW, $
      /VERBOSE $
  )

  changed = where(spectrum_clean NE spectrum_original, nchanged)

  print, 'Changed pixels = ', nchanged
  print, 'Changed ratio  = ', $
      100.0d * nchanged / n_elements(spectrum), ' %'

end