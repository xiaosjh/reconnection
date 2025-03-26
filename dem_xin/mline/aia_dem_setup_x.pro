;+ 
;===================================================================
; 
;     AIA_DEM_SETUP_X
;
;    constructs a temperature response structure for AIA
;    suitable for use in <xrt_dem_iterative2.pro>
;
;   INPUT:  none
;
;   OUTPUT:  tresp - temperature response structure for use
;                    with AIA data in <xrt_dem_iterative2.pro>
;
; NOTES:  be sure to :
; > setenv SSW_INSTR " $SSW_INSTR aia ontology "
; > source $SSW/gen/setup/setup.ssw
;
; before starting (ssw)IDL 
;
;  rebins responses to XRT temperature grid for potential 
;  use with joint data sets
;
;  written by SHS and Xin 07/11
;
;===================================================================
;-

pro aia_dem_setup_x,tresp

; time2='2023-31-Dec 21:40:00'     ; arbitrary obstime to get XRT tresp structure
time2='2010-Jun-22 08:00:00'     ; arbitrary obstime to get XRT tresp structure
w_resp2 = make_xrt_wave_resp(contam_time=time2)
trx = make_xrt_temp_resp(w_resp2, /apec)
tresp_old=trx 
nx=n_elements(trx.name)                ; 
logtx=trx.temp
lenx=26                                ; 


tresp0=aia_get_response(/temp, /dn,/chiantifix)


nc=n_elements(tresp0.channels)
logte=tresp0.logte
hist0=['Derived from aia_get_response','','']
units0='DN cm^5 s^-1 pix^-1'
descr0=' Provides temperature response function for an AIA channel, for a given plasma emission model.'
comm0='from aia_get_response; a copy of XRT tresp below length'
len0=101
typ0='AIA_TEMP_RESP'
chan0=tresp0.channels
dat0=tresp0.all

for j=0,nc-1 do begin
    trespj=tresp_old(0)
    trespj.type=typ0
    trespj.trs_str_descr=descr0
    trespj.name=chan0(j)
    trespj.temp=logtx(0:25)     ; was 10^logte
    tmp=reform(dat0(*,j))
    trespj.temp_resp=interpol(tmp,10.^logte,logtx(0:25))
    trespj.temp_resp_units=units0
    trespj.confg.name=chan0(j)
    trespj.length=lenx    
    trespj.history=hist0
    trespj.comments=[comm0,'','','',''] 
    if j eq 0 then tresp=trespj else tresp=[tresp,trespj]
endfor


return
end

