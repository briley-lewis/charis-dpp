pro charis_pdi_subtract_sky,pfname,medbox=medbox,prefname=prefname,suffname=suffname,pick=pick,scalesky=scalesky,scaleindchan=scaleindchan,skyrad=skyrad,skylim=skylim,ladder=ladder,verbose=verbose,help=help

if (N_PARAMS() eq 0 or keyword_set(help)) then begin
    print,"charis_pdi_subtract_sky,pfname,prefname=prefname,suffname=suffname,pick=pick,scalesky=scalesky,scaleindchan=scaleindchan,skyrad=skyrad,verbose=verbose,help=help"
    print," "
    print,"Example: charis_pdi_subtract_sky,'HR8799_low.info',/scalesky"
    print,""
    print,"***Important Keywords***"
    print,'*pfname - parameter file (e.g. HR8799_low.info)'
    print,"pick - Pick the sky frames from GUI"
    print,"scalesky - weight the sky frames by the relative median sky value of the red-most channel [outer regions]"
    print,"scaleindchan-weight the sky frames by the relative median sky value for *each* channel [outer regions]"
    print,"skyrad - radius beyond which to estimate sky contribution"
    print,"skylim - limit in cts/s to define 'channels with substantial sky emission'"
    print,"medbox [currently not included] - high-pass filter the data cube after subtracting sky"
    print,"***"
    goto,endofprogram
endif

; setupdir,reducdir=reducdir
datadir='./reduc/prep/'

if ~keyword_set(prefname) then prefname='n'

;charis_subtract_sky,'v1247ori_low.info', /scalesky, suffname='e_flat'
if ~keyword_set(suffname) then begin
    files=FILE_SEARCH(datadir+prefname+'*e_flat.fits', count=flat_count)
    if flat_count eq 0 then suffname='e' else suffname='e_flat'
endif

charis_subtract_sky,pfname,medbox=medbox,prefname=prefname,suffname=suffname,pick=pick,scalesky=scalesky,scaleindchan=scaleindchan,skyrad=skyrad,skylim=skylim,ladder=ladder,verbose=verbose,help=help

endofprogram:
end