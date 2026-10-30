pro charis_pdi_register_cube,pfname,prefname=prefname,$
suffname=suffname,$
polcent_offset=polcent_offset,$
astrogrid=astrogrid,$
method=method,$
rsub=rsub,$
refcube=refcube,$
medbox=medbox,$
ladder=ladder,$
xcorr=xcorr,$
rin_corr=rin_corr,$
rout_corr=rout_corr,$
revise=revise,$
splitpsf=splitpsf,$
keepprev=keepprev,$
nosmooth=nosmooth,$
fwhmlim=fwhmlim,$
smallsteps=smallsteps,$
smask=smask,$
nomask=nomask,$
fittrim=fittrim,$
checkquality=checkquality,$
wait=wait,$
verbose=verbose,$
help=help

if (N_PARAMS() eq 0 or keyword_set(help)) then begin
    print,'charis_pdi_register_cube.pro: Registers left and right pol cubes using charis_register_cube.'
    print,"Typically executed after 'charis_pdi_subtract_sky' and before 'charis_pdi_specphot_cal'."
    print,'Written by K. Lawson (2020)'
    print,''
    print,'**Calling Sequence**'
    print,"charis_pdi_register_cube,pfname,prefname=prefname,Lsuffname=Lsuffname,Rsuffname=Rsuffname,polcent_offset=polcent_offset,"
    print,"astrogrid=astrogrid,method=method,guessoffsets=guessoffsets,rsub=rsub,medbox=medbox,refcube=refcube,ladder=ladder,"
    print,"xcorr=xcorr,splitpsf=splitpsf,keepprev=keepprev,nosmooth=nosmooth,revise=revise,fwhmlim=fwhmlim,"
    print,"smallsteps=smallsteps,smask=smask,checkquality=checkquality,verbose=verbose,help=help"
    print,''
    print,'Example:'
    print,"charis_pdi_register_cube,'abaur_low.info',astrogrid=11.2"
    print,''
    print,"***Keywords***"
    print,'pfname=STRING - The parameter file (e.g. HR8799_low.info)'
    print,"prefname=STRING - The prefix of files to use as input. Defaults to 'n'."
    print,"polcent_offset=ARRAY - Two element array providing the pixel coordinate offset to the center of the right pol from the nominal center. Defaults to [31.758, 16.542]."
    print,"Other keywords as accepted by 'charis_register_cube'."
    goto,endofprogram
endif

if ~keyword_set(polcent_offset) then polcent_offset=[31.758, 16.542]

;* Figuring out what charis_register_cube is going to name our products:
datadir='./reduc/prep/'
reducdir='./reduc/reg/'

;******
;LOGGING
;the command line call

dpp_command=(recall_commands())[0]

;setting up log file
currenttime=systime()
ct=strsplit(currenttime,/extract,' ')
;format is month, day, hms, year
currenttime_log=ct[1]+ct[2]+'_'+ct[4]+'_'+ct[3]
trimmed_time_log=ct[1]+ct[2]+'_'+ct[4]
logfilename=currenttime_log+'.log'
logfilenametrimmed=trimmed_time_log

;now, look to see if logfile exists
checklog_list=file_search('./log/'+logfilenametrimmed+'*',count=count_checklog)
checklog=checklog_list[n_elements(checklog_list)-1]

if count_checklog eq 0 then begin
logname='./log/'+logfilename
openw,11,logname
printf,11,ct[3],dpp_command,format='(a,x,a)'
endif else begin
logname='./'+checklog
openw,11,logname,/append
printf,11,ct[3],dpp_command,format='(a,x,a)'
endelse
;*****

close,11

if ~keyword_set(ladder) then param,'fnum_sat',flist,/get,pfname=pfname else param,'fnum_lad',flist,/get,pfname=pfname
if ~keyword_set(prefname) then prefname='n'
; This could be used on any of four product types: 
; a) PDI flatfielded + skysubbed data
; b) Just PDI flatfielded
; c) Just sky subbed
; d) Neither
if ~keyword_set(suffname) then begin
    test=file_search(datadir+'*flat_skysub.fits')
    if n_elements(test) gt 1 then begin
        suffname='e_flat_skysub'
    endif else begin
        test=file_search(datadir+'*flat.fits')
        if n_elements(test) gt 1 then begin
            suffname='e_flat'
        endif else begin
            test=file_search(datadir+'*skysub.fits')
            if n_elements(test) gt 1 then begin
                suffname='e_skysub'
            endif else begin
                suffname='e'
            endelse
        endelse
    endelse
endif
print,'Using prefname: ',prefname
print,'Using suffname: ',suffname

filenum=nbrlist(flist)
filesin=filelist(filenum,nfiles,prefix=prefname,suffix=suffname)
filesout=filelist(filenum,prefix=prefname,suffix='reg')

;* Names to change defaults to for left and right products:
Lfilesout=filelist(filenum,prefix=prefname,suffix='leftreg')
Rfilesout=filelist(filenum,prefix=prefname,suffix='rightreg')

;testing for parity flip
testcube=readfits(test[0],h0)
htestcube=headfits(test[0])

;****parity flip for post-NBS data
parity_val=sxpar(htestcube,'PARITY',count=parity_count)

if parity_count gt 0 and parity_val eq 1 then begin
polcent_offset[1]*=-1.0
endif

print,'the left offsets are',-1*polcent_offset+[100,100]
print,'the right offsets are',polcent_offset+[100,100]
;stop

if ~keyword_set(astrogrid) then astrogrid=11.2

; Registration for left pol frames:
charis_register_cube,pfname,prefname=prefname,suffname=suffname,method=method,guessoffsets=-1*polcent_offset,$
  rsub=rsub, refcube=refcube, medbox=medbox, ladder=ladder, xcorr=xcorr, rin_corr=rin_corr,rout_corr=rout_corr,revise=revise, splitpsf=splitpsf,$
  keepprev=keepprev, nosmooth=nosmooth, fwhmlim=fwhmlim, smallsteps=smallsteps, smask=smask,fittrim=fittrim,$
  checkquality=checkquality, astrogrid=astrogrid, verbose=verbose, help=help

if ~keyword_set(ladder) then spawn,'mv reduc.log reduc_left.log'
for i=0L,nfiles-1 do spawn,'mv '+reducdir+filesout[i]+' '+reducdir+Lfilesout[i]

; Registration for right pol frames:
charis_register_cube,pfname,prefname=prefname,suffname=suffname,method=method,guessoffsets=polcent_offset,$
  rsub=rsub, refcube=refcube, medbox=medbox, ladder=ladder, xcorr=xcorr, rin_corr=rin_corr,rout_corr=rout_corr,revise=revise, splitpsf=splitpsf,$
  keepprev=keepprev, nosmooth=nosmooth, fwhmlim=fwhmlim, smallsteps=smallsteps, smask=smask,fittrim=fittrim,$
  checkquality=checkquality, astrogrid=astrogrid, verbose=verbose, help=help

if ~keyword_set(ladder) then spawn,'cp reduc.log reduc_right.log'; Just copy instead of moving this one to leave reduc.log for total intensity PSF sub later.
for i=0L,nfiles-1 do spawn,'mv '+reducdir+filesout[i]+' '+reducdir+Rfilesout[i]

; gcntrd does not like it if we replace zeros with nans before performing registration. Instead, we will have
; to load up the headers for the outputs of the above to get the shifts, load up the input images, replaces zeros
; with nans, reapply the shifts, and overwrite the leftreg and rightreg outputs

; In rare cases, skysubbed images manifest with exactly zero values INSIDE the FOV,
; which are propagated as nans if we don't treat them differently.
skysubbed = strmatch(suffname, '*skysub*')
for i=0L,nfiles-1 do begin
    im0 = readfits(datadir+filesin[i], /exten, /silent)
    sz=size(im0, /dim)
    if skysubbed then begin
        h1 = headfits(datadir+filesin[i], ext=1, /silent)
        imraw = readfits(('./data/raw/CRSA'+(sxpar(h1, 'origname'))+'_cube.fits'), ext=1)
        if parity_val eq 1 then outside_fov = where(reverse(imraw,1) eq 0) else outside_fov = where(imraw eq 0)
    endif else begin
        outside_fov = where(im0 eq 0)
    endelse
   
    if ~keyword_set(nomask) then im0[outside_fov] = !values.f_nan

    Lh0 = headfits(reducdir+Lfilesout[i], ext=0, /silent)
    Lh1 = headfits(reducdir+Lfilesout[i], ext=1, /silent)

    Rh0 = headfits(reducdir+Rfilesout[i], ext=0, /silent)
    Rh1 = headfits(reducdir+Rfilesout[i], ext=1, /silent)

    Lreg = fltarr(sz[0],sz[1],sz[2])
    Rreg = fltarr(sz[0],sz[1],sz[2])
    tot = fltarr(sz[0],sz[1],sz[2])
    for iL=0L,sz[2]-1 do begin
        Lshift = strsplit(sxpar(Lh1, 'shift_'+((string(iL)).trim())), /extract)
        Ldx = float(Lshift[0])
        Ldy = float(Lshift[1])
        Lreg[*,*,iL] = shift_sub(im0[*,*,iL],Ldx,Ldy)

        Rshift = strsplit(sxpar(Rh1, 'shift_'+((string(iL)).trim())), /extract)
        Rdx = float(Rshift[0])
        Rdy = float(Rshift[1])
        Rreg[*,*,iL] = shift_sub(im0[*,*,iL],Rdx,Rdy)
        ;tot[*,*,iL]=sqrt(Rreg[*,*,iL]^2.+Lreg[*,*,iL]^2)
        tot[*,*,iL]=Rreg[*,*,iL]+Lreg[*,*,iL]
    endfor
    ;Set any pixels where either L or R isn't defined to nan
    both_finite = finite(Lreg)*finite(Rreg)
    bad = where(~both_finite)

    if ~keyword_set(nomask) then begin
    Lreg[bad] = !values.f_nan
    Rreg[bad] = !values.f_nan
    tot[bad]=!values.f_nan
    endif

    writefits,reducdir+Lfilesout[i],0,Lh0
    writefits,reducdir+Lfilesout[i],Lreg,Lh1,/append

    writefits,reducdir+Rfilesout[i],0,Rh0
    writefits,reducdir+Rfilesout[i],Rreg,Rh1,/append 

    writefits,reducdir+filesout[i],0,Lh0
    writefits,reducdir+filesout[i],tot,Lh1,/append
endfor

close,/all

endofprogram:
end
