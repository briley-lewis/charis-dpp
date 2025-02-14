pro charis_pdi_hwp_match, pfname, use_headers=use_headers, hwpfilein=hwpfilein, prefname=prefname, suffname=suffname, hwpfileout=hwpfileout, logtype=logtype, help=help

if (N_PARAMS() eq 0 or keyword_set(help)) then begin
    print,'charis_pdi_hwp_match.pro: Uses HWP log to assign a HWP angle to each frame,'
    print,'then match to appropriate frames for PDI functionality.'
    print,'Current behavior: checks for the "RET-ANG1" entry in raw headers,'
    print,'which provides the HWP angle (but may be unreliable as of March 2022).'
    print,'If this is found, and a HWP log file is a) not specified by "hwpfilein",'
    print,'and b) not found in the current directory, then the header angles will be used.'
    print,'Otherwise, the HWP angles will be determined from the HWP log file.'
    print,''
    print,"Typically executed after 'charis_pdi_specphot_cal' and before 'charis_pdi'."
    print,'Written by K. Lawson (2021)'
    print,''
    print,'**Calling Sequence**'
    print,"charis_pdi_hwp_match,pfname,hwpfile=hwpfile, prefname=prefname, suffname=suffname, hwpfileout=hwpfileout"
    print,''
    print,'Example:'
    print,"charis_pdi_hwp_match,'abaur_low.info'"
    print,''
    print,"***Keywords***"
    print,'pfname=STRING - The parameter file (e.g. HR8799_low.info)'
    print,'/use_headers - Use HWP angle info from headers. Should be preferred for data after ~Summer 2022. Otherwise should NOT be used.'
    print,"hwpfilein=STRING - path to HWP log from Subaru. Default: searches for filename like 'vampHWPLog_*.csv'"
    print,"prefname=STRING - The prefix of files to use as input. Defaults to 'n'."
    print,"suffname=STRING - The suffix of files to use as input. Defaults to 'leftreg_cal' if present, or 'leftreg' otherwise."
    print,"hwpfileout=STRING - File into which to save HWP angles and matches. Defaults to 'hwp_info.txt'."
    print,"logtype=STRING - Specify the log type of 'hwpfilein' when an unusual naming convention is used: 'vampHWPLog' for a" 
    print,"                 VAMPIRES HWP Log file, or 'hwp_angles' for a csv file containing filenames (col header 'file') and "
    print,"                 corresponding HWP angles (col header 'hwp_ang')."
    goto,endofprogram
endif

reducdir='./reduc/'
datadir=reducdir+'reg/'

param,'fnum_sat',flist,/get,pfname=pfname
if ~keyword_set(prefname) then prefname='n'
if ~keyword_set(suffname) then begin
    test=file_search(datadir+'*leftreg_cal.fits')
    if n_elements(test) gt 1 then begin
        suffname='leftreg_cal'
    endif else begin
        test=file_search(datadir+'*leftreg.fits')
        if n_elements(test) gt 1 then begin
            suffname='leftreg'
        endif else begin
            datadir=reducdir+'prep/'
            suffname='left'
        endelse
    endelse
endif

print,'Using prefname: ',prefname
print,'Using suffname: ',suffname

if ~keyword_set(hwpfileout) then hwpfileout='hwp_info.txt'

filenum=nbrlist(flist)
files=filelist(filenum,nfiles,prefix=prefname,suffix=suffname)

; Check if headers contain the entry for the HWP switch angle (added ~ Feb 2022)
h1test = headfits(datadir+files[0], ext=1, /silent)
hrawtest = headfits(('./data/raw/CRSA'+(sxpar(h1test, 'origname'))+'_cube.fits'), ext=3, /silent)
hwptest = sxpar(hrawtest, 'RET-ANG1', count=hwp_ang_in_headers)

uni_angs = [0.d, 22.5d, 45.d, 67.5d]
if ~keyword_set(use_headers) then begin
    if ~keyword_set(hwpfilein) then begin
        hwp_search = file_search('./', 'vampHWPLog_*.csv', count=fcount)
        if fcount eq 1 then begin
            logtype = 'vampHWPLog'
            hwpfilein = hwp_search[0]
            print,'Using HWP log file:'
            print,hwpfilein
            use_headers = 0
        endif else begin
            if fcount eq 0 then begin
                hwp_search2 = file_search('./', 'hwp_angles*.csv', count=fcount2)
                if fcount2 eq 1 then begin
                    logtype = 'hwp_angles'
                    hwpfilein = hwp_search2[0]
                    print,'Using HWP log file:'
                    print,hwpfilein
                    use_headers = 0
                endif else begin ;fcount2 != 1
                    if fcount2 eq 0 then begin
                        if ~hwp_ang_in_headers then begin
                            print,'No HWP log files found'
                            print,'(e.g. "vampHWPLog_*.csv" or "keywords_log_*.csv"),'
                            print,'nor header entries providing the HWP angles.'
                            print,'Please add a HWP log file to your working directory,'
                            print,'or explicitly indicate the file to use by setting'
                            print,'the keyword "hwpfilein".'
                            goto,endofprogram
                        endif else begin
                            print,'No HWP log files found'
                            print,'(e.g. "vampHWPLog_*.csv" or "hwp_angles*.csv"),'
                            print,'using HWP angles from header entries.'
                            print,'Warning: header entries should *only* be used for data collected after ~Summer 2022. Prior to that,'
                            print,'any header info is likely incorrect, and an alternate method should be used for HWP matching.'
                            use_headers = 1
                            logtype='headers'
                        endelse
                    endif else begin ; fcount2 > 1
                        print,'Multiple HWP log files found:'
                        print, hwp_search2
                        print,'Please explicitly indicate the file to use by setting'
                        print,'the keyword "hwpfilein" to the file you wish to use.'
                        goto,endofprogram
                    endelse
                endelse
            endif else begin ; fcount > 1
                print,'Multiple HWP log files found:'
                print, hwp_search
                print,'Please explicitly indicate the file to use by setting'
                print,'the keyword "hwpfilein" to the file you wish to use.'
                goto,endofprogram
            endelse
        endelse
    endif else begin
        use_headers = 0
        if ~keyword_set(logtype) then begin
            if strpos(hwpfilein, 'vampHWPLog') ne -1 then begin
                logtype = 'vampHWPLog'
            endif else begin
                if strpos(hwpfilein, 'hwp_angles') ne -1 then begin
                    logtype = 'hwp_angles'
                endif else begin
                    print,'HWP log type not recognized from filename!'
                    print,'Please specify logtype="vampHWPLog" for a HWP'
                    print,'log produced by VAMPIRES, or logtype="hwp_angles"'
                    print,'for a csv log simply containing filenames (col header "file")
                    print,'with corresponding HWP angles (col header "hwp_angle").'
                    goto,endofprogram
                endelse
            endelse
        endif
    endelse
endif else begin
    if hwp_ang_in_headers then begin
        print,'Using HWP angles from header entries. Warning: Should *only* be used for data collected after ~Summer 2022. Prior to that,'
        print,'any header info is likely incorrect, and an alternate method should be used for HWP matching.'
        logtype='headers'
    endif else begin
        print,'No HWP info found in headers! /use_headers will not work for this data. Use a HWP log instead."
        goto,endofprogram
    endelse
endelse

if logtype eq 'vampHWPLog' or use_headers then begin
    if ~use_headers then begin
        readcol,hwpfilein,angs,utstart_string,utend_string,target,format='(f,a,a,a)'
        ; There are at least two utilized formats for PDI logs: 
        ; 1) start and end time columns are formatted as: hh:mm:ss.sss
        ; 2) yyyymmddThhmmss (where the T is literally the letter T)
        utchwp=fltarr(n_elements(angs))
        for i=0L,n_elements(angs)-1 do begin
            pos = STRPOS(utend_string[i], ':')
            if pos eq -1 then begin
                datetime = strsplit(utstart_string[i],/extract,'T')
                hour = strmid(datetime[1],0,2)
                minute = strmid(datetime[1],2,2)
                second = strmid(datetime[1],4,2)
                utchwp[i]=ten(hour,minute,second)
            endif else begin
                time=strsplit(utstart_string[i],/extract,':')   
                utchwp[i]=ten(time[0],time[1],time[2])
            endelse
        endfor
    endif

    hwp_ang=fltarr(nfiles)
    hwp_pos=intarr(nfiles)
    utcexp=fltarr(nfiles)
    utcexp_end=fltarr(nfiles)

    for i=0L,nfiles-1 do begin
        h=headfits(datadir+files[i], ext=1, /silent)
        exptime=sxpar(h, 'EXPTIME')

        origname = sxpar(h,'ORIGNAME')
        origfile = './data/raw/CRSA'+origname+'_cube.fits'
        hraw = headfits(('./data/raw/CRSA'+(sxpar(h, 'origname'))+'_cube.fits'), ext=3, /silent) ;This is the header of the raw image
        dit = (sxpar(hraw,'MJD-END') - sxpar(hraw,'MJD-STR'))*86400. ; In seconds
        timestring=sxpar(h,'UTC-TIME')
        time=strsplit(timestring,/extract,':')
        time[2] -= 0.5*dit

        utcexp[i]=ten(time[0],time[1],time[2])
        utcexp_end[i]=ten(time[0],time[1],time[2]+exptime)
        if use_headers then begin
            hwp_ang_i = sxpar(hraw, 'RET-ANG1')
            if hwp_ang_i eq -1 then hwp_ang_i = 0.d
            hwp_ang[i] = hwp_ang_i
        endif else begin
            dts = abs(utchwp-utcexp[i])
            min_dt = min(dts, argmin) ;Get index of nearest hwp move end (in time) to utcexp[i]
            if utchwp[argmin] lt utcexp[i] then hwpind = argmin ;If corresponding time is before utcexp[i], then that index gives the matching hwp pos
            if utchwp[argmin] gt utcexp[i] then hwpind = argmin-1 ;If corresponding time is after utcexp[i], then we want entry for one index prior
            hwp_ang[i] = angs[hwpind]
        endelse
        min_diff = min(abs(hwp_ang[i]-uni_angs), angpos)
        hwp_pos[i] = angpos
    endfor
endif else begin ;logtype eq 'hwp_angles'
    hwp_ang=fltarr(nfiles)
    hwp_pos=intarr(nfiles)
    utcexp = fltarr(nfiles)

    log = read_csv(hwpfilein, header = thead)
    thead = STRLOWCASE(thead)
    fileind = where(thead eq 'file')
    hwpangind = where(thead eq 'hwp_ang')
    if fileind eq -1 or hwpangind eq -1 then begin
        print,"Please ensure that your log file has column headers labeled:"
        print,'"file", indicating the name of the unprocessed file (CRSA*),'
        print,'and "hwp_ang", indicating the HWP angle for each file.'
    endif

    hwpanglog = log.(hwpangind)
    filenames = log.(fileind)
    filenames = filenames.replace('.fits', '')
    filenames = filenames.replace('_cube', '')
    filenames = filenames.replace('CRSA', '')
    for i=0L,nfiles-1 do begin
        h=headfits(datadir+files[i], ext=1, /silent)
        origname = sxpar(h,'ORIGNAME')
        logind = where(filenames eq origname)
        if logind eq -1 then begin
            hwp_ang[i] = -1
            hwp_pos[i] = -1
        endif else begin
            hwp_ang[i] = hwpanglog[logind]
            hwp_pos[i] = where(uni_angs eq hwp_ang[i])
        endelse
        timestring=sxpar(h,'UTC-TIME')
        time=strsplit(timestring,/extract,':')
        utcexp[i]=ten(time[0],time[1],time[2])
    endfor
endelse

nangs = n_elements(uni_angs)
ang_counts = intarr(nangs)
for i=0,nangs-1 do ang_counts[i] = total((hwp_pos eq i))

ncycles = min(ang_counts, min_occ_pos)
min_occ_msk = where(hwp_pos eq min_occ_pos)
inds = indgen(nfiles)
hwp_cyc1 = intarr(nfiles)
hwp_cyc1[*] = -1
used_inds = []
for i=0L,ncycles-1 do begin
    ind_min_occ = inds[min_occ_msk[i]]
    used_inds=[used_inds, ind_min_occ]
    hwp_cyc1[ind_min_occ] = i
    dts = abs(utcexp[ind_min_occ] - utcexp)
    for j=0,nangs-1 do begin
        if j ne min_occ_pos then begin
            ang_inds = inds[where(hwp_pos eq j, ang_count)]
            ang_inds_notused = []; I hope there's a better way to do this in IDL...
            for k=0,ang_count-1 do begin
                used = where(used_inds eq ang_inds[k], used_count)
                if used_count eq 0 then ang_inds_notused=[ang_inds_notused, ang_inds[k]]
            endfor
            if n_elements(ang_inds_notused) gt 0 then begin
                min_dt = min(dts[ang_inds_notused], dt_argmin)
                nearest_ind = ang_inds_notused[dt_argmin]
                used_inds = [used_inds, nearest_ind]
                hwp_cyc1[nearest_ind] = i
            endif
        endif
    endfor
endfor 

hwp_cyc2 = intarr(nfiles)
hwp_cyc2[*] = -1
used_inds = []
for ind=0L,ncycles-1 do begin
    i = ncycles-1-ind ; reverse order...
    ind_min_occ = inds[min_occ_msk[i]]
    used_inds=[used_inds, ind_min_occ]
    hwp_cyc2[ind_min_occ] = i
    dts = abs(utcexp[ind_min_occ] - utcexp)
    for j=0,nangs-1 do begin
        if j ne min_occ_pos then begin
            ang_inds = inds[where(hwp_pos eq j, ang_count)]
            ang_inds_notused = []; I hope there's a better way to do this in IDL...
            for k=0,ang_count-1 do begin
                used = where(used_inds eq ang_inds[k], used_count)
                if used_count eq 0 then ang_inds_notused=[ang_inds_notused, ang_inds[k]]
            endfor
            if n_elements(ang_inds_notused) gt 0 then begin
                min_dt = min(dts[ang_inds_notused], dt_argmin)
                nearest_ind = ang_inds_notused[dt_argmin]
                used_inds = [used_inds, nearest_ind]
                hwp_cyc2[nearest_ind] = i
            endif
        endif
    endfor
endfor 

tranges1 = []
tranges2 = []
for i=0L,ncycles-1 do begin
    times1 = utcexp[where(hwp_cyc1 eq i)]
    trange1 = max(times1) - min(times1)
    tranges1 = [tranges1, trange1]

    times2 = utcexp[where(hwp_cyc2 eq i)]
    trange2 = max(times2) - min(times2)
    tranges2 = [tranges2, trange2]
endfor

if mean(tranges1) lt mean(tranges2) then hwp_cyc = hwp_cyc1 else hwp_cyc = hwp_cyc2

openw,lunit,hwpfileout,/get_lun
printf, lunit, 'FILE, ', 'HWP_POS, ', 'HWP_ANG, ', 'CYCLE, '
for i=0,nfiles-1 do begin
    printf, lunit, filenum[i], hwp_pos[i], hwp_ang[i], hwp_cyc[i]
endfor
free_lun,lunit
endofprogram:
end