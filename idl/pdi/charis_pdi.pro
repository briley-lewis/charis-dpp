pro charis_pdi, pfname,$
 prefname=prefname, Lsuffname=Lsuffname, Rsuffname=Rsuffname, hwpinfo_name=hwpinfo_name, outfile=outfile,$
 meanadd=meanadd, collapse_type=collapse_type,$
 angoffset=angoffset, save_single_sums=save_single_sums, single_sum_zeros=single_sum_zeros,$
 no_mm_ip_corr=no_mm_ip_corr, no_mm_crosstalk_corr=no_mm_crosstalk_corr,$
 first_order_ip_corr=first_order_ip_corr, foip_rlims=foip_rlims,$
 bg_sub=bg_sub, bg_rlims = bg_rlims,$
 unres_pol_sub=unres_pol_sub, unres_pol_rlims=unres_pol_rlims, avg_unres_pol_coeffs=avg_unres_pol_coeffs, make_unres_pol_plot=make_unres_pol_plot,$
 collapse_cycles=collapse_cycles, save_full_seq=save_full_seq, save_derot_seq=save_derot_seq, save_uncorrected_derot=save_uncorrected_derot, no_adisync=no_adisync, $ 
 hwp_inds = hwp_inds, extra_pol_angle_corr=extra_pol_angle_corr, help=help

if (N_PARAMS() eq 0 or keyword_set(help)) then begin
    print,'charis_pdi.pro: Carries out single+double summing/differencing of PDI sequences and correction for instrumental effects, outputting Stokes parameter image cubes and wavelength-collapsed images.'
    print,"Typically executed after 'charis_pdi_hwp_match'."
    print,'Written by K. Lawson (2022)'
    print,''
    print,'**Calling Sequence**'
    print,'charis_pdi, pfname, prefname=prefname, Lsuffname=Lsuffname, Rsuffname=Rsuffname, hwpinfo_name=hwpinfo_name, outfile=outfile, meanadd=meanadd, collapse_type=collapse_type, angoffset=angoffset, single_sum_zeros=single_sum_zeros, no_mm_ip_corr=no_mm_ip_corr, no_mm_crosstalk_corr=no_mm_crosstalk_corr, first_order_ip_corr=first_order_ip_corr, foip_rlims=foip_rlims, bg_sub=bg_sub, bg_rlims = bg_rlims, unres_pol_sub=unres_pol_sub, unres_pol_rlims=unres_pol_rlims, avg_unres_pol_coeffs=avg_unres_pol_coeffs, collapse_cycles=collapse_cycles, save_full_seq=save_full_seq, save_derot_seq=save_derot_seq, save_uncorrected_derot=save_uncorrected_derot, no_adisync=no_adisync, hwp_inds = hwp_inds, extra_pol_angle_corr=extra_pol_angle_corr, help=help'
    print,''
    print,'Example:'
    print,"charis_pdi, 'ab_aur_low.info', /meanadd, /unres_pol_sub, /avg_unres_pol_coeffs"
    print,''
    print,"***Keywords***"
    print,"pfname=STRING - The parameter file (e.g. 'HR8799_low.info')"
    print,"prefname=STRING - The prefix of files to use as input. Defaults to 'n'."
    print,"Lsuffname=STRING - The suffix of left pol files to use as input. Defaults to 'leftreg_cal', or 'leftreg' if no images are found for 'leftreg_cal'."
    print,"Rsuffname=STRING - The suffix of right pol files to use as input. Defaults to 'rightreg_cal', or 'rightreg' if no images are found for 'rightreg_cal'."
    print,"hwpinfo_name=STRING - The path for the half-wave plate info file produced by charis_pdi_hwp_match. Defaults to 'hwp_info.txt'."
    print,"outfile=STRING - The prefix to use in naming final sequence-averaged PDI images and image cubes. Defaults to 'pdi'."
    print,"/meanadd - Use an outlier resistant mean rather than the median to combine exposures to final image cubes."
    print,"collapse_type=STRING{'median'|'mean'|'rmean'} - The method to use for calculating wavelength-collapsed products. 'median', 'mean', 'rmean' use the median, mean, or outlier resistant mean, respectively. If /mean_add is set, defaults to 'mean', else 'median'."
    print,"angoffset=FLOAT - Angle offset in degrees of the nominal image PA from true north. Defaults to the value for 'angoffset' from charis_get_constant."
    print,"/save_single_sums - Save single sum image cubes for use with total intensity PSF subtraction. Output: ./reduc/reg/n????{suffname}.fits"
    print,"/single_sum_zeros - Apply zeros instead of nans to pixels beyond the FOV for the output single sum image cubes."
    print,"/no_mm_ip_corr - Disables Mueller-matrix-based instrumental polarization corrections."
    print,"/no_mm_crosstalk_corr - Disables Mueller-matrix-based crosstalk and efficiency corrections."
    print,"/first_order_ip_corr - Use first order IP correction in place of MM corrections in atypical scenarios (e.g. significant optical path changes that invalidate the current MM model). /no_mm_ip_corr and /no_mm_crosstalk_corr must be set for this to work."
    print,"foip_rlims=ARRAY - A two element array giving the inner and outer pixel radii for the annulus in which to estimate instrumental pol when /first_order_ip_corr is set. Defaults to [10,25]."
    print,"/bg_sub - Perform background subtaction on Q, U, IQ, and IU. Not generally useful given CHARIS's narrow FOV, but provided as an option just in case."
    print,"bg_rlims=ARRAY - A two element array giving the inner and outer pixel radii for the annulus in which to compute the background signal when /bg_sub is set. Defaults to [55,60]."
    print,"/unres_pol_sub - Subtract unresolved polarization (i.e., originating from intervening dust or from unresolved disk signal)."
    print,"unres_pol_rlims=ARRAY - A two element array giving the inner and outer pixel radii for the annulus in which to compute the unresolved polarization when /unres_pol_sub is set. Defaults to [0,6], which uses the area beneath the standard 113mas coronagraph to estimate unres pol. Note: this should be set to avoid any significant resolved polarization and must be tuned for each target."
    print,"/avg_unres_pol_coeffs - When /unres_pol_sub is set: for each wavelength, use the median unresolved polarization across all cycles to perform unresolved polarization subtraction."
    print,"/make_unres_pol_plot - Create a figure showing calculated unres pol values per cycle and per wavelength, as well as average PA and altitude per Q and U cycle."
    print,"/collapse_cycles - Output wavelength-collapsed products for each corrected and derotated HWP cycle."
    print,"/save_full_seq - Save the Q,U,IQ,IU cubes resulting from double sums and differences BEFORE derotation and correction. Output: ./reduc/proc/{outfile}_c????_q.fits, etc."
    print,"/save_derot_seq - Save the derotated Q,U,IQ,IU cubes for each cycle after any corrections have been applied. Output: ./reduc/proc/{outfile}_c????_q_derot.fits, etc."
    print,"/save_uncorrected_derot - Save the derotated Q and U cubes for each cycle BEFORE any corrections have been applied. Output: ./reduc/proc/{outfile}_c????_q_derot_uncorr.fits, etc."
    print,"/no_adisync - Enable only if ADI synchronization was NOT used for your observations (this switch generally should not be used). If unsure: run the procedure with /collapse_cycles and check the ./reduc/proc/*c????_q_derot.fits images. If the pattern in Q is ~aligned in the derotated sequence, ADI sync was used and you do NOT need this option."
    print,"hwp_inds=ARRAY - Array of four integers indicating which values of the column 'HWP_POS' in the HWP info file correspond to Q+, U+, Q-, and U- respectively. Defaults to [0,1,2,3]."
    print,"extra_pol_angle_corr=FLOAT - Any additional offset (in degrees) by which to correct the polarization axes. Note: we already correct for the true north PA offset ('angoffset')."
    goto,endofprogram
endif

command = ((recall_commands())[0]).replace("'", '"')

reducdir='./reduc/'
datadir=reducdir+'reg/'
procdir=reducdir+'proc/'

if ~keyword_set(angoffset) then angoffset = charis_get_constant(name='angoffset') ;nominally 2.2 deg

if ~keyword_set(hwpinfo_name) then hwpinfo_name = 'hwp_info.txt'

readcol,hwpinfo_name,name,hwp_pos,hwp_ang,hwp_cycle,/silent,skipline=1,format="(A,I,F,I)"

if ~keyword_set(prefname) then prefname = 'n'

if ~keyword_set(Lsuffname) then begin
    files=FILE_SEARCH(datadir+prefname+'*leftreg_cal.fits', count=cal_count)
    if cal_count eq 0 then Lsuffname='leftreg' else Lsuffname='leftreg_cal'
endif

if ~keyword_set(Rsuffname) then begin
    files=FILE_SEARCH(datadir+prefname+'*rightreg_cal.fits', count=cal_count)
    if cal_count eq 0 then Rsuffname='rightreg' else Rsuffname='rightreg_cal'
endif

if ~keyword_set(outfile) then outfile='pdi'

if ~keyword_set(collapse_type) then begin
    if keyword_set(meanadd) then collapse_type='mean' else collapse_type='median'
endif

param,'fnum_sat',flist,/get,pfname=pfname

filenum_pol= nbrlist(flist)
files_polleft = datadir+filelist(filenum_pol,nfiles_polpair,prefix=prefname,suffix=Lsuffname)
files_polright = datadir+filelist(filenum_pol,nfiles_polpair,prefix=prefname,suffix=Rsuffname)

h0test = headfits(files_polleft[0],/silent)
imtest = readfits(files_polleft[0],h1test,/exten,/silent)
sz = size(imtest, /dim)
dimx = sz[0]
dimy = sz[1]

allpa = fltarr(nfiles_polpair)
for i=0,nfiles_polpair-1 do allpa[i] = sxpar((headfits(files_polleft[i], /ext)), 'parang')

get_charis_wvlh,h0test,lambda

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Get constant (per wavelength) values for Mueller matrix
delta_hwp = hwp_model(lambda, epsilon_hwp=epsilon_hwp)
delta_der = derotator_model(lambda, epsilon_der=epsilon_der)
delta_tel = m3_model(lambda, epsilon_m3=epsilon_tel)
epsilon_bs = 1.d
dtheta_hwp = -0.002d
dtheta_der = -0.0118d
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Remove old per-cycle products of same 'outfile' prefix. This avoids possible confusion; e.g., if you're using /collapse_cycles
; and want to rerun a reduction with the same name but with a different hwp_info file using fewer cycles.
previous_proc_files=FILE_SEARCH(procdir+outfile+'_c????_*.fits', count=fcount)
if fcount gt 0 then spawn,'rm '+procdir+outfile+'_c????_*.fits'

; Make filelists for input/output
valid_cycles = where(hwp_cycle ne -1)
cycles = (hwp_cycle[valid_cycles])[uniq(hwp_cycle[valid_cycles], sort(hwp_cycle[valid_cycles]))]
ncycles = n_elements(cycles) ; Exclude the unused cycle=-1 frames (if any)

IQout = procdir+filelist(cycles, ncyc, prefix=(outfile+'_c'), suffix='_iq'); Saving tot intensity products to 'reduc/reg/' so they can be used with PSFsub programs without alteration
IUout = procdir+filelist(cycles, ncyc, prefix=(outfile+'_c'), suffix='_iu')

Qout = procdir+filelist(cycles, ncyc, prefix=(outfile+'_c'), suffix='_q')
Uout = procdir+filelist(cycles, ncyc, prefix=(outfile+'_c'), suffix='_u')

; Prepare empty arrays for later use
q_dpas = fltarr(ncycles)
u_dpas = fltarr(ncycles)

qhcube = fltarr(sz[0],sz[1],sz[2],ncycles)
uhcube = fltarr(sz[0],sz[1],sz[2],ncycles)
iqhcube = fltarr(sz[0],sz[1],sz[2],ncycles)
iuhcube = fltarr(sz[0],sz[1],sz[2],ncycles)

q_ixs = fltarr(sz[2], ncycles)
q_qxs = fltarr(sz[2], ncycles)
q_uxs = fltarr(sz[2], ncycles)
u_ixs = fltarr(sz[2], ncycles)
u_qxs = fltarr(sz[2], ncycles)
u_uxs = fltarr(sz[2], ncycles)

; Setup for unres pol sub
if keyword_set(unres_pol_sub) then begin
    distarr = shift(dist(sz[0]), (sz[0])/2, (sz[0])/2)
    qstar = fltarr(sz[2], ncycles)
    ustar = fltarr(sz[2], ncycles)
    iqstar = fltarr(sz[2], ncycles)
    iustar = fltarr(sz[2], ncycles)
    if ~keyword_set(unres_pol_rlims) then unres_pol_rlims = [0, 6]
    unres_annulus = where((distarr ge unres_pol_rlims[0]) and (distarr le unres_pol_rlims[1]))

    telluric_channels = [5,6,14,21]
    qrel_fit = fltarr(sz[2], ncycles)
    urel_fit = fltarr(sz[2], ncycles)
    qalt = fltarr(ncycles)
    ualt = fltarr(ncycles)
    qparang = fltarr(ncycles)
    uparang = fltarr(ncycles)
endif

; Setup for background sub
if keyword_set(bg_sub) then begin
    if ~keyword_set(bg_rlims) then bg_rlims=[55,60]
    if ~keyword_set(distarr) then distarr = shift(dist(sz[0]), (sz[0])/2, (sz[0])/2)
    bg_annulus = where((distarr ge bg_rlims[0]) and (distarr le bg_rlims[1]))
endif

if keyword_set(first_order_ip_corr) then begin
    if keyword_set(no_mm_ip_corr) and keyword_set(no_mm_crosstalk_corr) then begin
        if ~keyword_set(foip_rlims) then foip_rlims=[10,25]
        if ~keyword_set(distarr) then distarr = shift(dist(sz[0]), (sz[0])/2, (sz[0])/2)
        foip_annulus = where((distarr ge foip_rlims[0]) and (distarr le foip_rlims[1]))
    endif else begin
        print,"Warning: to use first order instrumental pol. correction (/first_order_ip_corr), you must also set /no_mm_crosstalk_corr and /no_mm_ip_corr."
        print,"Proceeding without first order correction!"
        print,"Please restart the procedure with /no_mm_crosstalk_corr and /no_mm_ip_corr if you wish to use the first order correction instead."
        WAIT, 3. ; So this message doesn't get lost in the sea of imports.
        first_order_ip_corr = 0
    endelse
endif

if ~keyword_set(hwp_inds) then begin
    hwp_inds = [0,1,2,3]
endif

if keyword_set(save_single_sums) then begin
    ; Separately computing all single sums and writing to ./reduc/reg with standard total intensity filenames:
    for Ti=0,nfiles_polpair-1 do begin
        h0L = headfits(files_polleft[Ti], ext=0, /silent)
        leftcube= readfits(files_polleft[Ti],h1L,/exten,/silent)
        rightcube= readfits(files_polright[Ti],h1R,/exten,/silent)
        imcube = leftcube + rightcube
        if keyword_set(single_sum_zeros) then imcube[(where(~finite(imcube)))] = 0.
        file_singlesum = files_polleft[Ti].replace(Lsuffname, (Lsuffname.replace('left', '')))
        writefits,file_singlesum,0,h0L
        writefits,file_singlesum,imcube,h1L,/append
    endfor
endif

;;;;;;;;;;;;;;;    Main PDI loop    ;;;;;;;;;;;;;;;
for ci=0,ncycles-1 do begin
    cnum = cycles[ci]
    cycle_inds = where(hwp_cycle eq cnum, Ncyclefiles)

    Qp_ind= cycle_inds[where(hwp_pos[cycle_inds] eq hwp_inds[0])]
    Up_ind= cycle_inds[where(hwp_pos[cycle_inds] eq hwp_inds[1])]
    Qm_ind= cycle_inds[where(hwp_pos[cycle_inds] eq hwp_inds[2])]
    Um_ind= cycle_inds[where(hwp_pos[cycle_inds] eq hwp_inds[3])]

    if (n_elements(Qp_ind)+n_elements(Up_ind)+n_elements(Qm_ind)+n_elements(Um_ind)) ne 4 then begin
        print,'For cycle ',(string(cnum)).trim(),' more/fewer than 4 images are indicated in ',hwpinfo_name
        print,'Please ensure that each cycle index (except -1) appears exactly 4 times!'
        break
    endif
    cycle_pas = [allpa[Qp_ind], allpa[Up_ind], allpa[Qm_ind], allpa[Um_ind]]
 
    left_cycle = fltarr(sz[0],sz[1],sz[2],4)
    right_cycle = fltarr(sz[0],sz[1],sz[2],4)

    left_cycle[*,*,*,0] = readfits(files_polleft[Qp_ind],h1_qpL,/exten,/silent) ; IQ+L
    left_cycle[*,*,*,1] = readfits(files_polleft[Up_ind],h1_upL,/exten,/silent) ; IU+L
    left_cycle[*,*,*,2] = readfits(files_polleft[Qm_ind],h1_qmL,/exten,/silent) ; IQ-L
    left_cycle[*,*,*,3] = readfits(files_polleft[Um_ind],h1_umL,/exten,/silent) ; IU-L

    right_cycle[*,*,*,0] = readfits(files_polright[Qp_ind],h1_qpR,/exten,/silent) ; IQ+R
    right_cycle[*,*,*,1] = readfits(files_polright[Up_ind],h1_upR,/exten,/silent) ; IU+R
    right_cycle[*,*,*,2] = readfits(files_polright[Qm_ind],h1_qmR,/exten,/silent) ; IQ-R
    right_cycle[*,*,*,3] = readfits(files_polright[Um_ind],h1_umR,/exten,/silent) ; IU-R

    ; Single sums
    IQp = left_cycle[*,*,*,0] + right_cycle[*,*,*,0]
    IUp = left_cycle[*,*,*,1] + right_cycle[*,*,*,1]
    IQm = left_cycle[*,*,*,2] + right_cycle[*,*,*,2]
    IUm = left_cycle[*,*,*,3] + right_cycle[*,*,*,3]

    ; Double sums
    IQ = 0.5*(IQp + IQm)
    IU = 0.5*(IUp + IUm)

    ; Single diffs
    Qp = right_cycle[*,*,*,0] - left_cycle[*,*,*,0]
    Up = right_cycle[*,*,*,1] - left_cycle[*,*,*,1]
    Qm = right_cycle[*,*,*,2] - left_cycle[*,*,*,2]
    Um = right_cycle[*,*,*,3] - left_cycle[*,*,*,3]

    ; Double diffs
    PQ = 0.5*(Qp-Qm)
    PU = 0.5*(Up-Um)

    if keyword_set(no_mm_ip_corr) and keyword_set(no_mm_crosstalk_corr) then goto,skip_mm_calcs
    xqp = get_charis_mueller_matrix(h1_qpL, (hwp_ang[Qp_ind]), epsilon_bs, epsilon_der, delta_der, epsilon_hwp, delta_hwp, epsilon_tel, delta_tel, dtheta_der, dtheta_hwp, no_adisync=no_adisync)
    xup = get_charis_mueller_matrix(h1_upL, (hwp_ang[Up_ind]), epsilon_bs, epsilon_der, delta_der, epsilon_hwp, delta_hwp, epsilon_tel, delta_tel, dtheta_der, dtheta_hwp, no_adisync=no_adisync)
    xqm = get_charis_mueller_matrix(h1_qmL, (hwp_ang[Qm_ind]), epsilon_bs, epsilon_der, delta_der, epsilon_hwp, delta_hwp, epsilon_tel, delta_tel, dtheta_der, dtheta_hwp, no_adisync=no_adisync)
    xum = get_charis_mueller_matrix(h1_umL, (hwp_ang[Um_ind]), epsilon_bs, epsilon_der, delta_der, epsilon_hwp, delta_hwp, epsilon_tel, delta_tel, dtheta_der, dtheta_hwp, no_adisync=no_adisync)
    
    q_mr = 0.5*(xqp - xqm) ; Q double-diff Mueller row-vector, etc.
    q_ixs[*,ci] = q_mr[*, 0]
    q_qxs[*,ci] = q_mr[*, 1]
    q_uxs[*,ci] = q_mr[*, 2]

    u_mr = 0.5*(xup - xum)
    u_ixs[*,ci] = u_mr[*, 0]
    u_qxs[*,ci] = u_mr[*, 1]
    u_uxs[*,ci] = u_mr[*, 2]

    skip_mm_calcs:
    qpa = compute_average_angle([cycle_pas[0],cycle_pas[2]])
    upa = compute_average_angle([cycle_pas[1],cycle_pas[3]])

    h0 = headfits(files_polleft[Up_ind], ext=0, /silent)
    h1 = h1_upL

    q_dpas[ci] = qpa - allpa[Up_ind] ; Note: Up_ind here is not a typo -- we're using the WCS info from h1_upL, so we need angle offsets relative to this
    u_dpas[ci] = upa - allpa[Up_ind]
    
    ; Saving IQ, Q products...
    orignames=sxpar(h1_qpL,'origname')+', '+sxpar(h1_qmL,'origname')
    sxaddpar,h1,'cycle',orignames,'origname entries for frames in this cycle'
    writefits,IQout[ci],0,h0
    writefits,IQout[ci],IQ,h1,/append
    writefits,Qout[ci],0,h0
    writefits,Qout[ci],PQ,h1,/append

    ; Saving IU, U products...
    orignames=sxpar(h1_upL,'origname')+', '+sxpar(h1_umL,'origname')
    sxaddpar,h1,'cycle',orignames,'origname entries for frames in this cycle'
    writefits,IUout[ci],0,h0
    writefits,IUout[ci],IU,h1,/append
    writefits,Uout[ci],0,h0
    writefits,Uout[ci],PU,h1,/append
endfor

phi = reverse(transpose(angarr(sz[0])))
; Derotation, crosstalk corr, unres pol sub, and axis correction for no ADI sync
for ci=0,ncycles-1 do begin
    cnum = cycles[ci]

    qh0 = headfits(Qout[ci], ext=0, /silent)
    qcube = readfits(Qout[ci], qh1, /exten, /silent)
    qangoffset = angoffset + q_dpas[ci]
    charis_northup,qcube,qh0,qh1,angoffset=qangoffset,outputrot=qoutputrot

    uh0 = headfits(Uout[ci], ext=0, /silent)
    ucube = readfits(Uout[ci], uh1, /exten, /silent)
    uangoffset = angoffset + u_dpas[ci]
    charis_northup,ucube,uh0,uh1,angoffset=uangoffset,outputrot=uoutputrot

    ; loading and aligning total intensity results...
    iqh0 = headfits(IQout[ci], ext=0, /silent)
    iqcube = readfits(IQout[ci], iqh1, /exten, /silent)
    charis_northup,iqcube,iqh0,iqh1,angoffset=qangoffset,outputrot=iqoutputrot

    iuh0 = headfits(IUout[ci], ext=0, /silent)
    iucube = readfits(IUout[ci], iuh1, /exten, /silent)
    charis_northup,iucube,iuh0,iuh1,angoffset=uangoffset,outputrot=iuoutputrot

    if ~keyword_set(save_full_seq) then begin
        spawn,'rm ' + Qout[ci]
        spawn,'rm ' + Uout[ci]
        spawn,'rm ' + IQout[ci]
        spawn,'rm ' + IUout[ci]
    endif

    if keyword_set(save_uncorrected_derot) then begin
        Qout_derot = Qout[ci].replace('.fits', '_derot_uncorr.fits')
        writefits,Qout_derot,0,qh0
        writefits,Qout_derot,qcube,qh1,/append
        
        Uout_derot = Uout[ci].replace('.fits', '_derot_uncorr.fits')
        writefits,Uout_derot,0,uh0
        writefits,Uout_derot,ucube,uh1,/append
    endif
    
    if keyword_set(first_order_ip_corr) then begin
        for Li=0,sz[2]-1 do begin
            cq = median((qcube[*,*,Li]/iqcube[*,*,Li])[foip_annulus], /even)
            cu = median((ucube[*,*,Li]/iucube[*,*,Li])[foip_annulus], /even)
            qcube[*,*,Li] = (qcube[*,*,Li]) - (iqcube[*,*,Li])*cq
            ucube[*,*,Li] = (ucube[*,*,Li]) - (iucube[*,*,Li])*cu
        endfor
    endif

    ; IP CORRECTION:
    if ~keyword_set(no_mm_ip_corr) then begin
        for Li=0,sz[2]-1 do begin
            qcube[*,*,Li] = qcube[*,*,Li] - (iqcube[*,*,Li] * q_ixs[Li,ci])
            ucube[*,*,Li] = ucube[*,*,Li] - (iucube[*,*,Li] * u_ixs[Li,ci])
        endfor
    endif

    ; CROSSTALK CORRECTION:
    if ~keyword_set(no_mm_crosstalk_corr) then begin
        nx = long(sz[0])
        ny = long(sz[1])
        nL = long(sz[2])
        qinc = fltarr(nx, ny, nL)
        uinc = fltarr(nx, ny, nL)
        for Li=0,nL-1 do begin
            A = [[q_qxs[Li,ci], q_uxs[Li,ci]],$
                 [u_qxs[Li,ci], u_uxs[Li,ci]]]
            qrow = reform(qcube[*,*,Li], (nx*ny))
            urow = reform(ucube[*,*,Li], (nx*ny))
            B = [[qrow],$
                 [urow]]
            xrow = la_linear_equation(A, B)
            qinc[*,*,Li] = reform(xrow[*,0], nx, ny)
            uinc[*,*,Li] = reform(xrow[*,1], nx, ny)
        endfor
        qcube = qinc
        ucube = uinc
    endif
    ; CROSSTALK CORRECTION END

    ; Adjusting polarization axes so -Q is aligned with north-south and +Q is aligned with east-west (as necessary)
    pol_angle_corr = 0
    if keyword_set(no_adisync) then pol_angle_corr = (compute_average_angle([qoutputrot, uoutputrot])) else pol_angle_corr = angoffset
    if keyword_set(extra_pol_angle_corr) then pol_angle_corr += extra_pol_angle_corr
    if pol_angle_corr ne 0 then qcube = adjust_pol_axis(qcube, ucube, pol_angle_corr, phi, ucube_out=ucube)

    if keyword_set(bg_sub) then begin
        for Li=0,sz[2]-1 do begin
            qcube[*,*,Li] -= median((qcube[*,*,Li])[bg_annulus], /even)
            ucube[*,*,Li] -= median((ucube[*,*,Li])[bg_annulus], /even)
            iqcube[*,*,Li] -= median((iqcube[*,*,Li])[bg_annulus], /even)
            iucube[*,*,Li] -= median((iucube[*,*,Li])[bg_annulus], /even)
        endfor
    endif

    if keyword_set(unres_pol_sub) then begin
        qalt[ci] = sxpar(qh0, 'altitude')
        ualt[ci] = sxpar(uh0, 'altitude')
        qparang[ci] = (((sxpar(qh0, 'parang') - (angoffset + q_dpas[ci]))+180d) mod 360d) - 180d
        uparang[ci] = (((sxpar(uh0, 'parang') - (angoffset + u_dpas[ci]))+180d) mod 360d) - 180d
        for Li=0,sz[2]-1 do begin
            qslice = qcube[*,*,Li]
            iqslice = iqcube[*,*,Li]
            uslice = ucube[*,*,Li]
            iuslice = iucube[*,*,Li]
            qbad = where(~finite(qslice) or ~finite(iqslice))
            ubad = where(~finite(uslice) or ~finite(iuslice))
            qslice[qbad] = !values.f_nan
            iqslice[qbad] = !values.f_nan
            uslice[ubad] = !values.f_nan
            iuslice[ubad] = !values.f_nan
            qstar[Li,Ci] = total((qslice[unres_annulus]), /nan)
            iqstar[Li,Ci] = total((iqslice[unres_annulus]), /nan)
            ustar[Li,Ci] = total((uslice[unres_annulus]), /nan)
            iustar[Li,Ci] = total((iuslice[unres_annulus]), /nan)
        endfor
        ; Poly fit wavelength vs qrel and urel per cycle
        ; Fit only finite & non-telluric channels,
        ; then replace any non-finite values with the fit values 
        ; (leaving any finite telluric channels unchanged)
        qrel_ci = (qstar[*,ci]) / (iqstar[*,ci])
        qgood = finite(qrel_ci)
        qgood = where(qgood, complement=qbad)
        qgood_and_notellurics = finite(qrel_ci)
        qgood_and_notellurics[telluric_channels] = 0
        qgood_and_notellurics = where(qgood_and_notellurics)
        qfit = robust_poly_fit(lambda[qgood_and_notellurics], qrel_ci[qgood_and_notellurics], 2, numit=15)
        qrel_fit[*,ci] = qrel_ci
        qrel_fit[qbad,ci] = qfit[0] + qfit[1]*(lambda[qbad]) + qfit[2]*(lambda[qbad])^2.

        urel_ci = (ustar[*,ci]) / (iustar[*,ci])
        ugood = finite(urel_ci)
        ugood = where(ugood, complement=ubad)
        ugood_and_notellurics = finite(urel_ci)
        ugood_and_notellurics[telluric_channels] = 0
        ugood_and_notellurics = where(ugood_and_notellurics)
        ufit = robust_poly_fit(lambda[ugood_and_notellurics], urel_ci[ugood_and_notellurics], 2, numit=15)
        urel_fit[*,ci] = urel_ci
        urel_fit[ubad,ci] = ufit[0] + ufit[1]*(lambda[ubad]) + ufit[2]*(lambda[ubad])^2.
    endif

    qhcube[*,*,*,ci] = qcube
    uhcube[*,*,*,ci] = ucube
    iqhcube[*,*,*,ci] = iqcube
    iuhcube[*,*,*,ci] = iucube
    
    if keyword_set(save_derot_seq) then begin
        Qout_derot = Qout[ci].replace('.fits', '_derot.fits')
        writefits,Qout_derot,0,qh0
        writefits,Qout_derot,qcube,qh1,/append
        
        Uout_derot = Uout[ci].replace('.fits', '_derot.fits')
        writefits,Uout_derot,0,uh0
        writefits,Uout_derot,ucube,uh1,/append

        IQout_derot = IQout[ci].replace('.fits', '_derot.fits')
        writefits,IQout_derot,0,iqh0
        writefits,IQout_derot,iqcube,iqh1,/append
        
        IUout_derot = IUout[ci].replace('.fits', '_derot.fits')
        writefits,IUout_derot,0,iuh0
        writefits,IUout_derot,iucube,iuh1,/append
    endif
    
    if ci eq 0 then begin
        h0_main = qh0
        h1_main = qh1
        h1q_main = qh1
        h1u_main = uh1
    endif

    cycle_str = 'cyc_'+((string(cnum)).trim())
    qorignames = sxpar(qh1,'cycle')
    sxaddpar,h1q_main,cycle_str,qorignames
    uorignames = sxpar(uh1,'cycle')
    sxaddpar,h1u_main,cycle_str,uorignames
    orignames = qorignames+', '+uorignames
    sxaddpar,h1_main,cycle_str,orignames
endfor

if keyword_set(unres_pol_sub) then begin
    if keyword_set(avg_unres_pol_coeffs) then begin
        qrel_avg = median(qrel_fit, dimension=2, /even)
        urel_avg = median(urel_fit, dimension=2, /even)
    endif

    if keyword_set(make_unres_pol_plot) then begin
        ; Making a figure of the calculated unres pol per cycle for reference:
        loadct, 17, ncolors=ncycles, RGB_TABLE=colors
        pirel = sqrt(qrel_fit^2 + urel_fit^2)
        aolp_unres = (0.5*atan(urel_fit, qrel_fit))*!RADEG

        p_ymax = 100d*(max(pirel))
        p_yrange = [-p_ymax, p_ymax]
        for ci=0,ncycles-1 do begin
            color = reform(colors[ci,*],3)
            if ci eq 0 then begin
                p_q = PLOT(lambda, 100.d * qrel_fit[*, ci], symbol=0, color=color, DIM=[1600,800], linestyle='-', xtitle = 'Wavelength (nm)', ytitle='% q', layout=[4,2,1], margin=0.1, xrange = [1100,2420], yrange=p_yrange, /buffer)
                p_u = PLOT(lambda, 100.d * urel_fit[*, ci], symbol=0, color=color, linestyle='-', xtitle = 'Wavelength (nm)', ytitle='% u', layout=[4,2,2], margin=0.1, xrange = [1100,2420], yrange=p_yrange, /current)
                p_p = PLOT(lambda, 100.d * pirel[*, ci], symbol=0, color=color, linestyle='-', xtitle = 'Wavelength (nm)', ytitle='% pol', layout=[4,2,3], margin=0.1, xrange = [1100,2420], yrange=p_yrange, /current)
                p_a = PLOT(lambda, aolp_unres[*, ci], symbol=0, color=color, linestyle='-', xtitle = 'Wavelength (nm)', ytitle='AOLP', layout=[4,2,4], margin=0.1, xrange = [1100,2420], yrange=[-5,185], /current)

                p_qalt = PLOT([cycles[ci]], [qalt[ci]], symbol=24, sym_color='black', sym_fill_color=color, sym_filled=1, linestyle=6, xtitle = 'Cycle #', ytitle='Q altitudes', layout=[4,2,5], xrange = [cycles[0]-0.5, cycles[-1]+0.5], thick=5, margin=0.1, /current)
                p_ualt = PLOT([cycles[ci]], [ualt[ci]], symbol=24, sym_color='black', sym_fill_color=color, sym_filled=1, linestyle=6, xtitle = 'Cycle #', ytitle='U altitudes', layout=[4,2,6], xrange = [cycles[0]-0.5, cycles[-1]+0.5], thick=5, margin=0.1, /current)
                p_qpa = PLOT([cycles[ci]], [qparang[ci]], symbol=24, sym_color='black', sym_fill_color=color, sym_filled=1, linestyle=6, xtitle = 'Cycle #', ytitle='Q parangs', layout=[4,2,7], xrange = [cycles[0]-0.5, cycles[-1]+0.5], thick=5, margin=0.1, /current)
                p_upa = PLOT([cycles[ci]], [uparang[ci]], symbol=24, sym_color='black', sym_fill_color=color, sym_filled=1, linestyle=6, xtitle = 'Cycle #', ytitle='U parangs', layout=[4,2,8], xrange = [cycles[0]-0.5, cycles[-1]+0.5], thick=5, margin=0.1, /current)
            endif else begin
                p_q = PLOT(lambda, 100.d * qrel_fit[*, ci], symbol=0, color=color, linestyle='-', overplot=p_q)
                p_u = PLOT(lambda, 100.d * urel_fit[*, ci], symbol=0, color=color, linestyle='-', overplot=p_u)
                p_p = PLOT(lambda, 100.d * pirel[*, ci], symbol=0, color=color, linestyle='-', overplot=p_p)
                p_a = PLOT(lambda, aolp_unres[*, ci], symbol=0, color=color, linestyle='-', overplot=p_a)

                p_qalt = PLOT([cycles[ci]], [qalt[ci]], symbol=24, sym_color='black', sym_fill_color=color, sym_filled=1, thick=5, linestyle=6, overplot=p_qalt)
                p_ualt = PLOT([cycles[ci]], [ualt[ci]], symbol=24, sym_color='black', sym_fill_color=color, sym_filled=1, thick=5, linestyle=6, overplot=p_ualt)
                p_qpa = PLOT([cycles[ci]], [qparang[ci]], symbol=24, sym_color='black', sym_fill_color=color, sym_filled=1, thick=5, linestyle=6, overplot=p_qpa)
                p_upa = PLOT([cycles[ci]], [uparang[ci]], symbol=24, sym_color='black', sym_fill_color=color, sym_filled=1, thick=5, linestyle=6, overplot=p_upa)
            endelse
        endfor
        if keyword_set(avg_unres_pol_coeffs) then begin
            ; Plot the determined averages:
            pi_avg = (sqrt(qrel_avg^2 + urel_avg^2))
            aolp_avg = (0.5*atan(urel_avg, qrel_avg))*!RADEG
            p_q = PLOT(lambda, 100.d * qrel_avg, symbol=4, sym_filled=1, sym_color='white', sym_fill_color='black', linestyle=6, overplot=p_q, thick=5)
            p_u = PLOT(lambda, 100.d * urel_avg, symbol=4, sym_filled=1, sym_color='white', sym_fill_color='black', linestyle=6, overplot=p_u, thick=5)
            p_p = PLOT(lambda, 100.d * pi_avg, symbol=4, sym_filled=1, sym_color='white', sym_fill_color='black', linestyle=6, overplot=p_p, thick=5)
            p_a = PLOT(lambda, aolp_avg, symbol=4, sym_filled=1, sym_color='white', sym_fill_color='black', linestyle=6, overplot=p_a, thick=5)
        endif
        p_q.Save, outfile+"_unres_pol.png"
    endif   

    for ci=0,ncycles-1 do begin
        ; then replace the other values with the average for the actual subtraction
        if keyword_set(avg_unres_pol_coeffs) then begin
            qrel_fit[*,ci] = qrel_avg
            urel_fit[*,ci] = urel_avg
        endif
        for Li=0,sz[2]-1 do begin
            qslice = qhcube[*,*,Li,ci]
            iqslice = iqhcube[*,*,Li,ci]
            cq_Li = qrel_fit[Li,ci]

            uslice = uhcube[*,*,Li,ci]
            iuslice = iuhcube[*,*,Li,ci]
            cu_Li = urel_fit[Li,ci]

            qhcube[*,*,Li,ci] = qslice - (iqslice*cq_Li)
            uhcube[*,*,Li,ci] = uslice - (iuslice*cu_Li)
        endfor
        ; If saving derot stokes images, load and rewrite the Q and U images with additional header info
        if keyword_set(save_derot_seq) then begin
            Qout_derot = Qout[ci].replace('.fits', '_derot.fits')
            qh0 = headfits(Qout_derot, ext=0)
            qh1 = headfits(Qout_derot, ext=1)

            Uout_derot = Uout[ci].replace('.fits', '_derot.fits')
            uh0 = headfits(Uout_derot, ext=0)
            uh1 = headfits(Uout_derot, ext=1)

            ; Add the utilized coeffs to the headers:
            for Li=0,sz[2]-1 do begin
                sLi = (string(Li)).trim()
                sxaddpar, qh1, 'fcq_'+sLi, qrel_fit[Li,ci], 'Q unres pol coeff. for slice '+sLi
                sxaddpar, uh1, 'fcu_'+sLi, urel_fit[Li,ci], 'U unres pol coeff. for slice '+sLi
            end

            writefits,Qout_derot,0,qh0
            writefits,Qout_derot,(qhcube[*,*,*,ci]),qh1,/append
            writefits,Uout_derot,0,uh0
            writefits,Uout_derot,(uhcube[*,*,*,ci]),uh1,/append
        endif
    endfor    
endif
    
h1s = list(h1_main, h1q_main, h1u_main)
for i=0,2 do begin
    h1 = h1s[i]
    sxdelpar,h1,'cycle'; Removing "cycle" entry from original combo
    fxaddpar,h1,'command',command
    sxaddpar,h1,'prefname',prefname,'prefname'
    sxaddpar,h1,'Lsuffname',Lsuffname,'Lsuffname'
    sxaddpar,h1,'Rsuffname',Rsuffname,'Rsuffname'
    sxaddpar,h1,'hwpinfo_name',hwpinfo_name,'hwpinfo_name'
    sxaddpar,h1,'meanadd',keyword_set(meanadd),'meanadd'
    sxaddpar,h1,'collapse_type',collapse_type,'collapse_type'
    sxaddpar,h1,'angoffset',angoffset,'angoffset'
    sxaddpar,h1,'single_sum_zeros',keyword_set(single_sum_zeros),'single_sum_zeros'
    sxaddpar,h1,'no_mm_ip_corr',keyword_set(no_mm_ip_corr),'no_mm_ip_corr'
    sxaddpar,h1,'no_mm_crosstalk_corr',keyword_set(no_mm_crosstalk_corr),'no_mm_crosstalk_corr'
    if keyword_set(first_order_ip_corr) then sxaddpar,h1,'foip_rlims',strjoin((string(foip_rlims)).trim(),', '),'foip_rlims'
    if keyword_set(bg_sub) then sxaddpar,h1,'bg_rlims',strjoin((string(bg_rlims)).trim(),', '),'bg_rlims'
    if keyword_set(unres_pol_sub) then sxaddpar,h1,'unres_pol_rlims',strjoin((string(unres_pol_rlims)).trim(),', '),'unres_pol_rlims'
    sxaddpar,h1,'avg_unres_pol_coeffs',keyword_set(avg_unres_pol_coeffs),'avg_unres_pol_coeffs'
    sxaddpar,h1,'no_adisync',keyword_set(no_adisync),'no_adisync'
    sxaddpar,h1,'hwp_inds',strjoin((string(hwp_inds)).trim(),', '),'hwp_inds'
    if keyword_set(extra_pol_angle_offset) then sxaddpar,h1,'extra_pol_angle_offset',extra_pol_angle_offset,'extra_pol_angle_offset' else sxaddpar,h1,'extra_pol_angle_offset',0,'extra_pol_angle_offset'
    h1s[i] = h1
endfor

h1_main = h1s[0]
h1q_main = h1s[1]
h1u_main = h1s[2]

if keyword_set(collapse_cycles) then begin
    for ci=0,ncycles-1 do begin
        qcube_ci = qhcube[*,*,*,ci]
        ucube_ci = uhcube[*,*,*,ci]

        case collapse_type of
            'median': begin & qcol_ci=median(qcube_ci, dimension=3, /even) & ucol_ci=median(ucube_ci, dimension=3, /even) & end
            'rmean': begin & resistant_mean,qcube_ci,3.,qcol_ci,dimension=3 & resistant_mean,ucube_ci,3.,ucol_ci,dimension=3 & end
            'mean': begin & qcol_ci=mean(qcube_ci, dimension=3) & ucol_ci=mean(ucube_ci, dimension=3) & end
        endcase

        picol_ci = sqrt(qcol_ci^2 + ucol_ci^2)
        PIout_collapsed = Qout[ci].replace('_q.fits', '_pi_collapsed.fits')
        writefits,PIout_collapsed,0,h0_main
        writefits,PIout_collapsed,picol_ci,h1_main,/append     

        qout_collapsed = Qout[ci].replace('_q.fits', '_q_collapsed.fits')
        writefits,qout_collapsed,0,h0_main
        writefits,qout_collapsed,qcol_ci,h1_main,/append

        uout_collapsed = Qout[ci].replace('_q.fits', '_u_collapsed.fits')
        writefits,uout_collapsed,0,h0_main
        writefits,uout_collapsed,ucol_ci,h1_main,/append    
    endfor
endif

if ~keyword_set(meanadd) then begin
    qcube=median(qhcube, dimension=4, /even)
    ucube=median(uhcube, dimension=4, /even)
endif else begin
    resistant_mean,qhcube,3.,qcube,dimension=4
    resistant_mean,uhcube,3.,ucube,dimension=4
endelse

writefits,procdir+outfile+'_q.fits',0,h0_main
writefits,procdir+outfile+'_q.fits',qcube,h1q_main,/append

writefits,procdir+outfile+'_u.fits',0,h0_main
writefits,procdir+outfile+'_u.fits',ucube,h1u_main,/append

pi_cube = sqrt(qcube^2 + ucube^2)
writefits,procdir+outfile+'_pi.fits',0,h0_main
writefits,procdir+outfile+'_pi.fits',pi_cube,h1_main,/append

aolp = (0.5*atan(ucube, qcube))*!RADEG
writefits,procdir+outfile+'_aolp.fits',0,h0_main
writefits,procdir+outfile+'_aolp.fits',aolp,h1_main,/append

; Collapse Q and U wavelength cubes to a single broadband image
case collapse_type of
    'median': begin & qcol=median(qcube, dimension=3, /even) & ucol=median(ucube, dimension=3, /even) & end
    'rmean': begin & resistant_mean,qcube,3.,qcol,dimension=3 & resistant_mean,ucube,3.,ucol,dimension=3 & end
    'mean': begin & qcol=mean(qcube, dimension=3) & ucol=mean(ucube, dimension=3) & end
endcase

writefits,procdir+outfile+'_q_collapsed.fits',0,h0_main
writefits,procdir+outfile+'_q_collapsed.fits',qcol,h1q_main,/append
writefits,procdir+outfile+'_u_collapsed.fits',0,h0_main
writefits,procdir+outfile+'_u_collapsed.fits',ucol,h1u_main,/append

pi_col = sqrt(qcol^2 + ucol^2)
writefits,procdir+outfile+'_pi_collapsed.fits',0,h0_main
writefits,procdir+outfile+'_pi_collapsed.fits',pi_col,h1_main,/append

aolp_col = (0.5*atan(ucol, qcol))*!RADEG
writefits,procdir+outfile+'_aolp_collapsed.fits',0,h0_main
writefits,procdir+outfile+'_aolp_collapsed.fits',aolp_col,h1_main,/append

coverage_cube = product(finite((qhcube + uhcube)), 3)
if n_elements(size(coverage_cube, /dim)) lt 3 then coverage_map = coverage_cube else coverage_map = total(coverage_cube, 3) / double(ncycles)

writefits,'coverage_map.fits',0,h0_main
writefits,'coverage_map.fits',coverage_map,h1_main,/append

qphicube = get_azimuthal_stokes(qcube, ucube, phi, uphi=uphicube)
writefits,procdir+outfile+'_qphi.fits',0,h0_main
writefits,procdir+outfile+'_qphi.fits',qphicube,h1_main,/append
writefits,procdir+outfile+'_uphi.fits',0,h0_main
writefits,procdir+outfile+'_uphi.fits',uphicube,h1_main,/append

qphicol = get_azimuthal_stokes(qcol, ucol, phi, uphi=uphicol)
writefits,procdir+outfile+'_qphi_collapsed.fits',0,h0_main
writefits,procdir+outfile+'_qphi_collapsed.fits',qphicol,h1_main,/append
writefits,procdir+outfile+'_uphi_collapsed.fits',0,h0_main
writefits,procdir+outfile+'_uphi_collapsed.fits',uphicol,h1_main,/append

endofprogram:
end
