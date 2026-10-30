; charis_pdi_rlims_gridsearch.pro
; Grid search over unres_pol_rlims [rmin, rmax] combinations for charis_pdi.
; For each (rmin, rmax) pair, runs charis_pdi with a unique outfile prefix
; so that downstream Python plotting can load and compare all reductions.
;
; Usage:
;   .run pdi_rlims_gridsearch
;   pdi_rlims_gridsearch
;
; Outputs:  ./reduc/proc/pdi_rmin??_rmax??_*col*.fits  (wavelength-collapsed images)
;           ./pdi_rlims_gridsearch.log
;
; Written to accompany charis_pdi butterfly pattern grid check.

pro charis_pdi_rlims_gridsearch

  ; -------------------------------------------------------------------
  ; USER SETTINGS
  ; -------------------------------------------------------------------
  pfname = 'V2254Ori_low.info'

  ; rmin values to sweep (pixels)
  rmin_vals = [0, 2, 4, 6, 8, 10]

  ; rmax values to sweep (pixels)
  rmax_vals = [5, 10, 15, 20, 25, 30]

  ; Base flags (always applied)
  ; /meanadd, /bg_sub, /unres_pol_sub are always set
  ; outfile is varied per combo — see below

  logfile = 'pdi_rlims_gridsearch.log'
  ; -------------------------------------------------------------------

  ; Open log file
  openw, lun, logfile, /get_lun
  printf, lun, 'pdi_rlims_gridsearch.pro'
  printf, lun, 'Started: ' + systime()
  printf, lun, 'pfname  : ' + pfname
  printf, lun, ''
  printf, lun, string(format='(A6, A8, A8, A30)', 'Run', 'rmin', 'rmax', 'outfile')
  printf, lun, string(replicate('-', 60))

  n_rmin = n_elements(rmin_vals)
  n_rmax = n_elements(rmax_vals)
  n_total = n_rmin * n_rmax
  run_num = 0

  for i = 0, n_rmin - 1 do begin
    for j = 0, n_rmax - 1 do begin

      rmin = rmin_vals[i]
      rmax = rmax_vals[j]

      ; Skip invalid combos where rmin >= rmax
      if (rmin ge rmax) then begin
        printf, lun, string(format='(I6, I8, I8, A30)', run_num, rmin, rmax, 'SKIPPED (rmin>=rmax)')
        continue
      endif

      ; Build unique outfile name: pdi_rmin00_rmax00
      outfile = string(format='("pdi_rmin", I02, "_rmax", I02)', rmin, rmax)

      printf, lun, string(format='(I6, I8, I8, A30)', run_num, rmin, rmax, outfile)
      flush, lun

      print, ''
      print, '============================================================'
      print, string(format='("Run ", I0, "/", I0, " — rmin=", I0, "  rmax=", I0, "  outfile=", A)', $
        run_num, n_total, rmin, rmax, outfile)
      print, '============================================================'

      ; Call charis_pdi with this combo
      charis_pdi, pfname, $
        /meanadd, $
        /bg_sub, $
        /unres_pol_sub, $
        unres_pol_rlims = [rmin, rmax], $
        outfile = outfile

      run_num++

    endfor
  endfor

  printf, lun, ''
  printf, lun, 'Finished: ' + systime()
  printf, lun, string(format='("Completed ", I0, " runs.")', run_num)
  free_lun, lun

  print, ''
  print, 'Grid search complete.  Log written to: ' + logfile
  print, string(format='("Completed ", I0, " runs.")', run_num)

end