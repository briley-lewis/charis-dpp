pro charis_compile_check, root = root, quiet = quiet, nfail = nfail
  compile_opt idl2
  ; Compile every .pro under the CHARIS-DPP idl/ tree and report the ones
  ; that fail, then resolve the full call tree of the main entry points.
  ; Usage:  idl -e "charis_compile_check"      (idl/ must be on !path)
  ; Sets nfail= to the number of files that did not compile.
  if n_elements(root) eq 0 then root = file_dirname(file_dirname(routine_filepath('charis_compile_check')))
  files = file_search(root, '*.pro', count = nf)
  failed = list()
  for i = 0, nf - 1 do begin
    name = file_basename(files[i], '.pro')
    catch, err
    if err ne 0 then begin
      catch, /cancel
      failed.add, files[i] + '  :: ' + !error_state.msg
      continue
    endif
    resolve_routine, name, /either, /compile_full_file, /no_recompile
    catch, /cancel
  endfor
  nfail = failed.count()
  print, strtrim(nf, 2), ' files, ', strtrim(nfail, 2), ' failed to compile'
  foreach f, failed do print, '  FAIL ', f
  if ~keyword_set(quiet) then begin
    print, 'resolve_all from the PDI / registration entry points ...'
    resolve_all, /continue_on_error, /quiet, $
      resolve_procedure = ['charis_imprep', 'charis_pdi_flatfield', 'charis_pdi_subtract_sky', $
      'charis_pdi_register_cube', 'charis_register_cube', 'charis_pdi_specphot_cal', $
      'charis_pdi_hwp_match', 'charis_pdi', 'charis_set_epoch'], $
      resolve_function = ['charis_get_constant', 'm3_model', 'get_charis_mueller_matrix']
    print, 'resolve_all done'
  endif
end
