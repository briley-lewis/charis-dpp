function charis_compile_check_missing, path, names, duplicates = duplicates
  compile_opt idl2, hidden
  ; Names of the routines declared in `path` that are not compiled at all
  ; (empty when every declared routine compiled). A declared routine that is
  ; compiled from a *different* file (the same name defined twice in the
  ; tree, so /no_recompile kept the first one) is returned in duplicates=
  ; as 'NAME <- other/file.pro' instead; the tree still compiles.
  p = routine_info(/source)
  f = routine_info(/source, /functions)
  src = file_expand_path([p.path, f.path])
  rn = [p.name, f.name]
  missing = list()
  duplicates = list()
  foreach n, names do begin
    w = where(rn eq strupcase(n), nw)
    if nw eq 0 then begin
      missing.add, n
      continue
    endif
    if total(src[w] eq file_expand_path(path)) eq 0 then duplicates.add, strupcase(n) + ' <- ' + src[w[0]]
  endforeach
  duplicates = duplicates.toArray(type = 'string')
  return, missing.toArray(type = 'string')
end

pro charis_compile_check, root = root, quiet = quiet, nfail = nfail
  compile_opt idl2
  ; Compile every .pro under the CHARIS-DPP idl/ tree and report the ones
  ; that fail, then resolve the full call tree of the main entry points.
  ; Usage:  idl -e "charis_compile_check"      (idl/ must be on !path)
  ; Sets nfail= to the number of files that did not compile.
  ;
  ; A file counts as compiled when every pro/function it declares shows up
  ; in routine_info(/source) with that file as its source afterwards, so a
  ; file whose routine name differs from the file name (charis_rdiklip_2.pro
  ; defining charis_rdiklip) still passes, and a multi-routine file whose
  ; last routine has a syntax error still fails. Files with no pro/function header are
  ; @include fragments (craigm/fxfilter.pro and friends): they are skipped,
  ; since compiling them as $MAIN$ fails and poisons later includes.
  if n_elements(root) eq 0 then root = file_dirname(file_dirname(routine_filepath('charis_compile_check')))
  files = file_search(root, '*.pro', count = nf)
  failed = list()
  skipped = list()
  dups = list()
  for i = 0, nf - 1 do begin
    name = file_basename(files[i], '.pro')
    nl = file_lines(files[i])
    if nl eq 0 then begin
      skipped.add, files[i]
      continue
    endif
    txt = strarr(nl)
    openr, lun, files[i], /get_lun
    readf, lun, txt
    free_lun, lun
    hdr = stregex(txt, '^[[:space:]]*(pro|function)[[:space:]]+([A-Za-z0-9_$:]+)', /subexpr, /extract, /fold_case)
    w = where(hdr[0, *] ne '', nhdr)
    if nhdr eq 0 then begin
      skipped.add, files[i]
      continue
    endif
    declared = reform(hdr[2, w])
    catch, err
    if err ne 0 then begin
      catch, /cancel
      ; resolve_routine fails with "undefined procedure/function" when the
      ; routine name differs from the file name even though the file
      ; compiled; only an undefined declared routine is a real failure
      miss = charis_compile_check_missing(files[i], declared, duplicates = dd)
      if n_elements(dd) gt 0 then dups.add, files[i] + '  :: ' + strjoin(dd, ', ')
      if n_elements(miss) eq 0 then continue
      failed.add, files[i] + '  :: not compiled: ' + strjoin(miss, ', ') + '  (' + !error_state.msg + ')'
      continue
    endif
    resolve_routine, name, /either, /compile_full_file, /no_recompile
    catch, /cancel
    miss = charis_compile_check_missing(files[i], declared, duplicates = dd)
    if n_elements(dd) gt 0 then dups.add, files[i] + '  :: ' + strjoin(dd, ', ')
    if n_elements(miss) gt 0 then failed.add, files[i] + '  :: not compiled: ' + strjoin(miss, ', ')
  endfor
  nfail = failed.count()
  print, strtrim(nf, 2), ' files, ', strtrim(nfail, 2), ' failed to compile, ', $
    strtrim(skipped.count(), 2), ' skipped (no pro/function header)'
  foreach f, failed do print, '  FAIL ', f
  foreach f, skipped do print, '  SKIP ', f
  foreach f, dups do print, '  DUP  ', f
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
