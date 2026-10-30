function charis_path, pathname = pathname
  compile_opt idl2
  ; [BUGS.md #123, #319] paths come from the environment:
  ;   CHARIS_DPP_PATH    the idl/ directory of this repository
  ;   CHARIS_DPP_BTSETTL the BT-Settl model grid (lte*.dat.txt files)
  charispath = getenv('CHARIS_DPP_PATH')
  if charispath eq '' then begin
    ; fall back to the directory holding this file (idl/setup/ -> idl/)
    charispath = file_dirname(file_dirname(routine_filepath('charis_path', /is_function))) + '/'
  endif
  if strmid(charispath, strlen(charispath) - 1) ne '/' then charispath += '/'
  btsettldirectory = getenv('CHARIS_DPP_BTSETTL')
  if btsettldirectory eq '' then btsettldirectory = charispath + 'models/btsettl/'
  if strmid(btsettldirectory, strlen(btsettldirectory) - 1) ne '/' then btsettldirectory += '/'
  modeldirectory = charispath + 'models/'
  wavecaldirectory = charispath + 'cals/'
  filterresponsedirectory = charispath + 'tools/filters/filter_response/'
  planetdirectory = charispath + 'models/planetmodels/'
  modelspectrumdirectory = charispath + 'models/starmodels/'
  bonnefoydirectory = charispath + 'models/planetmodels/bonn_speclib/'
  montrealdirectory = charispath + 'models/planetmodels/Montreal/'
  case pathname of
    'charispath': pathdirectory = charispath
    'modeldir': pathdirectory = modeldirectory
    'wavecalpath': pathdirectory = wavecaldirectory
    'planetdir': pathdirectory = planetdirectory
    'filtresponsedir': pathdirectory = filterresponsedirectory
    'modelspecdir': pathdirectory = modelspectrumdirectory
    'montrealdir': pathdirectory = montrealdirectory
    'bonnefoydir': pathdirectory = bonnefoydirectory
    'btsettldir': pathdirectory = btsettldirectory
    else: message, 'charis_path: unknown pathname ' + string(pathname)
  endcase
  return, pathdirectory
end
