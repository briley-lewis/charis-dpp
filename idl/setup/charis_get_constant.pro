function charis_get_constant, name = name, mjd = mjd, label = label
  compile_opt idl2
  ; Instrument constants by epoch.  [BUGS.md #124, #302, #306, #307, #316]
  ;
  ; The source of truth is charis_dpp/data/instrument.toml in the Python
  ; package; this table mirrors its dated [[epoch]] blocks and a Python test
  ; (tests/unit/test_idl_parity.py) checks that the two agree.  Keep the
  ; `epochs.add, {...}` lines on one line each and the tag names unchanged:
  ; the test parses them.
  ;
  ; Usage: charis_get_constant(name='pixscale', mjd=sxpar(h0,'MJD'))
  ; Without mjd= the session default set by charis_set_epoch is used; if none
  ; has been set, the most recent epoch is used (the pre-existing behaviour
  ; of this routine, which always returned the latest values).  Later epochs
  ; only state what changed; earlier values carry over.  [BUGS.md #322]
  ;
  ; Spot guesses are absolute positions on a 201x201 frame (the Python
  ; package stores them relative to the frame centre).
  epochs = list()
  epochs.add, {mjd: 57388.0d, label: 'commissioning', dtel: 7.9d0, pixscale: 0.0162d0, angoffset: 2.2d0, satcon25: 7.5037d0, satcon50: 5.991d0, satconerr: 0.06d0, lat: 19.825504d0, lng: -155.47602d0, parity: 1, rotzp: 34.7d0, m3_m1: 2.104d0, m3_b1: 14.20d0, m3_m2: 2.100d0, m3_b2: 13.20d0, readtime: 1.475d0}
  epochs.add, {mjd: 57995.0d, label: 'satspot-recal-2017', satcon25: 6.412d0, satcon50: 4.92d0, satconerr: 0.05d0}
  epochs.add, {mjd: 60601.0d, label: 'Chen2023', dtel: 7.92d0, pixscale: 0.01615d0, angoffset: 2.03d0}
  epochs.add, {mjd: 60933.0d, label: 'post-NBS', parity: -1}
  epochs.add, {mjd: 61146.0d, label: 'M3-recoat-2026', m3_m1: 1.768d0, m3_b1: 12.77d0, m3_m2: 2.455d0, m3_b2: 15.67d0}

  common charis_epoch, charis_epoch_mjd
  if n_elements(mjd) gt 0 then begin
    mjd_use = double(mjd[0])
  endif else if n_elements(charis_epoch_mjd) gt 0 && charis_epoch_mjd gt 0 then begin
    mjd_use = double(charis_epoch_mjd)
  endif else begin
    mjd_use = 1d9
  endelse
  values = hash()
  label = ''
  foreach ep, epochs do begin
    if ep.mjd gt mjd_use then continue
    tags = tag_names(ep)
    for i = 0, n_elements(tags) - 1 do values[strlowcase(tags[i])] = ep.(i)
    label = ep.label
  endforeach

  case strlowcase(name) of
    'dtel': value = values['dtel']
    'pixscale': value = values['pixscale']
    'angoffset': value = values['angoffset']
    'lat': value = values['lat']
    'lng': value = values['lng']
    'parity': value = values['parity']
    'rotzp': value = values['rotzp']
    'readtime': value = values['readtime']
    'satcon25': value = values['satcon25']
    'satcon50': value = values['satcon50']
    'satconerr': value = values['satconerr']
    'm3_m1': value = values['m3_m1']
    'm3_b1': value = values['m3_b1']
    'm3_m2': value = values['m3_m2']
    'm3_b2': value = values['m3_b2']
    'epoch': value = label
    'lowspots_x': value = [59.2, 112.5, 141.75, 87]
    'lowspots_y': value = [113.5, 141, 86.5, 59]
    'jspots_x': value = [67.9, 110, 132., 89.9]
    'jspots_y': value = [110.7, 132.2, 89.3, 67.7]
    'hspots_x': value = [59.2, 112.9, 141.8, 86.96]
    'hspots_y': value = [113.5, 141.1, 86.5, 58.87]
    'kspots_x': value = [44.5, 117.5, 157.06, 82.5]
    'kspots_y': value = [118.5, 156, 81.5, 44]
    else: message, 'charis_get_constant: unknown constant ' + string(name)
  endcase
  return, value
end
