function derotator_model, lambda, delta_45_der = delta_45_der, delta_r_der = delta_r_der, epsilon_der = epsilon_der
  compile_opt idl2
  ; Linearly interpolates per-wavelength fits from McIntosh+2026

  ; Read the CSV file with wavelength-fit elliptical retardance params
  csv_path = file_dirname(routine_filepath()) + '/mcintosh_2026_spie_wavelength_fits.csv'
  data = read_csv(csv_path, header = headers)
  headers = strlowcase(headers)

  wavelength_idx = where(headers eq 'wavelength_bin', n_wavelength_idx)
  phi_45_idx = where(headers eq 'image_rotator_phi_45', n_phi_45_idx)
  phi_r_idx = where(headers eq 'image_rotator_phi_r', n_phi_r_idx)
  phi_h_idx = where(headers eq 'image_rotator_phi_h', n_phi_h_idx)

  wavelengths = data.(wavelength_idx[0])
  phi_45 = data.(phi_45_idx[0])
  phi_r = data.(phi_r_idx[0])
  phi_h = data.(phi_h_idx[0])
  
  ; Interpolate retardance values at the given wavelength
  delta_h_der = interpol(phi_h, wavelengths, lambda)
  delta_45_der = interpol(phi_45, wavelengths, lambda)
  delta_r_der = interpol(phi_r, wavelengths, lambda)
  epsilon_der = 0.d
  
  ; epsilon_der is just returned as input value
  
  return, delta_h_der
end