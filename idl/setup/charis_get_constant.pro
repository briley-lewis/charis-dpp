function charis_get_constant, name = name
  compile_opt idl2

  ; in case you need to change this ...
  Dtel = 7.92d0
  pixscale = 0.01615

  case name of
    'Dtel': begin
      value = Dtel 
    end

    'pixscale': begin
      value = pixscale 
    end

    'angoffset': begin
      value = 2.03
    end

    'lowspots_x': begin
      value = [59.2, 112.5, 141.75, 87]
    end
    'lowspots_y': begin
      value = [113.5, 141, 86.5, 59]
    end

    'jspots_x': begin
      value = [67.9, 110, 132., 89.9]
    end
    'jspots_y': begin
      value = [110.7, 132.2, 89.3, 67.7]
    end

    'hspots_x': begin
      value = [59.2, 112.9, 141.8, 86.96]
    end
    'hspots_y': begin
      value = [113.5, 141.1, 86.5, 58.87]
    end

    'kspots_x': begin
      value = [44.5, 117.5, 157.06, 82.5]
    end
    'kspots_y': begin
      value = [118.5, 156, 81.5, 44]
      ; value=[113.5,141,86.5,59]
    end

  endcase

  return, value
end
