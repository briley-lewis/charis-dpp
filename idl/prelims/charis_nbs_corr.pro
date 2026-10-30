function charis_nbs_corr,cube,h0,h1,rotzp=rotzp

;*****
;a simple function to correct the parity of post NBS CHARIS cubes and recalculate the WCS solution
;valid for CHARIS data taken after ~Sept 2025 (MJD 60933.0)
;we will assume north up, east-left

;also requires a fix to charis_register_cube to get guess positions of sat spots correct
;rotzp/north PA and pixel scale to be revised
;*****

if ~keyword_set(rotzp) then rotzp = 34.7
;if ~keyword_set(rotzp) then rotzp = 32.3
;if ~keyword_set(rotzp) then rotzp = 33.5
extrarot=113-rotzp


;correct the rotation offset

;left handed: impose this
totrot=-1*(sxpar(h1,'PARANG',/silent)+extrarot)

;right handed
;totrot=1*(sxpar(h1,'PARANG')+extrarot)
;sxaddpar,h1,'TOT_ROT',totrot



;left handed
xp=sxpar(h1,'XPIXSCAL',/silent)

;this is what you do if you want a right handed coordinate system
;xp=-1*sxpar(h1,'XPIXSCAL')
;sxaddpar,h1,'XPIXSCAL',xp

yp=sxpar(h1,'YPIXSCAL',/silent)

;****Delete Erroneous FITS header information****
;this will eliminate the WCS entirely

;goto,skipnowcs
sxdelpar,h1,'CD1_1'
sxdelpar,h1,'CD1_2'
sxdelpar,h1,'CD2_1'
sxdelpar,h1,'CD2_2'
sxdelpar,h1,'XPIXSCAL'
sxdelpar,h1,'YPIXSCAL'
sxdelpar,h1,'TOT_ROT'

;print,sxpar(h0,'CD1_1')

sxdelpar,h0,'CD1_1'
sxdelpar,h0,'CD1_2'
sxdelpar,h0,'CD2_1'
sxdelpar,h0,'CD2_2'
sxdelpar,h0,'XPIXSCAL'
sxdelpar,h0,'YPIXSCAL'
sxdelpar,h0,'TOT_ROT'
skipnowcs:

;*****
;Recalculate the WCS solution based off of the new calibration, save in FITS header
;*****

;angle on detector between Y axis and true north (estimate)
angle=!pi*totrot/180.

;populate primary FITS header and first extension

sxaddpar,h1,'COMMENT','post-NBS Parity Flip'
sxaddpar,h1,'parity',1
sxaddpar,h1,'CD1_1',cos(angle)*xp
sxaddpar,h1,'CD1_2',sin(angle)*xp
sxaddpar,h1,'CD2_1',-sin(angle)*yp
sxaddpar,h1,'CD2_2',cos(angle)*yp
sxaddpar,h1,'XPIXSCAL',xp
sxaddpar,h1,'YPIXSCAL',yp
sxaddpar,h1,'ROT_ZP',rotzp,'zp rot (degrees)'
sxaddpar,h1,'EX_ROT',extrarot,'extrarot:113-zp (degrees)'
sxaddpar,h1,'TOT_ROT',totrot,'Total rotation angle (degrees)'

sxaddpar,h0,'COMMENT','post-NBS Parity Flip'
sxaddpar,h0,'parity',1,'flipped from N-u,E-r to N-u, E-l'
sxaddpar,h0,'CD1_1',cos(angle)*xp
sxaddpar,h0,'CD1_2',sin(angle)*xp
sxaddpar,h0,'CD2_1',-sin(angle)*yp
sxaddpar,h0,'CD2_2',cos(angle)*yp
sxaddpar,h0,'XPIXSCAL',xp
sxaddpar,h0,'YPIXSCAL',yp
sxaddpar,h0,'ROT_ZP',rotzp,'zp rot (degrees)'
sxaddpar,h0,'EX_ROT',extrarot,'extrarot:113-zp (degrees)'
sxaddpar,h0,'TOT_ROT',totrot,'Total rotation angle (degrees)'

;****North-Up, East-Left correction

;flip parity for x axis
cube_corrected = reverse(cube,1)

;****

;****Return Results
; Return as an IDL structure (like a Python dict)
results = {cube: cube_corrected, h0: h0, h1: h1}

return, results


end




