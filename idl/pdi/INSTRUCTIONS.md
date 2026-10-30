# CHARIS PDI Instructions

*Originally written by Kellen Lawson (2022). Updated by Miles Lucas (2023).* Updated again by Thomas McIntosh (2026)

Assumes you're in a directory containing 
1. a sequence of extracted CHARIS PDI cubes (including sky frames)
2. a VAMPIRES half waveplate log file (e.g. vampHWPLog_*.csv)

For this example, we'll assume the following:
- Target is AB Aur
- Observations are from 2020/10/03.
- Used low-res broadband mode
- Used a coronagraph & satellite spots at 11.2 lambda/D spacing and 25nm amplitude throughout.

## Image preparation

```
charis_newobs, date='20201003', object='ab_aur', filter='low'
```
```
charis_imprep 
```
Open ab_aur_low.info, and correct:
1. fnum_sat and fnum_sky (the DPP often mistakes PDI images as sky frames. You should visually inspect your images to make sure these are correct.
2. sptype (to whatever spectral type you'd like to assume for specphotc/≥al
3. For pdi flat-fielding, `brightflatdir` should point to a directory containing extracted PDI flat images taken with your data. 
4. Optionally, if PDI dark flats were collected, you can point to them with the `darkflatdir` keyword.

## Image calibration

```
charis_pdi_flatfield, 'ab_aur_low.info', brightflatdir = '../pdi_flats/'
```
```
charis_pdi_subtract_sky, 'ab_aur_low.info', /scalesky
```

## Image registration
```
A, 'ab_aur_low.info', astrogrid=11.2
```

## Spectrophotometric calibration
```
charis_pdi_specphot_cal, 'ab_aur_low.info', modamp=25, ap_factor=3, starlib=1, av=av # where we basically account for reddening (read this in the related paper(s) of the target)
```

## PDI

First, match the HWPs using the VAMPIRES HWP log
```
charis_pdi_hwp_match, 'ab_aur_low.info'
```

Then, do PDI
```
charis_pdi, 'ab_aur_low.info', /meanadd
```