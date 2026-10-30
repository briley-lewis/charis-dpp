;+
; NAME:
;   SetFitParm.pro
;
; AUTHOR:
;	F.Bringezu, denet - Internetservice, Halle Germany,
;	bringezu@denet.de
;
; PURPOSE:
;   Provide a widget interface for creating a parinfo structure.
;   This parinfo structure can by used by mpfit routines of Craig B. Markwardt.
;
; MAJOR TOPICS:
;   Widget, mpfit.
;
; CALLING SEQUENCE:
;   parinfo=SetFitParm(used_parinfo)
;
; DESCRIPTION:
;
;   SetFitParm creates PARINFO using a widget interface.
;   PARINFO provides constraints for paramters used by the mpfit routines.
;
;   PARINFO is an array of structures, one for each parameter.
;
;   A detailed description can be found in the documentation of mpcurvefit.pro
;   This routine creates an array that contains a structure for each element.
;   The structure has the following entries.
;
;   - VALUE (DOUBLE): The starting parameter
;   - FIXED (BOOLEAN): 1 fix the parameter, 0 don't fix it at the
;     point given in VALUE.
;   - LIMITS (DBLARRAY(2)): Set upper and lower limit.
;   - LIMITED (BOOLEAN ARRAY 2):  Fix the limit.
;
;
;   The parameter OLDPARINFO is optional. OLDPARINFO is used to set
;   the default values in the widget.
;
;   You can simply run:
;   test=SetFitParm() to create the array for the first time.
;   Once the array is created it can be used to set the default values
;   in the widget by calling
;
;   test2=SetFitParm(test)
;
; INPUTS:
;
;
; OPTIONAL INPUTS:
;
;   OLDFITPARM - The default values of the new array
;
; INPUT KEYWORD PARAMETERS:
;
;   PARENT - if this widget is to be a child, set this keyword to the
;            parent widget ID.
;
; OUTPUT KEYWORD PARAMETERS:
;
;   CANCEL - if the user selected the cancel button on the SETFITPARM
;            widget, then this keyword will be set upon exit.
;
; OUTPUTS:
;   PARINFO array of structures
;
; SEE ALSO:
;   mpcurvefit
;
; MODIFICATION HISTORY:
;   Written, FB, 12/1999
;   Documented, FB, Jan 2000
;   Generalized positioning code, CM 01 Feb 2000
;
;-
; Copyright (C) 1999, F.Bringezu
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy and distribute unmodified copies for
; non-commercial purposes, and to modify and use for personal or
; internal use, is granted.  All other rights are reserved.
; -

pro SetFitParm_Events, event
  compile_opt idl2

  widget_control, event.id, get_value = buttonValue
  widget_control, event.id, get_uvalue = buttonUValue

  widget_control, event.top, get_uvalue = info, /no_copy

  case buttonUValue of
    'Cancel': widget_control, event.top, /destroy

    'FIX_HEIGHT': (*info.ptr).fparm(0).fixed = fix(event.select)
    'FIX_XMAX': (*info.ptr).fparm(1).fixed = fix(event.select)
    'FIX_WIDTH': (*info.ptr).fparm(2).fixed = fix(event.select)
    'FIX_OFFSET': (*info.ptr).fparm(3).fixed = fix(event.select)
    'FIX_SLOPE': (*info.ptr).fparm(4).fixed = fix(event.select)

    'LIMIT_HEIGHT_LOW': (*info.ptr).fparm(0).limited(0) = fix(event.select)
    'LIMIT_XMAX_LOW': (*info.ptr).fparm(1).limited(0) = fix(event.select)
    'LIMIT_WIDTH_LOW': (*info.ptr).fparm(2).limited(0) = fix(event.select)
    'LIMIT_OFFSET_LOW': (*info.ptr).fparm(3).limited(0) = fix(event.select)
    'LIMIT_SLOPE_LOW': (*info.ptr).fparm(4).limited(0) = fix(event.select)

    'LIMIT_HEIGHT_UP': (*info.ptr).fparm(0).limited(1) = fix(event.select)
    'LIMIT_XMAX_UP': (*info.ptr).fparm(1).limited(1) = fix(event.select)
    'LIMIT_WIDTH_UP': (*info.ptr).fparm(2).limited(1) = fix(event.select)
    'LIMIT_OFFSET_UP': (*info.ptr).fparm(3).limited(1) = fix(event.select)
    'LIMIT_SLOPE_UP': (*info.ptr).fparm(4).limited(1) = fix(event.select)

    'Accept': begin
      ; OK, get the information the user put into the form.
      ; Should do error checking, but...maybe later!

      widget_control, info.text(0), get_value = height
      widget_control, info.text(1), get_value = xmax
      widget_control, info.text(2), get_value = width
      widget_control, info.text(3), get_value = offset
      widget_control, info.text(4), get_value = slope

      widget_control, info.text(5), get_value = ll_height
      widget_control, info.text(6), get_value = ll_xmax
      widget_control, info.text(7), get_value = ll_width
      widget_control, info.text(8), get_value = ll_offset
      widget_control, info.text(9), get_value = ll_slope

      widget_control, info.text(10), get_value = ul_height
      widget_control, info.text(11), get_value = ul_xmax
      widget_control, info.text(12), get_value = ul_width
      widget_control, info.text(13), get_value = ul_offset
      widget_control, info.text(14), get_value = ul_slope

      ; Fill out the data structure with information
      ; collected from the form.

      (*info.ptr).fparm(0).value = height
      (*info.ptr).fparm(1).value = xmax
      (*info.ptr).fparm(2).value = width
      (*info.ptr).fparm(3).value = offset
      (*info.ptr).fparm(4).value = slope

      (*info.ptr).fparm(0).limits(0) = ll_height
      (*info.ptr).fparm(1).limits(0) = ll_xmax
      (*info.ptr).fparm(2).limits(0) = ll_width
      (*info.ptr).fparm(3).limits(0) = ll_offset
      (*info.ptr).fparm(4).limits(0) = ll_slope

      (*info.ptr).fparm(0).limits(1) = ul_height
      (*info.ptr).fparm(1).limits(1) = ul_xmax
      (*info.ptr).fparm(2).limits(1) = ul_width
      (*info.ptr).fparm(3).limits(1) = ul_offset
      (*info.ptr).fparm(4).limits(1) = ul_slope

      (*info.ptr).cancel = 0

      ; Destroy the widget program

      widget_control, event.top, /destroy
    end
    else:
  endcase

  if buttonValue ne 'Cancel' and buttonValue ne 'Accept' then begin
    widget_control, event.top, set_uvalue = info, /no_copy
  endif
end
; *******************************************************************

function SetFitParm, thisFParm, $
  parent = parent, $
  cancel = cancel
  compile_opt idl2

  on_error, 2

  if n_elements(thisFParm) eq 0 then $
    fparm = replicate({value: 0.d, fixed: 0, limited: [0, 0], limits: [0.d, 0]}, 5) $
  else fparm = thisFParm

  device, get_screen_size = screenSize
  xCenter = fix(screenSize[0] / 2.)
  yCenter = fix(screenSize[1] / 2.)
  xoff = xCenter - 150
  yoff = yCenter - 150

  ; Create a top-level base. Must have a Group Leader defined
  ; for Modal operation. If this widget is NOT modal, then it
  ; should only be called from the IDL command line as a blocking
  ; widget.

  if n_elements(parent) ne 0 then $
    tlb = widget_base(group_leader = parent $
    , uname = 'WID_BASE_0', xoffset = xoff, yoffset = yoff $
    , xsize = 380, ysize = 320, title = 'Parinfo Setup', space = 3, xpad = 3 $
    , ypad = 3, /floating, /modal, /base_align_center) else $
    tlb = widget_base(uname = 'WID_BASE_0', xoffset = xoff, yoffset = yoff $
    , xsize = 380, ysize = 320, title = 'Parinfo Setup', space = 3, xpad = 3 $
    , ypad = 3, /base_align_center)

  WID_BASE_1 = widget_base(tlb, uname = 'WID_BASE_1', frame = 1 $
  , xoffset = 10, yoffset = 55, xsize = 350, ysize = 200 $
  , title = 'IDL', space = 3, xpad = 3, ypad = 3)

  WID_LABEL_0 = widget_label(WID_BASE_1, uname = 'WID_LABEL_0' $
  , xoffset = 10, yoffset = 45, xsize = 42, ysize = 18 $
  , /align_left, value = 'Ampl.')

  WID_LABEL_1 = widget_label(WID_BASE_1, uname = 'WID_LABEL_1' $
  , xoffset = 10, yoffset = 75, xsize = 42, ysize = 18 $
  , /align_left, value = 'X(max)')

  WID_LABEL_2 = widget_label(WID_BASE_1, uname = 'WID_LABEL_2' $
  , xoffset = 10, yoffset = 105, xsize = 42, ysize = 18 $
  , /align_left, value = 'Width')

  WID_LABEL_3 = widget_label(WID_BASE_1, uname = 'WID_LABEL_3' $
  , xoffset = 10, yoffset = 135, xsize = 42, ysize = 18 $
  , /align_left, value = 'Offset')

  WID_LABEL_4 = widget_label(WID_BASE_1, uname = 'WID_LABEL_4' $
  , xoffset = 10, yoffset = 165, xsize = 42, ysize = 18 $
  , /align_left, value = 'Slope')

  ; ;;;;;;;;;;;;;;;;;;  table headline ;;;;;;;;;;;;;;;;;;;;;;;;
  ; ; with top labels VALUE AND LIMIT

  WID_LABEL_6 = widget_label(WID_BASE_1, uname = 'WID_LABEL_6' $
  , xoffset = 60, yoffset = 15, xsize = 42, ysize = 18 $
  , /align_left, value = 'Value')

  WID_LABEL_8 = widget_label(WID_BASE_1, uname = 'WID_LABEL_8' $
  , xoffset = 250, yoffset = 15, xsize = 40, ysize = 18 $
  , /align_left, value = 'Limits')

  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  WID_TEXT_0 = widget_text(WID_BASE_1, uname = 'WID_TEXT_0', frame = 1 $
  , xoffset = 60, yoffset = 40, value = [strtrim(string(fparm[0].value), 1)] $
  , xsize = 5, ysize = 1, /editable, /align_left) ; ;; Amplitude

  WID_TEXT_1 = widget_text(WID_BASE_1, uname = 'WID_TEXT_1', frame = 1 $
  , xoffset = 60, yoffset = 70, value = [strtrim(string(fparm[1].value), 1)] $
  , xsize = 5, ysize = 1, /editable) ; ;; X(max)

  WID_TEXT_2 = widget_text(WID_BASE_1, uname = 'WID_TEXT_2', frame = 1 $
  , xoffset = 60, yoffset = 100, value = [strtrim(string(fparm[2].value), 1)] $
  , xsize = 5, ysize = 1, /editable) ; ;; Width

  WID_TEXT_3 = widget_text(WID_BASE_1, uname = 'WID_TEXT_3', frame = 1 $
  , xoffset = 60, yoffset = 130, value = [strtrim(string(fparm[3].value), 1)] $
  , xsize = 5, ysize = 1, /editable) ; ;; Offset

  WID_TEXT_4 = widget_text(WID_BASE_1, uname = 'WID_TEXT_4', frame = 1 $
  , xoffset = 60, yoffset = 160, value = [strtrim(string(fparm[4].value), 1)] $
  , xsize = 5, ysize = 1, /editable) ; ;; Slope

  ; ;;;; Text widgets for lower limits

  WID_TEXT_5 = widget_text(WID_BASE_1, uname = 'WID_TEXT_5', frame = 1 $
  , xoffset = 210, yoffset = 40, value = [strtrim(string(fparm[0].limits(0)), 1)] $
  , xsize = 5, ysize = 1, /editable, /align_left) ; ;; Amplitude

  WID_TEXT_6 = widget_text(WID_BASE_1, uname = 'WID_TEXT_6', frame = 1 $
  , xoffset = 210, yoffset = 70, value = [strtrim(string(fparm[1].limits(0)), 1)] $
  , xsize = 5, ysize = 1, /editable) ; ;; X(max)

  WID_TEXT_7 = widget_text(WID_BASE_1, uname = 'WID_TEXT_7', frame = 1 $
  , xoffset = 210, yoffset = 100, value = [strtrim(string(fparm[2].limits(0)), 1)] $
  , xsize = 5, ysize = 1, /editable) ; ;; Width

  WID_TEXT_8 = widget_text(WID_BASE_1, uname = 'WID_TEXT_8', frame = 1 $
  , xoffset = 210, yoffset = 130, value = [strtrim(string(fparm[3].limits(0)), 1)] $
  , xsize = 5, ysize = 1, /editable) ; ;; Offset

  WID_TEXT_9 = widget_text(WID_BASE_1, uname = 'WID_TEXT_9', frame = 1 $
  , xoffset = 210, yoffset = 160, value = [strtrim(string(fparm[4].limits(0)), 1)] $
  , xsize = 5, ysize = 1, /editable) ; ;; Slope

  ; ;;;; Text widgets for upper limits

  WID_TEXT_10 = widget_text(WID_BASE_1, uname = 'WID_TEXT_10', frame = 1 $
  , xoffset = 290, yoffset = 40, value = [strtrim(string(fparm[0].limits(1)), 1)] $
  , xsize = 5, ysize = 1, /editable, /align_left) ; ;; Amplitude

  WID_TEXT_11 = widget_text(WID_BASE_1, uname = 'WID_TEXT_11', frame = 1 $
  , xoffset = 290, yoffset = 70, value = [strtrim(string(fparm[1].limits(1)), 1)] $
  , xsize = 5, ysize = 1, /editable) ; ;; X(max)

  WID_TEXT_12 = widget_text(WID_BASE_1, uname = 'WID_TEXT_12', frame = 1 $
  , xoffset = 290, yoffset = 100, value = [strtrim(string(fparm[2].limits(1)), 1)] $
  , xsize = 5, ysize = 1, /editable) ; ;; Width

  WID_TEXT_13 = widget_text(WID_BASE_1, uname = 'WID_TEXT_13', frame = 1 $
  , xoffset = 290, yoffset = 130, value = [strtrim(string(fparm[3].limits(1)), 1)] $
  , xsize = 5, ysize = 1, /editable) ; ;; Offset

  WID_TEXT_14 = widget_text(WID_BASE_1, uname = 'WID_TEXT_14', frame = 1 $
  , xoffset = 290, yoffset = 160, value = [strtrim(string(fparm[4].limits(1)), 1)] $
  , xsize = 5, ysize = 1, /editable) ; ;; Slope

  ; ;;;;;;;;;;; Container for checkboxes and checkboxes for FIXED ;;;;;;;;;;;;;;;;

  WID_BASE_2 = widget_base(WID_BASE_1, uname = 'WID_BASE_2' $
  , xoffset = 110, yoffset = 40, xsize = 20, ysize = 20 $
  , /align_top, /base_align_center, row = 1, /nonexclusive)

  WID_BUTTON_0 = widget_button(WID_BASE_2, /align_center, uvalue = 'FIX_HEIGHT' $
  , value = '')
  ; ;; Amplitude

  WID_BASE_3 = widget_base(WID_BASE_1, uname = 'WID_BASE_3' $
  , xoffset = 110, yoffset = 70, xsize = 20, ysize = 27 $
  , /align_top, /base_align_center, /nonexclusive)

  WID_BUTTON_1 = widget_button(WID_BASE_3, /align_left, uvalue = 'FIX_XMAX' $
  , value = '')
  ; ;; X(max)

  WID_BASE_4 = widget_base(WID_BASE_1 $
    , xoffset = 110, yoffset = 100, xsize = 20, ysize = 27, /align_top, /base_align_center, /nonexclusive)

  WID_BUTTON_2 = widget_button(WID_BASE_4, /align_left, uvalue = 'FIX_WIDTH' $
  , value = '')
  ; ;; Width

  WID_BASE_5 = widget_base(WID_BASE_1, uname = 'WID_BASE_5' $
  , xoffset = 110, yoffset = 130, xsize = 20, ysize = 27, /align_top, /base_align_center, /nonexclusive)

  WID_BUTTON_3 = widget_button(WID_BASE_5, /align_left, uvalue = 'FIX_OFFSET' $
  , value = '')
  ; ;; Slope

  WID_BASE_6 = widget_base(WID_BASE_1, uname = 'WID_BASE_6' $
  , xoffset = 110, yoffset = 160, xsize = 20, ysize = 27 $
  , /align_top, /base_align_center, /nonexclusive)

  WID_BUTTON_4 = widget_button(WID_BASE_6, /align_left, uvalue = 'FIX_SLOPE' $
  , value = '')
  ; ;; Slope

  ; ;;;;;;;;;;; Container for checkboxes and checkboxes for lower limited ;;;;;;;;;;;;;;;;

  WID_BASE_7 = widget_base(WID_BASE_1, uname = 'WID_BASE_7' $
  , xoffset = 180, yoffset = 40, xsize = 20, ysize = 27 $
  , /align_top, /base_align_center, title = 'IDL', space = 2, row = 1 $
  , /grid_layout, /nonexclusive)

  WID_BUTTON_5 = widget_button(WID_BASE_7 $
    , /align_left, uvalue = 'LIMIT_HEIGHT_LOW', value = '')
  ; ;; Height

  WID_BASE_8 = widget_base(WID_BASE_1, uname = 'WID_BASE_8' $
  , xoffset = 180, yoffset = 70, xsize = 20, ysize = 27 $
  , /align_top, /base_align_center, title = 'IDL', space = 2, row = 1 $
  , /grid_layout, /nonexclusive)

  WID_BUTTON_6 = widget_button(WID_BASE_8 $
    , /align_left, uvalue = 'LIMIT_XMAX_LOW', value = '')
  ; ;; Xmax

  WID_BASE_9 = widget_base(WID_BASE_1, uname = 'WID_BASE_9' $
  , xoffset = 180, yoffset = 100, xsize = 20, ysize = 27 $
  , /align_top, /base_align_center, title = 'IDL', space = 2, row = 1 $
  , /grid_layout, /nonexclusive)

  WID_BUTTON_7 = widget_button(WID_BASE_9 $
    , /align_left, uvalue = 'LIMIT_WIDTH_LOW', value = '')
  ; ;; Width

  WID_BASE_10 = widget_base(WID_BASE_1, uname = 'WID_BASE_10' $
  , xoffset = 180, yoffset = 130, xsize = 20, ysize = 27 $
  , /align_top, /base_align_center, title = 'IDL', space = 2, row = 1 $
  , /grid_layout, /nonexclusive)

  WID_BUTTON_8 = widget_button(WID_BASE_10, uname = 'WID_BUTTON_8' $
  , /align_left, uvalue = 'LIMIT_OFFSET_LOW', value = '')
  ; ;; Offset

  WID_BASE_11 = widget_base(WID_BASE_1, uname = 'WID_BASE_11' $
  , xoffset = 180, yoffset = 160, xsize = 20, ysize = 27 $
  , /align_top, /base_align_center, title = 'IDL', space = 2, row = 1 $
  , /grid_layout, /nonexclusive)

  WID_BUTTON_9 = widget_button(WID_BASE_11, uname = 'WID_BUTTON_9' $
  , /align_left, uvalue = 'LIMIT_SLOPE_LOW', value = '')
  ; ;; Offset

  ; ;;;;;;;;;;; Container for checkboxes and checkboxes for upper limited ;;;;;;;;;;;;;;;;

  WID_BASE_12 = widget_base(WID_BASE_1, uname = 'WID_BASE_12' $
  , xoffset = 265, yoffset = 40, xsize = 20, ysize = 27 $
  , /align_top, /base_align_center, title = 'IDL', space = 2, row = 1 $
  , /grid_layout, /nonexclusive)

  WID_BUTTON_10 = widget_button(WID_BASE_12, uname = 'WID_BUTTON_10' $
  , /align_left, uvalue = 'LIMIT_HEIGHT_UP', value = '')
  ; ;; Height

  WID_BASE_13 = widget_base(WID_BASE_1, uname = 'WID_BASE_13' $
  , xoffset = 265, yoffset = 70, xsize = 20, ysize = 27 $
  , /align_top, /base_align_center, title = 'IDL', space = 2, row = 1 $
  , /grid_layout, /nonexclusive)

  WID_BUTTON_11 = widget_button(WID_BASE_13, uname = 'WID_BUTTON_11' $
  , /align_left, uvalue = 'LIMIT_XMAX_UP', value = '')
  ; ;; Xmax

  WID_BASE_14 = widget_base(WID_BASE_1, uname = 'WID_BASE_14' $
  , xoffset = 265, yoffset = 100, xsize = 20, ysize = 27 $
  , /align_top, /base_align_center, title = 'IDL', space = 2, row = 1 $
  , /grid_layout, /nonexclusive)

  WID_BUTTON_12 = widget_button(WID_BASE_14, uname = 'WID_BUTTON_12' $
  , /align_left, uvalue = 'LIMIT_WIDTH_UP', value = '')
  ; ;; Width

  WID_BASE_15 = widget_base(WID_BASE_1, uname = 'WID_BASE_15' $
  , xoffset = 265, yoffset = 130, xsize = 20, ysize = 27 $
  , /align_top, /base_align_center, title = 'IDL', space = 2, row = 1 $
  , /grid_layout, /nonexclusive)

  WID_BUTTON_13 = widget_button(WID_BASE_15, uname = 'WID_BUTTON_13' $
  , /align_left, uvalue = 'LIMIT_OFFSET_UP', value = '')
  ; ;; Offset

  WID_BASE_16 = widget_base(WID_BASE_1, uname = 'WID_BASE_16' $
  , xoffset = 265, yoffset = 160, xsize = 20, ysize = 27 $
  , /align_top, /base_align_center, title = 'IDL', space = 2, row = 1 $
  , /grid_layout, /nonexclusive)

  WID_BUTTON_14 = widget_button(WID_BASE_16, uname = 'WID_BUTTON_14' $
  , /align_left, uvalue = 'LIMIT_SLOPE_UP', value = '')
  ; ;; Offset

  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  BUT_BASE = widget_base(tlb, row = 1, yoffset = 260, xoffset = 120, /align_center)
  cancel = widget_button(BUT_BASE, value = 'Cancel', uvalue = 'Cancel')
  ACCEPT = widget_button(BUT_BASE, value = 'Accept', uvalue = 'Accept')

  ; Create a pointer. This will point to the location where the
  ; information collected from the user will be stored.

  ptr = ptr_new({fparm: fparm, cancel: 1})

  ; Create info structure to hold information needed in event handler. This info structure containes also the parinfo
  ; that is used for mpfit.

  buttons = [WID_BUTTON_0, WID_BUTTON_1, WID_BUTTON_2, WID_BUTTON_3, WID_BUTTON_4, $
    WID_BUTTON_5, WID_BUTTON_6, WID_BUTTON_7, WID_BUTTON_8, WID_BUTTON_9, $
    WID_BUTTON_10, WID_BUTTON_11, WID_BUTTON_12, WID_BUTTON_13, WID_BUTTON_14, $
    cancel]
  text = [WID_TEXT_0, WID_TEXT_1, WID_TEXT_2, WID_TEXT_3, WID_TEXT_4, $
    WID_TEXT_5, WID_TEXT_6, WID_TEXT_7, WID_TEXT_8, WID_TEXT_9, $
    WID_TEXT_10, WID_TEXT_11, WID_TEXT_12, WID_TEXT_13, WID_TEXT_14]

  for i = 0, 4 do begin
    widget_control, buttons[i], set_button = fparm[i].fixed
  endfor

  for i = 5, 9 do begin
    widget_control, buttons[i], set_button = fparm[i - 5].limited(0)
  endfor

  for i = 10, 14 do begin
    widget_control, buttons[i], set_button = fparm[i - 10].limited(1)
  endfor

  info = {buttons: buttons, $ ; Identifier of widget holding buttons (checkboxes).
    text: text, $ ; Identifier of widget holding textfields.
    fparm: fparm, $ ; The actual parinfo
    ptr: ptr} ; The pointer

  ; Store the info structure in the top-level base

  widget_control, tlb, set_uvalue = info, /no_copy
  widget_control, /realize, tlb

  ; Register the program, set up event loop. Make this program a
  ; blocking widget. This will allow the program to also be called
  ; from IDL command line without a PARENT parameter. The program
  ; blocks here until the entire program is destroyed.

  xmanager, 'SetFitParm', tlb, event_handler = 'SetFitParm_Events'

  ; OK,

  newInfo = *ptr
  ptr_free, ptr

  ; All kinds of things can go wrong now. Let's CATCH them all.

  catch, error
  if error ne 0 then begin
    catch, /cancel

    ; If an error occurs, set the CANCEL flag and return -1.

    ok = dialog_message(!err_string)
    cancel = 1
    RETURN, -1
  endif

  ; If the error flag is set, let's disappear!

  cancel = newInfo.cancel
  if cancel then RETURN, fparm

  ; OK, try to read the data file. Watch out!

  RETURN, newInfo.fParm
end
; *******************************************************************