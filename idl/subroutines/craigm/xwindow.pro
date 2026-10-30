;+
; NAME:
;       XWINDOW
;
; PURPOSE:
;       This routine implements a "smart" resizeable graphics window.
;       It is used as a wrapper for built-in IDL graphics procedures
;       such as SURFACE, CONTOUR, PLOT, SHADE_SURF, etc. In additon,
;       it can be used to display any user-written graphics procedure
;       so long as that procedure follows three simple rules: (1) It
;       does not open it's own graphics windows, (2) It is defined with
;       no more than ten positional arguments (an unlimited number
;       of keyword arguments are allowed), and (3) It is defined
;       with an _EXTRA keyword.
;
;       Keyword arguments permit the window to have its own portion
;       of a color table and to be able to change the colors loaded in
;       that portion of the color table. Colors are updated
;       automatically on both 8-bit and 24-bit color displays. In
;       addition, the window colors can "protect" themselves. I mean
;       by this that the window can re-load its own colors into the
;       color table when the cursor is moved over the window. This
;       prevents other applications from changing the colors used to
;       display data in this window. (This is an issue mainly in
;       IDL 5.x applications where widget applications can run
;       concurrently with commands from the IDL command line.)
;
;       Keyword arguments also permit the window to create output
;       files of its contents. These files can be color and
;       gray-scale PostScript, GIF, TIFF, or JPEG files.
;
; AUTHOR:
;       ************* CM 19 Jan 1999 VERSION **********
;       Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;       craigm@lheamail.gsfc.nasa.gov
;
;       Originally by:
;       FANNING SOFTWARE CONSULTING
;       David Fanning, Ph.D.
;       2642 Bradbury Court
;       Fort Collins, CO 80521 USA
;       Phone: 970-221-0438
;       E-mail: davidf@dfanning.com
;       Coyote's Guide to IDL Programming: http://www.dfanning.com
;
; CATEGORY:
;       Widgets, Graphics.
;
; CALLING SEQUENCE:
;       XWINDOW, command, P1, P2, ... , keywords=..., ...
;
; REQUIRED INPUTS:
;       COMMAND: The graphics procedure command to be executed. This parameter
;       must be a STRING. Examples are 'SURFACE', 'CONTOUR', 'PLOT', etc.
;
; OPTIONAL INPUTS:
;       Pn: A positional parameter appropriate for the graphics command.
;           Any number of parameters between 0 and 10 may be given.
;
; INPUT KEYWORD PARAMETERS:
;       CPMENU: Setting this keyword adds a "Color Protection" button to the
;       "Controls" menu. Color protection can then be turned ON or OFF for the
;       window. Otherwise, the color protection scheme used to open the window
;       cannot be changed once the window is open. (See the PROTECT keyword.)
;       The default is to have this keyword OFF.
;
;       ERASE: Setting this keyword "erases" the contents of the current
;       graphics window before re-executing the graphics command. For example,
;       this keyword might need to be set if the graphics "command" is TVSCL.
;       The default is to NOT erase the display before reissuing the graphics
;       command.
;
;       _EXTRA: This keyword forms an anonymous structure of any unrecognized
;       keywords passed to the program. The keywords must be appropriate
;       for the graphics command being executed.
;
;       GROUP_LEADER: The group leader for this program. When the group leader
;       is destroyed, this program will be destroyed.
;
;       OUTPUT: Set this keyword if you want a "File Output" menu on
;       the menu bar. The default is to have this item turned OFF.
;
;       JUST_REGISTER: If this keyword is set, the XWINDOW program is just
;       registered with XMANAGER, but XMANAGER doesn't run. This is
;       useful, for example, if you want to open an XWINDOW window in
;       the widget definition module of another widget program.
;
;       NO_CHANGE_CONFIG: Normally as the XWINDOW graphics window is resized
;       the size (or aspect ratio, in the case of PostScript) of the
;       hardware configuration dialogs change to reflect the new size of
;       the graphics window. This results in file output that resembles
;       the current graphics window in size and aspect ratio. If you want
;       the file output dialogs to remember their current configuration
;       even as the window is resized, set this keyword.
;
;       NOMENU: Setting this keyword results in a graphics window without
;       menu items. The default is to have a "Controls" menu item in the
;       window menu bar with a "Quit" button. Setting this keyword
;       automatically turns of the COLORS, OUTPUT, and CPMENU menu
;       choices. (Note that the values specified by the COLORS keyword
;       will still be valid for color protection, but no "Change Colors..."
;       menu item will appear.)
;
;       PROTECT: If this keyword is set, color protection for the draw
;       widget is turned ON. What this means is that the window colors
;       (see the XCOLOR keyword) will be restored when the cursor enters
;       the draw widget window. This prevents someone at the IDL command
;       line in IDL 5.0 from changing the window display colors permanently.
;
;       WTITLE: This is the window title. It is the string "Resizeable
;       COMMAND Window (1)" by default, where COMMAND is the input
;       parameter. And the number (1 in this case) is the window
;       index number of the draw widget.
;
;       WXPOS: This is the initial X offset of the window. Default is to
;       position the window in the approximate middle of the display.
;
;       WYPOS: This is the initial Y offset of the window. Default is to
;       position the window in the approximate middle of the display.
;
;       WXSIZE: This is the initial X size of the window. Default is 400
;       pixels.
;
;       WYSIZE: This is the initial Y size of the window. Default is 400
;       pixels.
;
;       XCOLORS: Using this keyword adds a "Change Colors..." button to the
;       "Controls" menu. Set this keyword to the number of colors available
;       in the window and the starting index of the first color. For example,
;       to allow the window access to 100 colors, starting at color index 50
;       (i.e., color indices 50 to 149), use XColor=[100, 50]. If you use the
;       keyword syntax "/XColor", all the colors available will be used, not just
;       one color. If the keyword is set to a scalar value greater than 1, the
;       starting color index is set to 0. The default value for this keyword
;       is [(!D.N_COLORS < 256), 0]. Note that color "protection" may be
;       turned on (via the PROTECT keyword) even if this keyword is NOT used.
;
; OUTPUT KEYWORD PARAMETERS:
;       DRAWID: This keyword returns the draw widget identifier of the draw
;       widget created in XWINDOW.
;
;       TOP: This keyword returns the identifier of the top-level base widget
;       created by XWINDOW.
;
;       WID: This keyword returns the window index number of the draw widget
;       created in XWINDOW.
;
; COMMON BLOCKS:
;       None.
;
; SIDE EFFECTS:
;       If color protection is ON, the window colors are reloaded when the
;       cursor enters the XWINDOW graphics windows.
;
; RESTRICTIONS: This program requires three additional programs from
;       the Fanning Software Consulting library: PSWINDOW, PS_FORM (CM
;       version; available at
;       http://astrog.physics.wisc.edu/~craigm/idl), and XCOLORS. You
;       might also want to get the program TVIMAGE if you will be
;       displaying images in XWINDOW graphics windows.
;
;       If the "command" program requires keywords that are also keywords
;       to XWINDOW, then you must use the keyword twice on the command line.
;
; EXAMPLE:
;       To display a surface in the window, type:
;
;       XWINDOW, 'SURFACE', Dist(20), Charsize=1.5
;
;       To enable the Change Colors and File Output menu items, type:
;
;       XWINDOW, 'SHADE_SURF', Dist(30), /XColors, /Output
;
; MODIFICATION HISTORY:
;       Written by: David Fanning, October 96.
;       XSIZE and YSIZE keywords changed to WXSIZE and WYSIZE so as not to
;          conflict with these keywords on other programs. 14 April 1997, DWF.
;        Updated as non-blocking widget for IDL 5.0. 14 April 1997, DWF.
;        Extensively modified to work on either 8-bit or 24-bit displays,
;          to enable color protection schemes, to send the contents to a
;          number of different output files, and to give the user choices
;          about which menu items to enable. 21 April 1997, DWF.
;        Renamed COLORS keyword to XCOLORS and fixed a problem cancelling
;           out of File Configuration dialogs. 23 April 1997, DWF.
;        Modification: Craig Markwardt, 21 October 1997
;           Added capability for up to ten positional parameters
;        Modification: CM 15 May 1998
;           PS_FORM dependencies are not hardcoded now.  Requires
;           CM version of PS_FORM function.
;        Modification: CM 19 Jan 1999
;           Added Parent keyword to widget invocation of PS_FORM, and
;           make widgets MODAL-safe.
;-

pro NULL_EVENTS, event
  compile_opt idl2
end
; of NULL_EVENTS event handler *****************************************

function XWINDOW_ALERT, message, xoffset = xoff, yoffset = yoff
  compile_opt idl2

  ; Put up a message box

  if n_params() eq 0 then message = 'Please wait...'
  device, get_screen_size = screenSize
  if n_elements(xoff) eq 0 then xoff = (screenSize[0] / 2.0 - 100)
  if n_elements(yoff) eq 0 then yoff = (screenSize[1] / 2.0 - 75)

  tlb = widget_base(title = 'Writing a File...', xoffset = xoff, yoffset = yoff)
  label = widget_label(tlb, value = message)
  widget_control, tlb, /realize
  RETURN, tlb
end
; *******************************************************************

pro XWINDOW_COLOR_PROTECTION, event
  compile_opt idl2
  widget_control, event.top, get_uvalue = info, /no_copy
  widget_control, event.id, get_uvalue = buttonValue
  case buttonValue of
    'ON': begin
      info.protect = 1
      widget_control, info.cprotectOff, sensitive = 1
      widget_control, info.cprotectOn, sensitive = 0
      widget_control, info.drawId, tracking_events = 1
    end
    'OFF': begin
      info.protect = 0
      widget_control, info.cprotectOff, sensitive = 0
      widget_control, info.cprotectOn, sensitive = 1
      widget_control, info.drawId, tracking_events = 0
    end
  endcase
  widget_control, event.top, set_uvalue = info, /no_copy
end
; *******************************************************************

pro XWINDOW_CONFIGURATION_EVENTS, event
  compile_opt idl2

  widget_control, event.top, get_uvalue = info, /no_copy
  widget_control, event.id, get_uvalue = thisEvent
  case thisEvent of
    'SELECT_FILE': begin
      ; Start in the current directory.

      cd, current = startDirectory

      ; Use PICKFILE to pick a filename for writing.

      pick = Pickfile(path = startDirectory, /noconfirm, $
        get_path = path, /write)

      ; Make sure the user didn't cancel out of PICKFILE.

      if pick ne '' then widget_control, info.filenameId, set_value = pick
    end ; of the Select Filename button case

    'CANCEL': begin
      ; Have to exit here gracefully. Set CANCEL field in structure.

      formdata = {cancel: 1, create: 0}
      Handle_Value, info.ptr, formdata, /set

      ; Out of here!

      widget_control, event.top, /destroy
      RETURN
    end ; of the Cancel button case

    'ACCEPT': begin ; Gather the form information.

      ; Get the filename.

      widget_control, info.filenameId, get_value = filename

      filename = filename[0]

      ; Get the size info.

      widget_control, info.xsizeId, get_value = xsize
      widget_control, info.ysizeId, get_value = ysize

      ; Get the color info from the droplist widget.

      listIndex = widget_info(info.colordropId, /droplist_select)
      colortype = fix(abs(1 - listIndex))

      ; Get the order info from the droplist widget.

      order = widget_info(info.orderdropId, /droplist_select)
      order = fix(order)

      ; Get the quality fromt he slider widget, if needed

      if info.sliderId ne -1 then $
        widget_control, info.sliderId, get_value = quality else quality = -1

      ; Create the formdata structure from the information you collected.

      formdata = {filename: filename, xsize: xsize, ysize: ysize, $
        color: colortype, order: order, quality: quality, create: 0}

      ; Store the formdata in the pointer location.

      Handle_Value, info.ptr, formdata, /set

      ; Out of here!

      widget_control, event.top, /destroy
      RETURN
    end ; of the Accept button case

    'CREATE': begin ; Gather the form information.

      ; Get the filename.

      widget_control, info.filenameId, get_value = filename

      filename = filename[0]

      ; Get the size info.

      widget_control, info.xsizeId, get_value = xsize
      widget_control, info.ysizeId, get_value = ysize

      ; Get the color info from the droplist widget.

      listIndex = widget_info(info.colordropId, /droplist_select)
      colortype = fix(abs(1 - listIndex))

      ; Get the order info from the droplist widget.

      order = widget_info(info.orderdropId, /droplist_select)
      order = fix(order)

      ; Get the quality fromt he slider widget, if needed

      if info.sliderId ne -1 then $
        widget_control, info.sliderId, get_value = quality else quality = -1

      ; Create the formdata structure from the information you collected.

      formdata = {filename: filename, xsize: xsize, ysize: ysize, $
        color: colortype, order: order, quality: quality, create: 1}

      ; Store the formdata in the pointer location.

      Handle_Value, info.ptr, formdata, /set

      ; Out of here!

      widget_control, event.top, /destroy
      RETURN
    end ; of the Create button case

    else:
  endcase

  widget_control, event.top, set_uvalue = info, /no_copy
end
; of XWINDOW_CONFIGURATION_EVENTS event handler ************************

function XWINDOW_CONFIGURATION, filetype, config, title = title, $
  xoffset = xoffset, yoffset = yoffset, cancel = cancel, create = create, $
  parent = parent
  compile_opt idl2

  catch, error
  if error ne 0 then begin
    ok = WIDGET_MESSAGE(!err_string)
    RETURN, -1
  endif

  if n_elements(filetype) eq 0 then filetype = 'GIF'
  if n_elements(config) eq 0 then config = {xsize: 400, ysize: 400, $
    color: 1, filename: 'xwindow.gif', ncolors: (!d.n_colors < 256)}
  filetype = strupcase(filetype)
  if n_elements(title) eq 0 then title = 'Configure ' + $
    filetype + ' Output File'

  ; Check for placement offsets. Define defaults.

  if (n_elements(xoffset) eq 0) then begin
    device, get_screen_size = screenSize
    xoffset = (screenSize[0] - 200) / 2.
  endif
  if (n_elements(yoffset) eq 0) then begin
    device, get_screen_size = screenSize
    yoffset = (screenSize[1] - 100) / 2.
  endif

  ; Create widgets.

  thisRelease = strmid(!version.release, 0, 1)
  if thisRelease eq '5' and n_elements(parent) gt 0 then $
    extra_modal = {modal: 1, group_leader: parent[0]}
  tlb = widget_base(column = 1, title = title, xoffset = xoffset, $
    yoffset = yoffset, base_align_center = 1, _extra = extra_modal)

  bigbox = widget_base(tlb, column = 1, frame = 1, base_align_center = 1)

  ; Create the filename widgets.
  filebox = widget_base(bigbox, column = 1, base_align_center = 1)
  filename = config.filename
  filenamebase = widget_base(filebox, row = 1)
  filenamelabel = widget_label(filenamebase, value = 'Filename:')
  filenameID = widget_text(filenamebase, value = filename, /editable, $
    event_pro = 'NULL_EVENTS', scr_xsize = 320)

  ; Create a button to allow user to pick a filename.

  pickbutton = widget_button(filebox, value = 'Select Filename', $
    uvalue = 'SELECT_FILE')

  ; Create size widgets
  sizebox = widget_base(bigbox, column = 1, base_align_left = 1)
  sizebase = widget_base(sizebox, row = 1)
  xsizeID = cw_field(sizebase, value = config.xsize, title = 'XSize: ', $
    /integer)
  ysizeID = cw_field(sizebase, value = config.ysize, title = 'YSize: ', $
    /integer)

  ; File type and order.

  orderbase = widget_base(sizebox, row = 1)
  type = ['Color', 'Grayscale']
  order = ['0', '1']
  colordropID = widget_droplist(orderbase, value = type, $
    title = 'File Type: ', event_pro = 'NULL_EVENTS')
  orderdropID = widget_droplist(orderbase, value = order, $
    title = 'Display Order: ', event_pro = 'NULL_EVENTS')

  widget_control, colordropID, set_droplist_select = fix(abs(config.color - 1))
  widget_control, orderdropID, set_droplist_select = config.order

  ; Quality Slider if needed.

  if filetype eq 'JPEG' then $
    sliderID = widget_slider(bigbox, value = config.quality, max = 100, min = 0, $
      title = 'Compression Quality', event_pro = 'NULL_EVENTS', $
      scr_xsize = 350) else sliderID = -1

  ; Cancel and Accept buttons.

  buttonbase = widget_base(tlb, row = 1)
  cancelID = widget_button(buttonbase, value = 'Cancel', uvalue = 'CANCEL')
  createID = widget_button(buttonbase, value = 'Create File', uvalue = 'CREATE')
  ok = widget_button(buttonbase, value = 'Accept', uvalue = 'ACCEPT')

  widget_control, tlb, /realize

  ptr = HANDLE_CREATE()

  info = {filenameId: filenameID, xsizeId: xsizeID, $
    ysizeId: ysizeID, colordropId: colordropID, $
    orderdropId: orderdropID, ptr: ptr, sliderId: sliderID}

  widget_control, tlb, set_uvalue = info, /no_copy
  if thisRelease eq '4' then xmanager_modal = {modal: 1}
  xmanager, 'xwindow_configuration', tlb, $
    event_handler = 'XWINDOW_CONFIGURATION_EVENTS', _extra = xmanager_modal

  Handle_Value, ptr, formdata
  Handle_Free, ptr

  if n_elements(formdata) eq 0 then begin
    cancel = 1
    create = 0
    RETURN, -1
  endif

  fields = tag_names(formdata)
  create = formdata.create
  cancel = where(fields eq 'CANCEL')
  if cancel[0] eq -1 then begin
    cancel = 0
    newConfiguration = create_struct('XSIZE', formdata.xsize, $
      'YSIZE', formdata.ysize, 'COLOR', formdata.color, $
      'FILENAME', formdata.filename, 'ORDER', formdata.order, $
      'QUALITY', formdata.quality, name = 'XWINDOW_' + filetype)
    RETURN, newConfiguration
  endif else begin
    cancel = 1
    create = 0
    RETURN, -1
  endelse
end
; of XWINDOW_CONFIGURATION event handler *******************************

function XWindow_WhatTypeVariable, variable
  compile_opt idl2

  ; Use SIZE function to get variable info.

  varInfo = size(variable)

  ; The next to last element in varInfo has the data type.

  typeIndex = varInfo[varInfo[0] + 1]
  dataTypes = ['UNDEFINED', 'BYTE', 'INTEGER', 'LONG', 'FLOATING', $
    'DOUBLE', 'COMPLEX', 'STRING', 'STRUCTURE', 'DCOMPLEX']
  thisType = dataTypes[typeIndex]

  RETURN, thisType
end
; of XWindow_WhatTypeVariable utility routine **************************

pro XWINDOW_QUIT, event
  compile_opt idl2
  widget_control, event.top, /destroy
end
; of XWINDOW_QUIT procedure *********************************************

pro XWINDOW_CLEANUP, id
  compile_opt idl2
  widget_control, id, get_uvalue = info, /no_copy
  if n_elements(info) eq 0 then RETURN
  HANDLE_FREE, info.plotObjPtr
end
; of XWINDOW_CLEANUP cleanup procedure ***********************************

pro XWINDOW_CONFIGURE_FILES, event
  compile_opt idl2
  widget_control, event.top, get_uvalue = info, /no_copy

  ; What kind of file to configure?

  widget_control, event.id, get_value = whichFile
  case whichFile of
    'Configure PostScript File...': begin
      newkeywords = ps_form(defaults = info.ps, parent = info.top, $
        localdefaults = info.pslocal, cancel = cancel, create = create)
      if not cancel then info.ps = newkeywords
      if create then widget_control, info.psId, send_event = {id: info.psId, $
        top: event.top, handler: 0l}
    end

    'Configure GIF File...': begin
      config = info.gif
      newConfiguration = XWINDOW_CONFIGURATION('GIF', config, $
        cancel = cancel, create = create, parent = info.top)
      if not cancel then info.gif = newConfiguration
      if create then widget_control, info.gifId, send_event = {id: info.gifId, $
        top: event.top, handler: 0l}
    end

    'Configure TIFF File...': begin
      config = info.tiff
      newConfiguration = XWINDOW_CONFIGURATION('TIFF', config, $
        cancel = cancel, create = create, parent = info.top)
      if not cancel then info.tiff = newConfiguration
      if create then widget_control, info.tiffId, send_event = {id: info.tiffId, $
        top: event.top, handler: 0l}
    end

    'Configure JPEG File...': begin
      config = info.jpeg
      newConfiguration = XWINDOW_CONFIGURATION('JPEG', config, $
        cancel = cancel, create = create, parent = info.top)
      if not cancel then info.jpeg = newConfiguration
      if create then widget_control, info.jpegId, send_event = {id: info.jpegId, $
        top: event.top, handler: 0l}
    end
  endcase

  widget_control, event.top, set_uvalue = info, /no_copy
end
; of XWINDOW_CONFIGURE_FILES event handler ***********************************

pro XWINDOW_CREATE_FILES, event
  compile_opt idl2
  widget_control, event.top, get_uvalue = info, /no_copy

  ; There can be all kinds of problems writing a file.
  ; Trap errors here and try to get out of here.

  catch, error
  if error ne 0 then begin
    junk = widget_message(!err_string)
    ok = WIDGET_MESSAGE(['Problem writing file. The most', $
      'common problem is a mis-spelled filename.', 'Returning...'])
    if widget_info(id, /valid_id) then widget_control, id, /destroy
    if n_elements(thisDevice) gt 0 then set_plot, thisDevice
    if n_elements(info) ne 0 then widget_control, event.top, $
      set_uvalue = info, /no_copy
    RETURN
  endif

  id = XWINDOW_ALERT('Please be patient while writing a file...')

  ; Get the Plot Object.

  HANDLE_VALUE, info.plotObjPtr, plotObj, /no_copy

  ; What kind of file to create?

  widget_control, event.id, get_value = whichFile
  case whichFile of
    'Create PostScript File': begin
      keywords = info.ps
      thisDevice = !d.name
      tvlct, r, g, b, /get
      set_plot, 'PS'
      tvlct, info.r, info.g, info.b, info.bottom
      device, _extra = keywords
      ok = execute(plotObj.thisCommand)
      device, /close_file
      set_plot, thisDevice
      tvlct, r, g, b
    end

    'Create GIF File': begin
      config = info.gif

      ; Render graphic in Z-buffer.

      thisDevice = !d.name
      tvlct, rr, gg, bb, /get
      set_plot, 'Z'
      erase
      device, set_resolution = [config.xsize, config.ysize], $
        set_colors = info.ncolors
      ok = execute(plotObj.thisCommand)
      thisImage = tvrd()
      if config.color ne 1 then loadct, 0, ncolors = info.wcolors, $
        bottom = info.bottom else $
        tvlct, info.r, info.g, info.b, info.bottom
      tvlct, r, g, b, /get
      set_plot, thisDevice
      tvlct, rr, gg, bb

      ; Write GIF file.

      write_gif, config.filename, thisImage, r, g, b
    end ; of GIF file creation.

    'Create TIFF File': begin
      config = info.tiff

      ; Render graphic in Z-buffer.

      thisDevice = !d.name
      tvlct, rr, gg, bb, /get
      set_plot, 'Z'
      tvlct, info.r, info.g, info.b, info.bottom
      erase
      device, set_resolution = [config.xsize, config.ysize], $
        set_colors = info.ncolors
      ok = execute(plotObj.thisCommand)
      thisImage = tvrd()
      tvlct, r, g, b, /get
      set_plot, thisDevice
      tvlct, rr, gg, bb

      ; Write TIFF file. Use screen resolution.

      if config.color eq 1 then $
        TIFF_WRITE, config.filename, thisImage, config.order, $
        red = r, green = g, blue = b, xresol = round(!d.x_px_cm * 2.54), $
        yresol = round(!d.x_px_cm * 2.54) else $
        TIFF_WRITE, config.filename, thisImage, config.order, $
        xresol = round(!d.x_px_cm * 2.54), yresol = round(!d.x_px_cm * 2.54)
    end

    'Create JPEG File': begin
      config = info.jpeg

      ; Render graphic in Z-buffer.

      thisDevice = !d.name
      tvlct, rr, gg, bb, /get
      set_plot, 'Z'
      erase
      device, set_resolution = [config.xsize, config.ysize], $
        set_colors = info.ncolors
      tvlct, info.r, info.g, info.b, info.bottom
      ok = execute(plotObj.thisCommand)
      thisImage = tvrd()
      tvlct, r, g, b, /get
      set_plot, thisDevice
      tvlct, rr, gg, bb

      ; Write JPEG file.

      if config.color eq 1 then begin
        image24 = bytarr(3, config.xsize, config.ysize)
        image24[0, *, *] = r[thisImage]
        image24[1, *, *] = g[thisImage]
        image24[2, *, *] = b[thisImage]
        write_jpeg, config.filename, image24, true = 1, $
          quality = config.quality, order = config.order
      endif else $
        write_jpeg, config.filename, thisImage, $
        quality = config.quality, order = config.order
    end
  endcase

  ; Put the Plot Object back.

  HANDLE_VALUE, info.plotObjPtr, plotObj, /set, /no_copy
  widget_control, id, /destroy
  widget_control, event.top, set_uvalue = info, /no_copy
end
; of XWINDOW_CREATE_FILES event handler ***********************************

pro XWINDOW_COLORS, event
  compile_opt idl2
  widget_control, event.top, get_uvalue = info, /no_copy
  widget_control, event.id, get_uvalue = colors

  XCOLORS, group = event.top, ncolors = colors[0], bottom = colors[1], $
    title = 'Window ' + strtrim(info.wid, 2) + ' Colors', $
    notifyid = [info.drawId, event.top]

  widget_control, event.top, set_uvalue = info, /no_copy
end
; of XWINDOW_COLORS event handler ****************************************

pro XWINDOW_DRAW_EVENT, event
  compile_opt idl2
  widget_control, event.top, get_uvalue = info, /no_copy

  ; Get the Plot Object.

  HANDLE_VALUE, info.plotObjPtr, plotObj, /no_copy

  ; Need to respond to WIDGET_TRACKING events and
  ; XCOLORS_LOAD events.

  thisEvent = tag_names(event, /structure)

  if thisEvent eq 'WIDGET_TRACKING' then begin
    if event.enter eq 1 then $
      tvlct, info.r, info.g, info.b, info.bottom
  endif

  if thisEvent eq 'XCOLORS_LOAD' then begin
    info.r = event.r[info.bottom : info.bottom + info.wcolors - 1]
    info.g = event.g[info.bottom : info.bottom + info.wcolors - 1]
    info.b = event.b[info.bottom : info.bottom + info.wcolors - 1]
  endif

  ; Redisplay the command in the window if needed.

  ncolors = !d.n_colors
  if ncolors gt 256 then begin
    wset, info.wid
    if info.erase then erase
    ok = execute(plotObj.thisCommand)
  endif

  ; Put the Plot Object back.

  HANDLE_VALUE, info.plotObjPtr, plotObj, /set, /no_copy

  widget_control, event.top, set_uvalue = info, /no_copy
end
; of XWINDOW_DRAW_EVENT event handler **********************************

pro XWINDOW_RESIZE_EVENTS, event
  compile_opt idl2
  widget_control, event.top, get_uvalue = info, /no_copy

  ; Resize the draw widget.

  widget_control, info.drawid, xsize = event.x, ysize = event.y

  ; Get the Plot Object.

  HANDLE_VALUE, info.plotObjPtr, plotObj, /no_copy

  ; Redisplay the command in the window.

  widget_control, info.drawId, get_value = wid
  wset, wid
  if info.erase eq 1 then erase
  ok = execute(plotObj.thisCommand)

  ; Update file output configuration structures in necessary.

  if (info.output) and (not info.nochange) then begin
    info.gif.xsize = event.x
    info.gif.ysize = event.y
    info.tiff.xsize = event.x
    info.tiff.ysize = event.y
    info.jpeg.xsize = event.x
    info.jpeg.ysize = event.y
    if info.ps.inches eq 0 then newsizes = PSWINDOW(/cm) else $
      newsizes = PSWINDOW()
    info.ps.xsize = newsizes.xsize
    info.ps.ysize = newsizes.ysize
    info.ps.xoff = newsizes.xoffset
    info.ps.yoff = newsizes.yoffset
  endif

  ; Put the Plot Object back.

  HANDLE_VALUE, info.plotObjPtr, plotObj, /set, /no_copy

  widget_control, event.top, set_uvalue = info, /no_copy
end
; of XWINDOW_RESIZE_EVENTS event handler *********************************

pro XWINDOW, proName, $
  param1, param2, param3, param4, param5, $
  param6, param7, param8, param9, param10, $
  group_leader = group, $
  _extra = extra, wxsize = xsize, wysize = ysize, wid = wid, xcolors = colors, $
  drawid = drawid, wtitle = wtitle, just_register = justRegister, $
  output = output, no_change_config = nochange, erase = erase, $
  top = tlb, protect = protect, nomenu = nomenu, cpmenu = cpmenu, $
  wxpos = wxpos, wypos = wypos
  compile_opt idl2

  ; Return to main-level on error.

  on_error, 1

  ; Check keywords.

  if n_elements(xsize) eq 0 then xsize = 400
  if n_elements(ysize) eq 0 then ysize = 400
  nochange = keyword_set(nochange)
  extraFlag = keyword_set(extra)
  justRegister = keyword_set(justRegister)
  protect = keyword_set(protect)
  needColors = keyword_set(colors)
  needOutput = keyword_set(output)
  cpmenu = keyword_set(cpmenu)
  nomenu = keyword_set(nomenu)
  if nomenu then begin
    needOutput = 0
    cpmenu = 0
    cprotectON = -1l
    cprotectOFF = -1l
    psID = -1l
    gifID = -1l
    tiffID = -1l
    jpegID = -1l
  endif

  ; Make sure a window has been opened.

  thisWindowID = !d.window
  Window, xsize = 10, ysize = 10, /free, /pixmap
  wdelete, !d.window
  if thisWindowID ge 0 then wset, thisWindowID

  ; Set up color variables. If the user typed "/Colors"
  ; then use *all* colors!

  if needColors then begin
    if n_elements(colors) eq 1l and (colors[0] eq 1) then $
      colors = fix([(!d.n_colors < 256), 0])
    if (n_elements(colors) eq 1l) then colors = [colors[0] < 256l, 0]
  endif else colors = [(!d.n_colors < 256l), 0]

  colors = fix([(colors[0] < 256), colors[1]])
  wcolors = colors[0]
  bottom = colors[1]
  tvlct, r, g, b, /get
  r = r[bottom : bottom + wcolors - 1]
  g = g[bottom : bottom + wcolors - 1]
  b = b[bottom : bottom + wcolors - 1]
  if nomenu then needColors = 0

  ; Check for positional parameters. One parameter required.

  np = n_params()

  if np eq 0 then message, 'Sorry, at least one argument is required.'

  ; The first positional argument must be a string.

  if np gt 0 then begin
    thisType = XWindow_WhatTypeVariable(proName)
    if thisType ne 'STRING' then $
      message, 'First argument must be STRING type. Returning...'
  endif

  if n_elements(wtitle) eq 0 then $
    wtitle = 'Resizeable ' + strupcase(proName) + ' Window'

  ; Set up the Plot Object based on number of parameters and
  ; the extraFlag variable.

  if np eq 0 then $
    message, 'Must call XWINDOW with at least one argument.'

  if np gt 10 then $
    message, 'Must have fewer than 10 parameters in call to XWINDOW.'

  ; Compose the command
  thisCommand = proName
  for i = 2, np do $
    thisCommand = thisCommand + ', plotObj.param' + strtrim(i - 1, 2)

  ; If any extra keywords are present, put them on the command line
  if extraFlag then $
    thisCommand = thisCommand + ', _EXTRA=plotObj.extra' else $
    extra = 0

  ; We need to make sure that all of the parameters are set, at least to 0,
  ; so that they can be entered into the plotobj structure.
  if not keyword_set(param1) then param1 = 0
  if not keyword_set(param2) then param2 = 0
  if not keyword_set(param3) then param3 = 0
  if not keyword_set(param4) then param4 = 0
  if not keyword_set(param5) then param5 = 0
  if not keyword_set(param6) then param6 = 0
  if not keyword_set(param7) then param7 = 0
  if not keyword_set(param8) then param8 = 0
  if not keyword_set(param9) then param9 = 0
  if not keyword_set(param10) then param10 = 0

  plotObj = {thisCommand: thisCommand, extra: extra, $
    param1: param1, param2: param2, param3: param3, param4: param4, $
    param5: param5, param6: param6, param7: param7, param8: param8, $
    param9: param9, param10: param10}

  ; Store the Plot Object at a pointer location.

  plotObjPtr = HANDLE_CREATE()

  ; Create the widgets for this program.

  device, get_screen_size = screenSize
  if n_elements(wxpos) eq 0 then wxpos = (screenSize[0] - xsize) / 2.
  if n_elements(wypos) eq 0 then wypos = (screenSize[1] - ysize) / 2.
  if not nomenu then begin
    tlb = widget_base(tlb_size_events = 1, $
      xoffset = wxpos, yoffset = wypos, mbar = menubase)

    controls = widget_button(menubase, value = 'Controls', /menu)

    ; Need a COLORS button?

    if needColors then begin
      colorsID = widget_button(controls, value = 'Change Colors...', $
        event_pro = 'XWINDOW_COLORS', uvalue = colors)
    endif

    ; Need color protection buttons?

    if cpmenu then begin
      cprotect = widget_button(controls, value = 'Color Protection', $
        /menu, event_pro = 'XWindow_Color_Protection')
      cprotectON = widget_button(cprotect, value = 'ON', uvalue = 'ON')
      cprotectOFF = widget_button(cprotect, value = 'OFF', uvalue = 'OFF')
    endif else begin
      cprotectON = -1l
      cprotectOFF = -1l
    endelse

    ; Need FILE OUTPUT button?

    if needOutput then begin
      outputButton = widget_button(menubase, value = 'File Output', $
        /menu, event_pro = 'XWindow_Create_Files')
      psID = widget_button(outputButton, value = 'Create PostScript File')
      gifID = widget_button(outputButton, value = 'Create GIF File')
      tiffID = widget_button(outputButton, value = 'Create TIFF File')
      jpegID = widget_button(outputButton, value = 'Create JPEG File')
      configure = widget_button(outputButton, value = 'Configure Output File', $
        /menu, /separator, event_pro = 'XWindow_Configure_Files')
      ps_config = widget_button(configure, value = 'Configure PostScript File...')
      gif_config = widget_button(configure, value = 'Configure GIF File...')
      tiff_config = widget_button(configure, value = 'Configure TIFF File...')
      jpeg_config = widget_button(configure, value = 'Configure JPEG File...')
    endif else begin
      psID = -1l
      gifID = -1l
      tiffID = -1l
      jpegID = -1l
    endelse
    quit = widget_button(controls, value = 'Quit', event_pro = 'XWindow_Quit')
  endif else tlb = widget_base(tlb_size_events = 1, $
    xoffset = wxpos, yoffset = wypos)

  drawid = widget_draw(tlb, xsize = xsize, ysize = ysize, $
    event_pro = 'XWindow_Draw_Event', tracking_events = protect)

  widget_control, tlb, /realize
  widget_control, drawid, get_value = wid
  wset, wid

  ; Give each window a unique title.

  wtitle = wtitle + ' (' + strtrim(wid, 2) + ')'
  widget_control, tlb, tlb_set_title = wtitle

  ; Set color protection button sensitivity.

  if cpmenu then begin
    if protect then begin
      widget_control, cprotectON, sensitive = 0
      widget_control, cprotectOFF, sensitive = 1
    endif else begin
      widget_control, cprotectON, sensitive = 1
      widget_control, cprotectOFF, sensitive = 0
    endelse
  endif

  ; If something goes wrong executing the command, trap it.

  catch, error
  if error ne 0 then begin
    ok = WIDGET_MESSAGE(['There is a problem executing the command', $
      'string in XWINDOW. Please check keyword', $
      'spelling and command syntax. Returning...'])
    HANDLE_FREE, plotObjPtr
    widget_control, tlb, /destroy
    RETURN
  endif

  ok = execute(plotObj.thisCommand)
  if not ok then begin
    ok = WIDGET_MESSAGE(['There is a problem executing the command', $
      'string in XWINDOW. Please check keyword', $
      'spelling and command syntax. Returning...'])
    HANDLE_FREE, plotObjPtr
    widget_control, tlb, /destroy
    RETURN
  endif

  catch, /cancel

  ; Store the Plot Object in its pointer.

  HANDLE_VALUE, plotObjPtr, plotObj, /set, /no_copy

  ; Create an info structure.

  info = {top: tlb, $ ; Top level widget
    xsize: xsize, $ ; X size of window.
    ysize: ysize, $ ; Y size of window.
    wid: wid, $ ; Window index number.
    drawId: drawid, $ ; Draw widget identifier.
    cprotectOn: cprotectON, $ ; Color protection ON button.
    cprotectOff: cprotectOFF, $ ; Color protection OFF button.
    wtitle: wtitle, $ ; Window title.
    r: r, $ ; Red colors in window.
    g: g, $ ; Green colors in window.
    b: b, $ ; Blue colors in window.
    wcolors: wcolors, $ ; Number of window colors.
    gifId: gifID, $ ; ID of Create GIF file button.
    tiffId: tiffID, $ ; ID of Create TIFF file button.
    jpegId: jpegID, $ ; ID of Create JPEG file button.
    psId: psID, $ ; ID of Create PS file button.
    bottom: bottom, $ ; Starting color index.
    protect: protect, $ ; Protect colors flag.
    nomenu: nomenu, $ ; No menu flag.
    nochange: nochange, $ ; No change flag.
    erase: keyword_set(erase), $ ; Need erasure flag.
    ncolors: (!d.n_colors < 256), $ ; Size of color table.
    plotObjPtr: plotObjPtr, $ ; Pointer to plot object.
    output: needOutput} ; File Output menu flag.

  ; File Output configuration structures, if needed.

  if keyword_set(output) then begin
    cd, current = thisDir
    ps = ps_form(/init, filename = filepath(root_dir = thisDir, 'xwindow.ps'))
    pslocal = ps_form(/init, filename = filepath(root_dir = thisDir, 'xwindow.ps'), $
      xsize = 10., xoff = 0.5, ysize = 7.5, yoff = 0.5, color = 1, $
      landscape = 1)
    gif = {Xwindow_Gif, xsize: 400, ysize: 400, color: 1, $
      filename: filepath(root_dir = thisDir, 'xwindow.gif'), $
      order: 0, quality: -1}
    jpeg = {Xwindow_Jpeg, xsize: 400, ysize: 400, color: 1, $
      filename: filepath(root_dir = thisDir, 'xwindow.jpg'), $
      order: 0, quality: 75}
    tiff = {Xwindow_Tiff, xsize: 400, ysize: 400, color: 1, $
      filename: filepath(root_dir = thisDir, 'xwindow.tif'), $
      order: 1, quality: -1}
    info = create_struct(info, 'PS', ps, 'PSLOCAL', pslocal, 'GIF', gif, $
      'JPEG', jpeg, 'TIFF', tiff)
  endif

  ; Store the info structure in the TLB.

  widget_control, tlb, set_uvalue = info, /no_copy

  ; Register the program as on-blocking in 5.0.

  thisRelease = strmid(!version.release, 0, 1)
  if thisRelease eq '5' then $
    xmanager, 'xwindow', tlb, event_handler = 'XWINDOW_RESIZE_EVENTS', $
    cleanup = 'XWINDOW_CLEANUP', group_leader = group, $
    just_reg = justRegister, /no_block else $
    xmanager, 'xwindow', tlb, event_handler = 'XWINDOW_RESIZE_EVENTS', $
    cleanup = 'XWINDOW_CLEANUP', group_leader = group, just_reg = justRegister
end