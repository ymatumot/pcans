;###########################################################################
; AUTHOR
;
;   Yosuke Matsumoto
;   E-mail: ymatumot@ybb.ne.jp
;
;
; LICENSE
;
; This software is OSI Certified Open Source Software.
; OSI Certified is a certification mark of the Open Source Initiative.
;
; Copyright 2003 Yosuke Matsumoto
;
; This software is provided "as-is", without any express or
; implied warranty. In no event will the authors be held liable
; for any damages arising from the use of this software.
;
; Permission is granted to anyone to use this software for any
; purpose, including commercial applications, and to alter it and
; redistribute it freely, subject to the following restrictions:
;
; 1. The origin of this software must not be misrepresented; you must
;    not claim you wrote the original software. If you use this software
;    in a product, an acknowledgment in the product documentation
;    would be appreciated, but is not required.
;
; 2. Altered source versions must be plainly marked as such, and must
;    not be misrepresented as being the original software.
;
; 3. This notice may not be removed or altered from any source distribution.
;
; For more information on Open Source Software, visit the Open Source
; web site: http://www.opensource.org.
;
;###########################################################################


pro plot_clcnt, z,xax=xax,yax=yax,ct=ct,sym=sym,max=maxz,min=minz,$
                keep_aspect=keep_aspect,pos=pos,nocolbar=nocolbar,$
                noerase=noerase,_Extra=_Extra,charsize=charsize
common wdw, px,py,sx,sy

;; CHECK INPUT VARIABLE
if(not(keyword_set(z)))then begin
    z = file_read()
    if(n_elements(z) eq 1)then begin
        if(z eq 0)then return
    endif
endif


;; CHECK COLOR TABLE
if(not(keyword_set(ct)))then ct=0
if(ct gt 40) then begin
    print,'No such color table.'
    return
endif


;; READ THE SIZE OF THE MATRIX
z = reform(z)
info_z = size(z)
n_col = info_z(1)
n_line = info_z(2)
xsz = float(n_col)/float(n_col+n_line)*1280
ysz = float(n_line)/float(n_col+n_line)*1280
;; xsz = n_col*1.5
;; ysz = n_line*1.5
if(!d.window eq -1 and !d.name eq 'X')then window
if(!d.name ne 'PS'  and !d.name ne 'Z' and not keyword_set(keep_aspect))then $
  window,!d.window,xsize=xsz,ysize=ysz


;; SET X AND Y AXIS IF NOT
if(not(keyword_set(xax)))then xax = findgen(n_col)
if(not(keyword_set(yax)))then yax = findgen(n_line)
;; ADJUST X and Y AXIS TO PIXELS
xax2 = fltarr(n_col+1)
xax2(0) = 2.*xax(0)-xax(1)
xax2(1:n_col) = xax
yax2 = fltarr(n_line+1)
yax2(0) = 2.*yax(0)-yax(1)
yax2(1:n_line) = yax
xax2=xax2+(max(xax2)-min(xax2))/n_col*0.5
yax2=yax2+(max(yax2)-min(yax2))/n_line*0.5

if(not(keyword_set(noerase)))then erase

;; ADD COLOR BAR
if(not(keyword_set(nocolbar)))then begin

    if(keyword_set(pos) and n_elements(pos) eq 4)then begin
        xs=pos(0)
        xe=pos(1)
        ys=pos(2)-0.1
        ye=ys+0.01
    endif else begin
        xs = 0.15	
        xe = 0.92
        ys = 0.06
        ye = 0.08
    endelse

    if(keyword_set(sym) and min(z) lt 0.0) then begin
        color_bar,z,ct=ct,xs,xe,ys,ye,max=maxz,min=minz,/sym,charsize=charsize
    endif else begin
        color_bar,z,ct=ct,xs,xe,ys,ye,max=maxz,min=minz,charsize=charsize
    endelse

endif

;; SET THE SIZE OF THE WINDOW
if(keyword_set(pos) and n_elements(pos) eq 4)then begin
    xs=pos(0)
    xe=pos(1)
    ys=pos(2)
    ye=pos(3)
endif else begin
    if(keyword_set(nocolbar))then begin
        xs = 0.15
        xe = 0.92
        ys = 0.06
        ye = 0.95
    endif else begin
        xs = 0.15
        xe = 0.92
        ys = 0.18
        ye = 0.95
    endelse
endelse

set_window, xs, xe, ys, ye

;; DRAW
if( not(keyword_set(maxz)) )then maxz = max(z)
if( not(keyword_set(minz)) )then minz = min(z)
if (keyword_set(sym) and min(z) lt 0.0) then begin
    if (maxz gt abs(minz)) then minz = -maxz
    if (abs(minz) gt maxz) then maxz = -minz
endif

loadct, ct, /silent
if(!d.name eq 'PS') then begin
    tv,bytscl(z,max=maxz,min=minz),px(0)+1,py(0)+1,xsize=sx,ysize=sy
endif else begin
    tv, congrid(bytscl(z,max=maxz, min=minz),sx,sy), px(0)+1, py(0)+1
endelse

if( not(keyword_set(charsize)) )then charsize = 1.0
loadct,12,/silent
contour,congrid(z,n_col+1,n_line+1),xax2,yax2,$
  /nodata,position=[px(0),py(0),px(1),py(1)],$
  xs=1,ys=1,/noerase,/device,xticklen=-0.01,yticklen=-0.01,charsize=charsize,$
  _Extra=_Extra


;; RESET
loadct, 12, /silent


end





