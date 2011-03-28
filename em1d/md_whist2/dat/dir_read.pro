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


function dir_read, pat, gui

if(not(keyword_set(pat)) or keyword_set(gui))then begin
    if(not(keyword_set(pat)))then pat=''
    list = dialog_pickfile(filter=pat,/multiple_files,get_path=cd)
    if(list(0) eq '') then return,0
endif else begin
    pattern = '\ls '+pat
    spawn, pattern, list 
    if(list(0) eq '') then return,0
endelse

info = size(list)
n = info(1)

data = file_read(list(0))
info = size(data)
n_col = info(1)
n_line = info(2)
if(info(0) eq 1)then n_line=1

output = dblarr(n_col,n_line,n)


for i=0,n-1 do begin
    file = list(i)
    if(keyword_set(cd))then begin
        file=strmid(file,strlen(cd),strlen(file))
    endif

    print, i, '.  ', file, '  Reading......'
    output(*,*,i) = file_read(file)
endfor

return, reform(output)


end

