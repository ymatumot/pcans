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


function file_read, filename,format=format,string=string,silent=silent,compress=compress

if(not(keyword_set(filename)))then begin
    if(not(keyword_set(filename)))then filename=''
    filename = dialog_pickfile(filter=filename)
endif
if(filename eq '')then return,0

list = findfile(filename)

list = list(0)
if(strlen(list) eq 0) then begin
    print,'No such file'
    return,0
endif

if(keyword_set(compress))then begin
   spawn, 'gzip -cd '+filename+ '|wc -wl', input
endif else begin
   spawn, 'wc -wl '+filename, input
endelse

input = input(0)
data_info = strsplit(strcompress(strtrim(input, 2)), ' ',/extract)
if((double(data_info(1)) mod double(data_info(0))) ne 0.0 and not(keyword_set(format)))then begin
    print,'Not formatteded regularly'
    return,0
endif

if(double(data_info(0)) eq 0.0)then return,0

line = long(double(data_info(0)))
col = long(double(data_info(1))/line+0.5)

if(not(keyword_set(silent)))then begin
    print,'column:',col
    print,' line:',line
endif

if(keyword_set(string))then begin
    infile = strarr(col,line) 
endif else begin
    infile = dblarr(col,line) 
endelse

openr, /get_lun, unit, filename, compress=compress
readf, unit, infile, format=format
;if(col eq 1) then infile = reform(infile,/overwrite)


close,unit
free_lun,unit

return, infile

end
