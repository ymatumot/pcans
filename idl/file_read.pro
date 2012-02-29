;###########################################################################
; AUTHOR
;
;   Yosuke Matsumoto
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

if(not(keyword_set(filename)))then return,0
flist = file_search(filename,count=count)

if(count eq 0) then begin
    print,'No such file'
    return,0
endif

if(keyword_set(compress))then begin
   spawn, 'gzip -cd '+flist[0]+ '|wc -wl', input
endif else begin
   spawn, 'wc -wl '+flist[0], input
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
   tmp = strarr(col,line)
   infile = strarr(col,line,count)
endif else begin
   tmp = dblarr(col,line)
   infile = dblarr(col,line,count)
endelse


for l=0,count-1 do begin
   print, l, '.  ', flist[l], '  Reading......'
   openr, /get_lun, unit, flist[l], compress=compress
   readf, unit, tmp, format=format
   infile[*,*,l] = tmp

   close,unit
   free_lun,unit
endfor

return, reform(infile)

end
