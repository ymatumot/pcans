;d = file_read('mom/015000*.dat')

npos=50
info = size(d)
nx = info[1]
ny = info[2]

v0 = 0.5/4.0 ; typical Alfven speed
b0 = 28.025  ; reference mag field
e0 = v0*b0   ; typical electric field

;plot_clcnt,d[*,*,8]/28.025,ct=33
;plot_clcnt,d[*,*,8]/28.025,ct=33
ezi=d(*,*,13)+d(*,*,15)*d(*,*,7)-d(*,*,17)*d(*,*,6)
plot_clcnt,d(*,*,16)/v0,ct=33
;plot_clcnt,d(*,*,17)/v0,ct=33
;plot_clcnt,ezi/e0,ct=33

;; STARTPOINTS
r1 = fltarr(2,npos)
r1[0,*] = randomu(seed,npos)*(nx-1)
r1[1,*] = ny/2

fl1=field_lines_2d(d[*,*,6],d[*,*,7],r0=r1,npos=100,dir=+1)
fl2=field_lines_2d(d[*,*,6],d[*,*,7],r0=r1,npos=100,dir=-1)

for k=0,npos-1 do begin
   oplot,fl1[0,*,k],fl1[1,*,k],col=0
   oplot,fl2[0,*,k],fl2[1,*,k],col=0
endfor

end
