;; DATA INFORMATION
dx = 1.0
dt = 0.5*20.
c = 1.0D0
mr = 1.D0/16.D0
vai = 1.25D-1
vae = vai/sqrt(mr)
wgi = 2.21D-3
wge = wgi/mr
wpi = c/vai*wgi
wpe = c/vae*wge
beta = 0.04
rgi = c/wpi*sqrt(beta)
rge = rgi*sqrt(mr)
vte = rge*wge
ke = 1./(vte/wpe)

info = size(input)
col = info(1)
num = info(2)

output = dcomplexarr(col,num)

;; FFT
;; for k=0,num-1 do begin
;;     output(*,k) = fft(input(*,k),-1)
;; endfor
;; for i=0,col-1 do begin
;;     output(i,*) = fft(output(i,*),-1)
;; endfor
output = fft(input,-1)

;; DRAW W-K DIAGRAM
hcol = col/25
hnum = num/2.5

z = abs(output(0:hcol,0:hnum))
z = alog10(z)
;z = z^0.2

maxz = max(z)
minz = min(z)

kx = 2.*!pi*(findgen(hcol+1))/(col*dx)
w = 2.*!pi*(findgen(hnum+1))/(num*dt)

plot_clcnt,z[1:hcol,1:hnum],ct=33,xax=kx[1:hcol]*rge,yax=w[1:hnum]/wge,xtitle='!6k!Ix!N r!Ige!N',ytitle='!7x/x!I!6ge!N',title='!6Parallel waves',/keep,/ver_

;;******* add angular frequencies ********;;

loadct,0,/silent
;;electromagnetic wave
;; y = c*kx
;; oplot,kx/ke,y/wpe
;; y = 0.7*c*kx
;; oplot,kx/ke,y/wpe

;;Alfven wave
;y = kx*vai
;oplot,kx*rgi,y/wpe

;;Re-cut off frequency
y = fltarr(col)+1.0
y = y*0.5*(wge-wgi+sqrt((wge+wgi)^2+4.*(wpe^2+wpi^2)))
oplot,kx[1:hcol]*rge,y/wge,col=0,line=2

;;Le-cut off frequency
y = fltarr(col)+1.0
y = y*0.5*(-(wge-wgi)+sqrt((wge+wgi)^2+4.*(wpe^2+wpi^2)))
oplot,kx[1:hcol]*rge,y/wge,col=0,line=2

;;Upper hybrid frequency
;; y = kx
;; y[*] = 1.0
;; y = y*sqrt(wge^2+wpe^2)
;; oplot,kx*rge,y/wge

;;Lower hybrid frequency
;y[*] = 1.0
;y = y*sqrt( (wpi^2+wgi^2)/(1+wpe^2/wge^2) )
;oplot,kx*rgi,y/wpe

;;Electron gylo-frequency
y[*] = 1.0
oplot,kx[1:hcol]*rge,y,line=2,col=0

;;Ion gylo-frequency
;; y[*] = 1.0
;; y = y*wgi
;; oplot,kx,y

;;Electron plasma-frequency
;; y = sqrt(wpe^2*(1+mr)+3.*kx^2*vte^2/2.)
;; oplot,kx/ke,y/wpe

;;Ion plasma-frequency
;y[*] = 1.0
;y = x*vai*sqrt(beta)*sqrt(gam/2.)
;oplot,x,y
;y[*] = 1.0
;y = y*wpi
;oplot,x,y
;;******* end ********;;

;; RESET
loadct,12,/silent


end







