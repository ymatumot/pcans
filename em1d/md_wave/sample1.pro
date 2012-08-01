;; DATA INFORMATION
dx = 1.0
dt = 6.0 ;; TIME INTERVAL OF DATA OUTPUT
c = 1.0D0
mr = 1.D0/16.D0
vai = 1.25D-1
vae = vai/sqrt(mr)
wgi = 2.47D-3
wge = wgi/mr
wpi = c/vai*wgi
wpe = c/vae*wge
beta = 0.05
rgi = c/wpi*sqrt(beta)
rge = rgi*sqrt(mr)
vte = rge*wge
ke = 1./(vte/wpe)

info = size(input)
nx = info(1)
nt = info(2)
x = findgen(nx)*dx
time = findgen(nt)*dt
wnum = fltarr(nx)
freq = fltarr(nt)

; 2D-FFT
output = fft2d(input,x,time,-1,wnum=wnum,freq=freq)

;; SETUP FOR W-K DIAGRAM
z = alog10(abs(output))
kx = wnum*2.*!pi
w = freq*2.*!pi

;; DRAW RANGE
drx = [nx/2,0.55*nx]
drt = [nt/2,0.75*nt]

;; DRAW
plot_clcnt,z[drx[0]:drx[1],drt[0]:drt[1]],ct=33,$
           xax=kx[drx[0]:drx[1]]*rge,yax=w[drt[0]:drt[1]]/wge,$
           xtitle='!6k!Ix!N r!Ige!N',ytitle='!7x/x!I!6ge!N',$
           title='!6Parallel waves',/ver_,/keep,chars=1.5

;;OVERPLOT FREQUENCIES
loadct,0,/silent

;;electromagnetic wave
y = c*kx
oplot,kx*rge,y/wge

;;Alfven wave
;y = kx*vai
;oplot,kx*rgi,y/wpe

;;Re-cut off frequency
y = fltarr(nx)+1.0
y = y*0.5*(wge-wgi+sqrt((wge+wgi)^2+4.*(wpe^2+wpi^2)))
oplot,kx[drx[0]:drx[1]]*rge,y/wge,col=0,line=2

;;Le-cut off frequency
y = fltarr(nx)+1.0
y = y*0.5*(-(wge-wgi)+sqrt((wge+wgi)^2+4.*(wpe^2+wpi^2)))
oplot,kx[drx[0]:drx[1]]*rge,y/wge,col=0,line=2

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
oplot,kx[drx[0]:drx[1]]*rge,y,line=2,col=0

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







