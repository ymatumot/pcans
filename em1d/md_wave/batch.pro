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
img = image(z[drx[0]:drx[1],drt[0]:drt[1]],kx[drx[0]:drx[1]]*rge,w[drt[0]:drt[1]]/wge,$
           xtitle='$k_x r_{ge}$',ytitle='$\omega/\omega_{ge}$',$
           title='Parallel waves',aspect=0,rgb=33,$
           _extra=imgdefop)
img.xrange=[kx[drx[0]],kx[drx[1]]]*rge
img.yrange=[w[drt[0]],w[drt[1]]]/wge
cb = colorbar(target=img,/orientation,/textpos,font_size=16,font_name='Times',position=[0.925,0.3,0.95,0.7])

;;electromagnetic wave
y = c*kx
pl1 = plot(kx*rge,y/wge,over=img,_extra=pldefop)

;;Alfven wave
;y = kx*vai
;oplot,kx*rgi,y/wpe

;;Re-cut off frequency
y = fltarr(nx)+1.0
y = y*0.5*(wge-wgi+sqrt((wge+wgi)^2+4.*(wpe^2+wpi^2)))
pl2 = plot(kx[drx[0]:drx[1]]*rge,y/wge,line=2,over=img,_extra=pldefop)

;;Le-cut off frequency
y = fltarr(nx)+1.0
y = y*0.5*(-(wge-wgi)+sqrt((wge+wgi)^2+4.*(wpe^2+wpi^2)))
pl3 = plot(kx[drx[0]:drx[1]]*rge,y/wge,line=2,over=img,_extra=pldefop)

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
pl4 = plot(kx[drx[0]:drx[1]]*rge,y,line=2,over=img,_extra=pldefop)

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


end







