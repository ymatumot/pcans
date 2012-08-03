;; DATA INFORMATION
dx = 1.0
dt = 5.0 ;; TIME INTERVAL OF DATA OUTPUT
c = 1.0D0
q = 2.82d-3
mr = 1.D0/1837.D0
vai = 4.67D-3
vae = vai/sqrt(mr)
wgi = 1.54D-5
wce = wgi/mr
wpi = c/vai*wgi
wpe = c/vae*wce
beta = 1.0
rgi = c/wpi*sqrt(beta)
rge = rgi*sqrt(mr)
vte = rge*wce
b0 = wce/q
ke = 1./(vte/wpe)

;; READING WAVE DATA
data = file_read('dat/wk_by.dat')
info = size(data)
nx = info(1)
nt = info(2)
xax = findgen(nx)*dx
time = findgen(nt)*dt
wnum = fltarr(nx)
freq = fltarr(nt)
datakt = complexarr(nx,nt)

; READING TEMPS
t_para = file_read('mom/*_Txx_e.dat')
t_perp = file_read('mom/*_Tyy_e.dat')
info = size(t_para)
nt2 = info(2)
time2 = findgen(nt2)*250

; 1D-FFT IN X
for k=0,nt-1 do begin
   datakt[*,k] = fftf(data[*,k]/b0,xax,-1)
endfor
; 2D-FFT
datakw = fft2d(data/b0,xax,time,-1,wnum=wnum,freq=freq)

;; SETUP FOR W-K DIAGRAM
zkt = alog10(abs(datakt))
zkw = alog10(abs(datakw))

kx = wnum*2.*!pi
w = freq*2.*!pi

;; DRAW RANGES
drx = [0.5*nx,0.54*nx]
drt = [0,nt-1]
drw = [0.5*nt,0.525*nt]

;; DRAW T-Temp
window,0
plot,time2*wce, total(t_perp,1)/total(t_para,1),$
     xs=1,ys=2,$
     xtitle='!7X!I!6ce!NT',ytitle='!6T!Iperp!N/T!Ipara!N',$
     chars=2.0
write_png,'t_tratio.png',tvrd(true=1)

;; DRAW X-T
window,1
plot_clcnt,data/b0,xax=xax/(c/wpe),yax=time*wce,ct=33,$
           xtitle='!6x (c/!7x!6!Ipe!N)',ytitle='!7X!I!6ce!NT',$
           /ver_,chars=2.0,/keep
write_png,'x-t.png',tvrd(true=1)

;; DRAW K-T
window,2
plot_clcnt,zkt[drx[0]:drx[1],drt[0]:drt[1]],ct=33,$
           xax=kx[drx[0]:drx[1]]*(c/wpe),yax=time[drt[0]:drt[1]]*wce,$
           xtitle='!6k!Ix!N (c/!7x!6!Ipe!N)',ytitle='!7X!I!6ce!NT',$
           /ver_,/keep,chars=2.0
write_png,'k-t.png',tvrd(true=1)

;; DRAW T-FGM
window,3
plot, time[drt[0]:drt[1]/4.]*wce, abs(datakt[527,drt[0]:drt[1]/4.]), xs=1, ys=1, $
      xtitle='!7X!I!6ce!NT',ytitle='!6B!Iyk!N/B!I0!N',/ylog,chars=2.0,$
      title='!6FGM (k!Ix!Nc/!7x!6!Ipe!N)=0.7'
;; OVERPLOT LINEAR THEORY
oplot, time[drt[0]:drt[1]/4.]*wce, exp(0.27*time[drt[0]:drt[1]/4.]*wce),line=2
legend,['Simulation','Linear theory'],line=[0,2],position=[30,10],/data,chars=2.0
write_png,'k-t_fgm.png',tvrd(true=1)

;; DRAW W-K
window,5
plot_clcnt,zkw[drx[0]:drx[1],drw[0]:drw[1]],ct=33,$
           xax=kx[drx[0]:drx[1]]*(c/wpe),yax=w[drw[0]:drw[1]]/wce,$
           xtitle='!6k!Ix!N (c/!7x!6!Ipe!N)',ytitle='!7x/X!I!6ce!N',$
           /ver_,/keep,chars=2.0
;; OVERPLOT LINEAR THEORY
oplot,kx[drx[0]:drx[1]]*(c/wpe),$
      (kx[drx[0]:drx[1]]*(c/wpe))^2/(1.+(kx[drx[0]:drx[1]]*(c/wpe))^2),$
      line=2,thick=2
write_png,'k-w.png',tvrd(true=1)

;; RESET
loadct,12,/silent


end







