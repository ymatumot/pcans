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
pl1=plot(time2*wce, total(t_perp,1)/total(t_para,1),$
        xstyle=1,ystyle=2,$
        xtitle='$\Omega_{ce}T$',ytitle='$T_{\perp}/T_{\parallel}$',$
        _extra=pldefop)
pl1.save,'t_tratio.png',resolution=100

;; DRAW X-T
img1=image(data/b0,xax/(c/wpe),time*wce,rgb=33,$
           xtitle='$x (c/\omega_{pe})$',ytitle='$\Omega_{ce}T$',$
          _extra=imgdefop,aspect=0)
img1.save,'x-t.png',resolution=100

;; DRAW K-T
img2 = image(zkt[drx[0]:drx[1],drt[0]:drt[1]],kx[drx[0]:drx[1]]*(c/wpe),time[drt[0]:drt[1]]*wce,$
             rgb=33,xtitle='$k_x (c/\omega_{pe})$',ytitle='$\Omega_{ce}T$',$
             _extra=imgdefop,aspect=0)
img2.save,'k-t.png',resolution=100

; DRAW T-FGM
pl2 = plot(time[drt[0]:drt[1]/4.]*wce, abs(datakt[527,drt[0]:drt[1]/4.]),$
           xtitle='$\Omega_{ce}T$',ytitle='$B_{yk}/B_0$',/ylog,$
           title='FGM $(k_xc/\omega_{pe})=0.7$',name='Simulation',$
           _extra=pldefop)
pl2.xrange=[min(time[drt[0]:drt[1]/4.]),max(time[drt[0]:drt[1]/4.])]*wce
pl2.yrange=[min(abs(datakt[527,drt[0]:drt[1]/4.])),max(abs(datakt[527,drt[0]:drt[1]/4.]))]
;; OVERPLOT LINEAR THEORY
pl2o = plot(time[drt[0]:drt[1]/4.]*wce, exp(0.27*time[drt[0]:drt[1]/4.]*wce),line=2,over=pl2,$
            name='Linear theory')
;; ADDED LEGEND
leg = legend(target=[pl2,pl2o],position=[40,5],/data,font_size=18)
pl2.save,'k-t_fgm.png',resolution=100

;; DRAW W-K
img3 = image(zkw[drx[0]:drx[1],drw[0]:drw[1]],kx[drx[0]:drx[1]]*(c/wpe),w[drw[0]:drw[1]]/wce,$
             xtitle='$k_x (c/\omega_{pe})$',ytitle='$\omega/\Omega_{ce}$',rgb=33,$
             _extra=imgdefop,aspect=0)
;; OVERPLOT LINEAR THEORY
pl3 = plot(kx[drx[0]:drx[1]]*(c/wpe),(kx[drx[0]:drx[1]]*(c/wpe))^2/(1.+(kx[drx[0]:drx[1]]*(c/wpe))^2),$
           over=img3,line=2,$
           _extra=pldefop)
img3.save,'k-w.png',resolution=100


end







