;; paramters
n0 = 3000
vte = 2.83d-1
me = 1.0
mi = 1837.0
te = 0.5*me*n0*vte^2
wge = 4.00d-2
wpe = 2.00d-1
c = 1.0
q = 1.78d-3
b0 = wge*me*c/q
chars = 1.2

by = dir_read('../mom/*by.dat')
bz = dir_read('../mom/*bz.dat')
ve_para = dir_read('../mom/*Txx_e.dat')
ve_perp = dir_read('../mom/*Tyy_e.dat')

te_para = ve_para^2
te_perp = ve_perp^2
info=size(by)
nx = info(1)
nt = info(2)

xax = findgen(nx)/(c/wpe)
;time = findgen(nt)*50*0.5*wge
time = findgen(nt)*200*0.5*wge
be = (by^2+bz^2)/(8.*!pi)

;set_ps,'Bene-T.eps'
opps,file='Bene-T.eps'
plot,time[1:*],total(be[*,1:*],1)/(nx*te),/ylog,xs=2,ys=2,xrange=[0,200],$
     xtitle='!7X!I!6e!Nt',ytitle='!7d!NB!E2!N/(8!7p!6n!I0!Nk!IB!NT!i0!N)',$
     chars=chars,thick=2
clps

;set_ps,'Temp_ratio_T.eps'
opps,file='Temp_ratio_T.eps'
plot,time[1:*],total(Te_perp/Te_para,1)/nx,xs=2,ys=2,xrange=[0,200],yrange=[0,4],$
     xtitle='!7X!I!6e!Nt',ytitle='T!Ie_perp!N/T!Ie_para!N',chars=chars,thick=2
clps

bek = fltarr(nx,nt)
kx = 2.*!pi*findgen(nx/2+1)/nx
for l=0,nt-1 do begin
   bek[*,l] = fft(be[*,l]/te,-1)
endfor

;set_ps,'Bene-kx-T.eps'
opps,file='Bene-kx-T.eps'
plot_clcnt,abs(transpose(bek[1:nx/15,1:*],[1,0])),xrange=[0,400],$
           xax=time[1:*],yax=kx[1:nx/15]*c/wpe,ct=13,/keep,$
           xtitle='!7X!I!6e!Nt',ytitle='!6k!Ix!Nc/!7x!I!6e!N',chars=chars
clps

;set_ps,'x-bz.eps'
opps,file='x-bz.eps'
plot,xax,bz[*,20]/b0,xs=1,ys=2,thick=2,chars=chars,xrange=[0.0,400.],$
     xtitle='!6x!7x!Ie!N!6/c',ytitle='dB!Iz!N/B!I0!N',title='!7X!N!6!Ie!NT=20'
clps

end
