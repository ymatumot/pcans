window,xsize=600,ysize=900

c   = 1.0
va  = 2.0E-2
wgi = 2.0E-4
Mm  = 25
q   = 2.23d-3
v0  = va*10.
u0  = v0/sqrt(1.-v0^2)
L   = u0/wgi
n0  = 40
b0  = wgi*Mm/q
e0  = v0*b0/c
dt  = 2500.*0.5*wgi
ui_range = [-0.5*c,0.5*c]
ue_range = [-0.75*c,0.75*c]

ct  = 2
chars = 0.8 

pat = 'mom/*den_1.dat'
pattern = '\ls '+pat
spawn, pattern, list1
pat = 'mom/*den_2.dat'
pattern = '\ls '+pat
spawn, pattern, list2
pat = 'mom/*bz.dat'
pattern = '\ls '+pat
spawn, pattern, list3
pat = 'mom/*ex.dat'
pattern = '\ls '+pat
spawn, pattern, list4
pat = 'psd/*psd_i.dat'
pattern = '\ls '+pat
spawn, pattern, list5
pat = 'psd/*psd_e.dat'
pattern = '\ls '+pat
spawn, pattern, list6

info = size(list1)
nt   = info(1)
time = (findgen(nt)+1)*dt
dir   = 'img/'
fname = 'shock'

for i=23,nt-1 do begin

set_ps,dir+fname+strcompress(string(i,format='(I03)'),/remove)+'.eps'

print,'Reading... ',list1[i]
dni = file_read(list1[i])
print,'Reading... ',list2[i]
dne = file_read(list2[i])
print,'Reading... ',list3[i]
bz  = file_read(list3[i])
print,'Reading... ',list4[i]
ex  = file_read(list4[i])
print,'Reading... ',list5[i]
fi  = file_read(list5[i])
print,'Reading... ',list6[i]
fe  = file_read(list6[i])

info = size(dni)
nx = info(1)
x = (findgen(nx)+0.5)/L

max_n = max(dne[0:nx*9./10],xpos)
xpos = xpos/L
xrange = [xpos-0.8,xpos+0.4]

;; psd calculation
psd_calc,psdi,xax_i,yax_i,fi[0,*]/L,fi[1,*],$
         nbin_x=nx/15,nbin_y=101,max_vy=ui_range[1],min_vy=ui_range[0],max_vx=xrange[1],min_vx=xrange[0]
psd_calc,psde,xax_e,yax_e,fe[0,*]/L,fe[1,*],$
         nbin_x=nx/15,nbin_y=101,max_vy=ue_range[1],min_vy=ue_range[0],max_vx=xrange[1],min_vx=xrange[0]


plot,x,dni/n0,xs=1,ys=2,xtickname=replicate(' ',8),$
     ytitle='!6N/N!I0!N',chars=chars,pos=[0.125,0.85,0.85,0.99],xrange=xrange
oplot,x,dne/n0,col=50
legend,['!6ion','electron'],line=[0,0],color=[0,50]
xyouts,0.7,0.95,'!7X!6!Igi!NT='+strcompress(string(time[i],format='(f6.2)'),/remove),/normal
plot,x,bz/b0, xs=1,ys=2,xtickname=replicate(' ',6),$
     ytitle='!6B!Iz!N/B!I0!N',chars=chars,pos=[0.125,0.705,0.85,0.845],/noerase,xrange=xrange
plot,x,ex/e0,xs=1,ys=2,xtickname=replicate(' ',6),$
     ytitle='!6E!Ix!N/E!I0!N',chars=chars,pos=[0.125,0.56,0.85,0.700],/noerase,xrange=xrange
plot_clcnt,psdi,ct=ct,pos=[0.125,0.32,0.85,0.53],/nocolbar,/keep,/noerase,$
           xtickname=replicate(' ',6),chars=chars,xax=xax_i,yax=yax_i,$
           title='!6PSD-i',ytitle='!6u!Ix!N/c'
color_bar,psdi,ct=ct,0.88,0.32,0.90,0.53,maxz=maxz,minz=minz,charsize=charsize,/ver
plot_clcnt,psde,ct=ct,pos=[0.125,0.075,0.85,0.285],/nocolbar,/keep,/noerase,$
           chars=chars,xax=xax_i,yax=yax_e,$
           title='!6PSD-e',ytitle='!6u!Ix!N/c',$
           xtitle='!6x (U!I0!N/!7X!6!Igi!N)'
color_bar,psdi,ct=ct,0.88,0.075,0.90,0.285,maxz=maxz,minz=minz,charsize=charsize,/ver
device,/close

endfor

set_x

end

