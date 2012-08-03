window,xsize=600,ysize=900

c   = 1.0
va  = 2.0e-2
wgi = 2.0e-4
Mm  = 25
q   = 3.15d-3
v0  = va*10.
u0  = v0/sqrt(1.-v0^2)
L   = va/wgi
n0  = 20
b0  = wgi*Mm/q
e0  = v0*b0/c
dt  = 10000.*wgi

ct  = 13
chars = 0.75 

pat = 'mom/*den_i.dat'
pattern = '\ls '+pat
spawn, pattern, list1
pat = 'mom/*den_e.dat'
pattern = '\ls '+pat
spawn, pattern, list2
pat = 'mom/*bz.dat'
pattern = '\ls '+pat
spawn, pattern, list3
pat = 'mom/*ex.dat'
pattern = '\ls '+pat
spawn, pattern, list4
pat = 'mom/*ey.dat'
pattern = '\ls '+pat
spawn, pattern, list5

info = size(list1)
nt   = info(1)
time = (findgen(nt)+1)*dt
dir   = 'img/'
fname = 'shock2d'

;; for i=0,nt-1 do begin
for i=0,-1 do begin

print,'Reading... ',list1[i]
dni = file_read(list1[i])
print,'Reading... ',list2[i]
dne = file_read(list2[i])
print,'Reading... ',list3[i]
bz  = file_read(list3[i])
print,'Reading... ',list4[i]
ex  = file_read(list4[i])
print,'Reading... ',list5[i]
ey  = file_read(list5[i])

info = size(dni)
nx = info[1]
ny = info[2]
xax = (findgen(nx)+0.5)/L
yax = (findgen(ny)+0.5)/L

max_n = max(dne[0:nx*9/10,ny/2],xpos)
xpos = xpos
xrange = [max([xpos-8.*L,0]),min([xpos+5.*L,nx-1])]

plot_clcnt,dni[xrange[0]:xrange[1],*]/n0,title='!6N!Ii!N/N!I0!N',$
           xtitle='!6 X c/!7x!I!6pi!N',$
           ytitle='!6 Y c/!7x!I!6pi!N',chars=chars,$
           pos=[0.1,0.83,0.85,0.95],/keep,ct=ct,$
           xax=xax[xrange[0]:xrange[1]],yax=yax
color_bar,dni[xrange[0]:xrange[1],*]/n0,ct=ct,0.88,0.83,0.90,0.95,charsize=chars,/ver

plot_clcnt,dne[xrange[0]:xrange[1],*]/n0,title='!6N!Ie!N/N!I0!N',$
           xtitle='!6 X c/!7x!I!6pi!N',$
           ytitle='!6 Y c/!7x!I!6pi!N',chars=chars,$
           pos=[0.1,0.64,0.85,0.76],/keep,ct=ct,$
           xax=xax[xrange[0]:xrange[1]],yax=yax,/noerase
color_bar,dne[xrange[0]:xrange[1],*]/n0,ct=ct,0.88,0.64,0.90,0.76,charsize=chars,/ver

plot_clcnt,bz[xrange[0]:xrange[1],*]/b0,title='!6B!Iz!N/B!I0!N',$
           xtitle='!6 X c/!7x!I!6pi!N',$
           ytitle='!6 Y c/!7x!I!6pi!N',chars=chars,$
           pos=[0.1,0.45,0.85,0.57],/keep,ct=ct,$
           xax=xax[xrange[0]:xrange[1]],yax=yax,/noerase
color_bar,bz[xrange[0]:xrange[1],*]/b0,ct=ct,0.88,0.45,0.90,0.57,charsize=chars,/ver

plot_clcnt,ex[xrange[0]:xrange[1],*]/e0,title='!6E!Ix!N/E!I0!N',$
           xtitle='!6 X c/!7x!I!6pi!N',$
           ytitle='!6 Y c/!7x!I!6pi!N',chars=chars,$
           pos=[0.1,0.26,0.85,0.38],/keep,ct=ct,$
           xax=xax[xrange[0]:xrange[1]],yax=yax,/noerase
color_bar,ex[xrange[0]:xrange[1],*]/e0,ct=ct,0.88,0.26,0.90,0.38,charsize=chars,/ver

plot_clcnt,ey[xrange[0]:xrange[1],*]/e0,title='!6E!Iy!N/E!I0!N',$
           xtitle='!6 X c/!7x!I!6pi!N',$
           ytitle='!6 Y c/!7x!I!6pi!N',chars=chars,$
           pos=[0.1,0.07,0.85,0.19],/keep,ct=ct,$
           xax=xax[xrange[0]:xrange[1]],yax=yax,/noerase
color_bar,ey[xrange[0]:xrange[1],*]/e0,ct=ct,0.88,0.07,0.90,0.19,charsize=chars,/ver

write_png,dir+fname+strcompress(string(i,format='(I03)'),/remove)+'.png',tvrd(true=1)

endfor

;;PSD
fi = file_read('psd/030000_0600-0128_psd_i.dat')

psd_calc,psdi,xax,yax,fi[2,*],fi[3,*],nbin_x=51,nbin_y=51,$
         max_vx=0.3,min_vx=-0.3,max_vy=0.4,min_vy=-0.1

plot_clcnt,alog10(psdi),xax=xax,yax=yax,ct=13,min=-0.05,$
           xtitle='!6u!Ix!N/c!N',ytitle='!6u!Iy!N/c!N',/ver,$
           chars=1.75,title='!6(x,y) = (6.0,1.3) c/!7x!I!6pi!N'
write_png,'psdi_foot.png',tvrd(true=1)

end

