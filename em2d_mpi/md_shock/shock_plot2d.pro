window,xsize=600,ysize=900

c   = 1.0
va  = 1.0e-2
wgi = 5.0e-5
Mm  = 100
q   = 2.23d-3
v0  = va*10.
u0  = v0/sqrt(1.-v0^2)
L   = va/wgi
n0  = 40
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

;; for i=37,nt-1 do begin
for i=0,nt-1 do begin

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
           pos=[0.1,0.87,0.95,0.97],/keep,ct=ct,$
           xax=xax[xrange[0]:xrange[1]],yax=yax

plot_clcnt,dne[xrange[0]:xrange[1],*]/n0,title='!6N!Ie!N/N!I0!N',$
           xtitle='!6 X c/!7x!I!6pi!N',$
           ytitle='!6 Y c/!7x!I!6pi!N',chars=chars,$
           pos=[0.1,0.67,0.95,0.77],/keep,ct=ct,$
           xax=xax[xrange[0]:xrange[1]],yax=yax,/noerase

plot_clcnt,bz[xrange[0]:xrange[1],*]/b0,title='!6B!Iz!N/B!I0!N',$
           xtitle='!6 X c/!7x!I!6pi!N',$
           ytitle='!6 Y c/!7x!I!6pi!N',chars=chars,$
           pos=[0.1,0.47,0.95,0.57],/keep,ct=ct,$
           xax=xax[xrange[0]:xrange[1]],yax=yax,/noerase

plot_clcnt,ex[xrange[0]:xrange[1],*]/e0,title='!6E!Ix!N/E!I0!N',$
           xtitle='!6 X c/!7x!I!6pi!N',$
           ytitle='!6 Y c/!7x!I!6pi!N',chars=chars,$
           pos=[0.1,0.27,0.95,0.37],/keep,ct=ct,$
           xax=xax[xrange[0]:xrange[1]],yax=yax,/noerase

plot_clcnt,ey[xrange[0]:xrange[1],*]/e0,title='!6E!Iy!N/E!I0!N',$
           xtitle='!6 X c/!7x!I!6pi!N',$
           ytitle='!6 Y c/!7x!I!6pi!N',chars=chars,$
           pos=[0.1,0.07,0.95,0.17],/keep,ct=ct,$
           xax=xax[xrange[0]:xrange[1]],yax=yax,/noerase

device,/close

endfor

x

end

