c   = 1.0
delt = 0.5
va  = 2.0e-2
wgi = 2.0e-4
Mm  = 25
q   = 4.46e-3
ma  = 10.0
v0  = va*ma
u0  = v0/sqrt(1.-v0^2)
L   = va/wgi
n0  = 10
b0  = wgi*c*Mm/q
e0  = v0*b0/c
dt  = 2000.*delt*wgi

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
dir  = './fig/'
if(not(file_test(dir)))then spawn, 'mkdir '+dir

fname = 'shock2d'
ct  = 13

for i=0,nt-1 do begin

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

info = size(dni,/dimension)
nx = info[0]
ny = info[1]
xax = (findgen(nx)+0.5)/L
yax = (findgen(ny)+0.5)/L

max_n = max(total(dne[0:nx*9/10,*],2)/ny,xpos)
xrange = [max([xpos-8.*L,0]),min([xpos+4.*L,nx-1])]
margin = [0.1,0.1,0.15,0.01]
win = window(dimension=[1200,800],/no_toolbar)

img1 = image(dni[xrange[0]:xrange[1],*]/n0,xax[xrange[0]:xrange[1]],yax,$
             axis_style=2,xtickdir=1,ytickdir=1,font_size=14,$
             title='$N_i/N_0$',$
             xtitle='$X (c/\omega_{pi})$',$
             ytitle='$Y (c/\omega_{pi})$',$
             rgb=ct,layout=[1,5,1],margin=margin,/current)
pos = img1.position
pos[0] = pos[2]+0.05
pos[1] = pos[1]
pos[2] = pos[0]+0.025
pos[3] = pos[3]
cb = colorbar(target=img1,/orientation,/textpos,position=pos)

img2 = image(dne[xrange[0]:xrange[1],*]/n0,xax[xrange[0]:xrange[1]],yax,$
             axis_style=2,xtickdir=1,ytickdir=1,font_size=14,$
             title='$N_e/N_0$',$
             xtitle='$X (c/\omega_{pi})$',$
             ytitle='$Y (c/\omega_{pi})$',$
             rgb=ct,layout=[1,5,2],margin=margin,/current)
pos = img2.position
pos[0] = pos[2]+0.05
pos[1] = pos[1]
pos[2] = pos[0]+0.025
pos[3] = pos[3]
cb = colorbar(target=img2,/orientation,/textpos,position=pos)

img3 = image(bz[xrange[0]:xrange[1],*]/b0,xax[xrange[0]:xrange[1]],yax,$
             axis_style=2,xtickdir=1,ytickdir=1,font_size=14,$
             title='$B_z/B_0$',$
             xtitle='$X (c/\omega_{pi})$',$
             ytitle='$Y (c/\omega_{pi})$',$
             rgb=ct,layout=[1,5,3],margin=margin,/current)
pos = img3.position
pos[0] = pos[2]+0.05
pos[1] = pos[1]
pos[2] = pos[0]+0.025
pos[3] = pos[3]
cb = colorbar(target=img3,/orientation,/textpos,position=pos)

img4 = image(ex[xrange[0]:xrange[1],*]/e0,xax[xrange[0]:xrange[1]],yax,$
             axis_style=2,xtickdir=1,ytickdir=1,font_size=14,$
             title='$E_x/E_0$',$
             xtitle='$X (c/\omega_{pi})$',$
             ytitle='$Y (c/\omega_{pi})$',$
             rgb=ct,layout=[1,5,4],margin=margin,/current)
pos = img4.position
pos[0] = pos[2]+0.05
pos[1] = pos[1]
pos[2] = pos[0]+0.025
pos[3] = pos[3]
cb = colorbar(target=img4,/orientation,/textpos,position=pos)

img5 = image(ey[xrange[0]:xrange[1],*]/e0,xax[xrange[0]:xrange[1]],yax,$
             axis_style=2,xtickdir=1,ytickdir=1,font_size=14,$
             title='$E_y/E_0$',$
             xtitle='$X (c/\omega_{pi})$',$
             ytitle='$Y (c/\omega_{pi})$',$
             rgb=ct,layout=[1,5,5],margin=margin,/current)
pos = img5.position
pos[0] = pos[2]+0.05
pos[1] = pos[1]
pos[2] = pos[0]+0.025
pos[3] = pos[3]
cb = colorbar(target=img5,/orientation,/textpos,position=pos)
img1.save,dir+fname+strcompress(string(i,format='(I03)'),/remove)+'.png',resolution=100

endfor


;;PSD
fi = file_read('psd/030000_0500-0036_psd_i.dat')

psd_calc,psdi,xax,yax,fi[2,*],fi[3,*],nbin_x=51,nbin_y=51,$
         max_x=0.3,min_x=-0.3,max_y=0.55,min_y=-0.1

img6 = image(alog10(psdi>0),xax,yax,rgb=13,$
             xtitle='$u_x/c$',ytitle='$u_y/c$',$
             title='$(x,y) = ($'+strcompress(string(500/L,format='(f3.1)'),/remove)$
                    +'$,1.3) c/\omega_{pi}$',$
             _extra=imgdefop)
img6.save,'psdi_foot.png',resolution=100

end

