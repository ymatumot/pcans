;; SIMULATION PARAMETERS
c   = 1.0d1
va  = 4.0E-1
wgi = 8.0E-3
Mm  = 25
q   = 3.53d-2
v0  = va*3.
u0  = v0/sqrt(1.-v0^2)
L   = va/wgi
n0  = 64
b0  = wgi*Mm*c/q
e0  = v0*b0/c
t0  = 1.0/wgi
dt  = 2000.*5d-2
vte = 1.0d0
vti = 0.1d0
ui_range = [-1,1]/2.0
ue_range = [-1,1]

pat = 'mom/*_den_i.dat'
pattern = '\ls '+pat
spawn, pattern, list1
pat = 'mom/*_den_e.dat'
pattern = '\ls '+pat
spawn, pattern, list2
pat = 'mom/*_bz.dat'
pattern = '\ls '+pat
spawn, pattern, list3
pat = 'mom/*_ex.dat'
pattern = '\ls '+pat
spawn, pattern, list4
pat = 'psd/*_psd_i.dat'
pattern = '\ls '+pat
spawn, pattern, list5
pat = 'psd/*_psd_e.dat'
pattern = '\ls '+pat
spawn, pattern, list6

info = size(list1)
nt   = info(1)
time = (findgen(nt)+1)*dt/t0
dir = './fig/'
if(not(file_test(dir)))then spawn, 'mkdir '+dir
fname = 'shock'

win = window(dimension=[600,1200],/no_toolbar)
margin = [0.15,0.25,0.15,0.1]
font_size = 14

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
fi  = file_read(list5[i])
print,'Reading... ',list6[i]
fe  = file_read(list6[i])

info = size(dni,/dimension)
nx = info[0]
xax = (findgen(nx)+0.5)/L

max_n = max(dne[0:nx*9./10],xpos)
xpos = xpos/L
xrange = [min(fi[0,*]),max(fi[0,*])]/L

;; PSD CALCULATION
psd_calc,psdi,xax_i,yax_i,fi[0,*]/L,fi[1,*]/c,$
         nbin_x=nx/15,nbin_y=101,max_y=ui_range[1],min_y=ui_range[0],max_x=xrange[1],min_x=xrange[0]
psd_calc,psde,xax_e,yax_e,fe[0,*]/L,fe[1,*]/c,$
         nbin_x=nx/15,nbin_y=101,max_y=ue_range[1],min_y=ue_range[0],max_x=xrange[1],min_x=xrange[0]

;; PLOT & IMAGE
if(i eq 0)then begin

pl1 = plot(xax,dni/n0,xstyle=1,ystyle=2,$
           xtitle='$x (\lambda_i)$',ytitle='$N/N_0$',/current,margin=margin,font_size=font_size,$
           title='$\Omega_{gi}T=$'+strcompress(string(time[i],format='(f6.2)'),/remove),$
           layout=[1,5,1])

pl2 = plot(xax,bz/b0,xstyle=1,ystyle=2,$
           xtitle='$x (\lambda_i)$',ytitle='$B_z/B_0$',/current,margin=margin,font_size=font_size,$
           layout=[1,5,2])
pl3 = plot(xax,ex/e0,xstyle=1,ystyle=2,$
           xtitle='$x (\lambda_i)$',ytitle='$E_x/E_0$',/current,margin=margin,font_size=font_size,$
           layout=[1,5,3])
img1 = image(alog10(psdi>0),xax_i,yax_i,rgb=13,axis_style=2,xtickdir=1,ytickdir=1,$
             title='PSD-i',xtitle='$x (\lambda_i)$',ytitle='$u_x/c$',/current,margin=margin,$
             layout=[1,5,4],aspect=0,font_size=font_size)
pos = img1.position
cb1 = colorbar(target=img1,/textpos,/ori,font_size=font_size*0.75,$
               position=[pos[2]+0.05,pos[1],pos[2]+0.075,pos[3]])

img2 = image(alog10(psde>0),xax_e,yax_e,rgb=13,axis_style=2,xtickdir=1,ytickdir=1,$
             title='PSD-e',xtitle='$x (\lambda_i)$',ytitle='$u_x/c$',/current,margin=margin,$
             layout=[1,5,5],aspect=0,font_size=font_size)
pos = img2.position
cb2 = colorbar(target=img2,/textpos,/ori,font_size=font_size*0.75,$
               position=[pos[2]+0.05,pos[1],pos[2]+0.075,pos[3]])

endif else begin

pl1.setdata,xax,dni/n0
pl1.title='$\Omega_{gi}T=$'+strcompress(string(time[i],format='(f6.2)'),/remove)
pl2.setdata,xax,bz/b0
pl3.setdata,xax,ex/e0
img1.setdata,alog10(psdi>0),xax_i,yax_i
img2.setdata,alog10(psde>0),xax_e,yax_e

endelse

pl1.save,dir+'shock_plot'+strcompress(string(i,format='(i03)'),/remove)+'.eps'

endfor


end

