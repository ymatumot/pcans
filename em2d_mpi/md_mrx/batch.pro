d = file_read('mom/015000*.dat')
info = size(d,/dimension)
nx = info[0]
ny = info[1]

c = 1.0
wpi = 3.95d-2
l0 = c/wpi
va = 6.25d-1
wge = 7.91d-2
q = 2.82d-3
m = 1.0d0 
b0 = wge/q*m*c
e0 = va*b0

xax = (findgen(nx)-(nx-1)/2.)/l0
yax = findgen(ny)/l0
dx = xax[1]-xax[0]
dy = yax[1]-yax[0]

img=image(d[*,*,8]/b0,xax-0.5*dx,yax-0.5*dy,_extra=imgdefop,rgb=33,$
          xtitle='$x (\lambda_i)$',ytitle='$y (\lambda_i)$',title='$B_z/B_0$')
cb = colorbar(/textpos,/orientation,font_size=16,font_name='Times')

rotb_z = +0.5*(-shift(d[*,*,7],+1,0)+shift(d[*,*,7],-1,0)) $
         -0.5*(-shift(d[*,*,6],0,+1)+shift(d[*,*,6],0,-1))
rotb_z[0,*] = 0.0
rotb_z[nx-1,*] = 0.0

poisson_bp, -rotb_z, az
cnt=contour(float(az),xax-0.5*dx,yax-0.5*dy,c_label_show=0,c_color=0,c_thick=1.5,transparency=50,$
            c_value=min(az)+randomu(seed,50)*(max(az)-min(az)),/over)

end
