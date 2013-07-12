;d = file_read('mom/015000*.dat')

npos=50
info = size(d)
nx = info[1]
ny = info[2]

v0 = 0.5/4.0 ; typical Alfven speed : 1/(alpha*sqrt[me/mi])
b0 = 28.025  ; reference mag field  : manual input
e0 = v0*b0   ; typical electric field

;plot_clcnt,d[*,*,8]/28.025,ct=33
;plot_clcnt,d[*,*,8]/28.025,ct=33
ezi=d(*,*,13)+d(*,*,15)*d(*,*,7)-d(*,*,17)*d(*,*,6)
plot_clcnt,d(*,*,16)/v0,ct=33
;plot_clcnt,d(*,*,17)/v0,ct=33
;plot_clcnt,ezi/e0,ct=33

rotb_z = +0.5*(-shift(d[*,*,7],+1,0)+shift(d[*,*,7],-1,0)) $
         -0.5*(-shift(d[*,*,6],0,+1)+shift(d[*,*,6],0,-1))
rotb_z[0,*] = 0.0
rotb_z[nx-1,*] = 0.0

poisson_bp, -rotb_z, az
plot_cnt,float(az),25,c_labels=0

end
