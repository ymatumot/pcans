pro wave2

dat=file_read('wk_ey.dat')
x=dindgen(1001)+10000.0
ts=time_double(time_string(x))
wav=fltarr(2000)

for i=0,1000-1 do begin
wav[i]=dat[1000,i]
endfor

store_data,'wave',data={x:ts,y:wav}

stop
end
