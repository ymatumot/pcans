pro psd_calc,psd,xax,yax,data1,data2,$
             max_vx=max_vx,min_vx=min_vx,max_vy=max_vy,min_vy=min_vy,nbin_x=nbin_x,nbin_y=nbin_y

if(not(keyword_set(max_vx)))then max_vx = max(data1)
if(not(keyword_set(max_vy)))then max_vy = max(data2)
if(not(keyword_set(min_vx)))then min_vx = min(data1)
if(not(keyword_set(min_vy)))then min_vy = min(data2)
if(not(keyword_set(nbin_x)))then nbin_x = 20
if(not(keyword_set(nbin_y)))then nbin_y = 20

binsize_vx = (max_vx-min_vx)/double(nbin_x)
binsize_vy = (max_vy-min_vy)/double(nbin_y)

xax = min_vx+(findgen(nbin_x+1)+0.5)*binsize_vx
yax = min_vy+(findgen(nbin_y+1)+0.5)*binsize_vy

psd = hist_2d(data1,data2,bin1=binsize_vx,bin2=binsize_vy,$
              max1=max_vx,max2=max_vy,min1=min_vx,min2=min_vy)

end
