pro psd_calc,psd,xax,yax,data1,data2,$
             max_x=max_x,min_x=min_x,max_y=max_y,min_y=min_y,nbin_x=nbin_x,nbin_y=nbin_y

if(not(keyword_set(max_x)))then max_x = max(data1)
if(not(keyword_set(max_y)))then max_y = max(data2)
if(not(keyword_set(min_x)))then min_x = min(data1)
if(not(keyword_set(min_y)))then min_y = min(data2)
if(not(keyword_set(nbin_x)))then nbin_x = 20
if(not(keyword_set(nbin_y)))then nbin_y = 20

binsize_x = (max_x-min_x)/float(nbin_x)
binsize_y = (max_y-min_y)/float(nbin_y)

psd = hist_2d(data1,data2,bin1=binsize_x,bin2=binsize_y,$
              max1=max_x,max2=max_y,min1=min_x,min2=min_y)
xax = min_x+findgen(n_elements(psd[*,0]))*binsize_x
yax = min_y+findgen(n_elements(psd[0,*]))*binsize_y

end
