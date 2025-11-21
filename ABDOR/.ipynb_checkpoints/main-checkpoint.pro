files = file_search('Dat_LC_files/01*/pn*b80*.dat');all the lightcurves
nf    = n_elements(files); total number of the lightcurves

print,nf
print,files

n=0

file  = files[n]
data  = read_ascii(file,data_start = 5)
time  = reform(data.field1(0,*))
flux  = reform(data.field1(1,*))
error = reform(data.field1(2,*))

dim  = size(flux,/dimensions)
line = strmid(file,27,30)
dim  = dim[0]
time = time/60. ; time in minutes
print,line

;*******************************************************
; Selecting post-flare light curve
;*******************************************************
window,1
cgplot,time,flux,/xstyle,/ystyle,title=line,xtitle='Time [minutes]',ytitle='Counts $\tex s^{-1}$'

a     = 0  ;read,'From which exposure should I look for flare peak, ',a
a1    = where(flux(a:*) eq max(flux(a:*))) ; a1 index of peak in flare lc
;a1   = where(flux(a:200) eq max(flux(a:200)))
;a1   = 150
a1    = a1(0)+a
aa1   = 150 ;read,'number of exposure to be included ..',aa1
b1    = a1+aa1 ; b1 index of quiet part of lc after a1
time1 = time(a1:b1)
tim   = time1-time1(0)
flx1  = flux(a1:b1)

rnav  = 1 ;		;read,'Choose the smoothing window (ex: 30): ',rnav
flx   = smooth(flx1,rnav,/EDGE_TRUNCATE)

cgoplot,time(a1:b1),flux(a1:b1),color=cgcolor('green');Post-flare phase of lc

;*******************************************************
;Exponantial fit to calculate the exponential decay time
;*******************************************************
;err_in    = error(a1:b1)
;start_exp = [max(flux(a:*)),0.5]
;expr_exp  = 'P[0]*(exp(-(X)/P[1]))'
;b         = 15
;;read,'Where should I stop exponential fit..', b
;res_exp      = mpfitexpr(expr_exp,tim(0:b),flux(a1:a1+b),err_in,start_exp,perror = exp_err)
;res_exp_flux = res_exp[0]*exp(-(tim)/res_exp[1])
;print,res_exp
;cgplot,tim,flx,/xstyle,/ystyle,$
;title=line,xtitle='Time [minutes]',ytitle='Counts $\tex s^{-1}$',$
;charsize=3.,thick=4.,charthick=3
;cgoplot,tim,res_exp_flux,color=cgcolor('red'),thick=5

;*******************************************************
; Emd to find the IMFs and Back ground trend  
;*******************************************************
f_emd          = emd(flx)
dim_emd        = size(f_emd,/dimensions)
bg             = f_emd(*,dim_emd(1)-1);+f_emd(*,dim_emd(1)-2) ; Back ground trend
trend_sub_flux = flx-bg ; trend subtracted post-flare phase of lc

window,2

cgplot,tim,flx,/xstyle,/ystyle,$
title=line,xtitle='Time [minutes]',ytitle='Counts $\tex s^{-1}$',$
charsize=3.,thick=4.,charthick=3

cgoplot,tim,bg,color=cgcolor('gold'),thick = 5

;*******************************************************
; Find the best correlated IMF to trend subtracted post-flare lc
;*******************************************************

cv = fltarr(dim_emd(1))
;   .r 
for i = 0, dim_emd(1)-1 do begin
  cv(i) = correlate(trend_sub_flux,reform(f_emd(*,i)),/covariance) ; correlation coefficient finding
  print, i, cv(i)
endfor
;   end

imf_n = where(cv eq max(cv))
;   print,imf_n

window,3
cgplot,tim,trend_sub_flux,/xstyle,$
title=line,xtitle='Time [minutes]',ytitle='Counts $\tex s^{-1}$',$
charsize=3.,thick=4.,charthick=3

cgoplot,tim,f_emd(*,imf_n),color=cgcolor('red'),thick=5.

;*******************************************************
; Fitting Damping Sinusoidal Function to IMF that best correlates to trend subtracted lc
;*******************************************************

err_in       = .5
start        = [1.3,40.,0.,0.5]
expr_fit     = 'P[0]*sin(2*!pi*X/P[1]-P[2])*(exp(-(X)/P[3]))'
res_fit      = mpfitexpr(expr_fit,tim,f_emd(*,imf_n),err_in,start,perror                  = fit_err)
res_fit_flux = res_fit[0]*sin((2*!pi*tim)/res_fit[1]-res_fit[2])*(exp(-(tim)/res_fit[3]))
print,res_fit

cgoplot,tim,res_fit_flux,color=cgcolor('blue'),thick=5.

period       = res_fit[1]
damping_time = res_fit[3]


; ***************************************************************
; Caution
; only for cases where selecting all of the post-flare lc gives
; spurious results then we reduce the lc for the ananlysis
; ***************************************************************
dd = size(flx,/dimensions)

if (dd(0) le 100) then d=dd(0)-1 else d=100

f_emd_f          = emd(flx(0:d))
dim_emd_f        = size(f_emd_f,/dimensions)
trend_sub_flux_f = flx(0:d)-bg(0:d)
cv_f             = fltarr(dim_emd_f(1))

; .r
for i = 0, dim_emd_f(1)-1 do begin
  print,i, correlate(trend_sub_flux_f,reform(f_emd_f(*,i)),/covariance)
  cv_f(i) = correlate(trend_sub_flux_f,reform(f_emd_f(*,i)),/covariance)
endfor
; end

imf_n_f = where(cv_f eq max(cv_f))
ti      = tim(0:d)

window,4
cgplot,ti,trend_sub_flux_f,/xstyle,title=line,xtitle='time in minutes',ytitle='Counts/S'
cgoplot,ti,f_emd_f(*,imf_n_f),color=cgcolor('red')
err_in         = .5
start_f        = [1.3,40.,0.,0.5]
expr_fit_f     = 'P[0]*sin(2*!pi*X/P[1]-P[2])*(exp(-(X)/P[3]))'
res_fit_f      = mpfitexpr(expr_fit_f,ti,f_emd_f(*,imf_n_f),err_in,start_f)
res_fit_flux_f = res_fit_f[0]*sin((2*!pi*ti)/res_fit_f[1]-res_fit_f[2])*(exp(-(tim)/res_fit_f[3]))
print,res_fit_f
cgoplot,ti,res_fit_flux_f,color=cgcolor('blue');,linestyle=2

;'**********************************************************'
;'Making hard copies of the plots'
;'**********************************************************'

filename=strcompress('emd'+'.eps',/remove_all)   

cgPS_Open, filename
!P.Multi = [1, 1, 3, 2]
cgText, 0.5, 0.99, Alignment=0.5, /Normal, 'rgs-ccd lc for data set 27'

cgplot,time,flux,/xs,/ys,title='Original Light Curve',xtitle='time in minutes',ytitle='Counts/S',Position=[0.10, 0.10, 0.45, 0.90]
cgoplot,time(a1:b1),flux(a1:b1),color=cgcolor('green');LF

cgplot,tim,flx,/xs,/ys,title='Post-flare phase with back ground',xtitle='time in minutes',ytitle='Counts/S',Position=[0.55, 0.10, 0.95, 0.60], /NoErase
cgoplot,tim,bg,color=cgcolor('gold')
cgoplot,tim,res_exp_flux,color=cgcolor('Dark Red')
cglegend,colors=[cgcolor('gold'),cgcolor('Dark Red')],symsize=1.5,location=[0.64,0.5],titles=['Back ground', 'exp fit'],length=0.05,vspace=2.,/background,bg_color='white'

cgplot,tim,trend_sub_flux,/xs,/ys,title='Trend subtracted lc',xtitle='time in minutes',ytitle='Counts/S',Position=[0.55, 0.70, 0.95, 0.90], /NoErase
cgoplot,tim,f_emd(*,imf_n),color=cgcolor('red')
cgoplot,tim,res_fit_flux,color=cgcolor('blue');,linestyle=2
;*******************************************************************
cgplot,tim,trend_sub_flux,/xs,/ys,title='',xtitle='time in minutes',ytitle='Counts/S',Position=[0.10,0.10,0.45,0.9]
cgoplot,tim,f_emd(*,imf_n),color=cgcolor('red')
cgoplot,tim,res_fit_flux,color=cgcolor('blue');,linestyle=2

tau  = where(tim gt 2*fix(res_fit[3]))
ttau = tau(0)

if ttau eq -1 then ttau=n_elements(tim)-1
cgplot,tim(0:ttau)/res_fit[3],f_emd(0:ttau,imf_n),/xs,/ys,title='',xtitle='Damping time '+ cgsymbol('tau'),ytitle='Count Rate',Position=[0.55,0.10,0.9,0.9],color='red'
cgoplot,tim(0:ttau)/res_fit[3],res_fit_flux(0:ttau),color=cgcolor('blue')
;cglegend,colors=[cgcolor('red'),cgcolor('blue')],symsize=1.5,location=[0.64,0.5],titles=['best fit IMF', 'fit to IMF'],length=0.05,vspace=2.,/background,bg_color='white'
!P.Multi = 0    
cgPS_close

;'**********************************************************'
;'Wavelet analysis'
;'**********************************************************'

rep:
t = tim*60.; time in seconds
f = trend_sub_flux

for p = 0, 1 do begin

  date  = 'date'
  line  = 'flux'
  dtype = 'Flux in Normalized Units'
  dt    = dtype
  x_pos = 100
  y_pos = 100

  window,11,xs=600,ys=300
  plot,t/60.,f,title=line,xtitle='time(min)',ytitle=dtype,/xstyle

  wlet='n'
  read,'Like to see the wavelet analysis (y/n)', wlet
  ;if wlet eq 'y' then gen_wavelet,t,f,date,line,y_pos,dt
  if wlet eq 'y' then randomlet_1,t,f,date,line,y_pos

  loc='n'
  read,'Want to check again at the same location?<default n>:',loc

  if loc eq 'y' then goto,rep
  print,'Yloc of the chosen pix= ',y_pos
  read,'Like to choose other location then enter 0; ',p

endfor

if !d.window eq 20 then wdelete,20
if !d.window eq 19 then wdelete,19 
if !d.window eq 19 then wdelete,1 

read,'p1=',p1
read,'p2=',p2

;**********************************************************************
; saving the results of the analysis
;**********************************************************************
save,flux,time,a1,b1,res_fit_flux,res_fit,res_exp_flux,tim,flx,bg,trend_sub_flux,f_emd,res_exp,fit_err,exp_err,p1,p2,n,a,b,imf_n,aa1,rnav,filename='results.sav'


;**********************************************************************
; Printing the results of the analysis
;**********************************************************************
print,'**********************************************************'
print,'Period of QPP detected from IMF in minutes                    = ', res_fit[1]
print,'Error in Period of QPP detected from IMF in minutes           = ', fit_err[1]
print,'Damping Time of QPP detected from IMF in minutes              = ', res_fit[3]
print,'Error in Damping Time of QPP detected from IMF in minutes     = ', fit_err[3]
print,'Exponential Damping Time of post-flare LC in minutes          = ', res_exp[1]
print,'Error in Exponential Damping Time of post-flare LC in minutes = ', exp_err[1]
print,'p1                                                            = ', p1
print,'p2                                                            = ', p2
print,'**********************************************************'

;**********************************************************************
; renaming the save and eps files to represent data set used in the analysis 
;**********************************************************************
spawn,'mv emd.eps 27_rgs_emd.eps'
spawn,'mv flux.eps 27_rgs_wav.eps'
spawn,'mv results.sav 27_rgs_r.sav'

end
