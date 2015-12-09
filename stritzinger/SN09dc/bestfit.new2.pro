FORWARD_FUNCTION mpfitfun_eval, mpfitfun, mpfit

; This is the call-back function for MPFIT.  It evaluates the
; function, subtracts the data, and returns the residuals.
function mpfitfun_eval, p, dp, _EXTRA=extra

  common mpfitfun_common, fcn, x, y, err, wts, f, fcnargs

  ;; The function is evaluated here.  There are four choices,
  ;; depending on whether (a) FUNCTARGS was passed to MPFITFUN, which
  ;; is passed to this function as "hf"; or (b) the derivative
  ;; parameter "dp" is passed, meaning that derivatives should be
  ;; calculated analytically by the function itself.
  if n_elements(fcnargs) GT 0 then begin
      if n_params() GT 1 then f = call_function(fcn, x, p, dp, _EXTRA=fcnargs)$
      else                    f = call_function(fcn, x, p, _EXTRA=fcnargs)
  endif else begin
      if n_params() GT 1 then f = call_function(fcn, x, p, dp) $
      else                    f = call_function(fcn, x, p)
  endelse

  ;; Compute the deviates, applying either errors or weights
  if n_elements(err) GT 0 then begin
      result = (y-f)/err
  endif else if n_elements(wts) GT 0 then begin
      result = (y-f)*wts
  endif else begin
      result = (y-f)
  endelse
      
  ;; Make sure the returned result is one-dimensional.
  result = reform(result, n_elements(result), /overwrite)
  return, result
  
end

;; Implement residual and gradient scaling according to the
;; prescription of Cash (ApJ, 228, 939)
pro mpfitfun_cash, resid, dresid
  common mpfitfun_common, fcn, x, y, err, wts, f, fcnargs

  sz = size(dresid)
  m = sz(1)
  n = sz(2)

  ;; Do rudimentary dimensions checks, so we don't do something stupid
  if n_elements(y) NE m OR n_elements(f) NE m OR n_elements(resid) NE m then begin
      DIM_ERROR:
      message, 'ERROR: dimensions of Y, F, RESID or DRESID are not consistent'
  endif

  ;; Scale gradient by sqrt(y)/f
  gfact = temporary(dresid) * rebin(reform(sqrt(y)/f,m,1),m,n)
  dresid = reform(dresid, m, n, /overwrite)
  
  ;; Scale residuals by 1/sqrt(y)
  resid = temporary(resid)/sqrt(y)

  return
end

function mpfitfun, fcn, x, y, err, p, WEIGHTS=wts, FUNCTARGS=fa, $
                   BESTNORM=bestnorm, nfev=nfev, STATUS=status, $
                   parinfo=parinfo, query=query, CASH=cash, $
                   covar=covar, perror=perror, yfit=yfit, $
                   niter=niter, nfree=nfree, npegged=npegged, dof=dof, $
                   quiet=quiet, ERRMSG=errmsg, _EXTRA=extra

  status = 0L
  errmsg = ''

  ;; Detect MPFIT and crash if it was not found
  catch, catcherror
  if catcherror NE 0 then begin
      MPFIT_NOTFOUND:
      catch, /cancel
      message, 'ERROR: the required function MPFIT must be in your IDL path', /info
      return, !values.d_nan
  endif
  if mpfit(/query) NE 1 then goto, MPFIT_NOTFOUND
  catch, /cancel
  if keyword_set(query) then return, 1

  if n_params() EQ 0 then begin
      message, "USAGE: PARMS = MPFITFUN('MYFUNCT', X, Y, ERR, "+ $
        "START_PARAMS, ... )", /info
      return, !values.d_nan
  endif
  if n_elements(x) EQ 0 OR n_elements(y) EQ 0 then begin
      message, 'ERROR: X and Y must be defined', /info
      return, !values.d_nan
  endif

  if n_elements(err) GT 0 OR n_elements(wts) GT 0 AND keyword_set(cash) then begin
      message, 'ERROR: WEIGHTS or ERROR cannot be specified with CASH', /info
      return, !values.d_nan
  endif
  if keyword_set(cash) then begin
      scalfcn = 'mpfitfun_cash'
  endif

  ;; Use common block to pass data back and forth
  common mpfitfun_common, fc, xc, yc, ec, wc, mc, ac
  fc = fcn & xc = x & yc = y & mc = 0L
  ;; These optional parameters must be undefined first
  ac = 0 & dummy = size(temporary(ac))
  ec = 0 & dummy = size(temporary(ec))
  wc = 0 & dummy = size(temporary(wc))

  if n_elements(fa) GT 0 then ac = fa
  if n_elements(wts) GT 0 then begin
      wc = sqrt(abs(wts))
  endif else if n_elements(err) GT 0 then begin
      wh = where(err EQ 0, ct)
      if ct GT 0 then begin
          message, 'ERROR: ERROR value must not be zero.  Use WEIGHTS.', $
            /info
          return, !values.d_nan
      endif
      ec = err
  endif

  result = mpfit('mpfitfun_eval', p, SCALE_FCN=scalfcn, $
                 parinfo=parinfo, STATUS=status, nfev=nfev, BESTNORM=bestnorm,$
                 covar=covar, perror=perror, $
                 niter=niter, nfree=nfree, npegged=npegged, dof=dof, $
                 ERRMSG=errmsg, quiet=quiet, _EXTRA=extra)

  ;; Retrieve the fit value
  yfit = temporary(mc)
  ;; Some cleanup
  xc = 0 & yc = 0 & wc = 0 & ec = 0 & mc = 0 & ac = 0

  ;; Print error message if there is one.
  if NOT keyword_set(quiet) AND errmsg NE '' then $
    message, errmsg, /info

  return, result
end

FUNCTION MYFUNC, X, P 

common constants, Nimass

days2sec = 86400.0
MeV2erg = 1.6d-6 ;ergs/MeV

halfni = 6.077*days2sec ;s
lamni = alog(2)/halfni
eni = 1.75 * MeV2erg ;erg 

halfco = 77.12*days2sec ;s
lamco = alog(2)/halfco
eco = 3.76 * MeV2erg ;erg tot gamma and positron more or less

ECO_gam = 3.61 * MeV2erg
eco_pos = 0.12 * MeV2erg

;number of N(56) atoms
Nni0 = 2.151463d55 * NImass ;Nimass input by user
n = n_elements(X)

A = dblarr(n)
B = dblarr(n)
C = dblarr(n)
D = dblarr(n)
Edep = dblarr(n)

for i=0,n-1 do begin
    A[i] = lamni*Nni0*exp(-lamni*X[i])*eni
    B[i] = lamco*Nni0*(lamni / (lamni - lamco) )
    C[i] = (exp(-lamco*X[i]) - exp(-lamni*X[i]))
    D[i] = (eco_pos + eco_gam * (1.0 - exp(-((P*days2sec)/X[i])^2.d0)))
    Edep[i] =  A[i] + B[i]*C[i]*D[i]
endfor
;print, X
;print, "error=",alog10(Edep),P
;stop
RETURN, Edep
;RETURN,(Y-Edep) 
END

pro bestfit

common constants, Nimass

Nimass = 2.21
;Nimass = 0.84468407

days2sec = 86400.0
MeV2erg = 1.6d-6 ;ergs/Mev
riset = 28.0

;start = [50.0*days2sec, 100.0*days2sec, 10.0, 70.0]
start = [10.0]
addnum=10

;di everything in erg/s
readcol,'2009dc_bolo_uvoir.txt',a1,a2,a3,format='d,d,d'
;print,X,Y

a2 = 10.^a2
a3 = 0.1*a2

X1 = a1*days2sec + riset*days2sec ;day 0 = explosion day

a4 = a1[0] - addnum + dindgen(addnum)
x1_new=dblarr(n_elements(a1)+n_elements(a4))
nx1=n_elements(x1_new)

for i =0,addnum-1 do begin
x1_new[i] = a4[i]
endfor

for i = addnum, nx1-1 do begin
x1_new[i] = (a1[i-addnum])*days2sec + riset*days2sec ;day 0 = explosion day
endfor



if (a2[0] lt 100) then begin
    Y1 = 10.d0^(a2)
    Yerr1 = 10.d0^(a3)
endif else begin

Y1 = a2
Yerr1 = a3
endelse

n=n_elements(X1)
i=0

;plot, X1/days2sec, alog10(Y1)
;stop
;Starting value of fitting range
while (x1[i] lt ((50.+riset)*days2sec)) do begin
i=i+1
if (i gt n) then begin
    print,"ERROR in X: starting value in fitting range out of range"
    stop
endif
endwhile
start_fit = i

;Last value in fitting range
i=0
while (x1[i] lt ((100.+riset)*days2sec)) do begin
i=i+1
if (i gt n) then begin
    print,"ERROR in X: last value in fitting range too high"
    stop
endif
endwhile
stop_fit = i


X = dblarr(stop_fit - start_fit + 1)
Y = dblarr(stop_fit - start_fit + 1)
Yerr = dblarr(stop_fit - start_fit + 1)

for j=0,(stop_fit-start_fit) do begin
X[j] = X1[start_fit + j]
Y[j] = Y1[start_fit + j]
Yerr[j]=Yerr1[start_fit + j]
endfor

n = n_elements(X)
;stop
perror=1
result = 0.d0
result = mpfitfun('MYFUNC',X,Y,Yerr,start,perror=perror) ; result is in days now!
; stop
print, "to =", result, " in days" 
print,perror

;for k = 0, N-1 do begin
;days[k]=t[k]/days2sec
;endfor

!P.FONT=-1
!x.thick=6
!y.thick=6
!p.thick=6

MeV2erg = 1.6d-6 ;ergs/Mev

halfni = 6.077*days2sec ;s
lamni = alog(2)/halfni
eni = 1.75 * MeV2erg ;erg 

halfco = 77.12*days2sec ;s
lamco = alog(2)/halfco
eco = 3.76 * MeV2erg ;erg tot gamma and positron more or less

ECO_gam = 3.61 * MeV2erg
eco_pos = 0.12 * MeV2erg

;number of N(56) atoms
Nni0 = 2.151463d55 * NImass ;Nimass input by user

N = 3
m = nx1
A = dblarr(m)
B = dblarr(m)
C = dblarr(m)
D = dblarr(m,N)
Edep = dblarr(m,N)
to=[0.0001,result,9999.999]
for i=0,m-1 do begin
	for j = 0, N - 1 do begin
    A[i] = lamni*Nni0*exp(-lamni*x1_new[i])*eni
    B[i] = lamco*Nni0*(lamni / (lamni - lamco) )
    C[i] = (exp(-lamco*x1_new[i]) - exp(-lamni*x1_new[i]))
    D[i,j] = (eco_pos + eco_gam * (1.0 - exp(- ((to[j]*days2sec)/x1_new[i])^2.0)))
    Edep[i,j] =  A[i] + B[i]*C[i]*D[i,j]
endfor
	endfor

plot, (x1_new/days2sec)-riset,alog10(Edep[*,0]),xstyle=1,yrange = [40.d0,43.5],title='SN2009dc 1.78M_sun of Ni, rise_t=22 days to=50.5',xtitle='Days since max',ytitle='log Edep (ergs/s)' ,xrange=[-20.0,115],charthick=6

for j = 1, N-1 do begin
oplot, (x1_new/days2sec)-riset, alog10(Edep[*,j])
endfor
tek_color
oplot,(X1/days2sec)-riset,alog10(Y1),color=4,thick=1,psym=4

;n_data = n_elements(a3)

set_plot,'ps'
device,file='SN09dc_2.Edep.eps',/color,/encapsulated

plot, (x1_new/days2sec)-riset,alog10(Edep[*,0]),xstyle=1,yrange = [40.d0,43.5],title='SN2009dc 2.21 M_sun of Ni, rise_t=28 days, to=48.2 days',xtitle='Days since max',ytitle='log Edep (ergs/s)' ,xrange=[-20.0,165],thick=6,charthick=4


for j = 1, N-1 do begin
oplot, (x1_new/days2sec)-riset, alog10(Edep[*,j]),thick=6
endfor

oplot,(X1/days2sec)-riset,alog10(Y1),color=4,thick=6,psym=6


pp1=[50,50]
pp2=[44,40]
oplot,pp1,pp2,linestyle=2


;plot, (X1/days2sec)-riset,alog10(Edep[*,0]),psym=5,xstyle=1,yrange = [40.d0,43.5],title='SN1991T 0.40M_sun of Ni, to=26.584',xtitle='days since bol max',ytitle='log Edep (ergs/s)' ,xrange=[-20.0,115]

;for j = 1, N-1 do begin
;oplot, (X1/days2sec)-riset, alog10(Edep[*,j])
;endfor
;oplot,(X1/days2sec)-riset,alog10(Y1), color=4,thick=3

device, /close
set_plot,'X'

;spawn, 'gv out1.eps &'

 end


