subroutine ranlux(rvec,lenv)
!         Subtract-and-borrow random number generator proposed by
!         Marsaglia and Zaman, implemented by F. James with the name
!         RCARRY in 1991, and later improved by Martin Luescher
!         in 1993 to produce "Luxury Pseudorandom Numbers".
!     Fortran 77 coded by F. James, 1993
!
!       references:
!  M. Luscher, Computer Physics Communications  79 (1994) 100
!  F. James, Computer Physics Communications 79 (1994) 111
!
!   LUXURY LEVELS.
!   ------ ------      The available luxury levels are:
!
!  level 0  (p=24): equivalent to the original RCARRY of Marsaglia
!           and Zaman, very long period, but fails many tests.
!  level 1  (p=48): considerable improvement in quality over level 0,
!           now passes the gap test, but still fails spectral test.
!  level 2  (p=97): passes all known tests, but theoretically still
!           defective.
!  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible
!           correlations have very small chance of being observed.
!  level 4  (p=389): highest possible luxury, all 24 bits chaotic.
!
!!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!!!  Calling sequences for RANLUX:                                  ++
!!!!      CALL RANLUX (RVEC, LEN)   returns a vector RVEC of LEN     ++
!!!!                   32-bit random floating point numbers between  ++
!!!!                   zero (not included) and one (also not incl.). ++
!!!!      CALL RLUXGO(LUX,INT,K1,K2) initializes the generator from  ++
!!!!               one 32-bit integer INT and sets Luxury Level LUX  ++
!!!!               which is integer between zero and MAXLEV, or if   ++
!!!!               LUX > 24, it sets p=LUX directly.  K1 and K2   ++
!!!!               should be set to zero unless restarting at a break++
!!!!               point given by output of RLUXAT (see RLUXAT).     ++
!!!!      CALL RLUXAT(LUX,INT,K1,K2) gets the values of four integers++
!!!!               which can be used to restart the RANLUX generator ++
!!!!               at the current point by calling RLUXGO.  K1 and K2++
!!!!               specify how many numbers were generated since the ++
!!!!               initialization with LUX and INT.  The restarting  ++
!!!!               skips over  K1+K2*E9   numbers, so it can be long.++
!!!!   A more efficient but less convenient way of restarting is by: ++
!!!!      CALL RLUXIN(ISVEC)    restarts the generator from vector   ++
!!!!                   ISVEC of 25 32-bit integers (see RLUXUT)      ++
!!!!      CALL RLUXUT(ISVEC)    outputs the current values of the 25 ++
!!!!                 32-bit integer seeds, to be used for restarting ++
!!!!      ISVEC must be dimensioned 25 in the calling program        ++
!!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double precision rvec(lenv)
dimension seeds(24), iseeds(24), isdext(25)
parameter (maxlev=4, lxdflt=3)
dimension ndskip(0:maxlev)
dimension next(24)
parameter (twop12=4096., igiga=1000000000,jsdflt=314159265)
parameter (itwo24=2**24, icons=2147483563)
save notyet, i24, j24, carry, seeds, twom24, twom12, luxlev
save nskip, ndskip, in24, next, kount, mkount, inseed
integer luxlev
logical notyet
data notyet, luxlev, in24, kount, mkount /.true., lxdflt, 0,0,0/
data i24,j24,carry/24,10,0./
!                               default
!  Luxury Level   0     1     2   *3*    4
data ndskip/0,   24,   73,  199,  365 /
!orresponds to p=24    48    97   223   389
!     time factor 1     2     3     6    10   on slow workstation
!                 1    1.5    2     3     5   on fast mainframe
!
!  NOTYET is .TRUE. if no initialization has been performed yet.
!              Default Initialization by Multiplicative Congruential
if (notyet) then
   notyet = .false.
   jseed = jsdflt
   inseed = jseed
!         WRITE(6,'(A,I12)') ' RANLUX DEFAULT INITIALIZATION: ',JSEED
   luxlev = lxdflt
   nskip = ndskip(luxlev)
   lp = nskip + 24
   in24 = 0
   kount = 0
   mkount = 0
!        WRITE(6,'(A,I2,A,I4)')  ' RANLUX DEFAULT LUXURY LEVEL =  ',
!    +        LUXLEV,'      p =',LP
      twom24 = 1.
   do 25 i= 1, 24
      twom24 = twom24 * 0.5
   k = jseed/53668
   jseed = 40014*(jseed-k*53668) -k*12211
   if (jseed < 0)  jseed = jseed+icons
   iseeds(i) = mod(jseed,itwo24)
25    continue
   twom12 = twom24 * 4096.
   do 50 i= 1,24
   seeds(i) = real(iseeds(i))*twom24
   next(i) = i-1
50    continue
   next(1) = 24
   i24 = 24
   j24 = 10
   carry = 0.
   if (seeds(24) == 0.) carry = twom24
endif
!
!          The Generator proper: "Subtract-with-borrow",
!          as proposed by Marsaglia and Zaman,
!          Florida State University, March, 1989
!
do 100 ivec= 1, lenv
uni = seeds(j24) - seeds(i24) - carry
if (uni < 0.)  then
   uni = uni + 1.0
   carry = twom24
else
   carry = 0.
endif
seeds(i24) = uni
i24 = next(i24)
j24 = next(j24)
rvec(ivec) = uni
!  small numbers (with less than 12 "significant" bits) are "padded".
if (uni < twom12)  then
   rvec(ivec) = rvec(ivec) + twom24*seeds(j24)
!        and zero is forbidden in case someone takes a logarithm
   if (rvec(ivec) == 0.)  rvec(ivec) = twom24*twom24
endif
!        Skipping to luxury.  As proposed by Martin Luscher.
in24 = in24 + 1
if (in24 == 24)  then
   in24 = 0
   kount = kount + nskip
   do 90 isk= 1, nskip
   uni = seeds(j24) - seeds(i24) - carry
   if (uni < 0.)  then
      uni = uni + 1.0
      carry = twom24
   else
      carry = 0.
   endif
   seeds(i24) = uni
   i24 = next(i24)
   j24 = next(j24)
90    continue
endif
100 continue
kount = kount + lenv
if (kount >= igiga)  then
   mkount = mkount + 1
   kount = kount - igiga
endif
return
!
!           Entry to input and float integer seeds from previous run
entry rluxin(isdext)
!     The following IF block added by Phillip Helbig, based on conversation
!     with Fred James; an equivalent correction has been published by James.
if (notyet) then
   write(6,'(A)')  ' PROPER RESULTS ONLY WITH INITIALISATION FROM &
25 INTEGERS OBTAINED WITH RLUXUT'
   notyet = .false.
endif
   twom24 = 1.
   do 195 i= 1, 24
   next(i) = i-1
195    twom24 = twom24 * 0.5
   next(1) = 24
   twom12 = twom24 * 4096.
write(6,'(A)') ' FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:'
write(6,'(5X,5I12)') isdext
do 200 i= 1, 24
seeds(i) = real(isdext(i))*twom24
200 continue
carry = 0.
if (isdext(25) < 0)  carry = twom24
isd = iabs(isdext(25))
i24 = mod(isd,100)
isd = isd/100
j24 = mod(isd,100)
isd = isd/100
in24 = mod(isd,100)
isd = isd/100
luxlev = isd
  if (luxlev <= maxlev) then
    nskip = ndskip(luxlev)
!         WRITE (6,'(A,I2)') ' RANLUX LUXURY LEVEL SET BY RLUXIN TO: ',
!    +                         LUXLEV
  else  if (luxlev >= 24) then
    nskip = luxlev - 24
!         WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXIN TO:',LUXLEV
  else
    nskip = ndskip(maxlev)
!         WRITE (6,'(A,I5)') ' RANLUX ILLEGAL LUXURY RLUXIN: ',LUXLEV
    luxlev = maxlev
  endif
inseed = -1
return
!
!                    Entry to ouput seeds as integers
entry rluxut(isdext)
do 300 i= 1, 24
   isdext(i) = int(seeds(i)*twop12*twop12)
300 continue
isdext(25) = i24 + 100*j24 + 10000*in24 + 1000000*luxlev
if (carry > 0.)  isdext(25) = -isdext(25)
return
!
!                    Entry to output the "convenient" restart point
entry rluxat(lout,inout,k1,k2)
lout = luxlev
inout = inseed
k1 = kount
k2 = mkount
return
!
!                    Entry to initialize from one or three integers
entry rluxgo(lux,ins,k1,k2)
   if (lux < 0) then
      luxlev = lxdflt
   else if (lux <= maxlev) then
      luxlev = lux
   else if (lux < 24 .or. lux > 2000) then
      luxlev = maxlev
      write (6,'(A,I7)') ' RANLUX ILLEGAL LUXURY RLUXGO: ',lux
   else
      luxlev = lux
      do 310 ilx= 0, maxlev
        if (lux == ndskip(ilx)+24)  luxlev = ilx
310       continue
   endif
if (luxlev <= maxlev)  then
   nskip = ndskip(luxlev)
!        WRITE(6,'(A,I2,A,I4)') ' RANLUX LUXURY LEVEL SET BY RLUXGO :',
!    +        LUXLEV,'     P=', NSKIP+24
else
    nskip = luxlev - 24
!         WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXGO TO:',LUXLEV
endif
in24 = 0
if (ins < 0)  write (6,'(A)') &
   ' Illegal initialization by RLUXGO, negative input seed'
if (ins > 0)  then
  jseed = ins
!       WRITE(6,'(A,3I12)') ' RANLUX INITIALIZED BY RLUXGO FROM SEEDS',
!    +      JSEED, K1,K2
else
  jseed = jsdflt
!       WRITE(6,'(A)')' RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED'
endif
inseed = jseed
notyet = .false.
twom24 = 1.
   do 325 i= 1, 24
     twom24 = twom24 * 0.5
   k = jseed/53668
   jseed = 40014*(jseed-k*53668) -k*12211
   if (jseed < 0)  jseed = jseed+icons
   iseeds(i) = mod(jseed,itwo24)
325    continue
twom12 = twom24 * 4096.
   do 350 i= 1,24
   seeds(i) = real(iseeds(i))*twom24
   next(i) = i-1
350    continue
next(1) = 24
i24 = 24
j24 = 10
carry = 0.
if (seeds(24) == 0.) carry = twom24
!        If restarting at a break point, skip K1 + IGIGA*K2
!        Note that this is the number of numbers delivered to
!        the user PLUS the number skipped (if luxury > 0).
kount = k1
mkount = k2
if (k1+k2 /= 0)  then
  do 500 iouter= 1, k2+1
    inner = igiga
    if (iouter == k2+1)  inner = k1
    do 450 isk= 1, inner
      uni = seeds(j24) - seeds(i24) - carry
      if (uni < 0.)  then
         uni = uni + 1.0
         carry = twom24
      else
         carry = 0.
      endif
      seeds(i24) = uni
      i24 = next(i24)
      j24 = next(j24)
450     continue
500   continue
!         Get the right value of IN24 by direct calculation
  in24 = mod(kount, nskip+24)
  if (mkount > 0)  then
     izip = mod(igiga, nskip+24)
     izip2 = mkount*izip + in24
     in24 = mod(izip2, nskip+24)
  endif
!       Now IN24 had better be between zero and 23 inclusive
  if (in24 > 23) then
     write (6,'(A/A,3I11,A,I5)') &
    '  Error in RESTARTING with RLUXGO:','  The values', ins, &
     k1, k2, ' cannot occur at luxury level', luxlev
     in24 = 0
  endif
endif
return
end
