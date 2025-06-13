
      function elliptic_fm ( m )
    
      !*****************************************************************************80
      !
      !! ELLIPTIC_FM evaluates the complete elliptic integral F(M).
      !
      !  Discussion:
      !
      !    The value is computed using Carlson elliptic integrals:
      !
      !      F(m) = RF ( 0, 1-m, 1 ).
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    29 May 2018
      !
      !  Author:
      !
      !    John Burkardt
      !
      !  Parameters:
      !
      !    Input, real ( kind = rk ) M, the argument.
      !
      !    Output, real ( kind = rk ) ELLIPTIC_FM, the function value.
      !
      implicit none
    
      integer, parameter :: rk = kind ( 1.0D+00 )
    
      real ( kind = rk ) elliptic_fm
      real ( kind = rk ) errtol
      integer ierr
      real ( kind = rk ) m
      real ( kind = rk ) rf
      real ( kind = rk ) value
      real ( kind = rk ) x
      real ( kind = rk ) y
      real ( kind = rk ) z
    
      x = 0.0D+00
      y = 1.0D+00 - m
      z = 1.0D+00
      errtol = 1.0D-03
    
      value = rf ( x, y, z, errtol, ierr )
    
      elliptic_fm = value
    
      return
      end

      function elliptic_inc_fm ( phi, m )
    
      !*****************************************************************************80
      !
      !! ELLIPTIC_INC_FM evaluates the incomplete elliptic integral F(PHI,M).
      !
      !  Discussion:
      !
      !    The value is computed using Carlson elliptic integrals:
      !
      !      F(phi,m) = sin(phi) * RF ( cos^2 ( phi ), 1-m sin^2 ( phi ), 1 )
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    24 June 2018
      !
      !  Author:
      !
      !    John Burkardt
      !
      !  Parameters:
      !
      !    Input, real ( kind = rk ) PHI, M, the argument.
      !    0 <= PHI <= PI/2.
      !    0 <= M * sin^2(PHI) <= 1.
      !
      !    Output, real ( kind = rk ) ELLIPTIC_INC_FM, the function value.
      !
      implicit none
    
      integer, parameter :: rk = kind ( 1.0D+00 )
    
      real ( kind = rk ) cp
      real ( kind = rk ) elliptic_inc_fm
      real ( kind = rk ) errtol
      integer ierr
      real ( kind = rk ) m
      real ( kind = rk ) phi
      real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
      real ( kind = rk ) rf
      real ( kind = rk ) sp
      real ( kind = rk ) value
      real ( kind = rk ) x
      real ( kind = rk ) y
      real ( kind = rk ) z
    
      cp = cos ( phi )
      sp = sin ( phi )
      x = cp * cp
      y = 1.0D+00 - m * sp ** 2
      z = 1.0D+00
      errtol = 1.0D-03
    
      value = rf ( x, y, z, errtol, ierr )
    
      if ( ierr /= 0 ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'ELLIPTIC_INC_FM - Fatal error!'
        write ( *, '(a,i2)' ) '  RF returned IERR = ', ierr
        stop 1
      end if
    
      elliptic_inc_fm = sp * value
    
      return
      end


      function rf ( x, y, z, errtol, ierr )
    
    !*****************************************************************************80
    !
    !! RF computes an incomplete elliptic integral of the first kind, RF(X,Y,Z).
    !
    !  Discussion:
    !
    !    This function computes the incomplete elliptic integral of the first kind.
    !
    !    RF(X,Y,Z) = Integral ( 0 <= T < oo )
    !
    !                    -1/2     -1/2     -1/2
    !          (1/2)(T+X)    (T+Y)    (T+Z)    DT,
    !
    !    where X, Y, and Z are nonnegative and at most one of them is zero.
    !
    !    If X or Y or Z is zero, the integral is complete.
    !
    !    The duplication theorem is iterated until the variables are
    !    nearly equal, and the function is then expanded in Taylor
    !    series to fifth order.
    !
    !    Check by addition theorem:
    !
    !      RF(X,X+Z,X+W) + RF(Y,Y+Z,Y+W) = RF(0,Z,W),
    !      where X, Y, Z, W are positive and X * Y = Z * W.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    30 May 2018
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Bille Carlson, Elaine Notis.
    !    This FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Bille Carlson,
    !    Computing Elliptic Integrals by Duplication,
    !    Numerische Mathematik,
    !    Volume 33, 1979, pages 1-16.
    !
    !    Bille Carlson, Elaine Notis,
    !    Algorithm 577, Algorithms for Incomplete Elliptic Integrals,
    !    ACM Transactions on Mathematical Software,
    !    Volume 7, Number 3, pages 398-403, September 1981.
    !
    !  Parameters:
    !
    !    Input, real ( kind = rk ) X, Y, Z, the arguments in the integral.
    !
    !    Input, real ( kind = rk ) ERRTOL, the error tolerance.
    !    Relative error due to truncation is less than
    !      ERRTOL ^ 6 / (4 * (1 - ERRTOL)).
    !    Sample choices:
    !      ERRTOL   Relative truncation error less than
    !      1.D-3    3.D-19
    !      3.D-3    2.D-16
    !      1.D-2    3.D-13
    !      3.D-2    2.D-10
    !      1.D-1    3.D-7
    !
    !    Output, integer IERR, the error flag.
    !    0, no error occurred.
    !    1, abnormal termination.
    !
      implicit none
    
      integer, parameter :: rk = kind ( 1.0D+00 )
    
      real ( kind = rk ) c1
      real ( kind = rk ) c2
      real ( kind = rk ) c3
      real ( kind = rk ) e2
      real ( kind = rk ) e3
      real ( kind = rk ) epslon
      real ( kind = rk ) errtol
      integer ierr
      real ( kind = rk ) lamda
      real ( kind = rk ) lolim
      real ( kind = rk ) mu
      real ( kind = rk ) rf
      real ( kind = rk ) s
      real ( kind = rk ) uplim
      real ( kind = rk ) x
      real ( kind = rk ) xn
      real ( kind = rk ) xndev
      real ( kind = rk ) xnroot
      real ( kind = rk ) y
      real ( kind = rk ) yn
      real ( kind = rk ) yndev
      real ( kind = rk ) ynroot
      real ( kind = rk ) z
      real ( kind = rk ) zn
      real ( kind = rk ) zndev
      real ( kind = rk ) znroot
    !
    !  LOLIM AND UPLIM DETERMINE THE RANGE OF VALID ARGUMENTS.
    !  LOLIM IS NOT LESS THAN THE MACHINE MINIMUM MULTIPLIED BY 5.
    !  UPLIM IS NOT GREATER THAN THE MACHINE MAXIMUM DIVIDED BY 5.
    !
      save lolim
      save uplim
    
      data lolim /3.D-78/
      data uplim /1.D+75/
    
      if ( x < 0.0D+00 .or. y < 0.0D+00 .or. z < 0.0D+00 .or. x + y < lolim .or. x + z < lolim .or. y + z < lolim .or. uplim <= x .or. uplim <= y .or. uplim <= z ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'RF - Error!'
        write ( *, '(a)' ) '  Invalid input arguments.'
        write ( *, '(a,d23.16)' ) '  X = ', x
        write ( *, '(a,d23.16)' ) '  Y = ', y
        write ( *, '(a,d23.16)' ) '  Z = ', z
        write ( *, '(a)' ) ''
        ierr = 1
        rf = 0.0D+00
        return
      end if
    
      ierr = 0
      xn = x
      yn = y
      zn = z
    
      do
    
        mu = ( xn + yn + zn ) / 3.0d0
        xndev = 2.0d0 - ( mu + xn ) / mu
        yndev = 2.0d0 - ( mu + yn ) / mu
        zndev = 2.0d0 - ( mu + zn ) / mu
        epslon = max ( abs ( xndev ), abs ( yndev ), abs ( zndev ) )
    
        if ( epslon < errtol ) then
          c1 = 1.0d0 / 24.0d0
          c2 = 3.0d0 / 44.0d0
          c3 = 1.0d0 / 14.0d0
          e2 = xndev * yndev - zndev * zndev
          e3 = xndev * yndev * zndev
          s = 1.0d0 + ( c1 * e2 - 0.1d0 - c2 * e3 ) * e2 + c3 * e3
          rf = s / sqrt ( mu )
          return
        end if
    
        xnroot = sqrt ( xn )
        ynroot = sqrt ( yn )
        znroot = sqrt ( zn )
        lamda = xnroot * ( ynroot + znroot ) + ynroot * znroot
        xn = ( xn + lamda ) * 0.25d0
        yn = ( yn + lamda ) * 0.25d0
        zn = ( zn + lamda ) * 0.25d0
    
      end do
    
      end

      function EllipticPiDiff(r) result(output)
         real(kind=kind ( 1.0D+00 )), intent(in) :: r
         real(kind=kind ( 1.0D+00 )) :: output
         real(kind=kind ( 1.0D+00 )) :: sqrtr
         sqrtr = dSqrt(r)
         if (sqrtr < 0.0003) then
            print*, 'Warning generated kt2/sborn too low for good accuracy', r
         endif

         if (sqrtr <= 0.002) then
            output = (1.3206478192228157*(-0.000011496699697342602 + sqrtr*(-0.009143569959264996 + sqrtr*(-0.9776358452127787 + sqrtr*(-1.1944386536885805 + 1.*sqrtr))))*
     -    (1.1096190557369234 + sqrtr*(-366.45203594546075 + sqrtr*(519800.68722648243 + 
     -            sqrtr*(-4.035522349372441d8 + sqrtr*(1.8087340994427414d11 + sqrtr*(-4.422365781500885d13 + 4.59004629114671d15*sqrtr))))))*dATan(16.384611418022235*Sqrt(sqrtr))*
     -    (EXP(1.663043827638274*sqrtr) - 0.00003813365717666205*Cos(21.156592658523934*Sqrt(sqrtr)))*Log(2*sqrtr))/
     -  (EXP(1.663043827638274*sqrtr)*(-2.475163988222491d-6 + sqrtr*(-0.008358249399007839 + sqrtr*(-1.4555307600794787 + sqrtr*(-0.8003540934283068 + 1.*sqrtr)))))
         elseif (sqrtr <= 0.0045) then
            output = (1.3206478192228157*(-0.000011496699697342602 + sqrtr*(-0.009143569959264996 + sqrtr*(-0.9776358452127787 + sqrtr*(-1.1944386536885805 + 1.*sqrtr))))*
     -    (0.9984889856513807 + sqrtr*(-2.517948816848558 + sqrtr*(4906.532481831693 + 
     -            sqrtr*(-2.7555342655692487d6 + sqrtr*(7.354294337966862d8 + sqrtr*(-9.786970001383696d10 + 5.265654967993501d12*sqrtr))))))*dATan(16.384611418022235*Sqrt(sqrtr))*
     -    (EXP(1.663043827638274*sqrtr) - 0.00003813365717666205*Cos(21.156592658523934*Sqrt(sqrtr)))*Log(2*sqrtr))/
     -  (EXP(1.663043827638274*sqrtr)*(-2.475163988222491d-6 + sqrtr*(-0.008358249399007839 + sqrtr*(-1.4555307600794787 + sqrtr*(-0.8003540934283068 + 1.*sqrtr)))))
         elseif (sqrtr <= 0.008) then
            output = (-7.516144361969055d10*(-0.000011496699697342602 + sqrtr*(-0.009143569959264996 + sqrtr*(-0.9776358452127787 + sqrtr*(-1.1944386536885805 + 1.*sqrtr))))*
     -    (-1.751013982959634d-11 + sqrtr*(-7.382148542380893d-11 + sqrtr*(3.0630886296309246d-8 + 
     -            sqrtr*(-6.2149894088874274d-6 + sqrtr*(0.0006892943150668142 + sqrtr*(-0.04058356823292348 + 1.*sqrtr))))))*dATan(16.384611418022235*Sqrt(sqrtr))*
     -    (EXP(1.663043827638274*sqrtr) - 0.00003813365717666205*Cos(21.156592658523934*Sqrt(sqrtr)))*Log(2*sqrtr))/
     -  (EXP(1.663043827638274*sqrtr)*(-2.475163988222491d-6 + sqrtr*(-0.008358249399007839 + sqrtr*(-1.4555307600794787 + sqrtr*(-0.8003540934283068 + 1.*sqrtr)))))
         elseif (sqrtr <= 0.015) then
            output = (3.9539703374534434d8*(-0.000011496699697342602 + sqrtr*(-0.009143569959264996 + sqrtr*(-0.9776358452127787 + sqrtr*(-1.1944386536885805 + 1.*sqrtr))))*
     -    (3.352062466933086d-9 + sqrtr*(-5.19256736846263d-9 + sqrtr*(8.39256593190637d-7 + 
     -            sqrtr*(-0.00007069713049821346 + sqrtr*(0.0033955631379316543 + sqrtr*(-0.0891378118695529 + 1.*sqrtr))))))*dATan(16.384611418022235*Sqrt(sqrtr))*
     -    (EXP(1.663043827638274*sqrtr) - 0.00003813365717666205*Cos(21.156592658523934*Sqrt(sqrtr)))*Log(2*sqrtr))/
     -  (EXP(1.663043827638274*sqrtr)*(-2.475163988222491d-6 + sqrtr*(-0.008358249399007839 + sqrtr*(-1.4555307600794787 + sqrtr*(-0.8003540934283068 + 1.*sqrtr)))))
         elseif (sqrtr <= 0.03) then
            output = (7.932641779659658d6*(-0.000011496699697342602 + sqrtr*(-0.009143569959264996 + sqrtr*(-0.9776358452127787 + sqrtr*(-1.1944386536885805 + 1.*sqrtr))))*
     -    (1.663899296139464d-7 + sqrtr*(-8.03515583598002d-9 + sqrtr*(2.9397567590670356d-6 + 
     -            sqrtr*(-0.00022463011088694415 + sqrtr*(0.007948594759133629 + sqrtr*(-0.14000220116258497 + 1.*sqrtr))))))*dATan(16.384611418022235*Sqrt(sqrtr))*
     -    (EXP(1.663043827638274*sqrtr) - 0.00003813365717666205*Cos(21.156592658523934*Sqrt(sqrtr)))*Log(2*sqrtr))/
     -  (EXP(1.663043827638274*sqrtr)*(-2.475163988222491d-6 + sqrtr*(-0.008358249399007839 + sqrtr*(-1.4555307600794787 + sqrtr*(-0.8003540934283068 + 1.*sqrtr)))))
         elseif (sqrtr <= 0.1) then
            output = (-10494.732375926154*(-0.000011496699697342602 + sqrtr*(-0.009143569959264996 + sqrtr*(-0.9776358452127787 + sqrtr*(-1.1944386536885805 + 1.*sqrtr))))*
     -    (-0.00012585395451873276 + sqrtr*(-3.0785799530737665d-6 + sqrtr*(0.00022362182996384469 + 
     -            sqrtr*(-0.005543347193758585 + sqrtr*(0.06681152087762894 + sqrtr*(-0.405071352885517 + 1.*sqrtr))))))*dATan(16.384611418022235*Sqrt(sqrtr))*
     -    (EXP(1.663043827638274*sqrtr) - 0.00003813365717666205*Cos(21.156592658523934*Sqrt(sqrtr)))*Log(2*sqrtr))/
     -  (EXP(1.663043827638274*sqrtr)*(-2.475163988222491d-6 + sqrtr*(-0.008358249399007839 + sqrtr*(-1.4555307600794787 + sqrtr*(-0.8003540934283068 + 1.*sqrtr)))))
         elseif (sqrtr <= 0.2) then
            output = (77.50981175970828*(-0.000011496699697342602 + sqrtr*(-0.009143569959264996 + sqrtr*(-0.9776358452127787 + sqrtr*(-1.1944386536885805 + 1.*sqrtr))))*
     -    (0.01704421967489593 + sqrtr*(-0.00041549323618647545 + sqrtr*(0.00868830948785372 + sqrtr*(-0.08202033992424743 + sqrtr*(0.39897846700121553 + sqrtr*(-0.9890835747268558 + 1.*sqrtr))))))*
     -    dATan(16.384611418022235*Sqrt(sqrtr))*(EXP(1.663043827638274*sqrtr) - 0.00003813365717666205*Cos(21.156592658523934*Sqrt(sqrtr)))*Log(2*sqrtr))/
     -  (EXP(1.663043827638274*sqrtr)*(-2.475163988222491d-6 + sqrtr*(-0.008358249399007839 + sqrtr*(-1.4555307600794787 + sqrtr*(-0.8003540934283068 + 1.*sqrtr)))))
         else
            output = (-0.5697359385578563*(-0.000011496699697342602 + sqrtr*(-0.009143569959264996 + sqrtr*(-0.9776358452127787 + sqrtr*(-1.1944386536885805 + 1.*sqrtr))))*
     -    (-2.319675450777865 + sqrtr*(0.01482009761226051 + sqrtr*(-0.011804545319164346 + sqrtr*(-0.24601532828043546 + sqrtr*(1.0040463605662768 + sqrtr*(-1.6090052594727977 + 1.*sqrtr))))))*
     -    dATan(16.384611418022235*Sqrt(sqrtr))*(EXP(1.663043827638274*sqrtr) - 0.00003813365717666205*Cos(21.156592658523934*Sqrt(sqrtr)))*Log(2*sqrtr))/
     -  (EXP(1.663043827638274*sqrtr)*(-2.475163988222491d-6 + sqrtr*(-0.008358249399007839 + sqrtr*(-1.4555307600794787 + sqrtr*(-0.8003540934283068 + 1.*sqrtr)))))
         endif
         ! if (r<=11d0/128d0) then
         !    output = (0.17341683129703184d0 - 0.9825520079285872*(-0.0625d0 + r) + 4.116414898580441*(-0.0625 + r)**2 - 37.57897966023478*(-0.0625 + r)**3 + 354.61856844520617*(-0.0625 + r)**4)*Log(4d0*r)
         ! elseif (r<=21d0/128d0) then
         !    output = (0.12213591082749702d0 - 0.7229667725698756*(-0.125 + r) + 0.8149006661827478*(-0.125 + r)**2 - 9.41053563501984*(-0.125 + r)**3 + 17.835272551926323*(-0.125 + r)**4)*Log(4d0*r)
         ! else
         !    output = (0.07784562894186256d0 - 0.7358538011127749*(-0.1875 + r) - 1.2922720445194535*(-0.1875 + r)**2 - 19.609924430322337*(-0.1875 + r)**3 - 168.73839450285504*(-0.1875 + r)**4)*Log(4d0*r)
         ! endif
      end function EllipticPiDiff


      function EllipticPiDiff_new(r) result(output)
         real(kind=kind ( 1.0D+00 )), intent(in) :: r
         real(kind=kind ( 1.0D+00 )) :: output
         real(kind=kind ( 1.0D+00 )) :: sqrtr, epsilon, xmin, xmax
         real(kind=kind ( 1.0D+00 )), external :: dgauss

         if(sqrtr < 0.0003) then
            print*, 'Warning generated kt2/sborn too low for good accuracy', r
         endif
         
         xmin = 1d0/dsqrt(2d0)  ! Sin(Pi/4)
         xmax = 1d0
         sqrtr = dSqrt(r)
         epsilon = 1d-8*(0.25d0-r)**2   ! decrease as r->0.25

         output = dgauss(integrand, xmin, xmax, epsilon)
         output = output  *2*(2d0*sqrtr+1d0)/Sqrt(1d0-2d0*sqrtr)
         
         contains
         
         function integrand(t) result(out)
         real(kind=kind ( 1.0D+00 )), intent(in) :: t
         real(kind=kind ( 1.0D+00 )) :: out, n, k
         n = 1d0/(0.5d0-sqrtr)
         k = 2d0 - n
         out = 1d0/((1d0-n*t**2)*sqrt((1d0-t**2)*(1d0-k*t**2)))
         end
      
      end function EllipticPiDiff_new
