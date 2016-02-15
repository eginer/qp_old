 
      double precision function SABpartial_direction(l,zA,zB,A,B,nA,nB,gamA,gamB) 
      implicit double precision(a-h,o-z)
      dimension nA(3),nB(3)
      dimension A(3),B(3)
      gamtot=gamA+gamB
      SABpartial_direction=1.d0

      u=gamA/gamtot*A(l)+gamB/gamtot*B(l)
      arg=gamtot*u**2-gamA*A(l)**2-gamB*B(l)**2
      alpha=dexp(arg)
     &/gamtot**((1.d0+dfloat(nA(l))+dfloat(nB(l)))/2.d0)
      wA=dsqrt(gamtot)*(u-A(l))
      wB=dsqrt(gamtot)*(u-B(l))
      boundA=dsqrt(gamtot)*(zA-u)
      boundB=dsqrt(gamtot)*(zB-u)

      accu=0.d0
      do n=0,nA(l)
       do m=0,nB(l)
        integ=nA(l)+nB(l)-n-m
        accu=accu
     & +wA**n*wB**m*binom(nA(l),n)*binom(nB(l),m)
     & *(rinteg(integ,boundB)-rinteg(integ,boundA))
       enddo
      enddo
      SABpartial_direction=SABpartial_direction*accu*alpha
      end

