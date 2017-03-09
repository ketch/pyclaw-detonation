! =========================================================
subroutine rp1(maxmx,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =========================================================

! solve Riemann problems for the 1D Euler reactive Euler equations using Roe's
! approximate Riemann solver for the Euler system plus a passive
! tracer.

! waves: 3
! equations: 4

! Conserved quantities:
!       1 density
!       2 momentum
!       3 energy
!       4 lamda

! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell
! On output, wave contains the waves,
!            s the speeds,
!            amdq the  left-going flux difference  A^- \Delta q
!            apdq the right-going flux difference  A^+ \Delta q

! Note that the i'th Riemann problem has left state qr(:,i-1)
!                                    and right state ql(:,i)
! From the basic clawpack routine step1, rp is called with ql = qr = q.


    implicit double precision (a-h,o-z)
    dimension   ql(meqn,1-mbc:maxmx+mbc)
    dimension   qr(meqn,1-mbc:maxmx+mbc)
    dimension    s(mwaves,1-mbc:maxmx+mbc)
    dimension wave(meqn, mwaves,1-mbc:maxmx+mbc)
    dimension amdq(meqn,1-mbc:maxmx+mbc)
    dimension apdq(meqn,1-mbc:maxmx+mbc)

!     # local storage
!     ---------------
    dimension delta(meqn)
    dimension u(1-mbc:maxmx+mbc),enth(1-mbc:maxmx+mbc),Y(1-mbc:maxmx+mbc)
    dimension c(1-mbc:maxmx+mbc)
    logical :: efix
    common /cparam/  gamma,gamma1,qheat,xfspeed

    data efix /.true./    !# use entropy fix for transonic rarefactions

!     # Compute Roe-averaged quantities:

    do 20 i=2-mbc,mx+mbc
        rhsqrtl = dsqrt(qr(1,i-1))
        rhsqrtr = dsqrt(ql(1,i))
        pl = gamma1*(qr(3,i-1) &
        - 0.5d0*(qr(2,i-1)**2)/qr(1,i-1)-qr(4,i-1)*qheat)
        pr = gamma1*(ql(3,i) &
        - 0.5d0*(ql(2,i)**2)/ql(1,i)-qr(4,i)*qheat)
        rhsq2 = rhsqrtl + rhsqrtr
        u(i) = (qr(2,i-1)/rhsqrtl + ql(2,i)/rhsqrtr) / rhsq2
        enth(i) = (((qr(3,i-1)+pl)/rhsqrtl &
        + (ql(3,i)+pr)/rhsqrtr)) / rhsq2
        Y(i) = (qr(4,i-1)/rhsqrtl + ql(4,i)/rhsqrtr) / rhsq2
        c2 = gamma1*(enth(i) - .5d0*u(i)**2 - qheat*Y(i))
        c(i) = dsqrt(c2)
    20 END DO


    do 30 i=2-mbc,mx+mbc

    !        # find a1 thru a4, the coefficients of the 4 eigenvectors:

        delta(1) = ql(1,i) - qr(1,i-1)
        delta(2) = ql(2,i) - qr(2,i-1)
        delta(3) = ql(3,i) - qr(3,i-1)
        delta(4) = ql(4,i) - qr(4,i-1)
        a2 = gamma1/c(i)**2 * ((enth(i)-u(i)**2-qheat*Y(i))*delta(1) &
        + u(i)*delta(2) - delta(3) + qheat*delta(4))
        a4 = ( (c(i)-u(i))*delta(1) + delta(2) - c(i)*a2 ) / (2.d0*c(i))
        a1 = delta(1) - a2 - a4
        a3 = -Y(i)*delta(1) + delta(4) + Y(i)*a2

    !        # Compute the waves. These are simply the weights a1, a2, a3, a4 times eigvec

        wave(1,1,i) = a1
        wave(2,1,i) = a1*(u(i)-c(i))
        wave(3,1,i) = a1*(enth(i) - u(i)*c(i))
        wave(4,1,i) = a1*(Y(i))
        s(1,i) = u(i)-c(i)-xfspeed

        wave(1,2,i) = a2
        wave(2,2,i) = a2*u(i)
        wave(3,2,i) = a2*0.5d0*u(i)**2
        wave(4,2,i) = 0.d0
        s(2,i) = u(i)-xfspeed

        wave(1,3,i) = 0.d0
        wave(2,3,i) = 0.d0
        wave(3,3,i) = a3*qheat
        wave(4,3,i) = a3
        s(3,i) = u(i)-xfspeed

        wave(1,4,i) = a4
        wave(2,4,i) = a4*(u(i)+c(i))
        wave(3,4,i) = a4*(enth(i)+u(i)*c(i))
        wave(4,4,i) = a4*Y(i)
        s(4,i) = u(i)+c(i)-xfspeed
    30 END DO

!     # compute Godunov flux f0:
!     --------------------------


    if (efix) go to 110

!     # no entropy fix
!     ----------------

!     # amdq = SUM s*wave   over left-going waves
!     # apdq = SUM s*wave   over right-going waves

    do 100 m=1,4
        do 100 i=2-mbc, mx+mbc
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
            do 90 mw=1,mwaves
                if (s(mw,i) < 0.d0) then
                    amdq(m,i) = amdq(m,i) + s(mw,i)*wave(m,mw,i)
                else
                    apdq(m,i) = apdq(m,i) + s(mw,i)*wave(m,mw,i)
                endif
            90 END DO
    100 END DO
    go to 900

!-----------------------------------------------------

    110 continue

!     # With entropy fix
!     ------------------

!    # compute flux differences amdq and apdq.
!    # First compute amdq as sum of s*wave for left going waves.
!    # Incorporate entropy fix by adding a modified fraction of wave
!    # if s should change sign.

    do 200 i=2-mbc,mx+mbc

    !        # check 1-wave:
    !        ---------------

        rhoim1 = qr(1,i-1)
        pim1 = gamma1*(qr(3,i-1) - 0.5d0*qr(2,i-1)**2 / rhoim1 - &
        qr(4,i-1)*qheat)
        cim1 = dsqrt(gamma*pim1/rhoim1)
        s0 = qr(2,i-1)/rhoim1 - cim1 - xfspeed    !# u-c in left state (cell i-1)

    !        # check for fully supersonic case:
        if (s0 >= 0.d0 .AND. s(1,i) > 0.d0)  then
        !            # everything is right-going, so no left waves
            do 60 m=1,4
                amdq(m,i) = 0.d0
            60 END DO
            go to 200
        endif

        rho1 = qr(1,i-1) + wave(1,1,i)
        rhou1 = qr(2,i-1) + wave(2,1,i)
        rhoen1 = qr(3,i-1) + wave(3,1,i)
        rhoY1 = qr(4,i-1) + wave(4,1,i)

        p1 = gamma1*(rhoen1 - 0.5d0*rhou1**2/rho1 - rhoY1*qheat)
        c1 = dsqrt(gamma*p1/rho1)
        s1 = rhou1/rho1 - c1 - xfspeed !# u-c to right of 1-wave
        if (s0 < 0.d0 .AND. s1 > 0.d0) then
        !            # transonic rarefaction in the 1-wave
            sfract = s0 * (s1-s(1,i)) / (s1-s0)
        else if (s(1,i) < 0.d0) then
        !            # 1-wave is leftgoing
            sfract = s(1,i)
        else
        !            # 1-wave is rightgoing
            sfract = 0.d0   !# this shouldn't happen since s0 < 0
        endif
        do 120 m=1,4
            amdq(m,i) = sfract*wave(m,1,i)
        120 END DO

    !        # check 2-wave:wave corresponds u characteristic for entropy
    !        ---------------

        if (s(2,i) >= 0.d0) go to 200  !# 2-wave is rightgoing
        do 130 m=1,4
            amdq(m,i) = amdq(m,i) + s(2,i)*wave(m,2,i)
        130 END DO

    !        # check 3-wave: corresponds to u characteristic for Y
    !        ---------------

        if (s(3,i) >= 0.d0) go to 200  !# 3-wave is rightgoing
        do 140 m=1,4
            amdq(m,i) = amdq(m,i) + s(3,i)*wave(m,3,i)
        140 END DO

    !        # check 4-wave: corresponds to u+c
    !        ---------------

        rhoi = ql(1,i)
        pi = gamma1*(ql(3,i) - 0.5d0*ql(2,i)**2 / rhoi - ql(4,i)*qheat)
        ci = dsqrt(gamma*pi/rhoi)
        s3 = ql(2,i)/rhoi + ci - xfspeed     !# u+c in right state  (cell i)

        rho2 = ql(1,i) - wave(1,4,i)
        rhou2 = ql(2,i) - wave(2,4,i)
        rhoen2 = ql(3,i) - wave(3,4,i)
        rhoY2 = ql(4,i) - wave(4,4,i)

        p2 = gamma1*(en2 - 0.5d0*rhou2**2/rho2 - rhoY2*qheat)
        c2 = dsqrt(gamma*p2/rho2)
        s2 = rhou2/rho2 + c2 - xfspeed   !# u+c to left of 4-wave
        if (s2 < 0.d0 .AND. s3 > 0.d0) then
        !            # transonic rarefaction in the 4-wave
            sfract = s2 * (s3-s(4,i)) / (s3-s2)
        else if (s(4,i) < 0.d0) then
        !            # 3-wave is leftgoing
            sfract = s(4,i)
        else
        !            # 3-wave is rightgoing
            go to 200
        endif

        do 160 m=1,4
            amdq(m,i) = amdq(m,i) + sfract*wave(m,4,i)
        160 END DO
    200 END DO

!     # compute the rightgoing flux differences:
!     # df = SUM s*wave   is the total flux difference and apdq = df - amdq

    do 220 m=1,4
        do 220 i = 2-mbc, mx+mbc
            df = 0.d0
            do 210 mw=1,mwaves
                df = df + s(mw,i)*wave(m,mw,i)
            210 END DO
            apdq(m,i) = df - amdq(m,i)
    220 END DO

    900 continue
    return
    end subroutine rp1
