!
subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
!
! Roe-solver for the Euler equations with a tracer variable and separate shear
! and entropy waves.
!
! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell
!
! This data is along a slice in the x-direction if ixy=1
!                            or the y-direction if ixy=2.
! On output, wave contains the waves, s the speeds,
! and amdq, apdq the decomposition of the flux difference
!     f(qr(i-1)) - f(ql(i))
! into leftgoing and rightgoing parts respectively.  With the Roe solver we have
!      amdq  =  A^- \Delta q    and    apdq  =  A^+ \Delta q
! where A is the Roe matrix.  An entropy fix can also be incorporated into the
! flux differences.

! Note that the i'th Riemann problem has left state qr(:,i-1)
!                                   and right state ql(:,i)
! From the basic clawpack routines, this routine is called with ql = qr
!
! This routine has been made thread safe by removing the common block storage
! of the Roe-averages.
!


    implicit none

    ! Input
    integer, intent(in) :: ixy, maxm, meqn, mwaves, mbc, mx, maux
    real(kind=8), intent(in) :: ql(meqn, 1-mbc:maxm+mbc)
    real(kind=8), intent(in) :: qr(meqn, 1-mbc:maxm+mbc)
    real(kind=8), intent(in) :: auxl(maux, 1-mbc:maxm+mbc)
    real(kind=8), intent(in) :: auxr(maux, 1-mbc:maxm+mbc)

    ! Output
    real(kind=8), intent(in out) :: wave(meqn, mwaves, 1-mbc:maxm+mbc)
    real(kind=8), intent(in out) :: s(mwaves, 1-mbc:maxm+mbc)

    real(kind=8), intent(in out) :: apdq(meqn, 1-mbc:maxm+mbc)
    real(kind=8), intent(in out) :: amdq(meqn, 1-mbc:maxm+mbc)

    ! Local storage
    integer :: i, m, mw, mu, mv
    real(kind=8) :: rho_sqrtl, rho_sqrtr, rho_sq2, u, v, u2v2, pl, pr, enth, c, Y
    real(kind=8) :: rhol, rhor, ul, ur, vl, vr, cl, cr
    real(kind=8) :: csquare, euvY
    real(kind=8) :: delta(5), a(5)

    real(kind=8) :: sl(mwaves, 1-mbc:maxm+mbc)
    real(kind=8) :: sr(mwaves, 1-mbc:maxm+mbc)

    real(kind=8) :: rhoim1, pim1, cim1, s0, rho1, rhou1, rhov1, rhoY1, en1, p1, c1, s1
    real(kind=8) :: sfract, rhoi, pi, ci, s3, rho2, rhou2, rhov2, rhoY2, en2, p2, c2
    real(kind=8) :: s2, df

    ! Use entropy fix for transonic rarefactions
    logical, parameter :: efix = .true.

    ! Common block storage
    ! Ideal gas constant
    real(kind=8) :: gamma, gamma1, qheat, xfspeed, fspeed
    common /cparam/  gamma,gamma1,qheat, xfspeed

    ! Set mu to point to  the component of the system that corresponds to
    ! momentum in the direction of this slice, mv to the orthogonal momentum:
    if (ixy == 1) then
        mu = 2
        mv = 3
        fspeed = xfspeed
    else
        mu = 3
        mv = 2
        fspeed = 0.d0
    endif
    ! Note that notation for u and v reflects assumption that the Riemann
    ! problems  are in the x-direction with u in the normal direciton and v in
    ! the orthogonal direcion, but with the above definitions of mu and mv the
    ! routine also works  with ixy=2 and returns, for example, f0 as the Godunov
    ! flux g0 for the Riemann problems u_t + g(u)_y = 0 in the y-direction.

    ! Initialize waves and speeds
    wave = 0.d0
    s = 0.d0

    ! Primary loop over grid cell interfaces
    do i = 2-mbc, mx+mbc

       ! Compute left and right states
       rhol =qr(1,i-1)
       ul = qr(mu,i-1)/qr(1,i-1)
       vl = qr(mv,i-1)/qr(1,i-1)
       pl = gamma1 * (qr(4,i-1) - 0.5d0 * (qr(2,i-1)**2 + qr(3,i-1)**2) / qr(1,i-1) &
            - qheat*qr(5,i-1))
       cl = sqrt(gamma*pl/rhol)

       rhor =ql(1,i)
       ur = ql(mu,i)/ql(1,i)
       vr = ql(mv,i)/ql(1,i)
       pr = gamma1 * (ql(4,i) - 0.5d0 * (ql(2,i)**2 + ql(3,i)**2) / ql(1,i) &
            - qheat*ql(5,i))
       cr = sqrt(gamma*pr/rhor)

        ! Compute Roe-averaged quantities
        rho_sqrtl = sqrt(qr(1,i-1))
        rho_sqrtr = sqrt(ql(1,i))
        rho_sq2 = rho_sqrtl + rho_sqrtr

        u = (qr(mu,i-1) / rho_sqrtl + ql(mu,i) / rho_sqrtr) / rho_sq2
        v = (qr(mv,i-1) / rho_sqrtl + ql(mv,i) / rho_sqrtr) / rho_sq2
        u2v2 = u**2 + v**2

        pl = gamma1 * (qr(4,i-1) - 0.5d0 * (qr(2,i-1)**2 + qr(3,i-1)**2) / qr(1,i-1) &
             - qheat*qr(5,i-1))
        pr = gamma1 * (ql(4,i)   - 0.5d0 * (ql(2,i)**2   + ql(3,i)**2) / ql(1,i) &
             - qheat*ql(5,i))
        enth = (((qr(4,i-1) + pl) / rho_sqrtl + (ql(4,i) + pr) / rho_sqrtr)) / rho_sq2

        Y = (qr(5,i-1)/rho_sqrtl + ql(5,i)/rho_sqrtr) / rho_sq2

        csquare = gamma1 * (enth - 0.5d0 * u2v2 - qheat*Y)
        c = sqrt(csquare)
        euvY = enth - u2v2 - qheat*Y

        ! Now split the jump in q at each interface into waves and find a1 thru
        ! a5, the coefficients of the 5 eigenvectors:
        delta(1) = ql(1,i) - qr(1,i-1)
        delta(2) = ql(mu,i) - qr(mu,i-1)
        delta(3) = ql(mv,i) - qr(mv,i-1)
        delta(4) = ql(4,i) - qr(4,i-1)
        delta(5) = ql(5,i) - qr(5,i-1)

        a(2) = gamma1*(delta(1) * euvY + delta(2)*u + delta(3)*v - delta(4) + qheat*delta(5)) &
               / (csquare)
        a(4) = (a(2) - delta(1)) * Y + delta(5)
        a(1) = (-c*a(2) + c*delta(1) + u*delta(1) - delta(2)) / (2*c)
        a(3) = -v*delta(1) + delta(3)
        a(5) = -a(1) - a(2) + delta(1)


        ! Compute the waves
        ! Acoustic, left-going
        wave( 1,1,i) = a(1)
        wave(mu,1,i) = a(1) * (u - c)
        wave(mv,1,i) = a(1) * v
        wave( 4,1,i) = a(1) * (enth - u * c)
        wave( 5,1,i) = a(1) * Y

        s(1,i)  = u - c - fspeed
        sl(1,i) = ul - cl - fspeed
        sr(1,i) = ur - cr - fspeed

        ! Entropy
        wave( 1,2,i) = a(2)
        wave(mu,2,i) = a(2) * u
        wave(mv,2,i) = a(2) * v
        wave( 4,2,i) = a(2) * 0.5d0 * u2v2
        wave( 5,2,i) = 0.d0

        s(2,i) = u - fspeed
        sl(2,i) = ul - fspeed
        sr(2,i) = ur - fspeed


        ! Shear wave
        wave( 1,3,i) = 0.d0
        wave(mu,3,i) = 0.d0
        wave(mv,3,i) = a(3)
        wave( 4,3,i) = a(3) * v
        wave( 5,3,i) = 0.d0

        s(3,i) = u - fspeed
        sl(3,i) = ul - fspeed
        sr(3,i) = ur - fspeed


        ! Reactant wave
        wave( 1,4,i) = 0.d0
        wave(mu,4,i) = 0.d0
        wave(mv,4,i) = 0.d0
        wave( 4,4,i) = a(4) * qheat
        wave( 5,4,i) = a(4)

        s(4,i) = u - fspeed
        sl(4,i) = ul - fspeed
        sr(4,i) = ur - fspeed

        ! Acoustic, right-going
        wave( 1,5,i) = a(5)
        wave(mu,5,i) = a(5) * (u + c)
        wave(mv,5,i) = a(5) * v
        wave( 4,5,i) = a(5) * (enth + u * c)
        wave( 5,5,i) = a(5) * Y

        s(5,i) = u + c - fspeed
        sl(5,i) = ul + cl  - fspeed
        sr(5,i) = ur + cr - fspeed



    end do

    if (efix) then
        ! Compute flux differences amdq and apdq.
        ! Adds numerical viscosity
        !
        ! amdq = SUM s*wave   over left-going waves
        ! apdq = SUM s*wave   over right-going waves
        amdq = 0.d0
        apdq = 0.d0
        do i=2-mbc, mx+mbc
           do mw=1,mwaves
              a = max(abs(sl(mw,i)),abs(sr(mw,i)))
              amdq(:,i) = amdq(:,i) + 0.5d0*(s(mw,i) - a) * wave(:,mw,i)
              apdq(:,i) = apdq(:,i) + 0.5d0*(s(mw,i) + a) * wave(:,mw,i)
           enddo
        enddo
     else
        ! Compute flux differences amdq and apdq.
        ! No entropy fix
        !
        ! amdq = SUM s*wave   over left-going waves
        ! apdq = SUM s*wave   over right-going waves
        amdq = 0.d0
        apdq = 0.d0
        do i=2-mbc, mx+mbc
            do mw=1,mwaves
                if (s(mw,i) < 0.d0) then
                    amdq(:,i) = amdq(:,i) + s(mw,i) * wave(:,mw,i)
                else
                    apdq(:,i) = apdq(:,i) + s(mw,i) * wave(:,mw,i)
                endif
            enddo
        enddo

        ! End of entropy corrections
     endif

end subroutine rpn2
