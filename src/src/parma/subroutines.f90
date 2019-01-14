module parma
    use iso_c_binding

    interface
        subroutine setlocale_to_avoid_bug() bind(C, name = "setlocale_to_avoid_bug")
            use iso_c_binding
        end subroutine setlocale_to_avoid_bug
    end interface

contains

    ! *******************************************************
    function getSpec(ip, s, r, d, e, g)
        ! ip: particle ID
        ! s: W index
        ! r: cut-off rigidity in GV
        ! d: atmospheric depth in g/cm2
        ! e: energy in MeV/n
        ! g: local geometry effect
        ! *******************************************************
        implicit real(kind = 8) (a-h, o-z)

        getSpec = 0.0
        if(ip==0) then  ! neutron
            getSpec = getNeutSpec(s, r, d, g, e)
        elseif(ip>=1.and.ip<=28) then ! proton to Ni
            getSpec = getIonSpec(ip, s, r, d, e)
        elseif(ip==29.or.ip==30) then ! Muon
            iptmp = ip - 28
            getSpec = getMuonSpec(iptmp, s, r, d, e)
        else
            iptmp = ip - 28 ! 3:electron, 4:positron, 5:photon
            getSpec = getsecondary(iptmp, s, r, d, e)
        endif

        return

    end

    ! *******************************************************
    function getFl(ip, s, r, d)  ! get Fl value, s:solar modulation potential, r:Cut off rigidity, d:depth
        ! *******************************************************
        implicit real(kind = 8) (a-h, o-z)
        parameter (npart = 11) ! number of particle type, 0:neutron, 1:proton, 2:alpha, 3:electron, 4:positron, 5:photon, 6-11: Li-O
        parameter (nBdata = 4) ! B(1) - B(4) : Fl= B(1)*(exp(-B(2)*d)-B(3)*exp(-B(4)*d))
        parameter (nAdata = 10) ! A(1) - A(10) : Bmin = A(1)+A(2)*r+A(3)/(1+exp((r-A(4))/A(5))), Bmin = A(6)+A(7)*r+A(8)/(1+exp((r-A(9))/A(10)))
        parameter (nsor = 2) ! solar minimum & maximum
        character chatmp1*1, chatmp5*5
        character, save :: pname(0:npart)*6
        real(kind = 8), save :: A(0:npart, nBdata, nAdata), B(nBdata)
        real(kind = 8), save :: spot(nsor), FL(nsor)
        integer(kind = 4), save :: ifirst
        data ifirst/0/
        data spot/0.0, 150.0/
        data pname/'neutro', 'proton', 'alphaa', 'elemag', 'elemag', 'elemag', &
                & 'ions  ', 'ions  ', 'ions  ', 'ions  ', 'ions  ', 'ions  '/

        if(ifirst==0) then
            ifirst = 1
            do i = 0, 5
                if(i<=2) then ! neutron, proton, alpha
                    open(unit = 15, file = 'input/' // pname(i) // '/Rigid-Dep.inp', status = 'old')
                elseif(i==3) then
                    open(unit = 15, file = 'input/' // pname(i) // '/Rigid-Dep-EL.inp', status = 'old')
                elseif(i==4) then
                    open(unit = 15, file = 'input/' // pname(i) // '/Rigid-Dep-PO.inp', status = 'old')
                elseif(i==5) then
                    open(unit = 15, file = 'input/' // pname(i) // '/Rigid-Dep-PH.inp', status = 'old')
                endif
                read(15, '(a)') chatmp1
                do ib = 1, nBdata
                    read(15, 1010) chatmp5, (A(i, ib, ia), ia = 1, nAdata)
                enddo
                close(unit = 15)
            enddo
            open(unit = 15, file = 'input/ions/Rigid-Dep.inp', status = 'old')
            do i = npart - 5, npart
                read(15, '(a)') chatmp1
                do ib = 1, nBdata
                    read(15, 1010) chatmp5, (A(i, ib, ia), ia = 1, nAdata)
                enddo
            enddo
            close(unit = 15)
        endif
        1010 format(a5, 30es13.5)

        do is = 1, nsor ! solar minimum and maximum
            if(ip==1.or.ip==2.or.ip>=npart - 5) then ! need not correction
                r1 = r
            else ! for neutron, electron, positron and photon, Rc for high-altitude should be corrected
                if(ip==0) then
                    ipidx = ip ! for neutron
                else
                    ipidx = ip - 2 ! 1:electron, 2:positron, 3:photon
                endif
                r1 = r * getBestR(is, r, d, ipidx)
            endif
            do ib = 1, nBdata
                if(is==1) then ! solar minimum
                    B(ib) = A(ip, ib, 1) + A(ip, ib, 2) * r1 + A(ip, ib, 3) / (1 + exp((r1 - A(ip, ib, 4)) / A(ip, ib, 5)))
                else
                    B(ib) = A(ip, ib, 6) + A(ip, ib, 7) * r1 + A(ip, ib, 8) / (1 + exp((r1 - A(ip, ib, 9)) / A(ip, ib, 10)))
                endif
            enddo
            Fl(is) = B(1) * (exp(-B(2) * d) - B(3) * exp(-B(4) * d))
        enddo

        if(ip<=5) then
            pow = getPow(ip, d, r) ! get Power index
        else
            pow = getPow(ip + 2, d, r) ! in getPow, ip should be +2 for ions
        endif
        getFl = 0.0
        A2 = (Fl(1) - Fl(2)) / (getFFPfromW(spot(1))**pow - getFFPfromW(spot(2))**pow)
        A1 = Fl(1) - A2 * getFFPfromW(spot(1))**pow
        getFl = a1 + a2 * getFFPfromW(s)**pow

        return
    end

    ! **********************************************************
    function getFFPfromW(s) ! get FFP (MV) from W (sun spot number)
        ! **********************************************************
        implicit real(kind = 8) (a-h, o-z)
        if(s>=0) then
            getFFPfromW = 370.0 + 3.0e-1 * s**1.45 ! FFP in MV
        else
            getFFPfromW = 370.0 - 3.0e-1 * abs(s)**1.45 ! FFP in MV
        endif
        return
    end

    ! **********************************************************
    function getRfromE(iz, ia, Ek, Em) ! get Rigidity in MV from Kinetic Energy (MeV/n)
        ! **********************************************************
        implicit real(kind = 8) (a-h, o-z)
        getRfromE = sqrt((ia * Ek)**2 + 2 * Ek * ia * Em) / iz
        return
    end

    ! **********************************************************
    function getEfromR(iz, Rm, COR) ! get Kinetic Energy (MeV) from Rigidity (MV)
        ! **********************************************************
        implicit real(kind = 8) (a-h, o-z)
        getEfromR = sqrt((iz * COR)**2 + Rm**2) - Rm
        return
    end

    ! *******************************************************
    function getPow(ip, d, r)  ! get Power of solar modulation dependence, r:Cut off rigidity, d:depth
        ! *******************************************************
        implicit real(kind = 8) (a-h, o-z)
        parameter (npart = 13) ! number of particle type, neutron, proton, alpha, photon, electron positron, mu+, mu-, Li-O
        parameter (nBdata = 2) ! B(1) - B(4) : Pow = b1 + b2*d
        parameter (nAdata = 5) ! A(1) - A(5) : B = A(1)+A(2)*r+A(3)/(1+exp((r-A(4))/A(5)))
        character chatmp1*1, chatmp5*5
        character, save :: pname(0:npart)*6
        real(kind = 8), save :: A(0:npart, nBdata, nAdata), B(nBdata)
        integer(kind = 4), save :: ifirst
        data ifirst/0/
        data pname/'neutro', 'proton', 'alphaa', 'elemag', 'elemag', 'elemag', 'muon--', 'muon--', &
                & 'ions  ', 'ions  ', 'ions  ', 'ions  ', 'ions  ', 'ions  '/
        call setlocale_to_avoid_bug()

        if(ifirst==0) then
            ifirst = 1
            do i = 0, 7
                if(i<=2) then ! neutron, proton, alpha
                    open(unit = 15, file = 'input/' // pname(i) // '/solar-dep.inp', status = 'old')
                elseif(i==3) then
                    open(unit = 15, file = 'input/' // pname(i) // '/solar-dep-EL.inp', status = 'old')
                elseif(i==4) then
                    open(unit = 15, file = 'input/' // pname(i) // '/solar-dep-PO.inp', status = 'old')
                elseif(i==5) then
                    open(unit = 15, file = 'input/' // pname(i) // '/solar-dep-PH.inp', status = 'old')
                elseif(i==6) then
                    open(unit = 15, file = 'input/' // pname(i) // '/solar-dep.plus', status = 'old')
                elseif(i==7) then
                    open(unit = 15, file = 'input/' // pname(i) // '/solar-dep.mins', status = 'old')
                endif
                read(15, '(a)') chatmp1
                do ib = 1, nBdata
                    read(15, *) (A(i, ib, ia), ia = 1, nAdata)
                enddo
                close(unit = 15)
            enddo
            open(unit = 15, file = 'input/ions/solar-dep.inp', status = 'old')
            read(15, '(a)') chatmp1
            do i = npart - 5, npart
                do ib = 1, nBdata
                    read(15, *) (A(i, ib, ia), ia = 1, nAdata)
                enddo
            enddo
            close(15)
        endif

        do ib = 1, nBdata
            b(ib) = a(ip, ib, 1) + a(ip, ib, 2) * r + a(ip, ib, 3) / (1.0 + exp((r - a(ip, ib, 4)) / a(ip, ib, 5)))
        enddo

        getpow = b(1) + b(2) * d
        return
    end


    ! **********************************************************
    function getBestR(is, r, d, ip) ! get best Rc data for high-altitude correction
        ! **********************************************************
        parameter(nBdata = 6)
        parameter(ndep = 26)
        parameter(npart = 3) ! high-altitude correction is necessary only for electron, positron, photon
        implicit real(kind = 8) (a-h, o-z)
        real(kind = 8), save :: A(0:npart, nBdata, ndep), B(nBdata)
        real(kind = 8), save :: dep(ndep)
        character chatmp1*1, chatmp5*5
        character pname(npart)*2
        data ifirst/0/
        data pname/'EL', 'PO', 'PH'/
        call setlocale_to_avoid_bug()

        if(ifirst==0) then
            a(:, :, :) = 0.0
            do i = 0, npart
                if(i==0) then ! neutron
                    open(unit = 15, file = 'input/neutro/bestR.inp', status = 'old')
                else ! electron, positron, photon
                    open(unit = 15, file = 'input/elemag/bestR-' // pname(i) // '.inp', status = 'old')
                endif
                read(15, '(A)') chatmp1
                do id = 1, ndep
                    read(15, *) dep(id), (A(i, ib, id), ib = 1, nBdata)
                enddo
                close(unit = 15)
            enddo
            ifirst = 1
        endif

        do id = 1, ndep
            if(d<dep(id)) exit
        enddo
        if(id==1) then
            ratio = 0.0
            id = 2
        elseif(id==ndep + 1) then
            ratio = 1.0
            id = ndep
        else
            ratio = (d - dep(id - 1)) / (dep(id) - dep(id - 1))
        endif

        do ib = 1, nBdata
            B(ib) = A(ip, ib, id - 1) + (A(ip, ib, id) - A(ip, ib, id - 1)) * ratio
        enddo

        if(is==1) then ! solar minimum
            getBestR = 10**(b(1) + b(2) * r + b(3) / r)
        else
            getBestR = 10**(b(4) + b(5) * r + b(6) / r)
        endif

        return

    end

    ! *******************************************************
    function getNeutspec(s, r1, d1, g, e) ! get neutron flux
        !     s:Wolf number
        !     r:cut off rigidity (GV)
        !     d:air depth (g/cm^2)
        !     g:local geometry parameter, 0=< g =< 1: water weight fraction, 10:no-earth, 100:blackhole, -10< g < 0: pilot, g < -10: cabin
        !     e:neutron energy (MeV)
        ! *******************************************************
        implicit real(kind = 8) (a-h, o-z)
        parameter(nA = 12) ! number of basic spectrum parameter
        parameter(nG = 6) ! number of geometry parameter
        real(kind = 8), save :: A(nA) ! basic spectrum parameter
        real(kind = 8), save :: geo(nG) ! geometry parameter
        data ifirst/0/
        data airbus/2.45/  ! weight of airbus340 (100t)
        data Eth/2.5e-8/  ! themal energy
        call setlocale_to_avoid_bug()

        if(ifirst==0) then ! first time, get universal parameter (i.e: independent of all parameters)
            !      Read A parameter
            open(unit = 3, file = 'input/neutro/fitting-lowspec.inp', status = 'old')
            read(3, '(a)') chatmp1
            read(3, 1001) (A(ia), ia = 1, nA) !A(4)&A(12) is s,r,d-dependence, so will be changed
            close(unit = 3)
            ifirst = 1
        endif
        1001    format(30e13.5)

        r = max(1.0, r1) ! secondary particle fluxes are the same for Rc<1GV
        d = max(0.15, d1) ! secondary particle fluxes are the same for d < 0.15 g/cm2

        !     get condition dependent parameters
        Fl = getFl(0, s, r, d)
        A(12) = getA12(r, d)
        A(4) = getA4(r, d)
        call getGpara(g, geo) ! obtain G parameters
        !     calculate flux (/cm^2/s/lethargy)
        x = e  ! x is energy

        evap = a(1) * (x / a(2))**a(3) * exp(-x / a(2))
        gaus = a(4) * exp(-(log10(x) - log10(a(5)))**2 / (2 * log10(a(6))**2))
        conti = a(7) * log10(x / a(8)) * (1 + tanh(a(9) * log10(x / a(10)))) * (1 - tanh(a(11) * log10(x / a(12))))

        basic = conti + evap + gaus

        basic = basic * CorrNeut(s, r, d, e) ! correction for high altitude

        fG = geo(1) + geo(2) * log10(x / geo(3)) * (1 - tanh(geo(4) * log10(x / geo(5))))
        if(g<0.0) then
            if(g<-10.0) then ! cabin
                gtmp = g + 10.0
            else ! pilot
                gtmp = g
            endif
            fG = fG * (abs(gtmp) - 10 * int(abs(gtmp) / 10.0)) / airbus  ! consider size of aircraft
        endif
        geofactor = 10.0**fG
        ther = geo(6) * (x / Eth)**2 * exp(-(x / Eth))

        getNeutspec = Fl * (basic * geofactor + ther) / e

        return
    end

    ! *******************************************************
    subroutine getGpara(g, geo) ! get surroudning environment parameters
        ! *******************************************************
        implicit real(kind = 8) (a-h, o-z)
        parameter(nG = 6) ! number of geometry parameter
        dimension geo(nG) ! geometry parameter
        real(kind = 8), save :: P(24) ! Input parameters read from input file, P(1)-P(3) from Geo-Dep, P(4)-P(14) from Water-Dep, P(15)-P(24) from aircraft-dep
        character chatmp1*1, chatmp4*4
        call setlocale_to_avoid_bug()

        1001    format(a4, 1x, 30e13.5)

        if(P(1)==0.0d0) then
            open(unit = 3, file = 'input/neutro/Geo-Dep.inp', status = 'old')
            read(3, '(a)') chatmp1
            read(3, *) p(1), p(2), p(3)
            close(Unit = 3)
            open(unit = 3, file = 'input/neutro/Water-Dep.inp', status = 'old')
            read(3, '(a)') chatmp1
            read(3, 1001) chatmp4, p(4), p(5), p(6)
            read(3, 1001) chatmp4, p(7), p(8), p(9)
            read(3, 1001) chatmp4, p(10), p(11), p(12), p(13), p(14)
            close(unit = 3)
            open(unit = 3, file = 'input/neutro/Aircraft-Dep.inp', status = 'old')
            read(3, '(a)') chatmp1
            read(3, *) (p(ip), ip = 15, 19)  ! for pilot location
            read(3, *) (p(ip), ip = 20, 24)  ! for passenger & small aircraft configuration
            close(Unit = 3)
        endif

        if(g>=10.0) then ! in semi-infite atmosphere
            do ig = 1, nG
                geo(ig) = 0.0
            enddo
            geo(3) = 1.0  ! if geo(3)=0, the value should be in NaN
            geo(5) = 1.0  ! if geo(5)=0, the value should be in NaN
        elseif(g>=0.0) then ! for normal ground case
            geo(1) = p(1)
            geo(2) = p(2)
            geo(3) = 10.0**(p(4) + p(5) / (p(6) + g))
            geo(4) = p(3)
            geo(5) = p(7) + p(8) * g + p(9) * g**2
            geo(6) = (p(10) + p(11) * exp(-p(12) * g)) / (1 + p(13) * exp(-p(14) * g))
        else   ! pilot or cabin location
            is = 14
            if(g<=-10.0) is = is + 5 ! for passenger & small aircraft configuration, skip 5 more data
            do i = 1, 5
                geo(i) = p(is + i)
            enddo
            geo(6) = 0.0  ! no thermal component
        endif
        return
    end

    ! *******************************************************
    function getA4(r, d)  ! get A4 value, s:solar modulation potential, r:Cut off rigidity, d:depth
        ! *******************************************************
        implicit real(kind = 8) (a-h, o-z)
        parameter (nBdata = 4) ! B(1) - B(4) : Fl= B(1)*(exp(-B(2)*d)-B(3)*exp(-B(4)*d))
        parameter (nAdata = 6) ! A(1) - A(6) : B = A(1)+A(2)/(1+exp((r-A(3))/A(4))), A(5):A(1) for APmax, A(6):A(3) for APmax
        character chatmp1*1, chatmp5*5
        real(kind = 8), save :: A(nAdata), B(nBdata)
        call setlocale_to_avoid_bug()

        if(B(2)==0.0d0) then ! first time
            open(unit = 15, file = 'input/neutro/Depth-Dep-mid.out', status = 'old')  ! read B(2)-B(4) (independent of s,r,d)
            read(15, '(a)') chatmp1
            read(15, '(a)') chatmp1
            read(15, *) tmp1, tmp2, B(2), B(3), B(4)
            close(unit = 15)
            open(unit = 15, file = 'input/neutro/Rigid-Dep.inp', status = 'old')
            do i = 1, 5  ! skip 5 line, 1 header line + 4 Bdata
                read(15, '(a)') chatmp1
            enddo
            read(15, 1010) chatmp5, (A(ia), ia = 1, nAdata)
            close(unit = 15)
        endif
        1010    format(a5, 30e13.5)

        B(1) = A(1) + A(2) * r + A(3) / (1 + exp((r - A(4)) / A(5)))
        getA4 = B(1) + B(2) * d / (1 + B(3) * exp(B(4) * d))

        return
    end

    ! *******************************************************
    function getA12(r, d)  ! get Fl value, s:solar modulation potential, r:Cut off rigidity, d:depth
        ! *******************************************************
        implicit real(kind = 8) (a-h, o-z)
        parameter (nBdata = 4) ! B(1) - B(4) : Fl= B(1)*(exp(-B(2)*d)-B(3)*exp(-B(4)*d))
        parameter (nAdata = 6) ! A(1) - A(6) : B = A(1)+A(2)/(1+exp((r-A(3))/A(4))), A(5):A(1) for APmax, A(6):A(3) for APmax
        character chatmp1*1, chatmp5*5
        real(kind = 8), save :: A(nBdata, nAdata), B(nBdata)
        call setlocale_to_avoid_bug()

        if(B(4)==0.0d0) then ! first time
            open(unit = 15, file = 'input/neutro/Depth-Dep-hig.out', status = 'old')  ! read B(2),B(4) (independent of s,r,d)
            read(15, '(a)') chatmp1
            read(15, '(a)') chatmp1
            read(15, *) tmp0, tmp1, tmp2, tmp3, B(4)
            close(unit = 15)
            open(unit = 15, file = 'input/neutro/Rigid-Dep.inp', status = 'old')
            do i = 1, 6  ! skip 5 line, 1 header line + 5 Bdata
                read(15, '(a)') chatmp1
            enddo
            do ib = 1, 3 ! for B1 to B3
                read(15, 1010) chatmp5, (A(ib, ia), ia = 1, nAdata)
            enddo
            close(unit = 15)
        endif
        1010 format(a5, 30es13.5)

        do ib = 1, 3
            B(ib) = A(ib, 1) + A(ib, 2) * r + A(ib, 3) / (1 + exp((r - A(ib, 4)) / A(ib, 5)))
        enddo

        getA12 = B(1) * (exp(-B(2) * d) + B(3) * exp(-B(4) * d))

        return
    end

    ! *******************************************************
    function CorrNeut(s, r, d, e)  ! get correction factor for high altitude data, only solar minimum and maximum
        ! *******************************************************
        parameter(nhensu = 9) ! number of hensu
        parameter(mpara = 5) ! number of parameters to COR dependence
        parameter(ndep = 26) ! number of depth
        parameter(nsol = 2) ! solar minimum and maximum
        implicit real(kind = 8) (a-h, o-z)
        dimension ainp(mpara, nhensu, ndep, nsol)
        dimension dep(ndep)
        dimension b(nhensu), c(nsol) ! temporary dimension
        dimension spot(nsol)
        character chatmp1*1
        save ainp, spot, ifirst, dep
        data ifirst/0/
        data spot/0.0, 150.0/
        call setlocale_to_avoid_bug()

        if(ifirst==0) then ! first time call
            !     Read parameters for depth-independent parameters
            ! open(15,file='correction-depth-rigid-final.inp',status='old') ! only high altidue correction mode
            open(15, file = 'input/neutro/correction-depth-rigid.inp', status = 'old') ! all altitude correction mode
            read(15, '(a)') chatmp1
            do ih = 1, nhensu
                do id = 1, ndep
                    read(15, *) itmp, dep(id), ((ainp(ip, ih, id, is), ip = 1, mpara), is = 1, nsol)
                enddo
            enddo
            close(15)
            ifirst = 1
        endif

        ! find depth ID
        do id = 1, ndep
            if(d<dep(id)) exit
        enddo
        if(id==1) then
            ratio = 0.0
            id = 2
        elseif(id==ndep + 1) then
            ratio = 1.0
            id = ndep
        else
            ratio = (d - dep(id - 1)) / (dep(id) - dep(id - 1))
        endif

        rc = max(1.0, r)
        do is = 1, nsol
            do ih = 1, nhensu ! determine 9 hensu used in the correction equation
                d1 = ainp(1, ih, id - 1, is) + Ainp(2, ih, id - 1, is) * rc + Ainp(3, ih, id - 1, is) / (1 + exp((rc - Ainp(4, ih, id - 1, is)) / Ainp(5, ih, id - 1, is)))
                d2 = ainp(1, ih, id - 0, is) + Ainp(2, ih, id - 0, is) * rc + Ainp(3, ih, id - 0, is) / (1 + exp((rc - Ainp(4, ih, id - 0, is)) / Ainp(5, ih, id - 0, is)))
                b(ih) = d1 + (d2 - d1) * ratio
            enddo
            C(is) = 10**(b(1) + (b(2) * log10(e) + b(3)) * (1 - tanh(b(4) * log10(e / b(5)))) + (b(6) * log10(e) + b(7)) * (1 + tanh(b(8) * log10(e / b(9)))))
        enddo

        ip = 0 ! always neutron
        pow = getPow(ip, d, r) ! get Power index
        A2 = (c(1) - c(2)) / (getFFPfromW(spot(1))**pow - getFFPfromW(spot(2))**pow)
        A1 = c(1) - A2 * getFFPfromW(spot(1))**pow
        CorrNeut = a1 + a2 * getFFPfromW(s)**pow

        return
    end


    ! ******************************************************
    function getMuonSpec(ip, s, r1, d1, e) ! get muon spectrum, ip=30:mu+, =31:mu-
        ! ******************************************************
        implicit real(kind = 8) (a-h, o-z)
        parameter (nPart = 2) ! Muon+ or Muon-
        parameter (nsol = 2) ! solar minimum and maximum
        parameter (nAdata = 7) ! number of A parameter

        dimension Acurr(nAdata)
        dimension Fl(nsol), spot(nsol)

        data restmass/105.6/
        data spot/0.0, 150.0/
        data ethre/3.0e5/ ! threshold energy for high-energy muon correction

        if(e<1.0e-2) then
            getMuonSpec = 0.0
            return
        endif

        r = max(1.0, r1) ! muon fluxes are the same for Rc < 1 GV
        d = max(0.15, d1) ! secondary particle fluxes are the same for d < 0.15 g/cm2

        beta = sqrt(1.0 - (restmass / (restmass + e))**2)
        tmp = max(2.0, log10(e)) ! below 100 MeV, this value should be constant

        do is = 1, nsol
            iptmp = ip ! 1 for mu+, 2 for mu+
            call getAmuon(Acurr, iptmp, is, d, r)
            if(e>ethre) then
                acurr(5) = acurr(5) + 0.4 ! /(1+exp((5.8088849-tmp)/0.24867240)) ! high energy correction, see fit/muon/PowerCorrection
                acurr(1) = acurr(1) * ethre**0.4 ! /(1+exp((5.8088849-tmp)/0.24867240))
            endif
            Fl(is) = acurr(1) * (e + (acurr(2) + acurr(4) * tmp) / beta**acurr(3))**(-acurr(5)) * (1 + exp(-acurr(6) * (log(e) + acurr(7))))
        enddo

        iptmp = ip + 5 ! 6 for mu+, 7 for mu-
        pow = getPow(iptmp, d, r) ! get Power index
        A2 = (Fl(1) - Fl(2)) / (getFFPfromW(spot(1))**pow - getFFPfromW(spot(2))**pow)
        A1 = Fl(1) - A2 * getFFPfromW(spot(1))**pow
        getMuonSpec = a1 + a2 * getFFPfromW(s)**pow

        return
    end

    ! ******************************************************
    subroutine getAmuon(Acurr, ip, is, d, r) ! get A parameter
        ! ******************************************************
        implicit real(kind = 8) (a-h, o-z)
        parameter (npart = 2) ! mu+ and mu-
        parameter (nsol = 2) ! solar minimum and maximum
        parameter (nAdata = 7) ! number of A parameter A(1) to A(7)
        parameter (nBdata = 10) ! A_min = B1+B2*r+B3/(1+exp((r-B4)/B5), A_max = B6+B7*r+B8/(1+exp((r-B9)/B10)
        parameter (ndep = 26) ! number of depth
        dimension Acurr(nAdata)
        dimension Bdata(nAdata, nBdata, npart, ndep), dep(ndep)
        character chatmp1*1, chatmp4*4
        character charge(npart)*4
        character chanum(0:9)*1
        save Bdata, dep
        data chanum/'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/
        data charge/'plus', 'mins'/
        data ifirst/0/
        call setlocale_to_avoid_bug()

        if(ifirst==0) then  ! first time called this routine
            ifirst = 1
            do ip2 = 1, npart
                open(15, file = 'input/muon--/final135.' // charge(ip2), status = 'old') ! read A(1),A(3),A(5), they are solar independent
                read(15, *) chatmp1
                do id2 = 1, ndep
                    read(15, *) dep(id2), Bdata(1, 1, ip2, id2), Bdata(3, 1, ip2, id2), Bdata(5, 1, ip2, id2)
                enddo
                close(15)
                open(15, file = 'input/muon--/final2467.' // charge(ip2), status = 'old') ! read A(2),A(4),A(6),A(7)
                read(15, *) chatmp1
                do id2 = 1, ndep
                    read(15, *) dep(id2), (Bdata(2, ib, ip2, id2), ib = 1, nBdata)
                enddo
                read(15, *) chatmp1
                do id2 = 1, ndep
                    read(15, *) dep(id2), (Bdata(4, ib, ip2, id2), ib = 1, nBdata)
                enddo
                read(15, *) chatmp1
                do id2 = 1, ndep
                    read(15, *) dep(id2), (Bdata(6, ib, ip2, id2), ib = 1, nBdata)
                enddo
                read(15, *) chatmp1
                do id2 = 1, ndep
                    read(15, *) dep(id2), (Bdata(7, ib, ip2, id2), ib = 1, nBdata)
                enddo
                close(15)
            enddo
        endif

        ! find closest depth
        do id = 1, ndep
            if(d<dep(id)) exit
        enddo
        if(id==1) then
            ratio = 0.0
            id = 2
        elseif(id==ndep + 1) then
            ratio = 1.0
            id = ndep
        else
            ratio = (d - dep(id - 1)) / (dep(id) - dep(id - 1))
        endif

        rc = max(1.0, r) ! minimum Rc = 1.0GV

        do ia = 1, nAdata
            if(ia==1) then
                Acurr(ia) = log(Bdata(ia, 1, ip, id - 1)) + (log(Bdata(ia, 1, ip, id)) - log(Bdata(ia, 1, ip, id - 1))) * ratio ! log-interpolation
                Acurr(ia) = exp(Acurr(ia))
            elseif(ia==3.or.ia==5) then
                Acurr(ia) = Bdata(ia, 1, ip, id - 1) + (Bdata(ia, 1, ip, id) - Bdata(ia, 1, ip, id - 1)) * ratio
            else
                if(is==1) then
                    B1 = Bdata(ia, 1, ip, id - 1) + Bdata(ia, 2, ip, id - 1) * rc + Bdata(ia, 3, ip, id - 1) / (1 + exp((rc - Bdata(ia, 4, ip, id - 1)) / Bdata(ia, 5, ip, id - 1)))
                    B2 = Bdata(ia, 1, ip, id) + Bdata(ia, 2, ip, id) * rc + Bdata(ia, 3, ip, id) / (1 + exp((rc - Bdata(ia, 4, ip, id)) / Bdata(ia, 5, ip, id)))
                else
                    B1 = Bdata(ia, 6, ip, id - 1) + Bdata(ia, 7, ip, id - 1) * rc + Bdata(ia, 8, ip, id - 1) / (1 + exp((rc - Bdata(ia, 9, ip, id - 1)) / Bdata(ia, 10, ip, id - 1)))
                    B2 = Bdata(ia, 6, ip, id) + Bdata(ia, 7, ip, id) * rc + Bdata(ia, 8, ip, id) / (1 + exp((rc - Bdata(ia, 9, ip, id)) / Bdata(ia, 10, ip, id)))
                endif
                Acurr(ia) = B1 + (B2 - B1) * ratio
            endif
        enddo

        return
    end


    ! *******************************************************
    function getIonSpec(iz, s, r, d, e) ! get Ion flux
        !     s:Wolf number
        !     r:cut off rigidity (GV)
        !     d:air depth (g/cm^2)
        !     e:ion energy (MeV/n)
        ! *******************************************************
        implicit real(kind = 8) (a-h, o-z)
        parameter(nAdata = 6)
        parameter(npart = 28)
        parameter(ngroup = 6)

        character chatmp1*1

        real(kind = 8), save :: A(nAdata, ngroup) ! parameter used in combine.for
        real(kind = 8), save :: dEdx(npart)
        integer(kind = 4), save :: iAnum(npart)
        integer(kind = 4), save :: ifirst

        integer, save :: igidx(npart)    ! group index

        data igidx/1, 2, 3, 3, 3, 4, 4, 4, 4, 5 &
                &, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6 &
                &, 6, 6, 6, 6, 6, 6, 6, 6/

        data ianum/ 1, 4, 7, 9, 11, 12, 14, 16, 19, 20, 23, 24, 27, 28, 31, 32, 35, 40, 39, 40, 45, 48, 51, 52, 55, 56, 59, 59/
        data ifirst/0/
        data restmass/938.27d0/
        call setlocale_to_avoid_bug()

        if(ifirst==0) then  ! first time call
            ifirst = 1
            open(unit = 15, file = 'input/ions/Combine.inp', status = 'old')
            read(15, *) chatmp1
            do i = 1, ngroup
                read(15, *) (A(ia, i), ia = 1, nAdata)
            enddo
            close(unit = 15)
        endif

        if(e<1.0e-2) then  ! for lower energy, no output
            getIonSpec = 0.0
            return
        endif

        x = e ! x is energy
        ig = igidx(iz)

        tmp = restmass * iAnum(iz)

        Ecut = getEfromR(iZ, tmp, r * 1000.0) / iAnum(iz) - a(6, ig) * d  ! MeV/n
        EcPri = max(a(1, ig), Ecut * a(3, ig))
        EcSec = max(a(2, ig), Ecut * a(3, ig))

        if(iz<=2) then
            ip = iz
        else
            ip = iz + 3 ! in getsecondary, ip=iz+3 (electron, positron, photon) for ions
        endif

        getIonSpec = getPrimary(iz, iAnum(iz), s, d, x) * 0.5 * (tanh(a(4, ig) * (x / EcPri - 1)) + 1.0) &
                & + getsecondary(ip, s, r, d, x) * 0.5 * (tanh(a(5, ig) * (1 - x / EcSec)) + 1.0)

        return

    end

    ! **********************************************************
    function getsecondary(ip, s, r1, d1, e) ! get secondary particle flux (/cm^2/s/MeV)
        ! **********************************************************
        parameter(nBdata = 8)
        parameter(ndep = 26)
        parameter(npart = 11) ! proton, alpha, electron, positron, photon, Li, Be, B, C, N, O
        implicit real(kind = 8) (a-h, o-z)
        real(kind = 8), save :: A(nBdata, npart, ndep), B(nBdata)
        real(kind = 8), save :: dep(ndep)
        integer, save :: ifirst
        character chatmp1*1, chatmp5*5
        data ifirst/0/
        call setlocale_to_avoid_bug()

        r = max(1.0, r1) ! to get Fl, minimum Rc=1GV
        d = max(0.15, d1) ! secondary particle fluxes are the same for d < 0.15 g/cm2

        getsecondary = 0.0

        if(ifirst==0) then
            do i = 1, 5 ! proton, alpha, electron, positron, photon
                if(i==1) then
                    open(unit = 15, file = 'input/proton/fitting-lowspec.inp', status = 'old')
                elseif(i==2) then
                    open(unit = 15, file = 'input/alphaa/fitting-lowspec.inp', status = 'old')
                elseif(i==3) then
                    open(unit = 15, file = 'input/elemag/fitting-lowspec-EL.inp', status = 'old')
                elseif(i==4) then
                    open(unit = 15, file = 'input/elemag/fitting-lowspec-PO.inp', status = 'old')
                elseif(i==5) then
                    open(unit = 15, file = 'input/elemag/fitting-lowspec-PH.inp', status = 'old')
                endif
                read(15, '(A)') chatmp1
                do id = 1, ndep
                    read(15, *) dep(id), (A(ib, i, id), ib = 1, nBdata)
                enddo
                close(unit = 15)
            enddo
            open(15, file = 'input/ions/fitting-lowspec.inp', status = 'old')
            read(15, '(a)') chatmp1
            do i = npart - 5, npart
                read(15, *) (A(ib, i, 1), ib = 1, nBdata)
                do ib = 1, nBdata
                    do id = 1, ndep
                        A(ib, i, id) = A(ib, i, 1) ! for Li to O, depth independent
                    enddo
                enddo
            enddo
            close(15)
            ifirst = 1
        endif

        if(e<1.0e-2.or.ip>npart) then  ! for lower energy or heavier ions, no output
            getsecondary = 0.0
            return
        endif

        do id = 1, ndep
            if(d<dep(id)) exit
        enddo
        if(id==1) then
            ratio = 0.0
            id = 2
        elseif(id==ndep + 1) then
            ratio = 1.0
            id = ndep
        else
            ratio = (d - dep(id - 1)) / (dep(id) - dep(id - 1))
        endif

        do ib = 1, nBdata
            B(ib) = A(ib, ip, id - 1) + (A(ib, ip, id) - A(ib, ip, id - 1)) * ratio
        enddo

        if(ip==1.or.ip==2) then ! proton or alpha
            getsecondary = getFl(ip, s, r, d) * (b(1) * e**b(2)) / (1 + b(3) * e**b(4)) / (1 + b(5) * e**b(6)) * (1 + exp(-b(7) * (log(e) + b(8))))
        elseif(ip==3.or.ip==4) then ! electron or positron
            getsecondary = getFl(ip, s, r, d) * (b(1) * e**b(2)) / (1 + b(3) * e**b(4)) / (1 + b(5) * e**b(6))
        elseif(ip==5) then ! photon
            getsecondary = getFl(ip, s, r, d) * (b(1) * e**b(2)) * (1 + b(3) * e**b(4)) / (1 + b(5) * e**b(6)) / (1 + exp(-b(7) * (log(e) + b(8))))
        else ! Li,Be,B,C,N,O
            getsecondary = getFl(ip, s, r, d) * (b(1) * e**b(2)) / (1 + b(3) * e**b(4)) / (1 + b(5) * e**b(6)) / (1 + exp(-b(7) * (log(e) + b(8))))
        endif

        return

    end

    ! **********************************************************
    function getTOAspec(iZ, iA, Ek, Spot) ! get TOA spectrum in (/(MeV/n)/s/m^2/sr)
        !     Ek: Kinetic Energy in MeV/n
        !     Spot: Wolf number estimated from count rate of neutron monitors
        ! **********************************************************
        parameter(npart = 28)
        implicit real(kind = 8) (a-h, o-z)
        real(kind = 8), save :: Dpara(npart), alpha(npart), gamma(npart), bpara(npart) ! Data for Nymmik Model
        data Emp/938.27d0/ ! mass of proton, nucleus mass is simply assumed to be A*Emp

        data Dpara/1.85e4, 3.69e3, 19.50, 17.70, 49.20, 103.00, 36.70, 87.40, &
                &           3.19, 16.40, 4.43, 19.30, 4.17, 13.40, 1.15, 3.06, 1.30, &
                &           2.33, 1.87, 2.17, 0.74, 2.63, 1.23, 2.12, 1.14, 9.32, 0.10, 0.48/ ! ISO-Model taken from Matthia-ASR2013

        data alpha/2.85, 3.12, 3.41, 4.30, 3.93, 3.18, 3.77, 3.11, 4.05, 3.11, 3.14, &
                &           3.65, 3.46, 3.00, 4.04, 3.30, 4.40, 4.33, 4.49, 2.93, 3.78, 3.79, &
                &           3.50, 3.28, 3.29, 3.01, 4.25, 3.52/

        data gamma/2.74, 2.77, 2.82, 3.05, 2.96, 2.76, 2.89, 2.70, 2.82, 2.76, 2.84, &
                &           2.70, 2.77, 2.66, 2.89, 2.71, 3.00, 2.93, 3.05, 2.77, 2.97, 2.99, &
                &           2.94, 2.89, 2.74, 2.63, 2.63, 2.63/

        ! ***** Determine LIS spectra based on DLR Model *****************
        R = getRfromE(iz, ia, Ek, Emp * iA) * 0.001 ! rigidity in GV
        beta = sqrt(1 - (Emp * iA / (Emp * iA + Ek * iA))**2)
        dR2dE = 0.001 / iZ / beta * iA ! convert GV to MV, MeV to MeV/n
        SpecLIS = Dpara(iz) * beta**alpha(iz) / R**gamma(iz) * dR2dE
        ! *******************************************************************

        ! consider solar modulation
        R0 = getFFPfromW(Spot) * 0.001 ! FFP in GV from W, taken from Matthia-ASR2013
        delta = 0.02 * Spot + 4.7        ! taken from Matthia-ASR2013
        getTOAspec = SpecLIS * (R / (R + R0))**delta

        return
    end


    ! *************************************************************
    function getPrimary(iz, ia, s, d, e)  ! get Primary Flux
        ! *************************************************************
        parameter(nAdata = 3)
        parameter(npart = 28)
        parameter(nepoint = 6)
        parameter(ngroup = 6)
        parameter(ndep = 26)
        implicit real(kind = 8) (a-h, o-z)
        real(kind = 8), save :: A(nAdata, ngroup, ndep)
        real(kind = 8), save :: dEdxTable(npart, nepoint), epoint(nepoint) ! dE/dx for each particle
        real(kind = 8), save :: dep(ndep) ! depth (g/cm2)
        real(kind = 8), save :: down
        integer, save :: igidx(npart)    ! group index
        dimension B(nAdata) ! temporary used dimension

        character gname(ngroup)*2
        character chatmp1*1

        data down/0.0/

        data gname/'H-', 'He', 'Be', 'N-', 'Si', 'Fe'/

        data igidx/1, 2, 3, 3, 3, 4, 4, 4, 4, 5 &
                &, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6 &
                &, 6, 6, 6, 6, 6, 6, 6, 6/
        call setlocale_to_avoid_bug()

        if(down==0.0d0) then ! first time call this routine
            down = 1.720313d0
            do ig = 1, ngroup
                open(unit = 15, file = 'input/ions/primary-' // gname(ig) // '.inp', status = 'old')
                read(15, *) chatmp1
                do id = 1, ndep
                    read(15, *) dep(id), (a(i, ig, id), i = 1, nAdata)
                enddo
                close(15)
            enddo
            open(15, file = 'input/ions/dEdx-table.inp', status = 'old')
            read(15, '(a1)') chatmp1
            read(15, '(a1)') chatmp1
            read(15, *) (epoint(ie), ie = 1, nepoint)
            do ip = 1, npart
                read(15, *) itmp, itmp, (dEdxTable(ip, ie), ie = 1, nepoint)
            enddo
            close(15)
        endif

        ig = igidx(iz)

        do id = 1, ndep
            if(d<dep(id)) exit
        enddo
        if(id==1) then
            ratio = 0.0
            id = 2
        elseif(id==ndep + 1) then
            ratio = 1.0
            id = ndep
        else
            ratio = (d - dep(id - 1)) / (dep(id) - dep(id - 1))
        endif

        do i = 1, nAdata
            b(i) = a(i, ig, id - 1) + (a(i, ig, id) - a(i, ig, id - 1)) * ratio
        enddo

        ! find dEdx
        do ie = 2, nepoint - 1
            if(e<=epoint(ie)) exit
        enddo
        ratio = (log(e) - log(epoint(ie - 1))) / (log(epoint(ie)) - log(epoint(ie - 1))) ! log-logg interpolation
        ratio = max(0.0, min(1.0, ratio))
        tmp = log(dEdxTable(iz, ie - 1)) + (log(dEdxTable(iz, ie)) - log(dEdxTable(iz, ie - 1))) * ratio ! log-log interpolation
        dEdx = exp(tmp)

        Eini = e + dEdx * d  ! Energy at the TOA
        getPrimary = getTOAspec(iz, ia, Eini, s) * (b(1) * exp(-b(2) * d) + (1.0 - b(1)) * exp(-b(3) * d))
        getPrimary = getPrimary * 4.0 * acos(-1.0) * 1.0e-4 / down  ! convert (/(MeV/n)/s/m^2/sr) to (/(MeV/n)/s/cm^2)

        return

    end

    ! *************************************************************
    function get511flux(s, r, d)  ! get 511 keV photon flux in (/cm2/s)
        ! *************************************************************
        parameter(ndep = 26)
        implicit real(kind = 8) (a-h, o-z)
        real(kind = 8), save :: F511(ndep), dep(ndep)
        character chatmp1*1

        data ifirst/0/
        call setlocale_to_avoid_bug()

        if(ifirst==0) then
            ifirst = 1
            open(15, file = 'input/elemag/flux511keV.inp')
            read(15, '(a1)') chatmp1
            do id = 1, ndep
                read(15, *) dep(id), F511(id)
            enddo
            close(15)
        endif

        do id = 1, ndep
            if(d<dep(id)) exit
        enddo
        if(id==1) then
            ratio = 0.0
            id = 2
        elseif(id==ndep + 1) then
            ratio = 1.0
            id = ndep
        else
            ratio = (d - dep(id - 1)) / (dep(id) - dep(id - 1))
        endif

        Fratio = F511(id - 1) + (F511(id) - F511(id - 1)) * ratio

        iptmp = 5 ! photon index
        ene = 0.511 ! 511 keV
        get511flux = getsecondary(iptmp, s, r, d, ene) * Fratio

        return

    end


    ! ******************************************************
    function getd(alti, cido) ! getd in g/cm^2, alti in km
        ! if -90 < cido < 90, use MSIS database
        ! else, use US standard air 1976
        ! ******************************************************
        parameter(iMSIS = 0)  ! 0:US standard atmosphere, 1:NRLMSISE database
        parameter(maxUS = 75) ! number of altitude bin for US-Standard Air
        parameter(maxMSIS = 129) ! number of altitude bin for NRLMSISE-00
        parameter(maxlat = 36) ! number of latitude bin for NRLMSISE-00
        implicit real(kind = 8) (a-h, o-z)
        real(kind = 8), save :: altUS(maxUS) ! altitude data for US-Standard Air 1976
        real(kind = 8), save :: altMSIS(maxMSIS) ! altitude data for NRLMSISE-00
        real(kind = 8), save :: depUS(maxUS) ! atmospheric depth data for US-Standard Air 1976
        real(kind = 8), save :: depMSIS(maxMSIS, maxlat) ! atmospheric depth data for each altitude & latitude for NRLMSISE-00
        real(kind = 8), save :: glat(maxlat) ! latitude data
        integer(kind = 4), save :: ifirst
        character chatmp1*1
        data ifirst/0/
        call setlocale_to_avoid_bug()

        if(ifirst==0) then ! come to this routine first, so read input data
            ifirst = 1
            open(unit = 15, file = 'input/AtomDepth.inp', status = 'old')
            read(15, '(a)') chatmp1
            do ia = 1, maxUS ! read US standard air 1976 data
                read(15, *) altUS(ia), depUS(ia)
            enddo
            read(15, '(a)') chatmp1
            read(15, *) (glat(ido), ido = 1, maxlat)
            do ia = 1, maxMSIS ! read NRLMSISE-00 data
                read(15, *) altMSIS(ia), (depMSIS(ia, ido), ido = 1, maxlat)
            enddo
        endif

        if((cido<-90.01.or.cido>90.01).or.iMSIS==0) then ! US standard atmosphere 1976 mode
            do ia = 1, maxUS
                if(altUS(ia)>alti) exit
            enddo
            if(ia==1) then ! out of range
                write(6, *) 'Error in function getd'
                write(6, *) 'Altitude =', alti, ' (km) should be higher than', altUS(1), ' (km)'
                stop
            endif
            if(ia==maxUS + 1) then ! out of range
                write(6, *) 'Warning in function getd'
                write(6, *) 'Altitude =', alti, ' (km) is too high. It is assumed to be', altUS(maxUS), ' (km)'
                getd = depUS(maxUS)
                return
            endif
            ratio = (alti - altUS(ia - 1)) / (altUS(ia) - altUS(ia - 1))
            getd = depUS(ia - 1) + ratio * (depUS(ia) - depUS(ia - 1))
        else ! latitude is specified, so use NRLMSISE-00 data
            do ia = 1, maxMSIS
                if(altMSIS(ia)>alti) exit
            enddo
            if(ia==1) then ! out of range
                write(6, *) 'Error in function getd'
                write(6, *) 'Altitude =', alti, ' (km) should be higher than', altMSIS(1), ' (km)'
                stop
            endif
            if(ia==maxMSIS + 1) then ! out of range
                write(6, *) 'Error in function getd'
                write(6, *) 'Altitude =', alti, ' (km) is too high. It is assumed to be', altMSIS(maxMSIS), ' (km)'
                alti = altMSIS(maxMSIS)
                ia = maxMSIS
            endif
            ratio = (alti - altMSIS(ia - 1)) / (altMSIS(ia) - altMSIS(ia - 1))
            do ido = 2, maxlat - 1
                if(glat(ido)>cido) exit
            enddo
            ratio1 = min(1.0, max(0.0, (cido - glat(ido - 1)) / (glat(ido) - glat(ido - 1))))
            dep1 = depMSIS(ia - 1, ido - 1) + ratio * (depMSIS(ia, ido - 1) - depMSIS(ia - 1, ido - 1))
            dep2 = depMSIS(ia - 1, ido) + ratio * (depMSIS(ia, ido) - depMSIS(ia - 1, ido))
            getd = dep1 + ratio1 * (dep2 - dep1)
        endif
        return
    end

    ! ******************************************************
    function getr(cido, ckei) ! cido and ckei are center of ido&keido of each grid
        use iso_c_binding
        ! ******************************************************
        implicit none

        real(kind = 8) :: getr, cor1, cor2
        real(kind = 8), intent(in) :: cido, ckei
        real(kind = 8), save :: cordata(181, 361) ! maximum 1 deg step, +1 mean
        real(kind = 8), save :: dpido(181), dpkei(361) ! data point ido & keido
        real(c_double) :: temp1, temp2, temp3
        character(len = 12) chatmp1, chatmp2, chatmp3
        integer(kind = 4), save :: mkei, mido
        real(kind = 8), save :: skei, sido
        integer(kind = 4), save :: ifirst
        integer(kind = 4), save :: ifirst2
        integer(kind = 4) :: ik, id
        data ifirst/0/
        data ifirst2/0/

        call setlocale_to_avoid_bug() ! very important for some systems

        if(ifirst==0) then ! come to this routine first, so read input data
            ifirst = 1
            open(unit = 99, file = 'input/CORdata.inp', status = 'old')
            read(99, *) mkei, mido ! read step size
            !write(*,*), "read values are : ", mkei,mido
            read(99, *) chatmp1, chatmp2, chatmp3
            !write(*,*), "read value is : ", chatmp1,chatmp2,chatmp3
            skei = 360.0 / (mkei - 1)  ! keido step
            sido = 180.0 / (mido - 1)  ! ido step
            do id = mido, 1, -1     ! read from 90 to -90 deg
                do ik = mkei, 1, -1      ! read from 180 to -180 deg

                    read(99, FMT = *) dpkei(ik), dpido(id), cordata(id, ik)

                    if (ifirst2==0) then
                        write(*, *) "Test of potential BUG : "
                        write(*, *) "the next values should be : 180.0 90.0 0.0"
                        write(*, *) dpkei(ik), dpido(id), cordata(id, ik)
                        if (dpkei(ik)/=180.0) then
                            write(*, *) "first value is not 180.0 -> Aborting."
                            call abort()
                        endif
                        ifirst2 = 1
                    endif

                enddo
            enddo
            close(unit = 99)
        endif

        ! **** Determine Cut-off Rigidity ***********
        id = min(mido - 1, int((cido + 90.0) / sido) + 1)   ! lower ido bin
        ik = min(mkei - 1, int((ckei + 180.0) / skei) + 1)  ! lower keido bin

        if (id<1 .or. id>180) then
            write(*, *), "id index out of range"
            call abort();
        else
            !write(*,*), "id is ok"
        endif

        if (ik<1 .or. ik>360) then
            write(*, *), "ik index out of range"
            call abort()
        else
            !write(*,*), "ik is ok"
        endif

        cor1 = cordata(id, ik) * (dpido(id + 1) - cido) / (dpido(id + 1) - dpido(id)) + cordata(id + 1, ik) * (cido - dpido(id)) / (dpido(id + 1) - dpido(id))
        cor2 = cordata(id, ik + 1) * (dpido(id + 1) - cido) / (dpido(id + 1) - dpido(id)) + cordata(id + 1, ik + 1) * (cido - dpido(id)) / (dpido(id + 1) - dpido(id))

        if (isnan(cor1)) then
            write(*, *), "cor1 is nan"
            call abort()
        endif

        if (isnan(cor2)) then
            write(*, *), "cor2 is nan"
            call abort()
        endif

        getr = cor1 * (dpkei(ik + 1) - ckei) / (dpkei(ik + 1) - dpkei(ik)) + cor2 * (ckei - dpkei(ik)) / (dpkei(ik + 1) - dpkei(ik))

        return

    end


    ! subroutines for getting angular distribution
    function getSpecAngFinal(ip, s, r, d, e, g, ang)
        ! ip: particle ID '1:neutro','2:proton','3:he---4','4:muon--','5:elepos','6:photon'
        ! s: W index
        ! r: cut-off rigidity in GV
        ! d: atmospheric depth in g/cm2
        ! e: energy in MeV/n
        ! g: local geometry effect
        ! ang: cos(theta)
        implicit real(kind = 8) (a-h, o-z)
        if(ip==4.and.g>=0.0d0.and.g<=1.0d0) then ! ground level muon
            emin = 1.1535d4 ! minimum energy for correction data
            if(e>=emin) then ! full correction
                getSpecAngFinal = getSpecAng(ip, s, r, d, e, g, 1.0d0) * getGmuon(e, ang)
            else
                ratio = (getSpecAng(ip, s, r, d, emin, g, 1.0d0) * getGmuon(emin, ang)) / getSpecAng(ip, s, r, d, emin, g, ang)
                getSpecAngFinal = getSpecAng(ip, s, r, d, e, g, ang) * ratio
            endif
        else ! no correction
            getSpecAngFinal = getSpecAng(ip, s, r, d, e, g, ang)
        endif

        if(g>=100.0d0) getSpecAngFinal = getSpecAngFinal * BHfactor(ip, e, ang) ! Consider black hole

    end

    function getSpecAng(ip, s, r, d, e, g, ang)
        ! ip: particle ID
        ! s: W index
        ! r: cut-off rigidity in GV
        ! d: atmospheric depth in g/cm2
        ! e: energy in MeV/n
        ! g: local geometry effect
        ! ang: cos(theta)
        parameter(npart = 6) ! number of particle type (1:neutron, 2:proton, 3:heavy ion, 4:muon, 5:electron&positron, 6:photon
        parameter(nsur = 18)  ! number of surface, upto 52 km
        parameter(ncor = 7)   ! number of cut-off ridigity upto 20 GV
        parameter(maxfit = 8) ! number of parameters to express angular distribution
        parameter(mpeach = 10) ! number of parameters to express energy dependence of each parameter
        parameter(ifour = 4)  ! 4

        implicit real(kind = 8) (a-h, o-z)

        real(kind = 8), save :: ParaEdep(mpeach, maxfit, ncor, nsur, npart) ! parameter for expressing energy-differential energy dependence
        real(kind = 8), save :: ParaEint(maxfit, ncor, nsur, npart) ! parameter for expressing energy-integrated energy dependence
        real(kind = 8), save :: depth(nsur) ! depth (g/cm2) for parameters
        real(kind = 8), save :: cor(ncor)   ! Rc (GV) for parameters
        real(kind = 8), save :: ParaAdep(maxfit, ifour) ! parameter for expressing angular distribution, 1-4 is for each depth & Rc condition
        real(kind = 8), save :: ratio1, ratio2 ! ratio must be saved
        real(kind = 8), save :: emin(npart), emax(npart) ! maximum energy
        character chatmp1*1
        character pname(npart)*6

        data pname/'neutro', 'proton', 'he---4', 'muon--', 'elepos', 'photon'/
        data emin/   1.0d-7, 1.0d0, 1.0d0, 1.0d1, 1.0d-1, 1.0d-2/      ! minimum energy
        data emax/    1.0d4, 1.0d4, 1.0d4, 1.0d5, 1.0d4, 1.0d4/      ! maximum energy
        data phimin/1.0e-3/ ! minimum value of phi

        data ifirst/0/
        data ipold/0/
        data sold/0.0/
        data rold/0.0/
        data dold/0.0/
        data eold/0.0/
        data gold/0.0/

        call setlocale_to_avoid_bug()

        getSpecAng = 1.0

        ! Read parameters
        if(ifirst==0) then ! first time call this routine
            ifirst = 1
            do ip1 = 1, npart
                open(15, file = 'input/angle/' // pname(ip1) // '.out', status = 'old')
                read(15, '(a1)') chatmp1
                do i = 1, maxfit
                    do is = 1, nsur
                        do ic = 1, ncor
                            read(15, *) itmp, depth(is), cor(ic), (ParaEdep(ii, i, ic, is, ip1), ii = 1, mpeach)
                        enddo
                    enddo
                enddo
                close(15)
                open(15, file = 'input/angle/' // pname(ip1) // '-Eint.out', status = 'old')
                read(15, '(a1)') chatmp1
                do is = 1, nsur
                    do ic = 1, ncor
                        read(15, *) tmp, tmp, (ParaEint(i, ic, is, ip1), i = 1, maxfit)
                    enddo
                enddo
                close(15)
            enddo
        endif

        ! check previous condition
        if(ip==ipold.and.s==sold.and.r==rold.and.d==dold.and.e==eold.and.g==gold) goto 10 ! same condition, need not determine parameter again
        ipold = ip
        sold = s
        rold = r
        dold = d
        eold = e
        gold = g
        ! Find depth
        do is = 1, nsur
            if(d<depth(is)) exit
        enddo
        if(is==1) then
            ratio1 = 0.0
            is = 2
        elseif(is==nsur + 1) then
            ratio1 = 1.0
            is = nsur
        else
            ratio1 = (d - depth(is - 1)) / (depth(is) - depth(is - 1))
        endif
        ! Find COR
        do ic = 1, ncor
            if(r<cor(ic)) exit
        enddo
        if(ic==1) then
            ratio2 = 0.0
            ic = 2
        elseif(ic==ncor + 1) then
            ratio2 = 1.0
            ic = ncor
        else
            ratio2 = (r - cor(ic - 1)) / (cor(ic) - cor(ic - 1))
        endif

        idx = 0
        do is1 = is - 1, is
            do ic1 = ic - 1, ic
                idx = idx + 1
                do i = 1, maxfit
                    if(e==0.0d0) then ! energy integrated
                        ParaAdep(i, idx) = ParaEint(i, ic1, is1, ip)
                    else
                        ene = max(min(e, emax(ip)), emin(ip))
                        ParaAdep(i, idx) = getParaAdep(i, ene, ParaEdep(1, i, ic1, is1, ip))
                        if(ip==1.and.g>=0.0d0.and.g<=1.0d0) ParaAdep(i, idx) = ParaAdep(i, idx) + getGneut(e, i) ! Ground level neutron correction, e instead of ene is used because data are down to 1e-8
                        if(i==2.and.ParaAdep(1, idx) + ParaAdep(2, idx)<phimin) ParaAdep(2, idx) = -ParaAdep(1, idx) + min(ParaAdep(1, idx), phimin) ! avoid negative value at cos(theta)=-1.0
                        if(i==5.and.ParaAdep(4, idx) + ParaAdep(5, idx)<phimin) ParaAdep(5, idx) = -ParaAdep(4, idx) + min(ParaAdep(4, idx), phimin) ! avoid negative value at cos(theta)=1.0
                    endif
                enddo
            enddo
        enddo

        do idx = 1, ifour
            call adjustParaAdep(ParaAdep(1, idx)) ! integration value is adjusted to 1
        enddo

        10 continue ! if all conditions are same, jump to here

        B1 = funcAng(ang, ParaAdep(1, 1))
        B2 = funcAng(ang, ParaAdep(1, 2))
        B3 = funcAng(ang, ParaAdep(1, 3))
        B4 = funcAng(ang, ParaAdep(1, 4))

        C1 = B1 + (B2 - B1) * ratio2 ! Rc interpolation
        C2 = B3 + (B4 - B3) * ratio2 ! Rc interpolation

        getSpecAng = C1 + (C2 - C1) * ratio1 ! Depth interpolation

        if(ip==4.and.d>60.and.e>=7000.0) then ! for debugging
            continue
        endif

        return
    end

    function getParaAdep(i, x, A)
        ! i: parameter index
        ! x: energy
        ! A: ParaEdep
        parameter(mpeach = 10) ! number of parameters to express energy dependence of each parameter
        implicit real(kind = 8) (a-h, o-z)
        dimension A(mpeach)

        a4 = max(0.01, a(4))
        a7 = max(0.01, a(7))
        a10 = max(0.01, a(10))

        getParaAdep = a(1) + a(2) / (1 + exp((a(3) - log10(x)) / a4)) + a(5) / (1 + exp((a(6) - log10(x)) / a7)) + a(8) / (1 + exp((a(9) - log10(x)) / a10))

        if(i==3.or.i==6) getParaAdep = 10.0**getParaAdep       ! a(3) & a(6) are determined by log10 fit
        if(i==7) getParaAdep = max(0.1d0, min(1.0d0, getParaAdep)) ! a(7) should be 0.1 < a7 1.0
        if(i==1.or.i==4.or.i==8) getParaAdep = max(0.0d0, getParaAdep)  ! a(1), a(4) & a(8) should be positive

        return

    end

    subroutine AdjustParaAdep(A)
        parameter (maxfit = 8)
        implicit real(kind = 8) (a-h, o-z)
        dimension A(maxfit)
        data one/1.0d0/
        data atlow/-0.05/  ! lower threshold angle for 90 degree interpolation
        data athig/ 0.05/  ! higher threshold angle for 90 degree interpolation

        sum1 = max(0.0, getint(abs(atlow), one, a(1), a(2), a(3)))                        ! -1.0 to atlow
        sum2 = max(0.0, (funcAng(atlow, a(1)) + funcAng(athig, a(1))) * 0.5 * (athig - atlow)) ! atlow to athig
        sum3 = max(0.0, getint(athig, a(7), a(4), a(5), a(6)))                           ! athig to a(7)
        if(a(7)<one) then
            sum4 = max(0.0, (funcAng(a(7), a(1)) + funcAng(one, a(1))) * 0.5 * (one - a(7)))           ! a(7) to 1.0
        else
            sum4 = 0.0
        endif
        sum = max(1.0e-20, (sum1 + sum2 + sum3 + sum4) * 2.0 * acos(-1.0))

        A(1) = A(1) / sum
        A(2) = A(2) / sum
        A(4) = A(4) / sum
        A(5) = A(5) / sum
        A(8) = A(8) / sum

        return
    end

    function getint(x0, x1, a1, a2, a3)
        implicit real(kind = 8) (a-h, o-z)
        getint = a1 * (x1 - x0) + a2 * (x1**(a3 + 1) - x0**(a3 + 1)) / (a3 + 1)
        return
    end

    function funcAng(x, a)
        parameter (maxfit = 8)
        implicit real(kind = 8) (a-h, o-z)
        dimension A(maxfit)
        data atlow/-0.05/  ! lower threshold angle for 90 degree interpolation
        data athig/ 0.05/  ! higher threshold angle for 90 degree interpolation

        if(x<=atlow) then ! backward
            funcAng = a(1) + a(2) * abs(x)**a(3)
        elseif(x<athig) then ! intermediate
            tmp1 = a(1) + a(2) * abs(atlow)**a(3)
            tmp2 = a(4) + a(5) * athig**a(6)
            funcAng = tmp1 + (tmp2 - tmp1) * (x - atlow) / (athig - atlow) ! simply interpolate
        elseif(x<=a(7)) then ! forward
            funcAng = a(4) + a(5) * x**a(6)
        else ! after peak
            tmp1 = a(4) + a(5) * a(7)**a(6)
            tmp2 = a(8)
            funcAng = tmp1 + (tmp2 - tmp1) * (x - a(7)) / (1.0d0 - a(7))
        endif

        return
    end

    function getGmuon(emid, ang)
        implicit real(kind = 8) (a-h, o-z)
        call setlocale_to_avoid_bug()

        elog = min(6.0, max(4.062, log10(emid))) ! parameters are effective only between 10GeV to 1 TeV
        a2 = 9.9873006E-01 + 2.9141114E+00 / (1.0 + exp((5.8030900E+00 - elog) / 2.4585039E-01))
        a3 = 2.4226042e1 - 1.5142933e1 * elog + 3.2012346 * elog**2 - 2.2325286e-1 * elog**3
        a4 = 1.4970401e1 - 5.3110524 * elog + 4.7458357e-1 * elog**2

        if(ang<0.0d0) then
            getGmuon = 0.0d0
        else
            x = 1.0d0 / max(0.001, ang)
            Scal = a2 * sqrt(1 - EXP(-((1.0 / a2)**2)))
            getGmuon = (a2 * sqrt(1 - EXP(-((x / a2)**2))) - (1 - EXP(-a4 * (x - 1)**a3))) / Scal
        endif

        return

    end

    function getGneut(emid, ID)
        parameter(maxfit = 8)   ! number of fitting parameter for angular dependence
        parameter(maxEfit = 3)  ! number of fitting parameter for energy dependence
        implicit real(kind = 8) (a-h, o-z)
        real(kind = 8), save :: Gneut(maxEfit, maxfit)
        data ifirst/0/
        call setlocale_to_avoid_bug()

        if(ifirst==0) then ! first time call this routine
            ! Read Ground Correction Factor
            open(15, file = 'input/angle/NeutronGround.out', status = 'old')
            do i = 1, maxfit
                read(15, *) (Gneut(ie, i), ie = 1, maxEfit)
            enddo
            close(15)
            ifirst = 1
        endif

        elog = max(-8.0, log10(emid)) ! parameters are effective only above 0.01 eV
        getGneut = Gneut(1, id) / (1 + exp((elog - Gneut(2, id)) / Gneut(3, id)))

        return

    end

    function BHfactor(ip, e, ang) ! Black hole factor
        implicit real(kind = 8) (a-h, o-z)
        parameter(npart = 6) ! only neutron, elepos, and photon
        parameter(nBHpara = 3) ! number of parameter
        parameter(nBHeach = 7) ! number of parameter to represent each parameter
        real(kind = 8), save :: BHpara(nBHeach, nBHpara, npart)
        dimension dimtmp(nBHpara) ! temporary used dimension
        character chatmp1*1, chatmp40*40
        character pname(npart)*6
        data pname/'neutro', 'proton', 'he---4', 'muon--', 'elepos', 'photon'/
        data ifirst/0/
        call setlocale_to_avoid_bug()

        if(ifirst==0) then ! first time call this routine
            ! Read Black hole factor parameter (only for neutron, elepos, photon)
            do ip2 = 1, npart
                open(15, file = 'input/angle/BkH-' // pname(ip2) // '.inp', status = 'old')
                read(15, '(a1)') chatmp1
                do i = 1, nBHpara
                    read(15, *) (BHpara(ii, i, ip2), ii = 1, nBHeach)
                enddo
                close(15)
            enddo
            ifirst = 1
        endif

        if(ang<0.0) then
            BHfactor = 0.0 ! backward is always 0
        else
            do i = 1, nBHpara
                dimtmp(i) = doublesig(e, BHpara(1, i, ip))
            enddo
            BHfactor = dimtmp(1) + (dimtmp(2) - dimtmp(1)) * ang**dimtmp(3)
        endif

    end

    function doublesig(e, a) ! get double sigmoid
        implicit real(kind = 8) (a-h, o-z)
        parameter(nBHeach = 7) ! number of parameter to represent each parameter
        dimension a(nBHeach)
        doublesig = a(1) + a(2) / (1 + exp((a(3) - log10(e)) / a(4))) + a(5) / (1 + exp((a(6) - log10(e)) / a(7)))
        return
    end

    function getHP(iy0, im0, id0, ic) ! get FFP from FFP tables
        !   ic=1: obtain FFP from neutron monitor data
        !	ic=2: obtain FFP from Wolf number
        !	ic=3: suspected ground level event
        !	ic=4: Too long time ago or future
        !	ic=5: no such date
        parameter(nmonth = 12)
        parameter(nday = 31)
        parameter(iymax = 2020) ! maximum year
        parameter(iymin = 1614) ! earliest year
        implicit real(kind = 8) (a-h, o-z)
        real, save :: FFP(iymin:iymax, nmonth, nday)
        real, save :: FFPuso(iymin:iymax) ! 0 & nmonth+1 data are used for interpolation
        data ifirst/0/
        integer(kind = 4), save :: iystart, iyend, iysUs, iyeUs
        call setlocale_to_avoid_bug()

        if(ifirst==0) then ! first time call this subroutine
            do iy = iymin, iymax  ! Initialized
                FFPuso(iy) = 0.0
                do im = 1, nmonth
                    do id = 1, nday
                        FFP(iy, im, id) = 0.0
                    enddo
                enddo
            enddo
            open(unit = 15, file = 'input/FFPtable.day', status = 'old')
            read(15, *) iystart, iyend
            do im = 1, nmonth
                do id = 1, nday
                    read(15, *) itmp1, itmp2, (FFP(iy, im, id), iy = iystart, iyend)
                enddo
            enddo
            close(unit = 15)
            open(unit = 15, file = 'input/FFPtable.uso', status = 'old')
            read(15, *) iysUs, iyeUs
            do iy = iysUs, iyeUs
                read(15, *) itmp, FFPuso(iy)
            enddo
            close(unit = 15)
            ifirst = 1
        endif

        ! ****** Year, month, day Check **************
        if(iy0<iymin.or.iy0>iymax) then  ! out of range
            ic = 4
            getHP = 0.0
            return
        endif
        if(im0<1.or.im0>nmonth.or.id0<1.or.id0>nday) then
            ic = 5
            getHP = 0.0
            return
        endif

        ! ******Determine FFP from Neutron Monitor *************
        if(FFP(iy0, im0, id0)>-99.0) then ! data exist
            if(iy0>=iystart.and.FFP(iy0, im0, id0)==-1000.0d0) then ! no such date
                ic = 5
                getHP = 0.0
            else  ! neutron monitor data exist
                getHP = FFP(iy0, im0, id0)
                if(getHP>1000.0) then
                    ic = 3 ! GLE occurred
                    getHP = getHP - 10000.0
                else
                    ic = 1 ! normal data
                endif
            endif
            return
        endif

        !  ***** Determine FFP from Usoskin's data  *****
        if(FFPuso(iy0)/=0.0d0) then
            getHP = FFPuso(iy0)
            ic = 2  ! determine FFP from Usoskin's data
            return
        endif

        !  **** No data **************************
        ic = 4
        getHP = 0.0
        return

    end

end module parma
