!
!  $Author: pkubota $
!  $Date: 2009/03/03 16:36:38 $
!  $Revision: 1.20 $
!
MODULE PhysicsDriver

    USE Constants, ONLY : &
            ityp, imon, icg, iwv, idp, ibd, &
            tov, &  ! intent(in)
            tdelt, &  ! intent(in)
            tstrat, &  ! intent(in)
            tdampr, &  ! intent(in)
            tdampf, &  ! intent(in)
            h0, &  ! intent(in)
            rlaps, &  ! intent(in)
            pie, &
            pai12, &
            cp, &
            hl, &
            gasr, &
            rk, &
            grav, &
            solcon, &
            stefan, &
            tf, &
            epsfac, &
            tice, &
            oceald, &
            icealn, &
            icealv, &
            con_rd, &
            con_rv, &
            EPS, &
            EPSM1, &
            i8, &
            r8

    USE ModRadiationDriver, Only : &
            RadiationDriver, radtim

    USE Sizes, ONLY : &
            kMax, &
            jMax, &
            ibMax, &
            ibMaxPerJB, &
            sl, &
            si, &
            del, &
            cl

    USE Surface, ONLY : &
            surface_driver

    USE SFC_SSiB, ONLY : &
            Albedo, Phenology

    USE SFC_SiB2, ONLY : &
            Albedo_sib2, Phenology_sib2

    USE Sfc_Ibis_Interface, Only : Albedo_IBIS

    USE GwddDriver, ONLY : &
            Gwdd_Driver

    USE Diagnostics, ONLY : &
            updia, dodia, &
            StartStorDiag, &
            nDiag_txgwds, & ! gravity wave drag surface zonal stress
            nDiag_tygwds, & ! gravity wave drag surface meridional stress
            nDiag_gwduzc, & ! gravity wave drag zonal momentum change
            nDiag_gwdvmc    ! gravity wave drag meridional momentum change

    USE PblDriver, ONLY : &
            Pbl_Driver

    USE SfcPBLDriver, ONLY : &
            SfcPBL_Driver

    USE Options, ONLY : &
            varcut, &
            dogwd, &
            mxrdcc, &
            lcnvl, &
            lthncl, &
            cdhl, &
            istrt, &
            first, &
            co2val, &
            delt, &
            filta, &
            nfin0, &
            nfin1, &
            initlz, &
            nfcnv0, &
            nfcldr, &
            iswrad, &
            ilwrad, &
            iccon, &
            swint, &
            trint, &
            yrl, &
            monl, &
            dtc3x, &
            epsflt, &
            intg, &
            maxtid, &
            dt, &
            idate, &
            idatec, &
            kt, &
            ktm, &
            ktp, &
            jdt, &
            start, &
            ifilt, &
            Model1D, &
            dirfNameOutput, &
            schemes, &
            crdcld, &!hmjb
            nscalars, &
            microphys, &
            specSfc, & !for specSfc
            nClass, &
            nAeros, &
            OCFLUX, &
            omlmodel

    USE Init, ONLY : &
            nls, &
            nlcs

    USE FieldsPhysics, ONLY : &
            ! Coeficiente de transporte vertical para a turbulencia
            PBL_CoefKm, PBL_CoefKh, &
            ! Albedo
            AlbNirBeam, AlbNirDiff, AlbVisBeam, AlbVisDiff, &
            ! Radiation fields at last integer hour
            rSwToaDown, &
            rVisDiff, rNirDiff, rVisBeam, rNirBeam, &
            rVisDiffC, rNirDiffC, rVisBeamC, rNirBeamC, &
            rSwSfcNet, rSwSfcNetC, &
            ! Radiation fields at next integer hour
            ySwToaDown, &
            yVisDiff, yNirDiff, yVisBeam, yNirBeam, &
            yVisDiffC, yNirDiffC, yVisBeamC, yNirBeamC, &
            ySwHeatRate, ySwHeatRateC, &
            ySwSfcNet, ySwSfcNetC, &
            ! LW Radiation fields at last integer hour
            LwCoolRate, LwSfcDown, LwSfcNet, LwToaUp, &
            LwCoolRateC, LwSfcDownC, LwSfcNetC, LwToaUpC, &
            ! Cloud field
            cldsav, cldtot, iMask, &
            cldinv, cldsat, cldcon, cldson, &
            ! Microphysics
            clwd, emisd, taud, EFFCS, EFFIS, tg0, &
            botlv, convb, QSfc0, TSfc0, QSfcm, TSfcm, &
            convbs, convc, convcs, convt, convts, evap, geshem, &
            ppci, ppli, prcc, prcp1, prcp2, &
            prcp3, prcpt, prct, concld, cld, cu_hr, cu_kbot, cu_ktop, cu_Kuo, &
            !Planetary boundary layer
            gl0, Mmlen, tauresx, tauresy, &
            !Convection
            dudt, dvdt, &
            !Surface
            gtsea, gndvi, qliq, sens, sheleg, sigki, ssib, bstar, tc0, &
            toplv, tseam, ustr, var, tm0, qm0, tmm, qmm, tcm, tgm, tdm, capacm, wm, &
            vstr, zorl, xland, lowlyr, ustar, z0, tkemyj, &
            sflux_t, sflux_r, sflux_u, sflux_v, tstar, wstar, r_aer, &
            veg_type, frac_occ, npatches, npatches_actual, &
            thz0, qz0, uz0, ndvim, MskAnt, &
            vz0, pblh, akhs, akms, ct, o3mix, snow, htdisp, temp2m, umes2m, &
            tpert, qpert, HML, HUML, HVML, TSK, z0sea, & !hmjb!hmjb
            sm0, mlsi ! add solange 27-01-2012

    USE FieldsDynamics, ONLY : &
            fgpass_scalars, adr_scalars !,   &
    !fgTsfc, fgH, fgLE, fgTau       !for specSfc (these vars should come through FieldsPhysics)
    !

    USE Convection, ONLY : &
            cumulus_driver

    USE PhysicalFunctions, ONLY : &
            fpvs

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: DryPhysics
    PUBLIC :: HumidPhysics
    PUBLIC :: InitSimpPhys
    PUBLIC :: SimpPhys
    REAL(KIND = r8), ALLOCATABLE :: teq(:, :, :)
    REAL(KIND = r8), ALLOCATABLE :: tauri(:)
    REAL(KIND = r8), ALLOCATABLE :: alfa(:)

CONTAINS

    SUBROUTINE DryPhysics &
            (zs, gt, gq, gu, gv, gps, gyu, gyv, gtd, &
            gqd, colrad, ifday, tod, gtt, gqq, omg, latco, &
            lonrad, glnpm, cos2d, intcosz, kMax, ibMax, &
            sgTsfc, sgH, sgLE, sgTau, &
            gicem, gicet, gliqm, gliqt, gvarm, gvart)
        INTEGER, INTENT(in) :: ibMax
        INTEGER, INTENT(in) :: kMax
        REAL(KIND = r8), INTENT(in) :: zs    (ibMax)
        REAL(KIND = r8), INTENT(inout) :: gt    (ibMax, kMax)
        REAL(KIND = r8), INTENT(inout) :: gq    (ibMax, kMax)
        REAL(KIND = r8), INTENT(in) :: gu    (ibMax, kMax)
        REAL(KIND = r8), INTENT(in) :: gv    (ibMax, kMax)
        REAL(KIND = r8), INTENT(in) :: gps   (ibMax)
        REAL(KIND = r8), INTENT(inout) :: gyu   (ibMax, kMax)
        REAL(KIND = r8), INTENT(inout) :: gyv   (ibMax, kMax)
        REAL(KIND = r8), INTENT(inout) :: gtd   (ibMax, kMax)
        REAL(KIND = r8), INTENT(inout) :: gqd   (ibMax, kMax)
        REAL(KIND = r8), INTENT(in) :: colrad(ibMax)
        REAL(KIND = r8), INTENT(in) :: lonrad(ibMax)
        INTEGER, INTENT(in) :: ifday
        REAL(KIND = r8), INTENT(in) :: tod
        REAL(KIND = r8), INTENT(in) :: gtt   (ibMax, kMax)
        REAL(KIND = r8), INTENT(in) :: gqq   (ibMax, kMax)
        REAL(KIND = r8), INTENT(in) :: omg   (ibMax, kMax)
        REAL(KIND = r8), OPTIONAL, INTENT(inout) :: gicem (ibMax, kMax)
        REAL(KIND = r8), OPTIONAL, INTENT(inout) :: gicet (ibMax, kMax)
        REAL(KIND = r8), OPTIONAL, INTENT(inout) :: gliqm (ibMax, kMax)
        REAL(KIND = r8), OPTIONAL, INTENT(inout) :: gliqt (ibMax, kMax)
        REAL(KIND = r8), OPTIONAL, INTENT(inout) :: gvarm (ibMax, kMax, nClass + nAeros)
        REAL(KIND = r8), OPTIONAL, INTENT(inout) :: gvart (ibMax, kMax, nClass + nAeros)
        ! variables for specified surface (specSfc)
        REAL(KIND = r8), OPTIONAL, INTENT(inout) :: sgTsfc(ibMax)
        REAL(KIND = r8), OPTIONAL, INTENT(inout) :: sgH(ibMax)
        REAL(KIND = r8), OPTIONAL, INTENT(inout) :: sgLE(ibMax)
        REAL(KIND = r8), OPTIONAL, INTENT(inout) :: sgTau(ibMax)

        INTEGER, INTENT(in) :: latco
        REAL(KIND = r8), INTENT(inout) :: glnpm (ibMax)
        REAL(KIND = r8), INTENT(IN) :: cos2d (ibMax)
        LOGICAL, INTENT(IN) :: intcosz
        REAL(KIND = r8) :: ps    (ibMax)
        INTEGER :: ibLim, i
        REAL(KIND = r8) :: topog     (ibMax)

        topog = zs / grav
        ps = glnpm
        ibLim = ibMax!ibMaxPerJB(latco)
        IF(schemes==1)THEN
            DO i = 1, ibLim
                IF(iMask(i, latco)== 13_i8)THEN
                    gndvi  (i, latco) = 0.0_r8
                    ndvim  (i, latco) = 0.0_r8
                END IF
            END DO
        ELSE IF(schemes==2)THEN
            DO i = 1, ibLim
                IF(iMask(i, latco)==13_i8)THEN
                    gndvi (i, latco) = 0.0_r8
                    ndvim (i, latco) = 0.0_r8
                END IF
            END DO
            !    ELSE IF(schemes==3)THEN
            !        DO i=1,nCols
            !          IF(iMask(i)==15)THEN
            !             ndvi  (i) = 0.000
            !             ndvim (i) = 0.000
            !          END IF
            !       END DO
        END IF
        IF (microphys) THEN
            IF((nClass + nAeros)>0 .and. PRESENT(gvarm))THEN
                IF (specSfc) THEN
                    CALL physcs  (&
                            ! Coeficiente de transporte vertical para a turbulencia
                            PBL_CoefKm(1:ibLim, :, latco), PBL_CoefKh(1:ibLim, :, latco), &
                            gt(1:ibLim, :), gq(1:ibLim, :), gu(1:ibLim, :), gv(1:ibLim, :), &
                            gps(1:ibLim), &
                            ppli(1:ibLim, latco), &
                            ppci(1:ibLim, latco), gyu(1:ibLim, :), gyv(1:ibLim, :), gtd(1:ibLim, :), &
                            gqd(1:ibLim, :), ySwHeatRate(1:ibLim, :, latco), LwCoolRate(1:ibLim, :, latco), &
                            imask(1:ibLim, latco), rVisBeam(1:ibLim, latco), rVisDiff(1:ibLim, latco), &
                            rNirBeam(1:ibLim, latco), rNirDiff(1:ibLim, latco), LwSfcDown(1:ibLim, latco), &
                            gtsea(1:ibLim, latco), colrad(1:ibLim), sl, si(1:kMax + 1), del, ifday, &
                            tod, AlbVisBeam(1:ibLim, latco), AlbVisDiff(1:ibLim, latco), &
                            AlbNirBeam(1:ibLim, latco), AlbNirDiff(1:ibLim, latco), rSwToaDown(1:ibLim, latco), &
                            LwSfcNet(1:ibLim, latco), LwToaUp(1:ibLim, latco), gl0(1:ibLim, latco), &
                            zorl(1:ibLim, latco), gtt(1:ibLim, :), gqq(1:ibLim, :), &
                            sheleg(1:ibLim, latco), tseam(1:ibLim, latco), omg(1:ibLim, :), &
                            ySwHeatRateC(1:ibLim, :, latco), rVisBeamC(1:ibLim, latco), rVisDiffC(1:ibLim, latco), &
                            rNirBeamC(1:ibLim, latco), rNirDiffC(1:ibLim, latco), LwSfcDownC(1:ibLim, latco), &
                            LwSfcNetC(1:ibLim, latco), LwToaUpC(1:ibLim, latco), &
                            convts(1:ibLim, latco), convcs(1:ibLim, latco), convbs(1:ibLim, latco), &
                            convc(1:ibLim, latco), convt(1:ibLim, latco), convb(1:ibLim, latco), &
                            ustr(1:ibLim, latco), vstr(1:ibLim, latco), latco, &
                            yVisBeam(1:ibLim, latco), yVisDiff(1:ibLim, latco), yNirBeam(1:ibLim, latco), &
                            yNirDiff(1:ibLim, latco), ySwToaDown(1:ibLim, latco), yVisBeamC(1:ibLim, latco), &
                            yVisDiffC(1:ibLim, latco), yNirBeamC(1:ibLim, latco), yNirDiffC(1:ibLim, latco), &
                            cldsav(1:ibLim, latco), ssib(1:ibLim, latco), bstar(1:ibLim, latco), ibMaxPerJB(latco), kMax, &
                            sigki, lonrad(1:ibLim), ps(1:ibLim), &
                            var(1:ibLim, latco), sens(1:ibLim, latco), evap(1:ibLim, latco), &
                            cos2d(1:ibLim), intcosz, LwCoolRateC(1:ibLim, :, latco), topog(1:ibLim), &
                            o3mix(1:ibLim, :, latco), Mmlen(1:ibLim, latco), &
                            xland(1:ibLim, latco), lowlyr(1:ibLim, latco), &
                            ustar(1:ibLim, latco), z0(1:ibLim, latco), tkemyj(1:ibLim, :, latco), &
                            thz0(1:ibLim, latco), qz0 (1:ibLim, latco), &
                            uz0 (1:ibLim, latco), vz0 (1:ibLim, latco), &
                            pblh(1:ibLim, latco), akhs(1:ibLim, latco), akms(1:ibLim, latco), &
                            CT(1:ibLim, latco), snow(1:ibLim, latco), &
                            htdisp(1:ibLim, latco), temp2m(1:ibLim, latco), umes2m(1:ibLim, latco), &
                            cldtot(1:ibLim, :, latco), cldinv(1:ibLim, :, latco), cldsat(1:ibLim, :, latco), &
                            cldcon(1:ibLim, :, latco), cldson(1:ibLim, :, latco), clwd  (1:ibLim, :, latco), &
                            emisd (1:ibLim, :, latco), taud  (1:ibLim, :, latco), &
                            ySwSfcNet(1:ibLim, latco), ySwSfcNetC(1:ibLim, latco), &
                            rSwSfcNet(1:ibLim, latco), rSwSfcNetC(1:ibLim, latco), mskant(1:ibLim, latco), &
                            sm0(1:ibLim, :, latco), mlsi(1:ibLim, latco), & ! add solange 27-01-2012nClass
                            gndvi(1:ibLim, latco), ndvim(1:ibLim, latco), qliq(1:ibLim, :, latco), tpert(1:ibLim, latco), &
                            qpert(1:ibLim, latco), &
                            sflux_t(1:ibLim, latco), sflux_r(1:ibLim, latco), sflux_u(1:ibLim, latco), sflux_v(1:ibLim, latco), &
                            tstar(1:ibLim, latco), wstar(1:ibLim, latco), r_aer(1:ibLim, latco), veg_type(1:ibLim, latco, npatches), &
                            frac_occ(1:ibLim, latco, npatches), HML (1:ibLim, latco), HUML(1:ibLim, latco), &
                            HVML(1:ibLim, latco), TSK(1:ibLim, latco), z0sea(1:ibLim, latco), &
                            cu_hr (1:ibLim, :, latco), cu_kbot(1:ibLim, latco), cu_ktop(1:ibLim, latco), cu_Kuo(1:ibLim, latco), &
                            dudt     (1:ibLim, 1:kMax, latco), dvdt     (1:ibLim, 1:kMax, latco), &
                            EFFCS(1:ibLim, 1:kMax, latco), EFFIS (1:ibLim, 1:kMax, latco), tauresx(1:ibLim, latco), tauresy(1:ibLim, latco), &
                            sgTsfc(1:ibLim), sgH(1:ibLim), sgLE(1:ibLim), sgTau(1:ibLim), &
                            gicem(1:ibLim, :), gicet(1:ibLim, :), gliqm(1:ibLim, :), gliqt(1:ibLim, :), gvarm(1:ibLim, :, :), gvart(1:ibLim, :, :))
                ELSE !specSfc
                    CALL physcs  (&
                            ! Coeficiente de transporte vertical para a turbulencia
                            PBL_CoefKm(1:ibLim, :, latco), PBL_CoefKh(1:ibLim, :, latco), &
                            gt(1:ibLim, :), gq(1:ibLim, :), gu(1:ibLim, :), gv(1:ibLim, :), &
                            gps(1:ibLim), &
                            ppli(1:ibLim, latco), &
                            ppci(1:ibLim, latco), gyu(1:ibLim, :), gyv(1:ibLim, :), gtd(1:ibLim, :), &
                            gqd(1:ibLim, :), ySwHeatRate(1:ibLim, :, latco), LwCoolRate(1:ibLim, :, latco), &
                            imask(1:ibLim, latco), rVisBeam(1:ibLim, latco), rVisDiff(1:ibLim, latco), &
                            rNirBeam(1:ibLim, latco), rNirDiff(1:ibLim, latco), LwSfcDown(1:ibLim, latco), &
                            gtsea(1:ibLim, latco), colrad(1:ibLim), sl, si(1:kMax + 1), del, ifday, &
                            tod, AlbVisBeam(1:ibLim, latco), AlbVisDiff(1:ibLim, latco), &
                            AlbNirBeam(1:ibLim, latco), AlbNirDiff(1:ibLim, latco), rSwToaDown(1:ibLim, latco), &
                            LwSfcNet(1:ibLim, latco), LwToaUp(1:ibLim, latco), gl0(1:ibLim, latco), &
                            zorl(1:ibLim, latco), gtt(1:ibLim, :), gqq(1:ibLim, :), &
                            sheleg(1:ibLim, latco), tseam(1:ibLim, latco), omg(1:ibLim, :), &
                            ySwHeatRateC(1:ibLim, :, latco), rVisBeamC(1:ibLim, latco), rVisDiffC(1:ibLim, latco), &
                            rNirBeamC(1:ibLim, latco), rNirDiffC(1:ibLim, latco), LwSfcDownC(1:ibLim, latco), &
                            LwSfcNetC(1:ibLim, latco), LwToaUpC(1:ibLim, latco), &
                            convts(1:ibLim, latco), convcs(1:ibLim, latco), convbs(1:ibLim, latco), &
                            convc(1:ibLim, latco), convt(1:ibLim, latco), convb(1:ibLim, latco), &
                            ustr(1:ibLim, latco), vstr(1:ibLim, latco), latco, &
                            yVisBeam(1:ibLim, latco), yVisDiff(1:ibLim, latco), yNirBeam(1:ibLim, latco), &
                            yNirDiff(1:ibLim, latco), ySwToaDown(1:ibLim, latco), yVisBeamC(1:ibLim, latco), &
                            yVisDiffC(1:ibLim, latco), yNirBeamC(1:ibLim, latco), yNirDiffC(1:ibLim, latco), &
                            cldsav(1:ibLim, latco), ssib(1:ibLim, latco), bstar(1:ibLim, latco), ibMaxPerJB(latco), kMax, &
                            sigki, lonrad(1:ibLim), ps(1:ibLim), &
                            var(1:ibLim, latco), sens(1:ibLim, latco), evap(1:ibLim, latco), &
                            cos2d(1:ibLim), intcosz, LwCoolRateC(1:ibLim, :, latco), topog(1:ibLim), &
                            o3mix(1:ibLim, :, latco), Mmlen(1:ibLim, latco), &
                            xland(1:ibLim, latco), lowlyr(1:ibLim, latco), &
                            ustar(1:ibLim, latco), z0(1:ibLim, latco), tkemyj(1:ibLim, :, latco), &
                            thz0(1:ibLim, latco), qz0 (1:ibLim, latco), &
                            uz0 (1:ibLim, latco), vz0 (1:ibLim, latco), &
                            pblh(1:ibLim, latco), akhs(1:ibLim, latco), akms(1:ibLim, latco), &
                            CT(1:ibLim, latco), snow(1:ibLim, latco), &
                            htdisp(1:ibLim, latco), temp2m(1:ibLim, latco), umes2m(1:ibLim, latco), &
                            cldtot(1:ibLim, :, latco), cldinv(1:ibLim, :, latco), cldsat(1:ibLim, :, latco), &
                            cldcon(1:ibLim, :, latco), cldson(1:ibLim, :, latco), clwd  (1:ibLim, :, latco), &
                            emisd (1:ibLim, :, latco), taud  (1:ibLim, :, latco), &
                            ySwSfcNet(1:ibLim, latco), ySwSfcNetC(1:ibLim, latco), &
                            rSwSfcNet(1:ibLim, latco), rSwSfcNetC(1:ibLim, latco), mskant(1:ibLim, latco), &
                            sm0(1:ibLim, :, latco), mlsi(1:ibLim, latco), & ! add solange 27-01-2012nClass
                            gndvi(1:ibLim, latco), ndvim(1:ibLim, latco), qliq(1:ibLim, :, latco), tpert(1:ibLim, latco), &
                            qpert(1:ibLim, latco), &
                            sflux_t(1:ibLim, latco), sflux_r(1:ibLim, latco), sflux_u(1:ibLim, latco), sflux_v(1:ibLim, latco), &
                            tstar(1:ibLim, latco), wstar(1:ibLim, latco), r_aer(1:ibLim, latco), veg_type(1:ibLim, latco, npatches), &
                            frac_occ(1:ibLim, latco, npatches), HML (1:ibLim, latco), HUML(1:ibLim, latco), &
                            HVML(1:ibLim, latco), TSK(1:ibLim, latco), z0sea(1:ibLim, latco), &
                            cu_hr (1:ibLim, :, latco), cu_kbot(1:ibLim, latco), cu_ktop(1:ibLim, latco), cu_Kuo(1:ibLim, latco), &
                            dudt     (1:ibLim, 1:kMax, latco), dvdt     (1:ibLim, 1:kMax, latco), &
                            EFFCS(1:ibLim, 1:kMax, latco), EFFIS (1:ibLim, 1:kMax, latco), tauresx(1:ibLim, latco), tauresy(1:ibLim, latco), &
                            gicem(1:ibLim, :), gicet(1:ibLim, :), gliqm(1:ibLim, :), gliqt(1:ibLim, :), gvarm(1:ibLim, :, :), gvart(1:ibLim, :, :))
                ENDIF !specSfc
            ELSE
                IF (specSfc) THEN
                    CALL physcs  (&
                            ! Coeficiente de transporte vertical para a turbulencia
                            PBL_CoefKm(1:ibLim, :, latco), PBL_CoefKh(1:ibLim, :, latco), &
                            gt(1:ibLim, :), gq(1:ibLim, :), gu(1:ibLim, :), gv(1:ibLim, :), &
                            gps(1:ibLim), &
                            ppli(1:ibLim, latco), &
                            ppci(1:ibLim, latco), gyu(1:ibLim, :), gyv(1:ibLim, :), gtd(1:ibLim, :), &
                            gqd(1:ibLim, :), ySwHeatRate(1:ibLim, :, latco), LwCoolRate(1:ibLim, :, latco), &
                            imask(1:ibLim, latco), rVisBeam(1:ibLim, latco), rVisDiff(1:ibLim, latco), &
                            rNirBeam(1:ibLim, latco), rNirDiff(1:ibLim, latco), LwSfcDown(1:ibLim, latco), &
                            gtsea(1:ibLim, latco), colrad(1:ibLim), sl, si(1:kMax + 1), del, ifday, &
                            tod, AlbVisBeam(1:ibLim, latco), AlbVisDiff(1:ibLim, latco), &
                            AlbNirBeam(1:ibLim, latco), AlbNirDiff(1:ibLim, latco), rSwToaDown(1:ibLim, latco), &
                            LwSfcNet(1:ibLim, latco), LwToaUp(1:ibLim, latco), gl0(1:ibLim, latco), &
                            zorl(1:ibLim, latco), gtt(1:ibLim, :), gqq(1:ibLim, :), &
                            sheleg(1:ibLim, latco), tseam(1:ibLim, latco), omg(1:ibLim, :), &
                            ySwHeatRateC(1:ibLim, :, latco), rVisBeamC(1:ibLim, latco), rVisDiffC(1:ibLim, latco), &
                            rNirBeamC(1:ibLim, latco), rNirDiffC(1:ibLim, latco), LwSfcDownC(1:ibLim, latco), &
                            LwSfcNetC(1:ibLim, latco), LwToaUpC(1:ibLim, latco), &
                            convts(1:ibLim, latco), convcs(1:ibLim, latco), convbs(1:ibLim, latco), &
                            convc(1:ibLim, latco), convt(1:ibLim, latco), convb(1:ibLim, latco), &
                            ustr(1:ibLim, latco), vstr(1:ibLim, latco), latco, &
                            yVisBeam(1:ibLim, latco), yVisDiff(1:ibLim, latco), yNirBeam(1:ibLim, latco), &
                            yNirDiff(1:ibLim, latco), ySwToaDown(1:ibLim, latco), yVisBeamC(1:ibLim, latco), &
                            yVisDiffC(1:ibLim, latco), yNirBeamC(1:ibLim, latco), yNirDiffC(1:ibLim, latco), &
                            cldsav(1:ibLim, latco), ssib(1:ibLim, latco), bstar(1:ibLim, latco), ibMaxPerJB(latco), kMax, &
                            sigki, lonrad(1:ibLim), ps(1:ibLim), &
                            var(1:ibLim, latco), sens(1:ibLim, latco), evap(1:ibLim, latco), &
                            cos2d(1:ibLim), intcosz, LwCoolRateC(1:ibLim, :, latco), topog(1:ibLim), &
                            o3mix(1:ibLim, :, latco), Mmlen(1:ibLim, latco), &
                            xland(1:ibLim, latco), lowlyr(1:ibLim, latco), &
                            ustar(1:ibLim, latco), z0(1:ibLim, latco), tkemyj(1:ibLim, :, latco), &
                            thz0(1:ibLim, latco), qz0 (1:ibLim, latco), &
                            uz0 (1:ibLim, latco), vz0 (1:ibLim, latco), &
                            pblh(1:ibLim, latco), akhs(1:ibLim, latco), akms(1:ibLim, latco), &
                            CT(1:ibLim, latco), snow(1:ibLim, latco), &
                            htdisp(1:ibLim, latco), temp2m(1:ibLim, latco), umes2m(1:ibLim, latco), &
                            cldtot(1:ibLim, :, latco), cldinv(1:ibLim, :, latco), cldsat(1:ibLim, :, latco), &
                            cldcon(1:ibLim, :, latco), cldson(1:ibLim, :, latco), clwd  (1:ibLim, :, latco), &
                            emisd (1:ibLim, :, latco), taud  (1:ibLim, :, latco), &
                            ySwSfcNet(1:ibLim, latco), ySwSfcNetC(1:ibLim, latco), &
                            rSwSfcNet(1:ibLim, latco), rSwSfcNetC(1:ibLim, latco), mskant(1:ibLim, latco), &
                            sm0(1:ibLim, :, latco), mlsi(1:ibLim, latco), & ! add solange 27-01-2012nClass
                            gndvi(1:ibLim, latco), ndvim(1:ibLim, latco), qliq(1:ibLim, :, latco), tpert(1:ibLim, latco), &
                            qpert(1:ibLim, latco), &
                            sflux_t(1:ibLim, latco), sflux_r(1:ibLim, latco), sflux_u(1:ibLim, latco), sflux_v(1:ibLim, latco), &
                            tstar(1:ibLim, latco), wstar(1:ibLim, latco), r_aer(1:ibLim, latco), veg_type(1:ibLim, latco, npatches), &
                            frac_occ(1:ibLim, latco, npatches), HML (1:ibLim, latco), HUML(1:ibLim, latco), &
                            HVML(1:ibLim, latco), TSK(1:ibLim, latco), z0sea(1:ibLim, latco), &
                            cu_hr (1:ibLim, :, latco), cu_kbot(1:ibLim, latco), cu_ktop(1:ibLim, latco), cu_Kuo(1:ibLim, latco), &
                            dudt     (1:ibLim, 1:kMax, latco), dvdt     (1:ibLim, 1:kMax, latco), &
                            EFFCS(1:ibLim, 1:kMax, latco), EFFIS (1:ibLim, 1:kMax, latco), tauresx(1:ibLim, latco), tauresy(1:ibLim, latco), &
                            sgTsfc(1:ibLim), sgH(1:ibLim), sgLE(1:ibLim), sgTau(1:ibLim), &
                            gicem(1:ibLim, :), gicet(1:ibLim, :), gliqm(1:ibLim, :), gliqt(1:ibLim, :))
                ELSE  !for specSfc
                    CALL physcs  (&
                            ! Coeficiente de transporte vertical para a turbulencia
                            PBL_CoefKm(1:ibLim, :, latco), PBL_CoefKh(1:ibLim, :, latco), &
                            gt(1:ibLim, :), gq(1:ibLim, :), gu(1:ibLim, :), gv(1:ibLim, :), &
                            gps(1:ibLim), &
                            ppli(1:ibLim, latco), &
                            ppci(1:ibLim, latco), gyu(1:ibLim, :), gyv(1:ibLim, :), gtd(1:ibLim, :), &
                            gqd(1:ibLim, :), ySwHeatRate(1:ibLim, :, latco), LwCoolRate(1:ibLim, :, latco), &
                            imask(1:ibLim, latco), rVisBeam(1:ibLim, latco), rVisDiff(1:ibLim, latco), &
                            rNirBeam(1:ibLim, latco), rNirDiff(1:ibLim, latco), LwSfcDown(1:ibLim, latco), &
                            gtsea(1:ibLim, latco), colrad(1:ibLim), sl, si(1:kMax + 1), del, ifday, &
                            tod, AlbVisBeam(1:ibLim, latco), AlbVisDiff(1:ibLim, latco), &
                            AlbNirBeam(1:ibLim, latco), AlbNirDiff(1:ibLim, latco), rSwToaDown(1:ibLim, latco), &
                            LwSfcNet(1:ibLim, latco), LwToaUp(1:ibLim, latco), gl0(1:ibLim, latco), &
                            zorl(1:ibLim, latco), gtt(1:ibLim, :), gqq(1:ibLim, :), &
                            sheleg(1:ibLim, latco), tseam(1:ibLim, latco), omg(1:ibLim, :), &
                            ySwHeatRateC(1:ibLim, :, latco), rVisBeamC(1:ibLim, latco), rVisDiffC(1:ibLim, latco), &
                            rNirBeamC(1:ibLim, latco), rNirDiffC(1:ibLim, latco), LwSfcDownC(1:ibLim, latco), &
                            LwSfcNetC(1:ibLim, latco), LwToaUpC(1:ibLim, latco), &
                            convts(1:ibLim, latco), convcs(1:ibLim, latco), convbs(1:ibLim, latco), &
                            convc(1:ibLim, latco), convt(1:ibLim, latco), convb(1:ibLim, latco), &
                            ustr(1:ibLim, latco), vstr(1:ibLim, latco), latco, &
                            yVisBeam(1:ibLim, latco), yVisDiff(1:ibLim, latco), yNirBeam(1:ibLim, latco), &
                            yNirDiff(1:ibLim, latco), ySwToaDown(1:ibLim, latco), yVisBeamC(1:ibLim, latco), &
                            yVisDiffC(1:ibLim, latco), yNirBeamC(1:ibLim, latco), yNirDiffC(1:ibLim, latco), &
                            cldsav(1:ibLim, latco), ssib(1:ibLim, latco), bstar(1:ibLim, latco), ibMaxPerJB(latco), kMax, &
                            sigki, lonrad(1:ibLim), ps(1:ibLim), &
                            var(1:ibLim, latco), sens(1:ibLim, latco), evap(1:ibLim, latco), &
                            cos2d(1:ibLim), intcosz, LwCoolRateC(1:ibLim, :, latco), topog(1:ibLim), &
                            o3mix(1:ibLim, :, latco), Mmlen(1:ibLim, latco), &
                            xland(1:ibLim, latco), lowlyr(1:ibLim, latco), &
                            ustar(1:ibLim, latco), z0(1:ibLim, latco), tkemyj(1:ibLim, :, latco), &
                            thz0(1:ibLim, latco), qz0 (1:ibLim, latco), &
                            uz0 (1:ibLim, latco), vz0 (1:ibLim, latco), &
                            pblh(1:ibLim, latco), akhs(1:ibLim, latco), akms(1:ibLim, latco), &
                            CT(1:ibLim, latco), snow(1:ibLim, latco), &
                            htdisp(1:ibLim, latco), temp2m(1:ibLim, latco), umes2m(1:ibLim, latco), &
                            cldtot(1:ibLim, :, latco), cldinv(1:ibLim, :, latco), cldsat(1:ibLim, :, latco), &
                            cldcon(1:ibLim, :, latco), cldson(1:ibLim, :, latco), clwd  (1:ibLim, :, latco), &
                            emisd (1:ibLim, :, latco), taud  (1:ibLim, :, latco), &
                            ySwSfcNet(1:ibLim, latco), ySwSfcNetC(1:ibLim, latco), &
                            rSwSfcNet(1:ibLim, latco), rSwSfcNetC(1:ibLim, latco), mskant(1:ibLim, latco), &
                            sm0(1:ibLim, :, latco), mlsi(1:ibLim, latco), & ! add solange 27-01-2012nClass
                            gndvi(1:ibLim, latco), ndvim(1:ibLim, latco), qliq(1:ibLim, :, latco), tpert(1:ibLim, latco), &
                            qpert(1:ibLim, latco), &
                            sflux_t(1:ibLim, latco), sflux_r(1:ibLim, latco), sflux_u(1:ibLim, latco), sflux_v(1:ibLim, latco), &
                            tstar(1:ibLim, latco), wstar(1:ibLim, latco), r_aer(1:ibLim, latco), veg_type(1:ibLim, latco, npatches), &
                            frac_occ(1:ibLim, latco, npatches), HML (1:ibLim, latco), HUML(1:ibLim, latco), &
                            HVML(1:ibLim, latco), TSK(1:ibLim, latco), z0sea(1:ibLim, latco), &
                            cu_hr (1:ibLim, :, latco), cu_kbot(1:ibLim, latco), cu_ktop(1:ibLim, latco), cu_Kuo(1:ibLim, latco), &
                            dudt     (1:ibLim, 1:kMax, latco), dvdt     (1:ibLim, 1:kMax, latco), &
                            EFFCS(1:ibLim, 1:kMax, latco), EFFIS (1:ibLim, 1:kMax, latco), tauresx(1:ibLim, latco), tauresy(1:ibLim, latco), &
                            gicem(1:ibLim, :), gicet(1:ibLim, :), gliqm(1:ibLim, :), gliqt(1:ibLim, :))
                ENDIF !for specSfc
            END IF
        ELSE
            IF (specSfc) THEN
                CALL physcs  (&
                        ! Coeficiente de transporte vertical para a turbulencia
                        PBL_CoefKm(1:ibLim, :, latco), PBL_CoefKh(1:ibLim, :, latco), &
                        gt(1:ibLim, :), gq(1:ibLim, :), gu(1:ibLim, :), gv(1:ibLim, :), &
                        gps(1:ibLim), &
                        ppli(1:ibLim, latco), &
                        ppci(1:ibLim, latco), gyu(1:ibLim, :), gyv(1:ibLim, :), gtd(1:ibLim, :), &
                        gqd(1:ibLim, :), ySwHeatRate(1:ibLim, :, latco), LwCoolRate(1:ibLim, :, latco), &
                        imask(1:ibLim, latco), rVisBeam(1:ibLim, latco), rVisDiff(1:ibLim, latco), &
                        rNirBeam(1:ibLim, latco), rNirDiff(1:ibLim, latco), LwSfcDown(1:ibLim, latco), &
                        gtsea(1:ibLim, latco), colrad(1:ibLim), sl, si(1:kMax + 1), del, ifday, &
                        tod, AlbVisBeam(1:ibLim, latco), AlbVisDiff(1:ibLim, latco), &
                        AlbNirBeam(1:ibLim, latco), AlbNirDiff(1:ibLim, latco), rSwToaDown(1:ibLim, latco), &
                        LwSfcNet(1:ibLim, latco), LwToaUp(1:ibLim, latco), gl0(1:ibLim, latco), &
                        zorl(1:ibLim, latco), gtt(1:ibLim, :), gqq(1:ibLim, :), &
                        sheleg(1:ibLim, latco), tseam(1:ibLim, latco), omg(1:ibLim, :), &
                        ySwHeatRateC(1:ibLim, :, latco), rVisBeamC(1:ibLim, latco), rVisDiffC(1:ibLim, latco), &
                        rNirBeamC(1:ibLim, latco), rNirDiffC(1:ibLim, latco), LwSfcDownC(1:ibLim, latco), &
                        LwSfcNetC(1:ibLim, latco), LwToaUpC(1:ibLim, latco), &
                        convts(1:ibLim, latco), convcs(1:ibLim, latco), convbs(1:ibLim, latco), &
                        convc(1:ibLim, latco), convt(1:ibLim, latco), convb(1:ibLim, latco), &
                        ustr(1:ibLim, latco), vstr(1:ibLim, latco), latco, &
                        yVisBeam(1:ibLim, latco), yVisDiff(1:ibLim, latco), yNirBeam(1:ibLim, latco), &
                        yNirDiff(1:ibLim, latco), ySwToaDown(1:ibLim, latco), yVisBeamC(1:ibLim, latco), &
                        yVisDiffC(1:ibLim, latco), yNirBeamC(1:ibLim, latco), yNirDiffC(1:ibLim, latco), &
                        cldsav(1:ibLim, latco), ssib(1:ibLim, latco), bstar(1:ibLim, latco), ibMaxPerJB(latco), kMax, &
                        sigki, lonrad(1:ibLim), ps(1:ibLim), &
                        var(1:ibLim, latco), sens(1:ibLim, latco), evap(1:ibLim, latco), &
                        cos2d(1:ibLim), intcosz, LwCoolRateC(1:ibLim, :, latco), topog(1:ibLim), &
                        o3mix(1:ibLim, :, latco), Mmlen(1:ibLim, latco), &
                        xland(1:ibLim, latco), lowlyr(1:ibLim, latco), &
                        ustar(1:ibLim, latco), z0(1:ibLim, latco), tkemyj(1:ibLim, :, latco), &
                        thz0(1:ibLim, latco), qz0 (1:ibLim, latco), &
                        uz0 (1:ibLim, latco), vz0 (1:ibLim, latco), &
                        pblh(1:ibLim, latco), akhs(1:ibLim, latco), akms(1:ibLim, latco), &
                        CT(1:ibLim, latco), snow(1:ibLim, latco), &
                        htdisp(1:ibLim, latco), temp2m(1:ibLim, latco), umes2m(1:ibLim, latco), &
                        cldtot(1:ibLim, :, latco), cldinv(1:ibLim, :, latco), cldsat(1:ibLim, :, latco), &
                        cldcon(1:ibLim, :, latco), cldson(1:ibLim, :, latco), clwd  (1:ibLim, :, latco), &
                        emisd (1:ibLim, :, latco), taud  (1:ibLim, :, latco), &
                        ySwSfcNet(1:ibLim, latco), ySwSfcNetC(1:ibLim, latco), &
                        rSwSfcNet(1:ibLim, latco), rSwSfcNetC(1:ibLim, latco), mskant(1:ibLim, latco), &
                        sm0(1:ibLim, :, latco), mlsi(1:ibLim, latco), & ! add solange 27-01-2012
                        gndvi(1:ibLim, latco), ndvim(1:ibLim, latco), qliq(1:ibLim, :, latco), tpert(1:ibLim, latco), &
                        qpert(1:ibLim, latco), &
                        sflux_t(1:ibLim, latco), sflux_r(1:ibLim, latco), sflux_u(1:ibLim, latco), sflux_v(1:ibLim, latco), &
                        tstar(1:ibLim, latco), wstar(1:ibLim, latco), r_aer(1:ibLim, latco), veg_type(1:ibLim, latco, npatches), &
                        frac_occ(1:ibLim, latco, npatches), HML (1:ibLim, latco), HUML(1:ibLim, latco), &
                        HVML(1:ibLim, latco), TSK(1:ibLim, latco), z0sea(1:ibLim, latco), &
                        cu_hr (1:ibLim, :, latco), cu_kbot(1:ibLim, latco), cu_ktop(1:ibLim, latco), cu_Kuo(1:ibLim, latco), &
                        dudt     (1:ibLim, 1:kMax, latco), dvdt     (1:ibLim, 1:kMax, latco), &
                        EFFCS(1:ibLim, 1:kMax, latco), EFFIS (1:ibLim, 1:kMax, latco), tauresx(1:ibLim, latco), tauresy(1:ibLim, latco), &
                        sgTsfc(1:ibLim), sgH(1:ibLim), sgLE(1:ibLim), sgTau(1:ibLim))
            ELSE  !for specSfc
                CALL physcs  (&
                        ! Coeficiente de transporte vertical para a turbulencia
                        PBL_CoefKm(1:ibLim, :, latco), PBL_CoefKh(1:ibLim, :, latco), &
                        gt(1:ibLim, :), gq(1:ibLim, :), gu(1:ibLim, :), gv(1:ibLim, :), &
                        gps(1:ibLim), &
                        ppli(1:ibLim, latco), &
                        ppci(1:ibLim, latco), gyu(1:ibLim, :), gyv(1:ibLim, :), gtd(1:ibLim, :), &
                        gqd(1:ibLim, :), ySwHeatRate(1:ibLim, :, latco), LwCoolRate(1:ibLim, :, latco), &
                        imask(1:ibLim, latco), rVisBeam(1:ibLim, latco), rVisDiff(1:ibLim, latco), &
                        rNirBeam(1:ibLim, latco), rNirDiff(1:ibLim, latco), LwSfcDown(1:ibLim, latco), &
                        gtsea(1:ibLim, latco), colrad(1:ibLim), sl, si(1:kMax + 1), del, ifday, &
                        tod, AlbVisBeam(1:ibLim, latco), AlbVisDiff(1:ibLim, latco), &
                        AlbNirBeam(1:ibLim, latco), AlbNirDiff(1:ibLim, latco), rSwToaDown(1:ibLim, latco), &
                        LwSfcNet(1:ibLim, latco), LwToaUp(1:ibLim, latco), gl0(1:ibLim, latco), &
                        zorl(1:ibLim, latco), gtt(1:ibLim, :), gqq(1:ibLim, :), &
                        sheleg(1:ibLim, latco), tseam(1:ibLim, latco), omg(1:ibLim, :), &
                        ySwHeatRateC(1:ibLim, :, latco), rVisBeamC(1:ibLim, latco), rVisDiffC(1:ibLim, latco), &
                        rNirBeamC(1:ibLim, latco), rNirDiffC(1:ibLim, latco), LwSfcDownC(1:ibLim, latco), &
                        LwSfcNetC(1:ibLim, latco), LwToaUpC(1:ibLim, latco), &
                        convts(1:ibLim, latco), convcs(1:ibLim, latco), convbs(1:ibLim, latco), &
                        convc(1:ibLim, latco), convt(1:ibLim, latco), convb(1:ibLim, latco), &
                        ustr(1:ibLim, latco), vstr(1:ibLim, latco), latco, &
                        yVisBeam(1:ibLim, latco), yVisDiff(1:ibLim, latco), yNirBeam(1:ibLim, latco), &
                        yNirDiff(1:ibLim, latco), ySwToaDown(1:ibLim, latco), yVisBeamC(1:ibLim, latco), &
                        yVisDiffC(1:ibLim, latco), yNirBeamC(1:ibLim, latco), yNirDiffC(1:ibLim, latco), &
                        cldsav(1:ibLim, latco), ssib(1:ibLim, latco), bstar(1:ibLim, latco), ibMaxPerJB(latco), kMax, &
                        sigki, lonrad(1:ibLim), ps(1:ibLim), &
                        var(1:ibLim, latco), sens(1:ibLim, latco), evap(1:ibLim, latco), &
                        cos2d(1:ibLim), intcosz, LwCoolRateC(1:ibLim, :, latco), topog(1:ibLim), &
                        o3mix(1:ibLim, :, latco), Mmlen(1:ibLim, latco), &
                        xland(1:ibLim, latco), lowlyr(1:ibLim, latco), &
                        ustar(1:ibLim, latco), z0(1:ibLim, latco), tkemyj(1:ibLim, :, latco), &
                        thz0(1:ibLim, latco), qz0 (1:ibLim, latco), &
                        uz0 (1:ibLim, latco), vz0 (1:ibLim, latco), &
                        pblh(1:ibLim, latco), akhs(1:ibLim, latco), akms(1:ibLim, latco), &
                        CT(1:ibLim, latco), snow(1:ibLim, latco), &
                        htdisp(1:ibLim, latco), temp2m(1:ibLim, latco), umes2m(1:ibLim, latco), &
                        cldtot(1:ibLim, :, latco), cldinv(1:ibLim, :, latco), cldsat(1:ibLim, :, latco), &
                        cldcon(1:ibLim, :, latco), cldson(1:ibLim, :, latco), clwd  (1:ibLim, :, latco), &
                        emisd (1:ibLim, :, latco), taud  (1:ibLim, :, latco), &
                        ySwSfcNet(1:ibLim, latco), ySwSfcNetC(1:ibLim, latco), &
                        rSwSfcNet(1:ibLim, latco), rSwSfcNetC(1:ibLim, latco), mskant(1:ibLim, latco), &
                        sm0(1:ibLim, :, latco), mlsi(1:ibLim, latco), & ! add solange 27-01-2012
                        gndvi(1:ibLim, latco), ndvim(1:ibLim, latco), qliq(1:ibLim, :, latco), tpert(1:ibLim, latco), &
                        qpert(1:ibLim, latco), &
                        sflux_t(1:ibLim, latco), sflux_r(1:ibLim, latco), sflux_u(1:ibLim, latco), sflux_v(1:ibLim, latco), &
                        tstar(1:ibLim, latco), wstar(1:ibLim, latco), r_aer(1:ibLim, latco), veg_type(1:ibLim, latco, npatches), &
                        frac_occ(1:ibLim, latco, npatches), HML (1:ibLim, latco), HUML(1:ibLim, latco), &
                        HVML(1:ibLim, latco), TSK(1:ibLim, latco), z0sea(1:ibLim, latco), &
                        cu_hr (1:ibLim, :, latco), cu_kbot(1:ibLim, latco), cu_ktop(1:ibLim, latco), cu_Kuo(1:ibLim, latco), &
                        dudt     (1:ibLim, 1:kMax, latco), dvdt     (1:ibLim, 1:kMax, latco), &
                        EFFCS(1:ibLim, 1:kMax, latco), EFFIS (1:ibLim, 1:kMax, latco), tauresx(1:ibLim, latco), tauresy(1:ibLim, latco))
            ENDIF !for specSfc
        ENDIF

    END SUBROUTINE DryPhysics

    !          CALL HumidPhysics(1 , 1  ,kMax  , fgqmm(ib,1:kMax,jb), fgtmpp(ib,1:kMax,jb) , &
    !                           fgqp  (ib,1:kMax,jb), fgpsp   (ib,       jb), fgu   (ib,1:kMax,jb) , &
    !                           fgv   (ib,1:kMax,jb), omg     (ib,1:kMax,jb), fgtmpm(ib,1:kMax,jb) , &
    !                           fgtmp (ib,1:kMax,jb), fgq     (ib,1:kMax,jb), fgps  (ib       ,jb) , &
    !                           fgzs  (ib       ,jb), colrad2D(ib,       jb) )

    SUBROUTINE HumidPhysics(latco, ibMax, kMax, rqn, ftn, &
            fqn, fpn, gu, &
            gv, omg, gtmpm1, &
            gtmpm2, fgqm2, fgps2, &
            fgzs, colrad, lonrad, gicem, &
            gicep, gliqm, gliqp, gvarm, gvarp)

        ! fgqmm  -> fgqmm    time -> t-1
        ! fgtmp  -> fgtmp    time -> t+1
        ! fgq    -> fgq      time -> t+1
        ! fgps   -> fgps     time -> t+1
        ! fgumm  -> fgumm2   time -> t
        ! fgvmm  -> fgvmm2   time -> t
        ! omg    -> omg2     time -> t
        ! gtmpm  -> fgtmpmm  time -> t-1
        ! gtmpm2 -> fgtmpmm2 time -> t
        ! fgqm2  -> fgqmm2   time -> t
        ! fgps2  -> fgps     time -> t

        INTEGER, INTENT(IN) :: latco
        INTEGER, INTENT(IN) :: ibMax
        INTEGER, INTENT(IN) :: kMax
        REAL(KIND = r8), INTENT(INOUT) :: rqn  (ibMax, kMax)
        REAL(KIND = r8), INTENT(INOUT) :: ftn  (ibMax, kMax)
        REAL(KIND = r8), INTENT(INOUT) :: fqn  (ibMax, kMax)
        REAL(KIND = r8), INTENT(IN) :: fpn  (ibMax)
        REAL(KIND = r8), INTENT(IN) :: gu   (ibMax, kMax)
        REAL(KIND = r8), INTENT(IN) :: gv   (ibMax, kMax)
        REAL(KIND = r8), INTENT(IN) :: omg  (ibMax, kMax)

        !snf
        !
        REAL(KIND = r8), INTENT(IN) :: gtmpm1  (ibMax, kMax)
        REAL(KIND = r8), INTENT(IN) :: gtmpm2  (ibMax, kMax)
        REAL(KIND = r8), INTENT(IN) :: fgqm2   (ibMax, kMax)
        REAL(KIND = r8), INTENT(IN) :: fgps2   (ibMax)
        REAL(KIND = r8), INTENT(IN) :: fgzs    (ibMax)
        REAL(KIND = r8), INTENT(IN) :: colrad  (ibMax)
        REAL(KIND = r8), INTENT(IN) :: lonrad  (ibMax)
        REAL(KIND = r8), OPTIONAL, INTENT(INOUT) :: gicem (ibMax, kmax)
        REAL(KIND = r8), OPTIONAL, INTENT(INOUT) :: gicep (ibMax, kmax)
        REAL(KIND = r8), OPTIONAL, INTENT(INOUT) :: gliqm (ibMax, kmax)
        REAL(KIND = r8), OPTIONAL, INTENT(INOUT) :: gliqp (ibMax, kmax)
        REAL(KIND = r8), OPTIONAL, INTENT(INOUT) :: gvarm (ibMax, kmax, nClass + nAeros)
        REAL(KIND = r8), OPTIONAL, INTENT(INOUT) :: gvarp (ibMax, kmax, nClass + nAeros)

        INTEGER (KIND = i8) :: ibLim
        INTEGER (KIND = i8) :: i
        INTEGER (KIND = i8) :: k
        INTEGER :: ins
        INTEGER (KIND = i8) :: kprox
        REAL(KIND = r8) :: topog     (ibMax)
        REAL(KIND = r8) :: StpC_Temp (ibMax, kMax)
        REAL(KIND = r8) :: StpM_Temp (ibMax, kMax)
        REAL(KIND = r8) :: StpC_Umes (ibMax, kMax)
        REAL(KIND = r8) :: StpC_Uvel (ibMax, kMax)
        REAL(KIND = r8) :: StpC_Vvel (ibMax, kMax)

        REAL(KIND = r8) :: tracervar (ibMax, kMax, nClass + nAeros)
        REAL(KIND = r8) :: tracervarm(ibMax, kMax, nClass + nAeros)

        REAL(KIND = r8) :: tracerice (ibMax, kMax)
        REAL(KIND = r8) :: tracericem(ibMax, kMax)
        REAL(KIND = r8) :: tracerliq (ibMax, kMax)
        REAL(KIND = r8) :: tracerliqm (ibMax, kMax)

        REAL(KIND = r8) :: sst(ibMax)
        REAL(KIND = r8) :: fac
        REAL(KIND = r8) :: fac2
        REAL(KIND = r8) :: fac2x
        !snf
        !-------------------------------------------------
        !t-1     StpM_Temp
        !t       StpC_Temp ,StpC_Umes ,StpC_Uvel ,StpC_Vvel ,StpC_Pslc
        !t+1     StpP_Pslc
        !-----
        tracerice = 0.0_r8
        tracerliq = 0.0_r8
        tracericem = 0.0_r8
        tracerliqm = 0.0_r8
        tracervar = 0.0_r8
        tracervarm = 0.0_r8
        ibLim = ibMax!ibMaxPerJB(latco)
        topog(1:ibLim) = fgzs(1:ibLim) / (grav)
        fac = 0.5_r8
        IF(ifilt.EQ.0.AND.kt.EQ.0.AND.jdt.EQ.1) fac = 0.0_r8
        fac2 = 2.0_r8 * fac
        fac2x = 2.0_r8 * fac
        IF(ifilt.EQ.0.AND.kt.EQ.0.AND.jdt.EQ.2) fac2x = 2.0_r8
        DO i = 1, ibLim
            sst(i) = ABS(gtsea(i, latco))
        END DO
        DO k = 1, kMax
            DO i = 1, ibLim
                StpC_Uvel(i, k) = gu(i, k) / (SIN(colrad (i)))
                StpC_Vvel(i, k) = gv(i, k) / (SIN(colrad (i)))
                StpM_Temp(i, k) = gtmpm1(i, k) + tov(k)
                StpC_Temp(i, k) = gtmpm2(i, k) + tov(k)
                StpC_Umes(i, k) = fgqm2 (i, k)
            END DO
        END DO
        ! print*, 'gps=', fgps2, fpn; stop
        IF (microphys) THEN
            IF((nClass + nAeros)>0 .and. PRESENT(gvarm))THEN
                CALL cumulus_driver(&
                        ! Run Flags
                        ! Time info
                        fac, &          !REAL(KIND=r8), INTENT(IN   ) :: fac
                        fac2, &          !REAL(KIND=r8), INTENT(IN   ) :: fac2
                        fac2x, &          !REAL(KIND=r8), INTENT(IN   ) :: fac2x
                        ! Model Geometry
                        colrad   (1:ibLim), &          !REAL(KIND=r8), INTENT(IN   ) :: colrad(iMax)
                        lonrad   (1:ibLim), &          !REAL(KIND=r8), INTENT(IN   ) :: lonrad(iMax)
                        del      (1:kMax), &          !REAL(KIND=r8), INTENT(IN   ) :: del  (kMax)
                        sl       (1:kMax), &          !REAL(KIND=r8), INTENT(IN   ) :: sl   (kMax)
                        si       (1:kMax + 1), &          !REAL(KIND=r8), INTENT(IN   ) :: si   (kMax+1)
                        ! Model information
                        ibMaxPerJB(latco), &          !INTEGER      , INTENT(IN   ) :: iMax
                        kMax, &          !INTEGER      , INTENT(IN   ) :: kMax
                        latco, &          !INTEGER      , INTENT(IN   ) :: latco
                        imask    (1:ibLim, latco), &          !INTEGER(KIND=i8),INTENT(IN ) :: mask (iMax)
                        topog    (1:ibLim), &          !REAL(KIND=r8), INTENT(IN   ) :: zs   (iMax)
                        ! CONVECTION: convective clouds
                        convc    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: convc  (iMax)
                        convt    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: convt  (iMax)
                        convb    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: convb  (iMax)
                        toplv    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: toplv  (iMax)
                        botlv    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: botlv  (iMax)
                        convts   (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: convts (iMax)
                        convcs   (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: convcs (iMax)
                        convbs   (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: convbs (iMax)
                        concld   (1:ibLim, 1:kMax, latco), &    !REAL(KIND=r8), INTENT(INOUT) :: concld (iMax,kMax)
                        cld      (1:ibLim, 1:kMax, latco), &    !REAL(KIND=r8), INTENT(INOUT) :: cld  (iMax,kMax)
                        ! SURFACE:  Fields
                        sst      (1:ibLim), &          !REAL(KIND=r8), INTENT(IN   ) :: tsfc (iMax)
                        tpert    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(IN   ) :: tpert(iMax)
                        qpert    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(IN   ) :: qpert(iMax)
                        sens     (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(IN   ) :: sens (iMax)
                        evap     (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(IN   ) :: evap (iMax)
                        ustar    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(IN   ) :: ustar (iMax)
                        ! Precipitation Field
                        prcp1    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: prcp1  (iMax)
                        prcp2    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: prcp2  (iMax)
                        prcp3    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: prcp3  (iMax)
                        prcpt    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: prcpt  (iMax)
                        geshem   (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: geshem (iMax)
                        ppli     (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: ppli  (iMax)
                        ppci     (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: ppci  (iMax)
                        prct     (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: prct  (iMax)
                        prcc     (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: prcc  (iMax)
                        snow     (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: snowfl (iMax)
                        ! PBL:  Fields
                        PBL_CoefKh(1:ibLim, 1:kMax, latco), &   !REAL(KIND=r8), INTENT(IN   ) :: PBL_CoefKh(iMax,kMax)
                        tkemyj    (1:ibLim, 1:kMax, latco), &   !REAL(KIND=r8), INTENT(IN   ) :: tke  (iMax,kMax)
                        ! Microphysics
                        dudt     (1:ibLim, 1:kMax, latco), &
                        dvdt     (1:ibLim, 1:kMax, latco), &
                        qliq     (1:ibLim, 1:kMax, latco), &   !REAL(KIND=r8), INTENT(OUT  ) :: qliq(iMax,kMax)
                        EFFCS    (1:ibLim, 1:kMax, latco), &   !REAL(KIND=r8), INTENT(INOUT) :: EFFCS   (iMax,kMax)
                        EFFIS    (1:ibLim, 1:kMax, latco), &   !REAL(KIND=r8), INTENT(INOUT) :: EFFIS   (iMax,kMax)
                        ! Atmospheric fields
                        StpM_Temp(1:ibLim, :), &          !REAL(KIND=r8), INTENT(INOUT) :: ta (iMax,kMax)
                        StpC_Temp(1:ibLim, :), &          !REAL(KIND=r8), INTENT(INOUT) :: tb (iMax,kMax)
                        ftn      (1:ibLim, :), &          !REAL(KIND=r8), INTENT(INOUT) :: tc (iMax,kMax)
                        rqn      (1:ibLim, :), &          !REAL(KIND=r8), INTENT(INOUT) :: qa (iMax,kMax)
                        StpC_Umes(1:ibLim, :), &          !REAL(KIND=r8), INTENT(INOUT) :: qb (iMax,kMax)
                        fqn      (1:ibLim, :), &          !REAL(KIND=r8), INTENT(INOUT) :: qc (iMax,kMax)
                        StpC_Uvel(1:ibLim, :), &          !REAL(KIND=r8), INTENT(IN   ) :: ub   (iMax,kMax) ! (m/s)
                        StpC_Vvel(1:ibLim, :), &          !REAL(KIND=r8), INTENT(IN   ) :: vb   (iMax,kMax) ! (m/s)
                        omg      (1:ibLim, :), &          !REAL(KIND=r8), INTENT(IN   ) :: omgb (iMax,kMax) ! (Pa/s)
                        fpn      (1:ibLim), &          !REAL(KIND=r8), INTENT(IN   ) :: psb  (iMax)
                        fgps2    (1:ibLim), &          !REAL(KIND=r8), INTENT(IN   ) :: psb2 (iMax)
                        gicep(1:ibLim, 1:kMax), &          !REAL(KIND=r8), INTENT(INOUT) :: gicep  (iMax,kmax)
                        gicem(1:ibLim, 1:kMax), &          !REAL(KIND=r8), INTENT(INOUT) :: gicem  (iMax,kmax)
                        gliqp(1:ibLim, 1:kMax), &          !REAL(KIND=r8), INTENT(INOUT) :: gliqp  (iMax,kmax)
                        gliqm(1:ibLim, 1:kMax), &          !REAL(KIND=r8), INTENT(INOUT) :: gliqm  (iMax,kmax)
                        gvarp(1:ibLim, 1:kMax, 1:nClass + nAeros), &    !REAL(KIND=r8),OPTIONAL,   INTENT(INOUT) :: gvarp (iMax,kmax,nClass+nAeros)
                        gvarm(1:ibLim, 1:kMax, 1:nClass + nAeros))      !REAL(KIND=r8),OPTIONAL,   INTENT(INOUT) :: gvarm (iMax,kmax,nClass+nAeros)

            ELSE
                CALL cumulus_driver(&
                        ! Run Flags
                        ! Time info
                        fac, &          !REAL(KIND=r8), INTENT(IN   ) :: fac
                        fac2, &          !REAL(KIND=r8), INTENT(IN   ) :: fac2
                        fac2x, &          !REAL(KIND=r8), INTENT(IN   ) :: fac2x
                        ! Model Geometry
                        colrad   (1:ibLim), &          !REAL(KIND=r8), INTENT(IN   ) :: colrad(iMax)
                        lonrad   (1:ibLim), &          !REAL(KIND=r8), INTENT(IN   ) :: lonrad(iMax)
                        del      (1:kMax), &          !REAL(KIND=r8), INTENT(IN   ) :: del  (kMax)
                        sl       (1:kMax), &          !REAL(KIND=r8), INTENT(IN   ) :: sl   (kMax)
                        si       (1:kMax + 1), &          !REAL(KIND=r8), INTENT(IN   ) :: si   (kMax+1)
                        ! Model information
                        ibMaxPerJB(latco), &          !INTEGER      , INTENT(IN   ) :: iMax
                        kMax, &          !INTEGER      , INTENT(IN   ) :: kMax
                        latco, &          !INTEGER      , INTENT(IN   ) :: latco
                        imask    (1:ibLim, latco), &          !INTEGER(KIND=i8),INTENT(IN ) :: mask (iMax)
                        topog    (1:ibLim), &          !REAL(KIND=r8), INTENT(IN   ) :: zs   (iMax)
                        ! CONVECTION: convective clouds
                        convc    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: convc  (iMax)
                        convt    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: convt  (iMax)
                        convb    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: convb  (iMax)
                        toplv    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: toplv  (iMax)
                        botlv    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: botlv  (iMax)
                        convts   (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: convts (iMax)
                        convcs   (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: convcs (iMax)
                        convbs   (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: convbs (iMax)
                        concld   (1:ibLim, 1:kMax, latco), &    !REAL(KIND=r8), INTENT(INOUT) :: concld (iMax,kMax)
                        cld      (1:ibLim, 1:kMax, latco), &    !REAL(KIND=r8), INTENT(INOUT) :: cld   (iMax,kMax)
                        ! SURFACE:  Fields
                        sst      (1:ibLim), &          !REAL(KIND=r8), INTENT(IN   ) :: tsfc (iMax)
                        tpert    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(IN   ) :: tpert(iMax)
                        qpert    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(IN   ) :: qpert(iMax)
                        sens     (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(IN   ) :: sens (iMax)
                        evap     (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(IN   ) :: evap (iMax)
                        ustar    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(IN   ) :: ustar (iMax)
                        ! Precipitation Field
                        prcp1    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: prcp1  (iMax)
                        prcp2    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: prcp2  (iMax)
                        prcp3    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: prcp3  (iMax)
                        prcpt    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: prcpt  (iMax)
                        geshem   (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: geshem (iMax)
                        ppli     (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: ppli   (iMax)
                        ppci     (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: ppci   (iMax)
                        prct     (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: prct   (iMax)
                        prcc     (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: prcc   (iMax)
                        snow     (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: snowfl (iMax)
                        ! PBL:  Fields
                        PBL_CoefKh(1:ibLim, 1:kMax, latco), &   !REAL(KIND=r8), INTENT(IN   ) :: PBL_CoefKh(iMax,kMax)
                        tkemyj    (1:ibLim, 1:kMax, latco), &   !REAL(KIND=r8), INTENT(IN   ) :: tke  (iMax,kMax)
                        ! Microphysics
                        dudt     (1:ibLim, 1:kMax, latco), &
                        dvdt     (1:ibLim, 1:kMax, latco), &
                        qliq     (1:ibLim, 1:kMax, latco), &   !REAL(KIND=r8), INTENT(OUT  ) :: qliq(iMax,kMax)
                        EFFCS    (1:ibLim, 1:kMax, latco), &   !REAL(KIND=r8), INTENT(INOUT) :: EFFCS   (iMax,kMax)
                        EFFIS    (1:ibLim, 1:kMax, latco), &   !REAL(KIND=r8), INTENT(INOUT) :: EFFIS   (iMax,kMax)
                        ! Atmospheric fields
                        StpM_Temp(1:ibLim, :), &          !REAL(KIND=r8), INTENT(INOUT) :: ta (iMax,kMax)
                        StpC_Temp(1:ibLim, :), &          !REAL(KIND=r8), INTENT(INOUT) :: tb (iMax,kMax)
                        ftn      (1:ibLim, :), &          !REAL(KIND=r8), INTENT(INOUT) :: tc (iMax,kMax)
                        rqn      (1:ibLim, :), &          !REAL(KIND=r8), INTENT(INOUT) :: qa (iMax,kMax)
                        StpC_Umes(1:ibLim, :), &          !REAL(KIND=r8), INTENT(INOUT) :: qb (iMax,kMax)
                        fqn      (1:ibLim, :), &          !REAL(KIND=r8), INTENT(INOUT) :: qc (iMax,kMax)
                        StpC_Uvel(1:ibLim, :), &          !REAL(KIND=r8), INTENT(IN   ) :: ub   (iMax,kMax) ! (m/s)
                        StpC_Vvel(1:ibLim, :), &          !REAL(KIND=r8), INTENT(IN   ) :: vb   (iMax,kMax) ! (m/s)
                        omg      (1:ibLim, :), &          !REAL(KIND=r8), INTENT(IN   ) :: omgb (iMax,kMax) ! (Pa/s)
                        fpn      (1:ibLim), &          !REAL(KIND=r8), INTENT(IN   ) :: psb  (iMax)
                        fgps2    (1:ibLim), &          !REAL(KIND=r8), INTENT(IN   ) :: psb2 (iMax)
                        gicep(1:ibLim, 1:kMax), &          !REAL(KIND=r8), INTENT(INOUT) :: gicep  (iMax,kmax)
                        gicem(1:ibLim, 1:kMax), &          !REAL(KIND=r8), INTENT(INOUT) :: gicem  (iMax,kmax)
                        gliqp(1:ibLim, 1:kMax), &          !REAL(KIND=r8), INTENT(INOUT) :: gliqp  (iMax,kmax)
                        gliqm(1:ibLim, 1:kMax))          !REAL(KIND=r8), INTENT(INOUT) :: gliqm  (iMax,kmax)
            END IF
        ELSE

            CALL cumulus_driver(&
                    ! Run Flags
                    ! Time info
                    fac, &          !REAL(KIND=r8), INTENT(IN   ) :: fac
                    fac2, &          !REAL(KIND=r8), INTENT(IN   ) :: fac2
                    fac2x, &          !REAL(KIND=r8), INTENT(IN   ) :: fac2x
                    ! Model Geometry
                    colrad   (1:ibLim), &          !REAL(KIND=r8), INTENT(IN   ) :: colrad(iMax)
                    lonrad   (1:ibLim), &          !REAL(KIND=r8), INTENT(IN   ) :: lonrad(iMax)
                    del      (1:kMax), &          !REAL(KIND=r8), INTENT(IN   ) :: del  (kMax)
                    sl       (1:kMax), &          !REAL(KIND=r8), INTENT(IN   ) :: sl   (kMax)
                    si       (1:kMax + 1), &          !REAL(KIND=r8), INTENT(IN   ) :: si   (kMax+1)
                    ! Model information
                    ibMaxPerJB(latco), &          !INTEGER   , INTENT(IN   ) :: iMax
                    kMax, &          !INTEGER   , INTENT(IN   ) :: kMax
                    latco, &          !INTEGER   , INTENT(IN   ) :: latco
                    imask    (1:ibLim, latco), &          !INTEGER(KIND=i8),INTENT(IN ) :: mask (iMax)
                    topog    (1:ibLim), &          !REAL(KIND=r8), INTENT(IN   ) :: zs   (iMax)
                    ! CONVECTION: convective clouds
                    convc    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: convc  (iMax)
                    convt    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: convt  (iMax)
                    convb    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: convb  (iMax)
                    toplv    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: toplv  (iMax)
                    botlv    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: botlv  (iMax)
                    convts   (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: convts (iMax)
                    convcs   (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: convcs (iMax)
                    convbs   (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: convbs (iMax)
                    concld   (1:ibLim, 1:kMax, latco), &    !REAL(KIND=r8), INTENT(INOUT) :: concld (iMax,kMax)
                    cld      (1:ibLim, 1:kMax, latco), &    !REAL(KIND=r8), INTENT(INOUT) :: cld    (iMax,kMax)
                    ! SURFACE:  Fields
                    sst      (1:ibLim), &          !REAL(KIND=r8), INTENT(IN   ) :: tsfc (iMax)
                    tpert    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(IN   ) :: tpert(iMax)
                    qpert    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(IN   ) :: qpert(iMax)
                    sens     (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(IN   ) :: sens (iMax)
                    evap     (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(IN   ) :: evap (iMax)
                    ustar    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(IN   ) :: ustar (iMax)
                    ! Precipitation Field
                    prcp1    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: prcp1  (iMax)
                    prcp2    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: prcp2  (iMax)
                    prcp3    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: prcp3  (iMax)
                    prcpt    (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: prcpt  (iMax)
                    geshem   (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: geshem (iMax)
                    ppli     (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: ppli   (iMax)
                    ppci     (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: ppci   (iMax)
                    prct     (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: prct   (iMax)
                    prcc     (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: prcc   (iMax)
                    snow     (1:ibLim, latco), &          !REAL(KIND=r8), INTENT(INOUT) :: snowfl (iMax)
                    ! PBL:  Fields
                    PBL_CoefKh(1:ibLim, 1:kMax, latco), &   !REAL(KIND=r8), INTENT(IN   ) :: PBL_CoefKh(iMax,kMax)
                    tkemyj    (1:ibLim, 1:kMax, latco), &   !REAL(KIND=r8), INTENT(IN   ) :: tke  (iMax,kMax)
                    ! Microphysics
                    dudt     (1:ibLim, 1:kMax, latco), &
                    dvdt     (1:ibLim, 1:kMax, latco), &
                    qliq     (1:ibLim, 1:kMax, latco), &   !REAL(KIND=r8), INTENT(OUT  ) :: qliq(iMax,kMax)
                    EFFCS    (1:ibLim, 1:kMax, latco), &   !REAL(KIND=r8), INTENT(INOUT) :: EFFCS   (iMax,kMax)
                    EFFIS    (1:ibLim, 1:kMax, latco), &   !REAL(KIND=r8), INTENT(INOUT) :: EFFIS   (iMax,kMax)
                    ! Atmospheric fields
                    StpM_Temp(1:ibLim, :), &          !REAL(KIND=r8), INTENT(INOUT) :: ta (iMax,kMax)
                    StpC_Temp(1:ibLim, :), &          !REAL(KIND=r8), INTENT(INOUT) :: tb (iMax,kMax)
                    ftn      (1:ibLim, :), &          !REAL(KIND=r8), INTENT(INOUT) :: tc (iMax,kMax)
                    rqn      (1:ibLim, :), &          !REAL(KIND=r8), INTENT(INOUT) :: qa (iMax,kMax)
                    StpC_Umes(1:ibLim, :), &          !REAL(KIND=r8), INTENT(INOUT) :: qb (iMax,kMax)
                    fqn      (1:ibLim, :), &          !REAL(KIND=r8), INTENT(INOUT) :: qc (iMax,kMax)
                    StpC_Uvel(1:ibLim, :), &          !REAL(KIND=r8), INTENT(IN   ) :: ub   (iMax,kMax) ! (m/s)
                    StpC_Vvel(1:ibLim, :), &          !REAL(KIND=r8), INTENT(IN   ) :: vb   (iMax,kMax) ! (m/s)
                    omg      (1:ibLim, :), &          !REAL(KIND=r8), INTENT(IN   ) :: omgb (iMax,kMax) ! (Pa/s)
                    fpn      (1:ibLim), &          !REAL(KIND=r8), INTENT(IN   ) :: psb  (iMax)
                    fgps2    (1:ibLim), &          !REAL(KIND=r8), INTENT(IN   ) :: psb2 (iMax)
                    tracerice(1:ibLim, 1:kMax), &          !REAL(KIND=r8), INTENT(INOUT) :: gicep  (iMax,kmax)
                    tracericem(1:ibLim, 1:kMax), &          !REAL(KIND=r8), INTENT(INOUT) :: gicem  (iMax,kmax)
                    tracerliq(1:ibLim, 1:kMax), &          !REAL(KIND=r8), INTENT(INOUT) :: gliqp  (iMax,kmax)
                    tracerliqm(1:ibLim, 1:kMax))           !REAL(KIND=r8), INTENT(INOUT) :: gliqm  (iMax,kmax)
        END IF

    END SUBROUTINE HumidPhysics


    SUBROUTINE physcs(&
            ! Coeficiente de transporte vertical para a turbulencia
            PBL_CoefKm, PBL_CoefKh, gt, gq, gu, gv, gps, &
            ppli, ppci, gyu, gyv, gtd, gqd, ySwHeatRate, &
            LwCoolRate, imask, rVisBeam, rVisDiff, rNirBeam, rNirDiff, LwSfcDown, &
            tsea, colrad, sig, sigml, delsig, ifday, tod, &
            AlbVisBeam, AlbVisDiff, AlbNirBeam, AlbNirDiff, rSwToaDown, LwSfcNet, LwToaUp, &
            gl0, zorl, gtt, gqq, sheleg, tseam, omg, &
            ySwHeatRateC, rVisBeamC, rVisDiffC, rNirBeamC, rNirDiffC, LwSfcDownC, LwSfcNetC, &
            LwToaUpC, convts, convcs, convbs, convc, convt, convb, &
            ustr, vstr, latco, yVisBeam, yVisDiff, yNirBeam, yNirDiff, &
            ySwToaDown, yVisBeamC, yVisDiffC, yNirBeamC, yNirDiffC, cldsav, ssib, &
            bstar, ncols, kmax, sigki, lonrad, ps, var, &
            sens, evap, cos2d, intcosz, LwCoolRateC, topog, o3mix, &
            Mmlen, xland, lowlyr, ustar, z0, tkemyj, &
            thz0, qz0, uz0, vz0, pblh, akhs, akms, &
            ct, snow, htdisp, temp2m, umes2m, cldtot, cldinv, cldsat, cldcon, cldson, clwd, emisd, taud, &
            ySwSfcNet, ySwSfcNetC, rSwSfcNet, rSwSfcNetC, MskAnt, &
            sm0, mlsi, &  ! add solange 27-01-2012
            ndvi, ndvim, qliq, tpert, qpert, sflux_t, sflux_r, &
            sflux_u, sflux_v, tstar, wstar, r_aer, veg_type, frac_occ, HML, HUML, HVML, TSK, z0sea, &
            cu_hr, cu_kbot, cu_ktop, cu_Kuo, &
            dudt, dvdt, &
            EFFCS, EFFIS, tauresx, tauresy, &
            sgTsfc, sgH, sgLE, sgTau, &      ! variables for specified surface (specSfc)
            gicem, gicet, gliqm, gliqt, gvarm, gvart)

        !
        !
        ! physcs :main subroutine for turbulence closure
        !         hashvadahn radiation coupled 3-d model
        !         p.sellers  sib
        !         gps is in mb
        !==========================================================================
        ! ncols.....Number of grid points on a gaussian latitude circle
        ! jmax......Number of gaussian latitudes
        ! kmax......Number of sigma levels
        ! nls..... .Number of layers in the stratosphere.
        ! nlcs......nlcs =   30
        ! maxtid..../include/T062L28/restim.inc:
        !           constant integer, parameter maxtid=131760
        ! ityp......Numero das classes de solo 13
        ! imon......Max. number of month at year (12)
        ! icg.......Parameter of the vegetation  (icg=1 top e icg=2 bottom )
        ! iwv.......Compriment de onda iwv=1=visivel, iwv=2=infravermelho
        !           proximo, iwv=3 infravermelho termal
        ! idp.......Parameter to the layers of soils idp=1->3
        ! ibd.......Condiction of vegetation ibd=1 green / ibd=2
        ! gt........Temperature
        ! gq........Specific humidity
        ! gu........(zonal      velocity)*sin(colat)
        ! gv........(meridional velocity)*sin(colat)
        ! gps.......Surface pressure in mb
        ! tc0.......Temperatura da copa "dossel"(K)   modificada
        ! tg0.......Temperatura da superficie do solo  (K)   modificada
        ! td0.......Temperatura do solo profundo (K)   modificada
        ! w0(id)....Grau de saturacao de umidade do solo id=1 na camada superficial
        ! w0(id)....Grau de saturacao de umidade do solo id=2 na camada de raizes
        ! w0(id)....Grau de saturacao de umidade do solo id=3 na camada de drenagem
        ! sm0(id)...Conteudo de umidade do solo id=1 na camada superficial  (m3/m3)
        ! sm0(id)...Conteudo de umidade do solo id=2 na camada de raizes    (m3/m3)
        ! sm0(id)...Conteudo de umidade do solo id=3 na camada de drenagem  (m3/m3)
        ! capac0(iv).Agua interceptada iv=1 no dossel "water store capacity
        !             of leaves"(m)  modificada
        ! capac0(iv).Agua interceptada iv=2 na cobertura do solo (m)   modificada
        ! tcm........Temperatura da copa "dossel"(K)
        ! tgm........Temperatura da superficie do solo  (K)
        ! tdm........Temperatura do solo profundo (K)
        ! wm
        ! capacm.....Agua interceptada iv=2 na cobertura do solo (m)
        ! ppli.......Precipitation rate ( large scale )       (mm/s)
        ! ppci.......Precipitation rate ( cumulus )           (mm/s)
        ! gyu........-(dv/dt) negative of tendency of v*cos(latitude)
        ! gyv........(du/dt) tendency of zonal wind * cos(latitude)
        ! gtd
        ! gqd........Specific humidity
        ! ySwHeatRate........Heating rate due to shrt wave radiation in deg/sec
        ! LwCoolRate........Cooling rate due to long wave radiation in deg/sec
        ! imask......mascara continetal
        ! rVisBeam......visible beam cloudy skies (refer to downward surface fluxes)
        ! rVisDiff......visible diffuse cloudy skies (refer to downward surface fluxes)
        ! rNirBeam......near-ir beam cloudy skies (refer to downward surface fluxes)
        ! rNirDiff......near-ir diffuse cloudy skies (refer to downward surface fluxes)
        ! LwSfcDown.....downward longwave radiation at the bottom in w/m**2
        ! tsea.......effective surface radiative temperature ( tgeff )
        ! colrad.....colatitude  colrad=0-3.14 from np to sp in radians
        ! sig........sigma coordinate at middle of layer
        ! sigml......sigma coordinate at bottom of layer
        ! delsig      k=2  ****gu,gv,gt,gq,gyu,gyv,gtd,gqd,sig*** } delsig(2)
        !             k=3/2----sigml,ric,rf,km,kh,b,l -----------
        !             k=1  ****gu,gv,gt,gq,gyu,gyv,gtd,gqd,sig*** } delsig(1)
        !             k=1/2----sigml ----------------------------
        !
        ! istrt.......istrt = jdt =time step in getdia
        ! ifday.......model forecast day
        ! tod.........model forecast time of day in seconds
        ! AlbVisBeam.......visible beam surface albedo
        ! AlbVisDiff.......visible diffuse surface albedo
        ! AlbNirBeam.......near-ir beam surface albedo
        ! AlbNirDiff.......near-ir diffuse surface albedo
        ! uswtop......shortwave upward at top
        ! rSwToaDown......swinc....solar input at top of atmosphere
        ! LwSfcNet..........net surface ir radiation in w/m**2
        ! LwToaUp......long wave flux at top of atmosphere in w/m**2
        ! alon........define constant alon=0.0 at subroutine gfidi.f90:
        ! dt........time interval,usually =delt,but changes
        !             in nlnmi (dt=1.) and at dead start(delt/4,delt/2)
        ! intg........intg =2  time integration of surface physical variable
        !                      is done by leap-frog implicit scheme. this
        !                      conseves enegy and h2o.
        !             intg =1  time integration of surface physical variable
        !                      is done by backward implicit scheme.
        ! gl0.........maximum mixing length l0 in blackerdar's formula
        !             l=k0*z/(1+k0*z/l0)
        ! zorl........zorl (i)= 100.0 *zgrav*speedm(i)*rhi(i)
        !             zgrav =0.032 /grav and rhi   (i)=ct(i)*ustar(i), where
        !             ct is heat transfer coefficients.
        !             ustar is surface friction velocity  (m/s)
        !             speedm(i)=SQRT(gu(i)**2+gv(i)**2)*sincli, where
        !             sincli=1.0 /sinclt
        !
        ! gtt.........gtt =  gtmp(imx,kmax) input  : temperature (fourier).
        !                                   output : "s" as given by eq. (19)
        !                                            in noaa tech report nws 30
        ! gqq.........gqq = gq(imx,kmax)     input : specific humidity (fourier).
        !                                   output : tendency of specific humidity
        !                                            without sources and sinks
        !                                            (fourier).
        ! sheleg......snow amount in mm (equivalent water depth)
        ! tseam.......tseam = gtsea (ncols,jmax)  input(gsstcd) lt.0=sea surface temp.
        !                                                      gt.0=ground temp.
        ! omg.........omg   =  vertical velocity  (cb/sec)
        ! rVisBeamC......Visible beam clear sky (Refer to downward surface
        !             shortwave fluxe)
        ! rVisDiffC......Visible diffuse clear sky (Refer to downward surface
        !             shortwave fluxe)
        ! rNirBeamC......Near-IR beam clear skies (Refer to downward surface
        !             shortwave fluxe)
        ! rNirDiffC......Near-IR diffuse clear skies (Refer to downward surface
        !             shortwave fluxe)
        ! ySwHeatRateC........Heating rate due to shortwave (clear) (K/s)
        ! LwCoolRateC........Cooling rate due to longwave (clear) (K/s)   !hmjb
        ! LwSfcDownC......Downward longwave (clear) At the bottom (W/m2)
        ! LwSfcNetC.......net longwave at bottom (clear)
        ! LwToaUpC......longwave upward at top (clear)
        ! pie.........Constant pi=3.1415926e0
        ! stefan......Stefan Stefan Boltzman constant
        ! cpair.......Specific heat of air           (j/kg/k)
        ! hl..........heat of evaporation of water     (j/kg)
        ! grav........grav   gravity constant        (m/s**2)
        ! snomel......Calor latente de fusao is expressed in (j m-1)
        ! tf..........Temperatura de congelamento (K)=273.16e0
        ! clai........heat capacity of foliage
        ! cw..........liquid water heat capacity     (j/m**3)
        ! gasr........gas constant of dry air        (j/kg/k)
        ! epsfac......Constante 0.622 Razao entre as massas
        !             moleculares do vapor e do ar seco
        ! athird......Constant athird =                 1.0e0/3.0e0
        ! tice........tice=271.16 ice temperature ice
        ! oceald......oceald = 0.0419e0
        ! z0ice ......Roughness length of ice
        ! icealn......near-ir beam surface albedo or near-ir diffuse surface albedo
        ! icealv......visible beam surface albedo or visible diffuse surface albedo
        ! dtc3x.......time increment dt
        ! nmax........Number of point grid on continent
        ! nsx.........Phenology dates to fall within one year period
        ! itype.......Classe de textura do solo
        ! vcover(iv)..Fracao de cobertura de vegetacao iv=1 Top
        ! vcover(iv)..Fracao de cobertura de vegetacao iv=2 Bottom
        ! z0x.........Roughness length
        ! d...........Displacement height
        ! rdc.........constant related to aerodynamic resistance
        ! rbc.........Constant related to bulk boundary layer resistance
        ! z0..........Roughness length
        ! qm..........Reference specific humidity (fourier)
        ! tm..........Reference temperature    (fourier)                (k)
        ! um..........Razao entre zonal pseudo-wind (fourier) e seno da
        !             colatitude
        ! vm..........Razao entre meridional pseudo-wind (fourier) e seno da
        !             colatitude
        ! psur........Surface pressure in mb
        ! ppc.........Precipitation rate ( cumulus )           (mm/s)
        ! ppl.........Precipitation rate ( large scale )       (mm/s)
        ! radn........Downward sw/lw radiation at the surface
        ! tc..........Temperatura da copa "dossel"(K)
        ! tg .........Temperatura da superficie do solo (K)
        ! td..........Temperatura do solo profundo (K)
        ! capac(iv)...Agua interceptada iv=1 no dossel "water store capacity
        !             of leaves"(m)
        ! capac(iv)...Agua interceptada iv=2 na cobertura do solo (m)
        ! w(id).......Grau de saturacao de umidade do solo id=1 na camada superficial
        ! w(id).......Grau de saturacao de umidade do solo id=2 na camada de raizes
        ! w(id).......Grau de saturacao de umidade do solo id=3 na camada de drenagem
        ! ra..........Resistencia Aerodinamica (s/m)
        ! rb..........bulk boundary layer resistance
        ! rd..........Aerodynamic resistance between ground      (s/m)
        !             and canopy air space
        ! rc..........Resistencia do topo da copa
        ! rg..........Resistencia da base da copa
        ! tcta........Diferenca entre tc-ta                      (k)
        ! tgta........Diferenca entre tg-ta                      (k)
        ! ta..........Temperatura no nivel de fonte de calor do dossel (K)
        ! ea..........Pressure of vapor
        ! etc.........Pressure of vapor at top of the copa
        ! etg.........Pressao de vapor no base da copa
        ! btc.........btc(i)=EXP(30.25353  -5418.0  /tc(i))/(tc(i)*tc(i)).
        ! btg.........btg(i)=EXP(30.25353  -5418.0  /tg(i))/(tg(i)*tg(i))
        ! u2..........wind speed at top of canopy
        ! radt........net heat received by canopy/ground vegetation
        ! par.........par incident on canopy
        ! pd..........ratio of par beam to total par
        ! rst ........Resisttencia Estomatica "Stomatal resistence" (s/m)
        ! rsoil.......Resistencia do solo (s/m)
        ! phroot......Soil moisture potentials in root zone of each
        !             vegetation layer and summed soil+root resistance.
        ! hrr.........rel. humidity in top layer
        ! phsoil......soil moisture potential of the i-th soil layer
        ! cc..........heat capacity of the canopy
        ! cg..........heat capacity of the ground
        ! satcap......saturation liquid water capacity         (m)
        ! snow........snow amount
        ! dtc ........dtc(i)=pblsib(i,2,5)*dtc3x
        ! dtg.........dtg(i)=pblsib(i,1,5)*dtc3x
        ! dtm.........dtm(i)=pblsib(i,3,5)*dtc3x
        ! dqm ........dqm(i)=pblsib(i,4,5)*dtc3x
        ! stm.........Variavel utilizada mo cal. da Resisttencia
        ! closs.......Radiation loss from canopy
        ! gloss.......Radiation loss from ground
        ! ect.........Transpiracao no topo da copa (J/m*m)
        ! eci.........Evaporacao da agua interceptada no topo da copa (J/m*m)
        ! egt.........Transpiracao na base da copa (J/m*m)
        ! egi.........Evaporacao da neve (J/m*m)
        ! egs.........Evaporacao do solo arido (J/m*m)
        ! ec..........Soma da Transpiracao e Evaporacao da agua interceptada pelo
        !             topo da copa   ec   (i)=eci(i)+ect(i)
        ! eg..........Soma da transpiracao na base da copa +  Evaporacao do solo arido
        !             +  Evaporacao da neve  " eg   (i)=egt(i)+egs(i)+egi(i)"
        ! hc..........Total sensible heat lost of top from the veggies.
        ! hg..........Total sensible heat lost of base from the veggies.
        ! ecidif......check if interception loss term has exceeded canopy storage
        !             ecidif(i)=MAX(0.0   , eci(i)-capac(i,1)*hlat3 )
        ! egidif......check if interception loss term has exceeded canopy storage
        !             ecidif(i)=MAX(0.0   , egi(i)-capac(i,1)*hlat3 )
        ! ecmass......Mass of water lost of top from the veggies.
        ! egmass......Mass of water lost of base from the veggies.
        ! etmass......Total mass of water lost from the veggies.
        ! hflux.......Total sensible heat lost from the veggies
        ! chf.........Heat fluxes into the canopy  in w/m**2
        ! shf.........Heat fluxes into the ground, in w/m**2
        ! fluxef......Modified to use force-restore heat fluxes
        !             fluxef(i) = shf(i) - cg(i)*dtg(i)*dtc3xi " Garrat pg. 227"
        ! roff........runoff (escoamente superficial e drenagem)(m)
        ! drag........tensao superficial
        ! hgdtg.......n.b. fluxes expressed in joules m-2
        ! hgdtc.......n.b. fluxes expressed in joules m-2
        ! hgdtm.......n.b. fluxes expressed in joules m-2
        ! hcdtg.......n.b. fluxes expressed in joules m-2
        ! hcdtc.......n.b. fluxes expressed in joules m-2
        ! hcdtm.......n.b. fluxes expressed in joules m-2
        ! egdtg.......partial derivative calculation for latent heat
        ! egdtc.......partial derivative calculation for latent heat
        ! egdqm.......partial derivative calculation for latent heat
        ! ecdtg.......partial derivative calculation for latent heat
        ! ecdtc.......partial derivative calculation for latent heat
        ! ecdqm.......partial derivative calculation for latent heat
        ! deadtg
        ! deadtc
        ! deadqm
        ! bps
        ! psb
        ! dzm.........Altura media de referencia  para o vento para o calculo
        !             da estabilidade do escoamento
        ! em..........Pressao de vapor da agua
        ! gmt.........temperature related matrix virtual temperature tendency
        !             due to vertical diffusion
        ! gmq.........specific humidity related matrix specific humidity of
        !             reference (fourier)
        ! gmu.........wind related matrix
        ! cu..........Friction  transfer coefficients.
        ! cuni........Neutral friction transfer  coefficients.
        ! ctni........Neutral heat transfer coefficients.
        ! ustar.......Surface friction velocity  (m/s)
        ! tgeff.......effective ground temperature
        ! cosz........Cosine of zenith angle
        ! rhoair......Desnsidade do ar
        ! psy.........(cp/(hl*epsfac))*psur(i)
        ! rcp.........densidade do ar vezes o calor especifico do ar
        ! wc..........Minimo entre 1 e a razao entre a agua interceptada pelo
        !             indice de area foliar no topo da copa
        ! wg..........Minimo entre 1 e a razao entre a agua interceptada pelo
        !             indice de area foliar na base da copa
        ! fc..........Condicao de oravalho 0 ou 1 na topo da copa
        ! fg..........Condicao de oravalho 0 ou 1 na base da copa
        ! hrr.........rel. humidity in top layer
        ! ssib
        ! yVisBeam.......Downward Surface shortwave fluxe visible beam (cloudy)
        ! yVisDiff.......Downward Surface shortwave fluxe visible diffuse (cloudy)
        ! yNirBeam.......Downward Surface shortwave fluxe Near-IR beam (cloudy)
        ! yNirDiff.......Downward Surface shortwave fluxe Near-IR diffuse (cloudy)
        ! ySwToaDown......swinc....solar input at top of atmosphere
        ! yVisBeamC......Downward Surface shortwave fluxe visible beam (clear)
        ! yVisDiffC......Downward Surface shortwave fluxe visible diffuse (clear)
        ! yNirBeamC......Downward Surface shortwave fluxe Near-IR beam (clear)
        ! yNirDiffC......Downward Surface shortwave fluxe Near-IR diffuse (clear)
        ! cldsav......Cloud cover
        ! cp..........Specific heat of air           (j/kg/k)
        ! hl..........heat of evaporation of water     (j/kg)
        ! rgas........gas constant of dry air        (j/kg/k)
        ! g...........grav   gravity constant        (m/s**2)
        ! solcon......solar constant (wgne value)    (w/m**2)
        ! rmwmd.......fracao molar entre a agua e o ar
        ! swint.......sw subr. call interval in hours
        !             swint has to be less than or equal to trint
        !                              and mod(trint,swint)=0

        ! trint.......ir subr. call interval in hours
        ! yrl.........length of year in days
        ! idate(4)....output : idate(1) = initial hour of
        !                      idate(2) = day of month.
        !                      idate(3) = month of year.
        !                      idate(4) = year.

        ! idatec(4)...output : idatec(1)=current hour of
        !                   idatec(2)=current day of month.
        !                   idatec(3)=current month of year.
        !                   idatec(4)=current year.

        ! kt..........hour of present  time step
        ! ktm.........hour of previous time step
        ! jdt.........time step in getdia
        ! monl(12)....length of each month in days
        ! iswrad........shortwave radiation
        !               irad = NON: excluded
        !               irad = LCH: lacis & hansen
        !               irad = CRD: clirad (chou&lee, modified by tarasova&fomin)
        !               irad = UKM: ukmet office
        ! ilwrad........longwave radiation
        !               irad = NON: excluded
        !               irad = HRS: harshvardhan
        !               irad = CRD: clirad (chou&lee, modified by tarasova&fomin)
        !               irad = UKM: ukmet office
        ! iccon.......the physical process cumulus convection(kuo)
        !               iccon = NON: excluded
        !               iccon = KUO: kuo
        !               iccon = ARA: arakawa
        !               iccon = GRE: grell ensemble
        ! icld........>>> icld = 1    : old cloud emisivity (optical depth) setting
        !                   ccu :       0.05 *dp
        !                   css :       0.025*dp            for ice cloud t<253.0
        !                         0.05 *dp            for ice cloud t>253.0
        !            >>> icld = 2    : new cloud emisivity (optical depth) setting
        !                   ccu :       (0.16)*dp
        !                   css :        0.0                        t<-82.5c
        !                         (2.0e-6*(t-tcrit)**2)*dp    -82.5<t<-10.0c
        !                         (6.949e-3*(t-273)+.08)*dp   -10.0<t< 0.0c
        !                         (0.08)*dp                 -10.0<t< 0.0c
        !            >>> icld = 3    : ccm3 based cloud emisivity

        ! co2val......co2val is wgne standard value in ppm "co2val = /345.0/
        ! delt........time interval in sec (fixed throuh the integration)
        ! filta.......weight used on central time
        !              step of robert time filter.
        !              set in main routine "smf".filta=0.92e0
        ! nfin0.......input  file at time level t-dt
        ! nfin1.......input  file at time level t
        ! initlz......constant initlz=2.
        ! nfcnv0......initial information on convective clouds for int. radiation
        ! nfcldr......constant nfcldr = 74
        ! tbase.......constant tbase =  273.15e00
        ! latco.......latitude
        ! dodia.......Variable logical for search for combined field components.
        ! lvavl.......levels in available diagnostic (1 or kmax)
        ! nuavl.......unit code of available diagnostic
        ! itavl.......type of available diagnostic (1 gaussian, 2 spectral)
        ! iavrq.......Number of requested diagnostic
        ! ixavl.......Number available diagnostic components for combined fields
        ! inavl.......Number available diagnostic similar requested diagnostic
        ! iclcd.......requested diagnostic calculation code (0 direct
        !             calculation, > 0 add to requested field number iclcd,
        !             < 0 subtract from requested field number -iclcd )
        ! nucf........nurq  = unit code of requested diagnostic
        ! ixcf........Number of requested diagnostic
        ! incf........combined fields
        ! kravl.......Number of available diagnostics equivalence
        !              the desired diagnostic
        ! krcf........combined fields
        ! jrcf........combined fields
        ! jdt.........time step in getdia
        ! latco.......latco grid point reference the latitude
        ! cdhl........logical indicator for dhn output prognostics
        ! ustr........surface zonal stress umom(i)=fmom*um(ncount),
        !               where .fmom  momentum flux      in n/m**2
        !               fmom= rhoair(ncount)*cu(ncount)*ustar(ncount)
        !               um  (ncount)=gu (i,1)/sinclt
        !               gu          = (zonal velocity)*sin(colat)

        ! vstr........surface meridional stress.vmom(i)=rho(i)*gv(i)*rmi(i)
        !                           rho  (i)=gps(i)/(gr100*gt(i))
        !                           gr100 =gasr*0.01
        ! first.......control logical variable .true. or .false.
        ! mxrdcc......use maximum random converage for radiative conv. clouds
        !               constant logical mxrdcc = .true.
        ! lcnvl.......the lowest layer index where non-convective clouds can
        !               occur (ben says this should be 2 or more)
        !               constant lcnvl = 2
        ! lthncl......Minimum depth in mb of non-zero low level cloud
        !             consta lthncl=80
        ! convc.......ncols convective cloud cover in 3 hr. avrage
        ! convt.......ncols convective cloud top  (sigma layer)
        ! convb.......ncols convective cloud base (sigma layer)
        ! convts
        ! convcs
        ! convbs
        ! sigki ......sigki (k)=1.0e0/EXP(rk*LOG(sig(k))),  where "sig"
        !             sigma coordinate at middle of layer and rk=gasr/cp
        ! xVisBeam.......Downward Surface shortwave fluxe visible beam (cloudy)
        ! xVisDiff.......Downward Surface shortwave fluxe visible diffuse (cloudy)
        ! xNirBeam.......Downward Surface shortwave fluxe Near-IR beam (cloudy)
        ! xNirDiff.......Downward Surface shortwave fluxe Near-IR diffuse (cloudy)
        ! xswtop......shortwave upward at top  or  shortwave upward at top (clear)
        ! xVisBeamC......Downward Surface shortwave fluxe visible beam (clear)
        ! xVisDiffC......Downward Surface shortwave fluxe visible diffuse (clear)
        ! xNirBeamC......Downward Surface shortwave fluxe Near-IR beam (clear)
        ! xNirDiffC......Downward Surface shortwave fluxe Near-IR diffuse (clear)
        !==========================================================================
        !

        ! Time info

        real(KIND = r8) :: fator, pfator

        ! Model Geometry
        REAL(KIND = r8), INTENT(IN) :: colrad(ncols)
        REAL(KIND = r8), INTENT(IN) :: lonrad(ncols)
        REAL(KIND = r8), INTENT(IN) :: cos2d  (ncols)
        REAL(KIND = r8), INTENT(IN) :: sig   (kmax)
        REAL(KIND = r8), INTENT(IN) :: sigml (kmax + 1)
        REAL(KIND = r8), INTENT(IN) :: delsig(kmax)

        ! Model information
        INTEGER, INTENT(IN) :: latco
        INTEGER, INTENT(IN) :: ncols
        INTEGER, INTENT(IN) :: kmax
        INTEGER(KIND = i8), INTENT(INOUT) :: imask (ncols)

        ! Atmospheric fields
        REAL(KIND = r8), INTENT(IN) :: gps   (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: gt    (ncols, kmax)
        REAL(KIND = r8), INTENT(INOUT) :: gq    (ncols, kmax)
        REAL(KIND = r8), INTENT(IN) :: gu    (ncols, kmax)
        REAL(KIND = r8), INTENT(IN) :: gv    (ncols, kmax)
        REAL(KIND = r8), INTENT(IN) :: omg   (ncols, kmax)
        REAL(KIND = r8), INTENT(INOUT) :: tsea  (ncols)
        REAL(KIND = r8), OPTIONAL, INTENT(INOUT) :: gicem (ncols, kmax)
        REAL(KIND = r8), OPTIONAL, INTENT(INOUT) :: gicet (ncols, kmax)
        REAL(KIND = r8), OPTIONAL, INTENT(INOUT) :: gliqm (ncols, kmax)
        REAL(KIND = r8), OPTIONAL, INTENT(INOUT) :: gliqt (ncols, kmax)
        REAL(KIND = r8), OPTIONAL, INTENT(INOUT) :: gvarm (ncols, kmax, nClass + nAeros)
        REAL(KIND = r8), OPTIONAL, INTENT(INOUT) :: gvart (ncols, kmax, nClass + nAeros)
        ! Mass and energy turbulent diffusion coefficients
        REAL(KIND = r8), INTENT(INOUT) :: PBL_CoefKm(ncols, kmax + 1)
        REAL(KIND = r8), INTENT(INOUT) :: PBL_CoefKh(ncols, kmax + 1)

        ! variables for specified surface (specSfc)
        REAL(KIND = r8), OPTIONAL, INTENT(inout) :: sgTsfc(ncols)
        REAL(KIND = r8), OPTIONAL, INTENT(inout) :: sgH(ncols)
        REAL(KIND = r8), OPTIONAL, INTENT(inout) :: sgLE(ncols)
        REAL(KIND = r8), OPTIONAL, INTENT(inout) :: sgTau(ncols)

        ! SURFACE:  albedo
        REAL(KIND = r8), INTENT(INOUT) :: AlbVisBeam (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: AlbVisDiff (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: AlbNirBeam (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: AlbNirDiff (ncols)

        ! Radiation fields at last integer hour
        REAL(KIND = r8), INTENT(INOUT) :: rSwToaDown(ncols)
        REAL(KIND = r8), INTENT(INOUT) :: rVisBeam (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: rVisDiff (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: rNirBeam (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: rNirDiff (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: rVisBeamC(ncols)
        REAL(KIND = r8), INTENT(INOUT) :: rVisDiffC(ncols)
        REAL(KIND = r8), INTENT(INOUT) :: rNirBeamC(ncols)
        REAL(KIND = r8), INTENT(INOUT) :: rNirDiffC(ncols)
        REAL(KIND = r8), INTENT(INOUT) :: rSwSfcNet (ncols)! Abs Sfc SW
        REAL(KIND = r8), INTENT(INOUT) :: rSwSfcNetC(ncols)! Abs Sfc SW (clear)

        ! Radiation fields at next integer hour
        REAL(KIND = r8), INTENT(INOUT) :: ySwToaDown(ncols)
        REAL(KIND = r8), INTENT(INOUT) :: yVisBeam (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: yVisDiff (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: yNirBeam (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: yNirDiff (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: yVisBeamC(ncols)
        REAL(KIND = r8), INTENT(INOUT) :: yVisDiffC(ncols)
        REAL(KIND = r8), INTENT(INOUT) :: yNirBeamC(ncols)
        REAL(KIND = r8), INTENT(INOUT) :: yNirDiffC(ncols)
        REAL(KIND = r8), INTENT(INOUT) :: ySwHeatRate   (ncols, kmax)
        REAL(KIND = r8), INTENT(INOUT) :: ySwHeatRateC  (ncols, kmax)
        REAL(KIND = r8), INTENT(INOUT) :: ySwSfcNet (ncols)! Abs Sfc SW
        REAL(KIND = r8), INTENT(INOUT) :: ySwSfcNetC(ncols)! Abs Sfc SW (clear)

        ! LW Radiation fields at last integer hour
        REAL(KIND = r8), INTENT(INOUT) :: LwCoolRate (ncols, kmax)
        REAL(KIND = r8), INTENT(INOUT) :: LwSfcDown  (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: LwSfcNet   (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: LwToaUp    (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: LwCoolRateC(ncols, kmax)
        REAL(KIND = r8), INTENT(INOUT) :: LwSfcDownC (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: LwSfcNetC  (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: LwToaUpC   (ncols)

        ! CONVECTION: convective clouds
        REAL(KIND = r8), INTENT(IN) :: convc (ncols)
        REAL(KIND = r8), INTENT(IN) :: convt (ncols)
        REAL(KIND = r8), INTENT(IN) :: convb (ncols)
        REAL(KIND = r8), INTENT(IN) :: convts(ncols)
        REAL(KIND = r8), INTENT(IN) :: convcs(ncols)
        REAL(KIND = r8), INTENT(IN) :: convbs(ncols)

        ! Cloud field
        REAL(KIND = r8), INTENT(INOUT) :: cldsav(ncols)
        REAL(KIND = r8), INTENT(INOUT) :: cldtot(ncols, kmax)
        REAL(KIND = r8), INTENT(INOUT) :: cldinv(ncols, kmax)
        REAL(KIND = r8), INTENT(INOUT) :: cldsat(ncols, kmax)
        REAL(KIND = r8), INTENT(INOUT) :: cldcon(ncols, kmax)
        REAL(KIND = r8), INTENT(INOUT) :: cldson(ncols, kmax)

        ! Microphysics
        REAL(KIND = r8), INTENT(INOUT) :: clwd  (ncols, kmax)
        REAL(KIND = r8), INTENT(INOUT) :: emisd (ncols, kmax)
        REAL(KIND = r8), INTENT(INOUT) :: taud  (ncols, kmax)

        ! Chemistry
        REAL(KIND = r8), INTENT(INOUT) :: o3mix(ncols, kMax)

        REAL(KIND = r8), INTENT(INOUT) :: sm0   (ncols, 3) ! solange add 13-11-2012
        REAL(KIND = r8), INTENT(IN) :: ppli  (ncols)
        REAL(KIND = r8), INTENT(IN) :: ppci  (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: gyu   (ncols, kmax)
        REAL(KIND = r8), INTENT(INOUT) :: gyv   (ncols, kmax)
        REAL(KIND = r8), INTENT(INOUT) :: gtd   (ncols, kmax)
        REAL(KIND = r8), INTENT(INOUT) :: gqd   (ncols, kmax)
        INTEGER, INTENT(IN) :: ifday
        REAL(KIND = r8), INTENT(IN) :: tod
        REAL(KIND = r8), INTENT(INOUT) :: gl0   (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: zorl  (ncols)
        REAL(KIND = r8), INTENT(IN) :: gtt   (ncols, kmax)
        REAL(KIND = r8), INTENT(IN) :: gqq   (ncols, kmax)
        REAL(KIND = r8), INTENT(IN) :: sheleg(ncols)
        REAL(KIND = r8), INTENT(INOUT) :: tseam (ncols)
        REAL(KIND = r8), INTENT(IN) :: ssib(ncols)
        REAL(KIND = r8), INTENT(INOUT) :: bstar(ncols)
        INTEGER(KIND = i8), INTENT(INOUT) :: mlsi  (ncols) ! solange add 13-11-2012
        !
        !     this is for interpolating shortwave rad at ground

        REAL(KIND = r8), INTENT(INOUT) :: pblh (ncols)
        REAL(KIND = r8), INTENT(IN) :: cu_hr  (ncols, kmax)
        INTEGER, INTENT(IN) :: cu_kbot(ncols)
        INTEGER, INTENT(IN) :: cu_ktop(ncols)
        INTEGER, INTENT(IN) :: cu_Kuo (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: dudt (ncols, kmax)
        REAL(KIND = r8), INTENT(INOUT) :: dvdt (ncols, kmax)

        REAL(KIND = r8), INTENT(INOUT) :: EFFCS(ncols, kmax)
        REAL(KIND = r8), INTENT(INOUT) :: EFFIS(ncols, kmax)
        REAL(KIND = r8), INTENT(INOUT) :: tauresx(ncols)
        REAL(KIND = r8), INTENT(INOUT) :: tauresy(ncols)
        !
        !     these are for monitoring of gpv in gfidi.
        !

        REAL(KIND = r8), INTENT(INOUT) :: ustr(ncols)
        REAL(KIND = r8), INTENT(INOUT) :: vstr(ncols)

        REAL(KIND = r8), INTENT(IN) :: sigki (kmax)

        REAL(KIND = r8), INTENT(INOUT) :: ps    (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: var   (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: sens  (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: evap  (ncols)
        REAL(KIND = r8), INTENT(IN) :: topog(ncols)
        LOGICAL, INTENT(IN) :: intcosz

        REAL(KIND = r8), INTENT(INOUT) :: Mmlen(ncols)

        REAL(KIND = r8), PARAMETER :: alon = 0.0_r8
        REAL(KIND = r8), INTENT(IN) :: xland(ncols)
        INTEGER, INTENT(IN) :: lowlyr(ncols)
        REAL(KIND = r8), INTENT(INOUT) :: ustar(ncols)
        REAL(KIND = r8), INTENT(INOUT) :: z0(ncols)
        REAL(KIND = r8), INTENT(INOUT) :: tkemyj(ncols, kMax + 1)
        REAL(KIND = r8), INTENT(INOUT) :: snow (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: thz0 (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: qz0  (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: uz0  (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: vz0  (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: akhs (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: akms (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: ct   (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: htdisp(ncols)
        REAL(KIND = r8), INTENT(OUT) :: temp2m(ncols)
        REAL(KIND = r8), INTENT(OUT) :: umes2m(ncols)
        INTEGER(KIND = i8), INTENT(IN) :: mskant(ncols)
        REAL(KIND = r8), INTENT(INOUT) :: tpert(ncols)
        REAL(KIND = r8), INTENT(INOUT) :: qpert (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: sflux_t(ncols)
        REAL(KIND = r8), INTENT(INOUT) :: sflux_r(ncols)
        REAL(KIND = r8), INTENT(INOUT) :: sflux_u(ncols)
        REAL(KIND = r8), INTENT(INOUT) :: sflux_v(ncols)
        REAL(KIND = r8), INTENT(INOUT) :: tstar(ncols)
        REAL(KIND = r8), INTENT(INOUT) :: wstar(ncols)
        REAL(KIND = r8), INTENT(INOUT) :: r_aer(ncols)
        REAL(KIND = r8), INTENT(OUT) :: veg_type(ncols, npatches)
        REAL(KIND = r8), INTENT(OUT) :: frac_occ(ncols, npatches)

        REAL(KIND = r8), INTENT(INOUT) :: ndvi  (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: ndvim (ncols)
        REAL(KIND = r8), INTENT(IN) :: qliq  (ncols, kMax)

        REAL(KIND = r8), INTENT(INOUT) :: HML  (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: HUML (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: HVML (ncols)
        REAL(KIND = r8), INTENT(INOUT) :: z0sea(ncols)
        REAL(KIND = r8), INTENT(INOUT) :: TSK  (ncols)

        REAL(KIND = r8) :: QCF(ncols, kMax)
        REAL(KIND = r8) :: QCL(ncols, kMax)
        REAL(KIND = r8) :: QCR(ncols, kMax)
        ! Mass and energy turbulent diffusion coefficients
        REAL(KIND = r8) :: PBLSFC_CoefKm(ncols, kmax)
        REAL(KIND = r8) :: PBLSFC_CoefKh(ncols, kmax)
        REAL(KIND = r8) :: tkesfc(ncols, kMax + 1)

        REAL(KIND = r8) :: taux(1:nCols)
        REAL(KIND = r8) :: tauy(1:nCols)
        REAL(KIND = r8) :: SICE(ncols)
        INTEGER :: nmax
        INTEGER :: nsx(ncols)
        INTEGER :: itype(ncols)
        !
        !     prognostic variables
        !
        REAL(KIND = r8) :: qsfc  (ncols)
        REAL(KIND = r8) :: tsfc  (ncols)
        !
        !     variables calculated from above and ambient conditions
        !
        !
        !     heat fluxes : c-canopy, g-ground, t-trans, e-evap  in j m-2
        !
        !
        !     this is for coupling with closure turbulence model
        !
        REAL(KIND = r8) :: cosz   (ncols)
        REAL(KIND = r8) :: tmtx   (ncols, kmax, 3)
        REAL(KIND = r8) :: qmtx   (ncols, kmax, 5 + nClass + nAeros)
        REAL(KIND = r8) :: umtx   (ncols, kmax, 4)

        REAL(KIND = r8) :: tmsfc  (ncols, kmax, 3)
        REAL(KIND = r8) :: qmsfc  (ncols, kmax, 5 + nClass + nAeros)
        REAL(KIND = r8) :: umsfc  (ncols, kmax, 4)
        REAL(KIND = r8) :: gicem_local  (ncols, kmax)
        REAL(KIND = r8) :: gliqm_local  (ncols, kmax)
        REAL(KIND = r8) :: gvarm_local  (ncols, kmax, nClass + nAeros)
        REAL(KIND = r8) :: tsurf (ncols)
        REAL(KIND = r8) :: qsurf (ncols)
        REAL(KIND = r8) :: umom  (ncols)
        REAL(KIND = r8) :: vmom  (ncols)
        REAL(KIND = r8) :: slrad (ncols)
        REAL(KIND = r8) :: zenith(ncols)
        REAL(KIND = r8) :: QSS   (ncols)
        REAL(KIND = r8) :: sdelt
        REAL(KIND = r8) :: ratio
        REAL(KIND = r8) :: etime
        REAL(KIND = r8) :: xday

        ! Radiation field (Interpolated) at time = tod
        REAL(KIND = r8) :: xVisBeam (ncols)
        REAL(KIND = r8) :: xVisDiff (ncols)
        REAL(KIND = r8) :: xNirBeam (ncols)
        REAL(KIND = r8) :: xNirDiff (ncols)

        ! SSIB INIT: Solar radiation with cos2
        REAL(KIND = r8) :: ssib_VisBeam (ncols)
        REAL(KIND = r8) :: ssib_VisDiff (ncols)
        REAL(KIND = r8) :: ssib_NirBeam (ncols)
        REAL(KIND = r8) :: ssib_NirDiff (ncols)

        INTEGER :: k
        INTEGER :: kk, m
        INTEGER :: i
        INTEGER :: ncount
        REAL(KIND = r8) :: deltm
        REAL(KIND = r8) :: sindel
        REAL(KIND = r8) :: cosdel
        REAL(KIND = r8) :: fimxi
        REAL(KIND = r8) :: ctime
        REAL(KIND = r8) :: cos2  (ncols)
        REAL(KIND = r8) :: frh
        REAL(KIND = r8) :: btime
        REAL(KIND = r8) :: atime
        REAL(KIND = r8) :: tice01
        INTEGER :: month   (ncols)
        INTEGER :: month2  (ncols)
        REAL(KIND = r8) :: colrad2 (ncols)
        REAL(KIND = r8) :: zenith1 (ncols)
        REAL(KIND = r8) :: zenith2 (ncols)
        REAL(KIND = r8) :: sinclt2 (ncols)
        REAL(KIND = r8) :: chug  (ncols, kmax)
        REAL(KIND = r8) :: chvg  (ncols, kmax)
        REAL(KIND = r8) :: chtg  (ncols, kmax)
        REAL(KIND = r8) :: xdrag (ncols)
        REAL(KIND = r8) :: ydrag (ncols)
        REAL(KIND = r8) :: ELM(ncols, kMax)
        REAL(KIND = r8) :: topo (ncols)
        REAL(KIND = r8) :: wind (ncols)
        REAL(KIND = r8) :: rho (ncols)


        !------------------------------------------------------
        taux = 0.0_r8
        tauy = 0.0_r8
        SICE = 0.0_r8
        nsx = 0
        itype = 0
        !
        !     prognostic variables
        !
        qsfc = 0.0_r8
        tsfc = 0.0_r8
        !
        !     variables calculated from above and ambient conditions
        !
        !
        !     heat fluxes : c-canopy, g-ground, t-trans, e-evap  in j m-2
        !
        !
        !     this is for coupling with closure turbulence model
        !
        cosz = 0.0_r8
        tmtx = 0.0_r8
        qmtx = 0.0_r8
        umtx = 0.0_r8
        tmsfc = 0.0_r8
        qmsfc = 0.0_r8
        umsfc = 0.0_r8
        tsurf = 0.0_r8
        qsurf = 0.0_r8
        umom = 0.0_r8
        vmom = 0.0_r8
        slrad = 0.0_r8
        zenith = 0.0_r8
        QSS = 0.0_r8
        sdelt = 0.0_r8
        ratio = 0.0_r8
        etime = 0.0_r8
        xday = 0.0_r8

        ! Radiation field (Interpolated) at time = tod
        xVisBeam = 0.0_r8
        xVisDiff = 0.0_r8
        xNirBeam = 0.0_r8
        xNirDiff = 0.0_r8

        ! SSIB INIT: Solar radiation with cos2
        ssib_VisBeam = 0.0_r8
        ssib_VisDiff = 0.0_r8
        ssib_NirBeam = 0.0_r8
        ssib_NirDiff = 0.0_r8

        deltm = 0.0_r8
        sindel = 0.0_r8
        cosdel = 0.0_r8
        fimxi = 0.0_r8
        ctime = 0.0_r8
        cos2 = 0.0_r8
        frh = 0.0_r8
        btime = 0.0_r8
        atime = 0.0_r8
        tice01 = 0.0_r8
        month = 0
        month2 = 0
        colrad2 = 0.0_r8
        zenith1 = 0.0_r8
        zenith2 = 0.0_r8
        sinclt2 = 0.0_r8
        chug = 0.0_r8
        chvg = 0.0_r8
        chtg = 0.0_r8
        xdrag = 0.0_r8
        ydrag = 0.0_r8
        ELM = 0.0_r8
        topo = 0.0_r8

        DO  k = 1, kmax
            DO i = 1, ncols
                gyu(i, k) = gyu(i, k) + (dudt (i, k) * SIN(colrad(i)))
                gyv(i, k) = gyv(i, k) + (dvdt (i, k) * SIN(colrad(i)))
                dudt (i, k) = 0.0_r8
                dvdt (i, k) = 0.0_r8
            END DO
        END DO

        IF(dogwd.EQ.0)THEN

            CALL Gwdd_Driver(ps, gu, gv, gt, gq, chug, chvg, chtg, xdrag, ydrag, &
                    var, varcut, sigml, sig, delsig, ncols, kmax, latco, dt, imask, colrad, topog, &
                    pblh, cu_hr, cu_kbot, cu_ktop, cu_Kuo)

            DO  k = 1, kmax
                DO i = 1, ncols
                    gyu(i, k) = gyu(i, k) - chug(i, k)
                    gyv(i, k) = gyv(i, k) - chvg(i, k)
                    gtd(i, k) = gtd(i, k) + chtg(i, k)
                END DO
            END DO
        ENDIF
        !
        dtc3x = dt * REAL(intg, r8)
        !
        DO k = 1, kmax
            DO i = 1, ncols
                gq(i, k) = MAX(1.0e-12_r8, gq(i, k))
                gt(i, k) = gt(i, k) / (1.0e0_r8 + 0.608e0_r8 * gq(i, k))
            END DO
        END DO

        IF(initlz.GE.0.AND.kt.EQ.0.AND.jdt.EQ.1) THEN
            ncount = 0
            DO i = 1, ncols
                TSfc0(i, latco) = gt(i, 1)
                QSfc0(i, latco) = gq(i, 1)
                TSfcm(i, latco) = gt(i, 1)
                QSfcm(i, latco) = gq(i, 1)
                IF(imask(i).GE.1_i8) THEN
                    ncount = ncount + 1
                    tc0(ncount, latco) = gt(i, 1)
                    tg0(ncount, latco) = gt(i, 1)
                    tcm(ncount, latco) = gt(i, 1)
                    tgm(ncount, latco) = gt(i, 1)
                    tm0(ncount, latco) = gt(i, 1)
                    tmm(ncount, latco) = gt(i, 1)
                    qm0(ncount, latco) = gq(i, 1)
                    qmm(ncount, latco) = gq(i, 1)
                    IF(sheleg(i).GT.0.0e0_r8) THEN
                        tg0(ncount, latco) = MIN(tg0(ncount, latco), tf - 0.01e0_r8)
                        tgm(ncount, latco) = MIN(tgm(ncount, latco), tf - 0.01e0_r8)
                    END IF
                END IF
            END DO
        END IF
        DO k = 1, npatches
            DO i = 1, ncols
                veg_type(i, k) = REAL(imask(i), kind = r8)
                IF(k==1)THEN
                    IF(imask(i).GE.1_i8) THEN
                        frac_occ(i, k) = 0.0_r8
                    ELSE
                        frac_occ(i, k) = 1.0_r8
                    END IF
                ELSE
                    IF(imask(i).GE.1_i8) THEN
                        frac_occ(i, k) = 1.0_r8
                    ELSE
                        frac_occ(i, k) = 0.0_r8
                    END IF
                END IF
            END DO
        END DO
        ncount = 0
        DO i = 1, ncols
            IF(imask(i).GE.1_i8) THEN
                ncount = ncount + 1
                itype(ncount) = INT(imask(i))
            END IF
        END DO
        nmax = ncount
        !
        !     mon is the month used for vegetation data input
        !
        DO i = 1, ncols, 1
            month(i) = idatec(2)
            IF((((colrad(i) * 180.0_r8) / 3.1415926e0_r8) - 90.0_r8)  > 0.0_r8) THEN
                month(i) = month(i) + 6
                IF(month(i).GE.13) month(i) = month(i) - 12
            END IF
        END DO

        ncount = 0

        DO i = 1, ncols
            IF(imask(i).GE.1_i8) THEN
                ncount = ncount + 1
                month2 (ncount) = month(i)
                colrad2(ncount) = colrad(i)
                sinclt2(ncount) = SIN(colrad(i))
            END IF
        END DO
        !
        !     computation of astronomical parameters
        !     sdelt ;solar inclination
        !     etime ;correction factor to local time
        !     ratio ;factor relating to the distance between the earth and the sun
        !
        CALL radtim(idatec, sdelt, ratio, etime, tod, xday, yrl)

        sindel = SIN(sdelt)
        cosdel = COS(sdelt)
        fimxi = 24.0e0_r8 / 360.0_r8
        ctime = alon / 15.0e0_r8
        cos2 = 0.0e0_r8
        ncount = 0
        frh = (MOD(tod + 0.03125_r8, 3600.0_r8) - 0.03125_r8) / 3600.0_r8

        DO i = 1, ncols
            zenith1(i) = sindel * COS(colrad(i))
            wind  (i) = sqrt((gu (i, 1) / SIN(colrad(i)))**2 + (gv (i, 1) / SIN(colrad(i)))**2)
        ENDDO

        DO i = 1, ncols
            btime = fimxi * lonrad(i) + ctime
            atime = etime + pai12 * (12.0_r8 - idatec(1) - frh - btime)
            zenith2 (i) = cosdel * SIN(colrad(i)) * COS(atime)
            zenith  (i) = zenith1(i) + zenith2(i)
        END DO

        IF(ncount.EQ.0) ncount = 1
        IF(intcosz)THEN
            !cos2=cos2/REAL(ncount,r8)!!!!mudanca forcada
            cos2(1:ncols) = cos2d(1:ncols)
        ELSE
            cos2(1:ncols) = zenith(1:ncols)
        END IF
        ncount = 0
        DO i = 1, ncols
            !erg       wind(i)=sqrt((gu (i,1)/SIN( colrad(i)))**2 + (gv (i,1)/SIN( colrad(i)))**2)
            IF(imask(i).GE.1_i8) THEN
                ncount = ncount + 1
                cosz(ncount) = zenith(i)
            END IF
        END DO
        !
        !     sib setting  *phenology*
        !
        IF(nmax.GE.1) THEN
            IF(schemes==1)THEN
                nsx = 0
                CALL Phenology(latco, nCols, nmax, itype, colrad2, month2, xday, idatec, nsx)
            ELSE IF(schemes==2)THEN
                CALL Phenology_sib2(latco, nCols, idatec, ndvi, ndvim, colrad2, itype)
            END IF
        END IF
        !
        !     surface albedo (vis/nir and beam/diffuse)
        !     extinction coefficients
        !
        IF(schemes==1) THEN
            CALL Albedo(&
                    ! Model information
                    ncols, kMax, latco, &
                    nmax, nsx       (1:nCols), itype(1:nCols), &
                    imask     (1:nCols), &
                    ! Model Geometry
                    cosz      (1:nCols), zenith    (1:nCols), &
                    ! Time info
                    month2    (1:nCols), month     (1:nCols), &
                    ! Microphysics
                    taud      (1:ncols, 1:kMax), &
                    ! Atmospheric fields
                    wind      (1:nCols), tsea      (1:nCols), &
                    ! LW Radiation fields at last integer hour
                    LwSfcDown (1:nCols), &
                    ! SW Radiation field (Interpolated) at time = tod
                    xVisBeam  (1:nCols), xVisDiff   (1:nCols), xNirBeam  (1:nCols), &
                    xNirDiff  (1:nCols), &
                    ! Surface Albedo
                    AlbVisBeam(1:nCols), AlbVisDiff (1:nCols), AlbNirBeam(1:nCols), &
                    AlbNirDiff(1:nCols))
        ELSE IF(schemes==2) THEN
            CALL Albedo_sib2(&
                    ! Model information
                    ncols, kMax, latco, &
                    ! Model Geometry
                    cosz      (1:nCols), zenith      (1:nCols), &
                    ! Time info
                    month2    (1:nCols), month     (1:nCols), &
                    ! Atmospheric fields
                    wind      (1:nCols), tsea        (1:nCols), &
                    ! Microphysics
                    taud      (1:ncols, 1:kMax), &
                    ! LW Radiation fields at last integer hour
                    LwSfcDown (1:nCols), &
                    ! SW Radiation field (Interpolated) at time = tod
                    xVisBeam  (1:nCols), xVisDiff   (1:nCols), xNirBeam  (1:nCols), &
                    xNirDiff  (1:nCols), &
                    ! Surface Albedo
                    AlbVisBeam(1:nCols), AlbVisDiff(1:nCols), AlbNirBeam(1:nCols), &
                    AlbNirDiff(1:nCols))
        ELSE IF(schemes==3) THEN
            CALL Albedo_IBIS(&
                    ! Model information
                    latco, nCols, kMax, &
                    imask     (1:nCols), &
                    ! Model Geometry
                    zenith      (1:nCols), &
                    ! Time info
                    month2    (1:nCols), month     (1:nCols), &
                    ! Atmospheric fields
                    wind      (1:nCols), tsea        (1:nCols), &
                    ! Microphysics
                    taud      (1:ncols, 1:kMax), &
                    ! LW Radiation fields at last integer hour
                    LwSfcDown (1:nCols), &
                    ! SW Radiation field (Interpolated) at time = tod
                    xVisBeam  (1:nCols), xVisDiff   (1:nCols), xNirBeam  (1:nCols), &
                    xNirDiff  (1:nCols), &
                    ! Surface Albedo
                    AlbVisBeam(1:nCols), AlbVisDiff(1:nCols), AlbNirBeam(1:nCols), &
                    AlbNirDiff(1:nCols))

        END IF

        tice01 = tice + 0.01_r8
        DO i = 1, ncols
            IF(TRIM(OCFLUX) == 'WGFS')THEN
                IF(omlmodel)THEN
                    tsurf(i) = ABS(TSK(i))
                    IF (tsea(i) < 0.0_r8 .AND. ABS(tsea(i)) < tice + 0.01_r8) THEN
                        tsurf(i) = ABS(tsea(i))
                    END IF
                ELSE
                    tsurf(i) = ABS(tsea(i))
                END IF
            ELSE
                tsurf(i) = ABS(tsea(i))
            END IF
            qsurf(i) = 0.622e0_r8 * EXP(21.65605e0_r8 - 5418.0e0_r8 / tsurf(i)) / gps(i)
            !
            IF(initlz.GE.0.AND.kt.EQ.0.AND.ktm.EQ.-1) THEN
                IF(tsea(i).GT.0.0e0_r8.OR.(tsea(i).LT.0.0e0_r8.AND.ABS(tsea(i)).LT.tice01)) THEN
                    tsurf(i) = gt(i, 1)
                END IF
                IF(tsea(i).GT.0.0e0_r8.AND.sheleg(i).GT.0.0e0_r8) THEN
                    tsurf(i) = MIN(tf, tsurf(i))
                END IF
                IF(tsea(i).LT.0.0e0_r8.AND.ABS(tsea(i)).LT.tice01) THEN
                    tsurf(i) = MIN(tice, tsurf(i))
                END IF
            END IF
            qsurf(i) = 0.622e0_r8 * EXP(21.65605e0_r8 - 5418.0e0_r8 / tsurf(i)) / gps(i)
        END DO

        print*, 'shape qcf gicem', shape(QCF), shape(gicem), shape(gliqm) !, ncols, kmax
        print*, 'shape qcr gvarm', shape(QCR), shape(gvarm)
        IF (microphys) THEN
            QCF = gicem
            QCL = gliqm
            QCR = 0.0_r8
            IF((nClass + nAeros)>0 .and. PRESENT(gvarm))QCR(1:ncols, 1:kmax) = gvarm   (1:ncols, 1:kmax, 1)
        ELSE
            QCF = 0.0_r8
            QCL = 0.0_r8
            QCR = 0.0_r8
        END IF
        !
        !     radiation parameterization
        !
        CALL RadiationDriver (&
                ! Run Flags
                first, ifday, lcnvl, lthncl, nfin0, nfin1, nfcnv0, &
                intcosz, kt, mxrdcc, &
                ! Time info
                yrl, idatec, idate, tod, jdt, delt, &
                trint, swint, &
                ! Model Geometry
                colrad, lonrad, zenith, cos2d, &
                sig, sigml, delsig, &
                ! Model information
                latco, ncols, kmax, nls, nlcs, imask, &
                ! Atmospheric fields
                gps, gtt, gqq, tsurf, omg, tsea, &
                QCF, QCL, QCR, &
                ! CONVECTION: convective clouds
                convts, convcs, convbs, convc, convt, convb, &
                ! Surface Albedo
                AlbVisDiff, AlbNirDiff, AlbVisBeam, AlbNirBeam, &
                ! SW Radiation fields at last integer hour
                rSwToaDown, &
                rVisDiff, rNirDiff, rVisBeam, rNirBeam, &
                rVisDiffC, rNirDiffC, rVisBeamC, rNirBeamC, &
                rSwSfcNet, rSwSfcNetC, &
                ! SW Radiation fields at next integer hour
                ySwToaDown, &
                yVisDiff, yNirDiff, yVisBeam, yNirBeam, &
                yVisDiffC, yNirDiffC, yVisBeamC, yNirBeamC, &
                ySwHeatRate, ySwHeatRateC, &
                ySwSfcNet, ySwSfcNetC, &
                ! Radiation field (Interpolated) at time = tod
                xVisDiff, xNirDiff, xVisBeam, xNirBeam, &
                ! LW Radiation fields at last integer hour
                LwCoolRate, LwSfcDown, LwSfcNet, LwToaUp, &
                LwCoolRateC, LwSfcDownC, LwSfcNetC, LwToaUpC, &
                ! SSIB: Total radiation absorbed at ground
                slrad, &
                ! SSIB INIT: Solar radiation with cos2
                ssib_VisBeam, ssib_VisDiff, ssib_NirBeam, ssib_NirDiff, &
                ! Cloud field
                cldsav, cldtot, &
                cldinv, cldsat, cldcon, cldson, &
                ! Microphysics
                clwd, emisd, taud, EFFCS, EFFIS, &
                ! Chemistry
                co2val, o3mix)
        !
        !     yamada-mellor pbl surface parameterization ySwSfcNet,LwSfcNet
        !
        deltm = 0.5e0_r8 * dtc3x
        Mmlen = gl0
        tmsfc = 0.0_r8
        qmsfc = 0.0_r8
        umsfc = 0.0_r8
        !PRINT*,' SfcPBL_Drive'

        ! print specSfc values
        IF (specSfc) THEN
            !enable lines below to allow specSfc effects
            !tsurf(1:nCols)=sgTsfc
            !thz0 = sgTsfc
            !sens(1:nCols)=sgH
            !evap(1:nCols)=sgLE
            !!tsurf0(1:nCols,:)=sgTsfc
        ENDIF

        print*, 'before SfcPBL_Driver', shape(gicem), shape(gliqm), shape(gvarm)
        IF (microphys) THEN
            IF((nClass + nAeros)>0 .and. PRESENT(gvarm))THEN
                CALL SfcPBL_Driver(&
                        gu    (1:nCols, 1:kMax), gv  (1:nCols, 1:kMax), gt(1:nCols, 1:kMax), gq(1:nCols, 1:kMax), &
                        delsig(1:kmax), ncols, kmax, deltm, &
                        colrad(1:nCols), tmsfc(1:nCols, 1:kMax, 1:3), qmsfc(1:nCols, 1:kMax, 1:5 + nClass + nAeros), umsfc(1:nCols, 1:kMax, 1:4), &
                        Mmlen (1:ncols), dt, jdt, &
                        topog (1:ncols), xland (1:ncols), lowlyr(1:ncols), &
                        gps   (1:ncols), sig   (1:kmax), sigml(1:(kmax + 1)), sigki (1:kmax), &
                        ustar (1:ncols), tkesfc(1:nCols, 1:kMax + 1), QSfc0(1:nCols, latco), thz0(1:ncols), &
                        qz0   (1:ncols), uz0   (1:ncols), vz0  (1:ncols), z0    (1:ncols), &
                        pblh  (1:ncols), ELM(1:ncols, 1:kMax), akhs (1:ncols), akms  (1:ncols), &
                        latco, tsfc0(1:nCols, latco), CT   (1:ncols), htdisp(1:ncols), &
                        PBLSFC_CoefKm(1:ncols, 1:kmax), PBLSFC_CoefKh(1:ncols, 1:kmax), gicem(1:ncols, 1:kmax), &
                        gliqm(1:ncols, 1:kmax), gvarm   (1:ncols, 1:kmax, :))
            ELSE
                CALL SfcPBL_Driver(&
                        gu    (1:nCols, 1:kMax), gv  (1:nCols, 1:kMax), gt(1:nCols, 1:kMax), gq(1:nCols, 1:kMax), &
                        delsig(1:kmax), ncols, kmax, deltm, &
                        colrad(1:nCols), tmsfc(1:nCols, 1:kMax, 1:3), qmsfc(1:nCols, 1:kMax, 1:5 + nClass + nAeros), umsfc(1:nCols, 1:kMax, 1:4), &
                        Mmlen (1:ncols), dt, jdt, &
                        topog (1:ncols), xland (1:ncols), lowlyr(1:ncols), &
                        gps   (1:ncols), sig   (1:kmax), sigml(1:(kmax + 1)), sigki (1:kmax), &
                        ustar (1:ncols), tkesfc(1:nCols, 1:kMax + 1), QSfc0(1:nCols, latco), thz0(1:ncols), &
                        qz0   (1:ncols), uz0   (1:ncols), vz0  (1:ncols), z0    (1:ncols), &
                        pblh  (1:ncols), ELM(1:ncols, 1:kMax), akhs (1:ncols), akms  (1:ncols), &
                        latco, tsfc0(1:nCols, latco), CT   (1:ncols), htdisp(1:ncols), &
                        PBLSFC_CoefKm(1:ncols, 1:kmax), PBLSFC_CoefKh(1:ncols, 1:kmax), gicem(1:ncols, 1:kmax), &
                        gliqm(1:ncols, 1:kmax))

            END IF
        ELSE
            CALL SfcPBL_Driver(&
                    gu    (1:nCols, 1:kMax), gv  (1:nCols, 1:kMax), gt(1:nCols, 1:kMax), gq(1:nCols, 1:kMax), &
                    delsig(1:kmax), ncols, kmax, deltm, &
                    colrad(1:nCols), tmsfc(1:nCols, 1:kMax, 1:3), qmsfc(1:nCols, 1:kMax, 1:5 + nClass + nAeros), umsfc(1:nCols, 1:kMax, 1:4), &
                    Mmlen (1:ncols), dt, jdt, &
                    topog (1:ncols), xland (1:ncols), lowlyr(1:ncols), &
                    gps   (1:ncols), sig   (1:kmax), sigml(1:(kmax + 1)), sigki (1:kmax), &
                    ustar (1:ncols), tkesfc(1:nCols, 1:kMax + 1), QSfc0(1:nCols, latco), thz0(1:ncols), &
                    qz0   (1:ncols), uz0   (1:ncols), vz0  (1:ncols), z0    (1:ncols), &
                    pblh  (1:ncols), ELM(1:ncols, 1:kMax), akhs (1:ncols), akms  (1:ncols), &
                    latco, tsfc0(1:nCols, latco), CT   (1:ncols), htdisp(1:ncols), &
                    PBLSFC_CoefKm(1:ncols, 1:kmax), PBLSFC_CoefKh(1:ncols, 1:kmax))
        ENDIF
        !PRINT*,' surface_driver'
        !
        !     surface parameterization
        !
        IF (microphys) THEN
            QCF = gicem
            QCL = gliqm
        ELSE
            QCF = 0.0_r8
            QCL = 0.0_r8
        END IF

        !Test to reduce sfc temperature (start)
        !print*, '1-sindelcos=',1-sindel*COS(colrad(i))
        ! print*, 'zenith=',zenith
        fator = 0.0_r8 !=0 (cancel correction) 1.5_r8
        pfator = 0.5_r8 !0.7
        gt(1:nCols, 1) = gt(1:nCols, 1) + fator * (1 - pfator) * (1 - zenith)  !*(1-sindel*COS(colrad(i)))
        gt(1:nCols, 2) = gt(1:nCols, 2) + fator * (pfator) * (1 - zenith)  !*(1-sindel*COS(colrad(i)))
        gq(1:nCols, 1) = gq(1:nCols, 1) - cp / hl * fator * (zenith)  !*(1-sindel*COS(colrad(i)))
        !  print*, 'gt,gq', cp, hl, gt  (1:nCols,1:3) ,hl/cp*gq (1:nCols,1:3)
        ! fator=0.5_r8
        ! gt(1:nCols,1) =  gt(1:nCols,1) - fator vdcc
        ! gq(1:nCols,1) =  gq(1:nCols,1) * fator
        !Test to reduce sfc temperature (end)

        !print*, 'inside PhysicsDriver shape LwSfcNet', shape(LwSfcNet)

        !colocando tsurf queria mejorar el valor impreso en los archivos de salida. Y queria ver si hay un impacto en los flujos
        !de calor latente y sensible. Pero no funciona, pues dentro de surface_driver esta alterando el valor prescrito.
        !Necesito modificar esta parte.
        IF (specSfc) THEN
            !tsurf(1:nCols)=sgTsfc
            !TSK(1:nCols)=sgTsfc
            !tseam(1:nCols)=sgTsfc
            !sens(1:nCols)=sgH
            !evap(1:nCols)=sgLE
        ENDIF

        print*, 'before surface_driver'
        CALL surface_driver(&
                jdt, latco, dtc3x, nmax, &
                nCols, kMax, intg, nsx(1:ncols), &
                gt  (1:nCols, 1:kMax), gq (1:nCols, 1:kMax), gu(1:nCols, 1:kMax), gv (1:nCols, 1:kMax), &
                gps  (1:nCols), cosz(1:nCols), itype(1:nCols), month2(1:nCols), &
                ssib(1:nCols), imask(1:nCols), cos2(1:nCols), LwSfcDown (1:nCols), &
                ssib_VisBeam(1:nCols), ssib_VisDiff(1:nCols), ssib_NirBeam(1:nCols), ssib_NirDiff(1:nCols), &
                zenith(1:nCols), xVisBeam(1:nCols), xVisDiff(1:nCols), xNirBeam(1:nCols), &
                xNirDiff (1:nCols), ppli (1:nCols), ppci  (1:nCols), tmsfc(1:nCols, 1:kMax, 1:3), &
                qmsfc(1:nCols, 1:kMax, 1:5), umsfc(1:nCols, 1:kMax, 1:4), tsea(1:nCols), slrad(1:nCols), &
                tsurf (1:nCols), qsurf(1:nCols), colrad(1:nCols), sigki(1), &
                delsig(1:kMax), sens  (1:nCols), evap (1:nCols), umom (1:nCols), &
                vmom (1:nCols), zorl (1:nCols), tseam(1:nCols), sice(1:ncols), &
                ustar(1:nCols), qsfc  (1:nCols), tsfc(1:nCols), z0  (1:ncols), &
                htdisp(1:nCols), temp2m(1:nCols), umes2m(1:nCols), mskant(1:nCols), &
                taux(1:nCols), tauy(1:nCols), tkesfc(1:nCols, 1:1), ndvi(1:nCols), &
                ndvim(1:nCols), tod, bstar(1:nCols), sflux_t(1:nCols), &
                sflux_r (1:nCols), sflux_u(1:nCols), sflux_v(1:nCols), r_aer(1:nCols), &
                snow(1:nCols), cldtot(1:nCols, 1:kMax), HML    (1:nCols), HUML(1:nCols), &
                HVML(1:nCols), TSK(1:nCols), z0sea(1:nCols), ySwSfcNet(1:nCols), &
                LwSfcNet(1:nCols), pblh(1:nCols), QCF(1:nCols, 1:kMax), QCL(1:nCols, 1:kMax), &
                sm0(1:nCols, 1:3), mlsi(1:nCols), LwSfcDown (1:nCols), month    (1:nCols), &
                Mmlen (1:nCols))  !sm0 mlsi add solange 13-11-2012
        print*, 'after surface_driver'

        ! print specSfc values
        IF (specSfc) THEN
            !!print*, 'specSfc=', sgTsfc, sgH, sgLE, sgTau
            !!print*, 'check_specSfc', tsurf (1:nCols)-273.15, tsfc(1:nCols), sens(1:nCols), evap(1:nCols), sens+evap
            !!enable lines below to allow specSfc effects
            !tsurf(1:nCols)=sgTsfc
            !TSK(1:nCols)=sgTsfc
            !tseam(1:nCols)=sgTsfc
            !sens(1:nCols)=sgH
            !evap(1:nCols)=sgLE
            !thz0(1:nCols)=sgTsfc
            !sens(1:nCols)=sgH
            !evap(1:nCols)=sgLE
        ENDIF

        !
        !     yamada-mellor pbl parameterization
        !     ( solving the matrices from bottom to top )
        !
        print*, 'before pbl_driver'
        IF (microphys) THEN
            IF((nClass + nAeros)>0 .and. PRESENT(gvarm))THEN
                DO m = 1, nClass + nAeros
                    DO k = 1, kmax
                        DO i = 1, ncols
                            gicem_local  (i, k) = gicem  (i, k)
                            gliqm_local  (i, k) = gliqm  (i, k)
                            gvarm_local  (i, k, m) = gvarm  (i, k, m)
                        END DO
                    END DO
                END DO
                !print*, "Denis ************* nCols, kmax, nClass = ", nCols, kmax, nClass
                !print*, "Denis ************* nAeros, latco = ", nAeros, latco

                CALL pbl_driver(&
                        ncols, &
                        kmax, &
                        gu (1:nCols, 1:kmax), &
                        gv (1:nCols, 1:kmax), &
                        gt (1:nCols, 1:kmax), &
                        gq (1:nCols, 1:kmax), &
                        delsig(1:kmax), &
                        deltm, &
                        colrad (1:nCols), &
                        tmtx (1:ncols, 1:kmax, 1:3), &
                        qmtx (1:ncols, 1:kmax, 1:5 + nClass + nAeros), &
                        umtx (1:ncols, 1:kmax, 1:4), &
                        tmsfc (1:ncols, 1:kmax, 1:3), &
                        qmsfc(1:ncols, 1:kmax, 1:5 + nClass + nAeros), &
                        umsfc (1:ncols, 1:kmax, 1:4), &
                        gl0 (1:ncols), &
                        2, &
                        sigki (1:kmax), &
                        gps (1:ncols), &
                        sig (1:kmax), &
                        sigml (1:kmax + 1), &
                        z0 (1:ncols), &
                        xland (1:ncols), &
                        lowlyr (1:ncols), &
                        sice (1:ncols), &
                        snow (1:nCols), &
                        sens (1:nCols), &
                        evap (1:nCols), &
                        thz0 (1:nCols), &
                        qz0 (1:ncols), &
                        uz0 (1:ncols), &
                        vz0 (1:ncols), &
                        tkemyj (1:nCols, 1:kMax + 1), &
                        ustar (1:nCols), &
                        akhs (1:ncols), &
                        akms (1:ncols), &
                        latco, &
                        QSfc0 (1:nCols, latco), &
                        tsfc0 (1:nCols, latco), &
                        topog (1:ncols), &
                        ct (1:ncols), &
                        LwCoolRate (1:ncols, 1:kmax), &
                        LwCoolRateC (1:ncols, 1:kmax), &
                        cldtot (1:ncols, 1:kmax), &
                        cldinv (1:ncols, 1:kmax), &
                        cldsat (1:ncols, 1:kmax), &
                        cldson (1:ncols, 1:kmax), &
                        qliq (1:ncols, 1:kMax), &
                        bstar (1:ncols), &
                        htdisp (1:ncols), &
                        taux (1:ncols), &
                        tauy (1:ncols), &
                        kt, &
                        jdt, &
                        PBL_CoefKm (1:ncols, 1:kmax + 1), &
                        PBL_CoefKh (1:ncols, 1:kmax + 1), &
                        PBLSFC_CoefKm (1:ncols, 1:kmax), &
                        PBLSFC_CoefKh (1:ncols, 1:kmax), &
                        pblh (1:ncols), &
                        tpert (1:ncols), &
                        qpert (1:ncols), &
                        tstar (1:ncols), &
                        wstar (1:ncols), &
                        var (1:ncols), &
                        tauresx (1:ncols), &
                        tauresy (1:ncols), &
                        gicem_local (1:ncols, 1:kmax), &
                        gliqm_local (1:ncols, 1:kmax), &
                        gvarm_local (1:ncols, 1:kmax, 1:8 + nAeros))
            ELSE
                DO k = 1, kmax
                    DO i = 1, ncols
                        gicem_local  (i, k) = gicem  (i, k)
                        gliqm_local  (i, k) = gliqm  (i, k)
                    END DO
                END DO
                CALL pbl_driver(&
                        ncols, kmax, gu    (1:nCols, 1:kmax), gv    (1:nCols, 1:kmax), gt     (1:nCols, 1:kmax), &
                        gq    (1:nCols, 1:kmax), delsig(1:kmax), &
                        deltm, colrad (1:nCols), &
                        tmtx  (1:ncols, 1:kmax, 1:3), qmtx (1:ncols, 1:kmax, 1:5 + nClass + nAeros), umtx   (1:ncols, 1:kmax, 1:4), &
                        tmsfc (1:ncols, 1:kmax, 1:3), qmsfc(1:ncols, 1:kmax, 1:5 + nClass + nAeros), umsfc  (1:ncols, 1:kmax, 1:4), &
                        gl0   (1:ncols), 2, sigki  (1:kmax), &
                        gps   (1:ncols), sig   (1:kmax), sigml  (1:kmax + 1), &
                        z0    (1:ncols), xland (1:ncols), lowlyr (1:ncols), &
                        sice  (1:ncols), snow  (1:nCols), sens   (1:nCols), &
                        evap  (1:nCols), thz0  (1:nCols), qz0    (1:ncols), &
                        uz0   (1:ncols), vz0   (1:ncols), tkemyj (1:nCols, 1:kMax + 1), &
                        ustar (1:nCols), akhs  (1:ncols), akms   (1:ncols), &
                        latco, QSfc0 (1:nCols, latco), tsfc0  (1:nCols, latco), &
                        topog (1:ncols), ct    (1:ncols), LwCoolRate(1:ncols, 1:kmax), &
                        LwCoolRateC(1:ncols, 1:kmax), cldtot(1:ncols, 1:kmax), cldinv(1:ncols, 1:kmax), &
                        cldsat(1:ncols, 1:kmax), cldson(1:ncols, 1:kmax), qliq   (1:ncols, 1:kMax), &
                        bstar (1:ncols), htdisp(1:ncols), taux   (1:ncols), &
                        tauy  (1:ncols), kt, jdt, &
                        PBL_CoefKm(1:ncols, 1:kmax + 1), PBL_CoefKh(1:ncols, 1:kmax + 1), &
                        PBLSFC_CoefKm(1:ncols, 1:kmax), PBLSFC_CoefKh(1:ncols, 1:kmax), pblh   (1:ncols), &
                        tpert (1:ncols), qpert (1:ncols), tstar  (1:ncols), &
                        wstar (1:ncols), var(1:ncols), tauresx(1:ncols), &
                        tauresy(1:ncols), gicem_local (1:ncols, 1:kmax), gliqm_local  (1:ncols, 1:kmax))

            END IF
        ELSE
            DO k = 1, kmax
                DO i = 1, ncols
                    gicem_local  (i, k) = 0.0_r8
                    gliqm_local  (i, k) = 0.0_r8
                END DO
            END DO
            CALL pbl_driver(&
                    ncols, kmax, &
                    gu    (1:nCols, 1:kmax), gv    (1:nCols, 1:kmax), gt     (1:nCols, 1:kmax), &
                    gq    (1:nCols, 1:kmax), delsig(1:kmax), deltm, colrad (1:nCols), &
                    tmtx  (1:ncols, 1:kmax, 1:3), qmtx (1:ncols, 1:kmax, 1:5 + nClass + nAeros), umtx   (1:ncols, 1:kmax, 1:4), &
                    tmsfc (1:ncols, 1:kmax, 1:3), qmsfc(1:ncols, 1:kmax, 1:5 + nClass + nAeros), umsfc  (1:ncols, 1:kmax, 1:4), &
                    gl0   (1:ncols), 2, sigki  (1:kmax), &
                    gps   (1:ncols), sig   (1:kmax), sigml  (1:kmax + 1), &
                    z0    (1:ncols), xland (1:ncols), lowlyr (1:ncols), &
                    sice  (1:ncols), snow  (1:nCols), sens   (1:nCols), &
                    evap  (1:nCols), thz0  (1:nCols), qz0    (1:ncols), &
                    uz0   (1:ncols), vz0   (1:ncols), tkemyj (1:nCols, 1:kMax + 1), &
                    ustar (1:nCols), akhs  (1:ncols), akms   (1:ncols), &
                    latco, QSfc0 (1:nCols, latco), tsfc0  (1:nCols, latco), &
                    topog (1:ncols), ct    (1:ncols), LwCoolRate(1:ncols, 1:kmax), &
                    LwCoolRateC(1:ncols, 1:kmax), cldtot(1:ncols, 1:kmax), cldinv(1:ncols, 1:kmax), &
                    cldsat(1:ncols, 1:kmax), cldson(1:ncols, 1:kmax), qliq   (1:ncols, 1:kMax), &
                    bstar (1:ncols), htdisp(1:ncols), taux   (1:ncols), &
                    tauy  (1:ncols), kt, jdt, &
                    PBL_CoefKm(1:ncols, 1:kmax + 1), PBL_CoefKh(1:ncols, 1:kmax + 1), &
                    PBLSFC_CoefKm(1:ncols, 1:kmax), PBLSFC_CoefKh(1:ncols, 1:kmax), pblh   (1:ncols), &
                    tpert (1:ncols), qpert (1:ncols), tstar  (1:ncols), &
                    wstar (1:ncols), var(1:ncols), tauresx(1:ncols), &
                    tauresy(1:ncols), gicem_local (1:ncols, 1:kmax), gliqm_local  (1:ncols, 1:kmax))
        ENDIF
        !END IF
        DO I = 1, nCols
            rho(i) = 287.04_r8 * gt(i, 1) / (100.0_r8 * gps(i))
            ustar(i) = MAX(SQRT(SQRT(taux(i)**2 + tauy(i)**2) * rho(i)), 0.000001_r8) ! [m/s]
        END DO
        !
        !     umtx(i,k,3) zonal wind tendency due to vertical diffusion
        !     umtx(i,k,4) meridional wind tendency due to vertical diffusion
        !
        DO k = 1, kmax
            DO i = 1, ncols
                gyu(i, k) = gyu(i, k) + umtx(i, k, 3)
                gyv(i, k) = gyv(i, k) + umtx(i, k, 4)
            END DO
        END DO
        !
        !     qmtx(i,k,3) Specific humidity tendency due to vertical diffusion
        !
        DO k = 1, kmax
            DO i = 1, ncols
                gqd(i, k) = gqd(i, k) + qmtx(i, k, 3)
            END DO
        END DO
        !
        !     qmtx(i,k,4 and 5) ice and liquid water tendencies
        !                       due to vertical diffusion (if microphys = true)
        !
        IF (microphys) THEN
            DO k = 1, kmax
                DO i = 1, ncols
                    gicet(i, k) = gicet(i, k) + qmtx(i, k, 4)
                    gliqt(i, k) = gliqt(i, k) + qmtx(i, k, 5)
                END DO
            END DO
            IF((nClass + nAeros)>0 .and. PRESENT(gvarm))THEN
                DO kk = 1, nClass + nAeros
                    DO k = 1, kmax
                        DO i = 1, ncols
                            gvart(i, k, kk) = gvart(i, k, kk) + qmtx(i, k, 6 + kk - 1)
                        END DO
                    END DO
                END DO
            END IF
        ENDIF
        !
        !     tmtx(i,k,3) virtual temperature tendency due to vertical diffusion
        !
        DO k = 1, kmax
            DO i = 1, ncols
                ! Enver: Test with diffusion near surface (start)
                IF (k < 7) THEN
                    tmtx(i, k, 3) = (1.0_r8 + 0.608_r8 * gq(i, k)) * tmtx(i, k, 3) + &
                            0.608_r8 * gt(i, k) * qmtx(i, k, 3)
                    !tmtx(i,k,3) = 1.5_r8 * tmtx(i,k,3) !increasing vertical diffusion X%
                    tmtx(i, k, 3) = 0.5_r8 * tmtx(i, k, 3) !decresing vertical diffusion X%
                ELSEIF (k >= 7 .and. k < 10) THEN
                    tmtx(i, k, 3) = (1.0_r8 + 0.608_r8 * gq(i, k)) * tmtx(i, k, 3) + &
                            0.608_r8 * gt(i, k) * qmtx(i, k, 3)
                    !tmtx(i,k,3) = 1.5_r8 * tmtx(i,k,3) !increasing vertical diffusion X%
                    tmtx(i, k, 3) = 0.5_r8 * tmtx(i, k, 3) !decresing vertical diffusion X%
                ELSE
                    tmtx(i, k, 3) = (1.0_r8 + 0.608_r8 * gq(i, k)) * tmtx(i, k, 3) + &
                            0.608_r8 * gt(i, k) * qmtx(i, k, 3)
                ENDIF
                !Below are the older lines uncomment
                !tmtx(i,k,3) = (1.0_r8+ 0.608_r8 * gq(i,k)) * tmtx(i,k,3) + &
                !     0.608_r8 * gt(i,k) * qmtx(i,k,3)
                ! Enver: Test with diffusion near surface (end)
                gtd(i, k) = gtd(i, k) + tmtx(i, k, 3)
                gtd(i, k) = gtd(i, k) + (1.0_r8 + 0.608_r8 * gq(i, k)) * (ySwHeatRate(i, k) + LwCoolRate(i, k))
            END DO
        END DO

        IF((kt.NE.0) .OR. (jdt.NE.1)) THEN
            IF (cdhl(jdt)) THEN
                DO i = 1, ncols
                    ustr(i) = -umom(i)
                    vstr(i) = -vmom(i)
                END DO
            END IF
        END IF

        DO k = 1, kmax
            DO i = 1, ncols
                gt(i, k) = gt(i, k) * (1.0e0_r8 + 0.608e0_r8 * gq(i, k))
            END DO
        END DO
    END SUBROUTINE physcs


    SUBROUTINE InitSimpPhys(fgtmp, tequi, sl, dt)
        REAL(KIND = r8), INTENT(IN) :: fgtmp(:, :, :)
        REAL(KIND = r8), INTENT(INOUT) :: tequi(:, :, :)
        REAL(KIND = r8), INTENT(IN) :: sl(:)
        REAL(KIND = r8), INTENT(IN) :: dt
        INTEGER :: kMax, ibMax, jbMax
        INTEGER :: ib, jb, k
        INTEGER :: NstepSimpPhys = 1000
        REAL(KIND = r8) :: pi
        REAL(KIND = r8) :: beta
        REAL(KIND = r8) :: kappa
        REAL(KIND = r8) :: tf
        ibMax = SIZE(fgtmp, 1)
        jbMax = SIZE(fgtmp, 3)
        kMax = SIZE(sl)
        ALLOCATE(teq(ibMax, kMax, jbMax))
        ALLOCATE(tauri(kMax))
        ALLOCATE(alfa(kMax))
        pi = 4.0_r8 * ATAN(1.0_r8)
        !PKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPK
        ! DO jb = 1, jbMax
        !    DO ib = 1, ibMax
        !         tesurf(ib,jb) = tsfc0-tdelt*(SIN(.5*pi-colrad2D(ib,jb)))**2
        !    END DO
        ! END DO
        ! DO k = 1, kMax
        !    DO jb = 1, jbMax
        !       DO ib = 1, ibMax
        !       teq(ib,k,jb) = MAX(tstrat, fgtmp(ib,1,jb) + h0*rlaps*LOG(sl(k)))
        !       END DO
        !    END DO
        ! END DO
        !PKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPK
        !
        !     radiative damping inverse time
        !
        DO k = 1, kmax
            tauri(k) = tdampr(k) * 86400.0_r8
            tauri(k) = 1.0_r8 / tauri(k)
        END DO
        !
        !     frictional damping inverse time
        !
        DO k = 1, kmax
            alfa(k) = tdampf(k) * 86400.0_r8
            alfa(k) = 1.0_r8 / alfa(k)
        END DO

        IF(TRIM(start) == "cold") THEN
            beta = (1.0_r8 / 365.25_r8)          ! days/days
            kappa = 1.0_r8 / tdampr(1)        ! days^-1
            tf = (dt * REAL(NstepSimpPhys, r8)) / 86400.0_r8
            DO k = 1, kMax
                DO jb = 1, jbMax
                    DO ib = 1, ibMax
                        teq (ib, k, jb) = fgtmp(ib, k, jb) * (exp(-beta) * (1.0_r8 - exp(-kappa * tf)) + exp(-kappa * tf))
                        tequi(ib, k, jb) = teq (ib, k, jb)
                    END DO
                END DO
            END DO
        ELSE
            DO k = 1, kMax
                DO jb = 1, jbMax
                    DO ib = 1, ibMax
                        teq (ib, k, jb) = tequi(ib, k, jb)
                    END DO
                END DO
            END DO
        END IF
    END SUBROUTINE InitSimpPhys


    SUBROUTINE SimpPhys(gu, gv, gtmp, gyv, gyu, gtd, ibMax, ibLim, kMax, jb, jbMax, ibMaxPerJB, &
            iPerIJB, jPerIJB)
        !
        ! simfiz :simplified physics package
        !            newtonian cooling
        !            raleigh damping
        !
        ! arguments
        !     gu:     u*cos(latitude)
        !     gv:     v*cos(latitude)
        !     gtmp:   absolute temperature
        !     gyv:    (du/dt) tendency of zonal wind * cos(latitude)
        !     gyu:    -(dv/dt) negative of tendency of v*cos(latitude)
        !
        INTEGER, INTENT(IN) :: ibMax
        INTEGER, INTENT(IN) :: jbMax
        INTEGER, INTENT(IN) :: ibLim
        INTEGER, INTENT(IN) :: kMax
        INTEGER, INTENT(IN) :: jb
        REAL(KIND = r8), INTENT(IN) :: gu(ibMax, kMax)
        REAL(KIND = r8), INTENT(IN) :: gv(ibMax, kMax)
        REAL(KIND = r8), INTENT(IN) :: gtmp(ibMax, kMax)
        REAL(KIND = r8), INTENT(INOUT) :: gyv(ibMax, kMax)
        REAL(KIND = r8), INTENT(INOUT) :: gyu(ibMax, kMax)
        REAL(KIND = r8), INTENT(INOUT) :: gtd(ibMax, kmax)
        INTEGER, INTENT(IN) :: ibMaxPerJB (jbMax)
        INTEGER, INTENT(IN) :: iPerIJB    (ibMax, jbMax)
        INTEGER, INTENT(IN) :: jPerIJB    (ibMax, jbMax)
        INTEGER :: ib, k, i, j
        LOGICAL :: simfiz = .TRUE.
        REAL(KIND = r8), PARAMETER :: tendency_Temp(1:28) = (/ 1.0, 2.0, 3.0, 4.0, 5.0, &
                6.0, 7.0, 8.0, 9.0, 10.0, &
                11.0, 12.0, 13.0, 14.0, 15.0, &
                16.0, 17.0, 18.0, 19.0, 20.0, &
                21.0, 22.0, 23.0, 24.0, 25.0, &
                26.0, 27.0, 28.0/)
        INTEGER, PARAMETER :: IlocLon = 96 ! grid Position (longitude --> grads)
        INTEGER, PARAMETER :: jlocLat = 43 ! grid Position (latitude  --> grads)
        IF (simfiz) THEN
            !
            !     radiative damping
            !
            DO k = 1, kmax
                DO ib = 1, ibLim
                    gtd(ib, k) = gtd(ib, k) + tauri(k) * (teq(ib, k, jb) - (tov(k) + gtmp(ib, k)))
                END DO
            END DO
            !
            !     complete non-linear part of temperature tendency
            !     ------------------------------------------------

            DO k = 1, kmax
                DO ib = 1, ibMaxPerJB(jb)
                    j = jPerIJB(ib, jb)
                    i = iPerIJB(ib, jb)
                    IF(i == ilocLon .AND. j== jlocLat)THEN
                        gtd(ib, k) = gtd(ib, k) + tendency_Temp(k)
                    END IF
                END DO
            END DO
            !
            !     raleigh friction: zonal momentum equation
            !
            DO k = 1, kmax
                DO ib = 1, ibLim
                    gyu(ib, k) = gyu(ib, k) - alfa(k) * gu(ib, k)
                END DO
            END DO
            !
            !     raleigh friction: meridional momentum equation
            !
            DO k = 1, kmax
                DO ib = 1, ibLim
                    gyv(ib, k) = gyv(ib, k) - alfa(k) * gv(ib, k)
                END DO
            END DO

        END IF

    END SUBROUTINE SimpPhys
END MODULE PhysicsDriver
