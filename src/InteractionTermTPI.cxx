#include "InteractionTermTPI.hxx"

namespace JJCorrFitter
{
    InteractionTermTPI::InteractionTermTPI(SpinState state) :
    m_spinState(state),
    fD0s(0.,0.), 
    fF0s(0.,0.), 
    fD0t(0.,0.), 
    fF0t(0.,0.),
    fPionac(0),
    fOneoveracsq(0),
    fTwopioverac(0),
    fCoulqscpart(0),
    fEuler(0),
    fF0(0),
    fD0(0),
    fTwospin(0),
    fWritegrps(0),
    fPcount(0),
    fCoulombSteps(170),
    fRStarOutS(0),
    fRStarSideS(0),
    fRStarLongS(0),
    fRStarS(0),
    fRStarOut(0),
    fRStarSide(0),
    fRStarLong(0),
    fRStar(0),
    fKStarOut(0),
    fKStarSide(0),
    fKStarLong(0),
    fKStar(0) 
    {
        switch (state)
        {
            case SpinState::SpinAveraged:
                m_numberOfParams = 2;
                break;

            case SpinState::SpinSeparated:
                m_numberOfParams = 4;
                break;

            case SpinState::None:
                m_numberOfParams = 0;
                break;

            default:
                m_numberOfParams = 0;
                break;
        }

        InitializeGamow();
        m_InteractionTermName = "p-p Lednicky-Lyuboshitz (TPI)";
    }

    constexpr void InteractionTermTPI::InitializeGamow() 
    {
        using namespace std::complex_literals;

        fEuler = 0.577215665;
        fF0    = 7.77 / m_gevToFm;
        fD0    = 2.77 / m_gevToFm;

        fPionac  = 57.63975274 / m_gevToFm;
        fTwospin = 1;

        fD0t = {1.7 / m_gevToFm, 0.};
        fF0t = {-5.4 / m_gevToFm, 0.};

        // ESC08
        fF0s = {7.771 / m_gevToFm, 0.};
        fD0s = {2.754 / m_gevToFm, 0.};

        fOneoveracsq = 1.0 / (fPionac * fPionac);
        fTwopioverac = 2.0 * m_pi / fPionac;
    }

    void InteractionTermTPI::SetParameters(const std::vector<double> &pars)
    {
        if (pars.size() != m_numberOfParams)
        {
            throw std::length_error("InteractionTermTPI::SetParameters - Provided number of parameters does not match the expected number: " + std::to_string(m_numberOfParams));
        }

        switch (m_spinState)
        {
            case SpinState::SpinAveraged:
                fF0 = pars.at(0) / m_gevToFm;
                fD0 = pars.at(1) / m_gevToFm;
                break;

            case SpinState::SpinSeparated:
                fF0s = pars.at(0) / m_gevToFm;
                fD0s = pars.at(1) / m_gevToFm;
                fF0t = pars.at(2) / m_gevToFm;
                fD0t = pars.at(3) / m_gevToFm;
                break;
        
            default:
                break;
        } 
    }

    void InteractionTermTPI::SetMomentum(float kStar)
    {
        fKStar = kStar * m_mevToGev;
    }
    
    double InteractionTermTPI::GetValue(float rStar, float cosTheta)
    {
        return GetQuantumCoulombStrong(rStar / m_gevToFm,cosTheta);
    }

    constexpr double InteractionTermTPI::Gamow(double arg) const 
    {
        long double eta = fTwopioverac / arg;
        return (eta) *1.0 / (exp(eta) - 1.0);
    }

    constexpr double InteractionTermTPI::GetQuantumCoulombStrong(float rStar, float cosTheta) 
    {
        if (rStar < 0.0000000001)
            return 1.0;

        if (rStar < 1.0 / m_gevToFm) 
            return GetQuantumCoulomb(rStar,cosTheta);

        double tKstRst    = fKStar * rStar * cosTheta;
        long double kstar = fabs(fKStar);
        long double rho   = rStar * kstar;

        int ccase         = 0;
        int wavesign      = 1;

        // Classical limit - in the case of large k* we go to
        // classical coulomb interaction
        long double testp = rho * (1.0 + tKstRst / (rho));
        long double testm = rho * (1.0 - tKstRst / (rho));

        std::complex<long double> ffplus, ffminus;
        if ((testp > 15.0) && (testm > 15.0)) 
        {
            double asym = (1.0 - 1.0 / (rStar * (1.0 - tKstRst / rho) * fPionac * kstar * kstar)) / Gamow(kstar);
            asym = sqrt(asym);
            if (asym < 1.0)
                ffminus.real(1.0 + (asym - 1.0) * 2.0);
            else
                ffminus.real(1.0 + (asym - 1.0) / 2.0);

            ffminus.imag(sqrt(asym * asym - ffminus.real() * ffminus.real()));

            asym = (1.0 - 1.0 / (rStar * (1.0 + tKstRst / rho) * fPionac * kstar * kstar)) / Gamow(kstar);
            asym = sqrt(asym);
            if (asym < 1.0)
                ffplus.real(1.0 + (asym - 1.0) * 2.0);
            else
                ffplus.real(1.0 + (asym - 1.0) / 2.0);

            ffplus.imag(sqrt(asym * asym - ffplus.real() * ffplus.real()));
        }
        // Check for the classical limit in both functions separately
        else if (((testp < 15.0) && (testm < 15.0)))  // ||
        {
            // Calculate the F function
            GetFFdouble(rStar,cosTheta,ffplus, ffminus);
            ccase = 1;
        } 
        else if (testp < 15.0) 
        {
            GetFFsingle(rStar,cosTheta,ffplus, 1);
            GetFFsingle(rStar,cosTheta,ffminus, -1);
            if ((fabs(ffminus.real()) > 2.0) || fabs(ffminus.imag()) > 2.0) 
            {
                double asym = (1.0 - 1.0 / (rStar * (1.0 - tKstRst / (rho) *fPionac * kstar * kstar))) / Gamow(kstar);
                asym = sqrt(asym);
                if (asym < 1.0)
                    ffminus.real(1.0 + (asym - 1.0) * 2.0);
                else
                    ffminus.real(1.0 + (asym - 1.0) / 2.0);

                ffminus.imag(sqrt(asym * asym - ffminus.real() * ffminus.real()));
            }
            ccase = 2;
        }
        else 
        {
            GetFFsingle(rStar,cosTheta,ffminus, -1);
            GetFFsingle(rStar,cosTheta,ffplus, 1);
            if ((fabs(ffplus.real()) > 2.0) || fabs(ffplus.imag()) > 2.0) 
            {
                double asym = (1.0 - 1.0 / (rStar * (1.0 + tKstRst / (rho) *fPionac * kstar * kstar))) / Gamow(kstar);
                asym = sqrt(asym);
                if (asym < 1.0)
                    ffplus.real(1.0 + (asym - 1.0) * 2.0);
                else
                    ffplus.real(1.0 + (asym - 1.0) / 2.0);

                ffplus.imag(sqrt(asym * asym - ffplus.real() * ffplus.real()));
            }
            ccase = 3;
        }

        long double eta = 1.0 / (kstar * fPionac);
        long double hfun = GetH(eta);
        std::complex<long double> gtilde = GetG(eta, rho, hfun);
        std::complex<long double> gtilor = gtilde / static_cast<long double>(rStar);

        std::complex<long double> fcouls, fcoult;
        Getfc(kstar, eta, hfun, fcouls, fcoult);

        std::complex<long double> fgs = gtilor * fcouls;
        long double fgmods = norm(fgs);

        std::complex<long double> fgt = gtilor * fcoult;
        long double fgmodt = norm(fgt);

        std::complex<long double> expikr(cos(tKstRst),-sin(tKstRst));

        std::complex<long double> expikrc  = conj(expikr);
        std::complex<long double> ffplusc  = conj(ffplus);
        std::complex<long double> ffminusc = conj(ffminus);

        std::complex<long double> expikr2  = pow(expikr, 2);
        std::complex<long double> expikrc2 = conj(expikr2);
        std::complex<long double> sterm    = expikr2 * ffplus * ffminusc;
        std::complex<long double> tterm    = expikrc2 * ffminus * ffplusc;

        std::complex<long double> epfpc = expikr * ffplus;
        std::complex<long double> emfmc = expikrc * ffminus;

        long double fcgefhs = (fgs.real() * emfmc.real() + fgs.imag() * emfmc.imag());
        long double fcgefgs = (fgs.real() * epfpc.real() + fgs.imag() * epfpc.imag());

        long double fcgefht = (fgt.real() * emfmc.real() + fgt.imag() * emfmc.imag());
        long double fcgefgt = (fgt.real() * epfpc.real() + fgt.imag() * epfpc.imag());

        long double smult = 1 + wavesign;

        if (!finite(ffminus.real())) 
            std::cout << "FFMinus Re not a number ! " << testp << " " << testm << " " << ccase << std::endl;

        if (!finite(ffminus.imag())) 
            std::cout << "FFMinus Im not a number !" << testp << " " << testm << " " << ccase << std::endl;

        if ((ffplus.real() > 2.0) || (ffplus.real() < -2.0))
            std::cout << std::endl << "FFplus Re wild !" << ffplus.real() << " case " << ccase << " " << testp << " " << testm << std::endl;

        if ((ffplus.imag() > 2.0) || (ffplus.imag() < -2.0))
            std::cout << "FFplus Im wild !" << ffplus.imag() << " case " << ccase << " " << testp << " " << testm << std::endl;

        if ((ffminus.real() > 2.0) || (ffminus.real() < -2.0))
            std::cout << std::endl << "FFminus Re wild !" << ffminus.real() << " case " << ccase << std::endl;

        if ((ffminus.imag() > 2.0) || (ffminus.imag() < -2.0))
            std::cout << "FFminus Im wild !" << ffminus.imag() << " case " << ccase << std::endl;

        if (fTwospin == 1) 
        {
            wavesign            = 1;
            smult               = 2;
            long double singlet = (0.5 * Gamow(kstar) * (2.0 * fgmods * smult + norm(ffplus) + norm(ffminus) + wavesign * sterm.real() + wavesign * tterm.real() + smult * 2 * (fcgefhs + fcgefgs)));
            
            wavesign            = -1;
            smult               = 0;
            long double triplet = (0.5 * Gamow(kstar) * (2.0 * fgmodt * smult + norm(ffplus) + norm(ffminus) + wavesign * sterm.real() + wavesign * tterm.real() + smult * 2 * (fcgefht + fcgefgt)));

            return (0.25 * singlet + 0.75 * triplet);
        } 
        else
        {
            return (0.5 * Gamow(kstar) * (2.0 * fgmods * smult + norm(ffplus) + norm(ffminus) + wavesign * sterm.real() + wavesign * tterm.real() + smult * 2 * (fcgefhs + fcgefgs)));
        }
        
    }

    constexpr std::complex<long double> InteractionTermTPI::GetG(long double eta, long double rho, long double hfun) const 
    {
        std::complex<long double> gtemp;
        
        std::pair<long double,long double> tmp = Bfunpfun(eta, rho);
        long double bres = tmp.first, pres = tmp.second;

        long double bmult = 2.0 * eta * rho * bres;

        gtemp = {pres + bmult * (log(fabs(2.0 * eta * rho)) + 2.0 * fEuler - 1.0 + hfun), bmult * Chiim(eta)};

        return gtemp;
    }

    constexpr void InteractionTermTPI::Getfc(long double kstar, long double eta, long double hfun, std::complex<long double> &fcs, std::complex<long double> &fct) const 
    {
        std::complex<long double> ci;
        std::complex<long double> cia;

        std::complex<long double> fis;
        std::complex<long double> dis;
        std::complex<long double> fcinvs;

        std::complex<long double> fit;
        std::complex<long double> dit;
        std::complex<long double> fcinvt;

        ci = {hfun * 2.0 / fPionac,Chiim(eta)};
        //cia   = ci * static_cast<long double>(2.0 / fPionac);

        fis       = static_cast<long double>(1.) / fF0s;
        dis       = fD0s * static_cast<long double>(0.5 * kstar * kstar);
        fcinvs = fis + dis - cia;
        fcs       = static_cast<long double>(1.) / fcinvs;

        fit       = static_cast<long double>(1.) / fF0t;
        dit       = fD0t * static_cast<long double>(0.5 * kstar * kstar);
        fcinvt = fit + dit - cia;
        fct       = static_cast<long double>(1.) / fcinvt;
    }

    constexpr std::pair<long double,long double> InteractionTermTPI::Bfunpfun(long double eta, long double rho) const 
    {
        long double b0   = 1;
        long double b1   = eta * rho;
        long double bsum = b0 + b1;
        long double bnpu = 0;
        long double p0   = 1.0;
        long double p1   = 0.0;
        long double psum = p0 + p1;
        long double pnpu = 0;

        if (rho > m_pi * 4.0) 
        {
            return std::make_pair(sin(rho) / rho,cos(rho));
        }

        long double bn   = b1;
        long double bnmu = b0;
        long double pn   = p1;
        long double pnmu = p0;
        for (int iter = 1; iter < 100000; iter++) 
        {
            bnpu = 2 * eta * rho * bn - rho * rho * bnmu;
            bnpu /= (1.0 * iter + 1.0) * (1.0 * iter + 2.0);
            bsum += bnpu;

            pnpu = 2 * eta * rho * pn - rho * rho * pnmu - (2.0 * iter + 1.0) * 2.0 * eta * rho * bn;
            pnpu /= (1.0 * iter) * (1.0 * iter + 1.0);
            psum += pnpu;

            bnmu = bn;
            bn   = bnpu;

            pnmu = pn;
            pn   = pnpu;
            if ((fabs(pnmu) + fabs(bnmu) + fabs(bnpu) + fabs(pnpu)) < 1.0e-20) 
            {
                break;
            }
        }

        return std::make_pair(bsum,psum);
    }

    constexpr long double InteractionTermTPI::GetH(long double eta) const 
    {
        long double etasum = log(1.0 / eta) - fEuler;
        long double series = 0.0;
        long double x2inv  = (eta * eta);
        long double element;
        long double save;
        for (int iter = 1; iter < 1000000; iter++) 
        {
            element = ((1.0 * iter) * (1.0 * iter) + x2inv) * (1.0 * iter);
            element = 1.0 / element;
            if (iter == 1) 
                save = 1.0e-10 * element;

            series += element;
            if (element < save) 
                break;
        }
        series *= x2inv;
        etasum += series;

        return etasum;
    }

    constexpr double InteractionTermTPI::GetQuantumCoulomb(float rStar, float cosTheta) 
    {
        if (rStar < 0.0000000001) 
            return 1.0;

        double kstrst = fKStar * rStar * cosTheta;
        int ccase     = 0;
        int wavesign  = 1;

        if (fTwospin == 1) 
        {
            if (fPcount == 3)
                wavesign = 1;
            else
                wavesign = -1;
            fPcount++;
            if (fPcount == 4) 
                fPcount = 0;
        }

        // Classical limit - if distance is larger than Coulomb radius,
        // the interaction does not matter
        if (fabs(rStar) > fabs(fPionac)) 
            return (1.0 + wavesign * cos(2 * kstrst));

        // Classical limit - in the case of large k* we go to
        // classical coulomb interaction
        long double testp = fabs(fKStar) * rStar * (1.0 + kstrst / (rStar * fabs(fKStar)));
        long double testm = fabs(fKStar) * rStar * (1.0 - kstrst / (rStar * fabs(fKStar)));

        if ((testp > 15.0) && (testm > 15.0)) 
        {
            double fasymplus  = (1.0 - 1.0 / ((rStar + kstrst) * fPionac * fKStar * fKStar));
            double fasymminus = (1.0 - 1.0 / ((rStar - kstrst) * fPionac * fKStar * fKStar));
            return 0.5 * ((fasymplus + fasymminus) * cos(2 * kstrst) + (2.0 * sqrt(fasymplus * fasymminus)));
        }

        std::complex<long double> ffplus, ffminus;
        // Check for the classical limit in both functions separately
        if (((testp < 15.0) && (testm < 15.0)))
        {
            // Calculate the F function
            GetFFdouble(rStar,cosTheta,ffplus, ffminus);
            ccase = 1;
        } 
        else if (testp < 15.0) 
        {
            double asym;
            GetFFsingle(rStar,cosTheta,ffplus);
            asym =
                (1.0 - 1.0 / (rStar * (1.0 - kstrst / (rStar * fabs(fKStar)) * fPionac * fKStar * fKStar))) / Gamow(fabs(fKStar));
            asym = sqrt(asym);
            if (asym < 1.0)
                ffminus.real(1.0 + (asym - 1.0) * 2.0);
            else
                ffminus.real(1.0 + (asym - 1.0) / 2.0);

            ffminus.imag(sqrt(asym * asym - ffminus.real() * ffminus.real()));
            ccase      = 2;
        } 
        else 
        {
            double asym;
            GetFFsingle(rStar,cosTheta,ffminus, -1);
            asym = (1.0 - 1.0 / (rStar * (1.0 + kstrst / (rStar * fabs(fKStar)) * fPionac * fKStar * fKStar))) / Gamow(fabs(fKStar));
            asym = sqrt(asym);
            if (asym < 1.0)
                ffplus.real(1.0 + (asym - 1.0) * 2.0);
            else
                ffplus.real(1.0 + (asym - 1.0) / 2.0);

            ffplus.imag(sqrt(asym * asym - ffplus.real() * ffplus.real()));
            ccase = 3;
        }

        std::complex<long double> expikr(cos(kstrst),sin(kstrst));

        std::complex<long double> ffplusc = conj(ffplus);
        std::complex<long double> ffminusc = conj(ffminus);

        std::complex<long double> expikr2 = expikr * expikr;
        std::complex<long double> expikrc2 = conj(expikr2);
        std::complex<long double> sterm = expikr2 * ffplus * ffminusc;
        std::complex<long double> tterm = expikrc2 * ffminus * ffplusc;


        if (!finite(ffminus.real()))
        std::cout << "FFMinus Re not a number !" << " " << ccase << std::endl;

        if (!finite(ffminus.imag()))
        std::cout << "FFMinus Im not a number !" << " " << ccase << std::endl;

        if ((ffplus.real() > 2.0) || (ffplus.real() < -2.0)) 
            std::cout << std::endl << "FFplus Re wild !" << ffplus.real() << std::endl;

        if ((ffplus.imag() > 2.0) || (ffplus.imag() < -2.0)) 
            std::cout << "FFplus Im wild !" << ffplus.imag() << std::endl;

        if ((ffminus.real() > 2.0) || (ffminus.real() < -2.0))
            std::cout << std::endl << "FFminus Re wild !" << ffminus.real() << " " << ccase << std::endl;

        if ((ffminus.imag() > 2.0) || (ffminus.imag() < -2.0)) 
            std::cout << "FFminus Im wild !" << ffminus.imag() << " " << ccase << std::endl;

        fCoulqscpart = 0.5 * Gamow(fabs(fKStar)) * (norm(ffplus) + norm(ffminus));

        return (0.5 * Gamow(fabs(fKStar)) * (norm(ffplus) + wavesign * sterm.real() + wavesign * tterm.real() + norm(ffminus)));
    }

    constexpr void InteractionTermTPI::GetFFdouble(float rStar, float cosTheta, std::complex<long double> &ffp, std::complex<long double> &ffm) const 
    {
        std::vector<long double> comprep(fCoulombSteps,0);
        std::vector<long double> compimp(fCoulombSteps,0);
        std::vector<long double> comprem(fCoulombSteps,0);
        std::vector<long double> compimm(fCoulombSteps,0);
        long double eta, ksip, ksim;
        std::complex<long double> alfa, zetp, zetm;

        int nsteps;

        long double kstar  = fabs(fKStar);

        if ((kstar * rStar * (1 + cosTheta) < 5.0) && (kstar * rStar * (1 - cosTheta) < 5.0))
            nsteps = 25;
        else if ((kstar * rStar * (1 + cosTheta) < 10.0) && (kstar * rStar * (1 - cosTheta) < 10.0))
            nsteps = 45;
        else if ((kstar * rStar * (1 + cosTheta) < 15.0) && (kstar * rStar * (1 - cosTheta) < 15.0))
            nsteps = 110;
        else
            nsteps = 150;

        eta = 1.0 / (kstar * fPionac);
        alfa = {0.0,-eta};

        std::complex<long double> fcomp, scompp, scompm;
        long double tcomp;
        std::complex<long double> sump, summ;
        std::complex<long double> fcmult;

        long double rad = rStar;

        ksip = kstar * rad * (1 + cosTheta);
        ksim = kstar * rad * (1 - cosTheta);

        zetp = {0.0,ksip};

        zetm = {0.0,ksim};

        fcomp = 1.0;
        scompp = 1.0;
        scompm = 1.0;
        tcomp = 1.0;

        for (int istep = 0; istep < nsteps; istep++) 
        {
            sump = fcomp * scompp;
            summ = fcomp * scompm;

            sump = sump / (tcomp * tcomp);
            summ = summ / (tcomp * tcomp);


            if (istep == 0) 
            {
                comprep[istep] = sump.real();
                compimp[istep] = sump.imag();

                comprem[istep] = summ.real();
                compimm[istep] = summ.imag();
            } 
            else 
            {
                comprep[istep] = comprep[istep - 1] + sump.real();
                compimp[istep] = compimp[istep - 1] + sump.imag();

                comprem[istep] = comprem[istep - 1] + summ.real();
                compimm[istep] = compimm[istep - 1] + summ.imag();
            }

            fcmult = {alfa.real() + istep, alfa.imag()};

            fcomp  = fcomp * fcmult;
            scompp = scompp * zetp;
            scompm = scompm * zetm;
            tcomp *= (istep + 1);
        }

        ffp = {comprep[nsteps - 1], compimp[nsteps - 1]};
        ffm = {comprem[nsteps - 1], compimm[nsteps - 1]};
    }

    constexpr void InteractionTermTPI::GetFFsingle(float rStar, float cosTheta, std::complex<long double> &ffp, int sign) const 
    {
        std::vector<double> comprep(fCoulombSteps,0);
        std::vector<double> compimp(fCoulombSteps,0);
        double eta, ksip;
        std::complex<long double> alfa, zetp;

        int nsteps;

        double kstar  = fabs(fKStar);

        if (kstar * rStar * (1 + (sign*cosTheta)) > 15.0)
            nsteps = 170;
        else if (kstar * rStar * (1 + (sign*cosTheta)) > 10.0)
            nsteps = 45;
        else if (kstar * rStar * (1 + (sign*cosTheta)) > 5.0)
            nsteps = 35;
        else
            nsteps = 15;

        eta     = 1.0 / (kstar * fPionac);
        alfa = {0.0,-eta};

        std::complex<long double> fcomp, scompp;
        long double tcomp;
        std::complex<long double> sump;
        std::complex<long double> fcmult;

        double rad = rStar;

        ksip = kstar * rad * (1 + (sign*cosTheta));

        zetp = {0.0,ksip};

        fcomp = 1.0;
        scompp = 1.0;
        tcomp = 1.0;

        for (int istep = 0; istep < nsteps; istep++) 
        {
            sump = fcomp * scompp;

            sump = sump / (tcomp * tcomp);

            if (istep == 0) 
            {
                comprep[istep] = sump.real();
                compimp[istep] = sump.imag();
            } 
            else 
            {
                comprep[istep] = comprep[istep - 1] + sump.real();
                compimp[istep] = compimp[istep - 1] + sump.imag();
            }

            fcmult = {alfa.real() + istep,alfa.imag()};

            fcomp  = fcomp * fcmult;
            scompp = scompp * zetp;
            tcomp *= (istep + 1);

            if ((sump.real() * sump.real() + sump.imag() * sump.imag()) < 1.0e-14) // another arbitrary small number (at lest be consistent for God's sake) - JJ
            {
                nsteps = istep;
                break;
            }
        }

        ffp = {comprep[nsteps - 1],compimp[nsteps - 1]};
    }
} // namespace JJCorrFItter
