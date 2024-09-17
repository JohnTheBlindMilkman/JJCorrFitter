#include "Fitter1D.hxx"

namespace JJCorrFitter
{
    explicit Fitter1D::Fitter1D(std::unique_ptr<TH1> &&data) : 
        m_parameters{
            {Fitter1D::Parameter::Norm,{0,"N",0.1}},
            {Fitter1D::Parameter::Lambda,{1,"Lambda",0.1}},
            {Fitter1D::Parameter::Radius,{2,"Radius",0.1}}
        },
        m_corrFunction(
            std::string(m_dataToFit->GetName()) + "_fit",
            "Fit to histogram",
            m_dataToFit->GetXaxis()->GetFirst(),
            m_dataToFit->GetXaxis()->GetLast(),
            m_dataToFit->GetNbinsX()),
        m_minimiser(ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad")),
        m_dataToFit(std::move(data))
    {
        if (m_minimiser == nullptr)
            throw std::runtime_error("Fitter1D::Fitter1D - Error: cannot create minimizer");

        m_minimiser->SetMaxFunctionCalls(1000000);
        m_minimiser->SetTolerance(0.001);
        m_minimiser->SetPrintLevel(1);
    }

    Fitter1D::Fitter1D(Fitter1D &&other) noexcept :
        m_corrFunction(std::move(other.m_corrFunction))
    {
        m_minimiser = std::move(other.m_minimiser);
        other.m_minimiser = nullptr;

        m_dataToFit = std::move(other.m_dataToFit);
        other.m_dataToFit = nullptr;
    }

    Fitter1D& Fitter1D::operator=(Fitter1D &&other) noexcept
    {
        m_corrFunction = std::move(other.m_corrFunction);

        m_minimiser = std::move(other.m_minimiser);
        other.m_minimiser = nullptr;

        m_dataToFit = std::move(other.m_dataToFit);
        other.m_dataToFit = nullptr;
    }

    bool Fitter1D::Fit()
    {
        if (! AllParamsAreSet(m_parameters))
        {
            std::cout << "Fitter1D::Fit - Warining: not all praemeters were set before fitting. The minimizer may go haywire\n";
        }

        // do the fitting (easier said than done...)
        // ROOT::Math::Functor f(&myStupidMethod,m_parameters.size());
        // minimum->SetFunction(f);

        return m_minimiser->Minimize();
    }

    void Fitter1D::SetParameter(Parameter par, float start, float min, float max)
    {
        m_minimiser->SetLimitedVariable(m_parameters.at(par).number,m_parameters.at(par).name,start,m_parameters.at(par).step,min,max);
        m_parameters[par].isSet = true;
    }

    void Fitter1D::SetParameter(Parameter par, float start)
    {
        m_minimiser->SetFixedVariable(m_parameters.at(par).number,m_parameters.at(par).name,start);
        m_parameters[par].isSet = true;
    }

    bool Fitter1D::AllParamsAreSet(const std::map<Fitter1D::Parameter,Fitter1DHelper> &pars) noexcept
    {
        for (const auto & par : pars)
            if (!par.second.isSet)
                return false;

        return true;
    }

} // namespace JJCorrFitter
