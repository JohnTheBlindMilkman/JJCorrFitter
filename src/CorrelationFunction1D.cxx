#include "CorrelationFunction1D.hxx"

namespace JJCorrFitter
{
    CorrelationFunction1D::CorrelationFunction1D(std::unique_ptr<SourceFunctionImpl> &&source, std::unique_ptr<InteractionTermImpl> &&interact) :
    m_minKStar(0), m_maxKStar(500), m_nPoints(100), 
    m_histogramName("hCF_fit"), m_histogramTitle(""), m_kStarValues(m_nPoints)
    {
        m_correlationFunctionName = "1D";
        m_sourceFunction = std::exchange(source,nullptr);
        m_interactionTerm = std::exchange(interact,nullptr);
        m_correlationPoints = {};
        m_correlationErrors = {};
        m_kStarValues = SetKStarPoints(m_minKStar,m_maxKStar,m_nPoints);

        m_numberOfParams = 2;
        m_totalNumberOfParams = m_numberOfParams + m_sourceFunction->GetNParams() + m_interactionTerm->GetNParams();
    }

    std::unique_ptr<TH1D> CorrelationFunction1D::MakeHistogram(const std::vector<double> &points, const std::vector<double> &errors, const std::vector<double> &params)
    {
        std::unique_ptr<TH1D> tmp(new TH1D(m_histogramName.data(),m_histogramTitle.data(),m_nPoints,m_minKStar,m_maxKStar));
        for (int bin = 1; bin <= m_nPoints; ++bin)
        {
            tmp->SetBinContent(bin,params.at(0) * ((params.at(1) * (points.at(bin) - 1)) + 1));
            tmp->SetBinError(bin,errors.at(bin));
        }
            
        return tmp;
    }

    std::vector<double> CorrelationFunction1D::SetKStarPoints(double start, double stop, int nPoints)
    {
        std::vector<double> tmp(nPoints+1);
        const double step = (stop - start) / nPoints;

        for (int i = 0; i <= nPoints; ++i)
            tmp[i] = start + i*step;

        std::vector<double> tmp2 = tmp;
        tmp2.insert(tmp2.end(),m_normalisationPoints.begin(),m_normalisationPoints.end());
        std::sort(tmp2.begin(),tmp2.end());
        m_interactionTerm->SetMomentumBins(std::move(tmp2));
        m_interactionTerm->SetDistanceBins(CalculateIntSamplePoints(m_rIntRange,boost::math::quadrature::gauss_kronrod<double,61>::abscissa()));
        m_interactionTerm->SetCosThetaBins(CalculateIntSamplePoints(m_cosThetaIntRange,boost::math::quadrature::gauss<double,30>::abscissa()));
        m_interactionTerm->PopulateGrid();

        return tmp;
    }

    void CorrelationFunction1D::NormaliseFunction(std::vector<double> &points, std::vector<double> &errors)
    {
        double sum = 0;
        for (const double &bin : m_normalisationPoints)
        {
            if (bin < m_maxKStar)
            {
                sum += points.at(std::lower_bound(m_kStarValues.begin(),m_kStarValues.end(),bin) - m_kStarValues.begin());
            }
            else
            {
                double val, err;
                std::tie(val,err) = CalculatePoint(bin);
                sum += val;
            }
        }
        
        std::transform(points.begin(),points.end(),points.begin(), [&](auto &elem){return elem * (m_normalisationPoints.size()/sum);});
        std::transform(errors.begin(),errors.end(),errors.begin(), [&](auto &elem){return elem * (m_normalisationPoints.size()/sum);});
    }

    std::pair<double,double> CorrelationFunction1D::CalculatePoint()
    {
        std::cerr << "CorrelationFunction1D::CalculatePoint - This is a deleted function" << std::endl; 
        return std::make_pair(0,0);
    }

    std::pair<double,double> CorrelationFunction1D::CalculatePoint(double kStar)
    {
        using boost::math::quadrature::gauss;
        using boost::math::quadrature::gauss_kronrod;
        
        m_interactionTerm->SetMomentum(kStar);

        // TPI calculates wrong becasue of the 4 * pi * r * r factor!
        auto kooninPratt = [&](double r, double ctheta) -> double {return 4 * M_PI * r * r * m_sourceFunction->GetValue(r) * (m_interactionTerm->GetValue(r,ctheta));};

        auto f = [&](const double r) -> double
        {
            auto g = [&](const double ctheta) -> double
            {
                return kooninPratt(r,ctheta);
            };
            return gauss<double,30>::integrate(g,-1,1);
        };

        double error = 0.;
        double result = gauss_kronrod<double,61>::integrate(f,0.,50.,0,0,&error); // should be 0 to inf, but nothing happens above ~20 fm
        return std::make_pair(result,error);
    }

    void CorrelationFunction1D::SetBinning(const std::unique_ptr<TH1> &data, double minKstar, double maxKstar)
    {
        if (minKstar > maxKstar)
            throw std::logic_error("CorrelationFunction1D::SetBinning - minimal value of k* is greater than the maximal value");

        m_histogramName = std::string(data->GetName()) + "_fit";
        m_histogramTitle = std::string(data->GetTitle());

        const int minBin = (minKstar < 0) ? 0 : data->FindBin(minKstar);
        const int maxBin = (maxKstar < 0) ? data->GetNbinsX() : data->FindBin(maxKstar);

        m_minKStar = data->GetBinLowEdge(minBin);
        m_maxKStar = data->GetBinLowEdge(maxBin);
        m_nPoints = maxBin - minBin;

        m_kStarValues = SetKStarPoints(m_minKStar,m_maxKStar,m_nPoints);
    }

    void CorrelationFunction1D::SetBinning(const std::string &name, const std::string &title, int nPoints, double minKstar, double maxKstar)
    {
        if (minKstar > maxKstar)
            throw std::logic_error("CorrelationFunction1D::SetBinning - minimal value of k* is greater than the maximal value");

        m_histogramName = name;
        m_histogramTitle = title;
        m_nPoints = nPoints;
        m_minKStar = minKstar;
        m_maxKStar = maxKstar;

        m_correlationPoints.reserve(m_nPoints);
        m_correlationErrors.reserve(m_nPoints);
        m_kStarValues = SetKStarPoints(m_minKStar,m_maxKStar,m_nPoints);
    }

    std::unique_ptr<TH1> CorrelationFunction1D::Evaluate()
    {
        m_correlationPoints.clear();
        m_correlationPoints.resize(0);
        m_correlationErrors.clear();
        m_correlationErrors.resize(0);

        for (int bin = 0; bin <= m_nPoints; ++bin)
        {
            double val, err;
            std::tie(val,err) = CalculatePoint(m_kStarValues.at(bin));   
            m_correlationPoints.push_back(val);
            m_correlationErrors.push_back(err);
        }

        NormaliseFunction(m_correlationPoints,m_correlationErrors); 

        return MakeHistogram(m_correlationPoints,m_correlationErrors,m_corrFuncParams);
    }

    void CorrelationFunction1D::SetParameters(const std::vector<double> &generalPars,const std::vector<double> &srcPars,const std::vector<double> &psiPars)
    {
        if (generalPars.size() != m_numberOfParams)
        {
            throw std::length_error("CorrelationFunction1D::SetParameters - Provided number of parameters does not match the expected number: " + std::to_string(m_numberOfParams));
        }

        m_corrFuncParams = generalPars;
        m_sourceFunction->SetParameters(srcPars);
        m_interactionTerm->SetParameters(psiPars);
    }
} // namespace JJCorrFitter
