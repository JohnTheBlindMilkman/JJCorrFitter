#include "CorrelationFunction3D.hxx"

namespace JJCorrFitter
{
    CorrelationFunction3D::CorrelationFunction3D(std::unique_ptr<SourceFunctionImpl> &&source, std::unique_ptr<InteractionTermImpl> &&interact) :
    m_minKStar(0), m_maxKStar(500), m_nPoints(100), 
    m_histogramName("hCF_fit"), m_histogramTitle(""), m_kStarValues(m_nPoints)
    {
        m_correlationFunctionName = "1D";
        m_sourceFunction = std::exchange(source,nullptr);
        m_interactionTerm = std::exchange(interact,nullptr);
        m_correlationPoints = std::vector<std::vector<std::vector<double> > >(m_nPoints+1,std::vector<std::vector<double> >(m_nPoints+1,std::vector<double>(m_nPoints+1,0.)));
        m_correlationErrors = std::vector<std::vector<std::vector<double> > >(m_nPoints+1,std::vector<std::vector<double> >(m_nPoints+1,std::vector<double>(m_nPoints+1,0.)));
        m_kStarValues = SetKStarPoints(m_minKStar,m_maxKStar,m_nPoints);

        m_numberOfParams = 2;
        m_totalNumberOfParams = m_numberOfParams + m_sourceFunction->GetNParams() + m_interactionTerm->GetNParams();
    }

    std::unique_ptr<TH3D> CorrelationFunction3D::MakeHistogram(const std::vector<std::vector<std::vector<double> > > &points, const std::vector<std::vector<std::vector<double> > > &errors, const std::vector<double> &params)
    {
        std::unique_ptr<TH3D> tmp(new TH3D(m_histogramName.data(),m_histogramTitle.data(),m_nPoints,m_minKStar,m_maxKStar,m_nPoints,m_minKStar,m_maxKStar,m_nPoints,m_minKStar,m_maxKStar));
        for (int i = 1; i <= m_nPoints; ++i)
            for (int j = 1; j <= m_nPoints; ++j)
                for (int k = 1; k <= m_nPoints; ++k)
                {
                    if (std::isnan(points.at(i).at(j).at(k)) /* || std::isnan(errors.at(i).at(j).at(k)) */)
                    {
                        tmp->SetBinContent(i,j,k,0.);
                        tmp->SetBinError(i,j,k,0.);
                    }
                    else
                    {
                        tmp->SetBinContent(i,j,k,params.at(0) * ((params.at(1) * (points.at(i).at(j).at(k) - 1)) + 1));
                        tmp->SetBinError(i,j,k,errors.at(i).at(j).at(k));
                    }
                }
            
        return tmp;
    }

    std::vector<std::vector<std::vector<std::tuple<float,float,float> > > > CorrelationFunction3D::SetKStarPoints(float start, float stop, int nPoints)
    {
        std::vector<std::vector<std::vector<std::tuple<float,float,float> > > > tmp(nPoints+1,std::vector<std::vector<std::tuple<float,float,float> > >(nPoints+1,std::vector<std::tuple<float,float,float> >(nPoints+1,std::make_tuple(0.,0.,0.))));
        const float step = (stop - start) / nPoints;

        for (int i = 0; i <= nPoints; ++i)
            for (int j = 0; j <= nPoints; ++j)
                for (int k = 0; k <= nPoints; ++k)
                    tmp.at(i).at(j).at(k) = std::make_tuple(start + i*step,start + j*step,start + k*step);

        return tmp;
    }

    void CorrelationFunction3D::NormaliseFunction(std::vector<std::vector<std::vector<double> > > &points, std::vector<std::vector<std::vector<double> > > &errors)
    {
        double sum = 0;
        for (const double &binX : m_normalisationPoints)
            for (const double &binY : m_normalisationPoints)
                for (const double &binZ : m_normalisationPoints)
                {
                    //for now I don't know how to get around this
/*                     if (binX < m_maxKStar || binY < m_maxKStar || binZ < m_maxKStar)
                    {
                        sum += points.at(std::lower_bound(m_kStarValues.begin(),m_kStarValues.end(),binX) - m_kStarValues.begin())
                            .at(std::lower_bound(m_kStarValues.begin(),m_kStarValues.end(),binY) - m_kStarValues.begin())
                            .at(std::lower_bound(m_kStarValues.begin(),m_kStarValues.end(),binZ) - m_kStarValues.begin());
                    }
                    else
                    { */
                        double val, err;
                        std::tie(val,err) = CalculatePoint(binX,binY,binZ);
                        if (!std::isnan(val))
                            sum += val;
                    //}
                }
        for (int i = 0; i <= m_nPoints; ++i)
            for (int j = 0; j <= m_nPoints; ++j)
                for (int k = 0; k <= m_nPoints; ++k)
                {
                    points.at(i).at(j).at(k) *= (std::pow(m_normalisationPoints.size(),3)/sum);
                    errors.at(i).at(j).at(k) *= (std::pow(m_normalisationPoints.size(),3)/sum);
                }
    }

    std::pair<double,double> CorrelationFunction3D::CalculatePoint()
    {
        std::cerr << "CorrelationFunction3D::CalculatePoint - This is a deleted function" << std::endl; 
        return std::make_pair(0,0);
    }

    std::pair<double,double> CorrelationFunction3D::CalculatePoint(float kOut,float kSide,float kLong)
    {
        using boost::math::quadrature::gauss;
        using boost::math::quadrature::gauss_kronrod;
        
        m_interactionTerm->SetMomentum(CalculateModulus(kOut,kSide,kLong));

        auto kooninPratt = [&](double ro, double rs, double rl){return m_sourceFunction->GetValue(ro,rs,rl) * (m_interactionTerm->GetValue(CalculateModulus(ro,rs,rl),CalculateCosTheta(kOut,kSide,kLong,ro,rs,rl)));};

        auto f = [&](const double rOut)
        {
            auto g = [&](const double rSide)
            {
                auto h = [&](const double rLong)
                {
                    return kooninPratt(rOut,rSide,rLong);
                };
                return gauss<double,30>::integrate(h,0.,std::numeric_limits<double>::infinity());
            };
            return gauss<double,30>::integrate(g,0.,std::numeric_limits<double>::infinity());
        };

        double error = 0.;
        double result = gauss_kronrod<double,61>::integrate(f,0.,std::numeric_limits<double>::infinity(),0,0,&error);
        return std::make_pair(result,error);
    }

    constexpr double CorrelationFunction3D::CalculateCosTheta(float kOut,float kSide,float kLong, float rOut,float rSide,float rLong) noexcept
    {
        return (kOut*rOut + kSide*rSide + kLong*rLong) / (CalculateModulus(kOut,kSide,kLong) * CalculateModulus(rOut,rSide,rLong));
    }

    constexpr double CorrelationFunction3D::CalculateModulus(float xOut,float xSide,float xLong) noexcept
    {
        return std::sqrt(xOut*xOut + xSide*xSide + xLong*xLong);
    }

    void CorrelationFunction3D::SetBinning(const std::unique_ptr<TH1> &data, float minKstar, float maxKstar)
    {
        if (minKstar > maxKstar)
            throw std::logic_error("CorrelationFunction3D::SetBinning - minimal value of k* is greater than the maximal value");

        m_histogramName = std::string(data->GetName()) + "_fit";
        m_histogramTitle = std::string(data->GetTitle());

        const int minBin = (minKstar < 0) ? 0 : data->FindBin(minKstar);
        const int maxBin = (maxKstar < 0) ? data->GetNbinsX() : data->FindBin(maxKstar);

        m_minKStar = data->GetBinLowEdge(minBin);
        m_maxKStar = data->GetBinLowEdge(maxBin);
        m_nPoints = maxBin - minBin;

        m_kStarValues = SetKStarPoints(m_minKStar,m_maxKStar,m_nPoints);
        m_correlationPoints = std::vector<std::vector<std::vector<double> > >(m_nPoints+1,std::vector<std::vector<double> >(m_nPoints+1,std::vector<double>(m_nPoints+1,0.)));
        m_correlationErrors = std::vector<std::vector<std::vector<double> > >(m_nPoints+1,std::vector<std::vector<double> >(m_nPoints+1,std::vector<double>(m_nPoints+1,0.)));
    }

    void CorrelationFunction3D::SetBinning(const std::string &name, const std::string &title, int nPoints, float minKstar, float maxKstar)
    {
        if (minKstar > maxKstar)
            throw std::logic_error("CorrelationFunction3D::SetBinning - minimal value of k* is greater than the maximal value");

        m_histogramName = name;
        m_histogramTitle = title;
        m_nPoints = nPoints;
        m_minKStar = minKstar;
        m_maxKStar = maxKstar;

        m_kStarValues = SetKStarPoints(m_minKStar,m_maxKStar,m_nPoints);
        m_correlationPoints = std::vector<std::vector<std::vector<double> > >(m_nPoints+1,std::vector<std::vector<double> >(m_nPoints+1,std::vector<double>(m_nPoints+1,0.)));
        m_correlationErrors = std::vector<std::vector<std::vector<double> > >(m_nPoints+1,std::vector<std::vector<double> >(m_nPoints+1,std::vector<double>(m_nPoints+1,0.)));
    }

    std::unique_ptr<TH1> CorrelationFunction3D::Evaluate()
    {
        std::size_t counter = 0;
        double val, err, kx, ky, kz;
        /* for (int i = 0; i <= m_nPoints; ++i)
            for (int j = 0; j <= m_nPoints; ++j)
                for (int k = 0; k <= m_nPoints; ++k)
                {
                    std::tie(kx,ky,kz) = m_kStarValues.at(i).at(j).at(k);
                    std::tie(val,err) = CalculatePoint(kx,ky,kz);   
                    std::cout << counter << ": [" << kx << ',' << ky << ',' << kz << "] -> " << val << " +/- " << err << "\n";
                    m_correlationPoints.at(i).at(j).at(k) = val;
                    m_correlationErrors.at(i).at(j).at(k) = err;
                    ++counter;
                } */

        for (int i = 0; i <= m_nPoints; ++i)
        {
            std::tie(kx,ky,kz) = m_kStarValues.at(i).at(1).at(1);
            std::tie(val,err) = CalculatePoint(kx,ky,kz);
            std::cout << counter << ": [" << kx << ',' << ky << ',' << kz << "] -> " << val << " +/- " << err << "\n";
            m_correlationPoints.at(i).at(1).at(1) = val;
            m_correlationErrors.at(i).at(1).at(1) = err;
            ++counter;
        }

        for (int j = 0; j <= m_nPoints; ++j)
        {
            std::tie(kx,ky,kz) = m_kStarValues.at(1).at(j).at(1);
            std::tie(val,err) = CalculatePoint(kx,ky,kz);   
            std::cout << counter << ": [" << kx << ',' << ky << ',' << kz << "] -> " << val << " +/- " << err << "\n";
            m_correlationPoints.at(1).at(j).at(1) = val;
            m_correlationErrors.at(1).at(j).at(1) = err;
            ++counter;
        }

        for (int k = 0; k <= m_nPoints; ++k)
        {
            std::tie(kx,ky,kz) = m_kStarValues.at(1).at(1).at(k);
            std::tie(val,err) = CalculatePoint(kx,ky,kz);   
            std::cout << counter << ": [" << kx << ',' << ky << ',' << kz << "] -> " << val << " +/- " << err << "\n";
            m_correlationPoints.at(1).at(1).at(k) = val;
            m_correlationErrors.at(1).at(1).at(k) = err;
            ++counter;
        }

        NormaliseFunction(m_correlationPoints,m_correlationErrors); 

        return MakeHistogram(m_correlationPoints,m_correlationErrors,m_corrFuncParams);
    }

    void CorrelationFunction3D::SetParameters(const std::vector<double> &generalPars,const std::vector<double> &srcPars,const std::vector<double> &psiPars)
    {
        if (generalPars.size() != m_numberOfParams)
        {
            throw std::length_error("CorrelationFunction3D::SetParameters - Provided number of parameters does not match the expected number: " + std::to_string(m_numberOfParams));
        }

        m_corrFuncParams = generalPars;
        m_sourceFunction->SetParameters(srcPars);
        m_interactionTerm->SetParameters(psiPars);
    }
} // namespace JJCorrFitter
