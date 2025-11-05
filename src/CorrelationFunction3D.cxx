#include "CorrelationFunction3D.hxx"

namespace JJCorrFitter
{
    CorrelationFunction3D::CorrelationFunction3D(std::unique_ptr<SourceFunctionImpl> &&source, std::unique_ptr<InteractionTermImpl> &&interact) :
    m_minKStar(0), m_maxKStar(500), m_nPoints(100), 
    m_histogramName("hCF_fit"), m_histogramTitle("")
    {
        m_correlationFunctionName = "3D";
        m_sourceFunction = std::exchange(source,nullptr);
        m_interactionTerm = std::exchange(interact,nullptr);
        m_correlationPoints = Grid<double>(m_nPoints+1,m_nPoints+1,m_nPoints+1);
        m_correlationErrors = Grid<double>(m_nPoints+1,m_nPoints+1,m_nPoints+1);
        m_kStarValues = SetKStarPoints(m_minKStar,m_maxKStar,m_nPoints);

        m_numberOfParams = 2;
        m_totalNumberOfParams = m_numberOfParams + m_sourceFunction->GetNParams() + m_interactionTerm->GetNParams();
    }

    std::unique_ptr<TH3D> CorrelationFunction3D::MakeHistogram(const Grid<double> &points, const Grid<double> &errors, const std::vector<double> &params)
    {
        std::unique_ptr<TH3D> tmp(new TH3D(m_histogramName.data(),m_histogramTitle.data(),m_nPoints,m_minKStar,m_maxKStar,m_nPoints,m_minKStar,m_maxKStar,m_nPoints,m_minKStar,m_maxKStar));
        for (int i = 1; i <= m_nPoints; ++i)
            for (int j = 1; j <= m_nPoints; ++j)
                for (int k = 1; k <= m_nPoints; ++k)
                {
                    if (std::isnan(points(i - 1, j - 1, k - 1)))
                    {
                        tmp->SetBinContent(i,j,k,0.);
                        tmp->SetBinError(i,j,k,0.);
                    }
                    else
                    {
                        tmp->SetBinContent(i,j,k,params.at(0) * ((params.at(1) * (points(i - 1, j - 1, k - 1) - 1)) + 1));
                        tmp->SetBinError(i,j,k,errors(i - 1, j - 1, k - 1));
                    }
                }
            
        return tmp;
    }

    Grid<std::tuple<float,float,float> > CorrelationFunction3D::SetKStarPoints(float start, float stop, int nPoints)
    {
        Grid<std::tuple<float,float,float> > tmp(nPoints,nPoints,nPoints);
        const float step = (stop - start) / nPoints;

        constexpr double momBegin = 0.5;
        constexpr double momStep = 1;
        std::size_t momCounter = -1;
        std::vector<double> momBins(500,0.);
        std::generate(momBins.begin(),momBins.end(),[&momBegin,&momStep,&momCounter]{return momBegin + momStep * (++momCounter);});
        for (const auto &qOut : m_normalisationPoints)
            for (const auto &qSide : m_normalisationPoints)
                for (const auto &qLong : m_normalisationPoints)
                    momBins.push_back(CalculateModulus(qOut,qSide,qLong));

        std::sort(momBins.begin(),momBins.end());
        momBins.erase(
            std::unique(
                momBins.begin(),
                momBins.end(),
                [](const double &a, const double &b) -> bool
                {
                    return std::abs(a - b) <= std::sqrt(std::max(std::abs(a),std::abs(b)) * std::numeric_limits<double>::epsilon());
                }),
            momBins.end());

        for (int i = 0; i < nPoints; ++i)
            for (int j = 0; j < nPoints; ++j)
                for (int k = 0; k < nPoints; ++k)
                {
                    tmp(i,j,k) = std::make_tuple(start + i * step + 0.5 * step,start + j * step + 0.5 * step,start + k * step + 0.5 * step);
                }

        constexpr double rBegin = 0.5;
        constexpr double rStep = 0.5;
        std::size_t rCounter = -1;
        std::vector<double> rBins(200);
        std::generate(rBins.begin(),rBins.end(),[&rBegin,&rStep,&rCounter]{return rBegin + rStep * (++rCounter);});

        constexpr std::size_t elems = 200;
        constexpr double stepCt = 2. / elems;
        std::size_t counter = 0;
        std::vector<double> ctBins(elems,0.);
        std::generate(ctBins.begin(),ctBins.end(),[&stepCt,&counter]{return -1 + (++counter) * stepCt;});
        ctBins.push_back(1 + stepCt);

        m_interactionTerm->SetMomentumBins(std::move(momBins));
        m_interactionTerm->SetDistanceBins(std::move(rBins));
        m_interactionTerm->SetCosThetaBins(std::move(ctBins));
        m_interactionTerm->PopulateGrid();

        return tmp;
    }

    void CorrelationFunction3D::NormaliseFunction(Grid<double> &points, Grid<double> &errors)
    {
        std::size_t counter = 0;
        double sum = 0;
        for (const double &binX : m_normalisationPoints)
            for (const double &binY : m_normalisationPoints)
                for (const double &binZ : m_normalisationPoints)
                {
                    const auto [val,err] = CalculatePoint(binX,binY,binZ);
                    if (!std::isnan(val))
                    {
                        sum += val;
                        ++counter;
                    }
                }

        if (sum > 0)
        {
            for (auto &point : points)
                point *= counter / sum;
            for (auto &error : errors)
                error *= counter / sum;
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
                return gauss<double,30>::integrate(h,m_rIntRange.first,m_rIntRange.second);
            };
            return gauss<double,30>::integrate(g,m_rIntRange.first,m_rIntRange.second);
        };

        double error = 0.;
        double result = gauss_kronrod<double,61>::integrate(f,m_rIntRange.first,m_rIntRange.second,0,0,&error);
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

    void CorrelationFunction3D::SetBinning(const std::unique_ptr<TH1> &data, double minKstar, double maxKstar)
    {
        if (minKstar > maxKstar)
            throw std::logic_error("CorrelationFunction3D::SetBinning - minimal value of k* is greater than the maximal value");

        m_histogramName = std::string(data->GetName()) + "_fit";
        m_histogramTitle = std::string(data->GetTitle());

        const int minBin = (minKstar < 0) ? 0 : data->GetXaxis()->FindBin(minKstar);
        const int maxBin = (maxKstar < 0) ? data->GetNbinsX() : data->GetXaxis()->FindBin(maxKstar);

        m_minKStar = data->GetXaxis()->GetBinLowEdge(minBin);
        m_maxKStar = data->GetXaxis()->GetBinLowEdge(maxBin);
        m_nPoints = maxBin - minBin;

        m_kStarValues = SetKStarPoints(m_minKStar,m_maxKStar,m_nPoints);
        m_correlationPoints = Grid<double>(m_nPoints+1,m_nPoints+1,m_nPoints+1);
        m_correlationErrors = Grid<double>(m_nPoints+1,m_nPoints+1,m_nPoints+1);
    }

    void CorrelationFunction3D::SetBinning(const std::string &name, const std::string &title, int nPoints, double minKstar, double maxKstar)
    {
        if (minKstar > maxKstar)
            throw std::logic_error("CorrelationFunction3D::SetBinning - minimal value of k* is greater than the maximal value");

        m_histogramName = name;
        m_histogramTitle = title;
        m_nPoints = nPoints;
        m_minKStar = minKstar;
        m_maxKStar = maxKstar;

        m_kStarValues = SetKStarPoints(m_minKStar,m_maxKStar,m_nPoints);
        m_correlationPoints = Grid<double>(m_nPoints+1,m_nPoints+1,m_nPoints+1);
        m_correlationErrors = Grid<double>(m_nPoints+1,m_nPoints+1,m_nPoints+1);
    }

    std::unique_ptr<TH1> CorrelationFunction3D::Evaluate()
    {
        for (int i = 0; i < m_nPoints; ++i)
            for (int j = 0; j < m_nPoints; ++j)
                for (int k = 0; k < m_nPoints; ++k)
                {
                    const auto &[kx,ky,kz] = m_kStarValues(i,j,k);
                    const auto [val,err] = CalculatePoint(kx,ky,kz);   
                    m_correlationPoints(i,j,k) = val;
                    m_correlationErrors(i,j,k) = err;
                }

        NormaliseFunction(m_correlationPoints,m_correlationErrors); 

        return MakeHistogram(m_correlationPoints,m_correlationErrors,m_corrFuncParams);
    }

    std::unique_ptr<TH1> CorrelationFunction3D::EvaluateAtEdges(std::size_t nBins)
    {
        for (int i = 0; i < m_nPoints; ++i)
            for (int j = 0; j < nBins; ++j)
                for (int k = 0; k < nBins; ++k)
                {
                    const auto &[kx,ky,kz] = m_kStarValues(i,j,k);
                    const auto [val,err] = CalculatePoint(kx,ky,kz);
                    m_correlationPoints(i,j,k) = val;
                    m_correlationErrors(i,j,k) = err;
                }

        for (int i = 0; i < nBins; ++i)
            for (int j = 0; j < m_nPoints; ++j)
                for (int k = 0; k < nBins; ++k)
                {
                    const auto &[kx,ky,kz] = m_kStarValues(i,j,k);
                    const auto [val,err] = CalculatePoint(kx,ky,kz);   
                    m_correlationPoints(i,j,k) = val;
                    m_correlationErrors(i,j,k) = err;
                }

        for (int i = 0; i < nBins; ++i)
            for (int j = 0; j < nBins; ++j)
                for (int k = 0; k < m_nPoints; ++k)
                {
                    const auto &[kx,ky,kz] = m_kStarValues(i,j,k);
                    const auto [val,err] = CalculatePoint(kx,ky,kz);   
                    m_correlationPoints(i,j,k) = val;
                    m_correlationErrors(i,j,k) = err;
                }

        NormaliseFunction(m_correlationPoints,m_correlationErrors); 

        return MakeHistogram(m_correlationPoints,m_correlationErrors,m_corrFuncParams);
    }

    std::unique_ptr<TH1> CorrelationFunction3D::EvaluateAtPlanes(std::size_t nBins)
    {
        for (int i = 0; i < m_nPoints; ++i)
            for (int j = 0; j < m_nPoints; ++j)
                for (int k = 0; k < nBins; ++k)
                {
                    const auto &[kx,ky,kz] = m_kStarValues(i,j,k);
                    const auto [val,err] = CalculatePoint(kx,ky,kz);
                    m_correlationPoints(i,j,k) = val;
                    m_correlationErrors(i,j,k) = err;
                }

        for (int i = 0; i < nBins; ++i)
            for (int j = 0; j < m_nPoints; ++j)
                for (int k = 0; k < m_nPoints; ++k)
                {
                    const auto &[kx,ky,kz] = m_kStarValues(i,j,k);
                    const auto [val,err] = CalculatePoint(kx,ky,kz);   
                    m_correlationPoints(i,j,k) = val;
                    m_correlationErrors(i,j,k) = err;
                }

        for (int i = 0; i < m_nPoints; ++i)
            for (int j = 0; j < nBins; ++j)
                for (int k = 0; k < m_nPoints; ++k)
                {
                    const auto &[kx,ky,kz] = m_kStarValues(i,j,k);
                    const auto [val,err] = CalculatePoint(kx,ky,kz);   
                    m_correlationPoints(i,j,k) = val;
                    m_correlationErrors(i,j,k) = err;
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
