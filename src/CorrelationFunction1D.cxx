#include "CorrelationFunction1D.hxx"

namespace JJCorrFitter
{
    CorrelationFunction1D::CorrelationFunction1D()
    {
    }

    double CorrelationFunction1D::IntegralExpression(double x)
    {

    }

    double CorrelationFunction1D::Evaluate(int kStar, double rStar, double rInv)
    {
        ROOT::Math::WrappedFunction<> wf(IntegralExpression);
        ROOT::Math::Integrator integ(wf);
    }
} // namespace JJCorrFitter
