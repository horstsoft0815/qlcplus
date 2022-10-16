#include "CircularBuffer.h"

#include <numeric>

namespace CircularBufferUtils
{
std::vector<double> GetDifferenceVec(const CircularBuffer& buffer_p)
{
    std::size_t bufferSize = buffer_p.getCurrentSize();

    std::vector<double> output;
    output.reserve(bufferSize); // dont take -1, avoid special case bufferSize=0
    for(std::size_t i=0; i<bufferSize-1; ++i)
    {
        const double diff = buffer_p[i+1] - buffer_p[i];
        output.push_back(diff);
    }

    return output;
}

double operator*(const std::vector<double>& xVec_p, const std::vector<double>& yVec_p)
{
    return std::inner_product(xVec_p.begin(), xVec_p.end(), yVec_p.begin(), 0.);
}

LinearFit GetLinearFit(const std::vector<double>& xVec_p, const std::vector<double>& yVec_p)
{
    LinearFit output;
    std::vector<double> xVecm1 = xVec_p;

    for(double& entry : xVecm1)
    {
        entry -= 1.;
    }

    output.slope = (xVecm1*yVec_p)/(xVecm1*xVec_p);
    output.intercept = std::accumulate(yVec_p.begin(), yVec_p.end(), 0.) - output.slope * std::accumulate(xVec_p.begin(), xVec_p.end(), 0.);

    return output;
}

LinearFit GetLinearFit(const double timeStep_p, const CircularBuffer& buffer_p)
{
    std::size_t bufferSize = buffer_p.getCurrentSize();

    std::vector<double> xVec, yVec;
    xVec.reserve(bufferSize);
    yVec.reserve(bufferSize);
    for(std::size_t i=0; i<bufferSize; ++i)
    {
        xVec.push_back(i*timeStep_p);
        yVec.push_back(buffer_p[i]);
    }

    return GetLinearFit(xVec, yVec);
}
}
