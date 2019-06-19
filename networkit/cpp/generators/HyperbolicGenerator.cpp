/*
 * HyperbolicGenerator.cpp
 *
 *      Authors: Manuel Penschuck (networkit@manuel.jetzt)
 *
 */

#include <cmath>

#include <assert.h>
#include <algorithm>

#include <networkit/auxiliary/Parallel.hpp>

#include <networkit/geometric/HyperbolicSpace.hpp>

#include <networkit/generators/HyperbolicGenerator.hpp>
#include <networkit/generators/HyperbolicGeneratorBand.hpp>
#include <networkit/generators/HyperbolicGeneratorQuadTree.hpp>


namespace NetworKit {

namespace Hyperbolic {

// Constructors
GeneratorBase::GeneratorBase(count n, AverageDegree deg, PowerlawExponent ple, double T) :
    nodeCount(n),
    alpha(estimateAlpha(ple.get(), T)),
    temperature(T),
    R(HyperbolicSpace::getTargetRadius(n, 0.5 * deg.get() * n, alpha)),
    ple(ple.get()), preferPle(true),
    avgDegree(deg.get()), preferAvgDegree(true),
    angles(storageAngles),
    radii(storageRadii)
{
    checkParameters();
}

GeneratorBase::GeneratorBase(count n, AverageDegree deg, Alpha a, double T) :
    nodeCount(n),
    alpha(a.get()),
    temperature(T),
    R(HyperbolicSpace::getTargetRadius(n, 0.5 * deg.get() * n, alpha)),
    avgDegree(deg.get()), preferAvgDegree(true),
    angles(storageAngles),
    radii(storageRadii)
{
    checkParameters();
}

GeneratorBase::GeneratorBase(count n, Radius rad, PowerlawExponent ple, double T) :
    nodeCount(n),
    alpha(estimateAlpha(ple.get(), T)),
    temperature(T),
    R(rad.get()),
    ple(ple.get()), preferPle(true),
    angles(storageAngles),
    radii(storageRadii)
{
    checkParameters();
}

GeneratorBase::GeneratorBase(count n, Radius rad, Alpha a, double T) :
    nodeCount(n),
    alpha(a.get()),
    temperature(T),
    R(rad.get()),
    angles(storageAngles),
    radii(storageRadii)
{
    checkParameters();
}

GeneratorBase::GeneratorBase(const std::vector<double>& angles, const std::vector<double>& radii, Radius R, PowerlawExponent ple, double T) :
    nodeCount(angles.size()),
    alpha(estimateAlpha(ple.get(), T)),
    temperature(T),
    R(R.get()),
    ple(ple.get()), preferPle(true),
    angles(std::move(angles)),
    radii(std::move(radii))
{
    checkParameters();
}

GeneratorBase::GeneratorBase(const std::vector<double>& angles, const std::vector<double>& radii, Radius R, Alpha a, double T) :
    nodeCount(angles.size()),
    alpha(a.get()),
    temperature(T),
    R(R.get()),
    angles(std::move(angles)),
    radii(std::move(radii))
{
    checkParameters();
}

void GeneratorBase::checkParameters() {
    if (alpha <= 0.5)
        throw std::runtime_error("Alpha has to exceed 0.5");

    if (temperature < 0 || temperature == 1)
        throw std::runtime_error("Temperature must be non-negative and not 1.");

#ifndef NDEBUG
        if (!std::all_of(angles.cbegin(), angles.cend(), [] (double x) {return 0.0 <= x && x < 2 * PI;}))
            throw std::runtime_error("Angles must be in the interval [0: 2*Pi)");

        if (!std::all_of(radii.cbegin(), radii.cend(), [this] (double r) {return 0.0 <= r && r < R;})) // TODO: replace capture with R=R in C++14
            throw std::runtime_error("Radii must be in the interval [0: R)");
#endif
}

void GeneratorBase::samplePoints() {
    if (&radii != &storageRadii)
        throw std::runtime_error("samplePoints only allowed if no points were passed via constructor");

    storageAngles.resize(nodeCount);
    storageRadii.resize(nodeCount);

    // need to generate points
    HyperbolicSpace::fillPoints(storageAngles, storageRadii, R, alpha);
}

} // namespace Hyperbolic

HyperbolicGenerator::HyperbolicGenerator(count n, double deg, double ple, double T)
    : impl(T > 0 ? dynamic_cast<Hyperbolic::GeneratorBase*>(new HyperbolicGeneratorQuadTree(n, Hyperbolic::AverageDegree{deg}, Hyperbolic::PowerlawExponent{ple}, T))
                 : dynamic_cast<Hyperbolic::GeneratorBase*>(new HyperbolicGeneratorBand    (n, Hyperbolic::AverageDegree{deg}, Hyperbolic::PowerlawExponent{ple}, T)))
{}

HyperbolicGenerator::HyperbolicGenerator(count n, Hyperbolic::AverageDegree deg, Hyperbolic::PowerlawExponent exp, double T)
    : impl(T > 0 ? dynamic_cast<Hyperbolic::GeneratorBase*>(new HyperbolicGeneratorQuadTree(n, deg, exp, T))
                 : dynamic_cast<Hyperbolic::GeneratorBase*>(new HyperbolicGeneratorBand    (n, deg, exp, T)))
{}

HyperbolicGenerator::HyperbolicGenerator(count n, Hyperbolic::AverageDegree deg, Hyperbolic::Alpha alpha,          double T)
    : impl(T > 0 ? dynamic_cast<Hyperbolic::GeneratorBase*>(new HyperbolicGeneratorQuadTree(n, deg, alpha, T))
                 : dynamic_cast<Hyperbolic::GeneratorBase*>(new HyperbolicGeneratorBand    (n, deg, alpha, T)))
{}

HyperbolicGenerator::HyperbolicGenerator(count n, Hyperbolic::Radius rad,        Hyperbolic::PowerlawExponent exp, double T)
    : impl(T > 0 ? dynamic_cast<Hyperbolic::GeneratorBase*>(new HyperbolicGeneratorQuadTree(n, rad, exp, T))
                 : dynamic_cast<Hyperbolic::GeneratorBase*>(new HyperbolicGeneratorBand    (n, rad, exp, T)))
{}

HyperbolicGenerator::HyperbolicGenerator(count n, Hyperbolic::Radius rad,        Hyperbolic::Alpha alpha,          double T)
    : impl(T > 0 ? dynamic_cast<Hyperbolic::GeneratorBase*>(new HyperbolicGeneratorQuadTree(n, rad, alpha, T))
                 : dynamic_cast<Hyperbolic::GeneratorBase*>(new HyperbolicGeneratorBand    (n, rad, alpha, T)))
{}

HyperbolicGenerator::HyperbolicGenerator(const std::vector<double>& angles, const std::vector<double>& radii, Hyperbolic::Radius R, Hyperbolic::PowerlawExponent exp, double T)
    : impl(T > 0 ? dynamic_cast<Hyperbolic::GeneratorBase*>(new HyperbolicGeneratorQuadTree(angles, radii, R, exp, T))
                 : dynamic_cast<Hyperbolic::GeneratorBase*>(new HyperbolicGeneratorBand    (angles, radii, R, exp, T)))
{}

HyperbolicGenerator::HyperbolicGenerator(const std::vector<double>& angles, const std::vector<double>& radii, Hyperbolic::Radius R, Hyperbolic::Alpha alpha,          double T)
    : impl(T > 0 ? dynamic_cast<Hyperbolic::GeneratorBase*>(new HyperbolicGeneratorQuadTree(angles, radii, R, alpha, T))
                 : dynamic_cast<Hyperbolic::GeneratorBase*>(new HyperbolicGeneratorBand    (angles, radii, R, alpha, T)))
{}

std::vector<double> HyperbolicGenerator::getElapsedMilliseconds() const {
    dynamic_cast<HyperbolicGeneratorBand&>(*impl.get()).getElapsedMilliseconds();
}

void HyperbolicGenerator::setLeafCapacity(count capacity) {
    dynamic_cast<HyperbolicGeneratorQuadTree&>(*impl.get()).setLeafCapacity(capacity);
}

void HyperbolicGenerator::setTheoreticalSplit(bool split) {
    dynamic_cast<HyperbolicGeneratorQuadTree&>(*impl.get()).setTheoreticalSplit(split);
}

void HyperbolicGenerator::setBalance(double balance) {
    dynamic_cast<HyperbolicGeneratorQuadTree&>(*impl.get()).setBalance(balance);
}

} // namespace NetworKit
