/*
 * HyperbolicGeneratorBase.h
 *
 *  Created on: 20.05.2014
 *      Author: Manuel Penschuck (networkit@manuel.jetzt), Moritz v. Looz (moritz.looz-corswarem@kit.edu)
 */

#ifndef HYPERBOLICGENERATOR_H_
#define HYPERBOLICGENERATOR_H_

#include <vector>

#include <networkit/auxiliary/StrongType.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/generators/StaticGraphGenerator.hpp>

namespace NetworKit {

namespace Hyperbolic {

using PowerlawExponent = Aux::StrongType<double, struct PowerlawExponentTag>;
using Alpha = Aux::StrongType<double, struct AlphaTag>;
using AverageDegree = Aux::StrongType<double, struct AverageDegreeTag>;
using Radius = Aux::StrongType<double, struct RadiusTag>;

/**
 * Abstract base class for Hyperbolic Random Graph Generators.
 * Implementation should override generate(bool keepInput) and pull in the constructors provided
 * (see HyperbolicGenerator as an example).
 */
class GeneratorBase : public StaticGraphGenerator {
public:
    ///@{
    /**
     * Construct a random hyperbolic graph from scratch (i.e., first sample a set of points
     * and then generate a graph).
     *
     * Either the direct model parameters (@a alpha and @a R) can be provided, or their values
     * can be estimated as a function of the expected average degree and powerlaw exponent.
     * While this estimation works quite reliably for T=0, expect discrepancies for T > 0.
     *
     * @warning Parameter estimation may change in future releases to improve quality.
     *
     * @param n     Number of points/nodes to generate
     * @param T     Temperature with T >= 0 and T != 1.0
     * @param alpha Dispersion coefficient > 0.5
     * @param R     Radius of the hyperbolic disc
     * @param deg   Average degree to be generated (model is fitted, but actual degree might vary)
     * @param ple   Exponent of the powerlaw degree distribution to generate (ple > 2.0)
     */
    GeneratorBase(count n, AverageDegree deg, PowerlawExponent ple, double T = .0);
    GeneratorBase(count n, AverageDegree deg, Alpha alpa,           double T = .0);
    GeneratorBase(count n, Radius rad,        PowerlawExponent ple, double T = .0);
    GeneratorBase(count n, Radius rad,        Alpha alpha,          double T = .0);
    ///@}

    ///@{
    /**
     * Construct a random hyperbolic graph from a prescribed point set.
     *
     *
     * @warning Parameter estimation may change in future releases to improve quality.
     *
     * @param n     Number of points/nodes to generate
     * @param T     Temperature with T >= 0 and T != 1.0
     * @param alpha Dispersion coefficient > 0.5
     * @param R     Radius of the hyperbolic disc
     * @param ple   Exponent of the powerlaw degree distribution to generate (ple > 2.0)
     */
    GeneratorBase(const std::vector<double>& angles, const std::vector<double>& radii, Radius R, PowerlawExponent ple, double T = .0);
    GeneratorBase(const std::vector<double>& angles, const std::vector<double>& radii, Radius R, Alpha alpha,          double T = .0);
    ///@}

///@name Generation
///@{
    /// @return Graph to be generated according to parameters specified in constructor freeing memory by deleting the input point set.
    Graph generate() override {
        if (angles.empty())
            samplePoints();

        return generateGraph();
    }

    /**
     * Sample point coordinates.
     * When calling generate() or generateKeepingInput() this function is called automatically if
     * needed. Only call it explicitly if you want to access point coordinates using getAngles() or
     * getRadii() before generating the graph.
     */
    void samplePoints();
///@}

///@name Parameters and point coordinates
///@{
    /// @return Model Parameter Alpha (i.e., the displacement parameter)
    double getAlpha() const noexcept {
        return alpha;
    }

    /// @return Expected powerlaw exponent of the generated graph
    double getExpectedPowerlawExponent() const noexcept {
        return estimatePowerlawExponent(alpha, temperature);
    }

    /// @return Model Parameter T (i.e., the temperature)
    double getT() const noexcept {
        return temperature;
    }

    /// @return Model Paramter R (i.e., the target radius)
    double getR() const noexcept {
        return R;
    }

    /// @return Expected powerlaw exponent of the generated graph
    double getExpectedAverageDegree() const;

    ///@{
    /**
     * @return  Angles used to generate the graph
     */
    const std::vector<double>& getAngles() const noexcept {
        return angles;
    }
    ///@}

    ///@{
    /**
     * @return  Weights used to generate the graph
     */
    const std::vector<double>& getRadii() const noexcept {
        return radii;
    }
    ///@}
///@}

    /// @return the expected powerlaw exponent based on the alpha parameter and the temperature T
    static constexpr double estimatePowerlawExponent(double alpha, double T) {
        return (T < 1) ? 2.0 * alpha + 1   : 2.0 * T * alpha + 1;
    }

    /// @return the alpha parameter based on the requested powerlaw exponent ple and the temperature T
    static constexpr double estimateAlpha(double ple, double T) {
        return (T < 1) ? 0.5 * (ple - 1) : 0.5 * (ple - 1) / T;
    }

protected:
    virtual Graph generateGraph() = 0;

    count nodeCount;    ///< Number of nodes
    double alpha;       ///< Model parameter: Dispersion factor with alpha > 0.5
    double temperature; ///< Model parameter: Temperature >= 0
    double R;           ///< Model parameter: Radius of hyperbolic disc

    double ple;         ///< User requested powerlaw exponent rather then alpha
    bool preferPle{false};

    double avgDegree;   ///< User requested avgDegree rather than R
    bool preferAvgDegree{false};

    std::vector<double> storageAngles; ///< This vector will store the angles if they are generated by the generator itself (and angles will refer to it)
    std::vector<double> storageRadii;  ///< This vector will store the radii  if they are generated by the generator itself (and radii  will refer to it)

    const std::vector<double>& angles;
    const std::vector<double>& radii;

    /// Perform checks on the input parameters (constant time in release build)
    void checkParameters();

};

} // namespace Hyperbolic


// TODO: As soon as we remove the deprecated methods this instances, Hyperbolic can directly inherit from Hyperbolic::GeneratorBase and we can implement the delegation in generateGraph().
/**
 * @ingroup generators
 * Generator for Random Hyperbolic Graphs
 *
 * This class is a dispatch which tries to identify the most suited implementation
 * for the requested set of parameters.
 */
class HyperbolicGenerator final : public StaticGraphGenerator  {
public:
    /// Legacy constructor; same as HyperbolicGenerator(count, Hyperbolic::AverageDegree, Hyperbolic::PowerlawExponent, double T)
    explicit HyperbolicGenerator(count n=10000, double avgDegree=6, double ple=3, double T=0);

   ///@{
    /**
     * Construct a random hyperbolic graph from scratch (i.e., first sample a set of points
     * and then generate a graph).
     *
     * Either the direct model parameters (@a alpha and @a R) can be provided, or their values
     * can be estimated as a function of the expected average degree and powerlaw exponent.
     * While this estimation works quite reliably for T=0, expect discrepancies for T > 0.
     *
     * @warning Parameter estimation may change in future releases to improve quality.
     *
     * @param n     Number of points/nodes to generate
     * @param T     Temperature with T >= 0 and T != 1.0
     * @param alpha Dispersion coefficient > 0.5
     * @param R     Radius of the hyperbolic disc
     * @param deg   Average degree to be generated (model is fitted, but actual degree might vary)
     * @param ple   Exponent of the powerlaw degree distribution to generate (ple > 2.0)
     */
    HyperbolicGenerator(count n, Hyperbolic::AverageDegree deg, Hyperbolic::PowerlawExponent ple, double T = .0);
    HyperbolicGenerator(count n, Hyperbolic::AverageDegree deg, Hyperbolic::Alpha alpha,          double T = .0);
    HyperbolicGenerator(count n, Hyperbolic::Radius rad,        Hyperbolic::PowerlawExponent ple, double T = .0);
    HyperbolicGenerator(count n, Hyperbolic::Radius rad,        Hyperbolic::Alpha alpha,          double T = .0);
    ///@}

    ///@{
    /**
     * Construct a random hyperbolic graph from a prescribed point set.
     *
     *
     * @warning Parameter estimation may change in future releases to improve quality.
     *
     * @param n     Number of points/nodes to generate
     * @param T     Temperature with T >= 0 and T != 1.0
     * @param alpha Dispersion coefficient > 0.5
     * @param R     Radius of the hyperbolic disc
     * @param ple   Exponent of the powerlaw degree distribution to generate (ple > 2.0)
     */
    HyperbolicGenerator(const std::vector<double>& angles, const std::vector<double>& radii, Hyperbolic::Radius R, Hyperbolic::PowerlawExponent ple, double T = .0);
    HyperbolicGenerator(const std::vector<double>& angles, const std::vector<double>& radii, Hyperbolic::Radius R, Hyperbolic::Alpha alpha,          double T = .0);
    ///@}

    virtual ~HyperbolicGenerator() = default;

///@name Generation
///@{
    /// @return Graph to be generated according to parameters specified in constructor freeing memory by deleting the input point set.
    Graph generate() override {return impl->generate();}

    /**
     * Sample point coordinates.
     * When calling generate() or generateKeepingInput() this function is called automatically if
     * needed. Only call it explicitly if you want to access point coordinates using getAngles() or
     * getRadii() before generating the graph.
     */
    void samplePoints() {return impl->samplePoints();}
///@}

///@name Parameters and point coordinates
///@{
    /// @return Model Parameter Alpha (i.e., the displacement parameter)
    double getAlpha() const noexcept {return impl->getAlpha();}

    /// @return Expected powerlaw exponent of the generated graph
    double getExpectedPowerlawExponent() const noexcept {return impl->getExpectedPowerlawExponent();}

    /// @return Model Parameter T (i.e., the temperature)
    double getT() const noexcept {return impl->getT();}

    /// @return Model Paramter R (i.e., the target radius)
    double getR() const noexcept {return impl->getR();}

    /// @return Expected powerlaw exponent of the generated graph
    double getExpectedAverageDegree() const { return impl->getExpectedAverageDegree(); }

    ///@{
    /**
     * @return  Angles used to generate the graph
     */
    const std::vector<double>& getAngles() const noexcept {return impl->getAngles();}
    ///@}

    ///@{
    /**
     * @return  Weights used to generate the graph
     */
    const std::vector<double>& getRadii() const noexcept {return impl->getRadii();}


    /// @return the expected powerlaw exponent based on the alpha parameter and the temperature T
    static constexpr double estimatePowerlawExponent(double alpha, double T) {
        return Hyperbolic::GeneratorBase::estimatePowerlawExponent(alpha, T);
    }

    /// @return the alpha parameter based on the requested powerlaw exponent ple and the temperature T
    static constexpr double estimateAlpha(double ple, double T) {
        return Hyperbolic::GeneratorBase::estimateAlpha(ple, T);
    }

    ///! forward to HyperbolicGeneratorBand::getElapsedMilliseconds() if applicable; throw otherwise
    TLX_DEPRECATED(std::vector<double> getElapsedMilliseconds() const);

    ///! forward to HyperbolicGeneratorQuadTree::setLeafCapacity() if applicable; throw otherwise
    TLX_DEPRECATED(void setLeafCapacity(count capacity));

    ///! forward to HyperbolicGeneratorQuadTree::setTheoreticalSplit(bool) if applicable; throw otherwise
    TLX_DEPRECATED(void setTheoreticalSplit(bool split));

    ///! forward to HyperbolicGeneratorQuadTree::setBalance(double) if applicable; throw otherwise
    TLX_DEPRECATED(void setBalance(double balance));

protected:

    std::unique_ptr<Hyperbolic::GeneratorBase> impl;
};

} // namespace NetworKit

#endif /* HYPERBOLICGENERATOR_H_ */
