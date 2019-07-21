#!/usr/bin/env python3
import unittest
import os

from networkit import *

class TestDistance(unittest.TestCase):
	def setUp(self):
		self.L = readGraph("input/looptest1.gml", Format.GML) #without self-loops
		self.LL = readGraph("input/looptest2.gml", Format.GML) #with self-loops sprinkled in

	def test_Diameter(self):
		D = distance.Diameter(self.LL, distance.DiameterAlgo.EstimatedRange, error = 0.1)
		D.run()
		D = distance.Diameter(self.LL, distance.DiameterAlgo.EstimatedSamples, nSamples = 5)
		D.run()
		D = distance.Diameter(self.LL, distance.DiameterAlgo.Exact)
		D.run()

	def test_Eccentricity(self):
		E = distance.Eccentricity()
		E.getValue(self.LL, 0)


	def test_EffectiveDiameter(self):
		algo = distance.EffectiveDiameter(self.L)
		algo.run()
		algo = distance.EffectiveDiameter(self.LL)
		algo.run()


	def test_ApproxEffectiveDiameter(self):
		algo = distance.EffectiveDiameterApproximation(self.L)
		algo.run()
		algo = distance.EffectiveDiameterApproximation(self.LL)
		algo.run()


	def test_ApproxHopPlot(self):
		algo = distance.HopPlotApproximation(self.L)
		algo.run()
		algo = distance.HopPlotApproximation(self.LL)
		algo.run()


	def test_NeighborhoodFunction(self):
		algo = distance.NeighborhoodFunction(self.L)
		algo.run()
		algo = distance.NeighborhoodFunction(self.LL)
		algo.run()


	def test_ApproxNeighborhoodFunction(self):
		algo = distance.NeighborhoodFunctionApproximation(self.L)
		algo.run()
		algo = distance.NeighborhoodFunctionApproximation(self.LL)
		algo.run()

	def assertEdgeSetDistance(self, G1, G2, inter, union, jaccard, overlap):
		alg = distance.EdgeSetDistance(G1, G2).run()
		self.assertEqual(alg.sizeOfIntersection(), inter)
		self.assertEqual(alg.sizeOfUnion(), union)
		self.assertAlmostEqual(alg.jaccardIndex(), jaccard)
		self.assertAlmostEqual(alg.overlapCoefficient(), overlap)
		if G1.numberOfEdges() == G2.numberOfEdges():
			self.assertAlmostEqual(alg.perturbationScore(), overlap)

	def test_EdgeSetDistance(self):
		G1 = Graph(4)
		G2 = Graph(4)

		G1.addEdge(0, 1)
		G2.addEdge(1, 2)
		self.assertEdgeSetDistance(G1, G2, 0, 2, 0.0, 0.0)

		G1.addEdge(1, 2)
		self.assertEdgeSetDistance(G1, G2, 1, 2, 0.5, 1.0)

		G2.addEdge(0, 1)
		self.assertEdgeSetDistance(G1, G2, 2, 2, 1.0, 1.0)

		G1 = Graph(3, False, True)
		G2 = Graph(3, False, True)

		G1.addEdge(0, 1)
		G2.addEdge(1, 0)
		self.assertEdgeSetDistance(G1, G2, 0, 2, 0.0, 0.0)

		G1.addEdge(1, 0)
		self.assertEdgeSetDistance(G1, G2, 1, 2, 0.5, 1.0)

if __name__ == "__main__":
	unittest.main()
