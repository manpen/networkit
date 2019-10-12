#!/usr/bin/env python3
import unittest
import networkit as nk

class TestGraphTools(unittest.TestCase):
    def testMakeBipartite(self):
        # create line graph
        G = nk.Graph(100)
        for i in range(G.numberOfNodes() - 1):
            G.addEdge(i, i+1)

        bipariteG = nk.graph.GraphTools.makeBipartite(G, 50)

        self.assertEqual(bipariteG.numberOfEdges(), 1)

    def testMakeKPartite(self):
        # create line graph
        G = nk.Graph(100)
        part = nk.Partition(100)
        for i in range(G.numberOfNodes() - 1):
            G.addEdge(i, i+1)
            part[i] = 0 if i < 50 else 1
        part[99] = 1

        bipariteG = nk.graph.GraphTools.makeKPartite(G, part)
        self.assertEqual(bipariteG.numberOfEdges(), 1)

    def testMakeKPartiteIndependentSet(self):
        # create line graph and parition that is an independent set
        G = nk.Graph(100)
        part = nk.Partition(100)
        for i in range(G.numberOfNodes() - 1):
            G.addEdge(i, i+1)
            part[i] = i % 2
        part[99] = 1

        bipariteG = nk.graph.GraphTools.makeKPartite(G, part)
        self.assertEqual(bipariteG.numberOfEdges(), G.numberOfEdges())


if __name__ == "__main__":
    unittest.main()
