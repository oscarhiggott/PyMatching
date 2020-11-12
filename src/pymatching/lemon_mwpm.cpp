#include <iostream>
#include "lemon_mwpm.h"
#include <lemon/list_graph.h>
#include <lemon/matching.h>
using namespace lemon;
using namespace std;
int LemonTest()
{
  ListDigraph g;
  ListDigraph::Node u = g.addNode();
  ListDigraph::Node v = g.addNode();
  ListDigraph::Arc  a = g.addArc(u, v);
  cout << "Hello World! This is LEMON library here." << endl;
  cout << "We have a directed graph with " << countNodes(g) << " nodes "
       << "and " << countArcs(g) << " arc." << endl;
  return 0;
}

int MatchingTest()
{
  typedef lemon::ListGraph UGraph;
  typedef UGraph::EdgeMap<int> LengthMap;
  using lemon::INVALID;

  UGraph g;

  LengthMap length(g);

  typedef lemon::MaxWeightedPerfectMatching<UGraph,LengthMap> MWPM;
  // MWPM gg(g, length);
}


