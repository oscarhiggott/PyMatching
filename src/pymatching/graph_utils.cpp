  #include <iostream>
  #include <utility>
  #include <algorithm>
  #include <vector>
  #include "graph_utils.h"
  #include <boost/graph/graph_traits.hpp>
  #include <boost/graph/adjacency_list.hpp>
  #include <boost/graph/dijkstra_shortest_paths.hpp>

  using namespace boost;
  
  int TestGraph()
  {
    // create a typedef for the Graph type
    typedef adjacency_list<vecS, vecS, bidirectionalS> Graph;

    // Make convenient labels for the vertices
    enum { A, B, C, D, E, N };
    const int num_vertices = N;
    const char* name = "ABCDE";

    // writing out the edges in the graph
    typedef std::pair<int, int> Edge;
    Edge edge_array[] = 
    { Edge(A,B), Edge(A,D), Edge(C,A), Edge(D,C),
      Edge(C,E), Edge(B,D), Edge(D,E) };
    const int num_edges = sizeof(edge_array)/sizeof(edge_array[0]);

    // declare a graph object
    Graph g(num_vertices);

    // add the edges to the graph object
    for (int i = 0; i < num_edges; ++i)
      add_edge(edge_array[i].first, edge_array[i].second, g);
    return 0;
  }

struct EdgeData {
    int qubit_id;
    double weight;
};

  void QubitDijkstra()
  {
    typedef boost::adjacency_list < boost::listS, boost::vecS, boost::undirectedS,
    boost::no_property, EdgeData > graph_t;
    typedef graph_traits < graph_t >::vertex_descriptor vertex_descriptor;
    typedef graph_traits < graph_t >::edge_descriptor edge_descriptor;
    graph_t syndrome_graph(5);
    boost::add_edge(0, 1, {0, 3.0}, syndrome_graph);
    boost::add_edge(0, 2, {1, 1.0}, syndrome_graph);
    boost::add_edge(1, 2, {2, 7.0}, syndrome_graph);
    boost::add_edge(1, 3, {3, 5.0}, syndrome_graph);
    boost::add_edge(1, 4, {4, 1.0}, syndrome_graph);
    boost::add_edge(2, 3, {5, 2.0}, syndrome_graph);
    boost::add_edge(3, 4, {5, 7.0}, syndrome_graph);
    std::vector<double> distances(boost::num_vertices(syndrome_graph));
    std::vector<vertex_descriptor> p(num_vertices(syndrome_graph));
    vertex_descriptor from = 0;
    boost::dijkstra_shortest_paths(syndrome_graph, from,
        boost::weight_map(get(&EdgeData::weight, syndrome_graph))
        .distance_map(boost::make_iterator_property_map(distances.begin(),
                                                boost::get(boost::vertex_index, syndrome_graph)))
        .predecessor_map(&p[0]));
    
    std::cout << "distances and parents:" << std::endl;
    graph_traits < graph_t >::vertex_iterator vi, vend;
    for (tie(vi, vend) = vertices(syndrome_graph); vi != vend; ++vi) {
    std::cout << "distance(" << *vi << ") = " << distances[*vi] << ", ";
    std::cout << "parent(" << *vi << ") = " << p[*vi] << std::
        endl;
    }
    std::cout << std::endl;
    
  }

using namespace boost;

void DijkstraExample()
{
  typedef adjacency_list < listS, vecS, directedS,
    no_property, property < edge_weight_t, int > > graph_t;
  typedef graph_traits < graph_t >::vertex_descriptor vertex_descriptor;
  typedef graph_traits < graph_t >::edge_descriptor edge_descriptor;
  typedef std::pair<int, int> Edge;

  const int num_nodes = 5;
  enum nodes { A, B, C, D, E };
  char name[] = "ABCDE";
  Edge edge_array[] = { Edge(A, C), Edge(B, B), Edge(B, D), Edge(B, E),
    Edge(C, B), Edge(C, D), Edge(D, E), Edge(E, A), Edge(E, B)
  };
  int weights[] = { 1, 2, 1, 2, 7, 3, 1, 1, 1 };
  int num_arcs = sizeof(edge_array) / sizeof(Edge);
  graph_t g(edge_array, edge_array + num_arcs, weights, num_nodes);
  property_map<graph_t, edge_weight_t>::type weightmap = get(edge_weight, g);
  std::vector<vertex_descriptor> p(num_vertices(g));
  std::vector<int> d(num_vertices(g));
  vertex_descriptor s = vertex(A, g);

  dijkstra_shortest_paths(g, s, predecessor_map(&p[0]).distance_map(&d[0]));

  std::cout << "distances and parents:" << std::endl;
  graph_traits < graph_t >::vertex_iterator vi, vend;
  for (tie(vi, vend) = vertices(g); vi != vend; ++vi) {
    std::cout << "distance(" << name[*vi] << ") = " << d[*vi] << ", ";
    std::cout << "parent(" << name[*vi] << ") = " << name[p[*vi]] << std::
      endl;
  }
  std::cout << std::endl;
}