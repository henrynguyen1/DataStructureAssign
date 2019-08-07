/*
 * Notice that the list of included headers has
 * expanded a little. As before, you are not allowed
 * to add to this.
 */
#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <array>
#include <list>
#include <forward_list>
#include <deque>
#include <map>
#include <cstddef>
#include <string>
#include <utility>
#include <algorithm>
#include <limits>
#include <optional>
#include <exception>
#include <stdexcept>

#include "directed_graph.hpp"

using namespace std;




/*
 * Computes whether the input is a Directed Acyclic Graph (DAG).
 * A digraph is a DAG if there is no vertex that has a cycle.
 * A cycle is a non-empty set of [out-]edges that starts at one 
 * vertex, and returns to it.
 */
template <typename vertex>
bool is_dag(const directed_graph<vertex> & d) {
    // Initialise directed graph to test.
    directed_graph<vertex> graph(d);
   // iterate through the graphs vertex 
    for(auto i = graph.begin(); i != graph.end(); i++){
        //if graph has no vertices it has no cycle
        if(graph.num_vertices() <= 1){
          return true;
           } 
        // if any vertices have no out degree, the graph has no cycle
         if(graph.out_degree(*i) == 0){
             graph.remove_vertex(*i);
             return true;
           }
         }
   return false;
}
  
/*
 * Computes a topological ordering of the vertices.
 * For every vertex u in the order, and any of its
 * neighbours v, v appears later in the order than u.
 */
template <typename vertex>
std::list<vertex> topological_sort(const directed_graph<vertex> & d) {
    directed_graph<vertex> graph(d);
    list<vertex> topo_order;
  //vertices without incoming edges
  set<vertex> s;
  for (vertex v : graph.get_Vertices()){
     if (graph.in_degree(v)==0)
        s.insert(v);
  }
  // while ther are vertices in the set
  while (!s.empty()) {
     vertex v = *s.begin();
     s.erase(v);
     topo_order.push_back(v);
     // for each neighbour of u and v
     for (vertex u: graph.get_Neighbours(v)) {
     // remove (v,u) from graph;
     graph.remove_edge(v,u);
   //    if (u has no incoming edges) add u to s;
     if (graph.in_degree(u) ==0)
         s.insert(u);
        }
  }
      
    return topo_order;
}

/*
 * Given a DAG, computes whether there is a Hamiltonian path.
 * a Hamiltonian path is a path that visits every vertex
 * exactly once.
 */
template <typename vertex>
bool is_hamiltonian_dag(const directed_graph<vertex> & d) {
  queue <vertex> q;
  map<vertex, bool> m;
    
    
   while(!q.empty()){
       vertex u = q.front();
       q.pop();
       for (auto i = d.nbegin(u); i != d.nend(u); i++){
       if(d.adjacent(*i, u)) {
           return false;// not DAG
         }
       if(!m[*i]){
           m[*i] = true;
          q.push(*i);
         }
      } 
   }
    
return true;
}

/*
 * Computes the weakly connected components of the graph.
 * A [weak] component is the smallest subset of the vertices
 * such that the in and out neighbourhood of each vertex in
 * the set is also contained in the set.
 */
template <typename vertex>
std::vector<std::vector<vertex>> components(const directed_graph<vertex> & d) {
 //Prepare vector of vector of vertices for weak components and map of visited maps vertex
    vector <vector<vertex>> weakconnect;
    map<vertex, bool> visited;
    
    
//implement Breadth First Search
for(auto i = d.begin(); i != d.end(); i++){
   // intitalise a queue and add the inputted vertex mark all as unvisited
   if(!visited[*i]){
      queue<vertex> q;
      vector<vertex> vect;
       //add inputted vertex to queue and set to visited
       q.push(*i);
       visited[*i] = true;
       // while there are vertices in the queue
       while (!q.empty()){
       // initialise the top vertex and remove it from the queue
       vertex v = q.front();
       q.pop();
       // add vertex to the vectored for ordered vertices
       vect.push_back(v);
        //if vertex contains a neighbour then add it to the end of the queue .
       for (auto j = d.begin(); j != d.end(); j++){
           if(d.adjacent(*j, v) || d.adjacent(v, *j)){
               if(!visited[*j]){
               visited[*j] = true;
               q.push(*j);
              } 
            }
          }
        }
  weakconnect.push_back(vect);
       }
     }
return weakconnect;
}

/*
 * Computes the strongly connected components of the graph.
 * A strongly connected component is a subset of the vertices
 * such that for every pair u, v of vertices in the subset,
 * v is reachable from u and u is reachable from v.
 */

template <typename vertex>
std::vector<std::vector<vertex>> strongly_connected_components(const directed_graph<vertex> & d) {
    int index = 0;	
    vector <vector<vertex>> strongconnect;
    map<vertex, bool> visited;


// Implement Depth First Search
   for(auto i = d.begin(); i != d.end(); i++){
   if(!visited[*i]){
   stack<vertex> s;
   vector<vertex> vect;
       s.push(*i);
       visited[*i] = true;
       // while there are values in the stack
       while (!s.empty()){
      // initialise the top vertex and remove it from the stack
       vertex v = s.top();
       s.pop();
       vect.push_back(v);
        //if vertex contains a neighbour then set to visited.
       for (auto j = d.begin(); j != d.end(); j++){
           if(d.adjacent(*j, v) || d.adjacent(v, *j)){
               if(!visited[*j]){
               visited[*j] = true;
               s.push(*j);
               }
             }
           }
         }
    strongconnect.push_back(vect);
      }
    }
  return strongconnect;
}



/*
 * Computes the shortest distance from u to every other vertex
 * in the graph d. The shortest distance is the smallest number
 * of edges in any path from u to the other vertex.
 * If there is no path from u to a vertex, set the distance to
 * be the number of vertices in d plus 1.
 */
template <typename vertex>
std::unordered_map<vertex, std::size_t> shortest_distances(const directed_graph<vertex> & d, const vertex & u) {
  //Map of the distance - map each vertex to an integer distance from the source vertex u
  unordered_map<vertex, size_t> ordered;
  // keeps track of distance travelled from starting vertex 
  size_t distance=0;
  //Counter for distance of compared vertices from start of traversal
  size_t counter=1;
  queue <vertex> q;
  q.push(u);
  
  
  
	for (auto i = d.begin(); i != d.end(); i++){
		ordered[*i] = d.num_vertices()+1;
   }
  ordered[u] = 0;
  
  //Implement Dijkstra's algorithm 

  while(!q.empty()){
	  // initialise the top vertex and remove it from the queue
      vertex u = q.front();
    //Select the unvisited vertex with smallest distance and set it as the current vertex.
      if(counter <= 0){
        counter = q.size();
        distance++;
        }
        counter--;
        q.pop();
        ordered[u] = distance;
		 
	 //if vertex contains a neighbour then add it to the end of the queue .	  
     for (auto i = d.nbegin(u); i != d.nend(u); i++){
				if (ordered[*i] == d.num_vertices()+1){
  // distance of neighbour vertex is plus one of current position in traversal
          ordered[*i] = ordered[u] +1;
					q.push(*i);
          }
     }
 }
    return ordered;
}
