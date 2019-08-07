#ifndef DIRECTED_GRAPH_H
#define DIRECTED_GRAPH_H

//A large selection of data structures from the standard
//library. You need not feel compelled to use them all,
//but as you can't add any, they're all here just in case.
#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <iostream>
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
using namespace std;
//Forward declarations for classes below so they can be used below without worrying too much about the ordering.
template <typename vertex> class vertex_iterator;
template <typename vertex> class neighbour_iterator;
template <typename vertex> class directed_graph;


template <typename vertex>
class directed_graph {

private:
  vector<vector<bool>> adj_matrix;// stores adjacency matrix of the graph
  size_t size; // stores size of vector
  int edge_count; // stores total number of edges within the graph
  vector<vertex> vertices; // stores the vertices of the graph


public:


  directed_graph(); //A constructor for directed_graph. The graph should start empty.
  ~directed_graph(); //A destructor. Depending on how you do things, this may
  //not be necessary.
  
  bool contains(const vertex&) const; //Returns true if the given vertex is in the graph, false otherwise.

  bool adjacent(const vertex&, const vertex&) const; //Returns true if the first vertex is adjacent to the second, false otherwise.

  void add_vertex(const vertex&); //Adds the passed in vertex to the graph (with no edges).
  void add_edge(const vertex&, const vertex&); //Adds an edge from the first vertex to the second.

  void remove_vertex(const vertex&); //Removes the given vertex. Should also clear any incident edges.
  void remove_edge(const vertex&, const vertex&); //Removes the edge between the two vertices, if it exists.

  std::size_t in_degree(const vertex&) const; //Returns number of edges coming in to a vertex.
  std::size_t out_degree(const vertex&) const; //Returns the number of edges leaving a vertex.
  std::size_t degree(const vertex&) const; //Returns the degree of the vertex (both in and out edges).
  
  std::size_t num_vertices() const; //Returns the total number of vertices in the graph.
  std::size_t num_edges() const; //Returns the total number of edges in the graph.

  std::vector<vertex> get_vertices(); //Returns a vector containing all the vertices.
  std::vector<vertex> get_neighbours(const vertex&); //Returns a vector containing the neighbours of the given vertex.

  vertex_iterator<vertex> begin(); //Returns a graph_iterator pointing to the start of the vertex set.
  vertex_iterator<vertex> end(); //Returns a graph_iterator pointing to one-past-the-end of the vertex set.

  neighbour_iterator<vertex> nbegin(const vertex&); //Returns a neighbour_iterator pointing to the start of the neighbour set for the given vertex.
  neighbour_iterator<vertex> nend(const vertex&); //Returns a neighbour_iterator pointing to one-past-the-end of the neighbour set for the given vertex.

  std::vector<vertex> depth_first(const vertex&); //Returns the vertices of the graph in the order they are visited in by a depth-first traversal starting at the given vertex.
  std::vector<vertex> breadth_first(const vertex&); //Returns the vertices of the graph in the order they are visisted in by a breadth-first traversal starting at the given vertex.

  directed_graph<vertex> out_tree(const vertex&); //Returns a spanning tree of the graph starting at the given vertex using the out-edges.
  directed_graph<vertex> in_tree(const vertex&); //Returns a spanning tree of the graph starting at the given vertex using the in-edges.

  bool reachable(const vertex&, const vertex&) const; //Returns true if the second vertex is reachable from the first (can you follow a path of out-edges to get from the first to the second?). Returns false otherwise.

  int get_index(const vertex& u) const;
};

//The vertex_iterator class provides an iterator
//over the vertices of the graph.
//This is one of the harder parts, so if you're
//not too comfortable with C++ leave this for last.
//If you are, there are many ways of doing this,
//as long as it passes the tests, it's okay.
//You may want to watch the videos on iterators before starting.
template <typename vertex>
class vertex_iterator {

private:

  //You may need data members here.
  vector<vertex> v;
  size_t pos;

public:
  vertex_iterator(const vertex_iterator<vertex>&);
  vertex_iterator(const directed_graph<vertex>&, std::size_t);
  ~vertex_iterator();
  vertex_iterator<vertex> operator=(const vertex_iterator<vertex>&);
  bool operator==(const vertex_iterator<vertex>&) const;
  bool operator!=(const vertex_iterator<vertex>&) const;
  vertex_iterator<vertex> operator++();
  vertex_iterator<vertex> operator++(int);
  vertex operator*();
  vertex* operator->();
};

//The neighbour_iterator class provides an iterator
//over the neighbours of a given vertex. This is
//probably harder (conceptually) than the graph_iterator.
//Unless you know how iterators work.
template <typename vertex>
class neighbour_iterator {

private:

  //You may need data members here.
    vector<vertex> n;
    int row_index;
    size_t pos;

public:
  neighbour_iterator(const neighbour_iterator<vertex>&);
  neighbour_iterator(const directed_graph<vertex>&, const vertex&, std::size_t);
  ~neighbour_iterator();
  neighbour_iterator<vertex> operator=(const neighbour_iterator<vertex>&);
  bool operator==(const neighbour_iterator<vertex>&) const;
  bool operator!=(const neighbour_iterator<vertex>&) const;
  neighbour_iterator<vertex> operator++();
  neighbour_iterator<vertex> operator++(int);			
  vertex operator*();
  vertex* operator->();
};


//Define all your methods down here (or move them up into the header, but be careful you don't double up). If you want to move this into another file, you can, but you should #include the file here.
//Although these are just the same names copied from above, you may find a few more clues in the full
//method headers. Note also that C++ is sensitive to the order you declare and define things in - you
//have to have it available before you use it.



// Directed Graph


template <typename vertex> directed_graph<vertex>::directed_graph() {
	// constructor, set initial values
	 edge_count = 0;
}
template <typename vertex> directed_graph<vertex>::~directed_graph() {}
template <typename vertex> bool directed_graph<vertex>::contains(const vertex& u) const { 
   // for each vertex in vertices vector
	for (vertex v : vertices) {
		// if it is equal to the vertex we are looking for
        if (v == u) 
			// return true if found
         return true;
    }
	// else return false, if it has not been found
        return false;
}
	 
template <typename vertex> bool directed_graph<vertex>::adjacent(const vertex& u, const vertex& v) const {
	// if vertices are valid and the vertice inputed is not equal to what we are looking for and graph contains such vertices
	return adj_matrix[get_index(u)][get_index(v)] && u != v && contains(u) && contains(v);
}

template <typename vertex> int directed_graph<vertex>::get_index(const vertex& u) const {
	// for each vertex
	for (unsigned int i = 0; i < vertices.size(); i++) {
		// if the vertex equal to the searched vertex
		if (vertices[i] == u)
			//return index of the vertex position
			return i;
	}
	// else return if vertex does not exist
	return -1;
}
template <typename vertex> void directed_graph<vertex>::add_vertex(const vertex& u) {
	// if vertex inputed does not exist
	if(!contains(u)) {
		// add new vertex to vertices vector
		vertices.push_back(u);
		// add false value at the end of each row
		for (auto& row : adj_matrix){
			row.push_back(false);
		}
		// create another vector row with false as default values and add it to the graph
		adj_matrix.push_back(vector<bool>(vertices.size(), false));
	}
}
template <typename vertex> void directed_graph<vertex>::add_edge(const vertex& u, const vertex& v) {
	// get index position of vertices inputed and comparison too
	int u_pos = get_index(u),v_pos = get_index(v);
	// if graph doesnt contain vertices s or indexes are valid
	if (!contains(u) || !contains(v) || u == v)
	//terminate loop	
	return; 
	// set vertices index position to true in adjacency matrix and increment edge count by 1
	adj_matrix[u_pos][v_pos] = true;
	edge_count++;
	
	 
}
template <typename vertex> void directed_graph<vertex>::remove_vertex(const vertex& u) {
	// get vertex index
	int index = get_index(u);
	for (unsigned int i = 0; i < adj_matrix[index].size(); i++) {
		// decrease edge count if index is more than zero
			if (adj_matrix[index][i] > 0) 
				edge_count--;
		}
	    //remove vertex from vertices vector
		vertices.erase(vertices.begin() + index);
	   // remove vertex row from adj_matrix
		adj_matrix.erase(adj_matrix.begin() + index);
	    // remove vertex values from each row in adj_matrix
		for (unsigned int i = 0; i < adj_matrix.size(); i++) {
			adj_matrix[i].erase(adj_matrix[i].begin() + index);
		}
    
}
template <typename vertex> void directed_graph<vertex>::remove_edge(const vertex& u, const vertex& v) {
	// if graph doesnt vertices for such vertices and indexes are valid
	if (!contains(u) || !contains(v) || u == v)
    // terminate loop
	   return;
	//get index position of bother vertices
	int u_pos = get_index(u),v_pos = get_index(v);
	// set vertices index position to false in adjacency matrix and decrement edge count by 1
	adj_matrix[u_pos][v_pos] = false;
    edge_count--; 
}
template <typename vertex> std::size_t directed_graph<vertex>::in_degree(const vertex& u) const { 
	// initialise in degree value and set to zero
	int in_degree = 0;
	// get index for vertex
	int index = get_index(u);
	// if at u index position for v is equal to 1 there is an in-degree
	for (unsigned i =  0; i < adj_matrix.size(); i++){
		if(adj_matrix[i][index] == 1){
	    // increment in_degree
			in_degree++;
		}
	}
	return in_degree;
	
}
template <typename vertex> std::size_t directed_graph<vertex>::out_degree(const vertex& u) const { 
	// get index position
	int index = get_index(u);
	// initialise out degree and set to zero
	int out_degree = 0;
	// if true in the adjacency matrix means there is an out-degree
	for (bool out : adj_matrix[index]) {
		if (out)
			out_degree++;
	}
	return out_degree;
}
template <typename vertex> std::size_t directed_graph<vertex>::degree(const vertex& u) const {
    // add in and out degree values to find degree
	int degree = in_degree(u) + out_degree(u);
	return degree;
}
template <typename vertex> std::size_t directed_graph<vertex>::num_vertices() const { 
	// returns size of viertices vector
	return vertices.size();
}
template <typename vertex> std::size_t directed_graph<vertex>::num_edges() const { 
   // returns number of edges in graph
	return edge_count;
}
template <typename vertex> std::vector<vertex> directed_graph<vertex>::get_vertices() {
	// returns vertices vector values
	return vertices; 
}
template <typename vertex> std::vector<vertex> directed_graph<vertex>::get_neighbours(const vertex& u) {
	vector<vertex> neighbours;
	for (unsigned int i = 0; i < adj_matrix.size(); i++){
		// if index exists 
	    if(adj_matrix[u][i])  
			// add it to the neighbours vector
		neighbours.push_back(i);
	   }
	return neighbours;
}
template <typename vertex> vertex_iterator<vertex> directed_graph<vertex>::begin() {
	// build beginning vertex iterator
	return vertex_iterator<vertex>(*this,0);
}
template <typename vertex> vertex_iterator<vertex> directed_graph<vertex>::end() {
	// build ending vertex iterator
	return vertex_iterator<vertex>(*this, vertices.size());
}
template <typename vertex> neighbour_iterator<vertex> directed_graph<vertex>::nbegin(const vertex& u) {
	// build beginning neighbour iterator
	return neighbour_iterator<vertex>(*this, u, 0); 
}
template <typename vertex> neighbour_iterator<vertex> directed_graph<vertex>::nend(const vertex& u) { 
    // build ending neighbour iterator
	return neighbour_iterator<vertex>(*this, u, out_degree(u));	
}
template <typename vertex> std::vector<vertex> directed_graph<vertex>::depth_first(const vertex& u) { 
    vector<bool> visited (vertices.size(), false);
	for (unsigned i = 0; i < vertices.size(); i++){
	// Mark all index as not visited  
    visited[i] = false;
	}
	
    stack <vertex> s;
    s.push(u);
    vector<int> ordered;
	// while there are values in the stack
    while(!s.empty()){
	  // initialise the top vertex and remove it from the stack
      int n = s.top();
      s.pop();
	  // if top vertex not visited yet
      if(!visited[n]){
		  // set to visited
		  visited[n] = true;
		  // add vertex to the ordered vector
          ordered.push_back(n);
		  //if vertex contains a neighbour then add it to the stack.
	      for (unsigned i = vertices.size(); i != 0; i--){
				if (adj_matrix[n][i-1]){
					s.push(i-1);
				}
           }
        }
    }
    return ordered;	
 }
template <typename vertex> std::vector<vertex> directed_graph<vertex>::breadth_first(const vertex& u) { 
	 vector<bool> visited (vertices.size(), false);
	for (unsigned i = 0; i < vertices.size(); i++){
		// Mark all index as not visited 
		visited[i] = false;
	}
	// intitalise a queue and add the inputted vertex to the queue
    queue <vertex> q;
    q.push(u);
	vector<int> ordered;
    // while there are values in the queue
    while(!q.empty()){
	   // initialise the top vertex and remove it from the queue
      int n = q.front();
      q.pop();
		// if top vertex not visited yet
      if(!visited[n]){
		  // set to visited
		  visited[n] = true;
		  // add vertex to the ordered vector
		  ordered.push_back(n);
		  //if vertex contains a neighbour then add it to the end of the queue.
          for (unsigned i = 0; i <vertices.size(); i++){
				if (adj_matrix[n][i]){
					q.push(i);
				}
           }
        }
    }
    return ordered;
	 }
template <typename vertex> directed_graph<vertex> directed_graph<vertex>::out_tree(const vertex& u) { 
	directed_graph<vertex> out_graph;
	vector<int> parent(vertices.size());
	vector<int> key(vertices.size(), numeric_limits<int>::max());
	vector<bool> not_out_tree(vertices.size(), false);
	
	
	
	
	
	
	
	return directed_graph<vertex>(); }
template <typename vertex> directed_graph<vertex> directed_graph<vertex>::in_tree(const vertex& u) { return directed_graph<vertex>(); }
template <typename vertex> bool directed_graph<vertex>::reachable(const vertex& u, const vertex& v) const { return false; }


template <typename vertex> vertex_iterator<vertex>::vertex_iterator(const vertex_iterator<vertex>& other) {}
template <typename vertex> vertex_iterator<vertex>::vertex_iterator(const directed_graph<vertex>& graph, std::size_t position) {}
template <typename vertex> vertex_iterator<vertex>::~vertex_iterator() {}
template <typename vertex> vertex_iterator<vertex> vertex_iterator<vertex>::operator=(const vertex_iterator<vertex>& other) { return vertex_iterator<vertex>(directed_graph<vertex>(), std::size_t()); }
template <typename vertex> bool vertex_iterator<vertex>::operator==(const vertex_iterator<vertex>& other) const { return false; }
template <typename vertex> bool vertex_iterator<vertex>::operator!=(const vertex_iterator<vertex>& other) const { return false; }
template <typename vertex> vertex_iterator<vertex> vertex_iterator<vertex>::operator++() { return vertex_iterator<vertex>(directed_graph<vertex>(), std::size_t()); }
template <typename vertex> vertex_iterator<vertex> vertex_iterator<vertex>::operator++(int) { return vertex_iterator<vertex>(directed_graph<vertex>(), std::size_t()); }
template <typename vertex> vertex vertex_iterator<vertex>::operator*() { return vertex(); }
template <typename vertex> vertex* vertex_iterator<vertex>::operator->() { return &vertex(); }

template <typename vertex> neighbour_iterator<vertex>::neighbour_iterator(const neighbour_iterator<vertex>& other) {}
template <typename vertex> neighbour_iterator<vertex>::neighbour_iterator(const directed_graph<vertex>& graph, const vertex& u, std::size_t position) {}
template <typename vertex> neighbour_iterator<vertex>::~neighbour_iterator() {}
template <typename vertex> neighbour_iterator<vertex> neighbour_iterator<vertex>::operator=(const neighbour_iterator<vertex>& other) { return neighbour_iterator<vertex>(directed_graph<vertex>(), vertex(), std::size_t()); }
template <typename vertex> bool neighbour_iterator<vertex>::operator==(const neighbour_iterator<vertex>& other) const { return false; }
template <typename vertex> bool neighbour_iterator<vertex>::operator!=(const neighbour_iterator<vertex>& other) const { return false; }
template <typename vertex> neighbour_iterator<vertex> neighbour_iterator<vertex>::operator++() { return neighbour_iterator<vertex>(directed_graph<vertex>(), vertex(), std::size_t()); }
template <typename vertex> neighbour_iterator<vertex> neighbour_iterator<vertex>::operator++(int) { return neighbour_iterator<vertex>(directed_graph<vertex>(), vertex(), std::size_t()); }		
template <typename vertex> vertex neighbour_iterator<vertex>::operator*() { return vertex(); }
template <typename vertex> vertex* neighbour_iterator<vertex>::operator->() { return &vertex(); }


#endif