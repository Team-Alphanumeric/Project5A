
#include <iostream>
#include <limits.h>
#include "d_except.h"
#include <list>
#include <fstream>
#include "d_matrix.h"
#include <queue>
#include <vector>
#include <stack>

#include <boost/graph/adjacency_list.hpp>
#include "heapV.h"

#define LargeValue 99999999

using namespace std;
using namespace boost;

int const NONE = -1;  // Used to represent a node that does not exist

struct VertexProperties;
struct EdgeProperties;

typedef adjacency_list<vecS, vecS, bidirectionalS, VertexProperties, EdgeProperties> Graph;

struct VertexProperties
{
	pair<int, int> cell; // maze cell (x,y) value
	Graph::vertex_descriptor pred;
	bool visited;
	bool marked;
	int weight;
};

// Create a struct to hold properties for each edge
struct EdgeProperties
{
	int weight;
	bool visited;
	bool marked;
};

typedef adjacency_list<vecS, vecS, bidirectionalS, VertexProperties, EdgeProperties> Graph;

void initializeGraph(Graph &g,
	Graph::vertex_descriptor &start,
	Graph::vertex_descriptor &end, ifstream &fin)
	// Initialize g using data from fin.  Set start and end equal
	// to the start and end nodes.
{
	EdgeProperties e;

	int n, i, j;
	int startId, endId;
	fin >> n;
	fin >> startId >> endId;
	Graph::vertex_descriptor v;

	// Add nodes.
	for (int i = 0; i < n; i++)
	{
		v = add_vertex(g);
		if (i == startId)
			start = v;
		if (i == endId)
			end = v;
	}

	while (fin.peek() != '.')
	{
		fin >> i >> j >> e.weight;
		add_edge(i, j, e, g);
	}
}

// Mark all nodes in g as not visited.
void clearVisited(Graph &g)
{
	//loops through all the vertices
	pair<Graph::vertex_iterator, Graph::vertex_iterator> vItrRange = vertices(g);
	for (Graph::vertex_iterator v = vItrRange.first; v != vItrRange.second; ++v)
	{
		//marks each vertice as not visited		
		g[*v].visited = false;
	}
}

void findPathDFSRecursive(Graph &g, Graph::vertex_descriptor node)
{
	//mark the current node as visited
	g[node].visited = true;
	//as long as the end node isn't finsihed. a.k.a the path to the end isn't found
	// then keep searching
	
	//iterate through all the adjacnet nodes of the current node to visit them
	pair<Graph::adjacency_iterator, Graph::adjacency_iterator> vItrRange = adjacent_vertices(node, g);
	for (Graph::adjacency_iterator vItr = vItrRange.first; vItr != vItrRange.second; ++vItr)
	{
		//if the neigboring node has not been visited then go give them a visit
		//dont' ignore any of your neigbors neigbors!!!
		if (!g[*vItr].visited)
		{
			//recursively call the function so that the function works :|
			findPathDFSRecursive(g, *vItr);
		}

	}
}

bool BFSFindCycle(Graph &g, Graph::vertex_descriptor startNode)
{
	//declare all them variables
	Graph::adjacency_iterator vItr;
	queue<Graph::vertex_descriptor> path;
	pair<Graph::adjacency_iterator, Graph::adjacency_iterator> vItrRange;
	Graph::vertex_descriptor v;
	//push the startNode onto the queue because 
	//otherwise the queue would be empty which would 
	//mean the while loop would never execute
	path.push(startNode);
	//set the start node to true because we pushed it onto the queue
	//which means we visited it
	g[startNode].visited = true;
	//contiunue looping until the queue is empty which means that 
	//there was no path found
	while (!path.empty())
	{
		//set v equal to the node currently at the front of the queue
		v = path.front();
		//find the adajcent vertices for node v
		vItrRange = adjacent_vertices(v, g);
		//loop through all the adajcent nodes of v
		for (vItr = vItrRange.first; vItr != vItrRange.second; ++vItr)
		{
			//if an adajcent node has not been visited yet, then visit that node and
			//push it onto the back of the queue.
			//If an adjacent node is already visited, then that means that
			//there is a cycle in that graph because that node has already been visited
			if (g[*vItr].visited)
			{
				return true;
			}
			else 
			{
				//visit the current adajcent node
				g[*vItr].visited = true;
				//push the current adjacent node to the back of the queue
				path.push(*vItr);
				//set the position in the pred vector for the current node to the
				//value of the descriptor for its predecessor. This allows linking of
				//each node so the shortest path to a node can be recreated by finding
				//its predessesor nodes
				g[*vItr].pred = v;
			}
		}
		//if the end node hasn't been found then pop the front of the queue
		//so the process can keep going.
		path.pop();		
	}
	//if the graph traversal finishes without the graph trying to
	//visit an already visited node, that means there are no cycles
	// therefore, return false;
	return false;
}


//create a graph sf that contains a spanning forest on the graph g.
void findSpanningForest(Graph &g1, Graph &sf)
{
	//run BFS on the graph 
}

//Returns true if the graph g is connected. Otherwise false
bool isConnected(Graph &g)
{
	//mark all the nodes in the graph as not visited
	clearVisited(g);
	//find the start node for the graph
	pair<Graph::vertex_iterator, Graph::vertex_iterator> vItrRange = vertices(g);
	//do a graph traversal starting at any node.
	findPathDFSRecursive(g, *(vItrRange.first));
	// loop through all the nodes and if there is a node that wasn't visited, then
	// return false because the graph is not connected, otherwise if all the nodes
	// are visited, then return true because that means the graph is connected
	vItrRange = vertices(g);
	for (Graph::vertex_iterator vItr = vItrRange.first; vItr != vItrRange.second; ++vItr)
	{
		if (!g[*vItr].visited)
		{
			return false;
		}
	}
	return true;

}

//Returns true if the graph g contains a cycle. Otherwise, returns false
bool isCyclic(Graph &g)
{
	//Run Breadth first traversal, and if 
	// the traversal tries to push a node
	// that is already visited onto the queue, 
	pair<Graph::vertex_iterator, Graph::vertex_iterator> vItrRange = vertices(g);
	//loop through all the nodes and start the isCyclic on the first node that is 
	// not visited
	for (Graph::vertex_iterator vItr = vItrRange.first; vItr != vItrRange.second; ++vItr)
	{
		if (!g[*vItr].visited)
		{
			//return the result of BFSFindCycle
			return BFSFindCycle(g, *vItr);
		}
	}
	cout << "already visited" << endl;
	//there are no unvisted nodes, so
	//there are no cycles so return false
	return false;
}


int main()
{
	char x;
	ifstream fin;
	stack <int> moves;
	string fileName;

	// Read the name of the graph from the keyboard or
	// hard code it here for testing.

	fileName = "graph3.txt";

	//   cout << "Enter filename" << endl;
	//   cin >> fileName;

	fin.open(fileName.c_str());
	if (!fin)
	{
		cerr << "Cannot open " << fileName << endl;
		system("pause");
		exit(1);
	}

	try

	{
		cout << "Reading graph" << endl;
		Graph g;

		Graph::vertex_descriptor start, end;

		initializeGraph(g, start, end, fin);
		cout << "Num nodes: " << num_vertices(g) << endl;
		cout << "Num edges: " << num_edges(g) << endl;
		cout << endl;

		// cout << g;
		clearVisited(g);
		bool connected = false;
		bool cyclic = false;

		cout << "Calling isCyclic" << endl;
		cyclic = isCyclic(g); //not created yet

		if (cyclic)
			cout << "Graph contains a cycle" << endl;
		else
			cout << "Graph does not contain a cycle" << endl;

		cout << endl;

		cout << "Calling isConnected" << endl;
		connected = isConnected(g); // not created yet

		if (connected)
			cout << "Graph is connected" << endl;
		else
			cout << "Graph is not connected" << endl;

		cout << endl;
		cout << "Finding spanning forest" << endl;

		// Initialize an empty graph to contain the spanning forest
		Graph sf(num_vertices(g));

		//findSpanningForest(g, sf); //not created yet

		cout << endl;

		cout << "Calling isConnected" << endl;
		connected = isConnected(sf); //not created yet

		if (connected)
			cout << "Graph is connected" << endl;
		else
			cout << "Graph is not connected" << endl;
		cout << endl;

		cout << "Calling isCyclic" << endl;
		cyclic = isCyclic(sf); //not created yet

		if (cyclic)
			cout << "Graph contains a cycle" << endl;
		else
			cout << "Graph does not contain a cycle" << endl;
		cout << endl;
		system("pause");
	}
	catch (indexRangeError &ex)
	{
		cout << ex.what() << endl; exit(1);
	}
	catch (rangeError &ex)
	{
		cout << ex.what() << endl; exit(1);
	}
}
