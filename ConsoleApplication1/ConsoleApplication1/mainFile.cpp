#include <iostream>
#include <limits.h>
#include <vector>
#include <list>
#include <fstream>
#include <queue>
#include <boost/graph/adjacency_list.hpp>
#include <stack>
#include <queue>
#include <list>
#include <string>
#include "heapV.h"
#include "d_except.h"
#include "d_matrix.h"

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

// Modifies the new graph to be equivalent to the 
// old graph but with only the marked edges
void edgify(Graph &oldbj, Graph &newgj);

// Mark all nodes in g as not visited.
void setNodeWeights(int weight, Graph &g)
{
	//loops through all the vertices
	pair<Graph::vertex_iterator, Graph::vertex_iterator> vItrRange = vertices(g);
	for (Graph::vertex_iterator v = vItrRange.first; v != vItrRange.second; ++v)
	{
		//marks each vertex as not visited		
		g[*v].weight = weight;
	}
}

// Mark all nodes in g as not visited.
void clearVisited(Graph &g)
{
	//loops through all the vertices
	pair<Graph::vertex_iterator, Graph::vertex_iterator> vItrRange = vertices(g);
	for (Graph::vertex_iterator v = vItrRange.first; v != vItrRange.second; ++v)
	{
		//marks each vertex as not visited		
		g[*v].visited = false;
	}
}

// Unmark all edges in g 
void clearMarkedEdges(Graph &g)
{
	//loop through all the edges
	pair<Graph::edge_iterator, Graph::edge_iterator> eItrRange = edges(g);
	for (Graph::edge_iterator eItr = eItrRange.first; eItr != eItrRange.second; ++eItr)
	{
		g[*eItr].marked = false;
	}
}

// Unmark all vertices in g
void clearMarkedNodes(Graph &g)
{
	//loops through all the vertices
	pair<Graph::vertex_iterator, Graph::vertex_iterator> vItrRange = vertices(g);
	for (Graph::vertex_iterator v = vItrRange.first; v != vItrRange.second; ++v)
	{
		//marks each vertex as unmarked
		g[*v].marked = false;
	}
}

void clearPred(Graph &g)
{
	//loops through all the vertices
	pair<Graph::vertex_iterator, Graph::vertex_iterator> vItrRange = vertices(g);
	for (Graph::vertex_iterator v = vItrRange.first; v != vItrRange.second; ++v)
	{
		//marks each vertex's predecessor as undefined
		g[*v].pred = -1;
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

void findPathBFSSpanningTree(Graph &g, Graph::vertex_descriptor startNode)
{
	//declare all them vuriables
	Graph::adjacency_iterator vItr;
	pair<Graph::adjacency_iterator, Graph::adjacency_iterator> vItrRange;
	queue<Graph::vertex_descriptor> path;
	Graph::vertex_descriptor v;
	//push the startNode onto the queue because 
	//otherwise the queue would be empty which would 
	//mean the while loop would never execute
	path.push(startNode);
	//set the start node to true because we pushed it onto the queue
	//which means we visited it
	g[startNode].visited = true;
	//contiunue looping until the queue is empty 
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
			//push it onto the back of the queue
			if (!g[*vItr].visited)
			{
				//visit the current adajcent node
				g[*vItr].visited = true;	
				//mark the edge in between those nodes as true because
				// there should be an edge between them in the resulting graph
				pair<Graph::edge_descriptor, bool> checkEdge1 = edge(v, *vItr, g);
				pair<Graph::edge_descriptor, bool> checkEdge2 = edge(*vItr, v, g);
				g[checkEdge1.first].marked = true;
				g[checkEdge2.first].marked = true;				
				//push the current adjacent node to the back of the queue
				path.push(*vItr);
				//set the position in the pred vector for the current node to the
				//value of the descriptor for its predecessor. 
				g[*vItr].pred = v;
			}
		}		
		path.pop();		
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
			//If an adjacent node is already visited that is not
			//the predessor, then that means that

			//cout << "The current node is " << v << endl;
			//cout << "The predessor of node " << v << " is " << g[v].pred << endl;			
			//cout << "The adajcent node is " << *vItr << endl << endl;

			//there is a cycle in that graph because that node has already been visited
			if (g[*vItr].visited && (g[v].pred) != *vItr)
			{
				cout << "Cycle found " << endl;
				return true;
			}
			else if (!g[*vItr].visited)
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
void findSpanningForest(Graph &g, Graph &sf)
{
	// clear the given graph
	Graph temp; sf = temp;

	//all visited nodes and marked edges
	clearVisited(g);
	clearMarkedEdges(g);	
	//for every connected graph, find a spanning tree
	pair<Graph::vertex_iterator, Graph::vertex_iterator> vItrRange = vertices(g);
	for (Graph::vertex_iterator vItr = vItrRange.first; vItr != vItrRange.second; ++vItr)
	{
		if (!g[*vItr].visited)
		{
			//run BFS on the graph 
			//marks all the edges that are created for the spanning tree
			findPathBFSSpanningTree(g, *vItr);
		}
	}
	edgify(g,sf);
	return;
}

// adds in the edges for the new graph based on the marked edges from the old graph
void edgify(Graph &oldbj, Graph &newgj)
{
	// variable for storing edges
	pair<Graph::edge_descriptor, bool> newEdge;
	
	// create a graph with the same number of vertices
	newgj = Graph(num_vertices(oldbj));
	
	// add an edge to the new graph (newgj) corresponding to the marked edges in the old graph (g)
	pair<Graph::edge_iterator, Graph::edge_iterator> eItrRange = edges(oldbj);
	for (Graph::edge_iterator eItr = eItrRange.first; eItr != eItrRange.second; ++eItr)
	{
		// Returns the target vertex of edge e.
		Graph::vertex_descriptor targetVer = target(*eItr, oldbj);
		// Returns the source vertex of edge e.
		Graph::vertex_descriptor sourceVer = source(*eItr, oldbj);

		//if the edge is marked in the old graph, add it to the new graph
		if (oldbj[*eItr].marked)
		{
			// create the corresponding edge in the new graph
			newEdge = add_edge(targetVer, sourceVer, newgj);
			if(!newEdge.second)
				throw referenceError("In edgify: edge dupluicated in adding edges to the new graph");
			
			// preserve the weight of the old edge
			newgj[newEdge.first].weight = oldbj[*eItr].weight;
			
			//cout << "Keeping edge between vertices " << targetVer << " and " << sourceVer << endl;
		}
		else
		{
			//cout << "Removing edge between vertices " << targetVer << " and " << sourceVer << endl;
		}
	}
	return;
}


// returns the cumaltive cost of all of the edges in the graph
int edgeCost(Graph &g)
{
	// total cost variable
	int cost = 0;
	
	// Get a pair containing iterators pointing the beginning and end of the list of edges
	pair<Graph::edge_iterator, Graph::edge_iterator> eItrRange = edges(g);

	// Loop over all edges in the graph
	for (Graph::edge_iterator eItr= eItrRange.first; eItr != eItrRange.second; ++eItr)
    {
    	// add the cost of the edge to the total cost
    	cost += g[*eItr].weight;
    }  
    return cost;
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

// Returns the number of connected sub-graphs exist within the graph
int connections(Graph &g)
{
	// mark all the nodes in the graph as not visited
	clearVisited(g);
	
	// find bounds for iterating across all vertices
	pair<Graph::vertex_iterator, Graph::vertex_iterator> vItrRange = vertices(g);

	Graph::vertex_descriptor v;
	
	int numConnected = 0; // nuber of connected sub-graphs
	int i; // iteration variable
	int numV = num_vertices(g); // nubmer of vertices
	bool complete = false; // algorithm completion
	
	// while the algorithm hasn't terminated
	while(!complete)
	{
		i = 0;
		// search for an unvisited node
		for (Graph::vertex_iterator vItr = vItrRange.first; vItr != vItrRange.second; ++vItr)
		{
			// select node if unvisted, and break to traversal
			if (!g[*vItr].visited)
			{
				v = *vItr; break;
			}
		
			// if no node is unvisted after checking all vertices, algorithm is complete
			if(++i == numV)
				complete = true;
		}
		
		// if algorithm is not complete, traverse the sub-graph of the selected node
		if(!complete)
		{
			findPathDFSRecursive(g,v);	// traverse sub-graph
			numConnected++; // increment number of connected sub-graphs traverse
		}
	}
	
	return numConnected;
}

//Returns true if the graph g contains a cycle. Otherwise, returns false
bool isCyclic(Graph &g)
{
	bool result = false;
	//Run Breadth first traversal, and if 
	// the traversal tries to push a node
	// that is already visited onto the queue, 
	clearVisited(g);
	clearPred(g);
	pair<Graph::vertex_iterator, Graph::vertex_iterator> vItrRange = vertices(g);
	//loop through all the nodes and start the isCyclic on the first node that is 
	// not visited
	for (Graph::vertex_iterator vItr = vItrRange.first; vItr != vItrRange.second; ++vItr)
	{
		if (!g[*vItr].visited)
		{
			//accumlate whether or not cycles exists
			//for each connected graph
			result =  result || BFSFindCycle(g, *vItr);
		}
	}
	cout << endl;
	//there are no unvisted nodes, so
	//there are no cycles so return false
	return result;
}

// finds the min of definitively finitely exactly precisely only two numbers
int min(int a, int b) {return ((a<b) ? a : b);}

bool relaxTo0(int &nodeWeight, int edgeWeight)
{
	if(min(nodeWeight,edgeWeight) == edgeWeight)
	{	nodeWeight = edgeWeight; return true;	}
	else
	{	return false; }
}

void msfPrim(Graph &g, Graph &sf)
{
	// Unmark all nodes
	clearVisited(g); clearMarkedEdges(g); clearPred(g);
	
	// set all nodes as unreachable
	setNodeWeights(LargeValue,g);
	
	// Attempt to mark and clear weight of the start node
	try 
	{	g[0].weight = 0; g[0].pred = -1; }
	catch (...)
	{	throw rangeError("In msfPrim: Unable to mark and/or w8 the start node");	}
	
	
	// declare the heap for the nodes that we will use as a priority queue
	heapV<Graph::vertex_descriptor, Graph> nodes;
	
	// create vector for adding all of the nodes
	vector<Graph::vertex_descriptor> allNodes;
		
	// Get a pair containing vertex iterators pointing the beginning and end of the
	// list of nodes
	pair<Graph::vertex_iterator, Graph::vertex_iterator> vItrRange = vertices(g);

	
	// Loop over all nodes in the graph
	for (Graph::vertex_iterator vItr= vItrRange.first; vItr != vItrRange.second; ++vItr)
	{	allNodes.push_back(*vItr);	}
	
	// initilize the heap with all of the nodes
	nodes.initializeMinHeap(allNodes,g);
	
	// create graph with same number of vertices as the previous graph
	sf = Graph(num_vertices(g)); // set sf to a new graph initialized to same # of vertices 
	
	// declare algorithm variables: edge e, vertices u,v
	Graph::edge_descriptor e;
	Graph::vertex_descriptor u,v; 

	// Get a pair containing adjacent vertex iterators pointing the beginning and end of the
	// list of nodes
	pair<Graph::adjacency_iterator, Graph::adjacency_iterator>  vAdjRange;
	
	// Get the edge if it exists and whether it exists
	pair<Graph::edge_descriptor, bool> checkEdge, checkEdge2;
			
	// For(i=1:numNodes-1)
	for(int i=0; i<num_vertices(g); i++)
	{
		// get cheapest source node, u
		u = nodes.extractMinHeapMinimum(g);
		
		// if this node is unreachable, it is from a disconnected sub-graph
		// restart the algorithm by treating it as the next start node
		// by resetting the weight
		if(g[u].weight == LargeValue)
			g[u].weight = 0;
		
		// visit the current node
		g[u].visited = true;
		
		// mark the edge doublet from this node to its predecessor
		if(g[u].pred != -1)
		{
			checkEdge = edge(g[u].pred,u,g);
			checkEdge2 = edge(u,g[u].pred,g);
		
			if(checkEdge.second)
			{
				g[checkEdge.first].marked = true;
			}
			else
			{	throw expressionError("In msfPrim: Unable to find previous edge from source vertex"); }
		
			if(checkEdge2.second)
			{
				g[checkEdge2.first].marked = true;
			}
			else
			{	throw expressionError("In msfPrim: Unable to find previous edge from target vertex"); }
		
		}
		
		//// Find the/a cheapest edge e=(u,v) such that u is marked, v is unmarked and e is cheapest
		
		/// find cheapest edge of cheapest node
		
		// find adjacent vertices
		vAdjRange = adjacent_vertices(u,g);
		
		// Loop over all edges in the graph
		for (Graph::adjacency_iterator vAdj = vAdjRange.first; vAdj != vAdjRange.second; ++vAdj)
		{
			/* set each adjacent node's weight to the corresponding edge 
			 weight (edge from the current node to adj node) if cheaper than current weight
			*/
			
			// get the edge from current node, u, to adjacent node, v
			v = *vAdj;
			
			// only process if v is unvisited
			if( !g[v].visited )
			{
				// get the edge: throw error if it doesn't exist
				checkEdge = edge(u, v, g);
				
				if(checkEdge.second) { e = checkEdge.first;	}
				else { throw expressionError("In msfPrim: unable to find the edge from the current node to the adjacent node");	}
				
				// relax the adjacent vertex versus a 0 current node weight
				bool changed = relaxTo0(g[v].weight,g[e].weight);
				
				// if the weight changed, 
				if(changed)	
				{
					//set the predecessor of the adjacent to the current node
					g[v].pred = u;
								 
					// reset the heap (because v may or may not have been changed)
					nodes.minHeapDecreaseKey(v,g);
				}
			}
		}	
	}
	// place the marked edges into sf
	edgify(g,sf);
	return;
}

void analyzeGraph(Graph &g)
{
		Graph::vertex_descriptor u, v;
		Graph::edge_descriptor e;
		pair<Graph::vertex_iterator, Graph::vertex_iterator> vItrRange = vertices(g);
		pair<Graph::edge_iterator, Graph::edge_iterator> eItrRange = edges(g);

		
		cout << "The graph contains vertices ";
		int i=0; // manual counter
		for (Graph::vertex_iterator vItr= vItrRange.first; vItr != vItrRange.second; ++vItr)
		{
			cout << *vItr;
			cout << ( (i<num_vertices(g) - 1) ? ", " : "." );
			i++; // iteration counter
		} cout << endl;
		
		cout << "The graph contains edges ";

		i=0; // manual counter
		for (Graph::edge_iterator eItr= eItrRange.first; eItr != eItrRange.second; ++eItr)
		{
			cout << *eItr;
			cout << ( (i<num_vertices(g) - 1) ? ", " : "." );
			i++; // iteration counter
		} cout << endl;
 
		bool cyclic = isCyclic(g); 
		bool connected = isConnected(g); 
		int numCon = connections(g);
		int cost = edgeCost(g);
		
		
		cout << "The graph ";
		cout << ((cyclic) ? "contains" : "does not contain");
		cout << " a cycle, and ";
		cout << ((connected) ? "is" : "is not");
		cout << " connected with " << numCon << " sub-graph";
		cout << ((connected) ? "," : "s,");
		cout << " and a total edge cost of " << cost << "." << endl;
		
}


string int2string(int i)
{
	string txt="", temp="";
	while(i != 0)
	{
		txt += char((i%10)+48);
		i = i/10;
	}
	temp.resize(txt.size());
	for(i=0;i<txt.size();++i)
	{ temp[i] = txt[txt.size()-1-i]; }
	return temp;
}



int exploreGraph(int graphNum)
{
	ifstream fin;

	// Read the name of the graph from the keyboard or
	// hard code it here for testing.
	string path = "E:/Users/Thurston Brevett/Documents/Northeastern/Courses/Spring 2015/Algorithms/Project 5/";
	string file = "graph";	
	string fileName = path + file + int2string(graphNum) + ".txt";


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
		// initialize the graph	
		cout << "Reading graph " << graphNum << endl;
		Graph g;

		Graph::vertex_descriptor start, end;

		initializeGraph(g, start, end, fin);
		cout << "Num nodes: " << num_vertices(g) << endl;
		cout << "Num edges: " << num_edges(g) << endl;		
		
		cout << endl << endl << "Analyzing original graph:" << endl;
		analyzeGraph(g);
		
		clearVisited(g);

		// Initialize an empty graph to contain the spanning forest
		Graph sf;
		findSpanningForest(g, sf); 

		cout << endl << endl << "Analyzing default spanning forest graph:" << endl;
		analyzeGraph(sf);
		
		
		// Initialize an empty graph to contain the spanning forest
		Graph msf;
		msfPrim(g, msf); 

		cout << endl << endl << "Analyzing Prim's graph:" << endl;
		analyzeGraph(msf);
		
		
		system("pause");
		cout << endl << endl << endl << endl << endl;
	
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

int main()
{
	for(int i=1; i<=4; i++)
		exploreGraph(i);
}
