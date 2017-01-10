// Simple recursive dfs from notes in class -- uses adjacency list
public static void dfs(ArrayList[] graph, boolean[] visited, int v) {
    visited[v] = true;

    for (Integer next : ((ArrayList<Integer>)graph)[v])
        if (!visited[next])
            dfs(graph, visited, next);
}
// Simple bfs from class notes -- uses adjacency list
public static int[] bfs(ArrayList[] graph, int v) {
    int n = graph.length;
    int[] distance = new distance[n];
    Arrays.fill(distance, -1);
    visited[n] = true;
    ArrayDeque<Integer> q = new ArrayDeque<Integer>();
    q.offer(v);
    while (q.size() > 0) {
        int cur = q.poll();
        for (Integer next : ((ArrayList<Integer>)graph)[cur]) {
            if (distance[next] == -1) {
                distance[next] = distance[cur]+1;
                q.offer(next);
            }
        }
    }   
    return distance;
}

//------------------------------------------------------------------------------

// Modified dfs for topo sort with cycle detection
public static void dfs(List<Integer>[] graph, boolean[] used, 
    List<Integer> order, int u) {

		used[u] = true;
		for (int v : graph[u]) {			
            if (!used[v]) {
                // Edit this cycle detection if needed (right now toposort
                // continues after a cycle is detected)
                for (int i = 0; i < used.length; i++) {
                    if (used[i] == true && i == v) {
                        System.out.println("Cycle detected. Terminating.");
                        return;
                    }
                }
				
                dfs(graph, used, order, v);
            }
        }
		order.add(u);
}
// Only works with DAGs and doesn't work with cycles.
public static List<Integer> topologicalSort(List<Integer>[] graph) {
	int n = graph.length;
	boolean[] used = new boolean[n];
	List<Integer> order = new ArrayList<>();
	for (int i = 0; i < n; i++)
		if (!used[i])
			dfs(graph, used, order, i);
	Collections.reverse(order);
	return order;
}

//------------------------------------------------------------------------------

// Floyd-Warshall's DP algorithm and sample program
// For all-pairs shortest path.
// Works with negative edges, but not negative cycles.
public class Floyd {   
                                           // V == num vertices
    final static int INF = Integer.MAX_VALUE, V = 4;
 
    public int[][] floydWarshall(int graph[][]) {
        int dist[][] = new int[V][V];
 
        for (int i = 0; i < V; i++)
            for (int j = 0; j < V; j++)
                dist[i][j] = graph[i][j];
 
        for (int k = 0; k < V; k++) {
            for (int i = 0; i < V; i++) {
                for (int j = 0; j < V; j++) {
                    if (dist[i][k] + dist[k][j] < dist[i][j])
                        dist[i][j] = dist[i][k] + dist[k][j];
                }
            }
        }
 
        return dist;
    }
 
    // Example test
    public static void main (String[] args) {
        // Example weighted graph using adj matrix
        int graph[][] = { {0,   5,  INF, 10},
                          {INF, 0,   3, INF},
                          {INF, INF, 0,   1},
                          {INF, INF, INF, 0}
                        };
        Floyd fw = new Floyd();
        int[][] dist = fw.floydWarshall(graph);
    }
}

//------------------------------------------------------------------------------

// Cycle Detection in Directed Graph:
// To detect a back edge, we can keep track of vertices currently in
// recursion stack of function for DFS traversal. If we reach a vertex that
// is already in the recursion stack, then there is a cycle in the graph.
// ^ this also works with undirected graph that has no parallel edges.

// Cycle Detection for Undirected graph with no self-loops:
//    Use the union-find algorithm
// Find: Determine which subset a particular element is in. 
//    This can be used to check if two elements are in the same subset.
// Union: Join two subsets into a single subset.
// Algorithm: 
// Basically, it keeps track of parents. If two connected vertices have the
// same parent, then there's a cycle.

// This code is for an undirected graph with no self-loops...
public class Graph {
    int V, E;    // V == num vertices && E == num edges
    Edge edge[]; 
 
    Graph(int v, int e) {
        V = v;
        E = e;
        edge = new Edge[E];
        for (int i = 0; i < e; i++)
            edge[i] = new Edge();
    }
 
    int find(int parent[], int i) {
        if (parent[i] == -1)
            return i;
        return find(parent, parent[i]);
    }
 
    void union(int parent[], int x, int y) {
        int xset = find(parent, x);
        int yset = find(parent, y);
        parent[xset] = yset;
    }
 
    boolean isCyclic(Graph graph) {
        int parent[] = new int[graph.V];
 
        for (int i = 0; i < graph.V; i++)
            parent[i] = -1;
 
        // meat and potatoes
        for (int i = 0; i < graph.E; i++) {
            int x = graph.find(parent, graph.edge[i].src);
            int y = graph.find(parent, graph.edge[i].dest);
 
            if (x == y)
                return true;
 
            graph.union(parent, x, y);
        }
        return false;
    }
 
    // Example test
    public static void main (String[] args) {
        int numVert = 3, numEdge = 3;
        Graph graph = new Graph(numVert, numEdge);
 
        // Add example edges
        graph.edge[0].src = 0;
        graph.edge[0].dest = 1;
 
        graph.edge[1].src = 1;
        graph.edge[1].dest = 2;
 
        graph.edge[2].src = 0;
        graph.edge[2].dest = 2;
 
        if (graph.isCyclic(graph) == true)
            System.out.println("Graph contains cycle.");
        else
            System.out.println("Graph does not contain cycle.");
    }
}
class Edge {
    int src, dest;
}

//------------------------------------------------------------------------------

// Dijkstra's Algorithm using adj matrix
void dijkstra(int graph[][], int src) {
    int dist[] = new int[V]; // V == num vertices

    Boolean used[] = new Boolean[V];

    for (int i = 0; i < V; i++) {
        dist[i] = Integer.MAX_VALUE;
        used[i] = false;
    }

    dist[src] = 0;

    for (int count = 0; count < V-1; count++) {
        // Pick the minimum distance vertex from the set of unused vertices
        int u = minDistance(dist, used);

        used[u] = true;

        // Update dist value of the adjacent vertices for this vertex
        for (int v = 0; v < V; v++)
            if (!used[v] && graph[u][v] != 0 &&
                    dist[u] != Integer.MAX_VALUE &&
                    dist[u]+graph[u][v] < dist[v])
                dist[v] = dist[u] + graph[u][v];
    }

    return dist;
}

int minDistance(int dist[], Boolean used[]) {
    // Initialize min value
    int min = Integer.MAX_VALUE, minIndex = -1;

    for (int v = 0; v < V; v++)
        if (used[v] == false && dist[v] <= min) {
            min = dist[v];
            minIndex = v;
        }

    return minIndex;
}

//------------------------------------------------------------------------------

// Bellman-Ford is like Dikjstra's, but DOES work with negative edge weights
// still can't work with negative cycles though. 
static final int INF = Integer.MAX_VALUE;

	public static class Edge {
		int v, cost;

		public Edge(int v, int cost) {
			this.v = v;
			this.cost = cost;
		}
	}
    // Right now it returns true if no negative cycle,
    // can easily modify to return the dist array
	public static boolean bellmanFord(List<Edge>[] graph, int s, int[] dist, 
                                                                 int[] pred) {
		Arrays.fill(pred, -1);
		Arrays.fill(dist, INF);
		dist[s] = 0;
		int n = graph.length;
		boolean updated = false;
		for (int step = 0; step < n; step++) {
			updated = false;
			for (int u = 0; u < n; u++) {
				if (dist[u] == INF) continue;
				for (Edge e : graph[u]) {
					if (dist[e.v] > dist[u] + e.cost) {
						dist[e.v] = dist[u] + e.cost;
						dist[e.v] = Math.max(dist[e.v], -INF);
						pred[e.v] = u;
						updated = true;
					}
				}
			}
			if (!updated)
				break;
		}
		// if updated is true then a negative cycle exists
		return updated == false;
	}

//------------------------------------------------------------------------------



// In an undirected graph, a WIDEST path may be found as the path between
// the two vertices in the MAXIMUM spanning tree of the graph, and a minimax
// path may be found as the path between the two vertices in the minimum 
// spanning tree.
//  A MAXIMUM spanning tree is a spanning tree of a weighted graph
//  having maximum weight. It can be computed by negating the weights
//  for each edge and applying Kruskal's algorithm -- won't include here
//  because Ricardo implemented Kruskal's

//  Minimax path problem: asks for the path that minimizes the maximum 
//  weight of any of its edges


