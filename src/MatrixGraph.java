import java.util.Arrays;

/*
 *  A class for representing graphs via adjacency matrices.
 */
public class MatrixGraph extends Graph {
    private final boolean isDirected;

    /*
     *  The adjacency matrix is stored in the 2D array adj. If
     *  adj[x][y] stores Double.NaN, there is no edge (x,y); if
     *  it stores a number, there is an edge with that weight.
     *  Unweighted graphs are represented by giving all edges
     *  weight 1.
     */
    private double[][] adj;

    /*
     *  Create a new graph with the specified number of vertices
     *  (which will be denoted 0, 1, ..., numVertices-1) and
     *  either directed or undirected as appropriate. The graph
     *  is initialized to have no edges, i.e., every entry of
     *  the matrix is Double.NaN.
     */
    public MatrixGraph (int numVertices, boolean isDirected) {
        this.numVertices = numVertices;
        this.isDirected = isDirected;

        adj = new double[numVertices][numVertices];
        for (int i = 0; i < numVertices; i++)
            Arrays.fill(adj[i], Double.NaN);
    }

    /*
     *  Adds the edge (x,y) with the given weight. If the graph
     *  is undirected, we also add (y, x), since this is the same
     *  edge in an undirected graph.
     */
    @Override
    public void addEdge (int x, int y, double weight) {
        if (!isValidEdge (x, y, weight))
            throw new IllegalArgumentException();

        adj[x][y] = weight;
        if (!isDirected)
            adj[y][x] = weight;
    }

    /*
     *  The in-degree of a vertex in a directed graph is the
     *  number of incoming edges it has (i.e., the number of
     *  vertices y for which (y,x) is an edge). For an undir-
     *  ected graph, an incoming edge is just an edge so
     *  degree, in-degree and out-degree are the same thing.
     */
    @Override
    public int inDegree (int x) {
        if (!isValidVertex (x))
            throw new IllegalArgumentException ();

        int count = 0;
        for (int i = 0; i < numVertices; i++)
            if (!Double.isNaN (adj[i][x]))
                count++;
        return count;
    }

    /*
     *  See notes on in-degree, except now we're dealing with
     *  outgoing edges (x,y) instead of incoming edges (y,x).
     */
    @Override
    public int outDegree (int x) {
        if (!isValidVertex (x))
            throw new IllegalArgumentException ();

        int count = 0;
        for (int i = 0; i < numVertices; i++)
            if (!Double.isNaN (adj[x][i]))
                count++;
        return count;
    }

    /*
     *  Delete an edge by setting the relevant entry of the
     *  matrix (two entries if the graph is undirected) to
     *  Double.NaN.
     */
    @Override
    public void deleteEdge (int x, int y) {
        if (!isValidEdge (x, y))
            throw new IllegalArgumentException();

        adj[x][y] = Double.NaN;
        if (!isDirected)
            adj[y][x] = Double.NaN;
    }

    /*
     *  Return whether this graph is directed.
     */
    @Override
    public boolean isDirected () {
        return isDirected;
    }

    /*
     *  Return whether there's an edge between x and y.
     */
    @Override
    public boolean isEdge (int x, int y) {
        return isValidEdge (x,y) && !Double.isNaN (adj[x][y]);
    }

    /*
     *  Return the weight of the edge (x,y), or Double.NaN
     *  if there's no such edge.
     */
    @Override
    public double weight (int x, int y) {
        if (!isValidEdge (x,y))
            throw new IllegalArgumentException ();
        return adj[x][y];
    }

    /*
     *  Return the in-neighbours of a vertex, i.e., the set
     *  of all y such that (y,x) is an edge. In an undirected
     *  graph, in-neighbours, out-neighbours and neighbours
     *  are all the same thing.
     */
    @Override
    public int[] inNeighbours (int x) {
        int[] result = new int[inDegree (x)];
        int i = 0;
        for (int y = 0; y < numVertices; y++)
            if (!Double.isNaN (adj[y][x]))
                result[i++] = y;
        return result;
    }

    /*
     *  Return the out-neighbours of a vertex, i.e., the set
     *  of all y such that (x,y) is an edge.
     */
    @Override
    public int[] outNeighbours (int x) {
        int[] result = new int[outDegree (x)];
        int i = 0;
        for (int y = 0; y < numVertices; y++)
            if (!Double.isNaN (adj[x][y]))
                result[i++] = y;
        return result;
    }

    /*
     *  Test code.
     */
    public static void main (String[] args) {
        Graph g = new MatrixGraph(3, UNDIRECTED_GRAPH);
        g.addEdge(0,1);
        System.out.println("Edge (0,1) exists: " + g.isEdge (0,1));
        System.out.println("Edge (1,0) exists: " + g.isEdge (1,0));
        System.out.println("Edge (1,2) exists: " + g.isEdge (1,2));

        for (int i = 0; i < g.numVertices(); i++) {
            System.out.println("Degree of " + i + ": " + g.degree(i));
            System.out.println("Neighbours of " + i + ": " + Arrays.toString(g.neighbours(i)));
        }

        g.deleteEdge (1,0);
        System.out.println("After deleting (1,0):");
        for (int i = 0; i < g.numVertices(); i++) {
            System.out.println("Degree of " + i + ": " + g.degree(i));
            System.out.println("Neighbours of " + i + ": " + Arrays.toString(g.neighbours(i)));
        }
    }
}
