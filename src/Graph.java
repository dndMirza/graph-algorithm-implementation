/*
 *  Abstract class for representing graphs.  Implementations should
 *  have a constructor that takes as parameters an int indicating the
 *  number of vertices and a boolean indicating whether the graph is
 *  directed.  Implementations should allow weighted graphs, where an
 *  edge weight may be any double value except Double.NaN, which may
 *  be used to indicate the absence of an edge.
 */
public abstract class Graph {
    /*
     *  Constants for directed vs undirected.
     */
    public static final boolean DIRECTED_GRAPH = true;
    public static final boolean UNDIRECTED_GRAPH = false;

    /*
     *  The number of vertices in the graph.
     */
    protected int numVertices;

    /*
     *  Methods to add an edge from vertex x to vertex y. If no weight
     *  is specified, weight defaults to 1. To create your own implem-
     *  entation of Graphs, you'll need to implement the three-argu-
     *  ment version of the method.
     */
    public abstract void addEdge (int x, int y, double weight);
    public void addEdge (int x, int y) { addEdge (x, y, 1); }

    /*
     *  Delete an edge from the graph, if it exists. Does nothing if
     *  the edge does not exist.
     */
    public abstract void deleteEdge (int x, int y);

    /*
     *  Queries whether this graph is directed.
     */
    public abstract boolean isDirected ();

    /*
     *  Queries whether there is an edge from x to y in this graph.
     */
    public abstract boolean isEdge (int x, int y);

    /*
     *  Queries the weight of the edge (x,y) in this graph. May return
     *  Double.isNaN if the edge does not exist.
     */
    public abstract double weight (int x, int y);

    /*
     *  Methods relating to vertex degrees.  The out-degree of
     *  a vertex x is the number of vertices y such that (x,y) is an
     *  edge, and the in-degree is the number of vertices y such that
     *  (y,x) is an edge.  For an undirected graph, these two quan-
     *  tities are equal and are referred to as just the "degree"
     *  of the vertex.  For a directed graph, the degree of a vertex
     *  is the sum of its in-degree and out-degree.
     */
    public int degree (int x) {
        if (!isValidVertex(x))
            throw new IllegalArgumentException();

        if (isDirected())
            return inDegree (x) + outDegree (x);
        else
            return outDegree (x);
    }

    public abstract int inDegree (int x);
    public abstract int outDegree (int x);

    /*
     *  Query the number of vertices in this graph.
     */
    public int numVertices () {
        return numVertices;
    }

    /*
     *  Methods to return an array of a vertices neighbours in this
     *  graph. An out-neighbour of x is a vertex y such that (x,y) is
     *  an edge; an in-neighbour of x is a vertex such that (y,x) is
     *  and edge.  For an undirected graph, these are the same thing
     *  and are referred to as just the "neighbours" of x.  For a
     *  directed graph, the neighbours are all in-neighbours and out-
     *  neighbours.
     */
    public int[] neighbours (int x) {
        if (!isValidVertex (x))
            throw new IllegalArgumentException ();

        int[] result = new int[degree (x)];

        int[] outNeighbours = outNeighbours(x);
        System.arraycopy(outNeighbours, 0, result, 0, outNeighbours.length);
        if (isDirected()) {
            int[] inNeighbours = inNeighbours(x);
            System.arraycopy(inNeighbours, 0, result, outNeighbours.length,
                    inNeighbours.length);
        }
        return result;
    }
    public abstract int[] inNeighbours (int x);
    public abstract int[] outNeighbours (int x);

    /*
     *  Query whether x is a valid vertex in this graph, whether (x,y)
     *  is a valid possible edge (i.e., whether x and y are both valid
     *  vertices) and whether (x, y, weight) is a valid possible
     *  weighted edge (i.e., whether x and y are both valid vertices
     *  and weight is not Double.NaN). These methods are only intended
     *  to be used to check that arguments to other methods are valid.
     *  They should not be used outside the class.
     */
    protected boolean isValidVertex (int x) {
        return x >= 0 && x < numVertices;
    }

    protected boolean isValidEdge (int x, int y) {
        return x != y && isValidVertex (x) && isValidVertex(y);
    }

    protected boolean isValidEdge (int x, int y, double weight) {
         return isValidEdge (x,y) && !Double.isNaN (weight);
    }
}
