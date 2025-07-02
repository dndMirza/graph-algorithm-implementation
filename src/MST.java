import java.util.*;

public class MST {

    // Calculates the total weight of all edges in the graph. Each edge is only counted once for undirected graphs.
    public static double totalEdgeWeight(Graph g) {
        double totalWeight = 0;
        for (int i = 0; i < g.numVertices(); i++) {
            for (int j = i + 1; j < g.numVertices(); j++) {
                if (g.isEdge(i, j)) {
                    totalWeight += g.weight(i, j);
                }
            }
        }
        return totalWeight;
    }

    // Counts the number of connected components in the graph.
    public static int numberOfComponents(Graph g) {
        int numVertices = g.numVertices();
        boolean[] visited = new boolean[numVertices];
        int components = 0;

        for (int i = 0; i < numVertices; i++) {
            if (!visited[i]) {
                components++;
                bfs(g, i, visited); // go through component using bsf
            }
        }
        return components;
    }

    // Uses BFS to explore each unvisited component.
    private static void bfs(Graph g, int start, boolean[] visited) {
        Queue<Integer> queue = new LinkedList<>();
        queue.add(start);
        visited[start] = true;

        while (!queue.isEmpty()) {
            int current = queue.poll();
            for (int neighbor : g.neighbours(current)) {
                if (!visited[neighbor]) {
                    visited[neighbor] = true;
                    queue.add(neighbor);
                }
            }
        }
    }

    public static Graph getRandomGraph(int x, int y, double p) {
        int numVertices = x * y; // Total number of vertices in the grid
        Graph g = new MatrixGraph(numVertices, Graph.UNDIRECTED_GRAPH); // Create an undirected graph
        Random rand = new Random(); // Random generator for edges and weights

        int edgeCount = 0; // Counter for the number of edges added

        for (int i = 0; i < y; i++) {
            for (int j = 0; j < x; j++) {
                int vertex = getVertexId(i, j, x); // Get the vertex ID for this cell
                // Attempt to add an edge to the right neighbor
                if (j + 1 < x && rand.nextDouble() < p) {
                    int rightNeighbor = getVertexId(i, j + 1, x);
                    double weight = rand.nextDouble(); // Assign a random weight to the edge
                    g.addEdge(vertex, rightNeighbor, weight); // Add the edge
                    edgeCount++;
                }
                // Attempt to add an edge to the bottom neighbor
                if (i + 1 < y && rand.nextDouble() < p) {
                    int bottomNeighbor = getVertexId(i + 1, j, x);
                    double weight = rand.nextDouble(); // Assign a random weight to the edge
                    g.addEdge(vertex, bottomNeighbor, weight); // Add the edge
                    edgeCount++;
                }
                int[] coords = getCoordinates(vertex, x);
            }
        }
        return g; // Return the generated graph
    }

    // Converts grid coordinates (i, j) to a vertex ID
    private static int getVertexId(int i, int j, int x) {
        return i * x + j; // Calculate the vertex ID based on row and column
    }

    // Converts a vertex ID to grid coordinates (i, j)
    private static int[] getCoordinates(int vertexId, int x) {
        int i = vertexId / x; // Row number
        int j = vertexId % x; // Column number
        return new int[]{i, j}; // Return the coordinates as an array
    }

    public static class Edge {
        int vertex1; // One end of the edge
        int vertex2; // The other end of the edge
        double weight; // Weight of the edge

        public Edge(int vertex1, int vertex2, double weight) {
            this.vertex1 = vertex1;
            this.vertex2 = vertex2;
            this.weight = weight;
        }
    }

    private static List<Edge> preprocess(Graph g) {
        List<Edge> removedEdges = new ArrayList<>(); // store removed edges
        boolean hasDegreeOne; // check if any degree-1 vertex exists

        do {
            hasDegreeOne = false; // Reset for each iteration
            for (int i = 0; i < g.numVertices(); i++) {
                if (g.degree(i) == 1) { // Check if vertex has degree 1
                    int[] neighbors = g.neighbours(i); // Get neighbors of the vertex
                    if (neighbors.length > 0) {
                        int neighbor = neighbors[0]; // Get the only neighbor
                        double weight = g.weight(i, neighbor); // Get edge weight

                        // Remove edge and add it to the removed list
                        removedEdges.add(new Edge(i, neighbor, weight));
                        g.deleteEdge(i, neighbor); // Delete the edge from the graph
                        hasDegreeOne = true; // Found a degree-1 vertex, so repeat
                    }
                }
            }
        } while (hasDegreeOne); // Repeat while there are degree-1 vertices

        return removedEdges; // Return the list of removed edges
    }

    public static void spanningTree(Graph g) {
        // Preprocess the graph
        List<Edge> preprocessedEdges = preprocess(g);

        // Collect all remaining edges in the graph
        List<Edge> remainingEdges = new ArrayList<>();
        for (int i = 0; i < g.numVertices(); i++) {
            for (int j = i + 1; j < g.numVertices(); j++) {
                if (g.isEdge(i, j)) { // If an edge exists between i and j
                    remainingEdges.add(new Edge(i, j, g.weight(i, j)));
                }
            }
        }
        // Sort the edges in descending order of weight
        remainingEdges.sort(Comparator.comparingDouble((Edge e) -> e.weight).reversed());
        // Get the initial number of components in the graph
        int originalComponents = numberOfComponents(g);
        // Try removing edges while keeping the graph connected
        for (Edge edge : remainingEdges) {
            g.deleteEdge(edge.vertex1, edge.vertex2); // Remove the edge
            if (numberOfComponents(g) > originalComponents) {
                // If removing the edge disconnects the graph, re-add it
                g.addEdge(edge.vertex1, edge.vertex2, edge.weight);
            }
        }
        // Add back the preprocessed edges to the graph
        for (Edge edge : preprocessedEdges) {
            g.addEdge(edge.vertex1, edge.vertex2, edge.weight);
        }
    }

    public static void main(String[] args) {
        System.out.println("Part (a): Graph of Essex");
        Graph essexGraph = GraphOfEssex.getGraph();
        double totalWeightEssex = totalEdgeWeight(essexGraph);
        System.out.println("Total Edge Weight of Graph of Essex: " + totalWeightEssex);

        spanningTree(essexGraph);
        double mstWeightEssex = totalEdgeWeight(essexGraph);
        System.out.println("Total Weight of its MST: " + mstWeightEssex);

        System.out.println("\nPart (b): 1000 Random 10x10 Graphs");
        int numGraphs = 1000;
        double edgeProbability = 2.0 / 3;

        double totalEdgeWeightSum = 0;
        double totalPreprocessedWeightSum = 0;
        double totalWeightAfterPreprocessingSum = 0;
        double totalComponentsBeforeSum = 0;

        for (int i = 0; i < numGraphs; i++) {
            Graph randomGraph = getRandomGraph(10, 10, edgeProbability);

            // Compute total edge weight and components before preprocessing
            double totalEdgeWeightBefore = totalEdgeWeight(randomGraph);
            int componentsBeforePreprocessing = numberOfComponents(randomGraph);

            // Preprocessing
            List<Edge> preprocessedEdges = preprocess(randomGraph);
            double preprocessedWeight = preprocessedEdges.stream().mapToDouble(e -> e.weight).sum();

            // Calculate weight after preprocessing
            double weightAfterPreprocessing = totalEdgeWeightBefore - preprocessedWeight;

            // Add data for averages
            totalEdgeWeightSum += totalEdgeWeightBefore;
            totalPreprocessedWeightSum += preprocessedWeight;
            totalWeightAfterPreprocessingSum += weightAfterPreprocessing;
            totalComponentsBeforeSum += componentsBeforePreprocessing;
        }

        // Averages for part (b)
        System.out.printf("Average Number of Components: %.1f\n", totalComponentsBeforeSum / numGraphs);
        System.out.printf("Average Weight Before Preprocessing: %.6f\n", totalEdgeWeightSum / numGraphs);
        System.out.printf("Average Weight After Preprocessing: %.6f\n", totalWeightAfterPreprocessingSum / numGraphs);
        System.out.println("\nPart (c): 1000 Random 10x10 Graphs for MST Weight");
        double totalMstWeightSum = 0;

        for (int i = 0; i < numGraphs; i++) {
            Graph randomGraph = getRandomGraph(10, 10, edgeProbability);
            spanningTree(randomGraph);
            totalMstWeightSum += totalEdgeWeight(randomGraph);
        }
        // Average MST weight for part (c)
        System.out.printf("Average Total Edge Weight of MSTs: %.6f\n", totalMstWeightSum / numGraphs);
    }
}