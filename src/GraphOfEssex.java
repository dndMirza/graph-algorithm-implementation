public class GraphOfEssex {
    public static final int BASILDON   = 0;
    public static final int BRENTWOOD  = 1;
    public static final int BRAINTREE  = 2;
    public static final int CHELMSFORD = 3;
    public static final int CLACTON    = 4;
    public static final int COLCHESTER = 5;
    public static final int HARLOW     = 6;
    public static final int HARWICH    = 7;
    public static final int RAYLEIGH   = 8;
    public static final int SOUTHEND   = 9;
    public static final int WITHAM     = 10;

    public static final String[] TOWN = {"Basildon", "Brentwood", "Braintree", "Chelmsford",
                                         "Clacton", "Colchester", "Harlow", "Harwich", "Rayleigh",
                                         "Southend", "Witham"};

    public static Graph getGraph () {
        MatrixGraph g = new MatrixGraph (11, Graph.UNDIRECTED_GRAPH);

        g.addEdge (BASILDON, BRENTWOOD, 8);
        g.addEdge (BASILDON, CHELMSFORD, 12);
        g.addEdge (BASILDON, RAYLEIGH, 4);
        g.addEdge (BASILDON, SOUTHEND, 10);
        g.addEdge (BRENTWOOD, CHELMSFORD, 11);
        g.addEdge (BRENTWOOD, HARLOW, 14);
        g.addEdge (BRENTWOOD, RAYLEIGH, 13);
        g.addEdge (BRAINTREE, CHELMSFORD, 10);
        g.addEdge (BRAINTREE, COLCHESTER, 15);
        g.addEdge (BRAINTREE, HARLOW, 21);
        g.addEdge (BRAINTREE, WITHAM, 7);
        g.addEdge (CHELMSFORD, HARLOW, 17);
        g.addEdge (CHELMSFORD, RAYLEIGH, 13);
        g.addEdge (CHELMSFORD, SOUTHEND, 16);
        g.addEdge (CHELMSFORD, WITHAM, 9);
        g.addEdge (CLACTON, COLCHESTER, 13);
        g.addEdge (CLACTON, HARWICH, 10);
        g.addEdge (COLCHESTER, HARWICH, 17);
        g.addEdge (COLCHESTER, WITHAM, 14);
        g.addEdge (RAYLEIGH, SOUTHEND, 5);

        return g;
    }
}
