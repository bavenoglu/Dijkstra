import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import org.apache.commons.math4.legacy.linear.Array2DRowRealMatrix;
import org.apache.commons.math4.legacy.linear.EigenDecomposition;
import org.apache.commons.math4.legacy.linear.RealMatrix;
import org.apache.commons.math4.legacy.linear.RealVector;

public class Graph {
    public int V;
    public List<List<iPair>> adj;
    //public List<List<List<iPair>>> dijkstraAllPaths;
    public List<Integer> mnec;
    public List<Integer> k = new ArrayList();
    int dijkstraDistance[];
    int bellmanDistance[];
    int dijkstraArray[][];
    //int dijkstraArrayCheck[][];
    double dijkstraNodeDetailsArray[][];
    int dijkstraLoopArray[][];
    double dijkstraArray4Eigen[][];
    int bellmanLoopArray[][];
    int floydWarshallArray[][];
    int dijkstraPath[];
    int dijkstraPathHop[];
    int dijkstraArrayCount;
    int noOutgoingEdge[];
    int[] nodeMinList = new int[V];
    int weight = 10;
    int totalWeight = 0;
    int graphEdgeCount = 0;
    PriorityQueue<iPair> pq;
    PriorityQueue<iPair> pql;
    neo n = new neo();
    Random random = new Random();
    Graph(int V, int maxNodeEdgeCount) {
        this.V = V;
        adj = new ArrayList<>();
        //dijkstraAllPaths = new ArrayList<>();
        mnec = new ArrayList<>();
        dijkstraDistance = new int[V];
        bellmanDistance = new int[V];
        dijkstraArray = new int[V][V+1];
        //dijkstraArrayCheck = new int[V][V];
        dijkstraNodeDetailsArray = new double[V][8];
        dijkstraLoopArray = new int[V][V];
        dijkstraArray4Eigen = new double[V][V];
        bellmanLoopArray = new int[V][V];
        floydWarshallArray = new int[V][V];
        dijkstraPath = new int[V];
        dijkstraPathHop = new int[V];
        dijkstraArrayCount = V*V;
        noOutgoingEdge = new int[V];
        nodeMinList = new int[V];
        pq = new PriorityQueue<>(V, Comparator.comparingInt(o -> o.first));
        pql = new PriorityQueue<>(V, Comparator.comparingInt(o -> o.first));
        for (int i = 0; i < V; i++) {
            for (int j = 0; j < V; j++) {
                if (i == j) {
                    dijkstraArray[i][j] = 0;
                    dijkstraArrayCount--;
                    dijkstraArray4Eigen[i][j] = 0;
                } else {
                    dijkstraArray[i][j] = -1;
                    dijkstraArray4Eigen[i][j] = 0;
                }
            }
            dijkstraArray[i][V] = V-1;
        }
        //Arrays.stream(dijkstraArrayCheck).forEach(a -> Arrays.fill(a, 0));
        n.connect();
        n.deleteAll();
        for (int i = 0; i < V; i++) {
            adj.add(new ArrayList<>());
            //for random graphs
            mnec.add(1 + random.nextInt(maxNodeEdgeCount));
            n.createNode('n' + String.valueOf(i));
        }
    }
    void resetDijkstra(){
        dijkstraArrayCount = V*V;
        for (int i = 0; i < V; i++) {
            for (int j = 0; j < V; j++) {
                if (i == j) {
                    dijkstraArray[i][j] = 0;
                    dijkstraArrayCount--;
                } else {
                    dijkstraArray[i][j] = -1;
                }
            }
            dijkstraArray[i][V] = V-1;
        }
    }
    boolean addEdge(int u, int v, int w) {
        noOutgoingEdge[u]++;
        adj.get(u).add(new iPair(v, w));
        adj.get(v).add(new iPair(u, w));
        totalWeight = totalWeight + w;
        graphEdgeCount = graphEdgeCount + 1;
        n.createRelationship('n' + String.valueOf(u), 'n' + String.valueOf(v), String.valueOf(w), 'r' + String.valueOf(u) + String.valueOf(v) + String.valueOf(w));
        //n.createRelationship('n' + String.valueOf(v), 'n' + String.valueOf(u), String.valueOf(w), 'r' + String.valueOf(v) + String.valueOf(u) + String.valueOf(w));
        dijkstraArray4Eigen[u][v] = 1;
        dijkstraArray4Eigen[v][u] = 1;
        return true;
    }
    void createRandomGraph() {
        //int maxEdge = maxEdgeCount;
        for (int sourceNode = 0; sourceNode < V /*&& maxEdge > 0*/; sourceNode++) {
            int sourceNodeEdgeCount = getRandomEdgeCount(sourceNode);
            int sourceNodeActualEdgeCount = getEdgeCount(sourceNode);
            if (sourceNodeActualEdgeCount < sourceNodeEdgeCount) {
                sourceNodeEdgeCount = sourceNodeEdgeCount - sourceNodeActualEdgeCount;
                //if (maxEdge <= sourceNodeEdgeCount)
                //    sourceNodeEdgeCount = maxEdge;
                int edgeCount = 0;
                ArrayList lst = new ArrayList();
                lst = getRemainingDestinationList(sourceNode);
                while (edgeCount < sourceNodeEdgeCount && lst.size() > 0) {
                        int destinationIndex = random.nextInt(lst.size());
                        int destinationNode = (Integer) lst.get(destinationIndex);
                        lst.remove(destinationIndex);
                        int currentWeight = 1 + random.nextInt(weight);
                        if (addEdge(sourceNode, destinationNode, currentWeight)) {
                            //System.out.println("g.addEdge("+sourceNode+","+destinationNode+","+currentWeight+");");
                            edgeCount++;
                            //maxEdge--;
                        }
                }
            }
        }
    }
    void createBarabasiGraph(){
        createBarabasiDistribution();
        mnec.clear();
        for (int i = 0; i < V; i++) {
            mnec.add(k.get(i));
        }
        createRandomGraph();
    }
    void createBarabasiDistribution(){
        //creates a list of k values based on POWER-LAW DISTRIBUTION
        final double alpha = 2.5;
        Random random = new Random();
        List<Double> r = new ArrayList();
        List<Double> Dk = new ArrayList();
        double Dkcumulative;
        for(int i = 0; i < V; i++) {
            r.add(random.nextDouble());
        }
        Dk.add(1.0);
        Dkcumulative = 1;
        double step = 0.0;
        for(double i = 2; i < V; i++) {
            step = 1/Math.pow(i,alpha);
            Dkcumulative -= step;
            Dk.add(Dkcumulative);
        }
        Dk.add(0.0);
        for(int i = 0; i < V; i++) {
            double Dr = r.get(i);
            for(int j = 0; j < V; j++) {
                if(Dr < Dk.get(j) && Dr >= Dk.get(j+1)) {
                    k.add((int)Math.round(1/Dr));
                    break;
                }
            }
        }
    }
    ArrayList getRemainingDestinationList(int u) {
        boolean exists = false;
        ArrayList lst = new ArrayList();
        for (int i = 0; i < V; i++) {
            exists = false;
            if (i != u) {
                if (adj.get(i).size() < mnec.get(i)) {
                    for (int j = 0; j < adj.get(u).size(); j++)
                        if (i == ((iPair) adj.get(u).get(j)).first)
                            exists = true;
                    for (int j = 0; j < adj.get(i).size(); j++)
                        if (u == ((iPair) adj.get(i).get(j)).first)
                            exists = true;
                } else
                    exists = true;
                if (!exists)
                    lst.add(i);
            }
        }
        return lst;
    }
    int getNoOutgoingEdgeCount() {
        int noOutgoingEdgeCount = 0;
        for (int i = 0; i < V; i++) {
            if(noOutgoingEdge[i] == 0)
                noOutgoingEdgeCount++;
        }
        return noOutgoingEdgeCount;
    }
    int getOrphanedEdgeCount() {
        int totalOrphanedEdgeCount = 0;
        boolean exists = false;
        for (int i = 0; i < V; i++) {
            exists = false;
            for (int j = 0; j < V; j++)
                if (i != j) {
                    for (int k = 0; k < adj.get(j).size(); k++)
                        if (i == ((iPair) adj.get(j).get(k)).first)
                            exists = true;
                }
            if(adj.get(i).size() == 0 && !exists)
                totalOrphanedEdgeCount++;
        }
        return totalOrphanedEdgeCount;
    }
    int getEdgeCount(int u) {
        return adj.get(u).size();
    }
    double getGraphEdgeCount() {
        int totalEdgeCount = 0;
        for (int j = 0; j < V; j++) {
            totalEdgeCount = totalEdgeCount + this.getEdgeCount(j);
        }
        return (totalEdgeCount/2);
    }
    double getGraphMeanEdgeCount() {
        double meanEdgeCount = getGraphEdgeCount()/V;
        return meanEdgeCount;
    }
    double getGraphAverageWight(){
        return totalWeight/V;
    }
    int getNodeTotalWeight(int u){
        int total = 0;
        for (int i = 0; i < adj.get(u).size(); i++) {
            total = total + adj.get(u).get(i).second;
        }
        return total;
    }
    double getNodeStandardDeviation(int u){
        double sum = 0;
        double nodesize = adj.get(u).size();
        double avgweight = getNodeAverageWeight(u);
        for (int j = 0; j < nodesize; j++) {
            double nodeval = adj.get(u).get(j).second;
            sum += Math.pow(nodeval-avgweight,2);
        }
        return Math.sqrt(sum/nodesize);
    }
    double getNodeAverageWeight(int u){
        int nodeSize = adj.get(u).size();
        if(nodeSize == 0) nodeSize = 1;
        return (double) getNodeTotalWeight(u) / nodeSize;
    }
    double getTotalNodeAverageWeight(){
        double nodeAverageWeight = 0;
        for (int j = 0; j < V; j++) {
            nodeAverageWeight += getNodeAverageWeight(j);
        }
        return nodeAverageWeight;
    }
    int getNodeNeighbourCntByEdgeCnt(int u){
        int totalNeighbour = 0;
        for (iPair v : adj.get(u)) {
            totalNeighbour += adj.get(v.first).size() /*- 1*/;
        }
        return totalNeighbour;
    }
    int getTotalNodeNeighbourCntByEdgeCnt(){
        int totalNeighbour = 0;
        for (int j = 0; j < V; j++) {
            totalNeighbour += getNodeNeighbourCntByEdgeCnt(j);
        }
        return totalNeighbour;
    }
    double getNodeNeighboursAvgWeight(int u){
        int totalNeighbourWeight = 0;
        double avgNeighbourWeight = 0;
        //int totalNodeWeight = 0;
        for (iPair v : adj.get(u)) {
            totalNeighbourWeight = 0;
            for (iPair w : adj.get(v.first))
                //if(w.first != u)
                    totalNeighbourWeight += w.second;
            if(adj.get(v.first).size() > 1)
                avgNeighbourWeight += ((double) totalNeighbourWeight / (adj.get(v.first).size()/*-1*/));
            else avgNeighbourWeight += totalNeighbourWeight;
        }
        //avgNeighbourWeight = totalNeighbourWeight / getNodeNeighbourCntByEdgeCnt(u);
        return avgNeighbourWeight;
    }
    double getTotalNodeNeighboursWeight(){
        double totalNeighbourWeight = 0;
        for (int j = 0; j < V; j++) {
            totalNeighbourWeight = totalNeighbourWeight +  getNodeNeighboursAvgWeight(j);
        }
        return totalNeighbourWeight;
    }
    int getRandomEdgeCount(int u) {
        return mnec.get(u);
    }
    void printEdge(int u) {
        for (int i = 0; i < adj.get(u).size(); i++) {
            System.out.println(u + " - " + ((iPair) adj.get(u).get(i)).printiPair());
        }
    }
    static class iPair {
        int first, second;

        iPair(int first, int second) {
            this.first = first;
            this.second = second;
        }
        @Override
        public boolean equals(Object obj) {
            iPair ip = (iPair) obj;
            if (ip.first == this.first && ip.second == this.second)
                return true;
            else
                return false;
        }

        String printiPair() {
            return this.first + " - " + this.second;
        }
    }
    boolean checkNeighbour(int u, int v) {
        for (iPair x : adj.get(u)) {
            if(x.first == v)
                return true;
        }
        return false;
    }
    int getEdgeWeigth(int i, int j) {
        for (iPair x : adj.get(i)) {
            if(x.first == j)
                return x.second;
        }
        return -1;
    }
    int getDijkstraArrayNotFilledCount(){
        int differenceCount = 0;
        for (int i = 0; i < V; i++) {
            for (int j = 0; j < V; j++) {
                if (dijkstraArray[i][j] == -1) {
                    differenceCount++;
                }
            }
        }
        return differenceCount;
    }
    static class PairSize {
         int first;
         int second;

        PairSize(int first, int second) {
            this.first = first;
            this.second = second;
        }
    }
    static class dPairSize {
        int first;
        double second;

        dPairSize(int first, double second) {
            this.first = first;
            this.second = second;
        }
    }
    static class BestVertex {
        int vertexID;
        int arrayCount;
        String minusOneCount;
        int neighbourCount;
        double neighbourWeight;
        int eigenIndex;
        double clusCoefficient;
        int totalWeight;
        double averageWeight;
        double averageNeighbourWeight;

        BestVertex(int vertexID, int arrayCount, String minusOneCount, int neighbourCount,double neighbourWeight, int eigenIndex,double clusCoefficient, int totalWeight, double averageWeight, double averageNeighbourWeight) {
            this.vertexID = vertexID;
            this.arrayCount = arrayCount;
            this.minusOneCount = minusOneCount;
            this.neighbourCount = neighbourCount;
            this.neighbourWeight = neighbourWeight;
            this.eigenIndex = eigenIndex;
            this.clusCoefficient = clusCoefficient;
            this.totalWeight = totalWeight;
            this.averageWeight = averageWeight;
            this.averageNeighbourWeight = averageNeighbourWeight;
        }
    }
    void compareArrays(){
        System.out.println();
        int differenceCount = 0;
        for (int i = 0; i < V; i++) {
            for (int j = 0; j < V; j++) {
                if (dijkstraArray[i][j] != floydWarshallArray[i][j]) {
                    differenceCount++;
                }
            }
        }
        System.out.println("Dijkstra vs Floyd Difference Count: " + differenceCount);
        differenceCount = 0;
        for (int i = 0; i < V; i++) {
            for (int j = 0; j < V; j++) {
                if (dijkstraLoopArray[i][j] != floydWarshallArray[i][j]) {
                    differenceCount++;
                }
            }
        }
        System.out.println("Dijkstra Loop vs Floyd Difference Count: " + differenceCount);
        for (int i = 0; i < V; i++) {
            for (int j = 0; j < V; j++) {
                if (bellmanLoopArray[i][j] != floydWarshallArray[i][j]) {
                    differenceCount++;
                }
            }
        }
        System.out.println("Bellman-Ford vs Floyd Difference Count: " + differenceCount);
    }
    void printArrays(){
        System.out.println("----Dijkstra----------");
        for (int i = 0; i < V; i++) {
            for (int j = 0; j < V; j++) {
                System.out.print(dijkstraArray[i][j] + " ");
            }
            System.out.println();
        }
        System.out.println("----Dijkstra Loop----------");
        for (int i = 0; i < V; i++) {
            for (int j = 0; j < V; j++) {
                System.out.print(dijkstraLoopArray[i][j] + " ");
            }
            System.out.println();
        }
        System.out.println("----Floyd----------");
        for (int i = 0; i < V; i++) {
            for (int j = 0; j < V; j++) {
                System.out.print(floydWarshallArray[i][j] + " ");
            }
            System.out.println();
        }
    }
    void Dijkstra(int src) {
        pq.clear();
        Arrays.fill(dijkstraDistance, 999);
        Arrays.fill(dijkstraPath, 999);
        Arrays.fill(dijkstraPathHop, 0);
        pq.add(new iPair(0, src));
        dijkstraDistance[src] = 0;
        dijkstraPath[src] = 0;
        while (!pq.isEmpty()) {
            int u = pq.poll().second;
            for (iPair v : adj.get(u)) {
                int distSum = dijkstraDistance[u] + v.second;
                int dD = dijkstraDistance[v.first];
                if (dD >= distSum) {
                    dijkstraDistance[v.first] = distSum;
                    int temp = dijkstraArray[v.first][src];
                    dijkstraArray[v.first][src] = distSum;
                    dijkstraArray[src][v.first] = distSum;
                    if (temp == -1) {
                        dijkstraArrayCount--;dijkstraArrayCount--;
                        dijkstraArray[v.first][V] = dijkstraArray[v.first][V] - 1;
                        dijkstraArray[src][V] = dijkstraArray[src][V] - 1;
                    }
                    pq.add(new iPair(distSum, v.first));
                    dijkstraPath[v.first] = u;
                }
            }
        }

        for (int j = 0; j < V; j++) {
            if (dijkstraPath[j] == 999){
                for (int i = 0; i < V; i++)
                    if (dijkstraPath[i] != 999 && dijkstraArray[j][i] != 999){
                        dijkstraArray[j][i] = 999;
                        dijkstraArray[i][j] = 999;
                        dijkstraArrayCount--;dijkstraArrayCount--;
                        dijkstraArray[j][V] = dijkstraArray[j][V] - 1;
                        dijkstraArray[i][V] = dijkstraArray[i][V] - 1;
                    }
                dijkstraPathHop[j] = 0;
            }
            else {
                if(j != src)
                {
                    int k = j;
                    int i = dijkstraPath[j];
                    dijkstraPathHop[j] = 0;
                    if(i != src)
                    {
                        do{
                            do{
                                if(dijkstraArray[k][i] == -1)
                                {
                                    int weightSource = dijkstraDistance[k];
                                    int weightDest = dijkstraDistance[i];
                                    dijkstraArray[k][i] = Math.abs(weightDest - weightSource);
                                    dijkstraArray[i][k] = Math.abs(weightDest - weightSource);
                                    dijkstraArrayCount--;dijkstraArrayCount--;
                                    dijkstraArray[k][V] = dijkstraArray[k][V] - 1;
                                    dijkstraArray[i][V] = dijkstraArray[i][V] - 1;
                                }
                                i = dijkstraPath[i];
                            }while(i != src);
                            k = dijkstraPath[k];
                            i = dijkstraPath[k];
                            dijkstraPathHop[j] = ++dijkstraPathHop[j];
                        } while(k != src && i != src );
                    }
                }
            }
        }
    }
    void Dijkstra2(int src) {
        pq.clear();
        Arrays.fill(dijkstraDistance, 999);
        Arrays.fill(dijkstraPath, 999);
        pq.add(new iPair(0, src));
        dijkstraDistance[src] = 0;
        dijkstraPath[src] = src;
        Queue visitedVertices = new LinkedList();
        while (!pq.isEmpty()) {
            //if (dijkstraArrayCount == 0) break;
            iPair u = pq.poll();
            if (dijkstraArray[src][u.second] == -1) {
                dijkstraArray[src][u.second] = dijkstraDistance[u.second];
                dijkstraArray[u.second][src] = dijkstraDistance[u.second];
                dijkstraArrayCount -= 2;
                dijkstraArray[src][V] = dijkstraArray[src][V] - 1;
                dijkstraArray[u.second][V] = dijkstraArray[u.second][V] - 1;
                int val = dijkstraPath[u.second];
                while (val != src) {
                    if (dijkstraArray[u.second][val] == -1) {
                        int weightSource = dijkstraArray[src][u.second];
                        int weightDest = dijkstraArray[src][val];
                        int dif = Math.abs(weightDest - weightSource);
                        dijkstraArray[u.second][val] = dif;
                        dijkstraArray[val][u.second] = dif;
                        dijkstraArrayCount -= 2;
                        dijkstraArray[u.second][V] = dijkstraArray[u.second][V] - 1;
                        dijkstraArray[val][V] = dijkstraArray[val][V] - 1;
                    }
                    val = dijkstraPath[val];
                }
            }
            visitedVertices.add(u.second);
            for (iPair v : adj.get(u.second)) {
                if (!visitedVertices.contains(v.first)) {
                    int distSum = u.first + v.second;
                    int dD = dijkstraDistance[v.first];
                    if (dD >= distSum) {
                        dijkstraDistance[v.first] = distSum;
                        dijkstraPath[v.first] = u.second;
                        if (!pq.isEmpty()) {
                            iPair existiPair = new iPair(dD, v.first);
                            if (pq.contains(existiPair)) {
                                iPair newiPair = new iPair(dijkstraDistance[v.first], v.first);
                                if (!existiPair.equals(newiPair)) {
                                    pq.remove(existiPair);
                                    pq.add(newiPair);
                                }
                            } else pq.add(new iPair(dijkstraDistance[v.first], v.first));
                        } else
                            pq.add(new iPair(dijkstraDistance[v.first], v.first));
                    }
                }
            }
        }
        if(dijkstraArray[src][V] > 0) {
            for(int i = 0; i < V; i++){
                if(i != src && dijkstraArray[i][V] > 0){
                    if(dijkstraArray[src][i] == -1){
                        dijkstraArray[src][i] = 999;
                        dijkstraArray[i][src] = 999;
                        dijkstraArrayCount -= 2;
                        dijkstraArray[src][V] = dijkstraArray[src][V] - 1;
                        dijkstraArray[i][V] = dijkstraArray[i][V] - 1;
                    }
                    else {
                        for(int j = 0; j < V; j++){
                            if(i != src && j != src && i != j) {
                                if (dijkstraArray[i][j] == -1) {
                                    dijkstraArray[i][j] = 999;
                                    dijkstraArray[j][i] = 999;
                                    dijkstraArrayCount -= 2;
                                    dijkstraArray[i][V] = dijkstraArray[i][V] - 1;
                                    dijkstraArray[j][V] = dijkstraArray[j][V] - 1;
                                }
                            }
                        }
                    }
                    //if(dijkstraArrayCount == 0) break;
                }
            }
        }
    }
    void dijkstraArrayFillOrphanedValues(int orphanedNode){
        for(int i = 0; i < V; i++){
            if(dijkstraArray[orphanedNode][i] == -1){
                dijkstraArray[orphanedNode][i] = 999;
                dijkstraArray[i][orphanedNode] = 999;
                dijkstraArrayCount -= 2;
                dijkstraArray[orphanedNode][V] = dijkstraArray[orphanedNode][V] - 1;
                dijkstraArray[i][V] = dijkstraArray[i][V] - 1;
            }
        }
    };
    //Dijkstra algorithm loops for all vertices
    //Only the last vertex's distances are recorded to dijkstraDistance
    void printDijkstraPathandDistance(){
        for (int j = 0; j < V; j++) {
            System.out.println("Dijkstra Path-Distance: " + j + ": " + dijkstraPath[j] + " "+dijkstraDistance[j]);
        }
    }
    void printDijkstraPath(){
        for (int j = 0; j < V; j++) {
            System.out.print(dijkstraPath[j] + ",");
        }
    }
    void printDijkstraPathHop(){
        System.out.print(" %");
        for (int j = 0; j < V; j++) {
            System.out.print(dijkstraPathHop[j] + ",");
        }
    }
    void printDijkstraDistance(){
        System.out.print(" ||");
        for (int j = 0; j < V; j++) {
            System.out.print(dijkstraDistance[j] + ",");
        }
    }
    void printDijkstraArrayMinusOneCount(){
        //System.out.print(" ~");
        for (int j = 0; j < V; j++) {
            String frmt = String.format("%5.2f",(dijkstraArray[j][V]  * 10000) - (adj.get(j).size()*10000) - getNodeAverageWeight(j));// + (getNodeNeighboursAvgWeight(j) / adj.get(j).size() * (-1)));
            System.out.print(frmt + ", ");
        }
    }
    void dijkstraAllPairsWithLoops(int src) {
        pql.clear();
        Arrays.fill(dijkstraDistance, 999);
        pql.add(new iPair(0, src));
        dijkstraDistance[src] = 0;
        while (!pql.isEmpty()) {
            int u = pql.poll().second;
            for (iPair v : adj.get(u)) {
                int distSum = dijkstraDistance[u] + v.second;
                int dD = dijkstraDistance[v.first];
                if (dD >= distSum) {
                    dijkstraDistance[v.first] = distSum;
                    pql.add(new iPair(distSum, v.first));
                }
            }
        }
        for (int j = 0; j < V; j++) {
            dijkstraLoopArray[src][j] = dijkstraDistance[j];
        }
    }
    void Dijkstra3(int src) {
        pql.clear();
        Arrays.fill(dijkstraDistance, 999);
        Arrays.fill(dijkstraPath, 999);
        pql.add(new iPair(0, src));
        dijkstraDistance[src] = 0;
        while (!pql.isEmpty()) {
            int u = pql.poll().second;
            for (iPair v : adj.get(u)) {
                int distSum = dijkstraDistance[u] + v.second;
                int dD = dijkstraDistance[v.first];
                if (dD >= distSum) {
                    dijkstraDistance[v.first] = distSum;
                    //if(!pql.contains(v.first))
                    pql.add(new iPair(distSum,v.first));
                    dijkstraPath[v.first] = u;
                    int temp = dijkstraArray[v.first][src];
                    dijkstraArray[src][v.first] = dijkstraDistance[v.first];
                    dijkstraArray[v.first][src] = dijkstraDistance[v.first];
                    if (temp == -1) {
                        dijkstraArrayCount -= 2;
                        dijkstraArray[src][V] = dijkstraArray[src][V] - 1;
                        dijkstraArray[v.first][V] = dijkstraArray[v.first][V] - 1;
                    }
                }
            }
        }
        //Recursion
        for(int i = 0; i < V && dijkstraArrayCount > 0; i++){
            int val = dijkstraPath[i];
            while (val != src && val != 999 && dijkstraArrayCount > 0) {
                if (dijkstraArray[i][val] == -1) {
                    int weightSource = dijkstraArray[src][i];
                    int weightDest = dijkstraArray[src][val];
                    int dif = Math.abs(weightDest - weightSource);
                    dijkstraArray[i][val] = dif;
                    dijkstraArray[val][i] = dif;
                    dijkstraArrayCount -= 2;
                    dijkstraArray[i][V] = dijkstraArray[i][V] - 1;
                    dijkstraArray[val][V] = dijkstraArray[val][V] - 1;
                }
                val = dijkstraPath[val];
            }
        }
        //Disconnected
        if(dijkstraArrayCount > 0 && dijkstraArray[src][V] > 0) {
            for(int i = 0; i < V && dijkstraArrayCount > 0 && dijkstraArray[src][V] > 0; i++){
                if(i != src){
                    if(dijkstraArray[src][i] == -1) {
                        dijkstraArray[i][src] = 999;
                        dijkstraArray[src][i] = 999;
                        dijkstraArrayCount -= 2;
                        dijkstraArray[i][V] = dijkstraArray[i][V] - 1;
                        dijkstraArray[src][V] = dijkstraArray[src][V] - 1;
                        for (int j = 0; j < V && dijkstraArrayCount > 0; j++) {
                            if (j != src && i != j) {
                                if (dijkstraArray[src][j] != -1 && dijkstraArray[src][j] != 999) {
                                    dijkstraArray[i][j] = 999;
                                    dijkstraArray[j][i] = 999;
                                    dijkstraArrayCount -= 2;
                                    dijkstraArray[i][V] = dijkstraArray[i][V] - 1;
                                    dijkstraArray[j][V] = dijkstraArray[j][V] - 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    void Dijkstra4Nopq(int src) {
        Queue<Integer> pql = new LinkedList<>();
        Arrays.fill(dijkstraDistance, 999);
        Arrays.fill(dijkstraPath, 999);
        pql.add(src);
        dijkstraDistance[src] = 0;
        while (!pql.isEmpty()) {
            int u = pql.poll();
            for (iPair v : adj.get(u)) {
                int distSum = dijkstraDistance[u] + v.second;
                int dD = dijkstraDistance[v.first];
                if (dD >= distSum) {
                    dijkstraDistance[v.first] = distSum;
                    if(!pql.contains(v.first))
                        pql.add(v.first);
                    dijkstraPath[v.first] = u;
                    int temp = dijkstraArray[v.first][src];
                    dijkstraArray[src][v.first] = dijkstraDistance[v.first];
                    dijkstraArray[v.first][src] = dijkstraDistance[v.first];
                    if (temp == -1) {
                        dijkstraArrayCount -= 2;
                        dijkstraArray[src][V] = dijkstraArray[src][V] - 1;
                        dijkstraArray[v.first][V] = dijkstraArray[v.first][V] - 1;
                    }
                }
            }
        }
        //Recursion
        for(int i = 0; i < V && dijkstraArrayCount > 0; i++){
            int val = dijkstraPath[i];
            while (val != src && val != 999 && dijkstraArrayCount > 0) {
                if (dijkstraArray[i][val] == -1) {
                    int weightSource = dijkstraArray[src][i];
                    int weightDest = dijkstraArray[src][val];
                    int dif = Math.abs(weightDest - weightSource);
                    dijkstraArray[i][val] = dif;
                    dijkstraArray[val][i] = dif;
                    dijkstraArrayCount -= 2;
                    dijkstraArray[i][V] = dijkstraArray[i][V] - 1;
                    dijkstraArray[val][V] = dijkstraArray[val][V] - 1;
                }
                val = dijkstraPath[val];
            }
        }
        //Disconnected
        if(dijkstraArrayCount > 0 && dijkstraArray[src][V] > 0) {
            for(int i = 0; i < V && dijkstraArrayCount > 0 && dijkstraArray[src][V] > 0; i++){
                if(i != src){
                    if(dijkstraArray[src][i] == -1) {
                        dijkstraArray[i][src] = 999;
                        dijkstraArray[src][i] = 999;
                        dijkstraArrayCount -= 2;
                        dijkstraArray[i][V] = dijkstraArray[i][V] - 1;
                        dijkstraArray[src][V] = dijkstraArray[src][V] - 1;
                        for (int j = 0; j < V && dijkstraArrayCount > 0; j++) {
                            if (j != src && i != j) {
                                if (dijkstraArray[src][j] != -1 && dijkstraArray[src][j] != 999) {
                                    dijkstraArray[i][j] = 999;
                                    dijkstraArray[j][i] = 999;
                                    dijkstraArrayCount -= 2;
                                    dijkstraArray[i][V] = dijkstraArray[i][V] - 1;
                                    dijkstraArray[j][V] = dijkstraArray[j][V] - 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    void BellmanFord(int src) {
        Arrays.fill(bellmanDistance, 999);
        bellmanDistance[src] = 0;
        for (int i = 0; i < V; i++) {
            for (int j = 0; j < V; j++) {
                for (iPair p : adj.get(j)) {
                    int u = j, v = p.first, w = p.second;
                    if (bellmanDistance[u] != 999 && bellmanDistance[v] > bellmanDistance[u] + w) {
                        bellmanDistance[v] = bellmanDistance[u] + w;
                    }
                }
            }
        }
        for (int j = 0; j < V; j++) {
            bellmanLoopArray[src][j] = bellmanDistance[j];
        }
        //checks if there exist negative cycles in graph G
       /* for (int i = 0; i < V; i++) {
            for (iPair p : adj.get(i)) {
                int u = i, v = p.first, w = p.second;
                if (distance[u] != Integer.MAX_VALUE && distance[v] > distance[u] + w) {
                    neg = true;
                }
            }
        }*/
    }
    void fillDijkstraArray() {
        double totalNodeNeighbourCntByEdgeCnt = ((double) getTotalNodeNeighbourCntByEdgeCnt());
        if(totalNodeNeighbourCntByEdgeCnt == 0) totalNodeNeighbourCntByEdgeCnt = 1;
        double totalNodeNeighboursWeight = getTotalNodeNeighboursWeight();
        if(totalNodeNeighboursWeight == 0) totalNodeNeighboursWeight = 1;
        for (int j = 0; j < V; j++) {
            double neighbourNodeCount = ((V)*getNodeNeighbourCntByEdgeCnt(j))/totalNodeNeighbourCntByEdgeCnt;
            double edgeSize = (V*adj.get(j).size())/((double) graphEdgeCount*2);
            double neighbourAvgWeight = ((V)*getNodeNeighboursAvgWeight(j))/totalNodeNeighboursWeight;
            double nodeAverageWeight = (V*getNodeAverageWeight(j))/ getTotalNodeAverageWeight();
            dijkstraNodeDetailsArray[j][0] = edgeSize;
            dijkstraNodeDetailsArray[j][1] = neighbourNodeCount;
            dijkstraNodeDetailsArray[j][2] = nodeAverageWeight;
            dijkstraNodeDetailsArray[j][3] = dijkstraNodeDetailsArray[j][0] + dijkstraNodeDetailsArray[j][1];
            dijkstraNodeDetailsArray[j][4] = neighbourAvgWeight;
            dijkstraNodeDetailsArray[j][5] = dijkstraNodeDetailsArray[j][4] + dijkstraNodeDetailsArray[j][3];
            dijkstraNodeDetailsArray[j][6] = dijkstraNodeDetailsArray[j][0] + dijkstraNodeDetailsArray[j][1] - dijkstraNodeDetailsArray[j][2];
            dijkstraNodeDetailsArray[j][7] = dijkstraNodeDetailsArray[j][2] - dijkstraNodeDetailsArray[j][0] - dijkstraNodeDetailsArray[j][1];
//            dijkstraAllPaths.add(j,new ArrayList<>());
//            for (int i = 0; i < V; i++) {
//                    dijkstraAllPaths.get(j).add(new ArrayList<>());
//            }
        }

//        for (int i = 0; i < V; i++) {
//            for (iPair p : adj.get(i)) {
//                int u = i, v = p.first, w = p.second;
//                dijkstraArray[u][v] = w;
//                dijkstraArray[u][V] = dijkstraArray[u][V] - 1;
//                dijkstraArrayCount --;
//            }
//        }
    }
    void fillFloydWarshallArray() {
        //long GstartTime = System.nanoTime();
        for (int i = 0; i < V; i++) {
            for (int j = 0; j < V; j++) {
                if (i == j) {
                    floydWarshallArray[i][j] = 0;
                } else {
                    floydWarshallArray[i][j] = 999;
                }
            }
        }
        //long GendTime = System.nanoTime();
        //System.out.printf("Floyd Fill: %.3f", ((double)(GendTime - GstartTime)) / 1000000);
        //System.out.println();
        for (int i = 0; i < V; i++) {
            for (iPair p : adj.get(i)) {
                int u = i, v = p.first, w = p.second;
                floydWarshallArray[u][v] = w;
            }
        }
    }
    void floydWarshall(int src) {
        for (int k = 0; k < V; k++) {
            // Pick all vertices as source one by one
            for (int i = 0; i < V; i++) {
                // Pick all vertices as destination for the
                // above picked source
                for (int j = 0; j < V; j++) {
                    // If vertex k is on the shortest path
                    // from i to j, then update the value of
                    // dist[i][j]
                    if (floydWarshallArray[i][k] + floydWarshallArray[k][j]
                            < floydWarshallArray[i][j])
                        floydWarshallArray[i][j]
                                = floydWarshallArray[i][k] + floydWarshallArray[k][j];
                }
            }
        }
        /*for (int i = 0; i < V; i++) {
            for (int j = 0; j < V; j++) {
                System.out.print(floydWarshallArray[i][j]+" ");
            }
            System.out.println();
        }*/
    }
    void printToFile() throws IOException {
        FileWriter fw = new FileWriter("GraphOutput" + System.currentTimeMillis() + ".csv", true);
        BufferedWriter bw = new BufferedWriter(fw);
        bw.write("Dijkstra,Bellman,Floyd");
        bw.newLine();
        for (int i = 0; i < V; i++) {
            bw.write(dijkstraDistance[i]+","+bellmanDistance[i]+","+floydWarshallArray[0][i]);
            bw.newLine();
        }
        bw.close();
    }
    List<Integer> calculateEigenVector(){
        RealMatrix matrix = new Array2DRowRealMatrix(dijkstraArray4Eigen);
        EigenDecomposition e = new EigenDecomposition(matrix);
        RealVector arrayEigenVector = e.getEigenvector(0);
        double[] eigenVector = arrayEigenVector.toArray();
        //creates a list of indices based on the ordered values of eigen vector values
        List<Integer> indices = IntStream.range(0, eigenVector.length)
                .boxed().sorted(Comparator.comparingDouble(d -> eigenVector[d])).collect(Collectors.toList());
        return indices.reversed();
        /*
        //System.out.println("***************************************");
        double[] arrayEigenValue = e.getRealEigenvalues();
        for (int i = 0; i < e.getRealEigenvalues().length; i++) {
            System.out.println("eigenValue with index " + i + " " + arrayEigenValue[i]);
            System.out.println(e.getV());
            System.out.println(e.getD());
            System.out.println(e.getVT().getRowVector(0));
            RealVector arrayEigenVector = e.getEigenvector(i);
            for (int j = 0; j < arrayEigenVector.getDimension(); j++) {
                System.out.print(arrayEigenVector.getMaxIndex() + " ");
            }
            System.out.println();
            System.out.println();
        }
        System.out.println("***************************************");
        */
    }
    void calcNodeMinNeighborList(){
        PriorityQueue<PairSize> p = new PriorityQueue<>(Comparator.comparingInt(o -> o.second));
        int previousSecond = 999;
        for(int i = 0; i < V; i++){
            for (iPair v : adj.get(i))
                p.add(new PairSize(v.first,v.second));
            if(!p.isEmpty()){
                PairSize ps = p.poll();
                nodeMinList[ps.first]++;
                previousSecond = ps.second;
                do{
                    if(!p.isEmpty()){
                        ps = p.poll();
                        if(ps.second == previousSecond){
                            nodeMinList[ps.first]++;
                            //previousSecond = ps.second;
                        }
                    }
                }while(ps.second == previousSecond && !p.isEmpty());
                p.clear();
            }
        }
    }
}
