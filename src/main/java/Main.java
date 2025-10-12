import java.awt.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class Main {
    public static <maxEdgeCount> void main(String[] args) throws IOException {
        //System.out.println(Runtime.getRuntime().maxMemory());
        Random random = new Random();
        int V = 4, maxEdgeCount = 0, maxNodeEdgeCount = 0;
        int testloop = 1, avgloop = 0;
        int gtype = 1; // 1-Random, 2-Barabasi
        boolean closeAllPrintStatements = true;
        double result[][] = new double[testloop][26];
        int degreeResult[][] = new int[V][testloop]; for (int[] row : degreeResult) Arrays.fill(row, -1);
        //BigDecimal poissonResult[][] = new BigDecimal[V][testloop];
        double gc=0.0, dc=0.0, bfc=0.0, fwc=0.0, bcwc = 0.0, dapc=0.0, mindc=0.0, maxdc=0.0, negdc=0.0, cofmindc = 0.0, cofmaxdc = 0.0, countmaxdc = 0.0, countD2dc = 0.0, countmax2dc = 0.0, countmindc = 0.0, totalcofavgc = 0.0, cofavgc = 0.0, Totalcofavgc = 0.0, eigendc = 0.0, eigenMindc = 0.0, countmaxneighwdc = 0.0;
        double gt=0.0, dt=0.0, bft=0.0, fwt=0.0, bcwt = 0.0, dapt =0.0, mindt=0.0, maxdt=0.0, negdt=0.0, cofmindt=0.0, cofmaxdt=0.0, countmaxdt=0.0, countD2dt=0.0, countmax2dt=0.0, countmindt=0.0, cofavgt = 0.0, eigendt=0.0, eigenMindt=0.0, countmaxneighwdt=0.0, meedgcnt = 0.0, lnNt = 0.0, lnNlnKt = 0.0, nogedcnt = 0.0, orpht = 0.0;
        int TotaldfalseLoopCounter=0, TotalminfalseLoopCounter = 0, TotalmaxfalseLoopCounter = 0, TotalnegfalseLoopCounter = 0;
        int TotalcofMinfalseLoopCounter = 0, TotalcofMaxfalseLoopCounter = 0, TotalcountMaxfalseLoopCounter = 0,TotalcountD2falseLoopCounter = 0, TotalcountMax2falseLoopCounter = 0, TotalcountMinfalseLoopCounter = 0, TotaleigenfalseLoopCounter = 0, TotaleigenMinfalseLoopCounter = 0, TotalMaxNeighWfalseLoopCounter = 0;
        int dminLoopCounter = 0, maxminLoopCounter = 0, minminLoopCounter = 0, negminLoopCounter = 0, cofminminLoopCounter = 0, cofmaxminLoopCounter = 0, countmaxminLoopCounter = 0,countD2minLoopCounter = 0, countmax2minLoopCounter = 0, countminminLoopCounter = 0, eigenminLoopCounter = 0, eigenminminLoopCounter = 0, maxneighwminLoopCounter = 0;
        int totaldminLoopCounter = 0, totalMaxMinLoopCounter = 0, totalMinMinLoopCounter = 0, totalNegMinLoopCounter = 0, totalCofMinMinLoopCounter = 0, totalCofMaxMinLoopCounter = 0, totalCountMaxMinLoopCounter = 0,totalCountD2MinLoopCounter = 0, totalCountMax2MinLoopCounter = 0, totalCountMinMinLoopCounter = 0, totalEigenMinLoopCounter = 0, totalEigenMinMinLoopCounter = 0, totalMaxNeighWMinLoopCounter = 0;
        int permLoopCounter = 0, TotalPermLoopCounter = 0, TotalEdgeCount = 0;
        int dfalseLoopCounter = 0, minfalseLoopCounter = 0, maxfalseLoopCounter = 0, negfalseLoopCounter = 0, cofMinfalseLoopCounter = 0, cofMaxfalseLoopCounter = 0, countMaxfalseLoopCounter = 0, countD2falseLoopCounter = 0, countMax2falseLoopCounter = 0, countMinfalseLoopCounter = 0, eigenfalseLoopCounter = 0, eigenMinfalseLoopCounter = 0, maxNeighWfalseLoopCounter = 0;
        String fileName = "./results/" + gtype + "_" + V +"_" + testloop + "_" + System.currentTimeMillis() + ".csv";
//        String degreeFileName = "./results/Degree" + System.currentTimeMillis() + ".csv";
        FileWriter fw = new FileWriter(fileName, true);
//        FileWriter fwdegree = new FileWriter(degreeFileName, true);
        BufferedWriter bw = new BufferedWriter(fw);
//        BufferedWriter bwdegree = new BufferedWriter(fwdegree);
        bw.write(",V,maxEdgeCount,maxNodeEdgeCount,GraphT,DAll,DMin,DMax,DNeg,DCofMin,DCofMax,DCountMax,DCountD2,DCountMax2,DCountMin,DEigen,DEigenM,CntMaxW,DLoop,BellmanT,FloydT,Best,EdgeC,MeanEdgeC,lnN,lnN/ln<k>,NoOutgoingC,OrphanedC,ClusCofAvg,PermAvg");
        bw.newLine();
        System.out.println("                                               |  GraphT|  DAll  |   DMin |   DMax |  DNeg  | CofMin | CofMax | CntMax |   CntD2| CntMax2| CntMin |  Eigen | EigenM | CntMaxW|    Loop| Bellmn |  Floyd |    Best|  EdgeC |  MEdgeC|      lnN|  lnN/lnk| NoOutC | OrphanC|ClCofAvg|PermAvg");
        for(int loop=0; loop < testloop; loop++)
        {
            //V = 1 + random.nextInt(1024);
            if( V == 1) V = 2;
            maxEdgeCount = V * (V - 1) / 2; //(int) (2*V*Math.log10(V));  //1 + random.nextInt(V * (V - 1) / 2);
            maxNodeEdgeCount = V - 1; //1 + random.nextInt(V-1);
            long GstartTime = System.nanoTime();
            Graph g = new Graph(V, maxNodeEdgeCount);
            //if(gtype == 1) g.createRandomGraph();
            //else g.createBarabasiGraph();
            //g.addEdge(0,1,9);g.addEdge(0,2,6);g.addEdge(0,4,7);g.addEdge(0,5,4);g.addEdge(1,3,2);g.addEdge(1,4,1);g.addEdge(1,5,8);g.addEdge(3,4,1);
            g.addEdge(0,1,1);g.addEdge(0,2,2);g.addEdge(0,3,3);
            g.addEdge(1,2,3);g.addEdge(1,3,2);
            g.addEdge(2,3,1);
            long GendTime = System.nanoTime();
            gc = ((double) (GendTime - GstartTime)) / 1000000;
            result[loop][0] = gc;

            /*int totalSize = 0;
            for (int i = 0; i < V; i++) {
                    degreeResult[i][loop] = g.adj.get(i).size();
                    totalSize = totalSize + degreeResult[i][loop];
            }
            if(!closeAllPrintStatements) System.out.println("Total size: " + totalSize);*/
            /*Calculates poisson distribution. probability of k =  (-<k>th power of e) * (kth power of <k>) / k!
            closed it for performance
            float avgK = (float) totalSize / V;
            for (int i = 0; i < V; i++) {
                BigDecimal factorial = BigDecimal.valueOf(1);
                for(int s = 1; s <= degreeResult[i][loop]; s++){
                    factorial = factorial.multiply(BigDecimal.valueOf(s));
                }
                BigDecimal lhs = new BigDecimal(Math.exp(-1 * avgK));
                BigDecimal ths = BigDecimal.valueOf(Double.MAX_VALUE) ;
                if(!Double.isInfinite(Math.pow(avgK, degreeResult[i][loop])))
                    ths = BigDecimal.valueOf(Math.pow(avgK, degreeResult[i][loop]));
                poissonResult[i][loop] = lhs.multiply(ths).divide(factorial,50, RoundingMode.HALF_UP);
            }*/

            //Permutation
            perm permutation = new perm(g);
            permutation.generatePerm();
            for (int i = 0; i < permutation.correctPerm.length; i++) {
                if (permutation.correctPerm[i] > 0) {
                    permLoopCounter = i;
                    break;
                }
            }

            int returnVertex = -1;
            class theComparator implements Comparator<Graph.PairSize> {
                @Override
                public int compare(Graph.PairSize o1, Graph.PairSize o2) {
                    if (o1.second < o2.second)
                        return 1;
                    else if (o1.second > o2.second)
                        return -1;
                    return 0;
                }
            }
            class theDoubleComparator implements Comparator<Graph.dPairSize> {
                @Override
                public int compare(Graph.dPairSize o1, Graph.dPairSize o2) {
                    if (o1.second < o2.second)
                        return 1;
                    else if (o1.second > o2.second)
                        return -1;
                    return 0;
                }
            }
            class theBestVertexComparator implements Comparator<Graph.BestVertex> {
                @Override
                public int compare(Graph.BestVertex o1, Graph.BestVertex o2) {
                    if (o1.arrayCount > o2.arrayCount)
                        return 1;
                    else if (o1.arrayCount < o2.arrayCount)
                        return -1;
                    return 0;
                }
            }
            PriorityQueue<Graph.PairSize> dijkstraPriorityQueue = new PriorityQueue<>(V, Comparator.comparingInt(o -> o.second));
            //PriorityQueue<Graph.PairSize> dijkstraOnePriorityQueue = new PriorityQueue<>(V, Comparator.comparingInt(o -> o.second));
            PriorityQueue<Graph.dPairSize> dijkstraDoublePriorityQueue = new PriorityQueue<>(V, Comparator.comparingDouble(o -> o.second));
            PriorityQueue<Graph.dPairSize> dijkstraDoubleMaxPriorityQueue = new PriorityQueue<>(V, new theDoubleComparator());
            PriorityQueue<Graph.dPairSize> dijkstraNegDoubleMaxPriorityQueue = new PriorityQueue<>(V, new theDoubleComparator());
            //PriorityQueue<Graph.PairSize> dijkstraRepeatPriorityQueue = new PriorityQueue<>(V, Comparator.comparingInt(o -> o.second));
            PriorityQueue<Graph.PairSize> dijkstraMaxPriorityQueue = new PriorityQueue<Graph.PairSize>(V, new theComparator());
            PriorityQueue<Graph.PairSize> dijkstraMinPriorityQueue = new PriorityQueue<Graph.PairSize>(V, Comparator.comparingInt(o -> o.second));
            PriorityQueue<Graph.BestVertex> dijkstraBestVertexPriorityQueue = new PriorityQueue<>(V, new theBestVertexComparator());
            PriorityQueue<Graph.dPairSize> dijkstraDoubleDescPQCountMinNeighW = new PriorityQueue<>(V, Comparator.comparingDouble(o -> o.second));
            PriorityQueue<Graph.dPairSize> dijkstraDoubleDescPQCountMaxNeighW = new PriorityQueue<>(V, new theDoubleComparator());
            //Queue<Integer> remainingVertices = new LinkedList<>();

            /********** Clustering Coefficient Computation **********/
            {
                //The clustering coefficient captures the degree to which the neighbors of a given node link to each other.
                // Ci=2Li/ki(ki−1)
                g.resetDijkstra();
                dijkstraDoublePriorityQueue.clear();
                dijkstraDoubleMaxPriorityQueue.clear();
                for (int i = 0; i < V; i++) {
                    int coefficientEdge = 0;
                    for (int j = 0; j < g.adj.get(i).size(); j++) {
                        int nodef = g.adj.get(i).get(j).first;
                        for (int k = j + 1; k < g.adj.get(i).size(); k++) {
                            int nodes = g.adj.get(i).get(k).first;
                            boolean contains = g.adj.get(nodef).stream().anyMatch(x -> x.first == nodes);
                            if (contains)
                                coefficientEdge++;
                        }
                    }
                    if (g.adj.get(i).size() > 1) {
                        dijkstraDoublePriorityQueue.add(new Graph.dPairSize(i, (double) 2 * coefficientEdge / (g.adj.get(i).size() * (g.adj.get(i).size() - 1))));
                        dijkstraDoubleMaxPriorityQueue.add(new Graph.dPairSize(i, (double) 2 * coefficientEdge / (g.adj.get(i).size() * (g.adj.get(i).size() - 1))));
                    } else {
                        dijkstraDoublePriorityQueue.add(new Graph.dPairSize(i, 0));
                        dijkstraDoubleMaxPriorityQueue.add(new Graph.dPairSize(i, 0));
                    }
                }
                Iterator<Graph.dPairSize> iterator = dijkstraDoublePriorityQueue.iterator();
                totalcofavgc = 0.0;
                cofavgc = 0.0;
                while (iterator.hasNext()) {
                    Graph.dPairSize dps = iterator.next();
                    totalcofavgc += dps.second;
                }
                cofavgc = totalcofavgc / V;
                //System.out.println("ClusCofAvg finished!");
            }
            List<Graph.dPairSize> clusCoefList = dijkstraDoublePriorityQueue.stream()
                    .sorted(Comparator.comparingDouble(d -> d.first)).collect(Collectors.toList());
            List<Integer> indices = g.calculateEigenVector();

            //Finding best first vertex
            {
                g.resetDijkstra();
                dijkstraPriorityQueue.clear();
                for (int j = 0; j < V; j++) {
                    Graph.PairSize ps = new Graph.PairSize(j, g.adj.get(j).size());
                    dijkstraPriorityQueue.add(ps);
                }
                while (!dijkstraPriorityQueue.isEmpty()) {
                    g.resetDijkstra();
                    Graph.PairSize returnVertexPair = dijkstraPriorityQueue.poll();
                    g.Dijkstra3(returnVertexPair.first);
                    // node's neighbour's total edge count is divided by node's edge count - N/Edg
                    double neighbourWeight = g.getNodeNeighbourCntByEdgeCnt(returnVertexPair.first);
                    // clustering coefficient of the node - ClCo
                    double vertexclusCoef = clusCoefList.get(returnVertexPair.first).second;
                    String minusOneList = "";
                    for (int j = 0; j < V; j++) {
                        minusOneList = minusOneList + g.dijkstraArray[j][V]+",";
                    }
                    dijkstraBestVertexPriorityQueue.add(new Graph.BestVertex(returnVertexPair.first,g.dijkstraArrayCount,minusOneList,returnVertexPair.second,neighbourWeight,indices.indexOf(returnVertexPair.first),vertexclusCoef,g.getNodeTotalWeight(returnVertexPair.first),g.getNodeAverageWeight(returnVertexPair.first),g.getNodeNeighboursAvgWeight(returnVertexPair.first)));
                }
                if(!closeAllPrintStatements)
                {
                    System.out.println();
                    System.out.println("V\tAcnt\tNeig\tNeigEdg\tEig\tClCo\tTotW\tAvgW\tNeiW");
                }
                while (!dijkstraBestVertexPriorityQueue.isEmpty()){
                    Graph.BestVertex ps  = dijkstraBestVertexPriorityQueue.poll();
                    if(!closeAllPrintStatements)
                    {
                        String bestres = String.format("%3d\t%5d\t%5d\t%4.0f\t%3d\t%4.3f\t%3d\t%8.3f\t%4.3f",ps.vertexID,ps.arrayCount,ps.neighbourCount,ps.neighbourWeight,ps.eigenIndex,ps.clusCoefficient,ps.totalWeight,ps.averageWeight, ps.averageNeighbourWeight);
                        System.out.println(bestres) ;
                    }
                }
            }

            /********** Dijkstra with All Vertex Loops **********/
            // Apply Dijkstra algorithm to all vertices one by one.
            long DAPstartTime = System.nanoTime();
            for (int i = 0; i < V; i++)
                g.dijkstraAllPairsWithLoops(i);
            long DAPendTime = System.nanoTime();
            dapc = ((double) (DAPendTime - DAPstartTime)) / 1000000;
            /********** Bellman-Ford **********/
            // Apply Bellman-Ford algorithm to all vertices one by one.
            long BFstartTime = System.nanoTime();
            for (int i = 0; i < V; i++)
                g.BellmanFord(i);
            long BFendTime = System.nanoTime();
            bfc = ((double) (BFendTime - BFstartTime)) / 1000000;
            /********** Floyd-Warshall **********/
            // Floyd-Warshall algorithm is already an all-pairs algorithm.
            long FWstartTime = System.nanoTime();
            g.fillFloydWarshallArray();
            g.floydWarshall(0);
            long FWFendTime = System.nanoTime();
            fwc = ((double) (FWFendTime - FWstartTime)) / 1000000;
            /********** Best Case **********/
            // The best case by using the permutation
            g.resetDijkstra();
            long BCstartTime = System.nanoTime();
            for(int i : permutation.realList){
                g.Dijkstra3(i);
            }
            long BCendTime = System.nanoTime();
            bcwc = ((double) (BCendTime - BCstartTime)) / 1000000;
            /********** DAll  **********/
            /*{
                g.resetDijkstra();
                for (int j = 0; j < V; j++) {
                    Graph.PairSize ps = new Graph.PairSize(j, g.adj.get(j).size());
                    dijkstraPriorityQueue.add(ps);
                    remainingVertices.add(ps.first);
                }
                dminLoopCounter = 0;
                long DstartTime = System.nanoTime();
                if (g.dijkstraArrayCount > 0) {
                    returnVertex = dijkstraPriorityQueue.poll().first;
                    g.Dijkstra3(returnVertex);
                    remainingVertices.remove(returnVertex);
                    //g.printDijkstraPath();
                    //g.printDijkstraPathHop();
                    //g.printDijkstraDistance();
                    dminLoopCounter++;
                }
                int modNum = 0;
                PriorityQueue<Graph.PairSize> pq = null;
                while (g.dijkstraArrayCount > 0) {
                    if (modNum % 2 == 0) {
                        dijkstraRepeatPriorityQueue.clear();
                        for (int p = 0; p < V; p++) {
                            int finalP = p;
                            boolean contains = true;
                            if (returnVertex == 0 && p == 0)
                                contains = (Arrays.stream(g.dijkstraPath).filter(x -> x == finalP).count() > 1);
                            else if (returnVertex != p)
                                contains = IntStream.of(g.dijkstraPath).anyMatch(x -> x == finalP);
                            if (!contains) {
                                dijkstraRepeatPriorityQueue.add(new Graph.PairSize(p, g.adj.get(p).size()));
                                pq = dijkstraRepeatPriorityQueue;
                            }
                        }
                    } else {
                        dijkstraMaxPriorityQueue.clear();
                        for (int p = 0; p < V; p++)
                            dijkstraMaxPriorityQueue.add(new Graph.PairSize(p, g.dijkstraPathHop[p]));
                        pq = dijkstraMaxPriorityQueue;
                    }
                    while (g.dijkstraArrayCount > 0 && !pq.isEmpty()) {
                        returnVertex = pq.poll().first;
                        if (remainingVertices.contains(returnVertex)) {
                            g.Dijkstra3(returnVertex);
                            remainingVertices.remove(returnVertex);
                            //g.printDijkstraPath();
                            //g.printDijkstraPathHop();
                            //g.printDijkstraDistance();
                            dminLoopCounter++;
                            pq.clear();
                        }
                    }
                    modNum++;
                }
                long DendTime = System.nanoTime();
                dc = ((double) (DendTime - DstartTime)) / 1000000;
                //closed it for performance
                for(int i=0; i < dminLoopCounter; i++){
                    if(permutation.correctPerm[i] > 0){
                        dfalseLoopCounter++;
                        break;
                    }
                }
            }*/
            /********** min - Select the vertex with minimum neighbours **********/
            {
                g.resetDijkstra();
                dijkstraPriorityQueue.clear();
                minfalseLoopCounter = 0;
                minminLoopCounter = 0;
                long minDstartTime = System.nanoTime();
                for (int j = 0; j < V; j++) {
                    Graph.PairSize ps = new Graph.PairSize(j, g.adj.get(j).size());
                    dijkstraPriorityQueue.add(ps);
                }
                while (g.dijkstraArrayCount > 0) {
                    returnVertex = dijkstraPriorityQueue.poll().first;
                    g.Dijkstra3(returnVertex);
                    //g.printDijkstraPath();
                    //g.printDijkstraPathHop();
                    //g.printDijkstraDistance();
                    minminLoopCounter++;
                }
                long minDendTime = System.nanoTime();
                mindc = ((double) (minDendTime - minDstartTime)) / 1000000;
                for (int i = 0; i < minminLoopCounter; i++) {
                    if (permutation.correctPerm[i] > 0) {
                        minfalseLoopCounter++;
                        break;
                    }
                }
            }
            /********** max - Select the vertex with maximum neighbours **********/
            {
                g.resetDijkstra();
                dijkstraMaxPriorityQueue.clear();
                maxminLoopCounter = 0;
                maxfalseLoopCounter = 0;
                long maxDstartTime = System.nanoTime();
                for (int j = 0; j < V; j++) {
                    Graph.PairSize ps = new Graph.PairSize(j, g.adj.get(j).size());
                    dijkstraMaxPriorityQueue.add(ps);
                }
                while (g.dijkstraArrayCount > 0) {
                    returnVertex = dijkstraMaxPriorityQueue.poll().first;
                    g.Dijkstra3(returnVertex);
                    //g.printDijkstraPath();
                    //g.printDijkstraPathHop();
                    //g.printDijkstraDistance();
                    maxminLoopCounter++;
                }
                long maxDendTime = System.nanoTime();
                maxdc = ((double) (maxDendTime - maxDstartTime)) / 1000000;
                for(int i=0; i < maxminLoopCounter; i++){
                    if(permutation.correctPerm[i] > 0){
                        maxfalseLoopCounter++;
                        break;
                    }
                }
            }
            /********** neigAvg - Select the vertex with max of (total edge count of neighbours / edge count) **********/
            {
                g.resetDijkstra();
                negminLoopCounter = 0;
                negfalseLoopCounter = 0;
                long negDstartTime = System.nanoTime();
                for (int j = 0; j < V; j++) {
                    int totalNeighbour = 0;
                    for (Graph.iPair v : g.adj.get(j)) {
                        totalNeighbour += g.adj.get(v.first).size();
                    }
                    Graph.dPairSize ps = new Graph.dPairSize(j, (double) totalNeighbour / g.adj.get(j).size());
                    dijkstraNegDoubleMaxPriorityQueue.add(ps);
                }
                while (g.dijkstraArrayCount > 0) {
                    returnVertex = dijkstraNegDoubleMaxPriorityQueue.poll().first;
                    g.Dijkstra3(returnVertex);
                    //g.printDijkstraPath();
                    //g.printDijkstraPathHop();
                    //g.printDijkstraDistance();
                    negminLoopCounter++;
                }
                long negDendTime = System.nanoTime();
                negdc = ((double) (negDendTime - negDstartTime)) / 1000000;
                for(int i=0; i < negminLoopCounter; i++){
                    if(permutation.correctPerm[i] > 0){
                        negfalseLoopCounter++;
                        break;
                    }
                }
            }
            /********** cofMin - Select the vertex with minimum clustering coefficient **********/
            {
                g.resetDijkstra();
                cofminminLoopCounter = 0;
                cofMinfalseLoopCounter = 0;
                long cofMinDstartTime = System.nanoTime();
                while (g.dijkstraArrayCount > 0) {
                    returnVertex = dijkstraDoublePriorityQueue.poll().first;
                    g.Dijkstra3(returnVertex);
                    //g.printDijkstraPath();
                    //g.printDijkstraPathHop();
                    //g.printDijkstraDistance();
                    cofminminLoopCounter++;
                }
                long cofMinDendTime = System.nanoTime();
                cofmindc = ((double) (cofMinDendTime - cofMinDstartTime)) / 1000000;
                for(int i=0; i < cofminminLoopCounter; i++){
                    if(permutation.correctPerm[i] > 0){
                        cofMinfalseLoopCounter++;
                        break;
                    }
                }
            }
            /********** cofMax - Select the vertex with maximum clustering coefficient **********/
            {
                g.resetDijkstra();
                cofmaxminLoopCounter = 0;
                cofMaxfalseLoopCounter = 0;
                long cofMaxDstartTime = System.nanoTime();
                while (g.dijkstraArrayCount > 0) {
                    returnVertex = dijkstraDoubleMaxPriorityQueue.poll().first;
                    g.Dijkstra3(returnVertex);
                    //g.printDijkstraPath();
                    //g.printDijkstraPathHop();
                    //g.printDijkstraDistance();
                    cofmaxminLoopCounter++;
                }
                long cofMaxDendTime = System.nanoTime();
                cofmaxdc = ((double) (cofMaxDendTime - cofMaxDstartTime)) / 1000000;
                for(int i=0; i < cofmaxminLoopCounter; i++){
                    if(permutation.correctPerm[i] > 0){
                        cofMaxfalseLoopCounter++;
                        break;
                    }
                }
            }
            /********** countWeight - Select the vertex with maximum -1 values and minimum edge count **********/
            {
                g.resetDijkstra();
                g.fillDijkstraArray();
                dijkstraPriorityQueue.clear();
                countmaxminLoopCounter = 0;
                countMaxfalseLoopCounter = 0;
                long countMaxDstartTime = System.nanoTime();
                for (int j = 0; j < V; j++) {
                        Graph.PairSize ps = new Graph.PairSize(j, g.adj.get(j).size()*10000 - (int)g.dijkstraNodeDetailsArray[j][2]); //*1000+g.getNodeTotalWeight(j)
                        dijkstraPriorityQueue.add(ps);
                }
                if (g.dijkstraArrayCount > 0) {
                    returnVertex = dijkstraPriorityQueue.poll().first;
                    g.Dijkstra3(returnVertex);
                    //g.printDijkstraPath();
                    //g.printDijkstraPathHop();
                    //g.printDijkstraDistance();
                    countmaxminLoopCounter++;
                }
                while (g.dijkstraArrayCount > 0) {
                    dijkstraMaxPriorityQueue.clear();
                    for (int i = 0; i < V; i++) {
                        if (g.dijkstraArray[i][V] != 0) {
                            int second = g.dijkstraArray[i][V] * 10000 + (g.adj.get(i).size() * (-1));
                            dijkstraMaxPriorityQueue.add(new Graph.PairSize(i, second)); //
                        }
                    }
                    returnVertex = dijkstraMaxPriorityQueue.poll().first;
                    g.Dijkstra3(returnVertex);
                    //g.printDijkstraPath();
                    //g.printDijkstraPathHop();
                    //g.printDijkstraDistance();
                    countmaxminLoopCounter++;
                }
                long countMaxDendTime = System.nanoTime();
                countmaxdc = ((double) (countMaxDendTime - countMaxDstartTime)) / 1000000;
                for(int i=0; i < countmaxminLoopCounter; i++){
                    if(permutation.correctPerm[i] > 0){
                        countMaxfalseLoopCounter++;
                        break;
                    }
                }
            }
            /********** countDijkstra2 - modification within Dijkstra **********/
            /*{
                g.resetDijkstra();
                dijkstraDoubleDescPQCountMinNeighW.clear();
                countD2minLoopCounter = 0;
                g.fillDijkstraArray();
                long countD2DstartTime = System.nanoTime();
                for (int i = 0; i < V; i++) {
                    if(g.adj.get(i).size() > 0) {
                        Graph.dPairSize ps = new Graph.dPairSize(i,g.dijkstraNodeDetailsArray[i][6]);
                        dijkstraDoubleDescPQCountMinNeighW.add(ps);
                    }
                    else g.dijkstraArrayFillOrphanedValues(i);
                }
                returnVertex = dijkstraDoubleDescPQCountMinNeighW.poll().first;
                g.Dijkstra(returnVertex);
                countD2minLoopCounter++;
                while (g.dijkstraArrayCount > 0) {
                    dijkstraDoubleDescPQCountMaxNeighW.clear();
                    for (int i = 0; i < V; i++) {
                        if (g.dijkstraArray[i][V] > 0 && g.adj.get(i).size() > 0) {
                            dijkstraDoubleDescPQCountMaxNeighW.add(new Graph.dPairSize(i, (g.dijkstraArray[i][V]*10000) + g.dijkstraNodeDetailsArray[i][7]));
                        }
                    }
                    returnVertex = dijkstraDoubleDescPQCountMaxNeighW.poll().first;
                    g.Dijkstra(returnVertex);
                    countD2minLoopCounter++;
                }
                long countD2DendTime = System.nanoTime();
                countD2dc = ((double) (countD2DendTime - countD2DstartTime)) / 1000000;
                for(int i=0; i < countD2minLoopCounter; i++){
                    if(permutation.correctPerm[i] > 0){
                        countD2falseLoopCounter++;
                        break;
                    }
                }
                g.compareArrays();
            }*/
            /********** VerDetailPQ - with PriorityQueue **********/
            {
                g.resetDijkstra();
                dijkstraDoubleDescPQCountMinNeighW.clear();
                countmax2minLoopCounter = 0;
                countMax2falseLoopCounter = 0;
                g.fillDijkstraArray();
                long countMax2DstartTime = System.nanoTime();
                for (int i = 0; i < V; i++) {
                    if(g.adj.get(i).size() > 0) {
                        Graph.dPairSize ps = new Graph.dPairSize(i,g.dijkstraNodeDetailsArray[i][6]);
                        dijkstraDoubleDescPQCountMinNeighW.add(ps);
                    }
                    else g.dijkstraArrayFillOrphanedValues(i);
                }
                returnVertex = dijkstraDoubleDescPQCountMinNeighW.poll().first;
                g.Dijkstra3(returnVertex);
                countmax2minLoopCounter++;
                while (g.dijkstraArrayCount > 0) {
                    dijkstraDoubleDescPQCountMaxNeighW.clear();
                    for (int i = 0; i < V; i++) {
                        if (g.dijkstraArray[i][V] > 0 && g.adj.get(i).size() > 0) {
                            dijkstraDoubleDescPQCountMaxNeighW.add(new Graph.dPairSize(i, (g.dijkstraArray[i][V]*10000) + g.dijkstraNodeDetailsArray[i][7]));
                        }
                    }
                    returnVertex = dijkstraDoubleDescPQCountMaxNeighW.poll().first;
                    g.Dijkstra3(returnVertex);
                    countmax2minLoopCounter++;
                }
                long countMax2DendTime = System.nanoTime();
                countmax2dc = ((double) (countMax2DendTime - countMax2DstartTime)) / 1000000;
                for(int i=0; i < countmax2minLoopCounter; i++){
                    if(permutation.correctPerm[i] > 0){
                        countMax2falseLoopCounter++;
                        break;
                    }
                }
                //g.compareArrays();
            }
            /********** countMin - Select the vertex with maximum -1 values and maximum edge count **********/
            {
                g.resetDijkstra();
                dijkstraPriorityQueue.clear();
                countminminLoopCounter = 0;
                countMinfalseLoopCounter = 0;
                long countMinDstartTime = System.nanoTime();
                for (int j = 0; j < V; j++) {
                    Graph.PairSize ps = new Graph.PairSize(j, g.adj.get(j).size());
                    dijkstraPriorityQueue.add(ps);
                }
                if (g.dijkstraArrayCount > 0) {
                    returnVertex = dijkstraPriorityQueue.poll().first;
                    g.Dijkstra3(returnVertex);
                    //g.printDijkstraPath();
                    //g.printDijkstraPathHop();
                    //g.printDijkstraDistance();
                    countminminLoopCounter++;
                }
                while (g.dijkstraArrayCount > 0) {
                    dijkstraMinPriorityQueue.clear();
                    for (int i = 0; i < V; i++) {
                        if (g.dijkstraArray[i][V] != 0)
                            dijkstraMinPriorityQueue.add(new Graph.PairSize(i, (g.dijkstraArray[i][V] * (-10000)) - g.adj.get(i).size()));
                    }
                    returnVertex = dijkstraMinPriorityQueue.poll().first;
                    g.Dijkstra3(returnVertex);
                    //g.printDijkstraPath();
                    //g.printDijkstraPathHop();
                    //g.printDijkstraDistance();
                    countminminLoopCounter++;
                }
                long countMinDendTime = System.nanoTime();
                countmindc = ((double) (countMinDendTime - countMinDstartTime)) / 1000000;
                for(int i=0; i < countminminLoopCounter; i++){
                    if(permutation.correctPerm[i] > 0){
                        countMinfalseLoopCounter++;
                        break;
                    }
                }
            }
            /********** eigenMax - execute Dijkstra based on eigenvector **********/
            {
                g.resetDijkstra();
                dijkstraMaxPriorityQueue.clear();
                eigenminLoopCounter = 0;
                eigenfalseLoopCounter = 0;
                long eigenDstartTime = System.nanoTime();
                for (Integer element : indices) {
                    Graph.PairSize ps = new Graph.PairSize(element, g.adj.get(element).size());
                    dijkstraMaxPriorityQueue.add(ps);
                }
                while (g.dijkstraArrayCount > 0) {
                    returnVertex = dijkstraMaxPriorityQueue.poll().first;
                    g.Dijkstra3(returnVertex);
                    //g.printDijkstraPath();
                    //g.printDijkstraPathHop();
                    //g.printDijkstraDistance();
                    eigenminLoopCounter++;
                }
                long eigenDendTime = System.nanoTime();
                eigendc = ((double) (eigenDendTime - eigenDstartTime)) / 1000000;
                for (int i = 0; i < eigenminLoopCounter; i++) {
                    if (permutation.correctPerm[i] > 0) {
                        eigenfalseLoopCounter++;
                        break;
                    }
                }
            }
            /********** eigenMin - execute Dijkstra based on eigenvector **********/
            {
                g.resetDijkstra();
                dijkstraPriorityQueue.clear();
                eigenminminLoopCounter = 0;
                eigenMinfalseLoopCounter = 0;
                long eigenMinDstartTime = System.nanoTime();
                for (Integer element : indices) {
                    Graph.PairSize ps = new Graph.PairSize(element, g.adj.get(element).size());
                    dijkstraPriorityQueue.add(ps);
                }
                while (g.dijkstraArrayCount > 0) {
                    returnVertex = dijkstraPriorityQueue.poll().first;
                    g.Dijkstra3(returnVertex);
                    //g.printDijkstraPath();
                    //g.printDijkstraPathHop();
                    //g.printDijkstraDistance();
                    eigenminminLoopCounter++;
                }
                long eigenMinDendTime = System.nanoTime();
                eigenMindc = ((double) (eigenMinDendTime - eigenMinDstartTime)) / 1000000;
                for (int i = 0; i < eigenminminLoopCounter; i++) {
                    if (permutation.correctPerm[i] > 0) {
                        eigenMinfalseLoopCounter++;
                        break;
                    }
                }
            }
            /********** VerDetailNoPQ   **********/
            {
                g.resetDijkstra();
                dijkstraDoubleDescPQCountMinNeighW.clear();
                maxneighwminLoopCounter = 0;
                maxNeighWfalseLoopCounter = 0;
                g.fillDijkstraArray();
                long countMaxNeighWDstartTime = System.nanoTime();
                for (int j = 0; j < V; j++) {
                    if(g.adj.get(j).size() > 0) {
                        Graph.dPairSize ps = new Graph.dPairSize(j, g.dijkstraNodeDetailsArray[j][6]);
                        dijkstraDoubleDescPQCountMinNeighW.add(ps);
                    }
                    else g.dijkstraArrayFillOrphanedValues(j);
                }
                returnVertex = dijkstraDoubleDescPQCountMinNeighW.poll().first;
                g.Dijkstra4Nopq(returnVertex);
                maxneighwminLoopCounter++;
                while (g.dijkstraArrayCount > 0) {
                    dijkstraDoubleDescPQCountMaxNeighW.clear();
                    for (int i = 0; i < V; i++) {
                        if (g.dijkstraArray[i][V] > 0 && g.adj.get(i).size() > 0) {
                            dijkstraDoubleDescPQCountMaxNeighW.add(new Graph.dPairSize(i, (g.dijkstraArray[i][V]*10000) + g.dijkstraNodeDetailsArray[i][7]));
                        }
                    }
                    returnVertex = dijkstraDoubleDescPQCountMaxNeighW.poll().first;
                    g.Dijkstra4Nopq(returnVertex);
                    maxneighwminLoopCounter++;
                }
                long countMaxNeighWDendTime = System.nanoTime();
                countmaxneighwdc = ((double) (countMaxNeighWDendTime - countMaxNeighWDstartTime)) / 1000000;
                for(int i=0; i < maxneighwminLoopCounter; i++){
                    if(permutation.correctPerm[i] > 0){
                        maxNeighWfalseLoopCounter++;
                        break;
                    }
                }
                //g.compareArrays();
            }
             result[loop][1] = dc;
            result[loop][2] = mindc;
            result[loop][3] = maxdc;
            result[loop][4] = negdc;
            result[loop][5] = cofmindc;
            result[loop][6] = cofmaxdc;
            result[loop][7] = countmaxdc;
            result[loop][8] = countD2dc;
            result[loop][9] = countmax2dc;
            result[loop][10] = countmindc;
            result[loop][11] = eigendc;
            result[loop][12] = eigenMindc;
            result[loop][13] = countmaxneighwdc;
            result[loop][14] = dapc;
            result[loop][15] = bfc;
            result[loop][16] = fwc;
            result[loop][17] = bcwc;
            result[loop][18] = g.getGraphEdgeCount();
            result[loop][19] = g.getGraphMeanEdgeCount();
            result[loop][20] = Math.log(V);
            result[loop][21] = result[loop][18]/Math.log(result[loop][17]);
            result[loop][22] = g.getNoOutgoingEdgeCount();
            result[loop][23] = g.getOrphanedEdgeCount();
            result[loop][24] = cofavgc;
            result[loop][25] = permLoopCounter;
            if(loop > avgloop-1) {
                gt = gt + gc;
                dt = dt + dc;
                mindt = mindt + mindc;
                maxdt = maxdt + maxdc;
                negdt = negdt + negdc;
                cofmindt = cofmindt + cofmindc;
                cofmaxdt = cofmaxdt + cofmaxdc;
                countmaxdt = countmaxdt + countmaxdc;
                countD2dt = countD2dt + countD2dc;
                countmax2dt = countmax2dt + countmax2dc;
                countmindt = countmindt + countmindc;
                eigendt = eigendt + eigendc;
                eigenMindt = eigenMindt + eigenMindc;
                countmaxneighwdt = countmaxneighwdt + countmaxneighwdc;
                dapt = dapt + dapc;
                bft = bft + bfc;
                fwt = fwt + fwc;
                bcwt = bcwt + bcwc;
                TotalEdgeCount += g.graphEdgeCount;
                meedgcnt += result[loop][19];
                lnNt += result[loop][20];
                lnNlnKt += result[loop][21];
                nogedcnt += result[loop][22];
                orpht += result[loop][23];
                cofavgt = cofavgt + cofavgc;
                TotalPermLoopCounter += permLoopCounter;

                TotaldfalseLoopCounter += dfalseLoopCounter;
                TotalminfalseLoopCounter += minfalseLoopCounter;
                TotalmaxfalseLoopCounter += maxfalseLoopCounter;
                TotalnegfalseLoopCounter += negfalseLoopCounter;
                TotalcofMinfalseLoopCounter += cofMinfalseLoopCounter;
                TotalcofMaxfalseLoopCounter += cofMaxfalseLoopCounter;
                TotalcountMaxfalseLoopCounter += countMaxfalseLoopCounter;
                TotalcountD2falseLoopCounter += countD2falseLoopCounter;
                TotalcountMax2falseLoopCounter += countMax2falseLoopCounter;
                TotalcountMinfalseLoopCounter += countMinfalseLoopCounter;
                TotaleigenfalseLoopCounter += eigenfalseLoopCounter;
                TotaleigenMinfalseLoopCounter += eigenMinfalseLoopCounter;
                TotalMaxNeighWfalseLoopCounter += maxNeighWfalseLoopCounter;
                totaldminLoopCounter += dminLoopCounter;
                totalMinMinLoopCounter += minminLoopCounter;
                totalMaxMinLoopCounter += maxminLoopCounter;
                totalNegMinLoopCounter += negminLoopCounter;
                totalCofMinMinLoopCounter += cofminminLoopCounter;
                totalCofMaxMinLoopCounter += cofmaxminLoopCounter;
                totalCountMaxMinLoopCounter += countmaxminLoopCounter;
                totalCountD2MinLoopCounter += countD2minLoopCounter;
                totalCountMax2MinLoopCounter += countmax2minLoopCounter;
                totalCountMinMinLoopCounter += countminminLoopCounter;
                totalEigenMinLoopCounter += eigenminLoopCounter;
                totalEigenMinMinLoopCounter += eigenminminLoopCounter;
                totalMaxNeighWMinLoopCounter += maxneighwminLoopCounter;

            }
            String res = String.format("Lp: %3d, V: %6d, mEC: %6d, mNEC: %6d  ",loop,V,maxEdgeCount,maxNodeEdgeCount);
            System.out.print(res);
            bw.write(loop + "," + V +","+maxEdgeCount+","+maxNodeEdgeCount+",");
            for (int j = 0; j < 26; j++) {
                if(j < 25)
                    bw.write(String.format("%.3f",result[loop][j])+",");
                else
                    bw.write(String.format("%.3f",result[loop][j]));
                String res2 = String.format("%9.3f", result[loop][j]);
                System.out.printf(res2);
            }

            bw.newLine();
            System.out.println();
        }

        double ga = ((double) gt/(testloop-avgloop));
        double da = ((double) dt/(testloop-avgloop));
        double minda = ((double) mindt/(testloop-avgloop));
        double maxda = ((double) maxdt/(testloop-avgloop));
        double negda = ((double) negdt/(testloop-avgloop));
        double cofminda = ((double) cofmindt/(testloop-avgloop));
        double cofmaxda = ((double) cofmaxdt/(testloop-avgloop));
        double countmaxda = ((double) countmaxdt/(testloop-avgloop));
        double countD2da = ((double) countD2dt/(testloop-avgloop));
        double countmax2da = ((double) countmax2dt/(testloop-avgloop));
        double countminda = ((double) countmindt/(testloop-avgloop));
        double eigenda = ((double) eigendt/(testloop-avgloop));
        double eigenMda = ((double) eigenMindt/(testloop-avgloop));
        double countmaxneighwda = ((double) countmaxneighwdt/(testloop-avgloop));
        double dapa = ((double) dapt/(testloop-avgloop));
        double bfa = ((double) bft/(testloop-avgloop));
        double fwa = ((double) fwt/(testloop-avgloop));
        double bcwa = ((double) bcwt/(testloop-avgloop));
        double edgeca = (double) TotalEdgeCount/(testloop-avgloop);
        double meedgcna = (double) meedgcnt/(testloop-avgloop);
        double lnNa = (double) lnNt/(testloop-avgloop);
        double lnNlnKa = (double) lnNlnKt/(testloop-avgloop);
        double nogedcna = (double) nogedcnt/(testloop-avgloop);
        double orpha = (double) orpht/(testloop-avgloop);
        double cofavga = ((double) cofavgt/(testloop-avgloop));
        double prmlpcnt  =  ((double) TotalPermLoopCounter/(testloop-avgloop));
        bw.write(",V,maxEdgeCount,maxNodeEdgeCount,GraphT,DAll,DMin,DMax,DNeg,DCofMin,DCofMax,DCountMax,DCountD2,DCountMax2,DCountMin,DEigen,DEigenM,CntMaxW,DLoop,BellmanT,FloydT,Best,EdgeC,MeanEdgeC,lnN,lnN/ln<k>,NoOutgoingC,OrphanedC,ClusCofAvg,PermAvg");
        bw.newLine();
        bw.write(",");bw.write(",");bw.write(",");bw.write(",");
        bw.write(String.format("%.3f",ga));
        bw.write(",");
        bw.write(String.format("%.3f",da));
        bw.write(",");
        bw.write(String.format("%.3f",minda));
        bw.write(",");
        bw.write(String.format("%.3f",maxda));
        bw.write(",");
        bw.write(String.format("%.3f",negda));
        bw.write(",");
        bw.write(String.format("%.3f",cofminda));
        bw.write(",");
        bw.write(String.format("%.3f",cofmaxda));
        bw.write(",");
        bw.write(String.format("%.3f",countmaxda));
        bw.write(",");
        bw.write(String.format("%.3f",countD2da));
        bw.write(",");
        bw.write(String.format("%.3f",countmax2da));
        bw.write(",");
        bw.write(String.format("%.3f",countminda));
        bw.write(",");
        bw.write(String.format("%.3f",eigenda));
        bw.write(",");
        bw.write(String.format("%.3f",eigenMda));
        bw.write(",");
        bw.write(String.format("%.3f",countmaxneighwda));
        bw.write(",");
        bw.write(String.format("%.3f",dapa));
        bw.write(",");
        bw.write(String.format("%.3f",bfa));
        bw.write(",");
        bw.write(String.format("%.3f",fwa));
        bw.write(",");
        bw.write(String.format("%.3f",bcwa));
        bw.write(",");
        bw.write(String.format("%.3f",edgeca));
        bw.write(",");
        bw.write(String.format("%.3f",meedgcna));
        bw.write(",");
        bw.write(String.format("%.3f",lnNa));
        bw.write(",");
        bw.write(String.format("%.3f",lnNlnKa));
        bw.write(",");
        bw.write(String.format("%.3f",nogedcna));
        bw.write(",");
        bw.write(String.format("%.3f",orpha));
        bw.write(",");
        bw.write(String.format("%.3f",cofavga));
        bw.write(",");
        bw.write(String.format("%.3f",prmlpcnt));
        bw.newLine();
        bw.newLine();

        bw.write(String.format(",,Dall:, %d, %5.2f \n" , TotaldfalseLoopCounter, (double) totaldminLoopCounter/(testloop-avgloop)));
        bw.write(String.format(",,DMin:, %d, %5.2f \n" , TotalminfalseLoopCounter, (double) totalMinMinLoopCounter/(testloop-avgloop)));
        bw.write(String.format(",,DMax:, %d, %5.2f \n" , TotalmaxfalseLoopCounter, (double) totalMaxMinLoopCounter/(testloop-avgloop)));
        bw.write(String.format(",,DNeg:, %d, %5.2f \n" , TotalnegfalseLoopCounter, (double) totalNegMinLoopCounter/(testloop-avgloop)));
        bw.write(String.format(",,DCofMin:, %d, %5.2f \n" , TotalcofMinfalseLoopCounter, (double) totalCofMinMinLoopCounter/(testloop-avgloop)));
        bw.write(String.format(",,DCofMax:, %d, %5.2f \n" , TotalcofMaxfalseLoopCounter, (double) totalCofMaxMinLoopCounter/(testloop-avgloop)));
        bw.write(String.format(",,DCountMax:, %d, %5.2f \n" , TotalcountMaxfalseLoopCounter, (double) totalCountMaxMinLoopCounter/(testloop-avgloop)));
        bw.write(String.format(",,DCountD2:, %d, %5.2f \n" , TotalcountD2falseLoopCounter, (double) totalCountD2MinLoopCounter/(testloop-avgloop)));
        bw.write(String.format(",,DCountMax2:, %d, %5.2f \n" , TotalcountMax2falseLoopCounter, (double) totalCountMax2MinLoopCounter/(testloop-avgloop)));
        bw.write(String.format(",,DCountMin:, %d, %5.2f \n" , TotalcountMinfalseLoopCounter, (double) totalCountMinMinLoopCounter/(testloop-avgloop)));
        bw.write(String.format(",,DEigen:, %d, %5.2f \n" , TotaleigenfalseLoopCounter, (double) totalEigenMinLoopCounter/(testloop-avgloop)));
        bw.write(String.format(",,DEigenMin:, %d, %5.2f \n" , TotaleigenMinfalseLoopCounter, (double) totalEigenMinMinLoopCounter/(testloop-avgloop)));
        bw.write(String.format(",,DCountMaxNeighW:, %d, %5.2f \n" , TotalMaxNeighWfalseLoopCounter, (double) totalMaxNeighWMinLoopCounter/(testloop-avgloop)));
        System.out.println("                                               |  GraphT|  DAll  |   DMin |   DMax |  DNeg  | CofMin | CofMax | CntMax |   CntD2| CntMax2| CntMin |  Eigen | EigenM | CntMaxW|    Loop| Bellmn |  Floyd |    Best|  EdgeC |  MEdgeC|      lnN|  lnN/lnk| NoOutC | OrphanC|ClCofAvg|PermAvg");
        System.out.print("                                               ");
        System.out.printf("%9.3f",ga); //GraphT:
        System.out.printf("%9.3f",da); //DAll:
        System.out.printf("%9.3f",minda); //DMin:
        System.out.printf("%9.3f",maxda); // DMax:
        System.out.printf("%9.3f",negda); //DNeg:
        System.out.printf("%9.3f",cofminda); //CofMin:
        System.out.printf("%9.3f",cofmaxda); //CofMax:
        System.out.printf("%9.3f",countmaxda); //CntMax:
        System.out.printf("%9.3f",countD2da); //CntD2:
        System.out.printf("%9.3f",countmax2da); //CntMax2:
        System.out.printf("%9.3f",countminda); //CntMin:
        System.out.printf("%9.3f",eigenda); //Eigen:
        System.out.printf("%9.3f",eigenMda); //EigenM:
        System.out.printf("%9.3f",countmaxneighwda); //countmaxneighwda:
        System.out.printf("%9.3f",dapa); //Loop:
        System.out.printf("%9.3f",bfa); //Bellmn:
        System.out.printf("%9.3f",fwa); //Floyd:
        System.out.printf("%9.3f",bcwa); //Best
        System.out.printf("%9.3f",edgeca); //Edge Count
        System.out.printf("%9.3f",meedgcna); //Mean Edge Count
        System.out.printf("%9.3f",lnNa); //lnN
        System.out.printf("%9.3f",lnNlnKa); //lnN/ln<k>
        System.out.printf("%9.3f",nogedcna); //NoOutgoingC
        System.out.printf("%9.3f",orpha); //OrphanedC
        System.out.printf("%9.3f",cofavga); //Avg Coefficient
        System.out.printf("%9.3f",prmlpcnt); //Perm Cnt
        System.out.println();
        System.out.printf("Dall: %d %5.2f \n" , TotaldfalseLoopCounter, (double) totaldminLoopCounter/(testloop-avgloop));
        System.out.printf("DMin: %d %5.2f \n" , TotalminfalseLoopCounter, (double) totalMinMinLoopCounter/(testloop-avgloop));
        System.out.printf("DMax: %d %5.2f \n" , TotalmaxfalseLoopCounter, (double) totalMaxMinLoopCounter/(testloop-avgloop));
        System.out.printf("DNeg: %d %5.2f \n" , TotalnegfalseLoopCounter, (double) totalNegMinLoopCounter/(testloop-avgloop));
        System.out.printf("DCofMin: %d %5.2f \n" , TotalcofMinfalseLoopCounter, (double) totalCofMinMinLoopCounter/(testloop-avgloop));
        System.out.printf("DCofMax: %d %5.2f \n" , TotalcofMaxfalseLoopCounter, (double) totalCofMaxMinLoopCounter/(testloop-avgloop));
        System.out.printf("DCountMax: %d %5.2f \n" , TotalcountMaxfalseLoopCounter, (double) totalCountMaxMinLoopCounter/(testloop-avgloop));
        System.out.printf("DCountD2: %d %5.2f \n" , TotalcountD2falseLoopCounter, (double) totalCountD2MinLoopCounter/(testloop-avgloop));
        System.out.printf("DCountMax2: %d %5.2f \n" , TotalcountMax2falseLoopCounter, (double) totalCountMax2MinLoopCounter/(testloop-avgloop));
        System.out.printf("DCountMin: %d %5.2f \n" , TotalcountMinfalseLoopCounter, (double) totalCountMinMinLoopCounter/(testloop-avgloop));
        System.out.printf("DEigen: %d %5.2f \n" , TotaleigenfalseLoopCounter, (double) totalEigenMinLoopCounter/(testloop-avgloop));
        System.out.printf("DEigenMin: %d %5.2f \n" , TotaleigenMinfalseLoopCounter, (double) totalEigenMinMinLoopCounter/(testloop-avgloop));
        System.out.printf("DCountMaxNeighW: %d %5.2f \n" , TotalMaxNeighWfalseLoopCounter, (double) totalMaxNeighWMinLoopCounter/(testloop-avgloop));

        bw.close();
//        for (int i = 0; i < V; i++) {
//            for (int j = 0; j < testloop; j++) {
//                if(degreeResult[i][j] != -1) {
//                    bwdegree.write(degreeResult[i][j] + ",");
//                    bwdegree.write(poissonResult[i][j] + ",");
//                }
//                else bwdegree.write(",,");
//            }
//            bwdegree.newLine();
//        }
//        bwdegree.close();
        Desktop.getDesktop().open(new File(fileName));
        //Desktop.getDesktop().open(new File(degreeFileName));
    }
}

