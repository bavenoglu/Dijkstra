import java.util.Arrays;

public class Dijkstra extends Thread{
    Graph g;
    int src;
    Dijkstra(Graph g, int src){
        this.g = g;
        this.src = src;
    }
    public void run() {
        g.pq.clear();
        Arrays.fill(g.dijkstraDistance, 999);
        Arrays.fill(g.dijkstraPath, 999);
        Arrays.fill(g.dijkstraPathHop, 0);
        g.pq.add(new Graph.iPair(0, src));
        g.dijkstraDistance[src] = 0;
        g.dijkstraPath[src] = 0;
        while (!g.pq.isEmpty()) {
            int u = g.pq.poll().second;
            for (Graph.iPair v : g.adj.get(u)) {
                int distSum = g.dijkstraDistance[u] + v.second;
                int dD = g.dijkstraDistance[v.first];
                if (dD > distSum) {
                    g.dijkstraDistance[v.first] = distSum;
                    int temp = g.dijkstraArray[v.first][src];
                    g.dijkstraArray[v.first][src] = distSum;
                    g.dijkstraArray[src][v.first] = distSum;
                    if (temp == -1) {
                        g.dijkstraArrayCount--;g.dijkstraArrayCount--;
                    }
                    g.pq.add(new Graph.iPair(dD, v.first));
                    g.dijkstraPath[v.first] = u;
                }
            }
        }

        for (int j = 0; j < g.V; j++) {
            if (g.dijkstraPath[j] == 999){
                for (int i = 0; i < g.V; i++)
                    if (g.dijkstraPath[i] != 999 && g.dijkstraArray[j][i] != 999){
                        g.dijkstraArray[j][i] = 999;
                        g.dijkstraArray[i][j] = 999;
                        g.dijkstraArrayCount--;g.dijkstraArrayCount--;
                    }
                g.dijkstraPathHop[j] = 0;
            }
            else {
                if(j != src)
                {
                    int k = j;
                    int i = g.dijkstraPath[j];
                    g.dijkstraPathHop[j] = 0;
                    if(i != src)
                    {
                        do{
                            do{
                                if(g.dijkstraArray[k][i] == -1)
                                {
                                    int weightSource = g.dijkstraDistance[k];
                                    int weightDest = g.dijkstraDistance[i];
                                    g.dijkstraArray[k][i] = Math.abs(weightDest - weightSource);
                                    g.dijkstraArray[i][k] = Math.abs(weightDest - weightSource);
                                    g.dijkstraArrayCount--;g.dijkstraArrayCount--;
                                }
                                i = g.dijkstraPath[i];
                            }while(i != src);
                            k = g.dijkstraPath[k];
                            i = g.dijkstraPath[k];
                            g.dijkstraPathHop[j] = ++g.dijkstraPathHop[j];
                        } while(k != src && i != src );
                    }
                }
            }
        }
    }
}
