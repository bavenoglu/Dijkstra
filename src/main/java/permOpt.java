import java.util.Arrays;
import java.util.Vector;

public class permOpt {
    Graph g;
    Vector<Integer> tmp = new Vector<Integer>();
    Vector<Integer> realList = new Vector<Integer>();
    int[] correctPerm;
    boolean letExit = false;
    int dijkstraArrayPrevious[][];
    int dijkstraArrayCountPrevious = 0;
    permOpt(Graph g){
        this.g = g;
        correctPerm = new int[g.V];
    }
    void makeCombiUtil(int n, int left, int k) {
        if (k == 0) {
            //for(int i = 0; i < tmp.size(); i++)
            //{
            if(tmp.size() > 0){
                this.g.Dijkstra3(tmp.getLast());
                //System.out.print(tmp);System.out.println(tmp.getLast());
                if(g.dijkstraArrayCount == 0){
                    correctPerm[tmp.size()]++;
                    letExit = true;
                    realList = new Vector<>(tmp);
                    System.out.print(tmp); //to print the list of permutations open this
                    return;
                }
            }
            return;
        }
        else {
            for(int i = 0; i < tmp.size(); i++){
                    this.g.Dijkstra3(tmp.get(i));
                    //System.out.print(tmp);System.out.println(tmp.get(i));
                    if(g.dijkstraArrayCount == 0){
                        correctPerm[tmp.size()]++;
                        letExit = true;
                        realList = new Vector<>(tmp);
                        break;
                    }
                }
        }
        if(!letExit){
            dijkstraArrayPrevious = Arrays.stream(g.dijkstraArray)
                                    .map(int[]::clone).parallel()
                                    .toArray(int[][]::new);
            dijkstraArrayCountPrevious = g.dijkstraArrayCount;
            for (int i = left; i <= n; ++i) {
                tmp.add(i);
                if(!letExit)
                    makeCombiUtil(n, i + 1, k - 1);
                // Popping out last inserted element
                // from the vector
                tmp.remove(tmp.size() - 1);
                g.dijkstraArray = Arrays.stream(dijkstraArrayPrevious)
                                .map(int[]::clone).parallel()
                                .toArray(int[][]::new);
                g.dijkstraArrayCount = dijkstraArrayCountPrevious;
            }
        }
        g.resetDijkstra();
        dijkstraArrayPrevious = Arrays.stream(g.dijkstraArray)
                .map(int[]::clone).parallel()
                .toArray(int[][]::new);
        dijkstraArrayCountPrevious = g.dijkstraArrayCount;
    }
    void makeCombi(int n, int k) {
        makeCombiUtil(n-1, 0, k);
    }
    public void generatePerm()
    {
        g.resetDijkstra();
        for(int i = 0; i < g.V; i++){
            if(letExit)
                break;
            makeCombi(g.V, i);
        }
    }
}
