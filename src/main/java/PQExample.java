import java.util.PriorityQueue;

public class PQExample {
    public static void main(String[] args) {
        PriorityQueue<Integer> integerQueue = new PriorityQueue<>();
        PriorityQueue<Integer> integerQueueWithComparator = new PriorityQueue<>((Integer c1, Integer c2) -> Integer.compare(c1, c2));

        integerQueueWithComparator.add(5);
        integerQueue.add(5);

        integerQueueWithComparator.add(3);
        integerQueue.add(3);

        integerQueueWithComparator.add(2);
        integerQueue.add(2);

        integerQueueWithComparator.add(1);
        integerQueue.add(1);

        integerQueueWithComparator.add(4);
        integerQueue.add(4);

        System.out.println(integerQueue);
        System.out.println(integerQueueWithComparator);
    }
}
