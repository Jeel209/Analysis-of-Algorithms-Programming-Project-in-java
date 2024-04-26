import java.util.Arrays;
import java.util.Random;

class ClosestPairOfPoints {

    // Function to generate a random number 
    static int rand() {
        return new Random().nextInt(32768); 
    }


    // Function to generate the x and y coordinates for n points
    static double[][] generateRandomPoints(int n) {
        double[][] points = new double[n][2];
        Random randomGenerator = new Random();
        for (int i = 0; i < n; i++) {
            points[i][0] = randomGenerator.nextDouble();
            points[i][1] = randomGenerator.nextDouble();
        }
        return points;
    }

    // Function to compute the Euclidean distance between two points
    static double computeDistance(double x1, double y1, double x2, double y2) {
        return Math.sqrt(Math.pow((x1 - x2), 2) + Math.pow((y1 - y2), 2));
    }

    // Brute force algorithm to find the closest pair of points
    static void findClosestPairBruteForce(double[][] points, int n) {
        double minimumDistance = Double.POSITIVE_INFINITY;
        int point_Index1  = -1, point_Index2 = -1;

        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                double distance = computeDistance(points[i][0], points[i][1], points[j][0], points[j][1]);
                if (distance < minimumDistance) {
                    minimumDistance = distance;
                    point_Index1  = i;
                    point_Index2 = j;
                }
            }
        }
        System.out.printf("Brute Force: Closest pair is (%f, %f) and (%f, %f)\n",
                points[point_Index1 ][0], points[point_Index1 ][1], points[point_Index2][0], points[point_Index2][1]);
    }

    // Divide and conquer algorithm to find the closest pair of points
    static double[] findClosestPairDivideAndConquer(double[][] points, int left, int right) {
        if (left >= right)
            return new double[]{Double.POSITIVE_INFINITY, -1, -1};

        if (right - left == 1)
            return new double[]{computeDistance(points[left][0], points[left][1], points[right][0], points[right][1]), left, right};

        int mid = (left + right) / 2;
        double[] resultLeft = findClosestPairDivideAndConquer(points, left, mid);
        double[] resultRight =findClosestPairDivideAndConquer(points, mid + 1, right);

        double minimumDistance = Math.min(resultLeft[0], resultRight[0]);
        int point_Index1  = (resultLeft[0] < resultRight[0]) ? (int) resultLeft[1] : (int) resultRight[1];
        int point_Index2 = (resultLeft[0] < resultRight[0]) ? (int) resultLeft[2] : (int) resultRight[2];

        // Merge step
        double[][] strip = new double[right - left + 1][3];
        int k = 0;
        double midX = points[mid][0];
        for (int i = left; i <= right; i++) {
            if (Math.abs(points[i][0] - midX) < minimumDistance) {
                strip[k][0] = points[i][0];
                strip[k][1] = points[i][1];
                strip[k][2] = i;
                k++;
            }
        }

        Arrays.sort(strip, 0, k, (a, b) -> Double.compare(a[1], b[1]));

        for (int i = 0; i < k; i++) {
            for (int j = i + 1; j < Math.min(i + 8, k); j++) { 
                double dist = computeDistance(strip[i][0], strip[i][1], strip[j][0], strip[j][1]);
                if (dist < minimumDistance) {
                    minimumDistance = dist;
                    point_Index1  = (int) strip[i][2];
                    point_Index2 = (int) strip[j][2];
                }
            }
        }

        return new double[]{minimumDistance, point_Index1 , point_Index2};
    }

    // Main function to run the program
    public static void main(String[] args) {
        int[] ns = {10000,20000,30000,40000,50000,60000,70000,80000,90000,100000};
        int m = 10;

        for (int n : ns) {
            double totalBruteForceTime = 0.0;
            double totalDivideAndConquerTime = 0.0;

            for (int i = 0; i < m; i++) {
                double[][] points = generateRandomPoints(n);

                long startTime = System.nanoTime();
                findClosestPairBruteForce(points, n);
                long endTime = System.nanoTime();
                double bruteForceTime = (endTime - startTime) / 1e6; 
                totalBruteForceTime += bruteForceTime;

                startTime = System.nanoTime();
                double[] result = findClosestPairDivideAndConquer(points, 0, n - 1);
                endTime = System.nanoTime();
                double divideAndConquerTime = (endTime - startTime) / 1e6;
                totalDivideAndConquerTime += divideAndConquerTime;

                System.out.printf("n = %d, brute force time = %.2f ms, divide and conquer time = %.2f ms\n", n, bruteForceTime, divideAndConquerTime);
            }

            double averageBruteForceTime = totalBruteForceTime / m;
            double averageDivideAndConquerTime = totalDivideAndConquerTime / m;
            System.out.printf("Average running time for brute force algorithm with n = %d: %.2f ms\n", n, averageBruteForceTime);
            System.out.printf("Average running time for divide and conquer algorithm with n = %d: %.2f ms\n", n, averageDivideAndConquerTime);
}
}
}