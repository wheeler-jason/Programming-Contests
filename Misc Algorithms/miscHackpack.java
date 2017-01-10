public static long lcm(long a, long b) {
    return a * (b / gcd(a, b));
}
public static long lcm(long[] input) {
    long result = input[0];
    for(int i = 1; i < input.length; i++) 
        result = lcm(result, input[i]);
    return result;
}

//------------------------------------------------------------------------------

// Calculate the determinant of any size matrix
public static double det(double[][] a) {
    final double eps = 1e-9;
    int n = a.length;
    double res = 1;
    boolean[] used = new boolean[n];

    for (int i = 0; i < n; i++) {
      int p;
      for (p = 0; p < n; p++)
        if (!used[p] && Math.abs(a[p][i]) > eps)
          break;

      if (p >= n)
        return 0;

      res *= a[p][i];
      used[p] = true;

      double z = 1 / a[p][i];
      for (int j = 0; j < n; j++)
        a[p][j] *= z;
      for (int j = 0; j < n; ++j)
        if (j != p) {
          z = a[j][i];
          for (int k = 0; k < n; ++k)
            a[j][k] -= z * a[p][k];
        }
    }
    return res;
}

//------------------------------------------------------------------------------

public static int largestPrimeFactor(long number) {
    int i;

    for (i = 2; i <= number; i++) {
        if (number % i == 0) {
            number /= i;
            i--;
        }
    }

    return i;
}

public static void primeFactorize(int number) {
    int count;
    for (int i = 2; i<=(number); i++) {
        count = 0;
        while (number % i == 0) {
            number /= i;
            count++;
        }

        if (count != 0) {
            // Alternatively, append this to a StringBuilder or put
            // count into an int array indicating how many of that factor 
            System.out.println(i+ "**" + count);
        }
    }
}
// Finds num of prime factors for each number up to n
public static int[] numberOfPrimeDivisors(int n) {
    int[] divisors = new int[n + 1];
    Arrays.fill(divisors, 2, n + 1, 1);
    for (int i = 2; i * i <= n; ++i)
	    if (divisors[i] == 1)
		    for (int j = i; j * i <= n; j++)
			    divisors[j * i] = divisors[j] + 1;
    return divisors;
}

//------------------------------------------------------------------------------

public class MinCostFlow {

	public static int[] minCostFlow(int[][] cap, int[][] cost, int s, int t) {
		int n = cap.length;
		int[] d = new int[n];
		int[] p = new int[n];
		for (int flow = 0, flowCost = 0;; ++flow) {
			Arrays.fill(d, Integer.MAX_VALUE);
			d[s] = 0;
			for (int i = 0; i < n - 1; i++)
				for (int j = 0; j < n; j++)
					for (int k = 0; k < n; k++)
						if (cap[j][k] > 0 && d[j] < Integer.MAX_VALUE 
                                          && d[k] > d[j] + cost[j][k]) {
							d[k] = d[j] + cost[j][k];
							p[k] = j;
						}
			if (d[t] == Integer.MAX_VALUE)
				return new int[] { flowCost, flow };
			flowCost += d[t];
			for (int v = t; v != s; v = p[v]) {
				--cap[p[v]][v];
				++cap[v][p[v]];
			}
		}
	}

	public static void addEdge(int[][] cap, int[][] cost, int u, int v,
                                            int edgeCapacity, int edgeCost) {
		cap[u][v] = edgeCapacity;
		cost[u][v] = edgeCost;
		cost[v][u] = -edgeCost;
	}

	// Example test
	public static void main(String[] args) {
		int n = 3;
		int[][] capacity = new int[n][n];
		int[][] cost = new int[n][n];
		addEdge(capacity, cost, 0, 1, 3, 1);
		addEdge(capacity, cost, 0, 2, 2, 1);
		addEdge(capacity, cost, 1, 2, 2, 1);
		int[] costFlow = minCostFlow(capacity, cost, 0, 2);
		System.out.println(6 == costFlow[0]);
		System.out.println(4 == costFlow[1]);
	}
}

//------------------------------------------------------------------------------

// Longest increasing subsequence
public class Lis {

	public static int[] getLis(int[] x) {
		int n = x.length;
		int[] len = new int[n];
		Arrays.fill(len, 1);
		int[] pred = new int[n];
		Arrays.fill(pred, -1);
		int bi = 0;
		for (int i = 1; i < n; i++) {
			for (int j = 0; j < i; j++) {
				if (x[j] < x[i] && len[i] < len[j] + 1) {
					len[i] = len[j] + 1;
					pred[i] = j;
				}
			}
			if (len[bi] < len[i]) {
				bi = i;
			}
		}
		int cnt = len[bi];
		int[] res = new int[cnt];
		for (int i = bi; i != -1; i = pred[i]) {
			res[--cnt] = x[i];
		}
		return res;
	}

	// Example test
	public static void main(String[] args) {
		int[] a = { 1, 5, 4, 2, 3, 7, 6 };
		int[] lis = getLis(a);
		System.out.println(Arrays.toString(lis));
	}
}

//------------------------------------------------------------------------------



