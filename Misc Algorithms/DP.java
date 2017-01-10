// Longest Increasing Subsequence O(nlogn)
public static int LIS(int[] x) {
    int[] arr = new int[x.length];
    Arrays.fill(arr, Integer.MAX_VALUE);
    arr[0] = x[0];
    for(int i = 1; i < x.length; i++)
        arr[binSearch(arr, x[i])] = x[i];

    return binSearch(arr, Integer.MAX_VALUE / 2);
}

public static int binSearch(int[] x, int val) {
    int start = 0, end = x.length - 1, ret = -1;
    while(start <= end) {
        int m = (start+end)/2;
        if(x[m] > val)
            end = (ret = m) - 1;
        else
            start = m + 1;
    }

    return ret;
}

// Longest Common Subsequence utility
public static int[][] getLCSUtil(String a, String b) {
    int[][] table = new int[a.length() + 1][b.length() + 1];
    
    for(int i = 1; i <= a.length(); i++)
        for(int j = 1; j <= b.length(); j++)
            if(a.charAt(i - 1) == b.charAt(j - 1))
                table[i][j] = table[i - 1][j - 1] + 1;
            else if(table[i - 1][j] > table[i][j - 1])
                    table[i][j] = table[i - 1][j];
            else
                    table[i][j] = table[i][j - 1];
    
    return table;
}
// Longest Common Subsequence
// pre-condition: run the util first and pass this function the returned table
public static String getLCS(int[][] table, String a, int m, int n) {
    if(table[m][n] == 0)
        return "";
    
    if(table[m][n] == table[m - 1][n])
        return getSequence(table, a, m - 1, n);
    else if(table[m][n] == table[m][n - 1])
        return getSequence(table, a, m, n - 1);
    else
        return getSequence(table, a, m - 1, n - 1) + a.charAt(m - 1);
}

public static int knapsack(int[] w, int[] v, int cap) {
    int sack[] = new int[cap + 1];
    for(int i = 0; i < w.length; i++)
        for(int j = cap; j >= w[i]; j--) // 0-1 knapsack
            //for(int j = w[i]; j <= cap; j++) // unbounded knapsack
            sack[j] = Math.max(sack[j], sack[j - w[i]] + v[i]);
    return sack[cap];
}

public static boolean subsetSum(int[] set, int k, int sum) {
    
    if(k == set.length)
        return true;
    
    return subsetSum(set, k + 1, sum - set[k]) || subsetSum(set, k + 1, sum);
}