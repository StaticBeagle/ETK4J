package com.wildbitsfoundry.etk4j.math.linearalgebra;

import java.util.Arrays;

import static com.wildbitsfoundry.etk4j.math.linearalgebra.ColumnCounts.adjust;

class QrStructuralCounts {

    MatrixSparse A; // reference to input matrix
    int m, n; // short hand for number of rows and columns in A
    int[] leftmost = new int[0]; // left most column in each row
    int m2; // number of rows for QR after adding fictitious rows
    int[] pinv = new int[0]; // inverse permutation to ensure diagonal elements are all structurally nonzero
    // this is a row pivot
    int[] parent = new int[0]; // elimination tree
    int[] post = new int[0]; // post ordered tree
    IGrowArray gwork = new IGrowArray(); // generic work space
    int nz_in_V; // number of entries in V
    int nz_in_R; // number of entries in R
    int[] countsR = new int[0]; // column counts in R

    // ----- start location of different sections in work array inside of V
    int next;         // col=next[row] element in the linked list
    int head;         // row=head[col] first row in which col is the first
    int tail;         // row=tail[col] last row in which col is the first
    int nque;         // nque[col] number of elements in column linked list

    ColumnCounts columnCounts = new ColumnCounts(true);

    /**
     * Examins the structure of A for QR decomposition
     *
     * @param A matrix which is to be decomposed
     * @return true if the solution is valid or false if the decomposition can't be performed (i.e. requires column pivots)
     */
    public boolean process( MatrixSparse A ) {
        init(A);

        eliminationTree(A, true, parent, gwork);

        countNonZeroInR(parent);
        countNonZeroInV(parent);

        // if more columns than rows it's possible that Q*R != A. That's because a householder
        // would need to be created that's outside the  m by m Q matrix. In reality it has
        // a partial solution. Column pivot are needed.
        if (m < n) {
            for (int row = 0; row < m; row++) {
                if (gwork.data[head + row] < 0) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * Initializes data structures
     */
    void init( MatrixSparse A ) {
        this.A = A;
        this.m = A.rows;
        this.n = A.cols;

        this.next = 0;
        this.head = m;
        this.tail = m + n;
        this.nque = m + 2*n;

        if (parent.length < n || leftmost.length < m) {
            parent = new int[n];
            post = new int[n];
            pinv = new int[m + n];
            countsR = new int[n];
            leftmost = new int[m];
        }
        gwork.reshape(m + 3*n);
    }

    /**
     * Count the number of non-zero elements in R
     */
    void countNonZeroInR( int[] parent ) {
        postorder(parent, n, post, gwork);
        columnCounts.process(A, parent, post, countsR);
        nz_in_R = 0;
        for (int k = 0; k < n; k++) {
            nz_in_R += countsR[k];
        }
        if (nz_in_R < 0)
            throw new RuntimeException("Too many elements. Numerical overflow in R counts");
    }

    /**
     * <p>If ata=false then it computes the elimination tree for sparse lower triangular square matrix
     * generated from Cholesky decomposition. If ata=true then it computes the elimination tree of
     * A<sup>T</sup>A without forming A<sup>T</sup>A explicitly. In an elimination tree the parent of
     * node 'i' is 'j', where the first off-diagonal non-zero in column 'i' has row index 'j'; j &gt; i
     * for which l[k,i] != 0.</p>
     *
     * <p>This tree encodes the non-zero elements in L given A, e.g. L*L' = A, and enables faster to compute solvers
     * than the general purpose implementations.</p>
     *
     * <p>Functionally identical to cs_etree in csparse</p>
     *
     * @param A (Input) M by N sparse upper triangular matrix. If ata is false then M=N otherwise M &ge; N
     * @param ata If true then it computes elimination treee of A'A without forming A'A otherwise computes elimination
     * tree for cholesky factorization
     * @param parent (Output) Parent of each node in tree. This is the elimination tree. -1 if no parent. Size N.
     * @param gwork (Optional) Internal workspace. Can be null.
     */
    public static void eliminationTree( MatrixSparse A, boolean ata, int[] parent, IGrowArray gwork ) {
        final int m = A.rows;
        final int n = A.cols;

        if (parent.length < n)
            throw new IllegalArgumentException("parent must be of length N");

        final int[] work = adjust(gwork, n + (ata ? m : 0));

        final int ancestor = 0; // reference to index in work array
        final int previous = n; // reference to index in work array

        if (ata) {
            for (int i = 0; i < m; i++) {
                work[previous + i] = -1;
            }
        }

        // step through each column
        for (int k = 0; k < n; k++) {
            parent[k] = -1;
            work[ancestor + k] = -1;

            int idx0 = A.col_idx[k];   // node k has no parent
            int idx1 = A.col_idx[k + 1]; // node k has no ancestor

            for (int p = idx0; p < idx1; p++) {

                int nz_row_p = A.nz_rows[p];

                int i = ata ? work[previous + nz_row_p] : nz_row_p;

                int inext;
                while (i != -1 && i < k) {
                    inext = work[ancestor + i];
                    work[ancestor + i] = k;
                    if (inext == -1) {
                        parent[i] = k;
                        break;
                    } else {
                        i = inext;
                    }
                }

                if (ata) {
                    work[previous + nz_row_p] = k;
                }
            }
        }
    }


    /**
     * <p>Sorts an elimination tree {@link #eliminationTree} into postorder. In a postoredered tree, the d proper
     * descendants of any node k are numbered k-d through k-1. Non-recursive implementation for better performance.</p>
     *
     * <p>post[k] = i means node 'i' of the original tree is node 'k' in the postordered tree.</p>
     *
     * <p>See page 44</p>
     *
     * @param parent (Input) The elimination tree.
     * @param N Number of elements in parent
     * @param post (Output) Postordering permutation.
     * @param gwork (Optional) Internal workspace. Can be null
     */
    public static void postorder( int[] parent, int N, int[] post, IGrowArray gwork ) {
        if (parent.length < N)
            throw new IllegalArgumentException("parent must be at least of length N");
        if (post.length < N)
            throw new IllegalArgumentException("post must be at least of length N");

        int[] w = adjust(gwork, 3*N);

        // w[0] to w[N-1] is initialized to the youngest child of node 'j'
        // w[N] to w[2N-1] is initialized to the second youngest child of node 'j'
        // w[2N] to w[3N-1] is the stacked of nodes to be examined in the dfs
        final int next = N;

        // specify the linked list as being empty initially
        for (int j = 0; j < N; j++) {
            w[j] = -1;
        }
        // traverse nodes in reverse order
        for (int j = N - 1; j >= 0; j--) {
            // skip if j has no parent, i.e. is a root node
            if (parent[j] == -1)
                continue;
            // add j to the list of parents
            w[next + j] = w[parent[j]];
            w[parent[j]] = j;
        }

        // perform the DFS on each root node
        int k = 0;
        for (int j = 0; j < N; j++) {
            if (parent[j] != -1)
                continue;

            k = postorder_dfs(j, k, w, post, N);
        }
    }

    /**
     * Depth First Search used inside of {@link #postorder}.
     */
    protected static int postorder_dfs( int j, int k, int[] w, int[] post, int N ) {
        final int next = N;
        final int stack = 2*N;
        int top = 0; // top of the stack
        w[stack + top] = j;
        while (top >= 0) {
            int p = w[stack + top]; // next index in the stack to process
            int i = w[p];         // yongest child of p

            if (i == -1) {
                // p has no more unordered children left to process
                top--;
                post[k++] = p;
            } else {
                w[p] = w[next + i];
                top++;
                w[stack + top] = i;
            }
        }
        return k;
    }

    /**
     * Count the number of non-zero elements in V
     */
    void countNonZeroInV( int[] parent ) {
        int[] w = gwork.data;
        findMinElementIndexInRows(leftmost);
        createRowElementLinkedLists(leftmost, w);
        countNonZeroUsingLinkedList(parent, w);
    }

    /**
     * Non-zero counts of Householder vectors and computes a permutation
     * matrix that ensures diagonal entires are all structurally nonzero.
     *
     * @param parent elimination tree
     * @param ll linked list for each row that specifies elements that are not zero
     */
    void countNonZeroUsingLinkedList( int[] parent, int[] ll ) {

        Arrays.fill(pinv, 0, m, -1);
        nz_in_V = 0;
        m2 = m;

        for (int k = 0; k < n; k++) {
            int i = ll[head + k];           // remove row i from queue k
            nz_in_V++;                    // count V(k,k) as nonzero
            if (i < 0)                    // add a fictitious row since there are no nz elements
                i = m2++;
            pinv[i] = k;                  // associate row i with V(:,k)
            if (--ll[nque + k] <= 0)
                continue;
            nz_in_V += ll[nque + k];
            int pa;
            if ((pa = parent[k]) != -1) { // move all rows to parent of k
                if (ll[nque + pa] == 0)
                    ll[tail + pa] = ll[tail + k];
                ll[next + ll[tail + k]] = ll[head + pa];
                ll[head + pa] = ll[next + i];
                ll[nque + pa] += ll[nque + k];
            }
        }
        for (int i = 0, k = n; i < m; i++) {
            if (pinv[i] < 0)
                pinv[i] = k++;
        }

        if (nz_in_V < 0)
            throw new RuntimeException("Too many elements. Numerical overflow in V counts");
    }

    /**
     * Constructs a linked list in w that specifies which elements in each row are not zero (nz)
     *
     * @param leftmost index first elements in each row
     * @param w work space array
     */
    void createRowElementLinkedLists( int[] leftmost, int[] w ) {
        for (int k = 0; k < n; k++) {
            w[head + k] = -1;
            w[tail + k] = -1;
            w[nque + k] = 0;
        }

        // scan rows in reverse order creating a linked list of nz element indexes in each row
        for (int i = m - 1; i >= 0; i--) {
            int k = leftmost[i];      // 'k' = left most column in row 'i'
            if (k == -1)             // row 'i' is empty
                continue;
            if (w[nque + k]++ == 0)
                w[tail + k] = i;
            w[next + i] = w[head + k];
            w[head + k] = i;
        }
    }

//    private void printLinkedList(int w[]) {
//        System.out.print("head [");
//        for (int k = 0; k < n; k++) {
//            System.out.printf(" %2d",w[head+k]);
//        }
//        System.out.println(" ]");
//        System.out.print("tail [");
//        for (int k = 0; k < n; k++) {
//            System.out.printf(" %2d",w[tail+k]);
//        }
//        System.out.println(" ]");
//        System.out.print("nque [");
//        for (int k = 0; k < n; k++) {
//            System.out.printf(" %2d",w[nque+k]);
//        }
//        System.out.println(" ]");
//        System.out.print("next [");
//        for (int k = 0; k < m; k++) {
//            System.out.printf(" %2d",w[next+k]);
//        }
//        System.out.println(" ]");
//    }

    /**
     * Computes leftmost[i] =  min(find(A[i,:))
     * *
     *
     * @param leftmost (output) storage for left most elements
     */
    void findMinElementIndexInRows( int[] leftmost ) {
        Arrays.fill(leftmost, 0, m, -1);

        // leftmost[i] = min(find(A(i,:)))
        for (int k = n - 1; k >= 0; k--) {
            int idx0 = A.col_idx[k];
            int idx1 = A.col_idx[k + 1];

            for (int p = idx0; p < idx1; p++) {
                leftmost[A.nz_rows[p]] = k;
            }
        }
    }

    public void setGwork( IGrowArray gwork ) {
        this.gwork = gwork;
    }

    public int getFicticousRowCount() {
        return m2;
    }

    public int[] getLeftMost() {
        return leftmost;
    }

    public int[] getParent() {
        return parent;
    }

    public int[] getPinv() {
        return pinv;
    }

    public int getM2() {
        return m2;
    }
}
