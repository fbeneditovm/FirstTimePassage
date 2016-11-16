/******************************************************************************
 *  Robert Sedgewick and Kevin Wayne
 *  http://introcs.cs.princeton.edu/java/95linear/GaussianElimination.java.html 
 * 
 *  Gaussian elimination with partial pivoting.
 *
 *  % java GaussianElimination
 *  -1.0
 *  2.0
 *  2.0
 *
 ******************************************************************************/
package ufma.pe.firstpassagetime;

import java.math.BigDecimal;
import java.math.RoundingMode;

/**
 *
 * @author fbeneditovm
 */
public class GaussianElimitationBD {
    
    private static final BigDecimal EPSILON = new BigDecimal(1e-10);
    
    // Gaussian elimination with partial pivoting
    public static BigDecimal[] lsolve(BigDecimal[][] A, BigDecimal[] b) {
        int N  = b.length;

        for (int p = 0; p < N; p++) {

            // find pivot row and swap
            int max = p;
            for (int i = p + 1; i < N; i++) {
                if (A[i][p].abs().compareTo(A[max][p]) > 0) {
                    max = i;
                }
            }
            BigDecimal[] temp = A[p]; A[p] = A[max]; A[max] = temp;
            BigDecimal   t    = b[p]; b[p] = b[max]; b[max] = t;

            // singular or nearly singular
            if (A[p][p].abs().compareTo(EPSILON) <= 0) {
                throw new RuntimeException("Matrix is singular or nearly singular");
            }

            // pivot within A and b
            for (int i = p + 1; i < N; i++) {
                BigDecimal alpha = A[i][p].divide(A[p][p], 8, RoundingMode.HALF_UP);
                b[i] = b[i].subtract(alpha.multiply(b[p]));
                for (int j = p; j < N; j++) {
                    A[i][j] = A[i][j].subtract(alpha.multiply(A[p][j]));
                }
            }
        }

        // back substitution
        BigDecimal[] x = new BigDecimal[N];
        for (int i = N - 1; i >= 0; i--) {
            BigDecimal sum = new BigDecimal("0");
            for (int j = i + 1; j < N; j++) {
                sum = sum.add(A[i][j].multiply(x[j]));
            }
            x[i] = (b[i].subtract(sum)).divide(A[i][i], 8, RoundingMode.HALF_UP);
        }
        return x;
    }    
}
