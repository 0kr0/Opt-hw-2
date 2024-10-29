import java.util.Scanner;
import java.util.Arrays;

public class Interior_point{


     public static void main(String[] args) {
          // input
          Scanner in = new Scanner(System.in);

          System.out.println("The size of vector X/The number of variables:");
          int n = in.nextInt();

          System.out.println("The vector C/The coefficient of objective function:");
          double[][] c_arr = new double[n][1];
          for (int i = 0; i < n; i++) {
               c_arr[i][0] = in.nextDouble();
          }
          Matrix C;
          C = new Matrix(c_arr);

          System.out.println("The size of matrix A/The number of constraints:");
          int m = in.nextInt();

          System.out.println("The matrix A/The coefficient of constraint function:");
          double[][] A_arr = new double[m][n];
          for (int i = 0; i < m; i++) {
               for (int j = 0; j < n; j++) {
                    A_arr[i][j] = in.nextDouble();
               }
          }
          Matrix A = new Matrix(A_arr);

          System.out.println("Initial solution:");
          double[][] D_arr = new double[n][n];
          for (int i=0; i<n; i++) for (int j=0; j<n; j++) {
               if (i==j) D_arr[i][j] = in.nextDouble(); else D_arr[i][j]=0;
          }


          System.out.println("The vector b/The right-hand side numbers:");
          double[][] b_arr = new double[m][1];
          for (int i = 0; i < m; i++) {
               b_arr[i][0] = in.nextDouble();
          }
          Matrix b = new Matrix(b_arr);

          System.out.println("The approximation accuracy:");
          int accuracy = in.nextInt();
          if (accuracy < 0 ) {
               System.out.println("The approximation accuracy is negative");
               return;
          }

          // //check
          // if (!isApplicable(b)){
          //     System.out.println("The method is not applicable!");
          //     return;
          // }

          double alpha = 0.5;

          Matrix D = new Matrix(D_arr);
     }

}

class Matrix {
     private double[][] data;
     private final int rows;
     private final int cols;

     public Matrix(double[][] data) {
          this.rows = data.length;
          this.cols = data[0].length;
          this.data = new double[rows][cols];
          for (int i = 0; i < rows; i++) {
               System.arraycopy(data[i], 0, this.data[i], 0, cols);
          }
     }

     // Matrix addition
     public Matrix add(Matrix other) {
          checkDimensionCompatibility(other);
          double[][] result = new double[rows][cols];
          for (int i = 0; i < rows; i++) {
               for (int j = 0; j < cols; j++) {
                    result[i][j] = this.data[i][j] + other.data[i][j];
               }
          }
          return new Matrix(result);
     }

     // Matrix subtraction
     public Matrix subtract(Matrix other) {
          checkDimensionCompatibility(other);
          double[][] result = new double[rows][cols];
          for (int i = 0; i < rows; i++) {
               for (int j = 0; j < cols; j++) {
                    result[i][j] = this.data[i][j] - other.data[i][j];
               }
          }
          return new Matrix(result);
     }

     // Matrix multiplication
     public Matrix multiply(Matrix other) {
          if (this.cols != other.rows) {
               throw new IllegalArgumentException("Matrices cannot be multiplied: incompatible dimensions.");
          }
          double[][] result = new double[this.rows][other.cols];
          for (int i = 0; i < this.rows; i++) {
               for (int j = 0; j < other.cols; j++) {
                    for (int k = 0; k < this.cols; k++) {
                         result[i][j] += this.data[i][k] * other.data[k][j];
                    }
               }
          }
          return new Matrix(result);
     }

     // Matrix transposition
     public Matrix transpose() {
          double[][] result = new double[cols][rows];
          for (int i = 0; i < rows; i++) {
               for (int j = 0; j < cols; j++) {
                    result[j][i] = this.data[i][j];
               }
          }
          return new Matrix(result);
     }

     // Matrix inversion using Gaussian elimination
     public Matrix invert() {
          if (rows != cols) {
               throw new IllegalArgumentException("Only square matrices can be inverted.");
          }
          int n = rows;
          double[][] result = new double[n][n];
          double[][] augmented = new double[n][2 * n];

          // Prepare the augmented matrix
          for (int i = 0; i < n; i++) {
               System.arraycopy(this.data[i], 0, augmented[i], 0, n);
               augmented[i][i + n] = 1.0;
          }


          // Gaussian elimination
          for (int i = 0; i < n; i++) {
               double max = augmented[i][i];
               int maxRow = i;
               for (int k = i + 1; k < n; k++) {
                    if (Math.abs(augmented[k][i]) > max) {
                         max = Math.abs(augmented[k][i]);
                         maxRow = k;
                    }
               }

               double[] temp = augmented[i];
               augmented[i] = augmented[maxRow];
               augmented[maxRow] = temp;

               double divisor = augmented[i][i];
               for (int j = 0; j < 2 * n; j++) {
                    augmented[i][j] /= divisor;
               }

               for (int k = 0; k < n; k++) {
                    if (k != i) {
                         double factor = augmented[k][i];
                         for (int j = 0; j < 2 * n; j++) {
                              augmented[k][j] -= factor * augmented[i][j];
                         }
                    }
               }
          }

          // Extract the inverse matrix
          for (int i = 0; i < n; i++) {
               System.arraycopy(augmented[i], n, result[i], 0, n);
          }

          return new Matrix(result);
     }

     // Helper function to check matrix dimensions for addition and subtraction
     private void checkDimensionCompatibility(Matrix other) {
          if (this.rows != other.rows || this.cols != other.cols) {
               throw new IllegalArgumentException("Matrices must have the same dimensions.");
          }
     }

     // Print the matrix
     public void print() {
          for (double[] row : data) {
               System.out.println(Arrays.toString(row));
          }
     }
}
