import java.util.Scanner;

public class Interior_point{


     public static void main(String[] args) {
          // input
          Scanner in = new Scanner(System.in);

          System.out.println("The number of variables:");
          int n = in.nextInt();

          System.out.println("Vector C / The gradient / The coefficient of objective function:");
          double[][] c_arr = new double[n][1];
          for (int i = 0; i < n; i++) {
               c_arr[i][0] = in.nextDouble();
          }
          Matrix C;
          C = new Matrix(c_arr);

          System.out.println("The number of constraints:");
          int m = in.nextInt();

          System.out.println("The matrix A / The coefficient of constraint function:");
          double[][] A_arr = new double[m][n];
          for (int i = 0; i < m; i++) {
               for (int j = 0; j < n; j++) {
                    A_arr[i][j] = in.nextDouble();
               }
          }
          Matrix A = new Matrix(A_arr);

          System.out.println("Vector x / Initial solution:");
          double[][] x_arr = new double[n][1];
          boolean isnonnegative = true;
          for (int i=0; i<n; i++){
               x_arr[i][0] = in.nextDouble();
               if (x_arr[i][0] < 0){ isnonnegative = false;}
          }
          Matrix x = new Matrix(x_arr);
          if(!isnonnegative){
               System.out.println("Wrong input! Variables must be non-negative");
               return;
          }

          System.out.println("Vector b / The right-hand side numbers:");
          double[][] b_arr = new double[m][1];
          for (int i = 0; i < m; i++) {
               b_arr[i][0] = in.nextDouble();
          }
          Matrix b = new Matrix(b_arr);

          //check if A*x = b:
          if (!A.multiply(x).equals(b)){
               System.out.println("Wrong input! Ax must be equal b");
               return;
          }

          System.out.println("The approximation accuracy:");
          double accuracy = in.nextDouble();
          if (accuracy < 0 ) {
               System.out.println("The approximation accuracy is negative");
               return;
          }

          double alpha = 0.5;

          int counter=0;

          while (true) {

               counter++;
               if (counter>5000){
                    System.out.println("The method is not applicable!");
                    break;
               }

               double[][] D_arr = new double[n][n];
               for (int i=0; i<n; i++) for (int j=0; j<n; j++) {
                    if (i==j) D_arr[i][j] = x.data[i][0]; else D_arr[i][j]=0;
               }

               double[][] I_arr = new double[n][n];
               for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                         if (i != j) I_arr[i][j] = 0;
                         else I_arr[i][j] = 1;
               Matrix I = new Matrix(I_arr);

               Matrix D = new Matrix(D_arr);

               Matrix A_tilda = A.multiply(D);

               Matrix C_tilda = D.multiply(C);

               Matrix Z = A_tilda.multiply(A_tilda.transpose()).invert(); // Z = (A~ x A~T)^-1

               Matrix P = I.subtract(A_tilda.transpose().multiply(Z.multiply(A_tilda)));

               Matrix C_p = P.multiply(C_tilda);

               double minimal_value = Integer.MAX_VALUE;
               for (int i = 0; i < n; i++) {
                    if (C_p.data[i][0] < minimal_value) {
                         minimal_value = C_p.data[i][0];
                    }
               }

               double v = Math.abs(minimal_value);

               double[][] ones_arr = new double[n][1];
               for (int i=0; i<n; i++) ones_arr[i][0]=1;
               Matrix ones = new Matrix(ones_arr);

               Matrix C_p_v2 = C_p;

               for (int i = 0; i < n; i++) {
                    C_p_v2.data[i][0] = C_p.data[i][0] * alpha / v;
               }

               Matrix X_tilda = ones.add(C_p_v2);

               Matrix X_new = D.multiply(X_tilda);

               Matrix difference = X_new.subtract(x);

               double difference_length = 0;

               for (int i=0; i<n; i++) {
                    difference_length+=Math.pow(difference.data[i][0],2);
               }
               difference_length=Math.pow(difference_length, 0.5);

               if (difference_length<accuracy) {
                    System.out.println("Final vector X:");
                    for (int i=0; i<n; i++){
                         System.out.print(X_new.data[i][0]+" ");
                         System.out.println();
                    }
                    System.out.println("Result:");
                    double z = 0;
                    for (int i=0; i<n; i++){
                         z+=X_new.data[i][0]*C.data[i][0];
                    }
                    System.out.println(z);
                    break;
               }

               x = X_new;

          }
     }
}


class Matrix {
     public double[][] data;
     public int rows;
     public int cols;

     public boolean equals(Matrix other){
          if (this.rows == other.rows && this.cols == other.cols){
               for (int i = 0; i < rows; i ++){
                    for (int j = 0; j < cols; j++) {
                         if (this.data[i][j] != other.data[i][j]){
                              return false;
                         }
                    }
               }
               return true;
          } else{
               return false;
          }
     }

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
     public void checkDimensionCompatibility(Matrix other) {
          if (this.rows != other.rows || this.cols != other.cols) {
               throw new IllegalArgumentException("Matrices must have the same dimensions.");
          }
     }



}

