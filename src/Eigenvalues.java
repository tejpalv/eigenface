
import Jama.Matrix;
import Jama.EigenvalueDecomposition;

public class Eigenvalues {
   public static void main(String[] args) { 
      int N = 5;

      // create a symmetric positive definite matrix
      Matrix A = Matrix.random(N, N);
      A = A.transpose().times(A);

      // compute the spectral decomposition
      EigenvalueDecomposition e = A.eig();
      Matrix V = e.getV();
      Matrix D = e.getD();

      StdOut.print("A =");
      A.print(9, 6);
      StdOut.print("D =");
      D.print(9, 6);
      StdOut.print("V =");
      V.print(9, 6);

      // check that V is orthogonal
      StdOut.print("||V * V^T - I|| = ");
      StdOut.println(V.times(V.transpose()).minus(Matrix.identity(N, N)).normInf());

      // check that A V = D V
      StdOut.print("||AV - DV|| = ");
      StdOut.println(A.times(V).minus(V.times(D)).normInf());
   }

}
