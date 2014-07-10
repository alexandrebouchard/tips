package tips.pip;




import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;


import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import bayonet.distributions.Multinomial;
import briefj.Indexer;




public class PoissonParameters 
{
  
  public final double[][] subRateMtx;
  public final double insertRate,deleteRate;
  public final int numberOfCharacter, numberOfCharacterPlusGap;
  public final int gapIndex;
  public final double [][] Q;
  public final Indexer<Character> indexer;
  public final double [] quasiStatLogProbabilities, quasiStatProbs;

  public PoissonParameters(
      Indexer<Character> indexer, 
      double[][] subRateMtx,
      double insertRate, 
      double deletionRate)
  {
    if (insertRate <= 0.0 || deletionRate <= 0.0)
      throw new RuntimeException();
    if (indexer.size() != subRateMtx.length)
      throw new RuntimeException();
    if (subRateMtx.length != subRateMtx[0].length)
      throw new RuntimeException();
    
    this.indexer = indexer;
    this.insertRate = insertRate;
    this.deleteRate = deletionRate;
    this.subRateMtx = subRateMtx;
    this.numberOfCharacter = subRateMtx.length;
    this.numberOfCharacterPlusGap = numberOfCharacter + 1;
    this.gapIndex = numberOfCharacter;
    final double [] deleteRates = uniformDeletionRate(deletionRate, numberOfCharacter);
    this.Q = formQMtx(subRateMtx, deleteRates);
    this.quasiStatLogProbabilities = quasiStationaryDistributionFromRates(subRateMtx, deleteRates);
    this.quasiStatProbs = quasiStatLogProbabilities.clone();
    logInPlace(this.quasiStatLogProbabilities);
  }
  
  public static void logInPlace(double [] n)
  {
    for (int i = 0 ; i < n.length; i++)
      n[i] = Math.log(n[i]);
  }
  
  public static double [][] formQMtx(
      double [][] reversibleConditionalSubstitutionRates, 
      double [] absorptionsRates)
  {
    final int N = reversibleConditionalSubstitutionRates.length;
    double [][] result = new double[N+1][N+1];
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++)
        if (i != j)
          result[i][j] = reversibleConditionalSubstitutionRates[i][j];
    for (int i = 0; i < N; i++)
      result[i][N] = absorptionsRates[i];
    fillRateMatrixDiagonalEntries(result);
    return result;
  }
  
  // note: slightly different from quasiStationaryDistribution: returns an array of size Q while other was 1 smaller (extra comp is zero anyways)
  public static double [] quasiStationaryDistributionFromRates(
      double [][] subRates, 
      double [] absRates)
  {
    if (subRates.length !=  absRates.length)
      throw new RuntimeException();
    double [][] Q = formQMtx(subRates, absRates);
    double [][] marg = MatrixFunctions.expm(new DoubleMatrix(Q)).toArray2();
    
    double [][] reversibleCondPrMtx = extractReversibleCondPrMtx(marg);
    double [] sd = findStatDistn(reversibleCondPrMtx);
    double [] result = new double[Q.length];
    for (int i = 0; i < sd.length; i++)
      result[i] = sd[i];

    return result;
  }
  
  /**
   * The provided mtx should have:
   * @param _m: _m[i][j] = P(X_next = j | X_prev = i)
   * @return
   */
  public static double[] findStatDistn(double[][] _m)  // prev -> next
  {
    Matrix m = new Matrix(_m);
    double [] numbers = topEigenvector(m.transpose());
    Multinomial.normalize(numbers);
    return numbers;
  }
  
  /**
   * For mtx with positive entries: get the top eigenvector used in 
   * Perron-Frobenius theorem
   * @param M
   * @return
   */
  public static double [] topEigenvector(final Matrix M)
  {
    final int size = M.getColumnDimension();
    final EigenvalueDecomposition ed = M.eig();
    final Matrix D = ed.getD();
    final double [] imag = ed.getImagEigenvalues();
    int argmax = -1; 
    double max = Double.NEGATIVE_INFINITY;
    for (int i = 0; i < size; i++)
      if (Math.abs(D.get(i, i)) > max && imag[i] == 0.0)
        {
          max =Math.abs( D.get(i, i));
          argmax = i;
        }
    if (argmax == -1)
      throw new RuntimeException("Bad matrix");
    // little check to make sure it's positive for convenience
    final double sign = (ed.getV().get(0,argmax) < 0.0 ? -1.0 : 1.0);
    final double [] result = new double[size];
    for (int i = 0; i < size; i++)
      result[i] = ed.getV().get(i, argmax) * sign;
    return result;
  }
  
  public static double [][] extractReversibleCondPrMtx(double [][] fullIrreversiblePrMtx)
  {
    final int N = fullIrreversiblePrMtx.length - 1;
    double [][] result = new double[N][N];
    for (int i =0 ; i < N; i++)
    {
      for (int j = 0; j < N ; j ++)
      {
        result[i][j] = fullIrreversiblePrMtx [i][j];
        if (result[i][j] < 0) 
          throw new RuntimeException();
      }
      Multinomial.normalize(result[i]);
    }
    return result;
  }
  
  /**
   * Fill in place the diagonals to be equal for each row to the negative of the
   * off diagonals. 
   * @param rate
   */
  public static void fillRateMatrixDiagonalEntries(final double [][] rate)
  {
    int size = rate.length;
    for (int i = 0; i < size; i++)
    {
      double sum = 0.0;
      for (int j = 0; j < size; j++)
        if (i!=j)
          sum += rate[i][j];
      if (rate[i][i] != 0.0)
        throw new RuntimeException();
      rate[i][i] = -sum;
    }
  }
  
  private static double[] uniformDeletionRate(double deletionRate, int nChars)
  {
    double [] result = new double[nChars];
    for (int i = 0 ; i < nChars; i++)
      result[i] = deletionRate ; // / ((double) nChars);
    return result;
  }
 
}
