package tips.pip;

import java.util.Map;

import bayonet.math.NumericalUtils;



/**
 * Supporting functions for PIP likelihood computation
 * @author bouchard
 *
 */
public class PIPLikelihoodUtils
{
  public static double[][] standardLogScaleFelsensteinRecursion(
      final double[][] fels1,    final double[][] fels2, 
      final double[][] logPr1,   final double[][] logPr2)
  {
    final int ns = fels1.length;
    if (fels2.length != ns) throw new RuntimeException();
    final int ncs = logPr1.length;
    if (logPr1.length != ncs) throw new RuntimeException();
    for (double [][] item : new double[][][]{logPr1, logPr2})
      if (item[0].length != item.length)
        throw new RuntimeException();
    // resulting backward scores are indexed: site -> character, in log scale
    final double [][] result = new double[ns][ncs];
    
    for (int s = 0; s < ns; s++)
    {
      final double [] curFHat1 = fels1[s];
      final double [] curFHat2 = fels2[s];
      for (int v = 0; v < ncs; v++)
      {
        double prod = 0.0; 
        // take into account both left and right branches
        prod += sum(logPr1[v], curFHat1); // log scale!
        prod += sum(logPr2[v], curFHat2);
        result[s][v] = prod;
      }
    }
    return result;
  }
  
  private static double sum(final double [] logPrs, final double[] cacheAtLeaf)
  {
    double sum = Double.NEGATIVE_INFINITY;
    final int ncs = logPrs.length;
    for (int w = 0; w < ncs ; w++) 
      sum = NumericalUtils.logAdd(sum, logPrs[w] + cacheAtLeaf[w]);
    return sum;
  }

  public static double totalTreeLength(Map<?,Double> bls)
  {
    double sum = 0.0;
    for (double value : bls.values())
      sum += value;
    return sum;
  }
}
