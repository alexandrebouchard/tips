package tips.bd;

import bayonet.math.SpecialFunctions;
import briefj.collections.Counter;
import tips.StationaryProcess;


/**
 * A birth death process where the rate of insertion is constant,
 * and the rate of deletion is proportional to the size.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class ReversibleBDProcess implements StationaryProcess<Integer>
{
  private final double lambda, mu;

  /**
   * 
   * @param lambda Insertion rate.
   * @param mu Deletion rate.
   */
  public ReversibleBDProcess(double lambda, double mu)
  {
    this.lambda = lambda;
    this.mu = mu;
  }

  @Override
  public Counter<Integer> rates(Integer point)
  {
    Counter<Integer> result = new Counter<Integer>();
    result.setCount(point + 1, lambda);
    if (point > 0)
      result.setCount(point - 1, mu * point);
    return result;
  }

  @Override
  public double getStationaryProbability(Integer state)
  {
    final double rate = lambda / mu;
    return Math.exp(-rate + state * Math.log(rate) - SpecialFunctions.logFactorial(state));
  }

}
