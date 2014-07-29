package tips.bd;

import java.util.Random;

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
  
  public static ReversibleBDProcess normalizedIntensityWithExpectedLength(double expectedLen)
  {
    final double intensity = 1.0;
    return new ReversibleBDProcess(intensity/2.0, intensity/2.0/expectedLen);
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
  
  @Override
  public Integer sampleFromStationary(Random rand)
  {
    double uniformDraw = rand.nextDouble();
    
    double massSoFar = 0.0;
    for (int i = 0; i < MAX_STATIONARY_DRAW; i++)
    {
      massSoFar += getStationaryProbability(i);
      if (massSoFar >= uniformDraw)
        return i;
    }
    
    throw new RuntimeException();
  }
  private double MAX_STATIONARY_DRAW = 100000;

  /**
   * Required by TipsTreeLikelihood
   */
  @Override
  public int hashCode()
  {
    final int prime = 31;
    int result = 1;
    long temp;
    temp = Double.doubleToLongBits(lambda);
    result = prime * result + (int) (temp ^ (temp >>> 32));
    temp = Double.doubleToLongBits(mu);
    result = prime * result + (int) (temp ^ (temp >>> 32));
    return result;
  }

  /**
   * Required by TipsTreeLikelihood
   */
  @Override
  public boolean equals(Object obj)
  {
    if (this == obj)
      return true;
    if (obj == null)
      return false;
    if (getClass() != obj.getClass())
      return false;
    ReversibleBDProcess other = (ReversibleBDProcess) obj;
    if (Double.doubleToLongBits(lambda) != Double
        .doubleToLongBits(other.lambda))
      return false;
    if (Double.doubleToLongBits(mu) != Double.doubleToLongBits(other.mu))
      return false;
    return true;
  }

}
