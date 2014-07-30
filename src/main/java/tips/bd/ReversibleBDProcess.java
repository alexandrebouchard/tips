package tips.bd;

import java.util.Random;

import bayonet.math.SpecialFunctions;
import blang.annotations.FactorArgument;
import blang.annotations.FactorComponent;
import blang.variables.RealVariable;
import briefj.collections.Counter;
import tips.StationaryProcess;


/**
 * A birth death process where the rate of insertion is constant,
 * and the rate of deletion is proportional to the size.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class ReversibleBDProcess<P extends ReversibleBDProcess.Parameterization> implements StationaryProcess<Integer>
{
  @FactorComponent
  public final P parameters;
  
  public static interface Parameterization
  {
    public double getMu();
    public double getLambda();
  }
  
  public static class ExpectedLengthParameterization implements ReversibleBDProcess.Parameterization
  {
    @FactorArgument
    public final RealVariable expectedLength;
    
    private static final double intensity = 1.0;

    private ExpectedLengthParameterization(double expectedLength)
    {
      this.expectedLength = new RealVariable(expectedLength);
    }

    @Override
    public double getMu()
    {
      return intensity/2.0/expectedLength.getValue();
    }

    @Override
    public double getLambda()
    {
      return intensity/2.0;
    }

    @Override
    public String toString()
    {
      return "ExpectedLengthParameterization(expectedLength=" + expectedLength
          + ")";
    }
    
    
  }
  
  public static class FullParameterization implements ReversibleBDProcess.Parameterization
  {
    @FactorArgument
    public final RealVariable mu;

    @FactorArgument
    public final RealVariable lambda;
    
    private FullParameterization(double mu, double lambda)
    {
      this.mu = new RealVariable(mu);
      this.lambda = new RealVariable(lambda);
    }

    @Override
    public double getMu()
    {
      return mu.getValue();
    }

    @Override
    public double getLambda()
    {
      return lambda.getValue();
    }

    @Override
    public String toString()
    {
      return "FullParameterization(mu=" + mu + ", lambda=" + lambda + ")";
    }
  }

  /**
   * 
   * @param lambda Insertion rate.
   * @param mu Deletion rate.
   */
  public ReversibleBDProcess(P parameters)
  {
    this.parameters = parameters;
  }
  
  public static ReversibleBDProcess<ExpectedLengthParameterization> normalizedIntensityWithExpectedLength(double expectedLen)
  {
    return new ReversibleBDProcess<ExpectedLengthParameterization>(new ExpectedLengthParameterization(expectedLen));
  }
  
  public static ReversibleBDProcess<FullParameterization> muLambdaParameterized(double mu, double lambda)
  {
    return new ReversibleBDProcess<FullParameterization>(new FullParameterization(mu, lambda));
  }
  
  private double lambda()
  {
    return parameters.getLambda();
  }
  
  private double mu()
  {
    return parameters.getMu();
  }
  
  public static class InvalidParametersException extends RuntimeException
  {
    private static final long serialVersionUID = 1L;

    public InvalidParametersException(Object bad)
    {
      super("Invalid parameter: " + bad);
    }
  }

  @Override
  public Counter<Integer> rates(Integer point)
  {
    if (mu() < 0 || lambda() < 0)
      throw new InvalidParametersException(parameters);
    Counter<Integer> result = new Counter<Integer>();
    result.setCount(point + 1, lambda());
    if (point > 0)
      result.setCount(point - 1, mu() * point);
    return result;
  }

  @Override
  public double getStationaryProbability(Integer state)
  {
    final double rate = lambda() / mu();
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
  
  @Override
  public String toString()
  {
    return "ReversibleBDProcess(" + parameters + ")";
  }

  /**
   * Required by TipsTreeLikelihood
   */
  @Override
  public int hashCode()
  {
    final int prime = 31;
    int result = 1;
    long temp;
    temp = Double.doubleToLongBits(lambda());
    result = prime * result + (int) (temp ^ (temp >>> 32));
    temp = Double.doubleToLongBits(mu());
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
    if (Double.doubleToLongBits(lambda()) != Double
        .doubleToLongBits(other.lambda()))
      return false;
    if (Double.doubleToLongBits(mu()) != Double.doubleToLongBits(other.mu()))
      return false;
    return true;
  }

}
