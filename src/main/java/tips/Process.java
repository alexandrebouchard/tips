package tips;

import java.util.Random;



public interface Process<S>
{
  public double holdRate(S x);
  public double transitionProbability(S x, S y);
  public S sample(Random rand, S x);
}