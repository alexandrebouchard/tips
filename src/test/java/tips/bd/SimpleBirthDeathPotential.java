package tips.bd;

import tips.Potential;



public class SimpleBirthDeathPotential implements Potential<Integer>
{

  @Override
  public double get(Integer proposed, Integer target)
  {
    return Math.abs(proposed - target);
  }

}
