package tips.finite;

import tips.Potential;


/**
 * A finite example used for testing.
 *  
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class FiniteProcessPotential implements Potential<Integer>
{
  public double get(Integer proposed, Integer target)
  {
    return proposed.equals(target) ? 0 : 1;
  }

}
