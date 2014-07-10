package tips.bd;

import tips.Potential;
import tutorialj.Tutorial;


/**
 * A simple example used in the tutorial documentation.
 *  
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class SimpleBirthDeathPotential implements Potential<Integer>
{

  @Override
  @Tutorial(showSignature = true, showLink = true, linkPrefix = "src/test/java/")
  public double get(Integer proposed, Integer target)
  {
    // Just return how far we are from the target.
    return Math.abs(proposed - target);
  }

}
