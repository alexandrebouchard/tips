package tips.bd;

import java.util.Random;

import briefj.collections.Counter;
import tips.Process;
import tips.StationaryProcess;
import tutorialj.Tutorial;


/**
 * A simple example used in the tutorial documentation.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class SimpleBirthDeathProcess implements StationaryProcess<Integer>
{
  private double birthRate = 1.0;
  private double deathRate = 1.0;

  /**
   * (Note: using TIPS for a birth death process is overkill, as
   * specialized methods exist for numerical calculation of transition
   * probabilities of birth-death process, see 
   * FW Crawford and MA Suchard. Transition probabilities for general 
   * birth-death processes with applications in ecology, genetics, 
   * and evolution. J Math Biol, 65:553-580, 2012.)
   */
  @Override
  @Tutorial(showSignature = true, showLink = true, linkPrefix = "src/test/java/")
  public Counter<Integer> rates(Integer point)
  {
    Counter<Integer> result = new Counter<Integer>();
    result.setCount(point + 1, birthRate  * point + 1);
    if (point > 0)
      result.setCount(point - 1, deathRate * point);
    return result;
  }

  @Override
  public double getStationaryProbability(Integer state)
  {
    return 0.5;
  }

  @Override
  public Integer sampleFromStationary(Random rand)
  {
    return rand.nextBoolean() ? 0 : 1;
  }
}
