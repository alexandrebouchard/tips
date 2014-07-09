package tips.bd;

import briefj.collections.Counter;
import tips.Process;
import tutorialj.Tutorial;



public class SimpleBirthDeathProcess implements Process<Integer>
{
  private double birthRate = 1.0;
  private double deathRate = 1.0;

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
}
