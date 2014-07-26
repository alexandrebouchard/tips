package tips.finite;

import briefj.collections.Counter;
import tips.Process;
import tips.StationaryProcess;


/**
 * A simple example used for testing.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class FiniteProcess implements StationaryProcess<Integer>
{
  private final Counter<Integer> 
    rates0 = initRates(0),
    rates1 = initRates(1);

  public Counter<Integer> rates(Integer point)
  {
    if (point != 0 && point != 1)
      throw new RuntimeException();
    if (point == 0)
      return new Counter<Integer>(rates0);
    else
      return new Counter<Integer>(rates1);
  }

  private Counter<Integer> initRates(int start)
  {
    Counter<Integer> result = new Counter<Integer>();
    result.setCount(1 - start, 1.0);
    return result;
  }

  @Override
  public double getStationaryProbability(Integer state)
  {
    return 0.5;
  }
}
