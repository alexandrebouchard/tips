package tips;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

import bayonet.distributions.Exponential;
import briefj.BriefLists;

import tips.utils.ProcessUtils;




public class ForwardSampler<S>
{
  private final Process<S> process;
  
  public ForwardSampler(Process<S> process)
  {
    this.process = process;
  }

  private int MAX_SAMPLE_TRIES = 1000000;
  
  public S forwardSample(Random rand, S x, double t)
  {
    List<S> result = new ArrayList<S>();
    result.add(x);
    double s = 0.0;
    boolean success = false;
    S current = x;
    mainLoop:for  (int i = 0; i < MAX_SAMPLE_TRIES; i++)
    {
      Pair<S,Double> sample = sampleJumpWaitingTime(rand, current); 
      s += sample.getRight();
      if (s > t)
      {
        success = true;
        break mainLoop;
      }
      else
      {
        current = sample.getLeft();
        result.add(current);
      }
    }
    if (!success)
      throw new RuntimeException();
    return BriefLists.last(result);
  }
  
  private Pair<S,Double> sampleJumpWaitingTime(Random rand, S x)
  {
    final double holdRate = ProcessUtils.holdRate(process, x);
    final double sampledTime = Exponential.generate(rand, holdRate);
    S nextState = ProcessUtils.sample(process, rand, x);
    return Pair.of(nextState, sampledTime);
  }
}
