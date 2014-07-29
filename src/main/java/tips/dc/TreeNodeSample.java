package tips.dc;

import java.util.List;
import java.util.Random;

import bayonet.distributions.Multinomial;
import briefj.BriefMath;
import briefj.collections.Counter;

import com.google.common.collect.Lists;



public class TreeNodeSample<S>
{
  private final List<S> samples = Lists.newArrayList();
  private final double [] prs;
  final double logNormalization;
  TreeNodeSample(Counter<S> samples, double logNormalization)
  {
    BriefMath.checkCloseAndFinite(samples.totalCount(), 1.0);
    this.prs = new double[samples.size()];
    int i = 0;
    for (S key : samples.keySet())
    {
      prs[i++] = samples.getCount(key);
      this.samples.add(key);
    }
    this.logNormalization = logNormalization;
  }
  public S sample(Random rand)
  {
    return samples.get(Multinomial.sampleMultinomial(rand, prs));
  }
  public Counter<S> asCounter()
  {
    Counter<S> result = new Counter<S>();
    for (int i = 0; i < samples.size(); i++)
      result.incrementCount(samples.get(i), prs[i]);
    return result;
  }
  @Override
  public String toString()
  {
    return asCounter().toString();
  }
}