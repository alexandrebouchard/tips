package tips.utils;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import bayonet.distributions.Multinomial;
import briefj.collections.Counter;





public class ProposalRandom
{
  public double logProbability = 0.0;
  public double getLogProbability()
  {
    return logProbability;
  }

  public final Random rand;
  
  public ProposalRandom(Random rand)
  {
    this.rand = rand;
  }
  
  public <S> S sampleMultinomial(Counter<S> counter)
  {
    S item = sampleCounter(counter, rand);
    logProbability += Math.log(counter.getCount(item));
    return item;
  }
  
  public static <S> S sampleCounter(Counter<S> items, Random rand)
  {
    double [] prs = new double[items.size()];
    List<S> keys = new ArrayList<S>(items.size());
    int i =0 ;
    for (S key : items.keySet())
    {
      prs[i] = items.getCount(key);
      keys.add(key);
      i++;
    }
    int index = 
      Multinomial.sampleMultinomial(rand, prs);
    return keys.get(index);
  }

  public int sampleDiscreteUniform(int n)
  {
    logProbability = logProbability - Math.log(n);
    return rand.nextInt(n);
  }
  
  public boolean sampleBern(double p)
  {
    double u = rand.nextDouble();
    if (u < p)
    {
      logProbability += Math.log(p);
      return true;
    }
    else
    {
      logProbability += Math.log(1.0 - p);
      return false;
    }
  }
  
  public double sampleUniform(double len)
  {
    logProbability = logProbability - Math.log(len);
    return len * rand.nextDouble();
  }
}