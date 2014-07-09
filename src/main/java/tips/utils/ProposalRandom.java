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

//  public int sampleMultinomial(double[] ps)
//  {
//    int index = SampleUtils.sampleMultinomial(rand, ps);
//    logProbability += Math.log(ps[index]);
//    return index;
//  }
  
  public double sampleUniform(double len)
  {
    logProbability = logProbability - Math.log(len);
    return len * rand.nextDouble();
  }

//  /**
//   * Sample interval partitions with numberOfPartitions+1 blocks in it (think of numberOfPart as num of dividing lines to sample via unif order stat)
//   * @param numberOfPartitions
//   * @param len
//   * @return the length of the numberOfPartitions+1 blocks
//   */
//  public List<Double> samplePartitions(int numberOfPartitions, double len)
//  {
//    List<Double> pts = new ArrayList<Double>();
//    pts.add(0.0);
//    for (int i = 0; i < numberOfPartitions; i++)
//      pts.add(sampleUniform(len));
//    pts.add(len);   
//    
//    Collections.sort(pts);
//    
//    if (numberOfPartitions > 0)
//      logProbability += MathUtils.logFactorial(numberOfPartitions - 1);
//    
//    List<Double> result = new ArrayList<Double>();
//    for (int i = 1; i < pts.size(); i++)
//      result.add(pts.get(i) - pts.get(i-1));
//    
//    return result;
//  }

  
}