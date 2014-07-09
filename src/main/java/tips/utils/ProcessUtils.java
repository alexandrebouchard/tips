package tips.utils;

import java.util.Random;

import briefj.collections.Counter;

import tips.Process;

public class ProcessUtils
{
  public static <S> double holdRate(Process<S> process, S x)
  {
    return process.rates(x).totalCount();
  }
  
  public static <S> double transitionProbability(Process<S> process, S x, S y)
  {
    Counter<S> rates = process.rates(x);
    double norm = rates.totalCount();
    return rates.getCount(y) / norm;
  }
  
  public static <S> S sample(Process<S> process, Random rand, S x)
  {
    Counter<S> rates = new Counter<S>(process.rates(x));
    rates.normalize();
    return ProposalRandom.sampleCounter(rates, rand);
  }
}
