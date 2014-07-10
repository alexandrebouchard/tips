package tips.utils;



import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

import tips.Potential;
import tips.Proposal;
import tips.Process;

import bayonet.math.NumericalUtils;
import briefj.collections.Counter;




/**
 * Implements Equation (2) in Monir Hajiaghayi, 
 * Bonnie Kirkpatrick, Liangliang Wang and Alexandre Bouchard-Côté. 
 * (2014) Efficient Continuous-Time Markov Chain Estimation.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 * @param <S>
 */
public class PotProposal<S> implements Proposal<S>
{
  private final Process<S> process;
  private final Potential<S> potential;
  private final PotPropOptions options;
  
  public PotProposal(Process<S> process, Potential<S> potential,
      PotPropOptions options)
  {
    super();
    this.process = process;
    this.potential = potential;
    this.options = options;
  }
  
  public Pair<List<S>, Double> propose(Random rand, S x, S y, double t)
  {
    ProposalRandom pRand = new ProposalRandom(rand);
    double greed = options.greed;
    double stopPr = options.stopPr;
    if (options.automatic)
    {
      greed = options.randPr(rand);
      stopPr = options.randPr(rand);
    }

    List<S> proposed = propose(process, potential, pRand, x, y, options.specialSymbol, greed, stopPr);
    
    return Pair.of(proposed, Math.exp(pRand.getLogProbability()));
  }
  
  public static <S> List<S> propose(Process<S> process, Potential<S> pot, ProposalRandom pRand, S firstEndPoint, S lastEndPoint, Object specialSymbol, double greed, double stopPr)
  {
    List<S> result = new ArrayList<S>();
    
    S current = firstEndPoint;
    result.add(current);

    if (current.equals(lastEndPoint))
      if (pRand.sampleBern(stopPr))
        return result;
    
    boolean success = false;
    mainLoop:for (int i = 0; i < MAX_PROPOSE_ATTEMPTS; i++)
    {
      Counter<S> trPr = process.rates(current);
      trPr.normalize();
      buildPotProposal(pot, trPr, current, lastEndPoint, pRand.rand, specialSymbol, greed, stopPr);
      S next = pRand.sampleMultinomial(trPr);
      if (next.equals(specialSymbol))
      {
        result.add(lastEndPoint);
        success = true;
        break mainLoop;
      }
      else
      {
        result.add(next);
        current = next;
      }
    }
    if (!success)
      throw new RuntimeException();
    
    return result;
  }
  
  public static final String SPECIAL_SYMBOL = "___TARGET___";
  private static final int MAX_PROPOSE_ATTEMPTS = 1000000;

  public static <S> void buildPotProposal(Potential<S> pot, Counter<S> transitionPrs, S current, S target, Random rand, Object specialSymbol, double greed, double stopPr)
  {
    potPropDistort(transitionPrs, pot, current, target, greed);
    addTarget(transitionPrs, target, specialSymbol, stopPr);
  }
  
  @SuppressWarnings("unchecked")
  public static <S> void addTarget(@SuppressWarnings("rawtypes") Counter transitionPrs, S target, Object specialSymbol, double stopPr)
  {
    if (transitionPrs.keySet().contains(specialSymbol))
      throw new RuntimeException();
    
    final double curPr = transitionPrs.getCount(target);
    
    if (curPr == 0.0)
      return;
    
    transitionPrs.setCount(specialSymbol, curPr * stopPr);
    transitionPrs.setCount(target, curPr * (1.0 - stopPr));
  }
  
  public static <S> void potPropDistort(Counter<S> transitionPrs, Potential<S> pot, S current, S target, double greed)
  {
    NumericalUtils.checkIsClose(1.0, transitionPrs.totalCount());
    final double initialPotential = pot.get(current, target);
    if (current == null || target == null)
      throw new RuntimeException();
    if (greed < 0.5 || greed > 1.0)
      throw new RuntimeException();
    
    double pGood = 0.0, pBad = 0.0;
    for (S key : transitionPrs.keySet())
    {
      double delta = pot.get(key, target) - initialPotential;
      double pr = transitionPrs.getCount(key);
      if (delta == -1) // makes us closer to the target
        pGood += pr;
      else if (delta == +1 || delta == 0.0)
        pBad += pr;
      else if (Double.isInfinite(delta))
        ;
      else
        throw new RuntimeException(
            "\nd(" + key + "," + target + ") = " + pot.get(key, target) + "\n" +
            "d(" + current + "," + target + ") = " + pot.get(current, target));
      
      if (Double.isInfinite(delta))
        transitionPrs.setCount(key, 0.0);
    }
    
    if (pGood == 0.0 || pBad == 0.0)
    {
      transitionPrs.normalize();
      return;
    }
    
    double alpha = Math.max(greed, pGood);
    
    for (S key : transitionPrs.keySet())
    {
      double delta = pot.get(key, target) - initialPotential;
      double pr = transitionPrs.getCount(key);
      double newValue;
      if (Double.isInfinite(delta))
        newValue = 0.0;
      else
        newValue = (delta == -1 ? alpha/pGood : (1.0 - alpha)/pBad ) * pr;
      transitionPrs.setCount(key, newValue);
    }
    transitionPrs.normalize();
  }
}
