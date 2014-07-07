package tips.pip;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Random;

import muset.MSAPoset;
import muset.MSAPoset.Column;
import muset.util.Edge;
import muset.SequenceId;

import org.apache.commons.lang3.tuple.Pair;

import tips.Process;
import tips.ProposalRandom;
import tips.SparseProcess;
import bayonet.distributions.Exponential;
import bayonet.math.SpecialFunctions;
import briefj.collections.Counter;


public class PIPProcess implements SparseProcess<PIPString>, Process<PIPString>
{
  public final double lambda, mu;
  
  public PIPProcess(double lambda, double mu)
  {
    super();
    this.lambda = lambda;
    this.mu = mu;
  }

  @Override
  public double holdRate(PIPString x)
  {
    return rates(x).totalCount();
  }

  @Override
  public double transitionProbability(PIPString x, PIPString y)
  {
    Counter<PIPString> rates = rates(x);
    rates.normalize();
    return rates.getCount(y);
  }

  @Override
  public PIPString sample(Random rand, PIPString x)
  {
    Counter<PIPString> rates = rates(x);
    rates.normalize();
    return ProposalRandom.sampleCounter(rates, rand);
  }
  
  public MSAPoset sampleStationary(Random rand)
  {
    
    return createInitMSA(repeat(star, (int)samplePoisson(rand, lambda/mu)));
  }
  
  // Copied from numerical recipes 
  private static double oldm = -1, g, sq, alxm;
  public static double samplePoisson(Random random, double rate) {
    double xm = rate;
    double em, t, y;

    if (xm < 12.0) {
      if (xm != oldm) {
        oldm=xm;
        g=Math.exp(-xm);
      }
      em = -1;
      t=1.0;
      do {
        em += 1.0;
        t *= random.nextDouble();
      } while (t > g);
    } 
    else {
      if (xm != oldm) {
        oldm=xm;
        sq=Math.sqrt(2.0*xm);
        alxm=Math.log(xm);
        g=xm*alxm-SpecialFunctions.lnGamma(xm+1.0);
      }
      do 
      {
        do 
        {
          y=Math.tan(Math.PI*random.nextDouble());
          em=sq*y+xm;
        } while (em < 0.0);
        em=Math.floor(em);
        t=0.9*(1.0+y*y)*Math.exp(em*alxm-SpecialFunctions.lnGamma(em+1.0)-g);
      } while (random.nextDouble() > t);
    }
    return (int)em;
  }
  
  public MSAPoset createInitMSA(String str)
  {
    SequenceId firstId = indexedTaxon(0);
    Map<SequenceId,String> seqns = new LinkedHashMap<SequenceId, String>();
    seqns.put(firstId, str);
    return new MSAPoset(seqns);
  }
  
  public static char star = '*';
  public static String repeat(char star, int nRep)
  {
    StringBuilder result =new StringBuilder();
    for (int i = 0; i < nRep; i++)
      result.append(star);
    return result.toString();
  }
  
  public static SequenceId indexedTaxon(int idx)
  {
    String idxStr = "" + idx;
    int initLen = idxStr.length();
    for (int i = 0; i < 4 - initLen; i++)
      idxStr = '0' + idxStr;
    return new SequenceId("seq-" + idxStr);
  }
  
  // generate some
  // check their likelihood
  // handle empty case
  // test case
  
  public MSAPoset sample(Random rand, double branchLength)
  {
    MSAPoset init = sampleStationary(rand);
    return sample(rand, init, branchLength);
  }
  
  public static final int MAX_FWD_SAMPLE_STEPS  = 100000000;
  public MSAPoset sample(Random rand, MSAPoset init, double branchLength)
  {
    if (branchLength <= 0.0)
      throw new RuntimeException();
    MSAPoset curMSA = init;
    double lengthConsumed = 0.0;
    
    mainLoop:for (int i =0 ; i < MAX_FWD_SAMPLE_STEPS; i++)
    {
      Pair<MSAPoset,Double> sampled = sample(rand, curMSA);
      
      lengthConsumed += sampled.getRight();
      
      if (lengthConsumed > branchLength)
        break mainLoop;
      
      curMSA = sampled.getLeft();
    }
    
    return curMSA;
  }
  
  public static MSAPoset keepOnlyEndPts(MSAPoset msa, SequenceId first, SequenceId second)
  {
    SequenceId tF = indexedTaxon(0), tL = indexedTaxon(msa.sequences().size()-1);
    Map<SequenceId,String> newSeqns = new LinkedHashMap<SequenceId, String>();
    newSeqns.put(first,  msa.sequences().get(tF));
    newSeqns.put(second, msa.sequences().get(tL));
    MSAPoset result = new MSAPoset(newSeqns);
    for (Column c : msa.columns())
      if (c.getPoints().containsKey(tF) && c.getPoints().containsKey(tL))
        if (!result.tryAdding(new Edge(c.getPoints().get(tF), c.getPoints().get(tL), first, second)))
          throw new RuntimeException();
    return result;
  }
  
  public Pair<MSAPoset,Double> sample(Random rand, MSAPoset current)
  {
    SequenceId lastId = indexedTaxon(current.sequences().size()-1);
    String lastSeq = current.sequences().get(lastId);
    int len = lastSeq.length();
    double rate = lambda + mu * len;
    double insPr = lambda / rate;
    double time = Exponential.generate(rand, rate); 
    
    boolean isIns = sampleBern(insPr, rand);
    
    Map<SequenceId, String> seqns  = new LinkedHashMap<SequenceId, String>(current.sequences());
    String newSeq = repeat(star, len + (isIns ? +1 : -1));
    SequenceId newId = indexedTaxon(current.sequences().size());
    seqns.put(newId, newSeq);
    int pos = rand.nextInt(len + (isIns ? +1 : 0));
    
    MSAPoset newMSA = new MSAPoset(seqns);
    for (Column c : current.columns())
      if (!newMSA.tryAdding(c))
        throw new RuntimeException();
    
    for (int i = 0; i < pos; i++)
      if (!newMSA.tryAdding(new Edge(i, i, lastId, newId)))
        throw new RuntimeException();
    
    if (isIns)
    {
      for (int i = pos; i < len; i++)
        if (!newMSA.tryAdding(new Edge(i, i+1, lastId, newId)))
          throw new RuntimeException();
    }
    else
    {
      for (int i = pos+1; i < len; i++)
        if (!newMSA.tryAdding(new Edge(i, i-1, lastId, newId)))
          throw new RuntimeException();
    }
    
    return Pair.of(newMSA, time);
  }
  
  public static boolean sampleBern(final double prToBeTrue, Random rand)
  {
    return rand.nextDouble() < prToBeTrue;
  }
  

  @Override
  public Counter<PIPString> rates(PIPString point)
  {
    Counter<PIPString> rates = new Counter<PIPString>();
      
    // add ins
    double nInsPoints = point.characters.size()+1;
    for (int i = 0; i < nInsPoints; i++)
      rates.incrementCount(new PIPString(point.characters.subList(0, i), +1, point.characters.subList(i, point.characters.size())), lambda / nInsPoints);
      
    // add dels
    for (int i = 0; i < point.characters.size(); i++)
      rates.incrementCount(new PIPString(point.characters.subList(0, i), point.characters.subList(i+1, point.characters.size())), mu);   
      
    return rates;
  }
  
}