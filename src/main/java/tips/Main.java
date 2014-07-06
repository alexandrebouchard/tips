package tips;



import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import muset.MSAPoset;
import muset.SequenceId;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;




import briefj.collections.Counter;

import tips.pip.PIPPotential;
import tips.pip.PIPProcess;
import tips.pip.PIPString;




public class Main
{

  
  public static void main(String [] args)
  {
    double lambda =  2.0;
    double mu = 0.5;
    double bl =  0.3;
    Random genRand = new Random(10);
    SequenceId ta = new SequenceId("A"), tb = new SequenceId("B");
    
    PIPProcess process = new PIPProcess(lambda, mu);
    MSAPoset msa = PIPProcess.keepOnlyEndPts(process.sample(genRand, bl), ta, tb);
    System.out.println(msa);
    Pair<PIPString,PIPString> endPoints = getEndPoints(msa, ta, tb);
    PIPPotential pot = new PIPPotential();
    PotProposal<PIPString> prop = new PotProposal<PIPString>(process, pot, new PotPropOptions());
    ImportanceSampler<PIPString> is = new ImportanceSampler<PIPString>(prop, process);
    is.nParticles = 100;
    SummaryStatistics weightVariance = new SummaryStatistics();
    for (int i = 0; i < 10; i++)
    {
      Counter<List<PIPString>> samples = is.sample(endPoints.getLeft(), endPoints.getRight(), bl, weightVariance);
      double estimate = is.estimateZ(samples);
      System.out.println(estimate);
    }
  }
  
  public static Pair<PIPString, PIPString> getEndPoints(MSAPoset msa, SequenceId ta,
      SequenceId tb)
  {
    if (msa.sequences().size() != 2)
      throw new RuntimeException();
    Set<Integer> aligned1 = new HashSet<Integer>(), aligned2 = new HashSet<Integer>();
    for (muset.MSAPoset.Column c : msa.columns())
      if (c.getPoints().size() == 2)
      {
        if (c.getPoints().containsKey(ta)) aligned1.add(c.getPoints().get(ta));
        if (c.getPoints().containsKey(tb)) aligned2.add(c.getPoints().get(tb));
      }
    
    List<Integer> cur; PIPString pips1, pips2;
    
    cur = new ArrayList<Integer>(); for (int i = 0; i < msa.sequences().get(ta).length(); i++) cur.add(aligned1.contains(i) ? 0 : -1); pips1 = new PIPString(cur);
    cur = new ArrayList<Integer>(); for (int i = 0; i < msa.sequences().get(tb).length(); i++) cur.add(aligned2.contains(i) ? 0 : +1); pips2 = new PIPString(cur);
    
    return Pair.of(pips1, pips2);
  }

}
