package tips.utils;

import java.util.List;
import java.util.Random;
import java.util.Set;

import muset.MSAPoset;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import tips.Process;
import tips.Proposal;
import tips.TimeIntegratedPathSampler;
import tips.pip.PIPMain;
import tips.pip.PIPProcess;

import com.google.common.collect.Sets;


/**
 * Some baselines performing the same task as TIPS in more basic ways.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 */
public class Baselines
{
  /**
   * Simulate path and waiting times and see if the state at length
   * bl coincides with the target end point.
   * 
   * TODO: this is specialized to PIP for historical reasons, but 
   * could be made more general easily.
   */
  public static double standardIS(
      MSAPoset ref,
      double bl, 
      PIPProcess process, 
      int nPart, 
      Random rand, 
      SummaryStatistics weightStats)
  {
    double num = 0.0;
    
    for (int i = 0; i < nPart ; i++)
    {
      MSAPoset temp = process.createInitMSA(ref.sequences().get(PIPMain.ta));
      // simulate
      MSAPoset proposed = PIPProcess.keepOnlyEndPts(process.sample(rand, temp, bl), PIPMain.ta, PIPMain.tb);
      
      if (Sets.newLinkedHashSet(ref.edges()).equals(Sets.newLinkedHashSet(proposed.edges())) && ref.sequences().equals(proposed.sequences()))
      {
        num++;
        if (weightStats!=null) weightStats.addValue(1.0);
      }
      else
      {
        if (weightStats!=null) weightStats.addValue(0.0);
      }
    }
    
    return (num+1)/(nPart+2);
  }

  /**
   * Use the set of paths generated by the proposals and TIPS's marginalization 
   * machinery to sum the probabilities of a subset of paths, hence obtaining an 
   * approximation (from below) of the infinite set of paths.
   */
  public static <S> double exhaustiveSum(
      Random rand, 
      int nTries, 
      Process<S> process, 
      Proposal<S> proposal, 
      S x, S y, 
      double t)
  {
    double sum = 0.0;
    Set<List<S>> covered = Sets.newHashSet();
    
    for (int trial = 0; trial < nTries; trial++)
    {
      List<S> proposed = proposal.propose(rand, x, y, t).getLeft();
      if (!covered.contains(proposed))
      {
        covered.add(proposed);
        double integral = TimeIntegratedPathSampler.integral(process, proposed, t);
        for (int jIdx = 0; jIdx < proposed.size() - 1; jIdx++)
          integral *= ProcessUtils.transitionProbability(process, proposed.get(jIdx), proposed.get(jIdx+1));
        sum += integral;
      }
    }
    
    return sum;
  }

}
