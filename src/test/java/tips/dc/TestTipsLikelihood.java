package tips.dc;

import java.util.Random;

import org.junit.Test;

import tips.TimeIntegratedPathSampler;
import tips.finite.FiniteProcess;
import tips.finite.FiniteProcessPotential;

import conifer.TopologyUtils;
import conifer.TreeNode;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.DiscreteGammaMixture;
import conifer.models.MultiCategorySubstitutionModel;



public class TestTipsLikelihood
{
  private final double [][] rateMatrix = new double[][]{{-1,1},{1,-1}};
  
  /**
   * Use a finite rate matrix case where the true normalization is known.
   * 
   * Check that the likelihood converges to that.
   */
  @Test
  public void testAnalytical()
  {
    // generate data and compute its analytic logLikelihood
    Random genRandom = new Random(1);
    UnrootedTreeLikelihood<MultiCategorySubstitutionModel<DiscreteGammaMixture>> analytic = UnrootedTreeLikelihood.createEmpty(1, TopologyUtils.syntheticTaxaList(5)).withSingleRateMatrix(rateMatrix);
    analytic.generate(genRandom);
    System.out.println(analytic.observations);
    System.out.println("Analytic pr: " + analytic.logDensity());
    
    // prepare TIPS
    TipsTreeObservation<Integer> observations = new TipsTreeObservation<Integer>(1);
    for (TreeNode leaf : analytic.tree.leaves())
      observations.set(0, leaf, readObservation(analytic.observations.get(leaf)));
    
    TimeIntegratedPathSampler<Integer> sampler = new TimeIntegratedPathSampler<Integer>(new FiniteProcessPotential(), new FiniteProcess());
    sampler.nParticles = 100000;
    TipsTreeLikelihood<Integer> approximateLikelihood = new TipsTreeLikelihood<Integer>(analytic.tree, new FiniteProcess(), observations, sampler);
    
//    System.out.println(observations);
    System.out.println("nParticles: " + sampler.nParticles);
    System.out.println("Approximate pr: " + approximateLikelihood.recomputeLogDensity());
  }
  
  private int readObservation(Object o)
  {
    double [][] array = (double[][]) o;
    if (array.length != 1 || array[0].length != 2)
      throw new RuntimeException();
    return array[0][0] == 1.0 ? 0 : 1;
  }
}
