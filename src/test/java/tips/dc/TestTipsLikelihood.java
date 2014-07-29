package tips.dc;

import java.util.Random;


import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.junit.Assert;
import org.junit.Test;

import bayonet.distributions.Exponential.RateParameterization;
import blang.ForwardSampler;
import blang.ProbabilityModel;
import blang.annotations.DefineFactor;

import tips.Potential;
import tips.StationaryProcess;
import tips.bd.ReversibleBDProcess;
import tips.bd.SimpleBirthDeathPotential;
import tips.dc.TipsTreeLikelihood.LikelihoodCalculationMethod;
import tips.finite.FiniteProcess;
import tips.finite.FiniteProcessPotential;

import conifer.TopologyUtils;
import conifer.TreeNode;
import conifer.factors.NonClockTreePrior;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.DiscreteGammaMixture;
import conifer.models.MultiCategorySubstitutionModel;



public class TestTipsLikelihood
{
  
  /**
   * Use a finite rate matrix case where the true normalization is known.
   * 
   * Check that the likelihood converges to that.
   */
  @Test
  public void testAnalytical()
  {
    final int nTaxa = 5;
    final double [][] rateMatrix = new double[][]{{-1,1},{1,-1}};
    final Random genRandom = new Random(1);
    // generate data and compute its analytic logLikelihood

    UnrootedTreeLikelihood<MultiCategorySubstitutionModel<DiscreteGammaMixture>> analytic = UnrootedTreeLikelihood.createEmpty(1, TopologyUtils.syntheticTaxaList(nTaxa)).withSingleRateMatrix(rateMatrix);
    System.out.println(analytic.tree);
    analytic.generate(genRandom);
    System.out.println(analytic.observations);
    final double analyticLogPr = analytic.logDensity();
    System.out.println("Analytic pr: " + analyticLogPr);
    
    // prepare TIPS
    TipsTreeObservation<Integer> observations = new TipsTreeObservation<Integer>(1);
    for (TreeNode leaf : analytic.tree.leaves())
      observations.set(0, leaf, readObservation(analytic.observations.get(leaf)));
    
    TipsTreeLikelihood<Integer> approximateLikelihood = //new TipsTreeLikelihood<Integer>(analytic.tree, new FiniteProcess(), observations, new FiniteProcessPotential());
      TipsTreeLikelihood.fromObservations(observations).withTree(analytic.tree).withEvolutionaryProcess( new FiniteProcess(), new FiniteProcessPotential());
    
    SummaryStatistics mse = null, approximateLogLikelihood = null;
    
    for (LikelihoodCalculationMethod method : LikelihoodCalculationMethod.values())
    {
      System.out.println();
      System.out.println("========");
      System.out.println();
      approximateLikelihood.likelihoodCalculationMethod = method;
      System.out.println("Method: " + method);
      System.out.println();
      for (int nParticles : new int[]{1,10,100,1000})
      {
        System.out.println("nParticles: " + nParticles);
        mse = new SummaryStatistics();
        approximateLogLikelihood = new SummaryStatistics();
        
        for (int repeat = 0; repeat < 100; repeat++)
        {
          approximateLikelihood.setNParticles(nParticles);
          final double approx = approximateLikelihood.recomputeLogDensity();
          approximateLogLikelihood.addValue(approx);
          mse.addValue(Math.pow(approx - analyticLogPr, 2.0));
        }
        System.out.println("Approximate pr: " + approximateLogLikelihood.getMean());
        System.out.println("sd: " + approximateLogLikelihood.getStandardDeviation());
        System.out.println("rmse: " + Math.sqrt(mse.getMean()));
        System.out.println("---");
      }
      
      Assert.assertTrue(mse.getMean() < 0.05);
    }
  }
  
  public static void main(String [] args)
  {
    new TestTipsLikelihood().testBD();
  }
  
  public static class Model
  {
    private final StationaryProcess<Integer> process = ReversibleBDProcess.normalizedIntensityWithExpectedLength(10.0);
    private final Potential<Integer> potential = new SimpleBirthDeathPotential();
    private final int nTaxa = 5;
    
    @DefineFactor(onObservations = true)
    public final TipsTreeLikelihood<Integer> approximateLikelihood = new TipsTreeLikelihood<Integer>(1, TopologyUtils.syntheticTaxaList(nTaxa))
      .withEvolutionaryProcess(process, potential);
    
    @DefineFactor
    public final NonClockTreePrior<RateParameterization> treePrior = 
      NonClockTreePrior.on(approximateLikelihood.tree).withExponentialRate(0.5);
  }
  
  
  
  /**
   * Use a finite rate matrix case where the true normalization is known.
   * 
   * Check that the likelihood converges to that.
   */
  @Test
  public void testBD()
  {
    final Random genRandom = new Random(1);
    // prepare TIPS
    Model model = new Model();
    
    ProbabilityModel prModel = new ProbabilityModel(model);
    ForwardSampler forwardSampler = new ForwardSampler(prModel);
    forwardSampler.simulate(genRandom);
    
//    model.approximateLikelihood.generate(genRandom);
    System.out.println(model.treePrior.tree);
    System.out.println(model.approximateLikelihood.observations);
    
    SummaryStatistics approximateLogLikelihood = null;
    
    for (LikelihoodCalculationMethod method : LikelihoodCalculationMethod.values())
    {
      System.out.println();
      System.out.println("========");
      System.out.println();
      model.approximateLikelihood.likelihoodCalculationMethod = method;
      System.out.println("Method: " + method);
      System.out.println();
      for (int nParticles : new int[]{10,100,1000})
      {
        System.out.println("nParticles: " + nParticles);
        approximateLogLikelihood = new SummaryStatistics();
        
        for (int repeat = 0; repeat < 20; repeat++)
        {
          model.approximateLikelihood.setNParticles(nParticles);
          final double approx = model.approximateLikelihood.recomputeLogDensity();
//            useDC ? 
//              model.approximateLikelihood.recomputeLogDensity_dcSMC() :
//              model.approximateLikelihood.recomputeLogDensity_standardSMC(true);
          approximateLogLikelihood.addValue(approx);
        }
        System.out.println("Approximate pr: " + approximateLogLikelihood.getMean());
        System.out.println("sd: " + approximateLogLikelihood.getStandardDeviation());
        System.out.println("---");
      }
      
      Assert.assertTrue(approximateLogLikelihood.getMean() < 0.05);
    }
  }
  
  private int readObservation(Object o)
  {
    double [][] array = (double[][]) o;
    if (array.length != 1 || array[0].length != 2)
      throw new RuntimeException();
    return array[0][0] == 1.0 ? 0 : 1;
  }
}
