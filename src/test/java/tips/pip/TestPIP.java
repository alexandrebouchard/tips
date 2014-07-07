package tips.pip;



import java.util.List;
import java.util.Map;

import muset.LinearizedAlignment;
import muset.MSAPoset;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jgrapht.UndirectedGraph;
import org.junit.Test;

import tips.ImportanceSampler;
import bayonet.graphs.GraphUtils;
import briefj.Indexer;
import briefj.collections.Counter;
import briefj.collections.UnorderedPair;

import com.google.common.collect.Maps;



/**
 * Test the correctness of the method using the PIP evolutionary 
 * process, which has an analytical form for the transition 
 * probability. 
 * 
 * Check that we converge to that analytical value using the TIPS
 * approximation (throws if we do not get under a relative threshold
 * or if there is not an mse decrease twice in a row when increasing
 * the number of particles by a factor of 2).
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class TestPIP
{
  private static int nParticleIncreaseRounds = 10;
  private static int nTestInstances = 100;
  private static int nTestRepeats = 100;
  private double relativeThreshold = 1e-2;

  public static void main(String [] args)
  {
    new TestPIP().analyticTest();
  }
  
  @Test
  public void analyticTest()
  {
    PIPMain pipMain = new PIPMain();
    pipMain.potentialProposalOptions.automatic = false;
    
    for (int testInstance = 0; testInstance < nTestInstances ; testInstance++)
    {
      pipMain.generateNextData();
      System.out.println("Generated test case " + testInstance);
      
      MSAPoset msa = pipMain.getGeneratedEndPoints();
      
      System.out.println(pipMain.getFullGeneratedPath());
      
      double exact = Math.exp(exact(pipMain.mu, pipMain.lambda, pipMain.bl, msa));
      System.out.println("exact transition probability: " + exact);
      double threshold = exact * relativeThreshold;
      System.out.println("threshold: " + threshold);
      
      ImportanceSampler<PIPString> is = pipMain.buildImportanceSampler();
      is.nParticles = 1;
      SummaryStatistics weightVariance = new SummaryStatistics();
      double previousMSE = Double.POSITIVE_INFINITY;
      boolean previousBad = false;
      for (int particleIncreaseRound = 0; particleIncreaseRound < nParticleIncreaseRounds ; particleIncreaseRound++)
      {
        SummaryStatistics 
          mseStat = new SummaryStatistics(),
          meanStat =new SummaryStatistics();
        for (int testRepeat = 0; testRepeat < nTestRepeats ; testRepeat++)
        {
          Counter<List<PIPString>> samples = is.sample(pipMain.getStart(), pipMain.getEnd(), pipMain.bl, weightVariance);
          double estimate = is.estimateZ(samples);
          mseStat.addValue(Math.pow((estimate - exact), 2));
          meanStat.addValue(estimate);
        }
        double currentMSE = mseStat.getMean();
        System.out.println("mse (" + is.nParticles + " particles, " + nTestRepeats + " repeats): " + currentMSE + "\tmean: " + meanStat.getMean());
        
        if (currentMSE > previousMSE)
        {
          if (previousBad)
            throw new RuntimeException("Test suspicious: twice in a row, mse did not decrease when n particles " +
            		"was doubled (delta " + (previousMSE - currentMSE));
          previousBad = true;
        }
        else
          previousBad = false;
        
        is.nParticles *= 2;
        previousMSE = currentMSE;
      }
      
      if (previousMSE > threshold )
        throw new RuntimeException("Test suspicious: mse did no go under " + threshold);
      
      System.out.println();
    }
  }
  
  /**
   * Compute exact transition probabilities, which are available in the special
   * case of PIP thanks to analytic results (see Evolutionary Inference using 
   * the Poisson indel Process, Bouchard and Jordan, 2013).
   * @param mu
   * @param lambda
   * @param bl
   * @param msa
   * @return
   */
  public static double exact(double mu, double lambda, double bl, MSAPoset msa)
  { 
    LinearizedAlignment linearizedMSA = new LinearizedAlignment(msa);
    PIPTreeNode root = PIPTreeNode.nextUnlabelled();
    PIPTreeNode 
      n1 = PIPTreeNode.withLabel(PIPMain.ta.toString()),
      n2 = PIPTreeNode.withLabel(PIPMain.tb.toString());
    UndirectedGraph<PIPTreeNode, UnorderedPair<PIPTreeNode,PIPTreeNode>> topo = GraphUtils.newUndirectedGraph();
    topo.addVertex(n1);
    topo.addVertex(n2);
    topo.addVertex(root);
    topo.addEdge(n1, root);
    topo.addEdge(n2, root);
    Map<UnorderedPair<PIPTreeNode,PIPTreeNode>, Double> bls = Maps.newLinkedHashMap();
    bls.put(UnorderedPair.of(n1, root), fraction * bl);
    bls.put(UnorderedPair.of(n2, root), (1.0-fraction)*bl);
    
    double [][] trivialRateMtx = new double[][]{{1}};
    Indexer<Character> trivialIndex = new Indexer<Character>();
    trivialIndex.addToIndex(PIPProcess.star);
    
    PoissonParameters poissonParams = new PoissonParameters(trivialIndex, trivialRateMtx, lambda, mu);
    
    PIPLikelihoodCalculator calculator = new PIPLikelihoodCalculator(poissonParams, linearizedMSA, topo, bls, root);
    double logJoint = calculator.computeDataLogProbabilityGivenTree(); 
    
    double statRate = lambda / mu;
    PoissonDistribution pd = new PoissonDistribution(statRate);
    double logPrior = Math.log(pd.probability(msa.sequences().get(PIPMain.ta).length()));
    return logJoint - logPrior;
  }
  
  // used to test correctness
  public static double fraction = 0.5;
}
