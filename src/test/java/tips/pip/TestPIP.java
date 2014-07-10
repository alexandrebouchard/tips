package tips.pip;



import java.util.Map;

import muset.LinearizedAlignment;
import muset.MSAPoset;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jgrapht.UndirectedGraph;
import org.junit.Test;

import tips.TimeIntegratedPathSampler;
import bayonet.graphs.GraphUtils;
import briefj.Indexer;
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
 * the number of particles by a multiplicative factor).
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class TestPIP
{
  private static int nParticleIncreaseRounds = 3;
  private static int nTestInstances = 10;
  private static int nTestRepeats = 3;
  private double relativeThreshold = 1e-3;
  private double multiplicativeFactor = 1e2;

  public static void main(String [] args)
  {
    new TestPIP().analyticTest();
  }
  
  @Test
  public void analyticTest()
  {
    PIPMain pipMain = new PIPMain();
    pipMain.ensureLinearizationUnique = true;
    pipMain.lambda =  2.0;
    pipMain.mu = 0.5;
    pipMain.bl =  0.3;
    
    for (int testInstance = 0; testInstance < nTestInstances ; testInstance++)
    {
      pipMain.generateNextData();
      System.out.println("Generated test case " + testInstance);
      MSAPoset msa = (pipMain.getGeneratedEndPoints());
      
      double exact = Math.exp(exact(pipMain.mu, pipMain.lambda, pipMain.bl, msa));
      double threshold = relativeThreshold * exact;
      System.out.println("threshold: " + threshold);
      System.out.println("exact transition probability: " + exact);
      
      TimeIntegratedPathSampler<PIPString> is = pipMain.buildImportanceSampler();
      is.nParticles = 1;
      double previousMSE = Double.POSITIVE_INFINITY;
      boolean previousBad = false;
      for (int particleIncreaseRound = 0; particleIncreaseRound < nParticleIncreaseRounds ; particleIncreaseRound++)
      {
        SummaryStatistics mseStat = new SummaryStatistics();
        for (int testRepeat = 0; testRepeat < nTestRepeats ; testRepeat++)
        {
          double estimate = is.estimateTransitionPr(pipMain.getStart(), pipMain.getEnd(), pipMain.bl);
          mseStat.addValue(Math.pow((estimate - exact), 2));
        }
        double currentMSE = mseStat.getMean();
        System.out.println("mse (" + is.nParticles + " particles, " + nTestRepeats + " repeats): " + currentMSE);
        
        if (currentMSE > previousMSE)
        {
          if (previousBad)
            throw new RuntimeException("Test suspicious: twice in a row, mse did not decrease when n particles " +
            		"was doubled (delta " + (previousMSE - currentMSE));
          previousBad = true;
        }
        else
          previousBad = false;
        
        is.nParticles *= multiplicativeFactor;
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
  
  // used to test correctness in early tests (get same value for each value in (0,1) by reversibility
  private static final double fraction = 0.5;
}
