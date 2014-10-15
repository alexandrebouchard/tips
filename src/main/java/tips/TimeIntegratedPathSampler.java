package tips;

import java.util.Collections;
import java.util.List;
import java.util.Random;


import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;

import com.google.common.collect.Lists;

import tips.utils.PotPropOptions;
import tips.utils.PotProposal;
import tips.utils.ProcessUtils;
import tutorialj.Tutorial;


import briefj.collections.Counter;
 


/**
 * The main algorithm of this package: an IS algorithm for CTMCs where 
 * sojourn time are analytically marginalized.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 * @param <S>
 */
public class TimeIntegratedPathSampler<S>
{
  /**
   * 
   */
  public int nParticles = 1000;
  
  private final Proposal<S> proposal;
  private final Process<S> process;
  
  /**
   * Source of randomness used by the algorithm.
   * 
   * By default, the seed is fixed to 1.
   */
  public Random rand = new Random(1);
  
  /**
   * Unnormalized weight statistics
   * to the provided SummaryStatistics object so that diagnostic can be read off 
   * from that object (in particular, to get weights variance).
   * 
   * Leave null if not needed.
   */
  public SummaryStatistics unnormalizedWeightsStatistics = null;
  
  /**
   * Create a TIPS algorithm with a generic proposal mechanism.
   * 
   * @param proposal
   * @param process
   */
  public TimeIntegratedPathSampler(Proposal<S> proposal, Process<S> process)
  {
    this.proposal = proposal;
    this.process = process;
  }
  
  /**
   * Create a TIPS algorithm with a PotProposal, the provided potential, 
   * and default value for PotProposal's options.
   * 
   * (Syntactic sugar).
   * 
   * @param potential
   * @param process
   */
  public TimeIntegratedPathSampler(Potential<S> potential, Process<S> process)
  {
    this.proposal = new PotProposal<S>(process, potential, new PotPropOptions());
    this.process = process;
  }
  
  /**
   * 
   * @param x
   * @param y
   * @param t1
   * @param t2
   * @return
   */
  public Counter<S> sampleTreeCherries(S x, S y, double t1, double t2)
  {
    Counter<S> result = new Counter<S>();
    
    for (int particleIndex = 0; particleIndex < nParticles; particleIndex++)
    {
      Pair<S, Double> sampleTreeCherry = sampleTreeCherry(x, y, t1, t2);
      result.incrementCount(sampleTreeCherry.getLeft(), sampleTreeCherry.getRight());
    }
    
    return result;
  }
  
  public Pair<S,Double> sampleTreeCherry(S x, S y, double t1, double t2)
  {
    // sample a path between the two end points (start does not matter by reversibility)
    Pair<List<S>, Double> proposed = proposal.propose(rand, x, y, t1 + t2);
    final List<S> fullPath = proposed.getLeft();
    
    // sample a root 
    final int rootIndex = rand.nextInt(fullPath.size());
    
    // build the two paths
    List<S> list1 = Lists.newArrayList(fullPath.subList(0, rootIndex + 1)); Collections.reverse(list1);
    List<S> list2 = fullPath.subList(rootIndex, fullPath.size());
    
    // weight computation
    final S root = list2.get(0);
    final double weight = 
        (getStationaryPr(root) * computeUnnormalizedTargetPr(list1, t1) * computeUnnormalizedTargetPr(list2, t2))
      / (proposed.getRight() * (1.0 / fullPath.size()));
    
    return Pair.of(root, weight);
  }

  private StationaryProcess<S> stationary = null;
  public double getStationaryPr(S s)
  {
    if (stationary == null)
      stationary = (StationaryProcess<S>) process;
    return stationary.getStationaryProbability(s);
  }

//  public Pair<S,Double> sampleTreeCherry(S x, S y, double t1, double t2)
//  {
//    // sample a path between the two end points (start does not matter by reversibility)
//    Pair<List<S>, Double> proposed = proposal.propose(rand, x, y, t1 + t2);
//    final List<S> fullPath = proposed.getLeft();
//    
//    // sample a root 
//    final List<S> reversed = Lists.newArrayList(fullPath); Collections.reverse(reversed);
//    final double [] 
//      forwardJumpPrs = jumpPrs(fullPath),
//      backwardJumpPrs= jumpPrs(reversed);
//    final double [] prs = new double[fullPath.size()];
//    for (int i = 0; i < prs.length; i++)
//    {
//      
//    }
//    
//    double weight = 
//      proposed.getRight() +    // proposing the jump chain
//      xxx; // proposing a rooting
//  }
//  
//  private double[] jumpPrs(List<S> reversed)
//  {
//    double [] result = 
//    return null;
//  }

  /**
   * Estimate end-point transition probability.
   * 
   * @param x start point
   * @param y end point
   * @param t time between the two end points
   * @return Estimate for the end-point transition probability P(X_t = y|X_0 = x)
   */
  public double estimateTransitionPr(S x, S y, double t)
  {
    double sum = (Double) runTIPS(x, y, t, false);
    return sum/((double) nParticles);
  }
  
  /**
   * Extract samples from X_{0:t} | X_0, X_t.
   * 
   * @param x start point
   * @param y end point
   * @param t time between the two end points
   * @return An approximate unnormalized discrete measure built from the sampler. 
   *   The keys are paths, and the values, unnormalized weights. Note that identical 
   *   particles are grouped together by adding their weights.
   */
  @SuppressWarnings({ "unchecked", "rawtypes" })
  public Counter<List<S>> sampleEndPointConditionedPaths(S x, S y, double t)
  {
    return (Counter) runTIPS(x, y, t, true);
  }
  

  
  /**
   * First, runTIPS(), shown below, which is just an IS algorithm.
   * 
   * The argument ``keepPath`` controls whether sampled paths should be kept or not.
   * If true, then return a Counter over paths; if false, return 
   * just the sum of the weights.
   * The latter is useful because it runs in constant memory and still allows for
   * computing the transition probability estimate.
   */
  @Tutorial(showSignature = true, showLink = true)
  private Object runTIPS(S x, S y, double time, boolean keepPath)
  {
    Counter<List<S>> result = keepPath ? new Counter<List<S>>() : null;
    double sum = 0.0;
    
    for (int particleIndex = 0; particleIndex < nParticles; particleIndex++)
    {
      // propose
      Pair<List<S>, Double> proposed = proposal.propose(rand, x, y, time);
      
      // compute weight
      final double computeUnormalizedTargetPr = computeUnnormalizedTargetPr(proposed.getLeft(), time);
      
      if (unnormalizedWeightsStatistics != null) 
        unnormalizedWeightsStatistics.addValue(computeUnormalizedTargetPr/proposed.getRight());
      
      final double weight = computeUnormalizedTargetPr/proposed.getRight();
      
      if (keepPath)
        result.incrementCount(proposed.getLeft(), weight);
      else
        sum += weight;
    }
    
    return keepPath ? result : sum;
  }
  
  /**
   * Second, computeUnnormalizedTargetPr(), shown below, which compute the unnormalized target pr, 
   * i.e. the pr of the jumps time the 
   * integrated waiting times.
   */
  @Tutorial(showSignature = true, showLink = true)
  private double computeUnnormalizedTargetPr(List<S> proposedJumps, double time)
  {
    double weight = marginalizeSojournTimes(process, proposedJumps, time);
    for (int jumpIndex = 0; jumpIndex < proposedJumps.size() - 1; jumpIndex++)
      weight *= ProcessUtils.transitionProbability(process, proposedJumps.get(jumpIndex), proposedJumps.get(jumpIndex+1));
    return weight;
  }
  
  /**
   * And finally, marginalizeSojournTimes(), shown below, implementing Proposition 2 in
   * the paper.
   */
  @Tutorial(showSignature = true, showLink = true)
  public static <S> double marginalizeSojournTimes(Process<S> process, List<S> proposed, double t)
  {
    // build matrix
    final int size = proposed.size() + 1;
    double [][] mtx = new double[size][size];
    for (int i = 0; i < proposed.size(); i++)
    {
      final double curRate = ProcessUtils.holdRate(process, proposed.get(i));
      mtx[i][i] = -curRate * t;
      mtx[i][i+1] = curRate * t;
    }
    
    // return entry 0, size-2 of the matrix exponential
    return MatrixFunctions.expm(new DoubleMatrix(mtx)).get(0, size - 2);
  }

}