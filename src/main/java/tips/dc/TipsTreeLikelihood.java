package tips.dc;

import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

import com.google.common.collect.Lists;

import tips.Process;
import tips.TimeIntegratedPathSampler;
import bayonet.distributions.Multinomial;
import blang.annotations.FactorArgument;
import blang.annotations.FactorComponent;
import blang.factors.GenerativeFactor;
import briefj.BriefMath;
import briefj.collections.Counter;
import briefj.collections.Tree;
import conifer.TreeNode;
import conifer.UnrootedTree;
import conifer.UnrootedTreeUtils;


/**
 * - TODO
 *  - finish inference code
 *  - implement finite CTMC-based code
 *  - implement generation code
 *  - 
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 * @param <S>
 */
public class TipsTreeLikelihood<S> implements GenerativeFactor
{
  /**
   * 
   */
  @FactorArgument
  public final UnrootedTree tree;
  
  /**
   * 
   */
  @FactorComponent
  public final Process<S> evolutionaryProcess;
 
  /**
   * 
   */
  @FactorArgument(makeStochastic = true)
  public final TipsTreeObservation<S> observations;
  
  public final TimeIntegratedPathSampler<S> sampler;
  
  public TipsTreeLikelihood(UnrootedTree tree, Process<S> evolutionaryProcess,
      TipsTreeObservation<S> observations, TimeIntegratedPathSampler<S> sampler)
  {
    this.tree = tree;
    this.evolutionaryProcess = evolutionaryProcess;
    this.observations = observations;
    this.cached = null;
    this.sampler = sampler;
  }

  private class CachedDensity
  {
    private final double value;
    private final UnrootedTree tree;
    
    /**
     * Hack: assume that the process has a good hash function to detect changes in parameters values.
     */
    private final long processHash;
    
    /**
     * Hack: assume the observations have a good hash as well (or never change).
     */
    private final long observationsHash;
    
    private CachedDensity(TipsTreeLikelihood<S> currentConfig, double logDensity)
    {
      this.value = logDensity;
      this.tree = currentConfig.tree;
      this.processHash = currentConfig.evolutionaryProcess.hashCode();
      this.observationsHash = currentConfig.observations.hashCode();
    }

    private boolean valid()
    {
      if (UnrootedTreeUtils.computeTreeMetrics(this.tree, TipsTreeLikelihood.this.tree).get(UnrootedTreeUtils.TreeMetric.l1) != 0.0)
        return false;
      if (this.processHash != TipsTreeLikelihood.this.evolutionaryProcess.hashCode())
        return false;
      if (this.observationsHash != TipsTreeLikelihood.this.observations.hashCode())
        return false;
      return true;
    }
  }
  private CachedDensity cached = null;

  @Override
  public double logDensity()
  {
    if (cached == null || !cached.valid())
      refreshCache();
    return cached.value;
  }


  private void refreshCache()
  {
    double logDensity = recomputeLogDensity();
    cached = new CachedDensity(this, logDensity);
  }


  public double recomputeLogDensity()
  {
    Tree<Pair<TreeNode, Double>> topologicalCentroidRooting = TipsTreeUtils.topologicalCentroidRooting(tree);
    
    double result = 0.0;
    
    for (int site = 0; site < observations.nSites(); site++)
      result += sample(site, topologicalCentroidRooting).logNormalization;
    
    return result;
  }

  private static class TreeNodeSample<S>
  {
    private final List<S> samples = Lists.newArrayList();
    private final double [] prs;
    private final double logNormalization;
    private TreeNodeSample(Counter<S> samples, double logNormalization)
    {
      BriefMath.checkCloseAndFinite(samples.totalCount(), 1.0);
      this.prs = new double[samples.size()];
      int i = 0;
      for (S key : samples.keySet())
      {
        prs[i++] = samples.getCount(key);
        this.samples.add(key);
      }
      this.logNormalization = logNormalization;
    }
    public S sample(Random rand)
    {
      return samples.get(Multinomial.sampleMultinomial(rand, prs));
    }
  }

  private TreeNodeSample<S> sample(int siteIndex, Tree<Pair<TreeNode, Double>> node)
  {
    final Counter<S> samples = new Counter<S>(); 
    if (node.getChildren().size() == 0)
    {
      samples.incrementCount(observations.get(siteIndex, node.getLabel().getLeft()), 1.0);
      return new TreeNodeSample<S>(samples, 0.0);
    }
    else
    {
      if (node.getChildren().size() != 2)
        throw new RuntimeException("Only bifurcating trees supported (arity=" + node.getChildren().size() + ")");
      
      final Tree<Pair<TreeNode,Double>>
        child0 = node.getChildren().get(0),
        child1 = node.getChildren().get(1);
      
      final TreeNodeSample<S>
        childPop0 = sample(siteIndex, child0),
        childPop1 = sample(siteIndex, child1);
      
      final double 
        branchLength0 = child0.getLabel().getRight(),
        branchLength1 = child1.getLabel().getRight();
      
      final int nParticles = sampler.nParticles;
      for (int particleIndex = 0; particleIndex < nParticles; particleIndex++)
      {
        final S 
          childSample0 = childPop0.sample(sampler.rand),
          childSample1 = childPop1.sample(sampler.rand);
        
        final Pair<S, Double> sampleTreeCherry = sampler.sampleTreeCherry(
            childSample0, 
            childSample1, 
            branchLength0, 
            branchLength1);
        
        double weight = sampleTreeCherry.getRight();
        
        if (!child0.isLeaf())
          weight = weight / sampler.getStationaryPr(childSample0);
        
        if (!child1.isLeaf())
          weight = weight / sampler.getStationaryPr(childSample1);
        
        samples.incrementCount(sampleTreeCherry.getLeft(), weight);
      }
    
      final double logNormalization = 
        childPop0.logNormalization + 
        childPop1.logNormalization +
        Math.log(samples.totalCount()/((double) nParticles));
      
      samples.normalize();
      return new TreeNodeSample<S>(samples, logNormalization);
    }
  }

  @Override
  public void generate(Random random)
  {
    throw new RuntimeException();
  }
}
