package tips.dc;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;


import tips.ForwardSampler;
import tips.Potential;
import tips.Proposal;
import tips.StationaryProcess;
import tips.TimeIntegratedPathSampler;
import tips.bd.ReversibleBDProcess.InvalidParametersException;
import tips.utils.PotPropOptions;
import tips.utils.PotProposal;
import bayonet.distributions.Multinomial;
import bayonet.math.NumericalUtils;
import blang.annotations.FactorArgument;
import blang.annotations.FactorComponent;
import blang.factors.GenerativeFactor;
import briefj.collections.Counter;
import briefj.collections.Tree;
import conifer.TreeNode;
import conifer.UnrootedTree;
import conifer.UnrootedTreeUtils;
import conifer.factors.UnrootedTreeLikelihood;


  
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
  public final StationaryProcess<S> evolutionaryProcess;
 
  /**
   * 
   */
  @FactorArgument(makeStochastic = true)
  public final TipsTreeObservation<S> observations;
  
  /**
   * 
   */
  private int nParticles = 1000;
  
  /**
   * 
   */
  private final Proposal<S> proposal;
  
  /**
   * Source of randomness used by the algorithm.
   * 
   * By default, the seed is fixed to 1.
   */
  public Random rand = new Random(1);
  
  public void setNParticles(int nParticles)
  {
    this.nParticles = nParticles;
  }
  
  private TimeIntegratedPathSampler<S> getSampler()
  {
    TimeIntegratedPathSampler<S> result = new TimeIntegratedPathSampler<S>(proposal, evolutionaryProcess);
    result.nParticles = nParticles;
    result.rand = rand;
    return result;
  }
  
  public TipsTreeLikelihood(UnrootedTree tree, StationaryProcess<S> evolutionaryProcess,
      TipsTreeObservation<S> observations, Proposal<S> proposal)
  {
    this.tree = tree;
    this.evolutionaryProcess = evolutionaryProcess;
    this.observations = observations;
    this.cached = null;
    this.proposal = proposal;
  }
  
  public TipsTreeLikelihood(int nSites, List<TreeNode> leaves)
  {
    this(UnrootedTreeLikelihood.defaultTree(leaves), null, new TipsTreeObservation<S>(nSites), null);
  }
  
  public static <S> TipsTreeLikelihood<S> fromObservations(TipsTreeObservation<S> observations)
  {
    return new TipsTreeLikelihood<S>(null, null, observations, null);
  }
  
  public TipsTreeLikelihood<S> withEvolutionaryProcess(StationaryProcess<S> evolutionaryProcess, Potential<S> potential)
  {
    return new TipsTreeLikelihood<S>(this.tree, evolutionaryProcess, this.observations, new PotProposal<S>(evolutionaryProcess, potential, new PotPropOptions()));
  }
  
  public TipsTreeLikelihood<S> withTree(UnrootedTree tree)
  {
    return new TipsTreeLikelihood<S>(tree, this.evolutionaryProcess, this.observations, this.proposal);
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
      this.tree = new UnrootedTree(currentConfig.tree);
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
  
  private boolean testMode = false;
  
  void enableTestMode() 
  {
    testMode = true; 
  }

  @Override
  public double logDensity()
  {
    if (testMode)
      return testLogDensity();
    if (cached == null || !cached.valid())
      cached = recomputeCache();
    return cached.value;
  }

  private double testLogDensity()
  {
    this.rand = new Random(1);
    CachedDensity before = cached;
    CachedDensity after = recomputeCache();
    if (cached == null)
      cached = after;
    else
    {
      if (cached.valid())
        NumericalUtils.checkIsClose(before.value, after.value);
      else
      {
        if (NumericalUtils.isClose(before.value, after.value, NumericalUtils.THRESHOLD))
          throw new RuntimeException();
        cached = after;
      }
    }
    return cached.value;
  }

  private CachedDensity recomputeCache()
  {
    double logDensity = recomputeLogDensity(); //recomputeLogDensity_dcSMC();
    return new CachedDensity(this, logDensity);
  }
  
  public double recomputeLogDensity()
  {
    try 
    {
           if (likelihoodCalculationMethod == LikelihoodCalculationMethod.DC)
        return recomputeLogDensity_dcSMC();
      else if (likelihoodCalculationMethod == LikelihoodCalculationMethod.SMC_BOTTOM_UP)
        return recomputeLogDensity_standardSMC(true);
      else if (likelihoodCalculationMethod == LikelihoodCalculationMethod.SMC_DFS)
        return recomputeLogDensity_standardSMC(false);
      else
        throw new RuntimeException();
    }
    catch (InvalidParametersException e)
    {
      return Double.NEGATIVE_INFINITY;
    }
  }

  public static enum LikelihoodCalculationMethod { DC, SMC_BOTTOM_UP, SMC_DFS }
  public LikelihoodCalculationMethod likelihoodCalculationMethod = LikelihoodCalculationMethod.DC;

  private double recomputeLogDensity_dcSMC()
  {
    Tree<Pair<TreeNode, Double>> topologicalCentroidRooting = TipsTreeUtils.topologicalCentroidRooting(tree);
    
    double result = 0.0;
    
    for (int site = 0; site < observations.nSites(); site++)
    {
      TreeNodeSample<S> rootPopulation = divideAndConquer(site, topologicalCentroidRooting);
      result += rootPopulation.logNormalization;
    }
    
    return result;
  }
  
  private double recomputeLogDensity_standardSMC(boolean useBottomUpTraversal)
  {
    Tree<Pair<TreeNode, Double>> topologicalCentroidRooting = TipsTreeUtils.topologicalCentroidRooting(tree);
    
    double result = 0.0;
    
    for (int site = 0; site < observations.nSites(); site++)
    {
      standardSMCSamples = initStandardSMCSamples(site, topologicalCentroidRooting);
      standardSMCPrs = initStandardSMCPrs();
      standardSMCLogNorm = 0.0;
      standardSMC(site, topologicalCentroidRooting, useBottomUpTraversal);
      result += standardSMCLogNorm;
    }
    
    return result;
  }
  
  private double[] initStandardSMCPrs()
  {
    double [] result = new double[1];
    result[0] = 1.0;
    return result;
  }

  private List<Map<TreeNode, S>> initStandardSMCSamples(int site, Tree<Pair<TreeNode, Double>> topology)
  {
    Map<TreeNode,S> initialParticle = Maps.newLinkedHashMap();
    
    for (Tree<Pair<TreeNode, Double>> desc : topology.getPreOrderTraversal())
      if (desc.isLeaf())
      {
        TreeNode node = desc.getLabel().getLeft();
        initialParticle.put(node, observations.get(site, node));
      }
    
    return Collections.singletonList(initialParticle);
  }

  private List<Map<TreeNode,S>> standardSMCSamples = null;
  private double [] standardSMCPrs = null;
  private double standardSMCLogNorm = Double.NaN;
  private void standardSMC(int siteIndex, Tree<Pair<TreeNode, Double>> root, boolean useBottomUpTraversal)
  {
    List<Tree<Pair<TreeNode, Double>>> traversalOrder = useBottomUpTraversal ? bottomUpTraversal(root) : root.getPostOrderTraversal();
    nodeLoop : for (Tree<Pair<TreeNode, Double>> node : traversalOrder)
    {
      if (node.isLeaf())
        continue nodeLoop;

      TimeIntegratedPathSampler<S> sampler = getSampler();
      final int nParticles = sampler.nParticles;
      List<Map<TreeNode,S>> newSamples = Lists.newArrayList();
      double [] newWeights = new double[nParticles];
      
      final Tree<Pair<TreeNode,Double>>
        child0 = node.getChildren().get(0),
        child1 = node.getChildren().get(1);
      
      final double 
        branchLength0 = child0.getLabel().getRight(),
        branchLength1 = child1.getLabel().getRight();
      
      final TreeNode
        currentNode = node.getLabel().getLeft(),
        treeNode0 = child0.getLabel().getLeft(),
        treeNode1 = child1.getLabel().getLeft();
      
      for (int i = 0; i < nParticles; i++)
      {
        final Map<TreeNode,S> sample = standardSMCSamples.get(Multinomial.sampleMultinomial(sampler.rand, standardSMCPrs));
        
        final S 
          childSample0 = sample.get(treeNode0),
          childSample1 = sample.get(treeNode1);
        
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
        
        Map<TreeNode,S> newSample = Maps.newLinkedHashMap(sample);
        newSample.put(currentNode, sampleTreeCherry.getLeft());
        newSamples.add(newSample);
        newWeights[i] = weight;
      }
      
      double normalization = Multinomial.normalize(newWeights);
      
      standardSMCSamples = newSamples;
      standardSMCPrs = newWeights;
      standardSMCLogNorm += Math.log(normalization / nParticles);
    }
  }

  private List<Tree<Pair<TreeNode, Double>>> bottomUpTraversal(
      Tree<Pair<TreeNode, Double>> node)
  {
    List<List<Tree<Pair<TreeNode, Double>>>> recursions = Lists.newArrayList();
    
    for (Tree<Pair<TreeNode, Double>> child : node.getChildren())
      recursions.add(bottomUpTraversal(child));
    
    int nItems = 0;
    for (List<Tree<Pair<TreeNode, Double>>> recursion : recursions)
      if (recursion.size() > nItems)
        nItems = recursion.size();
    
    List<Tree<Pair<TreeNode, Double>>> result = Lists.newArrayList();
    for (int i = 0; i < nItems; i++)
      for (List<Tree<Pair<TreeNode, Double>>> recursion : recursions)
        if (i < recursion.size())
          result.add(recursion.get(i));
    
    result.add(node);
    
    return result;
  }

  private TreeNodeSample<S> divideAndConquer(int siteIndex, Tree<Pair<TreeNode, Double>> node)
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
        childPop0 = divideAndConquer(siteIndex, child0),
        childPop1 = divideAndConquer(siteIndex, child1);
      
      final double 
        branchLength0 = child0.getLabel().getRight(),
        branchLength1 = child1.getLabel().getRight();
      
      TimeIntegratedPathSampler<S> sampler = getSampler();
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
    Tree<Pair<TreeNode, Double>> aRooting = TipsTreeUtils.topologicalCentroidRooting(tree);
    
    for (int site = 0; site < observations.nSites(); site++)
    {
      // simulate root
      S rootState = evolutionaryProcess.sampleFromStationary(random);
      
      for (Tree<Pair<TreeNode, Double>> child : aRooting.getChildren())
        generate(random, rootState, child, site);
    }
    
  }

  private void generate(Random random, S parentState,
      Tree<Pair<TreeNode, Double>> node, int site)
  {
    // simulate
    final double branchLength = node.getLabel().getRight();
    ForwardSampler<S> forwardSampler = new ForwardSampler<S>(evolutionaryProcess);
    S currentState = forwardSampler.forwardSample(random, parentState, branchLength);
    
    if (node.isLeaf())
      observations.set(site, node.getLabel().getLeft(), currentState);
    
    for (Tree<Pair<TreeNode, Double>> child : node.getChildren())
      generate(random, currentState, child, site);
  }
}
