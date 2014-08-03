package tips.dc;

import java.util.Map;
import java.util.Random;

import tips.Potential;
import tips.bd.ReversibleBDProcess;
import tips.bd.SimpleBirthDeathPotential;
import tips.bd.SimpleBirthDeathProcess;
import bayonet.distributions.Exponential.RateParameterization;
import bayonet.distributions.Gamma;
import bayonet.distributions.Gamma.ScaleShapeParameterization;
import blang.ForwardSampler;
import blang.MCMCAlgorithm;
import blang.MCMCFactory;
import blang.ProbabilityModel;
import blang.annotations.DefineFactor;
import blang.processing.Processor;
import blang.processing.ProcessorContext;
import briefj.collections.UnorderedPair;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import conifer.TopologyUtils;
import conifer.TreeNode;
import conifer.UnrootedTree;
import conifer.UnrootedTreeUtils;
import conifer.UnrootedTreeUtils.TreeMetric;
import conifer.factors.NonClockTreePrior;
import conifer.processors.TreeDistanceProcessor;



public class RunTreeInference implements Processor, Runnable
{
  @Option public int nTaxa = 5;
  @Option public int nSites = 100;
  @Option public int nParticles = 10;
  @Option public Random generationRandom = new Random(1);
  @OptionSet(name = "mcmc") public final MCMCFactory mcmcFactory = new MCMCFactory();
  
  public class Model
  {
//    private final SimpleBirthDeathPotential potential = new SimpleBirthDeathPotential();
//    private final SimpleBirthDeathProcess process = new SimpleBirthDeathProcess();
    private final ReversibleBDProcess<ReversibleBDProcess.ExpectedLengthParameterization> process = ReversibleBDProcess.normalizedIntensityWithExpectedLength(10.0);
    private final Potential<Integer> potential = new SimpleBirthDeathPotential();
    
    @DefineFactor(onObservations = true)
    public final TipsTreeLikelihood<Integer> approximateLikelihood = new TipsTreeLikelihood<Integer>(nSites, TopologyUtils.syntheticTaxaList(nTaxa))
      .withEvolutionaryProcess(process, potential);
    
    @DefineFactor
    public final NonClockTreePrior<RateParameterization> treePrior = 
      NonClockTreePrior.on(approximateLikelihood.tree).withExponentialRate(5.0);
    
    @DefineFactor
    public final Gamma<ScaleShapeParameterization> processParameterPrior = Gamma.on(process.parameters.expectedLength).withScaleShape(0.5, 5);
  }
  
  public static void main(String [] args)
  {
    Mains.instrumentedRun(args, new RunTreeInference());
  }
  
  private Model modelSpec;
  private UnrootedTree trueTree;
  public void run()
  {
    modelSpec = new Model();
    
    // forward sample
    ProbabilityModel generatingModel = new ProbabilityModel(modelSpec);
    ForwardSampler forwardSampler = new ForwardSampler(generatingModel);
    forwardSampler.simulate(generationRandom);
    
    trueTree = new UnrootedTree(modelSpec.approximateLikelihood.tree);
    
    // make the true tree reasonably well behaved
    for (UnorderedPair<TreeNode, TreeNode> edge : modelSpec.treePrior.tree.getBranchLengths().keySet())
      modelSpec.treePrior.tree.updateBranchLength(edge, 0.1);
    System.out.println("True tree: " + trueTree);
    
    // set tree to wrong value
    modelSpec.treePrior.generate(generationRandom);
    for (UnorderedPair<TreeNode, TreeNode> edge : modelSpec.treePrior.tree.getBranchLengths().keySet())
      modelSpec.treePrior.tree.updateBranchLength(edge, 0.05);
    
    System.out.println("Init tree: " + modelSpec.approximateLikelihood.tree);
    Map<TreeMetric, Double> currentTreeMetrics = UnrootedTreeUtils.computeTreeMetrics(trueTree, modelSpec.approximateLikelihood.tree);
    System.out.println(currentTreeMetrics);
    
    
    // set param to some wrong value
//    System.out.println("True parameters: " + modelSpec.process.parameters);
//    modelSpec.process.parameters.expectedLength.setValue(1.0);
    
//    modelSpec.approximateLikelihood.enableTestMode();
    modelSpec.approximateLikelihood.setNParticles(nParticles);
    mcmcFactory.addProcessor(this);
    
    modelSpec.approximateLikelihood.rand = mcmcFactory.mcmcOptions.random;
    MCMCAlgorithm mcmc = mcmcFactory.build(modelSpec, false);
    TreeDistanceProcessor.setReferenceTree(mcmc, modelSpec.treePrior.tree, trueTree);
    System.out.println(mcmc.model);
    mcmc.run();
  }
  

  @Override
  public void process(ProcessorContext context)
  {
    // compute current distance to true tree
//    Map<TreeMetric, Double> currentTreeMetrics = UnrootedTreeUtils.computeTreeMetrics(trueTree, modelSpec.approximateLikelihood.tree);
//    System.out.println(currentTreeMetrics);
//    System.out.println(modelSpec.process.parameters.expectedLength);
  }
}
