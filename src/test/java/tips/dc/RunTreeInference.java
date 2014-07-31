package tips.dc;

import java.util.Random;

import org.junit.Test;

import tips.Potential;
import tips.bd.ReversibleBDProcess;
import tips.bd.SimpleBirthDeathPotential;
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
import briefj.run.Mains;
import conifer.TopologyUtils;
import conifer.factors.NonClockTreePrior;



public class RunTreeInference implements Processor, Runnable
{

  private final int nTaxa = 10;
  private Model modelSpec;
  
  public class Model
  {
    private final ReversibleBDProcess<ReversibleBDProcess.ExpectedLengthParameterization> process = ReversibleBDProcess.normalizedIntensityWithExpectedLength(10.0);
    private final Potential<Integer> potential = new SimpleBirthDeathPotential();
    
    @DefineFactor(onObservations = true)
    public final TipsTreeLikelihood<Integer> approximateLikelihood = new TipsTreeLikelihood<Integer>(10, TopologyUtils.syntheticTaxaList(nTaxa))
      .withEvolutionaryProcess(process, potential);
    
    @DefineFactor
    public final NonClockTreePrior<RateParameterization> treePrior = 
      NonClockTreePrior.on(approximateLikelihood.tree).withExponentialRate(0.5);
    
    @DefineFactor
    public final Gamma<ScaleShapeParameterization> processParameterPrior = Gamma.on(process.parameters.expectedLength).withScaleShape(0.5, 20);
  }
  
  public static void main(String [] args)
  {
    Mains.instrumentedRun(args, new RunTreeInference());
  }
  
  public void run()
  {
    modelSpec = new Model();
    
    // forward sample
    Random generationRandom = new Random(1);
    ProbabilityModel generatingModel = new ProbabilityModel(modelSpec);
    ForwardSampler forwardSampler = new ForwardSampler(generatingModel);
    forwardSampler.simulate(generationRandom);
    
    // set param to some wrong value
    System.out.println("True parameters: " + modelSpec.process.parameters);
//    modelSpec.process.parameters.expectedLength.setValue(1.0);
    
//    modelSpec.approximateLikelihood.enableTestMode();
    modelSpec.approximateLikelihood.setNParticles(100);
    MCMCFactory mcmcFactory = new MCMCFactory();
    mcmcFactory.mcmcOptions.thinningPeriod = 1;
    mcmcFactory.addProcessor(this);
    mcmcFactory.mcmcOptions.nMCMCSweeps = 100;
    MCMCAlgorithm mcmc = mcmcFactory.build(modelSpec, false);
    System.out.println(mcmc.model);
    mcmc.run();
  }

  @Override
  public void process(ProcessorContext context)
  {
    // TODO Auto-generated method stub
    
  }
}
