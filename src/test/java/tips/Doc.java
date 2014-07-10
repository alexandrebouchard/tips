package tips;


import java.util.Random;

import muset.MSAPoset;

import org.junit.Assert;
import org.junit.Test;


import tips.bd.SimpleBirthDeathPotential;
import tips.bd.SimpleBirthDeathProcess;
import tips.pip.PIPMain;
import tips.pip.PIPString;
import tips.pip.TestPIP;
import tips.utils.Baselines;
import tutorialj.Tutorial;




public class Doc
{
  /**
   * 
   * Summary
   * -------
   * 
   * TIPS is an approximate method to compute continuous time Markov chain (CTMC) end-point
   * probabilities. 
   * 
   * Please cite the following paper:
   * ``Monir Hajiaghayi, Bonnie Kirkpatrick, Liangliang Wang and Alexandre Bouchard-Côté. 
   * (2014) Efficient Continuous-Time Markov Chain Estimation. International Conference on 
   * Machine Learning (ICML).``
   * 
   * See: http://www.stat.ubc.ca/~bouchard/pub/icml2014.pdf
   * 
   * TIPS stands for Time Integrated Path Sampling.
   * 
   * 
   * Installation
   * ------------
   * 
   * There are several options available to install the package:
   * 
   * ### Integrate to a gradle script
   * 
   * Simply add the following lines (replacing 1.0.0 by the current version (see git tags)):
   * 
   * ```groovy
   * repositories {
   *  mavenCentral()
   *  jcenter()
   *  maven {
   *     url "http://www.stat.ubc.ca/~bouchard/maven/"
   *   }
   * }
   * 
   * dependencies {
   *   compile group: 'ca.ubc.stat', name: 'tips', version: '1.0.0'
   * }
   * ```
   * 
   * ### Compile using the provided gradle script
   * 
   * - Check out the source ``git clone git@github.com:alexandrebouchard/tips.git``
   * - Compile using ``gradle installApp``
   * - Add the jars in ``build/install/tips/lib/`` into your classpath
   * 
   * ### Use in eclipse
   * 
   * - Check out the source ``git clone git@github.com:alexandrebouchard/tips.git``
   * - Type ``gradle eclipse`` from the root of the repository
   * - From eclipse:
   *   - ``Import`` in ``File`` menu
   *   - ``Import existing projects into workspace``
   *   - Select the root
   *   - Deselect ``Copy projects into workspace`` to avoid having duplicates
   */
  @Tutorial(startTutorial = "README.md", showSource = false)
  public void installInstructions() {}
  
  /**
   * Example: PIP transition probabilities
   * -------------------------------------
   * 
   * Here is a simple example showing how to sample discrete paths and compute 
   * transition probabilities for the Poisson Indel Process (PIP) 
   * (http://www.pnas.org/content/early/2012/12/27/1220450110.full.pdf).
   * 
   * For more comprehensive tests, see TestPIP.
   */
  @Tutorial(showLink = true, linkPrefix = "src/test/java/")
  @Test
  public void pipExample()
  {
    System.out.println("PIP example");
    PIPMain pipMain = new PIPMain();
    
    // set some process parameters
    pipMain.bl = 0.5;
    pipMain.mu = 0.5;
    pipMain.lambda = 5;
    pipMain.ensureLinearizationUnique = true;
    
    // generate data
    MSAPoset align = pipMain.getGeneratedEndPoints();
    
    // compare to the exact transition probability
    double exact = Math.exp(TestPIP.exact(pipMain.mu, pipMain.lambda, pipMain.bl, align));
    System.out.println("exact = " + exact);
    
    // create a TIPS sampler
    TimeIntegratedPathSampler<PIPString> is = pipMain.buildImportanceSampler();
    is.nParticles = 10000;
    is.rand = new Random(1);
    pipMain.potentialProposalOptions.automatic = true;
      
    // sample
    double estimate = is.estimateZ(pipMain.getStart(), pipMain.getEnd(), pipMain.bl);
    System.out.println("TIPS = " + estimate);
    Assert.assertEquals(exact, estimate, 1e-6);
      
    // compare to some alternate methods
    System.out.println("approximateExhaustiveSum = " + 
        Baselines.exhaustiveSum(is.rand, is.nParticles, pipMain.getProcess(), pipMain.getProposal(), 
            pipMain.getStart(), pipMain.getEnd(), pipMain.bl));
    System.out.println("naiveIS = " + 
        Baselines.standardIS(pipMain.getGeneratedEndPoints(), pipMain.bl, pipMain.getProcess(), 
            is.nParticles, is.rand, null));
    System.out.println();
  }
  
  /**
   * Extending to other processes
   * ----------------------------
   * 
   * This is done in three steps:
   * 
   * ### Specifying the process
   * 
   * The first step consists in writing a class implementing 
   * the Process interface. This in turns means implementing 
   * a single method giving the rates of departure of a given 
   * state.
   * 
   * Here is a simple example:
   */
  @Tutorial(showSource = false, nextStep = SimpleBirthDeathProcess.class)
  public void extending1() {}
  
  /**
   * ### Creating a potential
   * 
   * The second step is to create a potential that will guide
   * the proposed paths towards the end point.
   * 
   * Again, here is a simple example:
   */
  @Tutorial(showSource = false, nextStep = SimpleBirthDeathPotential.class)
  public void extending2() {}
  
  /**
   * ### Calling TIPS
   * 
   * Here is an example of how to call TIPS with a custom process:
   * 
   */
  @Tutorial(showSource = true)
  @Test
  public void extending3() 
  {
    System.out.println("BD example");
    SimpleBirthDeathProcess process = new SimpleBirthDeathProcess();
    SimpleBirthDeathPotential potential = new SimpleBirthDeathPotential();
    TimeIntegratedPathSampler<Integer> sampler = new TimeIntegratedPathSampler<Integer>(potential, process);
    
    double t = 1;
    sampler.nParticles = 1000000;
    double estimate = sampler.estimateZ(1, 0, t);
    double exact = 0.25; // computed using FW Crawford and MA Suchard, 2012
    System.out.println("exact = " + exact);
    System.out.println("TIPS = " + estimate);
    Assert.assertEquals(0.25, estimate, 1e-3);
    System.out.println();
  }
  
  /**
   * Limitations
   * -----------
   * 
   * - The code currently does not support absorbing state, but this would be easy to fix.
   * - The method works best when the number of transitions between the end points is not too large.
   * - The code has not been optimized for speed. For example, the way rates are transmitted
   *   to the sampler is very wasteful (creating a Counter at each query)
   */
  @Tutorial(showSource = false)
  public void limitations() {}
}
