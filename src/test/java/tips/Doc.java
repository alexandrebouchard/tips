package tips;


import java.util.Random;

import muset.MSAPoset;

import org.junit.Test;

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
   * See http://www.stat.ubc.ca/~bouchard/pub/icml2014.pdf
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
    PIPMain pipMain = new PIPMain();
    
    // set some process parameters
    pipMain.bl = 0.5;
    pipMain.mu = 0.5;
    pipMain.lambda = 5;
    pipMain.ensureLinearizationUnique = true;
    
    // generate data
    MSAPoset align = pipMain.getGeneratedEndPoints();
    System.out.println(align);
    
    if (pipMain.ensureLinearizationUnique)
    {
      // compare to the exact transition probability
      double exact = Math.exp(TestPIP.exact(pipMain.mu, pipMain.lambda, pipMain.bl, align));
      System.out.println("exact = " + exact);
    }
    
    // create a sampler
    TimeIntegratedPathSampler<PIPString> is = pipMain.buildImportanceSampler();
    is.nParticles = 10000;
    is.rand = new Random(1);
    pipMain.potentialProposalOptions.automatic = true;
      
    // sample
    double estimate = is.estimateZ(pipMain.getStart(), pipMain.getEnd(), pipMain.bl);
    System.out.println("TIPS estimate = " + estimate);
      
    // compare to some alternate methods
    System.out.println("approximate exhaustive sum = " + Baselines.exhaustiveSum(is.rand, is.nParticles, pipMain.getProcess(), pipMain.getProposal(), pipMain.getStart(), pipMain.getEnd(), pipMain.bl));
    System.out.println("naive IS = " + Baselines.standardIS(pipMain.getGeneratedEndPoints(), pipMain.bl, pipMain.getProcess(), is.nParticles, is.rand, null));
  }
}
