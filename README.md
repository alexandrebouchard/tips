<!-- File generated by tutorialj -->

Summary
-------

TIPS is a method for approximate continuous time Markov chain (CTMC) end-point
probability computation and path sampling. 

If you use this library in your work, please cite the following paper:
``Monir Hajiaghayi, Bonnie Kirkpatrick, Liangliang Wang and Alexandre Bouchard-Côté. 
(2014) Efficient Continuous-Time Markov Chain Estimation. International Conference on 
Machine Learning (ICML).``

See: http://www.stat.ubc.ca/~bouchard/pub/icml2014.pdf

TIPS stands for Time Integrated Path Sampling.


Installation
------------

Prerequisite software:

- Java SDK 1.6+
- Gradle version 1.9+ (not tested on Gradle 2.0)

There are several options available to install the package:

### Integrate to a gradle script

Simply add the following lines (replacing 1.0.0 by the current version (see git tags)):

```groovy
repositories {
 mavenCentral()
 jcenter()
 maven {
    url "http://www.stat.ubc.ca/~bouchard/maven/"
  }
}

dependencies {
  compile group: 'ca.ubc.stat', name: 'tips', version: '1.0.0'
}
```

### Compile using the provided gradle script

- Check out the source ``git clone git@github.com:alexandrebouchard/tips.git``
- Compile using ``gradle installApp``
- Add the jars in ``build/install/tips/lib/`` into your classpath

### Use in eclipse

- Check out the source ``git clone git@github.com:alexandrebouchard/tips.git``
- Type ``gradle eclipse`` from the root of the repository
- From eclipse:
  - ``Import`` in ``File`` menu
  - ``Import existing projects into workspace``
  - Select the root
  - Deselect ``Copy projects into workspace`` to avoid having duplicates


Example: PIP transition probabilities
-------------------------------------

Here is a simple example showing how to sample discrete paths and compute 
transition probabilities for the Poisson Indel Process (PIP) 
(http://www.pnas.org/content/early/2012/12/27/1220450110.full.pdf).

For more comprehensive tests, see TestPIP.

```java
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
double estimate = is.estimateTransitionPr(pipMain.getStart(), pipMain.getEnd(), pipMain.bl);
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
```
<sub>From:[tips.Doc](src/test/java//tips/Doc.java)</sub>

Extending to other processes
----------------------------

This is done in three steps:

### Specifying the process

The first step consists in writing a class implementing 
the Process interface. This in turns means implementing 
a single method giving the rates of departure of a given 
state.

Here is a simple example:


(Note: using TIPS for a birth death process is overkill, as
specialized methods exist for numerical calculation of transition
probabilities of birth-death process, see 
FW Crawford and MA Suchard. Transition probabilities for general 
birth-death processes with applications in ecology, genetics, 
and evolution. J Math Biol, 65:553-580, 2012.)

```java
public briefj.collections.Counter rates(java.lang.Integer)
{
    Counter<Integer> result = new Counter<Integer>();
    result.setCount(point + 1, birthRate  * point + 1);
    if (point > 0)
      result.setCount(point - 1, deathRate * point);
    return result;
}
```
<sub>From:[tips.bd.SimpleBirthDeathProcess](src/test/java//tips/bd/SimpleBirthDeathProcess.java)</sub>

### Creating a potential

The second step is to create a potential that will guide
the proposed paths towards the end point.

Again, here is a simple example:


```java
public double get(java.lang.Integer,java.lang.Integer)
{
    // Just return how far we are from the target.
    return Math.abs(proposed - target);
}
```
<sub>From:[tips.bd.SimpleBirthDeathPotential](src/test/java//tips/bd/SimpleBirthDeathPotential.java)</sub>

### Calling TIPS

Here is an example of how to call TIPS with a custom process:


```java
System.out.println("BD example");
SimpleBirthDeathProcess process = new SimpleBirthDeathProcess();
SimpleBirthDeathPotential potential = new SimpleBirthDeathPotential();
TimeIntegratedPathSampler<Integer> sampler = new TimeIntegratedPathSampler<Integer>(potential, process);

double t = 1;
sampler.nParticles = 1000000;
double estimate = sampler.estimateTransitionPr(1, 0, t);
double exact = 0.25; // computed using FW Crawford and MA Suchard, 2012
System.out.println("exact = " + exact);
System.out.println("TIPS = " + estimate);
Assert.assertEquals(0.25, estimate, 1e-3);
System.out.println();
```

Limitations
-----------

- The code currently does not support absorbing state, but this would be easy to fix.
- The method works best when the number of transitions between the end points is not too large.
- The code has not been optimized for speed. For example, the way rates are transmitted
  to the sampler is very wasteful (creating a Counter at each query)
- Computations are not done in log scale nor with scalings, so numerical under/over flow could 
  occur when estimating small transition probabilities.


How it works
------------

The core of this package is in TimeIntegratedPathSampler, which in turns 
is based on two functions:


First, runTIPS(), shown below, which is just an IS algorithm.

The argument ``keepPath`` controls whether sampled paths should be kept or not.
If true, then return a Counter over paths; if false, return 
just the transition pr estimate (average of the weights).
The latter is useful because it runs in constant memory.

```java
private java.lang.Object runTIPS(S,S,double,boolean)
{
    Counter<List<S>> result = keepPath ? new Counter<List<S>>() : null;
    double sum = 0.0;

    for (int particleIndex = 0; particleIndex < nParticles; particleIndex++)
    {
      // propose
      Pair<List<S>, Double> proposed = proposal.propose(rand, x, y, t);

      // compute weight
      double weight = marginalizeSojournTimes(process, proposed.getLeft(), t);
      List<S> proposedJumps = proposed.getLeft();
      for (int jumpIndex = 0; jumpIndex < proposedJumps.size() - 1; jumpIndex++)
        weight *= ProcessUtils.transitionProbability(process, proposedJumps.get(jumpIndex), proposedJumps.get(jumpIndex+1));

      if (unnormalizedWeightsStatistics != null) 
        unnormalizedWeightsStatistics.addValue(weight/proposed.getRight());

      if (keepPath)
        result.incrementCount(proposed.getLeft(), weight/proposed.getRight());
      else
        sum += weight/proposed.getRight();
    }

    return keepPath ? result : sum;
}
```
<sub>From:[tips.TimeIntegratedPathSampler](src/main/java//tips/TimeIntegratedPathSampler.java)</sub>

Second, marginalizeSojournTimes(), shown below, implementing Proposition 2 in
the paper.

The argument ``keepPath`` controls whether sampled paths should be kept or not.
If true, then return a Counter over paths; if false, return 
just the transition pr estimate (average of the weights).
The latter is useful because it runs in constant memory.

```java
public static <S> double marginalizeSojournTimes(tips.Process,java.util.List,double)
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
```
<sub>From:[tips.TimeIntegratedPathSampler](src/main/java//tips/TimeIntegratedPathSampler.java)</sub>

