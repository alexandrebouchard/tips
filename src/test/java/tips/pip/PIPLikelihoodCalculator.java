package tips.pip;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;
import org.jgrapht.Graphs;
import org.jgrapht.UndirectedGraph;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import bayonet.graphs.GraphUtils;
import bayonet.math.NumericalUtils;
import bayonet.math.SpecialFunctions;
import briefj.collections.UnorderedPair;

import muset.LinearizedAlignment;
import muset.SequenceId;



/**
 * For some of the notation used in the inner workings of this class, see
 * The Poisson Indel Process, arXiv:1207.6327 
 * 
 * @author bouchard
 */
public class PIPLikelihoodCalculator
{
  
  public PIPLikelihoodCalculator(
      PoissonParameters pip,
      LinearizedAlignment linearizedMSA, 
      UndirectedGraph<PIPTreeNode, UnorderedPair<PIPTreeNode, PIPTreeNode>> topology,
      Map<UnorderedPair<PIPTreeNode, PIPTreeNode>,Double> branchLengths,
      PIPTreeNode root)
  {
    this.pip = pip;
    this.linearizedMSA = linearizedMSA;
    this.root = root;
    this.branchLengths = convert(branchLengths, topology, root);
    // derived fields
    this.postOrderTaxaTraversal = GraphUtils.postorder(topology, root); //new ArrayList<PIPTreeNode>();
//    PIPLikelihoodUtils.fillPostOrder(postOrderTaxaTraversal, tree.topology());
    this.nMSAColumns = linearizedMSA.nColumns();
    this.mu = pip.deleteRate;
    this.logMu = Math.log(mu);
    this.totalTreeLength = PIPLikelihoodUtils.totalTreeLength(branchLengths);
    this.nu = pip.insertRate * (totalTreeLength + 1.0/ mu);
    this.logNu = Math.log(nu);
    this.nSitesPlusFullGap = nMSAColumns + 1;
    this.fullGapIndex = nSitesPlusFullGap - 1;
    this.nCharacters = pip.numberOfCharacter;
    this.logPi = pip.quasiStatLogProbabilities;
    assert logPi.length == nCharacters;
    this.leaves = Sets.newLinkedHashSet(GraphUtils.leaves(topology)); //new HashSet<PIPTreeNode>(tree.topology().leaveContents());
    initLeftRight(topology, root);
  }
   
  private void initLeftRight(UndirectedGraph<PIPTreeNode, UnorderedPair<PIPTreeNode, PIPTreeNode>> topology, PIPTreeNode root2)
  {
    this.lefts = Maps.newLinkedHashMap();
    this.rights = Maps.newLinkedHashMap();
    Map<PIPTreeNode,PIPTreeNode> parentPtrs = GraphUtils.parentPointers(topology, root);
    for (PIPTreeNode internal : GraphUtils.internalNodes(topology))
    {
      List<PIPTreeNode> neighborListOf = Lists.newArrayList(Graphs.neighborListOf(topology, internal));
      if (internal.equals(root))
      {
        if (neighborListOf.size() != 2)
          throw new RuntimeException();
        
      }
      else
      {
        if (neighborListOf.size() != 3)
          throw new RuntimeException();
        neighborListOf.remove(parentPtrs.get(internal));
      }
      lefts.put(internal, neighborListOf.get(0));
      rights.put(internal, neighborListOf.get(1));
    }
  }

  private static Map<PIPTreeNode, Double> convert(
      Map<UnorderedPair<PIPTreeNode, PIPTreeNode>, Double> branchLengths2, UndirectedGraph<PIPTreeNode, UnorderedPair<PIPTreeNode, PIPTreeNode>> topology, PIPTreeNode root)
  {
    Map<PIPTreeNode,PIPTreeNode> parentPtrs = GraphUtils.parentPointers(topology, root);
    Map<PIPTreeNode,Double> result = Maps.newLinkedHashMap();
    for (UnorderedPair<PIPTreeNode, PIPTreeNode> key : branchLengths2.keySet())
    {
      double value = branchLengths2.get(key);
           if (parentPtrs.get(key.getFirst ()).equals(key.getSecond())) result.put(key.getFirst(), value);
      else if (parentPtrs.get(key.getSecond()).equals(key.getFirst()))  result.put(key.getSecond(),value);
      else throw new RuntimeException();
    }
    return result;
  }

  /**
   * If T is a tree, Y is the observed event, M is an hypothesized MSA event, this computes
   * log P(M, Y | T)
   * (note that by the way we define M, which includes the info of Y, this is the same as
   * log P(M | T)
   */
  public double computeDataLogProbabilityGivenTree()
  {
    computePeelingRecursions();
    prepareLogSurvivalPrs();
    prepareCommonAncestorIndicators();
    prepareLogInsertLocationPrs();
    double [] logZs = computeLogZ_phase1();
    return computeLogZ_phase2_nonGap(logZs) + logPhi(computeZ_phase2_gap(logZs));
  }
  
  // all private below
  
  // main vars
  private PoissonParameters pip;
  private LinearizedAlignment linearizedMSA;
  
  private final Map<PIPTreeNode, Double> branchLengths;
  private final PIPTreeNode root;
  
  // derived
  private List<PIPTreeNode> postOrderTaxaTraversal;
  private int nMSAColumns;
  
  /**
   * Deletion rate (per character)
   */
  private final double logMu, mu;
  
  /**
   * nu = lambda * (totalTreeLength + 1.0/mu)
   * where lambda is the insertion rate
   * 
   * This is the normalization of the intensity measure
   * in the Poisson characterization in arXiv:1207.6327
   */
  private final double nu, logNu;
  
  /**
   * Sum of the branch lengths of the tree
   */
  private final double totalTreeLength;
  
  /**
   * The number of character
   * e.g. this is equal to 4 if DNA is under study
   */
  private final int nCharacters;
  
  /**
   * The log of the quasi stationary distribution over the characters
   * Since the deletion rate is assumed constant, this is just the stat. dist
   * of the subsitution rate matrix
   * 
   * This should have the same dimensionality as nCharacters
   */
  private final double [] logPi;
  
  /**
   * When computing recursions, one pseudo-site is added at the end of the site indexed
   * recursions, corresponding to an empty MSA column
   * This is the index of this pseudo-site
   */
  private final int fullGapIndex;
  
  /**
   * = nColumns + 1 (see fullGapIndex)
   */
  private final int nSitesPlusFullGap;
  private final Set<PIPTreeNode> leaves;
  
  /**
   * We assume the rooted tree is binary, so this
   * returns the left and right edges
   * 
   * This could be relaxed without conceptual difficulties
   */
  private  Map<PIPTreeNode,PIPTreeNode> lefts, rights;
  
  // computational steps
  
  /**
   * For an internal nodes v, logSurvivalPrs.get(v) is a 
   * probability related to the edge e going from the parent
   * of v to v. It is the probability that the character 
   * survives at least to v given it was inserted at some uniform
   * location on e.
   * 
   * For v equal to the root, we set the value to log(1) = 0 for convenience
   * 
   * Denoted \beta in arXiv:1207.6327
   */
  private Map<PIPTreeNode, Double> logSurvivalPrs;
  
  /**
   * For each node and site, is this a common ancestor of all
   * the observed characters at that site?
   * Denoted 1[v \in A] in arXiv:1207.6327
   */
  private Map<PIPTreeNode, boolean[]> commonAncestorIndicators;
  
  /**
   * Arrays used to perform a standard peeling recursion
   * Given a taxon v, site s and character c,
   * logPeelingDetails.get(v)[s][c] returns the conditional probability of the subtree
   * rooted at v, the data under v, and for a single site s given that we have the character c at the root
   * 
   * Denoted \tilde f(\sigma) in arXiv:1207.6327
   */
  private Map<PIPTreeNode, double[][]> logPeelingDetails;
  
  /**
   * Obtained from logPeelingDetails by summing over the values at the root
   * under the quasi-stationary prior
   * 
   * Denoted \tilde f in arXiv:1207.6327
   */
  private Map<PIPTreeNode, double[]> logPeelingSummaries;
  
  /**
   * Given that there is an insertion, the log probabilities of where this insertion is
   * for taxon v not equal to the root, the corresponding edge is the one between v and its parent
   * For the root, the value is derived from the stationary le distribution (see arXiv:1207.6327)
   * 
   * Denoted \iota in arXiv:1207.6327
   */
  private Map<PIPTreeNode, Double> logInsertLocationPrs;
  
  /**
   * descCounts.get(v)[s] track of the number of non-gap symbol under each subtree rooted by v at
   * site s
   * 
   * An intermediate quantity used in the computation of commonAncestorIndicators
   */
  private Map<PIPTreeNode, int[]> descCounts;
  
  /**
   * nonGapCounts[s] gives the total number of nonGaps in the full tree
   */
  private int [] nonGapCounts;

  /**
   * Determine for each vertex and site if this point is a common ancestor to all
   * nucleotides observed in this column
   * @param t
   */
  private void prepareCommonAncestorIndicator(PIPTreeNode t) 
  {
    int [] newDescCount = new int[nMSAColumns];
    descCounts.put(t,newDescCount);
    
    if (isLeaf(t))
    {
      double [][] felArray = logPeelingDetails.get(t);
      for (int s = 0; s < nMSAColumns; s++)
        newDescCount[s] = isGap(felArray, s) ? 0 : 1;
    }
    else
    {
      int [] leftDescCount = descCounts.get(left(t)),
             rightDescCount= descCounts.get(right(t));
      for (int s = 0; s < nMSAColumns; s++)
        newDescCount[s] = leftDescCount[s] + rightDescCount[s];
    }
    
    boolean [] commonAncestorIndicator = new boolean[nMSAColumns];
    commonAncestorIndicators.put(t, commonAncestorIndicator);
    for (int s = 0; s < nMSAColumns; s++)
      commonAncestorIndicator[s] = nonGapCounts[s] == newDescCount[s];
  }

  /**
   * Summarizes an unknown and unbounded number of unobserved columns with
   * no modern descendants
   * 
   * See arXiv:1207.6327
   * 
   * @param z
   * @return
   */
  private double logPhi(double z)
  {
    return (z - 1.0) * nu + nMSAColumns * logNu - SpecialFunctions.logFactorial(nMSAColumns);
  }

  private void computePeelingRecursions()
  { 
    logPeelingDetails = new HashMap<PIPTreeNode, double[][]>();
    logPeelingSummaries = new HashMap<PIPTreeNode, double[]>();
    nonGapCounts = new int[nSitesPlusFullGap];
    
    // Felsenstein recursions themselves
    for (PIPTreeNode t : postOrderTaxaTraversal)
    {
      double [][] currentFelsensteinArray = isLeaf(t) ?
          initFelsensteinArray(t) :
          PIPLikelihoodUtils.standardLogScaleFelsensteinRecursion(logPeelingDetails.get(left(t)), logPeelingDetails.get(right(t)), marginalLogPr(left(t)), marginalLogPr(right(t)));
      logPeelingDetails.put(t, currentFelsensteinArray);
      logPeelingSummaries.put(t, logPeelingSummary(currentFelsensteinArray));
    } 
  }

  /**
   * Marginalizes the value of the nucleotide at the root of each of the intermediate bottom up recursions
   * @param currentFelsensteinArray
   * @return
   */
  private double[] logPeelingSummary(double[][] currentFelsensteinArray)
  {
    double [] result = new double[currentFelsensteinArray.length];
    for (int s = 0; s < nSitesPlusFullGap; s++)
    {
      double siteLogLikelihood = Double.NEGATIVE_INFINITY;
      double [] curCache = currentFelsensteinArray[s];
      for (int x = 0; x < nCharacters; x++)
        siteLogLikelihood = NumericalUtils.logAdd(siteLogLikelihood, curCache[x] + logPi[x]); // logAdd computes log(exp(x) + exp(y)) while avoiding underflows
      if (isSiteLogLikelihoodValid(siteLogLikelihood))
        result[s] = siteLogLikelihood;
      else
        throw new RuntimeException("Missing data not currently supported");
    }
    return result;
  }

  /**
   * Computes for each site the dot product over the branches of the insertionLocationPrs 
   * (denoted by iota in arXiv:1207.6327) with the column conditional probibilities 
   * (denoted f in arXiv:1207.6327)
   * This is the sum in p.12 of arXiv:1207.6327
   * @return site indexed array of probabilities
   */
  private double[] computeLogZ_phase1()
  {
    double [] workArray = new double[nSitesPlusFullGap];
    for (int s = 0; s < nSitesPlusFullGap; s++)
      workArray[s] = Double.NEGATIVE_INFINITY;
    for (PIPTreeNode t : postOrderTaxaTraversal)
      computeLogZ_phase1(workArray, t);
    return workArray;
  }
  
  /**
   * Takes the product over the columns (but not for the pseudo-site corresponding to 
   * the case where all nucleotides are deleted before reaching the leaves)
   * @param array
   * @return
   */
  private double computeLogZ_phase2_nonGap(double [] array)
  {
    double logProduct = 0.0;
    for (int s = 0; s < nMSAColumns; s++)
      logProduct += array[s];
    return logProduct;
  }
  
  private double computeZ_phase2_gap(double [] array)
  {
    return  Math.exp(array[fullGapIndex]);
  }

  /**
   * @param workArray keeps track of partial sums (in log space)
   * @param t
   */
  private void computeLogZ_phase1(double [] workArray, PIPTreeNode t)
  {
    boolean [] commonAncestorIndicator = commonAncestorIndicators.get(t);
    double logInsertPr = logInsertLocationPrs.get(t);
    double [] logPeelingSummary = logPeelingSummaries.get(t);
    double logSurvivalPr = logSurvivalPrs.get(t);
    for (int s = 0; s < nMSAColumns; s++)
      workArray[s] = NumericalUtils.logAdd( // logAdd computes log(exp(x) + exp(y)) while avoiding underflows
          commonAncestorIndicator[s] ? logInsertPr + logPeelingSummary[s] + logSurvivalPr : Double.NEGATIVE_INFINITY,
          workArray[s]);
    // logColumnConditionalPr is f_v in arXiv:1207.6327
    double logColumnConditionalPr = Math.log(1.0 + Math.exp(logSurvivalPr) * (Math.exp(logPeelingSummary[fullGapIndex])-1.0));
    workArray[fullGapIndex] = NumericalUtils.logAdd( // logAdd computes log(exp(x) + exp(y)) while avoiding underflows
        logInsertPr + logColumnConditionalPr, 
        workArray[fullGapIndex]);
  }
  
  private void prepareLogInsertLocationPrs()
  {
    logInsertLocationPrs = new HashMap<PIPTreeNode, Double>();
    double prefix = - Math.log(totalTreeLength + 1.0/mu);
    for (PIPTreeNode t : postOrderTaxaTraversal)
      logInsertLocationPrs.put(t, prefix +
          (isRoot(t) ?
              - logMu :
              Math.log(branchLength(t))));
  }

  private void prepareLogSurvivalPrs()
  {
    logSurvivalPrs = new HashMap<PIPTreeNode, Double>();
    for (PIPTreeNode t : postOrderTaxaTraversal)
    {
      final double bl = isRoot(t) ? Double.NaN : branchLength(t);
      logSurvivalPrs.put(t, 
          isRoot(t) ? 
              0.0 : 
              Math.log(1.0 - Math.exp(- mu * bl)) - logMu - Math.log(bl));
    }
  }
  
  private void prepareCommonAncestorIndicators()
  {
    descCounts = new HashMap<PIPTreeNode, int[]>();
    commonAncestorIndicators = new HashMap<PIPTreeNode, boolean[]>();
    for (PIPTreeNode t : postOrderTaxaTraversal)
      prepareCommonAncestorIndicator(t);
  }

  private double[][] marginalLogPr(PIPTreeNode t)
  {
    /*
     * NOTE: Q is not the substitution rate matrix itself, but it is easily constructed from it and the deletion rate
     * See arXiv:1207.6327
     */
    double [][] Q = pip.Q; 
    double [][] prs = MatrixFunctions.expm(new DoubleMatrix(Q).muli(branchLength(t))).toArray2();
    return log(prs);
  }
  
  public static double[][] log(double[][] ori)
  {
    double [][] processed = new double[ori.length][];
    for (int s = 0; s < processed.length; s++)
    {
      processed[s] = new double[ori[s].length];
      for (int c = 0; c < processed[s].length; c++)
        processed[s][c] = Math.log(ori[s][c]);
    }
    return processed;
  }
  
  private boolean isGap(double[][] felArray, int s)
  {
    return felArray[s][pip.gapIndex] == 0.0; // these are in log space!
  }
  
  private boolean isRoot(PIPTreeNode t)
  {
    return t.equals(root);
  }

  private double branchLength(PIPTreeNode t)
  {
    return branchLengths.get(t);
  }

  private PIPTreeNode right(PIPTreeNode t)
  {
    return lefts.get(t);
  }

  private PIPTreeNode left(PIPTreeNode t)
  {
    return rights.get(t);
  }

  private boolean isLeaf(PIPTreeNode t)
  {
    return leaves.contains(t);
  }

  private double [][] initFelsensteinArray(PIPTreeNode t) // note: this computes nonGapCounts at the same time
  {
    double [][] indicators = linearizedMSA.indicators(new SequenceId(t.toString()), pip.indexer, pip.gapIndex);
    double [][] result = new double[nSitesPlusFullGap][nCharacters+1];
    for (int s = 0; s < nMSAColumns; s++)
      for (int x = 0; x < nCharacters+1; x++)
      {
        double value = indicators[s][x];
        if (value == 1.0 && x < nCharacters)
          nonGapCounts[s]++;
        if (value != 0.0 && value != 1.0)
          throw new RuntimeException();
        result[s][x] = Math.log(value);
      }
    
    for (int x = 0; x < nCharacters+1; x++)
      result[fullGapIndex][x] =  Math.log(x == pip.gapIndex ? 1.0 : 0.0);
        
    return result;
  }

  private static boolean isSiteLogLikelihoodValid(double number)
  {
    return number <= 0.000001;
  }
  
  

}
