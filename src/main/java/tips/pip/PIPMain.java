package tips.pip;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.apache.commons.lang3.tuple.Pair;

import tips.Proposal;
import tips.TimeIntegratedPathSampler;
import tips.utils.PotPropOptions;
import tips.utils.PotProposal;

import briefj.BriefCollections;
import briefj.collections.Counter;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import muset.MSAPoset;
import muset.MSAPoset.Column;
import muset.SequenceId;
import muset.util.Edge;



public class PIPMain implements Runnable
{
  @Option public double lambda =  2.0;
  @Option public double mu = 0.5;
  @Option public double bl =  0.3;
  @Option public Random genRand = new Random(1);
  @Option public int nParticles = 1000;
  
  /**
   * Ensure that the generated alignment has a unique 
   * linearization. This is not required by the TIPS method,
   * but when comparing to the analytic PIP computation, the 
   * analytic formula computes the probability for one fixed
   * linearization, so this is only needed for these 
   * correctness checks.
   */
  @Option public boolean ensureLinearizationUnique = false;
  
  @OptionSet(name = "potential") 
  public PotPropOptions potentialProposalOptions = new PotPropOptions();

  private MSAPoset fullGeneratedPath;
  private MSAPoset generatedEndPoints;
  private Pair<PIPString,PIPString> endPoints;
  

  @Override
  public void run()
  {
    generateNextData();
    TimeIntegratedPathSampler<PIPString> sampler = buildImportanceSampler();
    sampler.nParticles = nParticles;
    Counter<List<PIPString>> sample = sampler.sample(getStart(), getEnd(), bl);
    System.out.println(sampler.estimateZ(sample));
  }
  
  public TimeIntegratedPathSampler<PIPString> buildImportanceSampler()
  {
    return new TimeIntegratedPathSampler<PIPString>(getProposal(), getProcess());
  }
  
  public Proposal<PIPString> getProposal()
  {
    PIPPotential pot = new PIPPotential();
    return new PotProposal<PIPString>(getProcess(), pot, potentialProposalOptions);
  }


  public void generateNextData()
  {
    PIPProcess process = getProcess();
    fullGeneratedPath = process.sample(genRand, bl);
    generatedEndPoints = PIPProcess.keepOnlyEndPts(fullGeneratedPath, PIPMain.ta, PIPMain.tb);
    if (ensureLinearizationUnique)
      generatedEndPoints = ensureUniqueLinearization(generatedEndPoints);
    endPoints = getEndPoints(generatedEndPoints, PIPMain.ta, PIPMain.tb);
  }

  private static MSAPoset ensureUniqueLinearization(MSAPoset endPts)
  {
    MSAPoset result = new MSAPoset(endPts);
    
    if (endPts.nSequences() != 2)
      throw new RuntimeException();
    
    Column previousIndel = null;
    for (Column c : endPts.linearizedColumns())
    {
      if (c.getPoints().size() == 2)
        previousIndel = null;
      else
      {
        if (previousIndel == null)
          previousIndel = c;
        else
        {
          Map<SequenceId, Integer>
            pts1 = previousIndel.getPoints(),
            pts2 = c.getPoints();
          SequenceId 
            s1 = BriefCollections.pick(pts1.keySet()),
            s2 = BriefCollections.pick(pts2.keySet());
          if (s1.equals(s2))
          {
            // both inserts or both deletes
            previousIndel = c;
          }
          else
          {
            // need to connect to aboid multiple linearizations
            if (!result.tryAdding(new Edge(pts1.get(s1), pts2.get(s2), s1, s2)))
              throw new RuntimeException();
            previousIndel = null;
          }
        }
      }
    }
    
    return result;
  }

  public PIPProcess getProcess()
  {
    return new PIPProcess(lambda, mu);
  }
  
  public static Pair<PIPString, PIPString> getEndPoints(MSAPoset msa, SequenceId ta,
      SequenceId tb)
  {
    if (msa.sequences().size() != 2)
      throw new RuntimeException();
    Set<Integer> aligned1 = new HashSet<Integer>(), aligned2 = new HashSet<Integer>();
    for (muset.MSAPoset.Column c : msa.columns())
      if (c.getPoints().size() == 2)
      {
        if (c.getPoints().containsKey(ta)) aligned1.add(c.getPoints().get(ta));
        if (c.getPoints().containsKey(tb)) aligned2.add(c.getPoints().get(tb));
      }
    
    List<Integer> cur; PIPString pips1, pips2;
    
    cur = new ArrayList<Integer>(); for (int i = 0; i < msa.sequences().get(ta).length(); i++) cur.add(aligned1.contains(i) ? 0 : -1); pips1 = new PIPString(cur);
    cur = new ArrayList<Integer>(); for (int i = 0; i < msa.sequences().get(tb).length(); i++) cur.add(aligned2.contains(i) ? 0 : +1); pips2 = new PIPString(cur);
    
    return Pair.of(pips1, pips2);
  }

  public static void main(String [] args)
  {
    Mains.instrumentedRun(args, new PIPMain());
  }

  public MSAPoset getGeneratedEndPoints()
  {
    ensureGenerated();
    return generatedEndPoints;
  }
  
  public PIPString getStart()
  {
    ensureGenerated();
    return endPoints.getLeft();
  }
  
  public PIPString getEnd()
  {
    ensureGenerated();
    return endPoints.getRight();
  }
  
  private void ensureGenerated()
  {
    if (generatedEndPoints == null)
      generateNextData();
  }
  
  public static SequenceId ta = new SequenceId("A"), tb = new SequenceId("B");
}
