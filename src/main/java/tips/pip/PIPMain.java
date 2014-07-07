package tips.pip;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.apache.commons.lang3.tuple.Pair;

import tips.ImportanceSampler;
import tips.PotPropOptions;
import tips.PotProposal;

import briefj.collections.Counter;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import muset.MSAPoset;
import muset.SequenceId;



public class PIPMain implements Runnable
{
  @Option public double lambda =  2.0;
  @Option public double mu = 0.5;
  @Option public double bl =  0.3;
  @Option public Random genRand = new Random(1);
  @Option public int nParticles = 1000;
  
  @OptionSet(name = "potential") 
  public PotPropOptions potentialProposalOptions = new PotPropOptions();

  private MSAPoset fullGeneratedPath;
  private MSAPoset generatedEndPoints;
  private Pair<PIPString,PIPString> endPoints;

  @Override
  public void run()
  {
    generateNextData();
    ImportanceSampler<PIPString> sampler = buildImportanceSampler();
    sampler.nParticles = nParticles;
    Counter<List<PIPString>> sample = sampler.sample(getStart(), getEnd(), bl);
    System.out.println(sampler.estimateZ(sample));
  }
  
  public ImportanceSampler<PIPString> buildImportanceSampler()
  {
    PIPProcess process = getProcess();
    PIPPotential pot = new PIPPotential();
    PotProposal<PIPString> prop = new PotProposal<PIPString>(process, pot, potentialProposalOptions);
    return new ImportanceSampler<PIPString>(prop, process);
  }

  public void generateNextData()
  {
    PIPProcess process = getProcess();
    fullGeneratedPath = process.sample(genRand, bl);
    generatedEndPoints = PIPProcess.keepOnlyEndPts(fullGeneratedPath, PIPMain.ta, PIPMain.tb);
    endPoints = getEndPoints(generatedEndPoints, PIPMain.ta, PIPMain.tb);
  }

  private PIPProcess getProcess()
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
