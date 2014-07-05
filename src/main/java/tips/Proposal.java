package tips;

import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;




public interface Proposal<S>
{
  public Pair<List<S>, Double> propose(Random rand, S x, S y, double t);
}
