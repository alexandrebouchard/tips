package tips;

import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;



/**
 * The implementation proposed in the paper is PotProposal, but 
 * other proposals could theoretically be used while still integrating 
 * over the sojourn times.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 * @param <S>
 */
public interface Proposal<S>
{
  /**
   * 
   * @param rand
   * @param x
   * @param y
   * @param t
   * @return A pair where the first item is the proposed path, 
   *   and the second, its proposal probability
   */
  public Pair<List<S>, Double> propose(Random rand, S x, S y, double t);
}
