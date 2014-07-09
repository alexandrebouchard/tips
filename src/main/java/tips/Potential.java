package tips;


/**
 * A potential used by PotProposal to guide proposed jump chains
 * towards the target.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 * @param <S>
 */
public interface Potential<S>
{
  /**
   * Lower value should roughly correspond to states closer to 
   * the target. 
   * 
   * Current limitation: we assume the value can only decrease or
   * increase by one. A value of POSITIVE_INFINITY can be used to 
   * denote going to an absorbing state/region that makes it 
   * impossible to reach the target.
   * 
   * @param proposed potential next state in the proposed sequence
   *   being built
   * @param target
   * @return potential value for proposed (relative to target)
   */
  public double get(S proposed, S target);
}
