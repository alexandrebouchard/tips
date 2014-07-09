package tips;


import briefj.collections.Counter;




/**
 * A Continuous time Markov chain.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 * @param <S> The type of the states
 */
public interface Process<S>
{
  public Counter<S> rates(S point);
}