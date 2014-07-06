package tips.pip.ref;





/**
 * The name of a node on a tree. Can be an informative string for the leaves of the 
 * tree, or some id for the internal nodes. We call the former case labelled,
 * and the latter case, unlabelled.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public final class PIPTreeNode
{
  /**
   * Create a labelled node with the given name.
   * Usually for the leaves only.
   * @param name
   * @return
   */
  public static PIPTreeNode withLabel(String name)
  {
    return new PIPTreeNode(name);
  }
  
  /**
   * Creates a unique internal node.
   * @return
   */
  public static PIPTreeNode nextUnlabelled()
  {
    synchronized (INTERNAL_PREFIX)
    {
      return unlabelled(nextId++);
    }
  }
  
  /**
   * Re-create an internal node with the given index
   * @param index
   * @return
   */
  public static PIPTreeNode unlabelled(int index)
  {
    return withLabel(INTERNAL_PREFIX + index);
  }
  
  @Override
  public String toString()
  {
    return description;
  }
  
  public boolean isLabelled() { return isLabelled; }
  
  private static int nextId = 0;
  
  private final String description;
  private final boolean isLabelled;
  
  private static final String INTERNAL_PREFIX = "unlabelled_";
  
  private PIPTreeNode(String description)
  {
    this.description = description;
    this.isLabelled = !description.startsWith(INTERNAL_PREFIX);
  }

  @Override
  public int hashCode()
  {
    final int prime = 31;
    int result = 1;
    result = prime * result
        + ((description == null) ? 0 : description.hashCode());
    result = prime * result + (isLabelled ? 1231 : 1237);
    return result;
  }

  @Override
  public boolean equals(Object obj)
  {
    if (this == obj)
      return true;
    if (obj == null)
      return false;
    if (getClass() != obj.getClass())
      return false;
    PIPTreeNode other = (PIPTreeNode) obj;
    if (description == null)
    {
      if (other.description != null)
        return false;
    } else if (!description.equals(other.description))
      return false;
    if (isLabelled != other.isLabelled)
      return false;
    return true;
  }
  
  
}
