import JaCoP.core.Var;
import JaCoP.search.SimpleSolutionListener;
import JaCoP.search.Search;
import JaCoP.search.SelectChoicePoint;

class SolutionListener<T extends Var> extends SimpleSolutionListener<T> {
  private long data;

  private static native boolean handleSolution(long data);

  SolutionListener(long data) { this.data = data; }
  
  public boolean executeAfterSolution(
      Search<T> search, SelectChoicePoint<T> select) {
    if (handleSolution(data))
      throw new InterruptSearch();
    return super.executeAfterSolution(search, select);
  }
}
