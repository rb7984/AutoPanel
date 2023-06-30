using System;
using System.Collections;
using System.Collections.Generic;

using Rhino;
using Rhino.Geometry;

using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;

using System.Linq;
using System.Windows.Forms;


/// <summary>
/// This class will be instantiated on demand by the Script component.
/// </summary>
public abstract class Script_Instance_ff93f : GH_ScriptInstance
{
  #region Utility functions
  /// <summary>Print a String to the [Out] Parameter of the Script component.</summary>
  /// <param name="text">String to print.</param>
  private void Print(string text) { /* Implementation hidden. */ }
  /// <summary>Print a formatted String to the [Out] Parameter of the Script component.</summary>
  /// <param name="format">String format.</param>
  /// <param name="args">Formatting parameters.</param>
  private void Print(string format, params object[] args) { /* Implementation hidden. */ }
  /// <summary>Print useful information about an object instance to the [Out] Parameter of the Script component. </summary>
  /// <param name="obj">Object instance to parse.</param>
  private void Reflect(object obj) { /* Implementation hidden. */ }
  /// <summary>Print the signatures of all the overloads of a specific method to the [Out] Parameter of the Script component. </summary>
  /// <param name="obj">Object instance to parse.</param>
  private void Reflect(object obj, string method_name) { /* Implementation hidden. */ }
  #endregion

  #region Members
  /// <summary>Gets the current Rhino document.</summary>
  private readonly RhinoDoc RhinoDocument;
  /// <summary>Gets the Grasshopper document that owns this script.</summary>
  private readonly GH_Document GrasshopperDocument;
  /// <summary>Gets the Grasshopper script component that owns this script.</summary>
  private readonly IGH_Component Component;
  /// <summary>
  /// Gets the current iteration count. The first call to RunScript() is associated with Iteration==0.
  /// Any subsequent call within the same solution will increment the Iteration count.
  /// </summary>
  private readonly int Iteration;
  #endregion
  /// <summary>
  /// This procedure contains the user code. Input parameters are provided as regular arguments,
  /// Output parameters as ref arguments. You don't have to assign output parameters,
  /// they will have a default value.
  /// </summary>
  #region Runscript
  private void RunScript(List<Polyline> x, List<string> y, ref object A)
  {
    alert.Clear();

    if (x.Count == y.Count)
    {
      for (int i = 0; i < x.Count; i++)
      {
        Check(x[i], y[i]);
      }
    }

    A = alert;
  }
  #endregion
  #region Additional
  // START CODE

  List<Point3d> alert = new List<Point3d>();

  public void Check(Polyline pl, string name)
  {
    string type = name.Split('|')[0];
    double len;
    double tolerance;

    if (type.Contains('A'))
    {
      tolerance = 90;
      len = pl[0].DistanceTo(pl[1]);
      double test = len + tolerance;

      if (test > 1251)
      {
        alert.Add(new Point3d(pl.BoundingBox.Center));
      }
    }
  }

  // END CODE
  #endregion
}