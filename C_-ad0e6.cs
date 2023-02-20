using System;
using System.Collections;
using System.Collections.Generic;

using Rhino;
using Rhino.Geometry;

using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;

using Rhino.Geometry.Intersect;
using System.Linq;

/// <summary>
/// This class will be instantiated on demand by the Script component.
/// </summary>
public abstract class Script_Instance_ad0e6 : GH_ScriptInstance
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
  private void RunScript(DataTree<object> x, object y, ref object A)
  {
    grid = new Grid(x);
    panelC41s = grid.panels;

  }
  #endregion
  #region Additional

  // fields
  public Grid grid;
  public List<PanelC41> panelC41s;

  public class Grid
  {
    public DataTree<object> param;
    public List<double> xCoordinates;
    public List<double> yCoordinates;

    public List<PanelC41> panels;

    public Grid(DataTree<object> x)
    {
      param = x;
      OrderLines(param);

      panels = Panels(xCoordinates, yCoordinates);
    }

    public void OrderLines(DataTree<object> param)
    {
      List<double> xtmp = new List<double>();
      List<double> ytmp = new List<double>();

      for (int i = 2; i < param.BranchCount; i++)
      {
        int fuga = (int)Convert.ToInt64(param.Branch(i)[0].ToString().Split('-')[2].ToString().Split('.')[0]);
        string dir = param.Branch(i)[0].ToString().Split('-')[1];

        List<double> tmp = new List<double>();

        if (dir == "h")
        {
          for (int j = 1; j < param.Branch(i).Count; j++)
          {
            Polyline tmpPl = (Polyline)param.Branch(i)[j];
            tmp.Add(tmpPl.PointAt(0).Y - fuga);
            tmp.Add(tmpPl.PointAt(0).Y + fuga);
          }

          ytmp.AddRange(tmp);
        }

        else if (dir == "v")
        {
          for (int j = 1; j < param.Branch(i).Count; j++)
          {
            Polyline tmpPl = (Polyline)param.Branch(i)[j];
            tmp.Add(tmpPl.PointAt(0).X - fuga);
            tmp.Add(tmpPl.PointAt(0).X + fuga);
          }

          xtmp.AddRange(tmp);
        }
      }

      xtmp.Sort();
      ytmp.Sort();

      xCoordinates = xtmp;
      yCoordinates = ytmp;

      xCoordinates.RemoveAt(xCoordinates.Count - 1);
      yCoordinates.RemoveAt(yCoordinates.Count - 1);
      xCoordinates.RemoveAt(0);
      yCoordinates.RemoveAt(0);

      //sortedList = pts.OrderBy(point => point.X).ToList();
    }

    public List<PanelC41> Panels(List<double> xCoordinates, List<double> yCoordinates)
    {
      List<PanelC41> panels = new List<PanelC41>();
      for (int i = 0; i < yCoordinates.Count; i= i+2)
      {
        for (int j = 0; j < xCoordinates.Count ;j = j + 2)
        {
          List<Point3d> points = new List<Point3d>
          {
            new Point3d(xCoordinates[j], yCoordinates[i],0),
            new Point3d(xCoordinates[j+1], yCoordinates[i],0),
            new Point3d(xCoordinates[j+1], yCoordinates[i+1],0),
            new Point3d(xCoordinates[j], yCoordinates[i + 1],0),
            new Point3d(xCoordinates[j], yCoordinates[i],0)
          };

          panels.Add(new PanelC41(points));
        }
      }

      return panels;
    }
  }

  public class PanelC41
  {
    //fields
    public double width;
    public double height;
    public double borderUp;
    public double borderDown;
    public double hookUp;
    public double hookDown;
    public Polyline pl;

    //Constructor
    public PanelC41(double width, double height, double borderUp, double borderDown, double hookUp, double hookDown, Polyline pl)
    {
      this.width = width;
      this.height = height;
      this.borderUp = borderUp;
      this.borderDown = borderDown;
      this.hookUp = hookUp;
      this.hookDown = hookDown;
      this.pl = pl;
    }
    public PanelC41(Polyline pl)
    {
      width = 0;
      height = 0;
      borderUp = 0;
      borderDown = 0;
      hookUp = 0;
      hookDown = 0;
      this.pl = pl;
    }
    public PanelC41(List<Point3d> list)
    {
      width = 0;
      height = 0;
      borderUp = 0;
      borderDown = 0;
      hookUp = 0;
      hookDown = 0;

      pl = new Polyline(list);
    }

  }

  #endregion
}