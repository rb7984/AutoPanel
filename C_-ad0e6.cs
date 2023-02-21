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
using Rhino.FileIO;


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
  private void RunScript(DataTree<object> x, List<Polyline> y, ref object A, ref object B, ref object C)
  {
    grid = new Grid(x, y);
    panelC41s = grid.panels;

    A = grid.Plines(panelC41s);
    B = grid.FreeTag(panelC41s);
    C = Archive(panelC41s);
  }
  #endregion
  #region Additional

  // fields
  public Grid grid;
  public List<PanelC41> panelC41s;
  public List<Polyline> plines;
  public List<bool> freeTag;

  public List<string> archive;

  public class Grid
  {
    public DataTree<object> param;
    List<Polyline> obs;
    public List<double> xCoordinates;
    public List<double> yCoordinates;

    public List<PanelC41> panels;
    public List<Point3d> pts;

    public Grid(DataTree<object> x, List<Polyline> y)
    {
      param = x;
      obs = y;

      OrderLines(param);

      pts = new List<Point3d>();

      panels = Panels(xCoordinates, yCoordinates);
    }

    public void OrderLines(DataTree<object> param)
    {
      List<double> xtmp = new List<double>();
      List<double> ytmp = new List<double>();

      for (int i = 2; i < param.BranchCount; i++)
      {
        int fuga = (int) Convert.ToInt64(param.Branch(i)[0].ToString().Split('-')[2].ToString().Split('.')[0]);
        string dir = param.Branch(i)[0].ToString().Split('-')[1];

        List<double> tmp = new List<double>();

        if (dir == "h")
        {
          for (int j = 1; j < param.Branch(i).Count; j++)
          {
            PolylineCurve tmpPl = (PolylineCurve) param.Branch(i)[j];
            tmp.Add(tmpPl.PointAt(0).Y - fuga);
            tmp.Add(tmpPl.PointAt(0).Y + fuga);
          }

          ytmp.AddRange(tmp);
        }

        else if (dir == "v")
        {
          for (int j = 1; j < param.Branch(i).Count; j++)
          {
            PolylineCurve tmpPl = (PolylineCurve) param.Branch(i)[j];
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

      for (int i = 0; i < yCoordinates.Count; i += 2)
      {
        for (int j = 0; j < xCoordinates.Count; j += 2)
        {
          var points = new List<Point3d>
            {
              new Point3d(xCoordinates[j], yCoordinates[i], 0),
              new Point3d(xCoordinates[j + 1], yCoordinates[i], 0),
              new Point3d(xCoordinates[j + 1], yCoordinates[i + 1], 0),
              new Point3d(xCoordinates[j], yCoordinates[i + 1], 0),
              new Point3d(xCoordinates[j], yCoordinates[i], 0)
              };

          pts.AddRange(points.GetRange(0, 4));
          panels.Add(new PanelC41(points, obs));
        }
      }

      return panels;
    }

    public List<Polyline> Plines(List<PanelC41> panels)
    {
      List<Polyline> polylines = panels.Select(i => i.pl).ToList();

      return polylines;
    }
    public List<bool> FreeTag(List<PanelC41> panels)
    {
      List<bool> freeTag = panels.Select(i => i.crossed).ToList();

      return freeTag;
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
    public bool crossed;
    public string name;

    //Constructor
    public PanelC41(List<Point3d> list, List<Polyline> obs)
    {
      borderUp = 0;
      borderDown = 0;
      hookUp = 0;
      hookDown = 0;

      pl = new Polyline(list);
      Intercept(obs);
      Name();
    }

    public void Intercept(List<Polyline> obs)
    {
      int a = 0;

      for (int i = 0; i < obs.Count; i++)
      {
        a += Intersection.CurveCurve(pl.ToPolylineCurve(), obs[i].ToPolylineCurve(), 0.1, 0.1).Count;
      }

      if (a == 0) crossed = false;
      else crossed = true;
    }

    public void Name()
    {
      width = Math.Round(pl[0].DistanceTo(pl[1]));
      height = Math.Round(pl[1].DistanceTo(pl[2]));
      name = width.ToString() + '-' + height.ToString();
    }
  }

  public List<string> Archive(List<PanelC41> panelC41s)
  {
    List<string> list;
    var orderedPanels = panelC41s.OrderBy(panelC41 => panelC41.width).ThenBy(panelC41 => panelC41.height);

    list = orderedPanels.Select(x => x.name).ToList();

    List<string> arc = new List<string> { list[0] };
    for (int i = 1; i < list.Count; i++)
    {
      if (list[i] != list[i - 1])
      {
        arc.Add(list[i]);
      }
    }

    return arc;
  }
  #endregion
}