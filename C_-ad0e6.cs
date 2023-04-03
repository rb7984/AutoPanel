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
using Rhino.DocObjects;
using Rhino.UI.Controls;
using Grasshopper.Kernel.Types.Transforms;
using System.Diagnostics;
using Rhino.UI;
using Rhino.Render.ChangeQueue;
using System.ComponentModel;
using System.Drawing.Printing;
using System.Text;
using System.Text.RegularExpressions;


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
  private void RunScript(DataTree<object> x, DataTree<object> y, List<Polyline> z, DataTree<object> u, ref object B, ref object A, ref object baselines, ref object pl, ref object layerNames, ref object names, ref object archive, ref object freeTag, ref object export, ref object normalVectors)
  {
    facade = new Facade(x, y, z, u);

    panelC41s = new List<PanelC41>();
    plines = new DataTree<Polyline>();
    this.archive = new List<string>();
    this.names = new DataTree<string>();
    types = new DataTree<string>();
    toExport = new List<string> { "Type,Width,Heigh,Marca,Facciata,Tin,TinH,Sec" };
    normals = new List<Vector3d>();

    foreach (Grid3d grid in facade.grids)
    {
      int i = facade.grids.IndexOf(grid);
      panelC41s.AddRange(grid.panels);
      plines.AddRange(grid.Plines(grid.panels), new GH_Path(i));
      types.AddRange(grid.TypeTags(grid.panels), new GH_Path(i));
      this.names.AddRange(grid.NameTags(grid.panels), new GH_Path(i));
      toExport.AddRange(grid.ToExport(grid.panels));
      normals.Add(grid.normal);
    }

    this.archive.AddRange(ArchiveTypes(panelC41s));

    B = facade.ints;
    A = facade.angles;
    baselines = facade.baseLines;
    pl = plines;
    layerNames = types;
    names = this.names;
    archive = this.archive;
    export = toExport;
    normalVectors = normals;
  }
  #endregion
  #region Additional

  // fields
  public Facade facade;
  public List<PanelC41> panelC41s;
  public DataTree<Polyline> plines;
  public DataTree<string> names;
  public DataTree<string> types;
  public List<string> toExport;
  public List<Vector3d> normals;

  public List<string> archive;

  public class Facade
  {
    //fields
    public List<Grid3d> grids;
    public List<Point3d> points;
    public List<double> angles;
    public List<Polyline> baseLines;

    public int cCount;

    //DELETE
    public List<int> ints;

    //contructor
    public Facade(DataTree<object> input, DataTree<object> obs, List<Polyline> baseline, DataTree<object> sec)
    {
      grids = new List<Grid3d>();
      cCount = 0;

      // counter = [3, 3, 3]
      var counter = FacadeCount(input);
      var pointer = 0;

      // Angles
      FillAngles(baseline);

      // Populate Datatrees for Grids
      for (int i = 0; i < counter.Count; i++)
      {
        DataTree<object> tmp = new DataTree<object>();
        DataTree<object> tmpObs = new DataTree<object>();
        var count = 0;

        for (int j = pointer; j < pointer + counter[i]; j++)
        {
          tmp.AddRange(input.Branch(j), new GH_Path(count)); count++;
        }
        pointer += counter[i];

        var k = i * 5;
        count = 0;
        for (int j = k; j < k + 5; j++)
        {
          tmpObs.AddRange(obs.Branch(j), new GH_Path(count));
          count++;
        }

        Grid3d tmpGrid = new Grid3d(tmp, tmpObs, baseLines[i]);
        grids.Add(tmpGrid);
      }

      // Left Border, Right Border, Name Uniformation and Sec
      PostPanels(grids, sec);

      // Uniform panelC41.ToExport
      CorrectExport();
    }

    public List<int> FacadeCount(DataTree<object> input)
    {
      List<int> k = new List<int> { 1 };
      List<double> compare = new List<double>();

      for (int i = 0; i < input.BranchCount; i++)
      {
        double a = char.GetNumericValue(input.Branch(i)[0].ToString().Split('-')[1][0]);
        compare.Add(a);
      }
      // k = [1]
      // compare = [0,0,0,1,1,1,2,2,2]

      for (int i = 1; i < compare.Count; i++)
      {
        if (compare[i] == compare[i - 1]) k[k.Count - 1]++;
        else k.Add(1);
      }

      // k = [3,3,3]
      ints = k;
      return k;
    }

    public void FillAngles(List<Polyline> baseline)
    {
      baseLines = new List<Polyline>();
      angles = new List<double>();

      foreach (Polyline p in baseline)
      {
        angles.Add(2 * Math.PI);

        for (int i = 1; i < p.Count - 1; i++)
        {
          baseLines.Add(new Polyline(new List<Point3d> { p[i - 1], p[i] }));
          Vector3d normal = Vector3d.CrossProduct(new Vector3d(p[i] - p[i - 1]), new Vector3d(0, 0, 1));

          if (Intersection.CurveCurve(new Line(p[i - 1], p[i + 1]).ToNurbsCurve(),
            new Line(p[i], normal, 100000).ToNurbsCurve(), 0.1, 0.1).Count > 0)
          {
            angles.AddRange(Enumerable.Repeat(Vector3d.VectorAngle(new Vector3d(p[i - 1] - p[i]), new Vector3d(p[i + 1] - p[i])), 2));
          }
          else
          {
            angles.AddRange(Enumerable.Repeat(2 * Math.PI - Vector3d.VectorAngle(new Vector3d(p[i - 1] - p[i]), new Vector3d(p[i + 1] - p[i])), 2));
          }
        }
        baseLines.Add(new Polyline(new List<Point3d> { p[p.Count - 2], p[p.Count - 1] }));

        angles.Add(2 * Math.PI);
      }
    }

    public void PostPanels(List<Grid3d> grid, DataTree<object> sec)
    {
      int fugaMax = 16;
      int tracker = 0;

      for (int j = 0; j < grid.Count; j++)
      {
        List<Polyline> secs = sec.Branch(j).Select(p => (PolylineCurve)p).ToList().Select(pp => pp.ToPolyline()).ToList();

        for (int i = 0; i < grid[j].panels.Count; i++)
        {
          double l = new Point3d(grid[j].panels[i].pl[0].X, grid[j].panels[i].pl[0].Y, 0).DistanceTo(new Point3d(grid[j].border[0].X, grid[j].border[0].Y, 0));
          double r = new Point3d(grid[j].panels[i].pl[1].X, grid[j].panels[i].pl[1].Y, 0).DistanceTo(new Point3d(grid[j].border[1].X, grid[j].border[1].Y, 0));

          // LEFT
          if (Math.Abs(l) < fugaMax)
          {
            var nn = BorderPanelName(grids.IndexOf(grid[j]), 0);
            if (grid[j].panels[i].type.Contains('F') || grid[j].panels[i].type.Contains('E'))
            {
              grid[j].panels[i].type = grid[j].panels[i].type.Replace("E", "K");
            }
            else { grid[j].panels[i].type = grid[j].panels[i].type.Replace("A", nn); }

            //BorderPanelReduction(grid[j].panels[i], j, 0);
          }

          // RIGHT
          if (Math.Abs(r) < fugaMax)
          {
            var nn = BorderPanelName(grids.IndexOf(grid[j]), 1);
            if (grid[j].panels[i].type.Contains('E') || grid[j].panels[i].type.Contains('F'))
            {
              grid[j].panels[i].type = grid[j].panels[i].type.Replace("E", "K");
            }
            else { grid[j].panels[i].type = grid[j].panels[i].type.Replace("A", nn); }

            //BorderPanelReduction(grid[j].panels[i], j, 1);
          }

          // C Tracker
          if (grid[j].panels[i].type.Contains("C"))
          {
            grid[j].panels[i].cCount = tracker;
            tracker++;
          }

          // Name correction (B*C)
          if (grid[j].panels[i].type.Contains("C*B"))
          {
            grid[j].panels[i].type = grid[j].panels[i].type.Replace("C*B", "B*C");
          }

          // SEC
          PanelSec(grid[j].panels[i], secs);

          // REDUCTION
          BorderPanelReduction(grid[j].panels[i], j);
        }
      }
    }

    public string BorderPanelName(int i, int a)
    {
      double angle = angles[(2 * i) + a];

      if (angle == 2 * Math.PI)
      {
        if (a == 0) return "E";
        else return "F";
      }
      // I and J Panels
      else if (grids[i].monopanel)
      {
        angle = (180 / Math.PI) * angle;
        angle = (360 - angle) / 2;
        if (a == 0) return "J" + "." + Math.Round(angle).ToString();
        else return "I" + "." + Math.Round(angle).ToString();
      }
      else
      {
        angle = (180 / Math.PI) * angle;
        angle = (360 - angle) / 2;
        if (a == 0) return "H" + "." + Math.Round(angle).ToString();
        else return "G" + "." + Math.Round(angle).ToString();
      }
    }

    public void BorderPanelReduction(PanelC41 panel, int gridIndex)
    {
      Point3d p = new Point3d();

      var reduce = new Action<PanelC41, int, int>((panelPass, reduction, leftRight) =>
      {
        Point3d tmpPoint = panel.pl[leftRight];
        Circle c = new Circle(tmpPoint, reduction);

        var events = Intersection.CurveCurve(panel.pl.ToNurbsCurve(), c.ToNurbsCurve(), 0.01, 0.01);

        p = new Point3d(events[0].PointA);

      });

      if (panel.type.Contains('.'))
      {
        string[] keys = new string[] { "G", "H", "I", "J" };

        string sKeyResult = keys.FirstOrDefault<string>(s => panel.type.Contains(s));

        switch (sKeyResult)
        {
          case "G":
          case "I":

            double angle = angles[gridIndex * 2 + 1];
            angle = (180 / Math.PI) * angle;
            angle = (360 - angle) / 2;

            if (Math.Abs(angle - 45) > 0.1 && angle > 0)
            {
              reduce(panel, 12, 1);

              List<Point3d> list = new List<Point3d>
              {
                panel.pl[0],
                p,
                new Point3d(p.X,p.Y, panel.pl[2].Z),
                panel.pl[3],
                panel.pl[0]
              };

              panel.pl = new Polyline(list);
            }
            else if (Math.Abs(angle - 45) < 0.1)
            {
              reduce(panel, 10, 1);

              List<Point3d> list = new List<Point3d>
              {
                panel.pl[0],
                p,
                new Point3d(p.X,p.Y, panel.pl[2].Z),
                panel.pl[3],
                panel.pl[0]
              };

              panel.pl = new Polyline(list);
            }

            break;

          case "H":
          case "J":

            double angle2 = angles[gridIndex * 2];
            angle2 = (180 / Math.PI) * angle2;
            angle2 = (360 - angle2) / 2;

            if (Math.Abs(angle2 - 45) > 0.1 && angle2 > 0)
            {
              reduce(panel, 12, 0);

              List<Point3d> list = new List<Point3d>
              {
                p,
                panel.pl[1],
                panel.pl[2],
                new Point3d(p.X, p.Y, panel.pl[3].Z),
                p
              };

              panel.pl = new Polyline(list);
            }
            else if (Math.Abs(angle2 - 45) < 0.1)
            {
              reduce(panel, 10, 0);

              List<Point3d> list = new List<Point3d>
              {
                p,
                panel.pl[1],
                panel.pl[2],
                new Point3d(p.X, p.Y, panel.pl[3].Z),
                p
              };

              panel.pl = new Polyline(list);
            }

            break;

          default:
            break;
        }

      }
    }

    //public Point3d PanelReduction(PanelC41 panel, int reduction, int leftRight)
    //{
    //  Point3d tmpPoint = panel.pl[leftRight];
    //  Circle c = new Circle(tmpPoint, reduction);

    //  var events = Intersection.CurveCurve(panel.pl.ToNurbsCurve(), c.ToNurbsCurve(), 0.01, 0.01);

    //  return new Point3d(events[0].PointA);
    //}

    public void PanelSec(PanelC41 panel, List<Polyline> pls)
    {
      foreach (Polyline pl in pls)
      {
        var tmp = Intersection.CurveCurve(pl.ToPolylineCurve(), panel.pl.ToPolylineCurve(), 0.1, 0.1);

        if (tmp.Count != 0)
        {
          panel.sec = Math.Round(panel.pl[0].DistanceTo(tmp[0].PointA));
          panel.tinH = panel.height - 6;
        }
      }
    }

    public void CorrectExport()
    {
      foreach (Grid3d grid in grids)
      {
        int a = grids.IndexOf(grid);
        foreach (PanelC41 p in grid.panels)
        {
          p.width = Math.Round(p.pl[0].DistanceTo(p.pl[1]));
          p.height = Math.Round(p.pl[1].DistanceTo(p.pl[2]));
          p.name = p.width.ToString() + '-' + p.height.ToString();

          p.toExcel = p.type + ","
            + p.width.ToString() + ","
            + p.height.ToString() + ","
            + p.type + "." + p.width.ToString() + "." + p.height.ToString() + ","
            + a.ToString() + ","
            + p.tin.ToString() + ","
            + p.tinH.ToString() + ","
            + p.sec.ToString();
        }
      }
    }
  }

  public class Grid3d
  {
    DataTree<object> param;
    public DataTree<Polyline> obs;

    public Polyline border;
    public Polyline height;
    public Vector3d normal;
    public Vector3d ax;

    // tutte le coordinate delle fughe in X e Z
    List<Point3d> pointCoordinates;
    List<double> zCoordinates;

    public List<PanelC41> panels;
    public List<Point3d> pts;

    // DELETE
    public bool monopanel;

    public Grid3d(DataTree<object> x, DataTree<object> y, Polyline border)
    {
      param = x;
      obs = new DataTree<Polyline>();
      ax = new Vector3d(border[1] - border[0]);
      this.border = border;

      for (int i = 0; i < y.BranchCount; i++)
      { obs.AddRange(y.Branch(i).Select(pl => (PolylineCurve)pl).ToList().Select(pl => pl.ToPolyline()).ToList(), new GH_Path(i)); }

      Borders(x.Branch(0)[1], ax);
      OrderLines(param);

      pts = new List<Point3d>();

      panels = Panels(pointCoordinates, zCoordinates);

      //TEST
      OrderOutputH(panels);
      PostPanelsWindows(panels);
      //TEST

      OrderOutputV(panels);

      PostPanels(panels);
    }

    // Assegna height e normal
    public void Borders(object a, Vector3d vector)
    {
      PolylineCurve p = (PolylineCurve)a;
      Polyline tmpPl = p.ToPolyline();
      List<Point3d> list = tmpPl.Select(point => point).ToList();
      list.RemoveAt(list.Count - 1);
      list = list.OrderBy(point => point.Z).ToList();

      height = new Polyline(new List<Point3d> { border[0], new Point3d(border[0].X, border[0].Y, list[list.Count - 1].Z) });
      normal = Vector3d.CrossProduct(vector, new Vector3d(height[1] - height[0]));
    }

    // Crea Liste di coordinate dalle quali fare i pannelli
    public void OrderLines(DataTree<object> param)
    {
      int verticalPanelCount = 0;

      List<Point3d> ptmp = new List<Point3d>();
      List<double> ztmp = new List<double>();

      // i = 1 perché il primo è il bordo
      for (int i = 1; i < param.BranchCount; i++)
      {

        double fuga = Convert.ToInt64(param.Branch(i)[0].ToString().Split('-')[2].ToString().Split('.')[0]);
        if (fuga == 0) { fuga += 0.1; }
        char dir = param.Branch(i)[0].ToString().Split('-')[1][1];

        List<Point3d> tmpP = new List<Point3d>();

        if (dir == 'h')
        {
          // start at j = 1 perchè la prima entry è il nome del layer 
          for (int j = 1; j < param.Branch(i).Count; j++)
          {
            PolylineCurve tmpPl = (PolylineCurve)param.Branch(i)[j];
            ztmp.Add(tmpPl.PointAt(0).Z - fuga);
            ztmp.Add(tmpPl.PointAt(0).Z + fuga);
          }
        }

        else if (dir == 'v')
        {
          verticalPanelCount += param.Branch(i).Count;
          // start at j = 1 perchè la prima entry è il nome del layer
          for (int j = 1; j < param.Branch(i).Count; j++)
          {
            PolylineCurve tmpPl = (PolylineCurve)param.Branch(i)[j];
            Point3d tmpPoint = tmpPl.PointAtStart;
            Point3d tp = new Point3d(tmpPl.PointAt(0).X, tmpPl.PointAt(0).Y, border.PointAt(0).Z);
            Circle c = new Circle(tp, fuga);
            var events = Intersection.CurveCurve(border.ToNurbsCurve(), c.ToNurbsCurve(), 0.01, 0.01);
            foreach (IntersectionEvent ie in events) { tmpP.Add(new Point3d(ie.PointA)); }
          }

          ptmp.AddRange(tmpP);
        }
      }

      // Sorting of Z coordinates
      ztmp.Sort();
      ztmp.RemoveAt(0);
      ztmp.RemoveAt(0);
      ztmp.RemoveAt(ztmp.Count - 1);
      ztmp.RemoveAt(ztmp.Count - 1);
      double[] t = new double[] { height[0].Z, height[1].Z };
      ztmp.Insert(0, t.Min());
      ztmp.Add(t.Max());

      zCoordinates = ztmp;

      ptmp = ptmp.OrderBy(point => point.DistanceTo(border[0])).ToList();
      ptmp.RemoveAt(0);
      ptmp.RemoveAt(ptmp.Count - 1);
      ptmp.AddRange(new List<Point3d> { border.PointAt(0), border.PointAt(1) });
      pointCoordinates = ptmp.OrderBy(point => point.DistanceTo(border[0])).ToList();

      if (verticalPanelCount > 2) monopanel = false;
      else monopanel = true;
    }

    public List<PanelC41> Panels(List<Point3d> pointCoordinates, List<double> zCoordinates)
    {
      List<PanelC41> panels = new List<PanelC41>();

      for (int i = 0; i < zCoordinates.Count; i += 2)
      {
        for (int j = 0; j < pointCoordinates.Count; j += 2)
        {
          var points = new List<Point3d> {
                new Point3d(pointCoordinates[j].X, pointCoordinates[j].Y, zCoordinates[i]),
                new Point3d(pointCoordinates[j + 1].X, pointCoordinates[j + 1].Y, zCoordinates[i]),
                new Point3d(pointCoordinates[j + 1].X, pointCoordinates[j + 1].Y, zCoordinates[i + 1]),
                new Point3d(pointCoordinates[j].X, pointCoordinates[j].Y, zCoordinates[i + 1]),
                new Point3d(pointCoordinates[j].X, pointCoordinates[j].Y, zCoordinates[i])
          };

          pts.AddRange(points.GetRange(0, 4));

          PanelC41 tmpPanel = new PanelC41(points, obs, border);
          //EDIT
          //if (!tmpPanel.crossed[1]) panels.Add(tmpPanel);
          if (!(tmpPanel.crossed[0] && tmpPanel.crossed[1])) panels.Add(tmpPanel);
        }
      }

      return panels;
    }

    public void OrderOutputV(List<PanelC41> panels)
    {
      List<PanelC41> tmp;

      tmp = panels.OrderBy(a => new Point3d(a.firstCorner.X, a.firstCorner.Y, 0).DistanceTo(border[0])).ThenBy(b => b.firstCorner.Z).ToList();

      this.panels.Clear();

      this.panels.AddRange(tmp);
    }

    public void OrderOutputH(List<PanelC41> panels)
    {
      List<PanelC41> tmp;

      tmp = panels.OrderBy(a => a.firstCorner.Z).ThenBy(b => new Point3d(b.firstCorner.X, b.firstCorner.Y, 0).DistanceTo(border[0])).ToList();

      this.panels.Clear();

      this.panels.AddRange(tmp);
    }

    void PostPanelsWindows(List<PanelC41> panels)
    {
      List<int> tmp = new List<int>();

      // panels[i] = Panel inside of windows
      for (int i = 1; i < panels.Count - 1; i++)
      {
        if (panels[i].crossed[0] == false && panels[i].crossed[1] == true)
        {
          if (!panels[i - 1].type.Contains('E'))
          {
            panels[i - 1].type = panels[i - 1].type.Replace('A', 'F');
          }
          else
          {
            panels[i - 1].type = panels[i - 1].type.Replace('E', 'K');
          }
          if (!panels[i + 1].type.Contains('F'))
          {
            panels[i + 1].type = panels[i + 1].type.Replace('A', 'E');
          }
          else
          {
            panels[i + 1].type = panels[i + 1].type.Replace('F', 'K');
          }
        }
      }
    }

    void PostPanels(List<PanelC41> panels)
    {
      List<int> tmp = new List<int>();

      // Panel with horizontal interceptions
      for (int i = 1; i < panels.Count - 1; i++)
      {
        if (panels[i].crossed[0] == true && panels[i].crossed[1] == false)
        {
          panels[i - 1].type += "*W";
          panels[i + 1].type += "*B";

          if (panels[i].type.Contains("E") || panels[i].type.Contains("F") || panels[i].type.Contains("J") || panels[i].type.Contains("I") || panels[i].type.Contains("K"))
          {
            panels[i].crossed[0] = false;
          }
          else
          {
            tmp.Add(i);
          }
        }
        if (panels[i].crossed[0] == false && panels[i].crossed[1] == true)
        {
          panels[i - 1].type = panels[i - 1].type.Replace("A", "KZ");
          panels[i + 1].type = panels[i + 1].type.Replace("A", "KZ");

          tmp.Add(i);
        }
      }

      tmp = tmp.OrderByDescending(x => x).ToList();
      for (int i = 0; i < tmp.Count; i++)
      {
        panels.RemoveAt(tmp[i]);
      }

      // All panels
      for (int i = 0; i < panels.Count - 1; i++)
      {
        if (panels[i].type.Contains('B'))
        {
          panels[i + 1].type += "*D";
          panels[i + 1].tin = panels[i + 1].width;
          panels[i + 1].toExcel += "," + panels[i + 1].width.ToString();
        }
      }

      for (int i = 0; i < panels.Count; i++)
      {
        if (panels[i].type.Contains('W'))
        {
          panels[i].type = panels[i].type.Replace('W', 'B');
        }
        if (panels[i].type.Contains('D'))
        {
          panels[i].type = panels[i].type.Replace('D', 'B');
        }
        if (panels[i].type.Contains('Z'))
        {
          panels[i].type = panels[i].type.Replace("Z", "*B");
        }
      }

      // BorderUp and Down
      for (int i = 0; i < panels.Count; i++)
      {
        if (Math.Abs(panels[i].pl[0].Z - height[0].Z) < 0.01 && !panels[i].type.Contains('B'))
        {
          panels[i].type += "*B";
          if (!panels[i + 1].type.Contains('B'))
          {
            panels[i + 1].type += "*B";
            panels[i + 1].tin = panels[i + 1].width;
            panels[i + 1].toExcel += "," + panels[i + 1].width.ToString();
          }
        }
        else if (Math.Abs(panels[i].pl[3].Z - height[1].Z) < 0.01 && !panels[i].type.Contains('B'))
        {
          panels[i].type += "*B";
        }
      }
    }

    public List<Polyline> Plines(List<PanelC41> panels)
    {
      List<Polyline> polylines = panels.Select(i => i.pl).ToList();

      return polylines;
    }

    public List<string> NameTags(List<PanelC41> panels)
    {
      List<string> names = panels.Select(i => i.type).ToList();
      for (int i = 0; i < panels.Count; i++)
      {
        names[i] = names[i] + "|" + panels[i].name;
        if (panels[i].cCount != -1) names[i] += "/" + panels[i].cCount.ToString();
      }

      return names;
    }

    public List<string> TypeTags(List<PanelC41> panels)
    {
      //List<string> types = panels.Select(i => i.type.Split('.')[0]).ToList();
      List<string> types = new List<string>();
      foreach (PanelC41 panel in panels)
      {
        string[] type = panel.type.Split('.');
        var adding = type[0];
        if (type.Length > 1 && type[1].Contains('B'))
        {
          adding += "*B";
        }
        if (type.Length > 1 && type[1].Contains('C'))
        {
          adding += "*C";
        }
        types.Add(adding);
      }

      return types;
    }

    public List<string> ToExport(List<PanelC41> panels)
    {
      List<string> export = panels.Select(i => i.toExcel).ToList();

      return export;
    }

  }

  public class PanelC41
  {
    //fields
    public double width;
    public double height;
    public double tin;
    public double tinH;
    public double sec;
    public int cCount;

    public Polyline pl;
    public bool[] crossed;
    public string type;
    public string name;
    public Point3d firstCorner;
    public Polyline border;
    public Vector3d ax;
    public Point3d gridOrigin;
    public Point3d gridOriginZ0;

    public string toExcel;

    //constructor
    public PanelC41(List<Point3d> list, DataTree<Polyline> obs, Polyline border)
    {
      type = "A";
      tin = 0;
      tinH = 0;
      sec = 0;
      cCount = -1;

      crossed = new bool[2] { false, false };

      pl = new Polyline(list);

      firstCorner = pl[0];

      ax = new Vector3d(border[1] - border[0]);
      gridOrigin = border[0];
      gridOriginZ0 = new Point3d(border[0].X, border[0].Y, 0);
      this.border = border;

      List<int> wall = InterceptedPanelWall(obs.Branch(4));

      // [0 = internal,
      // 1 = verticalInterception,
      // 2 = horizontalInteception,
      // 3 = corner]

      if (!wall.Contains(0))
      {
        List<int> i = Intercept(obs);

        if (i[0] != 0)
        {
          for (int j = 1; j < i[0] + 1; j++)
          {
            crossed = InterceptedPanelVertical(obs.Branch(0), i[j], 8);
          }
        }
        if (i[i[0] + 1] != 0)
        {
          for (int j = i[i[0] + 1] + 1; j < i.Count; j++)
          {
            crossed = InterceptedPanelHorizontal(obs.Branch(2), i[j], 8, false);
          }
        }
        for (int j = 0; j < wall.Count; j++)
        {
          if (wall[j] == 1)
          {
            crossed = InterceptedPanelVertical(obs.Branch(4), j, 0);
            type = type.Replace('W', 'B');
          }
          if (wall[j] == 2) crossed = InterceptedPanelHorizontal(obs.Branch(4), j, 8, true);
        }
      }

      Name();
      Detached(obs.Branch(1));
      InterceptedPanelVarious(obs.Branch(3));

      toExcel = type + "," + width.ToString() + "," + height.ToString() +
        "," + type + "." + width.ToString() + "." + height.ToString();
    }

    public void Name()
    {
      width = Math.Round(pl[0].DistanceTo(pl[1]));
      height = Math.Round(pl[1].DistanceTo(pl[2]));
      name = width.ToString() + '-' + height.ToString();
    }

    public List<int> Intercept(DataTree<Polyline> obs)
    {
      // Intercept():
      // > plines that intesect the panel ==> var tmp
      // > plines that completely include the panel ==> var tmp2

      List<int> intercepted = new List<int>();
      intercepted.Add(0);

      // Branch(0) windows. Polylines centrate sul pannello verticalmente
      for (int i = 0; i < obs.Branch(0).Count; i++)
      {
        //var tmp = Intersection.CurveCurve(pl.ToPolylineCurve(), obs.Branch(0)[i].ToPolylineCurve(), 0.1, 0.1).Count;
        var tmp2 = GenericInterception(obs.Branch(0)[i]);

        //if (tmp != 0)
        if (tmp2)
        {
          intercepted.Add(i);
          intercepted[0]++;
        }
      }

      // Branch(2) pavements. Polylines centrate sul pannello orizzontalmente
      intercepted.Add(0);

      for (int i = 0; i < obs.Branch(2).Count; i++)
      {
        var tmp = Intersection.CurveCurve(pl.ToPolylineCurve(), obs.Branch(2)[i].ToPolylineCurve(), 0.1, 0.1).Count;

        if (tmp != 0)
        {
          intercepted.Add(i);
          intercepted[intercepted[0] + 1]++;
        }
      }

      // intercepted structure:
      //
      //   0  1  2  3  4  5  6
      //   |-----|  |--------|
      // 
      // [ 2, 2, 5, 3, 0, 7, 9]
      //   ^  '  '  ^  '  '  '
      //
      //   pointer  pointer

      return intercepted;
    }

    public bool[] InterceptedPanelVertical(List<Polyline> obs, int i, double fuga)
    {
      // Funziona solo per finestre centrate sui pannelli
      // Viene considerata solo la Z quindi

      double[] zObs = new double[] { obs[i][0].Z, obs[i][1].Z, obs[i][2].Z };
      double[] z = new double[] { pl[0].Z, pl[1].Z, pl[2].Z };

      // pannelli interamente contenuti
      if (zObs.Max() >= z.Max() && zObs.Min() <= z.Min())
      {
        return new bool[] { false, true };
      }
      //pannelli intersecati sotto
      else if (z.Max() > zObs.Max())
      {
        double newZ = zObs.Max();
        List<Point3d> tmp = new List<Point3d>()
        {
          new Point3d(pl[0].X, pl[0].Y, newZ + fuga),
          new Point3d(pl[1].X, pl[1].Y, newZ + fuga),
          new Point3d(pl[2]),
          new Point3d(pl[3]),
          new Point3d(pl[4].X, pl[4].Y, newZ + fuga),
        };

        pl = new Polyline(tmp);
        if (!type.Contains('W')) type += "*W";
        return new bool[] { false, false };
      }
      //pannelli intersecati sopra
      else if (z.Min() < zObs.Min())
      {
        double newZ = zObs.Min();
        List<Point3d> tmp = new List<Point3d>()
        {
          new Point3d(pl[0]),
          new Point3d(pl[1]),
          new Point3d(pl[2].X, pl[2].Y, newZ - fuga),
          new Point3d(pl[3].X, pl[3].Y, newZ - fuga),
          new Point3d(pl[4]),
        };

        pl = new Polyline(tmp);
        if (!type.Contains('W')) type += "*W";
        return new bool[] { false, false };
      }

      else { return new bool[] { false, false }; }
    }

    public bool[] InterceptedPanelHorizontal(List<Polyline> obs, int i, double fuga, bool wall)
    {
      // wall = false: per rettangoli orizzontali centrati sui pannelli
      // wall = true: rettangoli alti più di un pannello
      // Viene considerata solo distanza X-Y

      List<Point3d> p = new List<Point3d> { pl[0], pl[1] };

      List<Point3d> pObs = obs[i].Select(point => point).ToList();
      pObs.RemoveAt(pObs.Count - 1);
      pObs = pObs.OrderBy(point => point.Z).ToList();
      double[] zObs = new double[] { pObs[0].Z, pObs[3].Z };
      pObs.RemoveRange(2, 2);
      pObs = pObs.OrderBy(point => new Point3d(point.X, point.Y, 0).DistanceTo(gridOriginZ0)).ToList();

      bool leftObs = new Point3d(pObs[0].X, pObs[0].Y, 0).DistanceTo(gridOriginZ0) < 0.01;
      bool a = new Point3d(pl[1].X, pl[1].Y, 0).DistanceTo(gridOriginZ0) > new Point3d(pObs[1].X, pObs[1].Y, 0).DistanceTo(gridOriginZ0);
      bool b = new Point3d(pObs[1].X, pObs[1].Y, 0).DistanceTo(gridOriginZ0) > new Point3d(pl[0].X, pl[0].Y, 0).DistanceTo(gridOriginZ0);
      bool c = new Point3d(pl[1].X, pl[1].Y, 0).DistanceTo(gridOriginZ0) > new Point3d(pObs[0].X, pObs[0].Y, 0).DistanceTo(gridOriginZ0);
      bool d = new Point3d(pObs[0].X, pObs[0].Y, 0).DistanceTo(gridOriginZ0) > new Point3d(pl[0].X, pl[0].Y, 0).DistanceTo(gridOriginZ0);

      if (!leftObs && (c && d))
      {
        List<Point3d> tmp = new List<Point3d>()
        {
          new Point3d(pl[0]),
          new Point3d(pObs[0].X, pObs[0].Y, pl[1].Z),
          new Point3d(pObs[0].X, pObs[0].Y, pl[2].Z),
          new Point3d(pl[3]),
          new Point3d(pl[4])
        };

        pl = new Polyline(tmp);

        if (new Point3d(p[0].X, p[0].Y, 0).DistanceTo(gridOriginZ0) <= 0.01) { type = type.Replace("A", "J"); }
        else type = type.Replace("A", "F");
      }
      else if (leftObs && a && b)
      {
        List<Point3d> tmp = new List<Point3d>()
        {
          new Point3d(pObs[1].X, pObs[1].Y, pl[0].Z),
          new Point3d(pl[1]),
          new Point3d(pl[2]),
          new Point3d(pObs[1].X, pObs[1].Y, pl[3].Z),
          new Point3d(pObs[1].X, pObs[1].Y, pl[4].Z)
        };

        pl = new Polyline(tmp);

        if (new Point3d(pl[1].X, pl[1].Y, 0).DistanceTo(new Point3d(border[1].X, border[1].Y, 0)) < 0.01) { type = type.Replace("A", "I"); }
        else type = type.Replace("A", "E");
      }

      // wall = t: return [F,F] and keep the panel
      // wall = f: return [T,F] to name B panel (i+1), (i-1)
      if (wall) return new bool[] { false, false };
      else return new bool[] { true, false };
    }

    public void InterceptedPanelVarious(List<Polyline> obs)
    {
      double[] z = new double[] { pl[1].Z, pl[2].Z };
      int counter = 0;

      for (int i = 0; i < obs.Count; i++)
      {
        List<Point3d> pObs = new List<Point3d> { obs[i][0], obs[i][1], obs[i][2], obs[i][3] };
        pObs = pObs.OrderBy(point => point.Z).ToList();
        double[] zObs = pObs.GetRange(1, 2).Select(point => point.Z).ToArray();

        pObs.RemoveRange(2, 2);
        pObs = pObs.OrderBy(point => new Point3d(point.X, point.Y, 0).DistanceTo(gridOriginZ0)).ToList();

        // Pannelli intersecati da un ostacolo
        var tmp = Intersection.CurveCurve(pl.ToPolylineCurve(), obs[i].ToPolylineCurve(), 0.1, 0.1).Count;
        counter += tmp;

        // Pannelli C circondati
        bool panelIn = new Point3d(pObs[1].X, pObs[1].Y, 0).DistanceTo(gridOriginZ0) > new Point3d(pl[1].X, pl[1].Y, 0).DistanceTo(gridOriginZ0)
          && new Point3d(pObs[0].X, pObs[0].Y, 0).DistanceTo(gridOriginZ0) < new Point3d(pl[0].X, pl[0].Y, 0).DistanceTo(gridOriginZ0)
          && zObs[1] > z[1] && z[0] > zObs[0];
        if (panelIn)
        {
          counter++;
          crossed = new bool[] { true, true };
        }

        // Pannelli C con un foro
        bool panelOut = new Point3d(pObs[1].X, pObs[1].Y, 0).DistanceTo(gridOriginZ0) < new Point3d(pl[1].X, pl[1].Y, 0).DistanceTo(gridOriginZ0)
          && new Point3d(pObs[0].X, pObs[0].Y, 0).DistanceTo(gridOriginZ0) > new Point3d(pl[0].X, pl[0].Y, 0).DistanceTo(gridOriginZ0)
          && zObs[1] < z[1] && z[0] < zObs[0];
        if (panelOut) counter++;
      }

      if (counter > 0)
      {
        type += "*C";
      }
    }

    public List<int> InterceptedPanelWall(List<Polyline> obs)
    {
      List<int> result = new List<int>();

      foreach (Polyline polyline in obs)
      {
        int i = obs.IndexOf(polyline);
        result.Add(-1);

        List<Point3d> pObs = polyline.Select(point => point).ToList();
        pObs.RemoveAt(pObs.Count - 1);
        pObs = pObs.OrderBy(point => point.Z).ToList();
        double[] zObs = new double[] { pObs[0].Z, pObs[3].Z };
        pObs.RemoveRange(2, 2);
        pObs = pObs.OrderBy(point => new Point3d(point.X, point.Y, 0).DistanceTo(gridOriginZ0)).ToList();

        bool horizontalIn = new Point3d(pObs[1].X, pObs[1].Y, 0).DistanceTo(gridOriginZ0) >= new Point3d(pl[1].X, pl[1].Y, 0).DistanceTo(gridOriginZ0)
  && new Point3d(pObs[0].X, pObs[0].Y, 0).DistanceTo(gridOriginZ0) <= new Point3d(pl[0].X, pl[0].Y, 0).DistanceTo(gridOriginZ0);
        bool verticalIn = zObs[1] >= pl[2].Z && zObs[0] <= pl[0].Z;

        bool verticalIntersection = (pl[2].Z > zObs[1] && zObs[1] > pl[0].Z) | (pl[2].Z > zObs[0] && zObs[0] > pl[0].Z);

        bool a = new Point3d(pl[1].X, pl[1].Y, 0).DistanceTo(gridOriginZ0) > new Point3d(pObs[1].X, pObs[1].Y, 0).DistanceTo(gridOriginZ0);
        bool b = new Point3d(pObs[1].X, pObs[1].Y, 0).DistanceTo(gridOriginZ0) > new Point3d(pl[0].X, pl[0].Y, 0).DistanceTo(gridOriginZ0);
        bool c = new Point3d(pl[1].X, pl[1].Y, 0).DistanceTo(gridOriginZ0) > new Point3d(pObs[0].X, pObs[0].Y, 0).DistanceTo(gridOriginZ0);
        bool d = new Point3d(pObs[0].X, pObs[0].Y, 0).DistanceTo(gridOriginZ0) > new Point3d(pl[0].X, pl[0].Y, 0).DistanceTo(gridOriginZ0);
        bool horizontalIntersection = (a && b) || (c && d);

        if (horizontalIn && verticalIn)
        {
          result[i] = 0;
          crossed = new bool[] { true, true };
        }
        else if (verticalIntersection && horizontalIn)
        {
          result[i] = 1;
        }
        else if (horizontalIntersection && verticalIn)
        {
          if (result[i] == -1)
          {
            result[i] = 2;
          }
        }
        else if (horizontalIntersection && verticalIntersection)
        {
          result[i] = 3;
          type += "*B*C";
        }
      }

      return result;
    }

    public bool GenericInterception(Polyline obs)
    {
      bool ret = false;

      double[] z = new double[] { pl[1].Z, pl[2].Z };
      int counter = 0;

      List<Point3d> pObs = new List<Point3d> { obs[0], obs[1], obs[2], obs[3] };
      pObs = pObs.OrderBy(point => point.Z).ToList();
      double[] zObs = pObs.GetRange(1, 2).Select(point => point.Z).ToArray();

      pObs.RemoveRange(2, 2);
      pObs = pObs.OrderBy(point => new Point3d(point.X, point.Y, 0).DistanceTo(gridOriginZ0)).ToList();

      // Pannelli circondati
      bool panelIn = new Point3d(pObs[1].X, pObs[1].Y, 0).DistanceTo(gridOriginZ0) > new Point3d(pl[1].X, pl[1].Y, 0).DistanceTo(gridOriginZ0)
        && new Point3d(pObs[0].X, pObs[0].Y, 0).DistanceTo(gridOriginZ0) < new Point3d(pl[0].X, pl[0].Y, 0).DistanceTo(gridOriginZ0)
        && zObs[1] > z[1] && z[0] > zObs[0];

      if (panelIn)
      {
        counter++;
        ret = true;
      }

      return ret;
    }

    public void Detached(List<Polyline> det)
    {
      // Per ogni Line di intersezione
      int counter = 0;

      for (int i = 0; i < det.Count; i++)
      {
        var tmp = Intersection.CurveCurve(pl.ToPolylineCurve(), det[i].ToPolylineCurve(), 0.1, 0.1).Count;
        counter += tmp;
      }

      if (counter != 0)
      {
        type += "*B";
      }
    }

  }

  public List<string> ArchiveTypesOld(List<PanelC41> panelC41s)
  {
    List<string> list;
    var orderedPanels = panelC41s.OrderBy(panelC41 => panelC41.type);

    list = orderedPanels.Select(x => x.type.Split('.')[0]).ToList();

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

  public List<string> ArchiveTypes(List<PanelC41> panelC41s)
  {
    List<string> list = panelC41s.Select(x => Regex.Replace(x.type, @"[\d.]", String.Empty)).ToList();

    list = list.Distinct().ToList();

    return list;
  }
  #endregion
}