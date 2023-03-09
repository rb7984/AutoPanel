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
  private void RunScript(DataTree<object> x, DataTree<object> y, List<Polyline> z, ref object A, ref object B, ref object pl, ref object layerNames, ref object names, ref object archive, ref object freeTag, ref object export, ref object PanelC41)
  {
    facade = new Facade(x, y, z);

    panelC41s = new List<PanelC41>();
    plines = new DataTree<Polyline>();
    this.archive = new List<string>();
    this.names = new DataTree<string>();
    types = new DataTree<string>();
    toExport = new List<string> { "Type,Width,Heigh,Marca,Facciata" };
    normals = new List<Vector3d>();

    foreach (Grid3d grid in facade.grids)
    {
      int i = facade.grids.IndexOf(grid);
      panelC41s.AddRange(grid.panels);
      plines.AddRange(grid.Plines(grid.panels), new GH_Path(i));
      this.names.AddRange(grid.NameTags(grid.panels), new GH_Path(i));
      types.AddRange(grid.TypeTags(grid.panels), new GH_Path(i));
      toExport.AddRange(grid.ToExport(grid.panels));
      normals.Add(grid.normal);
    }

    this.archive.AddRange(ArchiveTypes(panelC41s));

    A = facade.angles;
    B = facade.baseLines;
    pl = plines;
    layerNames = types;
    names = this.names;
    archive = this.archive;
    export = toExport;
    PanelC41 = normals;
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

    //contructor
    public Facade(DataTree<object> input, DataTree<object> obs, List<Polyline> baseline)
    {
      grids = new List<Grid3d>();

      // counter = [3, 3, 3]
      var counter = FacadeCount(input);
      var pointer = 0;

      FillAngles(baseline);

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

      for (int i = 0; i < grids.Count; i++)
      {
        PostPanels(grids[i]);
      }

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

    public void PostPanels(Grid3d grid)
    {
      int fugaMax = 16;

      for (int i = 0; i < grid.panels.Count; i++)
      {
        double l = new Point3d(grid.panels[i].pl[0].X, grid.panels[i].pl[0].Y, 0).DistanceTo(new Point3d(grid.border[0].X, grid.border[0].Y, 0));
        double r = new Point3d(grid.panels[i].pl[1].X, grid.panels[i].pl[1].Y, 0).DistanceTo(new Point3d(grid.border[1].X, grid.border[1].Y, 0));

        //LEFT
        if (Math.Abs(l) < fugaMax)
        {
          if (grid.panels[i].type == "A*C") grid.panels[i].type = "C*" + BorderPanel(grids.IndexOf(grid), 0);
          else grid.panels[i].type = BorderPanel(grids.IndexOf(grid), 0);
        }
        //RIGHT
        if (Math.Abs(r) < fugaMax)
        {
          if (grid.panels[i].type == "A*C") grid.panels[i].type = "C*" + BorderPanel(grids.IndexOf(grid), 1);
          else grid.panels[i].type = BorderPanel(grids.IndexOf(grid), 1);
        }
      }
    }

    public string BorderPanel(int i, int a)
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

    public void CorrectExport()
    {
      foreach (Grid3d grid in grids)
      {
        int a = grids.IndexOf(grid);
        foreach (PanelC41 p in grid.panels)
        {
          p.toExcel = p.type + "," + p.width.ToString() + "," + p.height.ToString() +
        "," + p.type + "." + p.width.ToString() + "." + p.height.ToString() + "," + a.ToString();
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

      OrderOutput(panels);

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
        int fuga = (int)Convert.ToInt64(param.Branch(i)[0].ToString().Split('-')[2].ToString().Split('.')[0]);
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
            var events = Intersection.CurveCurve(border.ToNurbsCurve(), c.ToNurbsCurve(), 0.1, 0.1);
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
          if (!tmpPanel.crossed) panels.Add(tmpPanel);
        }
      }

      return panels;
    }

    public void OrderOutput(List<PanelC41> panels)
    {
      List<PanelC41> tmp;

      tmp = panels.OrderBy(a => new Point3d(a.firstCorner.X, a.firstCorner.Y, 0).DistanceTo(border[0])).ThenBy(b => b.firstCorner.Z).ToList();

      this.panels.Clear();

      this.panels.AddRange(tmp);
    }

    void PostPanels(List<PanelC41> panels)
    {
      for (int i = 0; i < panels.Count - 1; i++)
      {
        if (panels[i].type == "B" || panels[i].type == "B*C")
        {
          if (panels[i + 1].type == "C" || panels[i + 1].type == "A*C") panels[i + 1].type = "D*C";
          else panels[i + 1].type = "D";
        }
      }
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

    public List<string> DimensionsTags(List<PanelC41> panels)
    {
      List<string> nameTags = panels.Select(i => i.name).ToList();

      return nameTags;
    }

    public List<string> NameTags(List<PanelC41> panels)
    {
      List<string> names = panels.Select(i => i.type).ToList();
      for (int i = 0; i < panels.Count; i++)
      {
        names[i] = names[i] + "|" + panels[i].name;
      }

      return names;
    }

    public List<string> TypeTags(List<PanelC41> panels)
    {
      List<string> types = panels.Select(i => i.type.Split('.')[0]).ToList();

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
    public double borderUp;
    public double borderDown;
    public double hookUp;
    public double hookDown;

    public Polyline pl;
    public bool crossed;
    public string type;
    public string name;
    public Point3d firstCorner;
    public Vector3d ax;
    public Point3d gridOrigin;
    public Point3d gridOriginZ0;

    public string toExcel;

    //constructor
    public PanelC41(List<Point3d> list, DataTree<Polyline> obs, Polyline border)
    {
      borderUp = 1;
      borderDown = 2;
      hookUp = 3;
      hookDown = 4;
      type = "A";

      crossed = false;

      pl = new Polyline(list);

      firstCorner = pl[0];

      ax = new Vector3d(border[1] - border[0]);
      gridOrigin = border[0];
      gridOriginZ0 = new Point3d(border[0].X, border[0].Y, 0);

      int[] wall = InterceptedPanelWall(obs.Branch(4));

      // [0 = internal,
      // 1 = verticalInterception,
      // 2 = horizontalInteception,
      // 3 = corner]

      if (wall[0] != 0)
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
            crossed = InterceptedPanelHorizontal(obs.Branch(2), i[j], 8);
          }
        }
        if (wall[0] == 1) crossed = InterceptedPanelVertical(obs.Branch(4), wall[1], 8);
        if (wall[0] == 2) crossed = InterceptedPanelHorizontal(obs.Branch(4), wall[1], 8);
      }

      Name();
      Detached(obs.Branch(1));
      InterceptedPanelVarious(obs.Branch(3));

      toExcel = type + "," + width.ToString() + "," + height.ToString() +
        "," + type + "." + width.ToString() + "." + height.ToString();
    }

    public PanelC41(List<Point3d> list)
    {
      borderUp = 1;
      borderDown = 2;
      hookUp = 3;
      hookDown = 4;
      type = "A";

      crossed = false;

      pl = new Polyline(list);

      firstCorner = pl[0];

      Name();

      toExcel = type + "," + width.ToString() + "," + height.ToString() +
        "," + type + "." + width.ToString() + "." + height.ToString();
    }

    public List<int> Intercept(DataTree<Polyline> obs)
    {
      List<int> intercepted = new List<int>();
      intercepted.Add(0);

      // Branch(0) windows. Polylines centrate sul pannello verticalmente
      for (int i = 0; i < obs.Branch(0).Count; i++)
      {
        var tmp = Intersection.CurveCurve(pl.ToPolylineCurve(), obs.Branch(0)[i].ToPolylineCurve(), 0.1, 0.1).Count;

        if (tmp != 0)
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

      if (intercepted[0] == 0 && intercepted[intercepted[0] + 1] == 0) crossed = false;
      else crossed = true;
      return intercepted;
    }

    public void Name()
    {
      width = Math.Round(pl[0].DistanceTo(pl[1]));
      height = Math.Round(pl[1].DistanceTo(pl[2]));
      name = width.ToString() + '-' + height.ToString();
    }

    public bool InterceptedPanelVertical(List<Polyline> obs, int i, double fuga)
    {
      // Funziona solo per finestre centrate sui pannelli
      // Viene considerata solo la Z quindi

      double[] zObs = new double[] { obs[i][0].Z, obs[i][1].Z, obs[i][2].Z };
      double[] z = new double[] { pl[0].Z, pl[1].Z, pl[2].Z };

      if (zObs.Max() >= z.Max() && zObs.Min() <= z.Min())
      {
        return true;
      }
      else if (z.Max() > zObs.Max())
      {
        double newZ = zObs.Max();
        List<Point3d> tmp = new List<Point3d>()
        {
          new Point3d(pl[0].X, pl[0].Y, newZ+ fuga),
          new Point3d(pl[1].X, pl[1].Y, newZ+ fuga),
          new Point3d(pl[2]),
          new Point3d(pl[3]),
          new Point3d(pl[4].X, pl[4].Y, newZ+ fuga),
        };

        pl = new Polyline(tmp);
        return false;
      }
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

        return false;
      }

      else { return false; }
    }

    public bool InterceptedPanelHorizontal(List<Polyline> obs, int i, double fuga)
    {
      // Funziona solo per rettangoli orizzontali centrati sui pannelli
      // Viene considerata solo distanza X-Y

      List<Point3d> p = new List<Point3d> { pl[0], pl[1] };
      List<Point3d> pObs = new List<Point3d> { obs[i][0], obs[i][1], obs[i][2], obs[i][3] };
      pObs.OrderBy(point => point.Z);
      pObs.RemoveRange(2, 2);

      pObs = pObs.OrderBy(point => new Point3d(point.X, point.Y, 0).DistanceTo(gridOrigin)).ToList();

      if (new Point3d(pObs[0].X, pObs[0].Y, 0).DistanceTo(gridOrigin) > 0.01)
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

        type = "F";

        return false;
      }
      else if (pl[1].DistanceTo(gridOrigin) > pObs[1].DistanceTo(gridOrigin))
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

        type = "E";

        return false;
      }
      else
      {
        return true;
      }
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
          crossed = true;
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

    public int[] InterceptedPanelWall(List<Polyline> obs)
    {
      int[] result = new int[2] { -1, -1 };

      foreach (Polyline polyline in obs)
      {
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
        bool horizontalIntersection = (a && b) | (c && d);

        if (horizontalIn && verticalIn)
        {
          result.SetValue(0, 0);
          crossed = true;
        }
        else if (verticalIntersection && horizontalIn)
        {
          result.SetValue(1, 0);
          result.SetValue(obs.IndexOf(polyline), 1);
        }
        else if (horizontalIntersection && verticalIn)
        {
          if (result[0] == -1)
          {
            result.SetValue(2, 0);
            result.SetValue(obs.IndexOf(polyline), 1);
          }
        }
        else if (horizontalIntersection && verticalIntersection)
        {
          result.SetValue(3, 0);
          type += "*C";
        }
      }

      return result;
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
        type = "B";
      }
    }

    public bool SameSign(double num1, double num2)
    {
      return num1 >= 0 && num2 >= 0 || num1 < 0 && num2 < 0;
    }
  }

  public List<string> ArchiveDimensions(List<PanelC41> panelC41s)
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

  public List<string> ArchiveTypes(List<PanelC41> panelC41s)
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
  #endregion
}