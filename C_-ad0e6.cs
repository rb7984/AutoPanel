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
  private void RunScript(DataTree<object> x, DataTree<object> y, ref object A, ref object pl, ref object layerNames, ref object names, ref object archive, ref object freeTag, ref object export, ref object PanelC41)
  {
    facade = new Facade(x, y);

    A = facade.gridCount;
    //pl = facade.grids[0].Plines(panelC41s);
    //layerNames = facade.grids[0].TypeTags(panelC41s);
    //names = facade.grids[0].NameTags(panelC41s);
    //archive = ArchiveTypes(panelC41s);
    //freeTag = facade.grids[0].FreeTag(panelC41s);
    //export = facade.grids[0].toExport(panelC41s);
  }
  #endregion
  #region Additional

  // fields
  public Facade facade;
  public List<PanelC41> panelC41s;

  public List<string> archive;

  public class Facade
  {
    //fields
    public List<Grid3d> grids;
    public DataTree<PanelC41> panelsTree;
    public List<int> gridCount;

    //contructor
    public Facade(DataTree<object> input, DataTree<object> obs)
    {
      grids = new List<Grid3d>();
      panelsTree = new DataTree<PanelC41>();
      gridCount = new List<int>();

      for (int i = 1; i < input.Branch(0).Count; i++)
      {
        int count = 0;
        DataTree<object> tmp = new DataTree<object>();
        tmp.AddRange(new List<object> { input.Branch(0)[0], input.Branch(0)[i] }); count++;
        tmp.AddRange(new List<object> { input.Branch(1)[0], input.Branch(1)[i] }, new GH_Path(count++)); count++;

        for (int j = 2; j < input.BranchCount - 1; j++)
        {
          var a = (int)input.Branch(j)[0].ToString().Split('-')[1][1];
          if (a == i - 1) { tmp.AddRange(input.Branch(j), new GH_Path(count++)); count++; }
        }
        tmp.AddRange(new List<object> { input.Branch(input.BranchCount - 1)[0], input.Branch(input.BranchCount - 1)[i] }, new GH_Path(count++));

        DataTree<object> tmpObs = new DataTree<object>();

        for (int j = i - 1; j < obs.BranchCount; j = j + input.Branch(0).Count - 1)
        {
          tmpObs.AddRange(obs.Branch(j), new GH_Path());
        }

        Grid3d tmpGrid = new Grid3d(tmp, tmpObs);
        //grids.Add(tmpGrid);
        //panelsTree.AddRange(tmpGrid.panels, new GH_Path(i - 1));
        gridCount.Add(0);
      }
    }

    public double FacadeCount(DataTree<object> input)
    {
      double a = char.GetNumericValue(input.Branch(input.BranchCount - 1)[0].ToString().Split('-')[1][0]);

      return a;
    }
  }

  public class Grid
  {
    public DataTree<object> param;
    public DataTree<Polyline> obs;

    // tutte le coordinate delle fughe in X e Y
    public List<double> xCoordinates;
    public List<double> yCoordinates;

    public List<PanelC41> panels;
    public List<Point3d> pts;

    public Grid(DataTree<object> x, DataTree<object> y)
    {
      param = x;
      obs = new DataTree<Polyline>();

      for (int i = 0; i < y.BranchCount; i++)
      { obs.AddRange(y.Branch(i).Select(pl => (PolylineCurve)pl).ToList().Select(pl => pl.ToPolyline()).ToList(), new GH_Path(i)); }

      OrderLines(param);

      pts = new List<Point3d>();

      panels = Panels(xCoordinates, yCoordinates);

      OrderOutput(panels);

      PostPanels(panels);
    }

    public void OrderLines(DataTree<object> param)
    {
      List<double> xtmp = new List<double>();
      List<double> ytmp = new List<double>();

      // param.BranchCount - 1 perchè l'ultimo branch sono i bordi esterni
      for (int i = 2; i < param.BranchCount - 1; i++)
      {
        int fuga = (int)Convert.ToInt64(param.Branch(i)[0].ToString().Split('-')[2].ToString().Split('.')[0]);
        string dir = param.Branch(i)[0].ToString().Split('-')[1];

        List<double> tmp = new List<double>();

        if (dir == "h")
        {
          // start at j = 1 perchè la prima linea è il nome del layer 
          for (int j = 1; j < param.Branch(i).Count; j++)
          {
            PolylineCurve tmpPl = (PolylineCurve)param.Branch(i)[j];
            tmp.Add(tmpPl.PointAt(0).Y - fuga);
            tmp.Add(tmpPl.PointAt(0).Y + fuga);
          }

          ytmp.AddRange(tmp);
        }

        else if (dir == "v")
        {
          // start at j = 1 perchè la prima linea è il nome del layer
          for (int j = 1; j < param.Branch(i).Count; j++)
          {
            PolylineCurve tmpPl = (PolylineCurve)param.Branch(i)[j];
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

      Borders(param);
    }

    public void OrderOutput(List<PanelC41> panels)
    {
      List<PanelC41> tmp = panels.OrderBy(a => a.firstCorner[0]).ThenBy(b => b.firstCorner[1]).ToList();

      this.panels.Clear();

      this.panels.AddRange(tmp);
    }

    void PostPanels(List<PanelC41> panels)
    {
      for (int i = 0; i < panels.Count - 1; i++)
      {
        if (panels[i].type == "B" || panels[i].type == "B*C")
        {
          if (panels[i + 1].type == "C") panels[i + 1].type = "D*C";
          else panels[i + 1].type = "D";
        }
      }
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

          PanelC41 tmpPanel = new PanelC41(points, obs);
          if (!tmpPanel.crossed) panels.Add(tmpPanel);
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
      List<string> types = panels.Select(i => i.type).ToList();

      return types;
    }

    public List<string> toExport(List<PanelC41> panels)
    {
      List<string> export = panels.Select(i => i.toExcel).ToList();
      export.Insert(0, "Type,Width,Heigh,Marca");

      return export;
    }

    void Borders(DataTree<object> param)
    {
      int[] fughe = param.Branch(param.BranchCount - 1)[0].ToString().Split('-')[2].Split('.').Select(element => Convert.ToInt32(element)).ToArray();

      PolylineCurve curve = (PolylineCurve)param.Branch(param.BranchCount - 1)[1];

      double[] x = new double[4];
      double[] y = new double[4];

      for (int i = 0; i < 4; i++)
      {
        x[i] = curve.Point(i).X;
        y[i] = curve.Point(i).Y;
      }

      yCoordinates.Insert(0, Math.Round(y.Min() + fughe[0]));
      xCoordinates.Insert(xCoordinates.Count, Math.Round(x.Max() - fughe[1]));
      yCoordinates.Insert(yCoordinates.Count, Math.Round(y.Max() - fughe[2]));
      xCoordinates.Insert(0, Math.Round(x.Min() + fughe[3]));
    }
  }

  public class Grid3d
  {
    public DataTree<object> param;
    public DataTree<Polyline> obs;
    double[] rotation;

    // tutte le coordinate delle fughe in X e Z
    public List<double> xCoordinates;
    public List<double> zCoordinates;

    public List<PanelC41> panels;
    public List<Point3d> pts;

    public Grid3d(DataTree<object> x, DataTree<object> y)
    {
      param = x;
      obs = new DataTree<Polyline>();

      // [rad, X,Y,Z]
      rotation = Rotation(x.Branch(0)[1]);

      for (int i = 0; i < y.BranchCount; i++)
      { obs.AddRange(y.Branch(i).Select(pl => (PolylineCurve)pl).ToList().Select(pl => pl.ToPolyline()).ToList(), new GH_Path(i)); }

      OrderLines(param);

      //pts = new List<Point3d>();

      //panels = Panels(xCoordinates, zCoordinates);



      //OrderOutput(panels);

      //PostPanels(panels);
    }

    public double[] Rotation(object a)
    {
      PolylineCurve tmpPl = (PolylineCurve)a;

      var vec = tmpPl.PointAtEnd - tmpPl.PointAtStart;

      double angle = Vector3d.VectorAngle(vec, new Vector3d(1, 0, 0));

      double[] coordinates = new double[]
      {
        angle,
        tmpPl.PointAtStart.X,
        tmpPl.PointAtStart.Y,
        tmpPl.PointAtStart.Z
      };

      return coordinates;
    }

    public void OrderLines(DataTree<object> param)
    {
      List<double> xtmp = new List<double>();
      List<double> ztmp = new List<double>();

      // i = 2 perché i primi due sono i due bordi singoli
      // param.BranchCount - 1 perchè l'ultimo branch sono i bordi esterni
      for (int i = 2; i < param.BranchCount - 1; i++)
      {
        int fuga = (int)Convert.ToInt64(param.Branch(i)[0].ToString().Split('-')[2].ToString().Split('.')[0]);
        char dir = param.Branch(i)[0].ToString().Split('-')[1][0];

        List<double> tmp = new List<double>();

        if (dir == 'h')
        {
          // start at j = 1 perchè la prima linea è il nome del layer 
          for (int j = 1; j < param.Branch(i).Count; j++)
          {
            PolylineCurve tmpPl = (PolylineCurve)param.Branch(i)[j];
            tmp.Add(tmpPl.PointAt(0).Z - fuga);
            tmp.Add(tmpPl.PointAt(0).Z + fuga);
          }

          ztmp.AddRange(tmp);
        }

        else if (dir == 'v')
        {
          // start at j = 1 perchè la prima linea è il nome del layer
          for (int j = 1; j < param.Branch(i).Count; j++)
          {
            PolylineCurve tmpPl = (PolylineCurve)param.Branch(i)[j];
            Point3d tmpPoint = tmpPl.PointAtStart;
            tmpPoint.Transform(Transform.Rotation(rotation[0], new Point3d(rotation[1], rotation[2], 0)));
            tmp.Add(tmpPoint.X - fuga);
            tmp.Add(tmpPoint.X + fuga);
          }

          xtmp.AddRange(tmp);
        }
      }

      xtmp.Sort();
      ztmp.Sort();

      xCoordinates = xtmp;
      zCoordinates = ztmp;

      Borders(param);
    }

    void Borders(DataTree<object> param)
    {
      int[] fughe = param.Branch(param.BranchCount - 1)[0].ToString().Split('-')[2].Split('.').Select(element => Convert.ToInt32(element)).ToArray();
      //int a = Convert.ToInt32(param.Branch(param.BranchCount - 1)[0].ToString().Split('-')[2].Split('.')[1]);
      //PolylineCurve curve = (PolylineCurve)param.Branch(param.BranchCount - 1)[1];
      //Point3d tmpPoint = curve.PointAtStart;
      //tmpPoint.Transform(Transform.Rotation(rotation[0], new Point3d(rotation[1], rotation[2], 0)));

      //double[] x = new double[4];
      //double[] z = new double[4];

      //for (int i = 0; i < 4; i++)
      //{
      //  x[i] = curve.Point(i).X;
      //  z[i] = curve.Point(i).Z;
      //}

      //zCoordinates.Insert(0, Math.Round(z.Min() + fughe[0]));
      //xCoordinates.Insert(xCoordinates.Count, Math.Round(x.Max() - fughe[1]));
      //zCoordinates.Insert(zCoordinates.Count, Math.Round(z.Max() - fughe[2]));
      //xCoordinates.Insert(0, Math.Round(x.Min() + fughe[3]));
    }

    public void OrderOutput(List<PanelC41> panels)
    {
      List<PanelC41> tmp = panels.OrderBy(a => a.firstCorner[0]).ThenBy(b => b.firstCorner[1]).ToList();

      this.panels.Clear();

      this.panels.AddRange(tmp);
    }

    void PostPanels(List<PanelC41> panels)
    {
      for (int i = 0; i < panels.Count - 1; i++)
      {
        if (panels[i].type == "B" || panels[i].type == "B*C")
        {
          if (panels[i + 1].type == "C") panels[i + 1].type = "D*C";
          else panels[i + 1].type = "D";
        }
      }
    }

    public List<PanelC41> Panels(List<double> xCoordinates, List<double> zCoordinates)
    {
      List<PanelC41> panels = new List<PanelC41>();

      for (int i = 0; i < zCoordinates.Count; i += 2)
      {
        for (int j = 0; j < xCoordinates.Count; j += 2)
        {
          var points = new List<Point3d>
            {
              new Point3d(xCoordinates[j], 0, zCoordinates[i]),
              new Point3d(xCoordinates[j + 1], 0, zCoordinates[i]),
              new Point3d(xCoordinates[j + 1], 0, zCoordinates[i + 1]),
              new Point3d(xCoordinates[j], 0, zCoordinates[i + 1]),
              new Point3d(xCoordinates[j], 0, zCoordinates[i])
              };

          pts.AddRange(points.GetRange(0, 4));

          PanelC41 tmpPanel = new PanelC41(points, obs);
          if (!tmpPanel.crossed) panels.Add(tmpPanel);
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
      List<string> types = panels.Select(i => i.type).ToList();

      return types;
    }

    public List<string> toExport(List<PanelC41> panels)
    {
      List<string> export = panels.Select(i => i.toExcel).ToList();
      export.Insert(0, "Type,Width,Heigh,Marca");

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
    public double[] firstCorner = new double[2];

    public string toExcel;

    //Constructor
    public PanelC41(List<Point3d> list, DataTree<Polyline> obs)
    {
      borderUp = 1;
      borderDown = 2;
      hookUp = 3;
      hookDown = 4;
      type = "A";

      crossed = false;

      pl = new Polyline(list);
      firstCorner[0] = pl[0].X;
      firstCorner[1] = pl[0].Y;

      int[] i = Intercept(obs);

      if (i[0] >= 0)
      {
        crossed = InterceptedPanelVertical(obs.Branch(0), i[0], 8);
      }
      if (i[1] >= 0)
      {
        crossed = InterceptedPanelHorizontal(obs.Branch(2), i[1], 8);
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

    public int[] Intercept(DataTree<Polyline> obs)
    {
      // ora non funziona per ostacoli multipli che intersecano oggetti

      int counter = 0;
      int b = -1;
      int c = -1;

      for (int i = 0; i < obs.Branch(0).Count; i++)
      {
        var tmp = Intersection.CurveCurve(pl.ToPolylineCurve(), obs.Branch(0)[i].ToPolylineCurve(), 0.1, 0.1).Count;
        counter += tmp;
        if (tmp != 0) b = i;
      }

      for (int i = 0; i < obs.Branch(2).Count; i++)
      {
        var tmp = Intersection.CurveCurve(pl.ToPolylineCurve(), obs.Branch(2)[i].ToPolylineCurve(), 0.1, 0.1).Count;
        counter += tmp;
        if (tmp != 0) c = i;
      }

      if (counter == 0) crossed = false;
      else crossed = true;
      return new int[] { b, c };
    }

    public bool InterceptedPanelVertical(List<Polyline> obs, int i, double fuga)
    {
      // Funziona solo per finestre centrate sui pannelli
      // Viene considerata solo la Y quindi

      double[] yObs = new double[] { obs[i][0].Y, obs[i][1].Y, obs[i][2].Y };
      double[] y = new double[] { pl[0].Y, pl[1].Y, pl[2].Y };

      if (yObs.Max() >= y.Max() && yObs.Min() <= y.Min())
      {
        return true;
      }
      else if (y.Max() > yObs.Max())
      {
        double newY = yObs.Max();
        List<Point3d> tmp = new List<Point3d>()
        {
          new Point3d(pl[0].X, newY + fuga,0),
          new Point3d(pl[1].X, newY+ fuga,0),
          new Point3d(pl[2]),
          new Point3d(pl[3]),
          new Point3d(pl[4].X, newY+ fuga,0),
        };

        pl = new Polyline(tmp);
        return false;
      }
      else if (y.Min() < yObs.Min())
      {
        double newY = yObs.Min();
        List<Point3d> tmp = new List<Point3d>()
        {
          new Point3d(pl[0]),
          new Point3d(pl[1]),
          new Point3d(pl[2].X, newY - fuga,0),
          new Point3d(pl[3].X, newY - fuga,0),
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
      // Viene considerata solo la X quindi

      double[] xObs = new double[] { obs[i][0].X, obs[i][1].X, obs[i][2].X };
      double[] x = new double[] { pl[0].X, pl[1].X, pl[2].X };

      if (xObs.Max() >= x.Max() && xObs.Min() <= x.Min())
      {
        return true;
      }

      else if (x.Max() > xObs.Max())
      {
        double newX = xObs.Max();
        List<Point3d> tmp = new List<Point3d>()
        {
          new Point3d(newX + fuga, pl[0].Y, 0),
          new Point3d(pl[1]),
          new Point3d(pl[2]),
          new Point3d(newX + fuga, pl[3].Y,0 ),
          new Point3d(newX + fuga, pl[0].Y, 0)
        };

        pl = new Polyline(tmp);
        return false;
      }

      else if (x.Min() < xObs.Min())
      {
        double newX = xObs.Min();
        List<Point3d> tmp = new List<Point3d>()
        {
          new Point3d(pl[0]),
          new Point3d(newX - fuga, pl[1].Y, 0),
          new Point3d(newX - fuga, pl[2].Y, 0),
          new Point3d(pl[3]),
          new Point3d(pl[0])
        };

        pl = new Polyline(tmp);

        return false;
      }

      else { return false; }
    }

    public void InterceptedPanelVarious(List<Polyline> obs)
    {
      int counter = 0;

      double[] x = new double[] { pl[0].X, pl[1].X };
      double[] y = new double[] { pl[1].Y, pl[2].Y };

      for (int i = 0; i < obs.Count; i++)
      {
        // Pannelli intersecati da un ostacolo
        var tmp = Intersection.CurveCurve(pl.ToPolylineCurve(), obs[i].ToPolylineCurve(), 0.1, 0.1).Count;
        counter += tmp;

        double[] xObs = new double[] { obs[i][0].X, obs[i][1].X };
        double[] yObs = new double[] { obs[i][1].Y, obs[i][2].Y };

        // Pannelli C circondati
        if (xObs[1] > x[1] && x[0] > xObs[0] && yObs[1] > y[1] && y[0] > yObs[0])
        {
          counter++;
          crossed = true;
        }

        // Pannelli C con un foro
        if (xObs[1] < x[1] && x[0] < xObs[0] && yObs[1] < y[1] && y[0] < yObs[0]) counter++;
      }

      if (counter > 0)
      {
        if (type == "A") type = "C";
        else type += "*C";
      }
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

    list = orderedPanels.Select(x => x.type).ToList();

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