using System;
using System.Collections;
using System.Collections.Generic;

using Rhino;
using Rhino.Geometry;

using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;

using Rhino.Commands;
using Rhino.Display;
using System.Drawing;
using Rhino.DocObjects;


/// <summary>
/// This class will be instantiated on demand by the Script component.
/// </summary>
public abstract class Script_Instance_ba9bf : GH_ScriptInstance
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
  private void RunScript(bool bake, DataTree<GeometryBase> plines, DataTree<GeometryBase> hatches, DataTree<string> names, DataTree<Plane> textLocations, string savingPath, List<Polyline> baseLines, ref object A)
  {
    RhinoDoc doc = RhinoDoc.ActiveDoc;

    // TEXT ATTRIBUTE
    ObjectAttributes textAttribute = new ObjectAttributes();
    textAttribute.LayerIndex = doc.Layers.FindName("text").Index;

    // GREY HATCH ATTRIBUTE
    ObjectAttributes greyHatchesAtt = new ObjectAttributes();
    greyHatchesAtt.LayerIndex = doc.Layers.FindName("hatch_grey").Index;
    Color colorGray = Color.Gray;
    greyHatchesAtt.ObjectColor = colorGray;
    greyHatchesAtt.ColorSource = ObjectColorSource.ColorFromObject;

    if (bake)
    {
      for (int i = 0; i < plines.BranchCount; i++)
      {
        ConstructionPlane cp = new ConstructionPlane();
        Plane p = new Plane(baseLines[i][0], new Vector3d(baseLines[i][1] - baseLines[i][0]), new Vector3d(0, 0, 1));

        cp.Plane = p;

        doc.Views.ActiveView.ActiveViewport.SetConstructionPlane(cp);

        for (int j = 0; j < plines.Branch(i).Count; j++)
        {
          ObjectAttributes colorHatchesAtt = new ObjectAttributes();
          colorHatchesAtt.LayerIndex = doc.Layers.FindName("hatch_colors").Index;
          int[] tmp = ColorRB(names.Branch(i)[j]);
          Color color = Color.FromArgb(tmp[0], tmp[1], tmp[2]);
          colorHatchesAtt.ObjectColor = color;
          colorHatchesAtt.ColorSource = ObjectColorSource.ColorFromObject;

          ObjectAttributes plineAtt = new ObjectAttributes();
          plineAtt.LayerIndex = doc.Layers.FindName(names.Branch(i)[j]).Index;

          TextEntity t = new TextEntity
          {
            PlainText = names.Branch(i)[j],
            Plane = textLocations.Branch(i)[j]
          };

          RhinoDocument.Objects.Add(plines.Branch(i)[j], plineAtt);
          RhinoDocument.Objects.Add(t, textAttribute);
          RhinoDocument.Objects.Add(hatches.Branch(i)[j], greyHatchesAtt);
          RhinoDocument.Objects.Add(hatches.Branch(i)[j], colorHatchesAtt);
        }

        RhinoApp.RunScript("-SelAll", false);
        string path = savingPath + "/" + i.ToString() + ".dwg";

        RhinoApp.RunScript("-Export " + path + " " + "Enter ", false);
        RhinoApp.RunScript("-SelAll", false);
        RhinoApp.RunScript("-Delete", false);
      }
    }
  }
  #endregion
  #region Additional

  public int[] ColorRB(string layerName)
  {
    int[] k;

    if (layerName == "A")
    {
      k = new int[] { 242, 68, 68, 255 };
      return k;
    }
    else if (layerName == "B")
    {
      k = new int[] { 82, 137, 227, 255 };
      return k;
    }
    else if (layerName == "D")
    {
      k = new int[] { 150, 92, 250, 255 };
      return k;
    }
    else if (layerName == "E")
    {
      k = new int[] { 77, 200, 219, 255 };
      return k;
    }
    else if (layerName == "F")
    {
      k = new int[] { 174, 233, 242, 255 };
      return k;
    }
    else if (layerName == "G")
    {
      k = new int[] { 52, 140, 153, 255 };
      return k;
    }
    else if (layerName == "H")
    {
      k = new int[] { 10, 169, 194, 255 };
      return k;
    }
    else
    {
      k = new int[] { 227, 208, 61, 255 };
      return k;
    }
  }
  #endregion
}