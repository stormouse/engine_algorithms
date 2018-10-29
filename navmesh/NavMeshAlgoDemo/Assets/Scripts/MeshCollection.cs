using System.Collections;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;

public class MeshCollection : MonoBehaviour {

    public GameObject[] targets;
    public Material lineMaterial;

    private void Start()
    {
        PostRenderer.Instance.OnPostRenderCallbacks = DrawTriangles;
    }

    private void OnGUI()
    {
        if (GUILayout.Button("Add Triangles to Navigation Volume"))
        {
            foreach (var target in targets)
            {
                Mesh mesh = target.GetComponent<MeshFilter>().mesh;
                Transform t = target.transform;
                Matrix4x4 modelToWorld = t.localToWorldMatrix;
                Vector3[] triangle = new Vector3[3];
                for (int i = 0; i < mesh.triangles.Length; i += 3)
                {
                    triangle[0] = modelToWorld * mesh.vertices[mesh.triangles[i + 0]]; triangle[0] += t.position;
                    triangle[1] = modelToWorld * mesh.vertices[mesh.triangles[i + 1]]; triangle[1] += t.position;
                    triangle[2] = modelToWorld * mesh.vertices[mesh.triangles[i + 2]]; triangle[2] += t.position;
                    FindObjectOfType<NavigationVolume>().MarkAsWalkable(triangle);
                }
            }
        }

        if(GUILayout.Button("Generate Navmesh Areas"))
        {
            var areas = FindObjectOfType<NavigationVolume>().GenerateNavmeshAreas();
            FindObjectOfType<NavigationVolume>().areas = areas;
            Debug.Log(string.Format("#areas: {0}", areas.Count));
        }

        /*

        if (GUILayout.Button("Generate Simplified Contour"))
        {
            // FindObjectOfType<NavigationVolume>().GenerateContour();
            FindObjectOfType<NavigationVolume>().FindHolesAndWalkables();
        }

        if (GUILayout.Button("Shrink Walkable Contour"))
        {
            FindObjectOfType<NavigationVolume>().ShrinkContour(0.5f);
        }

        if (GUILayout.Button("Clip Walkable Area"))
        {
            FindObjectOfType<NavigationVolume>().TestClipWalkableArea();
        }*/

        if (GUILayout.Button("Do partition"))
        {
            FindObjectOfType<NavigationVolume>().DoTestPartition();
        }

        if (GUILayout.Button("Do triangulation"))
        {
            FindObjectOfType<NavigationVolume>().DoTestTriangulation();
        }

        if (GUILayout.Button("Do pathFinding"))
        {
            FindObjectOfType<NavigationVolume>().DoTestPathFinding();
        }
    }

    private void DrawTriangles()
    {
        foreach (var target in targets)
        {
            Mesh mesh = target.GetComponent<MeshFilter>().mesh;
            Transform t = target.transform;
            Matrix4x4 modelToWorld = t.localToWorldMatrix;

            GL.PushMatrix();
            lineMaterial.SetPass(0);
            GL.MultMatrix(modelToWorld);
            for (int i = 0; i < mesh.triangles.Length; i += 3)
            {
                GL.Begin(GL.LINES);
                GL.Vertex(mesh.vertices[mesh.triangles[i + 0]]);
                GL.Vertex(mesh.vertices[mesh.triangles[i + 1]]);
                GL.Vertex(mesh.vertices[mesh.triangles[i + 2]]);
                GL.Vertex(mesh.vertices[mesh.triangles[i + 0]]);
                GL.End();
            }
            GL.PopMatrix();
        }
    }

}
