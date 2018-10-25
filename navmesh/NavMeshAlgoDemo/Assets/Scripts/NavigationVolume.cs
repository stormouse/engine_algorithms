using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;


public class NavmeshArea
{
    public List<Vector3> contour;
    public List<List<Vector3>> holes;
    public NavmeshArea() { holes = new List<List<Vector3>>(); }
}

public class Partition
{
    public List<Vector3> points;
}

public class NavmeshTrangular
{
    public List<Vector3> points;
}

public class NavigationVolume : MonoBehaviour {


    class VertexNode
    {
        static int latestId = 0;
        public int id;
        public Vector3 v;
        public VertexNode next;
        public VertexNode prev;
        public VertexNode(Vector3 _v)
        {
            id = latestId++;
            v = _v;
            next = prev = null;
        }

        public static void ResetStaticId()
        {
            latestId = 0;
        }
    }

    class VertexWANode
    {
        public static int latestId = 0;

        public enum NodeType { Endpoint, Intersection, };
        public enum NodeDirection { In, Out, };

        public int Id;
        public Vector3 Vertex;
        public NodeType Type;
        public VertexWANode ClipperPrev;
        public VertexWANode ClipperNext;
        public VertexWANode PolygonPrev;
        public VertexWANode PolygonNext;

        // we hold the polygon edges and the intersections belongs to the clipper
        public VertexWANode Next { get { return ClipperNext; } set { ClipperNext = value; } }
        public VertexWANode Prev { get { return ClipperPrev; } set { ClipperPrev = value; } }

        // shouldn't use when Type == Endpoint
        public NodeDirection Direction;

        public VertexWANode(Vector3 v, NodeType type)
        {
            Id = latestId++;
            Vertex = v;
            Type = type;
            ClipperPrev = ClipperNext = PolygonPrev = PolygonNext = null;
        }
    }


    public float gridSize;
    public Bounds bounds;

    private int numRows;
    private int numCols;

    private int[,] data;

    private List<Vector3> m_contour = new List<Vector3>();
    private List<List<Vector3>> m_holes = new List<List<Vector3>>();
    private List<List<Vector3>> walkables = new List<List<Vector3>>();
    public List<NavmeshArea> areas = null;
    public List<Partition> partitions = new List<Partition>();
    public List<NavmeshTrangular> triangulars = new List<NavmeshTrangular>();
    

    private void Start()
    {
        Debug.Log(Vector3.Angle(new Vector3(1, 0, 0), new Vector3(-1, 0, 0)));


        Debug.Log(Vector3.Angle(new Vector3(0, 0, 1), new Vector3(1, 0, 0)));
        numRows = Mathf.CeilToInt(bounds.extents.z * 2 / gridSize);
        numCols = Mathf.CeilToInt(bounds.extents.x * 2 / gridSize);
        bounds.extents = new Vector3(numCols * gridSize / 2, bounds.extents.y, numRows * gridSize / 2);
        data = new int[numRows, numCols];
        for (int i = 0; i < numRows; i++) for (int j = 0; j < numCols; j++) data[i, j] = 0;
    }


    public void MarkAsWalkable(Vector3[] triangle)
    {
        Vector3[] flatTriangle = new Vector3[triangle.Length];
        for(int i = 0; i < triangle.Length; i++)
        {
            flatTriangle[i] = Vector3.Scale(triangle[i], new Vector3(1, 0, 1));
        }
        
        for(int i = 0; i < numRows; i++)
        {
            for(int j = 0; j < numCols; j++)
            {
                Vector3[] verts = new Vector3[4];
                Vector3 center = bounds.min + Vector3.right * j * gridSize + Vector3.forward * i * gridSize + Vector3.up * bounds.extents.y;
                verts[0] = new Vector3(center.x - 0.5f * gridSize, 0, center.z - 0.5f * gridSize);
                verts[1] = new Vector3(center.x + 0.5f * gridSize, 0, center.z - 0.5f * gridSize);
                verts[2] = new Vector3(center.x + 0.5f * gridSize, 0, center.z + 0.5f * gridSize);
                verts[3] = new Vector3(center.x - 0.5f * gridSize, 0, center.z + 0.5f * gridSize);
                bool overlapTest = false;
                for(int k = 0; k < 4; k++)
                {
                    float s0 = Vector3.Cross(flatTriangle[1] - flatTriangle[0], verts[k] - flatTriangle[1]).y;
                    float s1 = Vector3.Cross(flatTriangle[2] - flatTriangle[1], verts[k] - flatTriangle[2]).y;
                    float s2 = Vector3.Cross(flatTriangle[0] - flatTriangle[2], verts[k] - flatTriangle[0]).y;
                    if(s0*s1 > 0f && s1*s2 > 0f && s2*s0 > 0f)
                    {
                        overlapTest = true;
                        break;
                    }
                }
                if (overlapTest) data[i, j] = 1;
            }
        }
    }


    private Vector3 IndexToVector(int x, int z)
    {
        return bounds.min + Vector3.right * x * gridSize + Vector3.forward * z * gridSize + Vector3.up * bounds.extents.y;
    }

    private int CalcDir(int dx, int dz)
    {
        if (dx != 0) return dx > 0 ? 0 : 2;
        return dz > 0 ? 3 : 1;
    }

    public void GenerateContour()
    {
        bool foundTopLeft = false;
        for(int i = 0;!foundTopLeft && i < numRows; i++)
        {
            for(int j = 0; j < numCols; j++)
            {
                if(data[i, j] == 1)
                {
                    foundTopLeft = true;
                    m_contour = ContourWalk(j, i, 1);
                    break;
                } 
            }
        }
    }

    private List<Vector3> ContourWalk(int startX, int startZ, int areaId)
    {
        List<Vector3> ct = new List<Vector3>();
        // 0:RIGHT, 1:UP, 2:LEFT, 3:DOWN
        int dir = 2;
        int[] delta = { 1, 0, 0, -1, -1, 0, 0, 1 };
        ct.Add(IndexToVector(startX, startZ));

        int x = startX, z = startZ;
        bool loop = false;
        while(!loop)
        {
            for(int i = 0; i < 4; i++)
            {
                int nx = x + delta[(dir + 3 + i) % 4 * 2 + 0], nz = z + delta[(dir + 3 + i) % 4 * 2 + 1];
                if (nx < 0 || nx >= numCols || nz < 0 || nz >= numRows) continue;
                if (data[nz, nx] == areaId)
                {
                    dir = CalcDir(nx - x, nz - z);
                    x = nx; z = nz;
                    break;
                }
            }
            if (x == startX && z == startZ) loop = true;
            else ct.Add(IndexToVector(x, z));
        }

        return ct;
    }

    private float DistanceToSegment(Vector3 t, Vector3 a, Vector3 b)
    {
        Vector3 ab = b - a;
        Vector3 bt = t - b;
        Vector3 n = Vector3.Dot(ab, ab) * bt - Vector3.Dot(bt, ab) * ab; // triple cross
        n.Normalize();

        return Vector3.Dot(t - a, n);
    }

    public List<Vector3> SimplifyContour(List<Vector3> contour, float threshold = 0.0f)
    {
        if (threshold < 1e-3f) threshold = 2.0f * gridSize;
        contour.Add(contour[0]);
        int i0 = 0;
        int i1 = contour.Count / 4;
        int i2 = contour.Count / 2;
        int i3 = contour.Count * 3 / 4;
        int i4 = contour.Count - 1;
        var p1 = SimplifyContourRecursive(contour, i0, i1, threshold);
        var p2 = SimplifyContourRecursive(contour, i1, i2, threshold);
        var p3 = SimplifyContourRecursive(contour, i2, i3, threshold);
        var p4 = SimplifyContourRecursive(contour, i3, i4, threshold);

        // contour.Clear();
        List<Vector3> simplified = new List<Vector3>();
        simplified.Add(contour[i0]);
        foreach (var p in p1) simplified.Add(contour[p]);
        simplified.Add(contour[i1]);
        foreach (var p in p2) simplified.Add(contour[p]);
        simplified.Add(contour[i2]);
        foreach (var p in p3) simplified.Add(contour[p]);
        simplified.Add(contour[i3]);
        foreach (var p in p4) simplified.Add(contour[p]);
        // simplified.Add(contour[i4]);

        return simplified;
    }

    private List<int> SimplifyContourRecursive(List<Vector3> contour, int s, int t, float threshold)
    {
        List<int> anchor = new List<int>();
        float maxDistance = 0.0f;
        int maxIndex = -1;
        for(int i = s + 1; i < t; i++)
        {
            float distance = DistanceToSegment(contour[i], contour[s], contour[t]);
            if (distance < threshold) continue;
            if(maxIndex == -1 || maxDistance < distance)
            {
                maxDistance = distance;
                maxIndex = i;
            }
        }
        if (maxIndex == -1) return anchor;
        else
        {
            var pre = SimplifyContourRecursive(contour, s, maxIndex, threshold);
            var post = SimplifyContourRecursive(contour, maxIndex, t, threshold);
            foreach (var element in pre) anchor.Add(element);
            anchor.Add(maxIndex);
            foreach (var element in post) anchor.Add(element);
        }
        return anchor;
    }

    public float targetRadius = 0.5f;

    public List<NavmeshArea> GenerateNavmeshAreas()
    {
        List<NavmeshArea> rawAreas = FindAreas();
        List<NavmeshArea> realAreas = new List<NavmeshArea>();
        foreach(var a in rawAreas)
        {
            realAreas.AddRange(PurifyArea(a));
        }

        return realAreas;
    }

    private List<NavmeshArea> PurifyArea(NavmeshArea area)
    {
        List<NavmeshArea> result = new List<NavmeshArea>();
        List<List<Vector3>> contours = OffsetByRadius(area.contour, targetRadius, Vector3.up);
        List<List<Vector3>> holes = new List<List<Vector3>>();
        foreach (var h in area.holes)
            holes.AddRange(OffsetByRadius(h, targetRadius, Vector3.down));

        MergePolygons(holes);

        bool[] shouldRemove = new bool[holes.Count];
        
        for (int i = 0; i < contours.Count; i++)
        {
            for (int j = 0; j < holes.Count; j++)
            {
                List<List<Vector3>> tmp;
                if (PolygonSubtract(contours[i], holes[j], out tmp))
                {
                    contours.RemoveAt(i--);
                    foreach (var tmpList in tmp) contours.Add(tmpList);
                    shouldRemove[j] = true; // the hole cut through the contour, not a hole anymore
                    if (i < 0) break;
                } // else nothing changed
            }
        }
        
        foreach (var c in contours)
        {
            NavmeshArea a = new NavmeshArea();
            a.contour = c;
            for (int i = 0; i < holes.Count; i++)
            {
                if (!shouldRemove[i] && GetInsideness(holes[i][0], c))
                {
                    a.holes.Add(holes[i]);
                }
            }
            result.Add(a);
        }

        return result;
    }


    public List<NavmeshArea> FindAreas()
    {
        List<NavmeshArea> result = new List<NavmeshArea>();
        Dictionary<int, NavmeshArea> dict = new Dictionary<int, NavmeshArea>();
        NavmeshArea extArea = null;
        int[,] visited = new int[numRows, numCols]; // 0 by default
        int[] delta = { 1, 0, 0, -1, -1, 0, 0, 1 };
        for (int i = 0; i < numRows; i++)
        {
            for (int j = 0; j < numCols; j++)
            {
                if (visited[i, j] == 0)
                {
                    int areaId = data[i, j] == 1 ? result.Count + 1 : -1;
                    visited[i, j] = areaId;
                    Queue<int> q = new Queue<int>();
                    q.Enqueue(i * numCols + j);
                    while (q.Count > 0)
                    {
                        int front = q.Dequeue();
                        int x = front % numCols, z = front / numCols;
                        for (int k = 0; k < 4; k++)
                        {
                            int nx = x + delta[k * 2 + 0],
                                nz = z + delta[k * 2 + 1];
                            if (nx < 0 || nx >= numCols || nz < 0 || nz >= numRows) continue;
                            if (visited[nz, nx] != 0) continue;
                            if (data[nz, nx] == data[i, j])
                            {
                                q.Enqueue(nz * numCols + nx);
                                visited[nz, nx] = areaId;
                            }
                        }
                    }
                    if (data[i, j] == 1)
                    {
                        NavmeshArea newArea = new NavmeshArea();
                        newArea.contour = SimplifyContour(ContourWalk(j, i, 1));
                        result.Add(newArea);
                        extArea = newArea;
                    }
                    else if (extArea != null) extArea.holes.Add(SimplifyContour(ContourWalk(j, i, 0)));
                }
                else if(visited[i, j] > 0)
                {
                    extArea = result[visited[i, j] - 1];
                }
            }
        }
        Debug.Log(string.Format("# Navmesh areas: {0}", result.Count));
        return result;
    }



    private void MergePolygons(List<List<Vector3>> polygons)
    {
        int n = polygons.Count;
        int i = n - 1;
        while(i > 0)
        {
            int innerCount = 0;
            for(int j = i - 1; j >= 0; j--)
            {
                innerCount++;
                if (innerCount > 100) throw new Exception("innerCount > 100");
                List<Vector3> newPolygon;
                if(PolygonUnion(polygons[j], polygons[i], out newPolygon))
                {
                    polygons[i] = newPolygon;
                    polygons.RemoveAt(j);
                    break;
                }
            }
            i--;
            if (i < 0) throw new Exception("i < 0");
        }
    }

    public void FindHolesAndWalkables()
    {
        bool[,] visited = new bool[numRows, numCols]; // false by default
        int[] delta = { 1, 0, 0, -1, -1, 0, 0, 1 };
        for(int i = 0; i < numRows; i++)
        {
            for(int j = 0; j < numCols; j++)
            {
                if(!visited[i, j])
                {
                    Queue<int> q = new Queue<int>();
                    q.Enqueue(i * numCols + j);
                    while(q.Count > 0)
                    {
                        int front = q.Dequeue();
                        int x = front % numCols, z = front / numCols;
                        for(int k = 0; k < 4; k++)
                        {
                            int nx = x + delta[k * 2 + 0],
                                nz = z + delta[k * 2 + 1];
                            if (nx < 0 || nx >= numCols || nz < 0 || nz >= numRows) continue;
                            if (visited[nz, nx]) continue;
                            if (data[nz, nx] == data[i, j])
                            {
                                q.Enqueue(nz * numCols + nx);
                                visited[nz, nx] = true;
                            }
                        }
                    }
                    if (data[i, j] == 0) m_holes.Add(SimplifyContour(ContourWalk(j, i, 0)));
                    else if (data[i, j] == 1) walkables.Add(SimplifyContour(ContourWalk(j, i, 1)));
                }
            }
        }
        Debug.Log(string.Format("#Areas: {0}, #Holes: {1}", walkables.Count, m_holes.Count));
        m_holes.Remove(m_holes[0]); // remove the exterior space
    }

    public void ShrinkContour(float radius = 0.5f)
    {
        List<List<Vector3>> shrinkedWalkables = new List<List<Vector3>>();
        foreach (var w in walkables)
        {
            shrinkedWalkables.AddRange(OffsetByRadius(w, radius, Vector3.up));
        }
        walkables = shrinkedWalkables;

        List<List<Vector3>> expandedHoles = new List<List<Vector3>>();
        foreach(var h in m_holes)
        {
            expandedHoles.AddRange(OffsetByRadius(h, radius, Vector3.down));
        }
        m_holes = expandedHoles;
    }

    List<List<Vector3>> clippingResults = new List<List<Vector3>>();
    public void TestClipWalkableArea()
    {
        bool changed = PolygonSubtract(walkables[0], m_holes[0], out clippingResults);
    }



    List<Vector3> debugIntersections = new List<Vector3>();
    
    private bool OnSegment(ref Vector3 p, ref Vector3 q, ref Vector3 r)
    {
        return Mathf.Min(p.x, r.x) - 1e-3f <= q.x && q.x <= Mathf.Max(p.x, r.x) + 1e-3f &&
            Mathf.Min(p.z, r.z) - 1e-3f <= q.z && q.z <= Mathf.Max(p.z, r.z) + 1e-3f;
    }

    private List<List<Vector3>> OffsetByRadius(List<Vector3> contour, float r, Vector3 up)
    {
        VertexNode head = null, tail = null;
        for(int i=0;i<contour.Count;i++)
        {
            int i1 = (i + 1) % contour.Count;
            var n = Vector3.Cross(up, contour[i1] - contour[i]).normalized;
            if (head == null)
            {
                head = new VertexNode(contour[i] + n * r);
                tail = head;
            }
            else
            {
                tail.next = new VertexNode(contour[i] + n * r);
                tail.next.prev = tail;
                tail = tail.next;
            }
            tail.next = new VertexNode(contour[i1] + n * r);
            tail.next.prev = tail;
            tail = tail.next;
        }
        tail.next = new VertexNode(head.v); // form a loop
        tail.next.prev = tail;
        tail = tail.next;

        List<List<Vector3>> polygons = new List<List<Vector3>>();
        Queue<VertexNode> queue = new Queue<VertexNode>();
        queue.Enqueue(head);
        
        while (queue.Count > 0)
        {
            head = queue.Dequeue();
            var p = head;
            while (p.next != null)
            {
                Vector3 v1 = p.v, v2 = p.next.v;  // current segment
                var q = p.next.next;              // skip next segment (they intersect head to tail)
                while (q != null && q.next != null)
                {
                    Vector3 u1 = q.v, u2 = q.next.v;  // segment to detect intersection
                    Vector3 t;
                    if(p == head && q.next.next == null)
                    {
                        // doomed intersection, do nothing
                        break;
                    }
                    if (SegmentIntersection(v1, v2, u1, u2, out t))
                    {
                        if (OnSegment(ref v1, ref t, ref v2) && OnSegment(ref u1, ref t, ref u2))
                        {
                            t.y = 1;
                            debugIntersections.Add(t);

                            var pnext = p.next;
                            var qnext = q.next;

                            p.next = new VertexNode(t);
                            p.next.next = qnext;
                            qnext.prev = p.next;
                            p.next.prev = p;

                            var newHead = new VertexNode(t);
                            newHead.next = pnext;
                            pnext.prev = newHead;

                            q.next = new VertexNode(t);
                            q.next.prev = q;

                            queue.Enqueue(newHead);

                            // queue.Enqueue(head); // keep on running on this contour
                            q = qnext.prev;
                        }
                    }
                    q = q.next;
                }
                p = p.next;
            }
            // if (!splitted) // always simple polygon
            {
                float area = 0;
                p = head;
                List<Vector3> ct = new List<Vector3>();
                while (p.next != null)
                {
                    area += (p.next.v.x - p.v.x) * (p.next.v.z + p.v.z);
                    ct.Add(p.v);
                    p = p.next;
                }
                if (area > 1.0f) // sign same as contour
                {
                    polygons.Add(ct);
                }
            }
        }

        return polygons;
    }

    private bool _IsLeft(Vector3 p, Vector3 v1, Vector3 v2)
    {
        return Vector3.Cross(v2 - v1, p - v2).y > 0;
    }

    private bool _Upward(Vector3 v0, Vector3 v1)
    {
        return v1.z > v0.z;
    }

    private bool _throughZ(Vector3 v0, Vector3 v1, float z)
    {
        return Mathf.Min(v0.z, v1.z) < z && z <= Mathf.Max(v0.z, v1.z);
    }




    private bool GetInsideness(Vector3 point, List<Vector3> v)
    {
        int wn = 0;
        for(int i = 0; i < v.Count; i++)
        {
            int j = (i + 1 == v.Count) ? 0 : i + 1;
            if (_throughZ(v[i], v[j], point.z) && _Upward(v[i], v[j]) && _IsLeft(point, v[i], v[j]))
                wn++;
            else if (_throughZ(v[i], v[j], point.z) && !_Upward(v[i], v[j]) && !_IsLeft(point, v[i], v[j]))
                wn--;
        }
        return wn != 0; // true <--> inside
    }



    private List<VertexWANode> GetWAIntersections(VertexWANode contourHead, VertexWANode clipperHead, bool headInside)
    {
        List<VertexWANode> intersections = new List<VertexWANode>();
        var p = contourHead;
        do
        {
            var pnext = p.Next;
            {
                var q = clipperHead;
                do
                {
                    var qnext = q.Next;
                    Vector3 t;
                    if (SegmentIntersection(p.Vertex, pnext.Vertex, q.Vertex, qnext.Vertex, out t))
                    {
                        if (OnSegment(ref p.Vertex, ref t, ref pnext.Vertex) && OnSegment(ref q.Vertex, ref t, ref qnext.Vertex))
                        {
                            t.y = 1;
                            VertexWANode inter = new VertexWANode(t, VertexWANode.NodeType.Intersection);
                            inter.PolygonPrev = p;
                            inter.PolygonNext = p.Next;     // p.Next !== pnext, if there're multiple insertions
                            inter.ClipperPrev = q;
                            inter.ClipperNext = q.Next;

                            p.Next.Prev = inter;
                            p.Next = inter;
                            q.Next.Prev = inter;
                            q.Next = inter;
                            intersections.Add(inter);       // direction is determined later
                        }
                    }
                    q = qnext;
                } while (q != clipperHead);  // end of clipper do-while
            }
            p = pnext;
        } while (p != contourHead); // end of contour do-while

        var dir = headInside ? VertexWANode.NodeDirection.Out : VertexWANode.NodeDirection.In;
        p = contourHead;
        do
        {
            if (p.Type == VertexWANode.NodeType.Intersection)
            {
                p.Direction = dir;
                dir = (dir == VertexWANode.NodeDirection.In) ? VertexWANode.NodeDirection.Out : VertexWANode.NodeDirection.In;
                p = p.PolygonNext;
            }
            else p = p.Next;
        } while (p != contourHead);
        
        return intersections;
    }


    private VertexWANode MakeLinkedLoop(List<Vector3> polygon, bool reverse = false)
    {
        int n = polygon.Count;
        int d = reverse ? -1 : 1;
        int s = reverse ? n-1 : 0;
        int t = reverse ? -1 : n;

        VertexWANode head = new VertexWANode(polygon[s], VertexWANode.NodeType.Endpoint);
        VertexWANode tail = head;
        for(int i = s + d; i != t; i += d)
        {
            int j = (i + n) % n;
            tail.Next = new VertexWANode(polygon[j], VertexWANode.NodeType.Endpoint);
            tail.Next.Prev = tail;
            tail = tail.Next;
        }
        tail.Next = head;

        return head;
    }


    private bool PolygonSubtract(List<Vector3> contour, List<Vector3> clipper, out List<List<Vector3>> result)
    {
        result = new List<List<Vector3>>();

        bool clipperInsideContour = true;

        /* clipper must be cw to get exterior polygons (contour - hole) */
        VertexWANode clpHead = new VertexWANode(clipper[clipper.Count-1], VertexWANode.NodeType.Endpoint);
        VertexWANode clpTail = clpHead;
        if (!GetInsideness(clipper[clipper.Count - 1], contour)) clipperInsideContour = false;
        for (int i = clipper.Count - 2; i >= 0; i--)
        {
            if (!GetInsideness(clipper[i], contour)) clipperInsideContour = false;
            clpTail.Next = new VertexWANode(clipper[i], VertexWANode.NodeType.Endpoint);
            clpTail.Next.Prev = clpTail;
            clpTail = clpTail.Next;
        }
        clpTail.Next = clpHead;

        if (clipperInsideContour)
        {
            List<Vector3> tmp = new List<Vector3>(); // copy contour
            tmp.AddRange(contour);
            result.Add(tmp);
            return false; // not changed
        }
        
        bool contourInsideClipper = true;
        Dictionary<int, bool> insideClipper = new Dictionary<int, bool>();
        VertexWANode ctrHead = new VertexWANode(contour[0], VertexWANode.NodeType.Endpoint);
        VertexWANode ctrTail = ctrHead;
        insideClipper[ctrHead.Id] = GetInsideness(contour[0], clipper);
        if (!insideClipper[ctrHead.Id]) contourInsideClipper = false;
        for (int i = 1; i < contour.Count; i++)
        {
            var node = new VertexWANode(contour[i], VertexWANode.NodeType.Endpoint);
            insideClipper[node.Id] = GetInsideness(contour[i], clipper);
            if (!insideClipper[node.Id]) contourInsideClipper = false;
            ctrTail.Next = node;
            ctrTail.Next.Prev = ctrTail;
            ctrTail = ctrTail.Next;
        }
        ctrTail.Next = ctrHead; // form a loop

        if (contourInsideClipper) return true; // return empty list

        var intersections = GetWAIntersections(ctrHead, clpHead, insideClipper[ctrHead.Id]);
        if (intersections.Count == 0) return false;
        if (intersections.Count > 8) return false; // unlikely, should be collinear boundaries

        // start from outward intersections to get exterior polygons (contour - hole)
        HashSet<int> visited = new HashSet<int>();
        foreach(var it in intersections)
        {
            if (!visited.Contains(it.Id) && it.Direction == VertexWANode.NodeDirection.Out)
            {
                var q = it;
                List<Vector3> polygon = new List<Vector3>();
                do
                {
                    visited.Add(q.Id);
                    polygon.Add(new Vector3(q.Vertex.x, q.Vertex.y, q.Vertex.z)); // copy
                    if (q.Type == VertexWANode.NodeType.Intersection)
                    {
                        if (q.Direction == VertexWANode.NodeDirection.Out)
                            q = q.PolygonNext;
                        else
                            q = q.ClipperNext;
                    }
                    else { q = q.Next; }
                } while (q != it);
                result.Add(polygon);
            }
        }

        return true;
    }


    private bool PolygonUnion(List<Vector3> ct1, List<Vector3> ct2, out List<Vector3> result)
    {
        result = new List<Vector3>();

        VertexWANode head1 = MakeLinkedLoop(ct1);
        VertexWANode head2 = MakeLinkedLoop(ct2);

        var intersections = GetWAIntersections(head1, head2, GetInsideness(head1.Vertex, ct2));
        if (intersections.Count == 0) return false;
        
        foreach(var it in intersections)
        {
            if(it.Direction == VertexWANode.NodeDirection.Out)
            {
                var p = it;
                do
                {
                    result.Add(p.Vertex);
                    if (p.Type == VertexWANode.NodeType.Intersection)
                    {
                        if (p.Direction == VertexWANode.NodeDirection.Out)
                            p = p.PolygonNext;
                        else
                            p = p.ClipperNext;
                    }
                    else p = p.Next;
                } while (p != it);
                break; // only one union polygon can be made
            }
        }
        return true;
    }

    private bool SegmentIntersection(Vector3 v1, Vector3 v2, Vector3 u1, Vector3 u2, out Vector3 intersection)
    {
        float a1 = v2.z - v1.z, b1 = v1.x - v2.x, c1 = a1 * v1.x + b1 * v1.z;
        float a2 = u2.z - u1.z, b2 = u1.x - u2.x, c2 = a2 * u1.x + b2 * u1.z;
        float det = a1 * b2 - a2 * b1;
        if(det == 0.0f)
        {
            intersection = Vector3.zero;
            return false;
        }
        else
        {
            intersection = new Vector3((b2 * c1 - b1 * c2) / det, 0, (a1 * c2 - a2 * c1) / det);
            return true;
        }
    }

    public void DoTestPartition()
    {
        NavmeshArea a = areas[0];
        
        // partition among points on contour and holes
        for (int i = 0; i < a.contour.Count - 2; i ++)
        {
            for (int j = i + 2; j < a.contour.Count; j++)
            {
                // intersect test
                Partition p = new Partition
                {
                    points = new List<Vector3>()
                };
                p.points.Add(a.contour[i]);
                p.points.Add(a.contour[j]);

                bool intersect = false;

                intersect = IsIntersectWithOtherPartitions(p);
                if (!intersect)
                {
                    intersect = IsIntersectWithContourEdges(p, a.contour);
                }
                if (!intersect)
                {
                    foreach (var h in a.holes)
                    {
                        intersect = IsIntersectWithContourEdges(p, h);
                        
                        if (intersect)
                            break;
                    }
                }
                if (intersect)
                {
                    continue;
                }
                // angle test 重合点会有bug
                {
                    Vector3 partiVec = a.contour[i] - a.contour[j];
                    Vector3 farVec = a.contour[(j + 1) % a.contour.Count] - a.contour[j];
                    Vector3 closeVec = a.contour[(j - 1) % a.contour.Count] - a.contour[j];

                    float angleJ = Vector3.Angle(farVec, closeVec);
                    if (Vector3.Cross(closeVec, farVec).y > 0f)
                        angleJ = 360 - angleJ;
                    float angleI = Vector3.Angle(partiVec, closeVec);
                    if (Vector3.Cross(closeVec, partiVec).y > 0f)
                        angleI = 360 - angleI;

                    if (angleI > angleJ)
                    {
                        continue;
                    }
                }

                partitions.Add(p);
            }
            foreach (var hole in a.holes)
            {
                for (int j = 0; j < hole.Count; j ++)
                {
                    // intersect test
                    Partition p = new Partition
                    {
                        points = new List<Vector3>()
                    };
                    p.points.Add(a.contour[i]);
                    p.points.Add(hole[j]);

                    bool intersect = false;

                    intersect = IsIntersectWithOtherPartitions(p);
                    if (!intersect)
                    {
                        intersect = IsIntersectWithContourEdges(p, a.contour);
                    }
                    if (!intersect)
                    {
                        foreach (var h in a.holes)
                        {
                            intersect = IsIntersectWithContourEdges(p, h);

                            if (intersect)
                                break;
                        }
                    }
                    if (intersect)
                    {
                        continue;
                    }

                    partitions.Add(p);
                }
            }
        }
        // partition among points on holes
        for (int t = 0; t < a.holes.Count; t ++)
        {
            for (int i = 0; i < a.holes[t].Count - 2; i ++)
            {
                for (int j = i + 2; j < a.holes[t].Count; j++)
                {
                    // intersect test
                    Partition p = new Partition
                    {
                        points = new List<Vector3>()
                    };
                    p.points.Add(a.holes[t][i]);
                    p.points.Add(a.holes[t][j]);

                    bool intersect = false;

                    intersect = IsIntersectWithOtherPartitions(p);
                    if (!intersect)
                    {

                        intersect = IsIntersectWithContourEdges(p, a.contour);
                    }
                    if (!intersect)
                    {
                        foreach (var h in a.holes)
                        {
                            intersect = IsIntersectWithContourEdges(p, h);

                            if (intersect)
                                break;
                        }
                    }
                    if (intersect)
                    {
                        continue;
                    }
                    // angle test
                    {
                        Vector3 partiVec = a.holes[t][i] - a.holes[t][j];
                        Vector3 farVec = a.holes[t][(j + 1) % a.holes[t].Count] - a.holes[t][j];
                        Vector3 closeVec = a.holes[t][(j - 1) % a.holes[t].Count] - a.holes[t][j];

                        float angleJ = Vector3.Angle(farVec, closeVec);
                        if (Vector3.Cross(closeVec, farVec).y > 0f)
                            angleJ = 360 - angleJ;
                        float angleI = Vector3.Angle(partiVec, closeVec);
                        if (Vector3.Cross(closeVec, partiVec).y > 0f)
                            angleI = 360 - angleI;

                        if (angleI < angleJ)
                        {
                            continue;
                        }
                    }
                    partitions.Add(p);
                }
            }

            // to points on other holes
            if (t != a.holes.Count - 1)
            {
                for (int s = t + 1; s < a.holes.Count; s++)
                {
                    for (int i = 0; i < a.holes[t].Count; i ++)
                    {
                        for (int j = 0; j < a.holes[s].Count; j ++)
                        {
                            Partition p = new Partition
                            {
                                points = new List<Vector3>()
                            };
                            p.points.Add(a.holes[t][i]);
                            p.points.Add(a.holes[s][j]);

                            bool intersect = false;

                            intersect = IsIntersectWithOtherPartitions(p);
                            if (!intersect)
                            {
                                intersect = IsIntersectWithContourEdges(p, a.contour);
                            }
                            if (!intersect)
                            {
                                foreach (var h in a.holes)
                                {
                                    intersect = IsIntersectWithContourEdges(p, h);

                                    if (intersect)
                                        break;
                                }
                            }
                            if (intersect)
                            {
                                continue;
                            }

                            partitions.Add(p);

                        }
                    }
                }
            }
        }

        Debug.Log(string.Format("# Partitions : {0}", partitions.Count));
    }

    private bool IsIntersectWithOtherPartitions(Partition p)
    {
        for (int k = 0; k < partitions.Count; k++)
        {

            Vector3 left = partitions[k].points[0] - p.points[0];
            Vector3 right = partitions[k].points[1] - p.points[0];
            Vector3 middle = p.points[1] - p.points[0];
            Vector3 left0 = p.points[0] - partitions[k].points[0];
            Vector3 right0 = p.points[1] - partitions[k].points[0];
            Vector3 middle0 = partitions[k].points[1] - partitions[k].points[0];
            if (Vector3.Cross(left, middle).y * Vector3.Cross(right, middle).y < 0f && Vector3.Cross(left0, middle0).y * Vector3.Cross(right0, middle0).y < 0f)
            {
                return true;
            }
        }
        return false;
    }

    private bool IsIntersectWithContourEdges(Partition p, List<Vector3> contour)
    {
        for (int k = 0; k < contour.Count; k++)
        {
            Vector3 left = contour[k] -p.points[0];
            Vector3 right = contour[(k + 1) % contour.Count] -p.points[0];
            Vector3 middle = p.points[1] -p.points[0];
            Vector3 left0 =p.points[0] - contour[k];
            Vector3 right0 = p.points[1] - contour[k];
            Vector3 middle0 = contour[(k + 1) % contour.Count] - contour[k];
            if (Vector3.Cross(left, middle).y * Vector3.Cross(right, middle).y < 0f && Vector3.Cross(left0, middle0).y * Vector3.Cross(right0, middle0).y < 0f)
            {
                return true;
            }
        }
        return false;
    }

    private void OnDrawGizmosSelected()
    {
        // Draw volume
        Gizmos.color = Color.green;
        Gizmos.DrawWireCube(bounds.center, bounds.extents * 2);
    }

    private void OnDrawGizmos()
    {

        // Draw walkable grids
        //Gizmos.color = Color.yellow;
        //for (int i = 0; i < numRows; i++)
        //{
        //    for (int j = 0; j < numCols; j++)
        //    {
        //        if (data[i, j] == 1)
        //        {
        //            Gizmos.DrawWireCube(bounds.min + Vector3.right * j * gridSize + Vector3.forward * i * gridSize + Vector3.up * bounds.extents.y, Vector3.one * gridSize);
        //        }
        //    }
        //}

        //// Draw walkable contour
        //Gizmos.color = Color.red;
        //foreach (var contour in walkables)
        //{
        //    if (contour.Count > 0)
        //    {
        //        for (int i = 0; i < contour.Count - 1; i++)
        //        {
        //            Gizmos.DrawLine(contour[i], contour[i + 1]);
        //        }
        //        Gizmos.DrawLine(contour[contour.Count - 1], contour[0]);
        //    }
        //}


        //Gizmos.color = Color.green;
        //foreach(var t in debugIntersections)
        //{
        //    Gizmos.DrawWireSphere(t, 0.1f);
        //}


        //// Draw holes
        //Gizmos.color = Color.magenta;
        //foreach(var hole in holes)
        //{
        //    for (int i = 0; i < hole.Count - 1; i++)
        //    {
        //        Gizmos.DrawLine(hole[i], hole[i + 1]);
        //    }
        //    Gizmos.DrawLine(hole[hole.Count - 1], hole[0]);
        //}


        //// Draw test polygon subtraction
        //Gizmos.color = new Color(0.5f, 0.8f, 0.8f);
        //foreach(var p in clippingResults)
        //{
        //    for(int i = 0; i < p.Count - 1; i++)
        //    {
        //        Gizmos.DrawLine(p[i], p[i + 1]);
        //    }
        //    Gizmos.DrawLine(p[p.Count - 1], p[0]);
        //}

        if (partitions != null)
        {
            foreach (var p in partitions)
            {
                Gizmos.color = Color.yellow;
                Gizmos.DrawLine(p.points[0], p.points[1]);
            }
        }

        if(areas != null)
        {
            foreach(var a in areas)
            {
                Gizmos.color = new Color(1f, 0f, 0f);
                Gizmos.DrawLine(a.contour[0], a.contour[1]);
                Gizmos.color = new Color(0.5f, 0.8f, 0.8f);

                for (int i = 1; i < a.contour.Count - 1; i++)
                {
                    Gizmos.DrawLine(a.contour[i], a.contour[i + 1]);
                }
                Gizmos.color = new Color(0f, 0f, 1f);
                Gizmos.DrawLine(a.contour[a.contour.Count - 1], a.contour[0]);

                
                foreach (var hole in a.holes)
                {
                    Gizmos.color = new Color(1f, 0f, 0f);
                    Gizmos.DrawLine(hole[0], hole[1]);
                    Gizmos.color = Color.magenta;
                    for (int i = 1; i < hole.Count - 1; i++)
                    {
                        Gizmos.DrawLine(hole[i], hole[i + 1]);
                    }
                    Gizmos.color = new Color(0f, 0f, 1f);

                    Gizmos.DrawLine(hole[hole.Count - 1], hole[0]);
                }
            }
        }
    }
}
