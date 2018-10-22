using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class NavigationVolume : MonoBehaviour {

    public float gridSize;
    public Bounds bounds;

    private int numRows;
    private int numCols;

    private int[,] data;

    private List<Vector3> m_contour = new List<Vector3>();
    private List<List<Vector3>> holes = new List<List<Vector3>>();
    private List<List<Vector3>> walkables = new List<List<Vector3>>();

    private void Start()
    {
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
                    if (data[i, j] == 0) holes.Add(SimplifyContour(ContourWalk(j, i, 0)));
                    else if (data[i, j] == 1) walkables.Add(SimplifyContour(ContourWalk(j, i, 1)));
                }
            }
        }
        Debug.Log(string.Format("#Areas: {0}, #Holes: {1}", walkables.Count, holes.Count));
        holes.Remove(holes[0]); // remove the exterior space
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
        foreach(var h in holes)
        {
            expandedHoles.AddRange(OffsetByRadius(h, radius, Vector3.down));
        }
        holes = expandedHoles;
    }

    List<List<Vector3>> clippingResults = new List<List<Vector3>>();
    public void TestClipWalkableArea()
    {
        clippingResults = PolygonSubtract(walkables[0], holes[0]);
    }

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

            string debugOutput = "";
            float a = 0.0f;
            while (p.next != null)
            {
                a += (p.next.v.x - p.v.x) * (p.next.v.z + p.v.z);
                debugOutput += p.id.ToString() + " ";
                p = p.next;
            }

            p = head;
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


    private List<List<Vector3>> PolygonSubtract(List<Vector3> contour, List<Vector3> clipper)
    {
        List<List<Vector3>> result = new List<List<Vector3>>();

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
            return result;
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

        if (contourInsideClipper) return result; // return empty list

        List<VertexWANode> intersections = new List<VertexWANode>();
        
        var p = ctrHead;
        do
        {
            var pnext = p.Next;
            if(insideClipper[p.Id] != insideClipper[pnext.Id])
            {
                var q = clpHead;
                do
                {
                    var qnext = q.Next;
                    Vector3 t;
                    if(SegmentIntersection(p.Vertex, pnext.Vertex, q.Vertex, qnext.Vertex, out t))
                    {
                        if(OnSegment(ref p.Vertex, ref t, ref pnext.Vertex) && OnSegment(ref q.Vertex, ref t, ref qnext.Vertex))
                        {
                            t.y = 1;
                            VertexWANode inter = new VertexWANode(t, VertexWANode.NodeType.Intersection);
                            inter.PolygonPrev = p;
                            inter.PolygonNext = pnext;
                            inter.ClipperPrev = q;
                            inter.ClipperNext = qnext;
                            inter.Direction = (insideClipper[p.Id] == false) ? VertexWANode.NodeDirection.In : VertexWANode.NodeDirection.Out;

                            p.Next = inter;
                            pnext.Prev = inter;
                            q.Next = inter;
                            qnext.Prev = inter;
                            intersections.Add(inter);
                        }
                    }
                    q = qnext;
                } while (q != clpHead);  // end of clipper do-while
            }
            p = pnext;
        }
        while (p != ctrHead); // end of contour do-while

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
                    polygon.Add(q.Vertex);
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

        return result;
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

    private void OnDrawGizmos()
    {
        //// Draw volume
        //Gizmos.color = Color.green;
        //Gizmos.DrawWireCube(bounds.center, bounds.extents * 2);

        //// Draw walkable grids
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

        // Draw walkable contour
        Gizmos.color = Color.red;
        foreach (var contour in walkables)
        {
            if (contour.Count > 0)
            {
                for (int i = 0; i < contour.Count - 1; i++)
                {
                    Gizmos.DrawLine(contour[i], contour[i + 1]);
                }
                Gizmos.DrawLine(contour[contour.Count - 1], contour[0]);
            }
        }


        Gizmos.color = Color.green;
        foreach(var t in debugIntersections)
        {
            Gizmos.DrawWireSphere(t, 0.1f);
        }


        // Draw holes
        Gizmos.color = Color.magenta;
        foreach(var hole in holes)
        {
            for (int i = 0; i < hole.Count - 1; i++)
            {
                Gizmos.DrawLine(hole[i], hole[i + 1]);
            }
            Gizmos.DrawLine(hole[hole.Count - 1], hole[0]);
        }


        // Draw test polygon subtraction
        Gizmos.color = new Color(0.5f, 0.8f, 0.8f);
        foreach(var p in clippingResults)
        {
            for(int i = 0; i < p.Count - 1; i++)
            {
                Gizmos.DrawLine(p[i], p[i + 1]);
            }
            Gizmos.DrawLine(p[p.Count - 1], p[0]);
        }
    }
}
