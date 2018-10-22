using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class PostRenderer : MonoBehaviour {

    public static PostRenderer Instance { get; private set; }
    private Material lineMat;

    public delegate void PostRenderFunc();
    public PostRenderFunc OnPostRenderCallbacks = null;

    private void Start()
    {
        Instance = this;
        lineMat = new Material(Shader.Find("Unlit/Color"));
        lineMat.color = Color.green;
    }


    private void OnPostRender()
    {
        if(OnPostRenderCallbacks != null)
            OnPostRenderCallbacks.Invoke();
    }
}
