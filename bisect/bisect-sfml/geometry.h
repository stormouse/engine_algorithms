#pragma once

#include "linalg.h"

using float2 = linalg::float2;


static float cross2d(float2 a, float2 b)
{
	return a.x * b.y - a.y * b.x;
}


static bool ray_intersection(float2 p, float2 dir, float2 p0, float2 p1, float2& intersection)
{
	float2 v1 = p - p0;
	float2 v2 = p1 - p0;
	float2 v3(-dir.y, dir.x);

	float dot = linalg::dot(v2, v3);
	if (abs(dot) < 0.000001)
		return false;

	float t1 = cross2d(v2, v1) / dot;
	float t2 = linalg::dot(v1, v3) / dot;

	if (t1 >= 0.0f && (t2 >= 0.0f && t2 <= 1.0f))
	{
		intersection = p + dir * t1;
		return true;
	}
	
	return false;
}


static float2 inner_normal(float2 p0, float2 p1, float2 p2)
{
	float2 d = p1 - p0;
	float2 normal(-d.y, d.x);
	return linalg::dot(normal, p2 - p1) < 0 ? -normal : normal;
}
