#pragma once

#include <random>
#include <math.h>

extern std::default_random_engine _engine = std::default_random_engine(1);
extern std::uniform_real_distribution<float> _r_distribution = std::uniform_real_distribution<float>();

static float random ()
{
	return _r_distribution ( _engine );
}

static int random_int ( int min, int max )
{
	float r = _r_distribution ( _engine );
	return (int)floorf ( ( max - min ) * r + min );
}

static int random_int ( int max )
{
	return random_int ( 0, max );
}
