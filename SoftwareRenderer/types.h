#pragma once

#define min3(a, b, c) min(a, min(b, c))
#define max3(a, b, c) max(a, max(b, c))
#define clamp(n, a, b) min(max(n, a), b)

typedef unsigned int uint;
typedef uint8_t u8;
typedef uint32_t u32;
typedef uint64_t u64;

struct vec2
{
	float x;
	float y;
};

struct vec3
{
	float x;
	float y;
	float z;
};

struct vec4
{
	float x;
	float y;
	float z;
	float w;


	vec4() :
		x(0),
		y(0),
		z(0),
		w(0)
	{}
	vec4(float _x, float _y, float _z, float _w) :
		x(_x),
		y(_y),
		z(_z),
		w(_w)
	{}
};

struct mat4x4
{
	float a0; float a1; float a2; float a3;
	float b0; float b1; float b2; float b3;
	float c0; float c1; float c2; float c3;
	float d0; float d1; float d2; float d3;
};

struct quaternion
{
	float w;
	float x;
	float y;
	float z;
};

vec4 operator/(vec4 v, float f)
{
	vec4 rv = vec4(
		v.x / f,
		v.y / f,
		v.z / f, 
		v.w / f
	);

	return rv;
}

vec4 operator-(vec4 v1, vec4 v2)
{
	return vec4(v1.x-v2.x, v1.y-v2.y, v1.z-v2.z, v1.w-v2.w);
}

vec4 normalized(vec4 v)
{
	float magnitude = sqrtf(v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w);
	return v / magnitude; 
}

float magnitude(vec4 v)
{
	return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w);
}

float dot(vec3 a, vec3 b)
{
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

float magnitude(vec3 v)
{
	return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
}

vec3 normalized(vec3 v)
{
	float mag = magnitude(v);
	vec3 rv = { v.x / mag, v.y / mag, v.z / mag };
	return rv;
}

quaternion RotationAroundAxis(float theta, vec4 axis)
{
	float halfTheta = theta / 2.0f;
	float sinHalfTheta = sin(halfTheta);

	axis = normalized(axis);

	quaternion q = {
		cos(halfTheta),
		axis.x * sinHalfTheta,
		axis.y * sinHalfTheta,
		axis.z * sinHalfTheta
	};

	return q;
}

quaternion QuaternionInverse(quaternion q)
{
	float sqrMagnitude = q.w*q.w + q.x*q.x + q.y*q.y + q.z*q.z;
	quaternion inverse = { 
		q.w / sqrMagnitude, 
		-q.x / sqrMagnitude, 
		-q.y / sqrMagnitude, 
		-q.z / sqrMagnitude 
	};

	return inverse;
}

mat4x4 RotationFromQuaternion(const quaternion& q)
{
	mat4x4 rv = {
		1 - 2*q.y*q.y - 2*q.z*q.z, 2*q.x*q.y - 2*q.w*q.z,     2*q.x*q.z + 2*q.w*q.y,     0.0f,
		2*q.x*q.y + 2*q.w*q.z,     1 - 2*q.x*q.x - 2*q.z*q.z, 2*q.y*q.z - 2*q.w*q.x,     0.0f,
		2*q.x*q.z - 2*q.w*q.y,     2*q.y*q.z + 2*q.w*q.x,     1 - 2*q.x*q.x - 2*q.y*q.y, 0.0f,
		0.0f,                      0.0f,                      0.0f,                      1.0f
	};

	return rv;
}

quaternion operator*(const quaternion& q1, const quaternion& q2)
{
	quaternion rv = {
		q1.w*q2.w - q1.x*q2.x - q1.y*q2.y - q1.z*q2.z,
		q1.w*q2.x + q1.x*q2.w + q1.y*q2.z - q1.z*q2.y,
		q1.w*q2.y - q1.x*q2.z + q1.y*q2.w + q1.z*q2.x,
		q1.w*q2.z + q1.z*q2.y - q1.y*q2.x + q1.z*q2.w
	};
	return rv;
}

/*
vec4 operator*(const quaternion& q, const vec4& v)
{
	assert(0 && "Incomplete");
	// TODO: this math is wrong.  Simply upconverting the
	// vector to a quaternion is not the right approach.
	// Most engines appear to just convert to a rotation matrix
	// to apply to a vector.
	//
	// The only example I've found of actual vector / quaternion
	// multiplication is here:
	// https://github.com/okamstudio/godot/blob/master/core/math/quat.h
	//
	// Need to analyse what, exactly, is going on there.  Perhaps 
	// start with the equivlent matrix product and see if it can be
	// reasonably simplified?
	quaternion inverse = QuaternionInverse(q);
	quaternion vecAsQuat = { 1.0f, v.x, v.y, v.z };

	quaternion result = q * vecAsQuat * inverse;

	vec4 rv = { result.x, result.y, result.z, v.w };

	return rv;
}*/

mat4x4 operator*(const mat4x4& mA, const mat4x4& mB)
{
	mat4x4 rv;

	rv.a0 = mA.a0 * mB.a0 + mA.a1 * mB.b0 + mA.a2 * mB.c0 + mA.a3 * mB.d0;
	rv.a1 = mA.a0 * mB.a1 + mA.a1 * mB.b1 + mA.a2 * mB.c1 + mA.a3 * mB.d1;
	rv.a2 = mA.a0 * mB.a2 + mA.a1 * mB.b2 + mA.a2 * mB.c2 + mA.a3 * mB.d2;
	rv.a3 = mA.a0 * mB.a3 + mA.a1 * mB.b3 + mA.a2 * mB.c3 + mA.a3 * mB.d3;
    
	rv.b0 = mA.b0 * mB.a0 + mA.b1 * mB.b0 + mA.b2 * mB.c0 + mA.b3 * mB.d0;
	rv.b1 = mA.b0 * mB.a1 + mA.b1 * mB.b1 + mA.b2 * mB.c1 + mA.b3 * mB.d1;
	rv.b2 = mA.b0 * mB.a2 + mA.b1 * mB.b2 + mA.b2 * mB.c2 + mA.b3 * mB.d2;
	rv.b3 = mA.b0 * mB.a3 + mA.b1 * mB.b3 + mA.b2 * mB.c3 + mA.b3 * mB.d3;
    
	rv.c0 = mA.c0 * mB.a0 + mA.c1 * mB.b0 + mA.c2 * mB.c0 + mA.c3 * mB.d0;
	rv.c1 = mA.c0 * mB.a1 + mA.c1 * mB.b1 + mA.c2 * mB.c1 + mA.c3 * mB.d1;
	rv.c2 = mA.c0 * mB.a2 + mA.c1 * mB.b2 + mA.c2 * mB.c2 + mA.c3 * mB.d2;
	rv.c3 = mA.c0 * mB.a3 + mA.c1 * mB.b3 + mA.c2 * mB.c3 + mA.c3 * mB.d3;

	rv.d0 = mA.d0 * mB.a0 + mA.d1 * mB.b0 + mA.d2 * mB.c0 + mA.d3 * mB.d0;
	rv.d1 = mA.d0 * mB.a1 + mA.d1 * mB.b1 + mA.d2 * mB.c1 + mA.d3 * mB.d1;
	rv.d2 = mA.d0 * mB.a2 + mA.d1 * mB.b2 + mA.d2 * mB.c2 + mA.d3 * mB.d2;
	rv.d3 = mA.d0 * mB.a3 + mA.d1 * mB.b3 + mA.d2 * mB.c3 + mA.d3 * mB.d3;

	return rv;
}

vec4 operator*(const mat4x4& m, const vec4& v)
{
	vec4 rv;

    rv.x = v.x * m.a0 + v.y * m.a1 + v.z * m.a2 + v.w * m.a3;
    rv.y = v.x * m.b0 + v.y * m.b1 + v.z * m.b2 + v.w * m.b3;
    rv.z = v.x * m.c0 + v.y * m.c1 + v.z * m.c2 + v.w * m.c3;
    rv.w = v.x * m.d0 + v.y * m.d1 + v.z * m.d2 + v.w * m.d3;

	return rv;
}

mat4x4 FrustumMatrix(float r, float l, float t, float b, float n, float f)
{
	mat4x4 rv = { 0 };

	float width = r - l;
	float height = t - b;
	float depth = f - n;

	rv.a0 = (2 * n) / width;
	rv.a2 = (r + l) / width;
	rv.b1 = (2 * n) / height;
	rv.b2 = (t + b) / height;
	rv.c2 = -(f + n) / depth;
	rv.c3 = -(2 * n * f) / depth;
	rv.d2 = -1;

	return rv;
}

mat4x4 TranslationMatrix(float dx, float dy, float dz)
{
	mat4x4 rv = { 1.0f, 0.0f, 0.0f, dx,
				  0.0f, 1.0f, 0.0f, dy,
				  0.0f, 0.0f, 1.0f, dz,
				  0.0f, 0.0f, 0.0f, 1.0f };
	return rv;
}


mat4x4 MakeTransformMatrix(quaternion rotation, vec3 scale, vec3 position)
{
	mat4x4 quatMatrix = RotationFromQuaternion(rotation);
	mat4x4 transform = {
		scale.x * quatMatrix.a0, scale.y * quatMatrix.a1, scale.z * quatMatrix.a2, position.x,
		scale.x * quatMatrix.b0, scale.y * quatMatrix.b1, scale.z * quatMatrix.b2, position.y,
		scale.x * quatMatrix.c0, scale.y * quatMatrix.c1, scale.z * quatMatrix.c2, position.z,
		0.0f,                    0.0f,                    0.0f,                    1.0f
	};
	return transform;
}

/*
            [cos(theta) -sin(theta) 0]
Rz(theta) = |sin(theta)  cos(theta) 0|
            [0           0          1]

            [1 0          0          ]
Rx(theta) = |0 cos(theta) -sin(theta)|
            [0 sin(theta)  cos(theta)]


            [ cos(theta)  0 sin(theta)]
Ry(theta) = | 0           1 0         |
            [-sin(theta)  0 cos(theta)]
*/

enum RotationAxis
{
	RA_X,
	RA_Y,
	RA_Z,
};

mat4x4 RotationMatrix(float theta, RotationAxis axis)
{
	mat4x4 rv = { 0.0f };

	float c = cos(theta);
	float s = sin(theta);
	switch (axis)
	{
	case RA_X:
		rv.a0 = 1.0f;
		rv.b1 = c;
		rv.b2 = -s;
		rv.c1 = s;
		rv.c2 = c;
		rv.d3 = 1.0f;
		break;
	case RA_Y:
		rv.a0 = c;
		rv.a2 = s;
		rv.b1 = 1.0f;
		rv.c0 = -s;
		rv.c2 = c;
		rv.d3 = 1.0f;
		break;
	case RA_Z:
		rv.a0 = c;
		rv.a1 = -s;
		rv.b0 = s;
		rv.b1 = c;
		rv.c2 = 1.0f;
		rv.d3 = 1.0f;
		break;
	default:
		assert(0);
	}
	return rv;
}

#define RGBA32(r, g, b, a) (((a & 0xff) << 24) | ((r & 0xff) << 16) | ((g & 0xff) << 8) | (b & 0xff))

struct Color
{
	float r;
	float g;
	float b;
	float a;

	Color() :
		r(0.0f),
		g(0.0f),
		b(0.0f),
		a(1.0f)
	{
	}

	Color(float _r, float _g, float _b, float _a) :
		r(clamp(_r, 0.0f, 1.0f)),
		g(clamp(_g, 0.0f, 1.0f)),
		b(clamp(_b, 0.0f, 1.0f)),
		a(clamp(_a, 0.0f, 1.0f))
	{
	}

	Color(u32 rgba)
	{
		a = ((rgba & 0xff000000) >> 24) / 255.0f;
		r = ((rgba & 0x00ff0000) >> 16) / 255.0f;
		g = ((rgba & 0x0000ff00) >> 8 ) / 255.0f;
		b = ((rgba & 0x000000ff)      ) / 255.0f;
	}

	u32 rgba() const
	{
		u8 alpha = (u8)(a * 255.0f);
		u8 red   = (u8)(r * 255.0f);
		u8 green = (u8)(g * 255.0f);
		u8 blue  = (u8)(b * 255.0f);

		return RGBA32(red, green, blue, alpha);
	}

	Color& operator*=(float f)
	{
		f = clamp(f, 0.0f, 1.0f);
		r*=f;
		g*=f;
		b*=f;

		return *this;
	}

	Color& operator*=(const Color& c2)
	{
		r*=c2.r;
		g*=c2.g;
		b*=c2.b;
		a*=c2.a;

		return *this;
	}
};

Color operator*(const Color& c1, float f)
{
	f = clamp(f, 0.0f, 1.0f);
	return Color(c1.r*f, c1.g*f, c1.b*f, c1.a);
}

Color operator*(const Color& c1, const Color& c2)
{
	return Color(c1.r*c2.r, c1.g*c2.g, c1.b*c2.b, c1.a*c2.a);
}

Color operator+(const Color& c1, const Color& c2)
{
	float r = clamp((c1.r + c2.r), 0.0f, 1.0f);
	float g = clamp((c1.g + c2.g), 0.0f, 1.0f);
	float b = clamp((c1.b + c2.b), 0.0f, 1.0f);
	float a = clamp((c1.a + c2.a), 0.0f, 1.0f);

	return Color(r, g, b, a);
}

struct Texture
{
	uint width;
	uint height;
	char * data;
};

struct Material
{
	Color ambientColor;
	Color diffuseColor;
	Color specularColor;
	float specularIntensity;
	u8 illumType;

	Texture * ambientTexture;
	Texture * diffuseTexture;
	char * name;
};

struct Mesh
{
	u32 vertexCount;
	vec4 * vertices;
	u32 uvwCount;
	vec3 * uvws;
	u32 normalCount;
	vec3 * normals;

	u32 indexCount;
	uint * indices;

	Material * material;
	char * name;
};
