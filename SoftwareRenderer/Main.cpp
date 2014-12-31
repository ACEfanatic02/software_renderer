#include <Windows.h>
#include <cassert>

#define _USE_MATH_DEFINES 1
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cfloat>

// For std::sort
#include <algorithm> 

// Disable warnings about C runtime functions
// In a production setting this is a bad idea, but this is just proof-of-concept.
#pragma warning(disable: 4996)

#define min3(a, b, c) min(a, min(b, c))
#define max3(a, b, c) max(a, max(b, c))
#define clamp(n, a, b) min(max(n, a), b)

typedef unsigned int uint;
typedef uint8_t u8;
typedef uint32_t u32;
typedef uint64_t u64;

static const int gScreenWidth = 1280;
static const int gScreenHeight = 720;
static bool gRunning;

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

vec4 normalized(vec4 v)
{
	float magnitude = sqrtf(v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w);
	return v / magnitude; 
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

void DEBUGPrintMat4x4(const mat4x4& m)
{
	wchar_t buffer[256] = {0};

	swprintf(buffer, _countof(buffer), L"[ %.4f, %.4f, %.4f, %.4f ]\n", m.a0, m.a1, m.a2, m.a3);
	OutputDebugStringW(buffer);

	swprintf(buffer, _countof(buffer), L"[ %.4f, %.4f, %.4f, %.4f ]\n", m.b0, m.b1, m.b2, m.b3);
	OutputDebugStringW(buffer);

	swprintf(buffer, _countof(buffer), L"[ %.4f, %.4f, %.4f, %.4f ]\n", m.c0, m.c1, m.c2, m.c3);
	OutputDebugStringW(buffer);

	swprintf(buffer, _countof(buffer), L"[ %.4f, %.4f, %.4f, %.4f ]\n", m.d0, m.d1, m.d2, m.d3);
	OutputDebugStringW(buffer);
}

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

void DEBUGPrintVec4(const vec4& v)
{
	wchar_t buffer[256] = {0};
	swprintf(buffer, _countof(buffer), L"< %.4f, %.4f, %.4f, %.4f >\n", v.x, v.y, v.z, v.w);
	OutputDebugStringW(buffer);
}

void TestMatrixMultiply()
{
	float r = (gScreenWidth / 2.0f);
	float l = -(gScreenWidth / 2.0f);
	float t = (gScreenHeight / 2.0f);
	float b = -(gScreenHeight / 2.0f);
	float n = 1.0f;
	float f = 100.0f;
	mat4x4 test = FrustumMatrix(r, l, t, b, n, f);
	mat4x4 ident = { 0.0f };

	ident.a0 = 1.0f;
	ident.b1 = 1.0f;
	ident.c2 = 1.0f;
	ident.d3 = 1.0f;

	vec4 testVec( -100.0f, -100.0f, 100.0f, 1.0f );

	OutputDebugStringW(L"Original:\n");

	DEBUGPrintMat4x4(test);

	OutputDebugStringW(L"\nAfter multiply by identity:\n");

	DEBUGPrintMat4x4(test * ident);

	OutputDebugStringW(L"\nTest vector:\n");

	DEBUGPrintVec4(testVec);

	OutputDebugStringW(L"\nTransformed vector:\n");

	vec4 transformed = test * testVec;
	DEBUGPrintVec4(transformed);

	OutputDebugStringW(L"\n90degree rotation matrix around Y:\n");

	mat4x4 rotMat = RotationMatrix((float)M_PI_2, RA_Y);
	DEBUGPrintMat4x4(rotMat);
	DEBUGPrintVec4(rotMat * testVec);

	OutputDebugStringW(L"\n90 degree quaternion rotation around Y:\n");

	vec4 yAxis = vec4( 0.0f, 1.0f, 0.0f, 0.0f );
	quaternion quat = RotationAroundAxis((float)M_PI_2, yAxis);
//	DEBUGPrintVec4(quat * testVec);
}

LRESULT CALLBACK
Win32MainWindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
	switch (uMsg)
	{
	case WM_CLOSE:
	case WM_DESTROY:
		{
			gRunning = false;
			PostQuitMessage(0);
		} break;
	default:
		return DefWindowProc(hwnd, uMsg, wParam, lParam);
	}

	return 0;
}

struct win32_backbuffer
{
	HDC bmpHdc;
	BITMAPINFO bmpInfo;
	void * bmpMemory;
	HBITMAP bmpHandle;
	float * depthBuffer;
};

static void 
Win32AllocateDIBSection(win32_backbuffer * backbuffer, int width, int height)
{
	backbuffer->bmpHdc = CreateCompatibleDC(0);

	backbuffer->bmpInfo.bmiHeader.biSize = sizeof(backbuffer->bmpInfo.bmiHeader);
	backbuffer->bmpInfo.bmiHeader.biWidth = width;
	backbuffer->bmpInfo.bmiHeader.biHeight = height;
	backbuffer->bmpInfo.bmiHeader.biPlanes = 1;
	backbuffer->bmpInfo.bmiHeader.biBitCount = 32;
	backbuffer->bmpInfo.bmiHeader.biCompression = BI_RGB;

	backbuffer->bmpHandle = CreateDIBSection(backbuffer->bmpHdc, &backbuffer->bmpInfo, DIB_RGB_COLORS, &backbuffer->bmpMemory, NULL, 0);

	assert(backbuffer->bmpHandle);

	backbuffer->depthBuffer = (float *)calloc(sizeof(float), width * height);
}

static void
Win32ClearBackbuffer(win32_backbuffer * backbuffer, u32 color)
{
	int width = backbuffer->bmpInfo.bmiHeader.biWidth;
	int height = backbuffer->bmpInfo.bmiHeader.biHeight;

	u32 * pixels = (u32 *)backbuffer->bmpMemory;
	for (int i = 0; i < width * height; ++i)
	{
		pixels[i] = color;
		backbuffer->depthBuffer[i] = -FLT_MAX;
	}
}

static void 
Win32RedrawWindow(HWND window, int x, int y, int width, int height, win32_backbuffer * backbuffer)
{
	HDC windowDC = GetDC(window);
	StretchDIBits(windowDC, x, y, width, height, 0, 0, 
		backbuffer->bmpInfo.bmiHeader.biWidth, 
		backbuffer->bmpInfo.bmiHeader.biHeight,
		backbuffer->bmpMemory,
		&backbuffer->bmpInfo,
		DIB_RGB_COLORS,
		SRCCOPY);

	ReleaseDC(window, windowDC);
}

#define RGBA32(r, g, b, a) (((a & 0xff) << 24) | ((r & 0xff) << 16) | ((g & 0xff) << 8) | (b & 0xff))

struct Color
{
	float r;
	float g;
	float b;
	float a;

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

struct Texture;

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

static bool ReadLine(const char ** lineStart, const char ** lineEnd, const char ** fileCur) {
    const char *cur = *fileCur;
    if (!*cur) {
        return false; // EOF
    }
 
    while (*cur == ' ' || *cur == '\t') {
        ++cur;
    }
 
    *lineStart = cur;
 
    const char * end = cur;
    while (true) {
        char ch = *cur++;
        if (!ch) {
            --cur;  // Preserve terminating NUL
            break;
        }
        else if (ch == '\n') {
            break;
        }
        else if (!isspace(ch)) {
            end = cur;
        }
    }
 
    *lineEnd = end;
    *fileCur = cur;
 
    return true;
}

int ReadUint8(u8 * u, const char * c)
{
	const char * start = c;
	const char * end;

	while (*c && !isdigit(*c)) ++c;

	*u = (u8)atoi(c);

	while (*c && isdigit(*c)) ++c;
	end = c;

	return end - start;
}

int ReadFloat32(float * f, const char * c)
{
	const char * start = c;
	const char * end;
	while (*c && !isdigit(*c)) ++c;
	
	*f = (float)atof(c);

	while (*c && (isdigit(*c) || *c == '.')) ++c;
	end = c;

	return end - start;
}

int ReadColor(Color * color, const char * c)
{
	color->a = 1.0f; // Colors are given as RGB, no alpha channel.

	const char * start = c;
	const char * end;

	while (*c && !(isdigit(*c))) ++c;

	c += ReadFloat32(&color->r, c);
	c += ReadFloat32(&color->g, c);
	c += ReadFloat32(&color->b, c);

	end = c;

	return end - start;
}

// NOTE: Only loads the first material in a .mtl file.
// TODO: Figure out a proper API for returning multiple materials.  
// We would rather avoid returning allocated memory if possible.
// Store all results in a resource cache and return a pointer into that?
// Then need to allow queries by name for a particular material.
void LoadMaterial(char * filename, Material * material)
{
	FILE * file = fopen(filename, "r");
	fseek(file, 0, SEEK_END);
	u32 fileLength = ftell(file);  // This is a minimum, may overestimate length of text files.
	fseek(file, 0, SEEK_SET);

	char * bytes = (char *)calloc(1, fileLength);
	fileLength = fread(bytes, 1, fileLength, file);

	// Clear material struct
	memset(material, sizeof(Material), 0);

	const char * cur;
	const char * lineStart;
	const char * lineEnd;

	while (ReadLine(&lineStart, &lineEnd, &cur))
	{
		if (lineStart[0] == '#') continue;  // Comment

		if (!strncmp(lineStart, "newmtl", 6))
		{
			if (material->name) 
			{
				// We've already loaded a material.
				// TODO:  Loading for multiple materials.
				break;
			}
			const char * nameStart = lineStart + 6;
			const char * nameEnd = lineEnd;

			int nameLength = (nameEnd - nameStart) + 1;
			material->name = (char *)calloc(1, nameLength);
			strncpy(material->name, nameStart, nameLength);
		}
		else if (!strncmp(lineStart, "Ns", 2))
		{
			ReadFloat32(&material->specularIntensity, lineStart);
		}
		else if (!strncmp(lineStart, "illum", 5))
		{
			ReadUint8(&material->illumType, lineStart);
		}
		else if (*lineStart == 'K')
		{
			lineStart++;
			if      (*lineStart == 'a') ReadColor(&material->ambientColor, lineStart);
			else if (*lineStart == 'd') ReadColor(&material->diffuseColor, lineStart);
			else if (*lineStart == 's') ReadColor(&material->specularColor, lineStart);
		}
		else if (!strncmp(lineStart, "map_", 4))
		{
			// TODO:  Texture map.
		}
	}

	free(bytes);
}

void LoadMesh(char * filename, Mesh * mesh)
{
}

static void 
RenderTestGradient(win32_backbuffer * backbuffer)
{
	int width = backbuffer->bmpInfo.bmiHeader.biWidth;
	int height = backbuffer->bmpInfo.bmiHeader.biHeight;

	u32 * pixels = (u32*)backbuffer->bmpMemory;

	for (int x = 0; x < width; ++x)
	{
		for (int y = 0; y < height; ++y) 
		{
			pixels[y * width + x] = RGBA32(x, y, x-y, 0xff);
		}
	}
}

// Rasterizer.
//
// Based on the triangle rasterizer described by Fabian Giesen here:
// https://fgiesen.wordpress.com/2013/02/08/triangle-rasterization-in-practice/
//

static void
SetPixel(win32_backbuffer * backbuffer, int x, int y, u32 color, float depth)
{
	int width = backbuffer->bmpInfo.bmiHeader.biWidth;
	int height = backbuffer->bmpInfo.bmiHeader.biHeight;

	assert(x >= 0 && x < width);
	assert(y >= 0 && y < height);

	int idx = y * width + x;

	if (depth > backbuffer->depthBuffer[idx])
	{
		u32 * pixels = (u32 *)backbuffer->bmpMemory;
		pixels[idx] = color;
		backbuffer->depthBuffer[idx] = depth;
	}
}

static void
SetPixel(win32_backbuffer * backbuffer, vec4 p, Color c0, Color c1, Color c2, float d0, float d1, float d2, float l0, float l1, float l2)
{
	int width = backbuffer->bmpInfo.bmiHeader.biWidth;
	int height = backbuffer->bmpInfo.bmiHeader.biHeight;

	int x = (int)(p.x + 0.5f);
	int y = (int)(p.y + 0.5f);

	int idx = y * width + x;

	float depth = d0 + l1*(d1-d0) + l2*(d2-d0);
	if (depth > backbuffer->depthBuffer[idx])
	{
		u32 * pixels = (u32 *)backbuffer->bmpMemory;

		float r = c0.r + l1*(c1.r - c0.r) + l2*(c2.r - c0.r);
		float g = c0.g + l1*(c1.g - c0.g) + l2*(c2.g - c0.g);
		float b = c0.b + l1*(c1.b - c0.b) + l2*(c2.b - c0.b);
		float a = c0.a + l1*(c1.a - c0.a) + l2*(c2.a - c0.a);
		Color color(r, g, b, a);

		pixels[idx] = color.rgba();
		backbuffer->depthBuffer[idx] = depth;
	}
}

// Compute (twice) the area of the triangle abc.
static float
Orient2D(const vec4& a, const vec4& b, const vec4& c)
{
	return (b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x);
}

static bool IsTopLeft(const vec4& a, const vec4& b)
{
	return (abs(a.y - b.y) <= 0.5f) // Top edge
		|| (b.x > a.x);             // Left edge
}

static void
Rasterize(win32_backbuffer * backbuffer, vec4 v0, vec4 v1, vec4 v2, Color c0, Color c1, Color c2)
{
	float screenHalfWidth = gScreenWidth / 2.0f;
	float screenHalfHeight = gScreenHeight / 2.0f;

	// Map x and y to screen coordinates.
	v0.x = ((v0.x / v0.w) + 1.0f) * screenHalfWidth;
	v0.y = ((v0.y / v0.w) + 1.0f) * screenHalfHeight;
	v1.x = ((v1.x / v1.w) + 1.0f) * screenHalfWidth;
	v1.y = ((v1.y / v1.w) + 1.0f) * screenHalfHeight;
	v2.x = ((v2.x / v2.w) + 1.0f) * screenHalfWidth;
	v2.y = ((v2.y / v2.w) + 1.0f) * screenHalfHeight;

	// Near-clip test.
	if (v0.w > 0.0f || v1.w > 0.0f || v2.w > 0.0f) return;

	// Calculate triangle bounding box
	float minX = min3(v0.x, v1.x, v2.x);
	float minY = min3(v0.y, v1.y, v2.y);
	float maxX = max3(v0.x, v1.x, v2.x);
	float maxY = max3(v0.y, v1.y, v2.y);

	// Clip bounding box to screen.
	minX = max(minX, 0);
	minY = max(minY, 0);
	maxX = min(maxX, gScreenWidth - 1.0f);
	maxY = min(maxY, gScreenHeight - 1.0f);

	float stepSize = 0.5f; // Pixels to step in each direction.

	// Per-step deltas for barycentric weights.
	float a01 = (v0.y - v1.y) * stepSize;
	float b01 = (v1.x - v0.x) * stepSize;
	float a12 = (v1.y - v2.y) * stepSize;
	float b12 = (v2.x - v1.x) * stepSize;
	float a20 = (v2.y - v0.y) * stepSize;
	float b20 = (v0.x - v2.x) * stepSize;

	vec4 p( minX, minY, 0.0f, 1.0f );

	// Fill-rule bias
	float bias0 = IsTopLeft(v1, v2) ? 0.0f : -1.0f;
	float bias1 = IsTopLeft(v2, v0) ? 0.0f : -1.0f;
	float bias2 = IsTopLeft(v0, v1) ? 0.0f : -1.0f;
	
	// Calculate barycentric coordinates for first pixel.
	float w0_row = Orient2D(v1, v2, p) + bias0;
	float w1_row = Orient2D(v2, v0, p) + bias1;
	float w2_row = Orient2D(v0, v1, p) + bias2;

	// Backface culling.
	if (w0_row + w1_row + w2_row < 0.0f) return;

	// Normalized barycentric coordinates (for interpolation)
	float twice_area = w0_row + w1_row + w2_row;
	float lambda0_row = w0_row / twice_area;
	float lambda1_row = w1_row / twice_area;
	float lambda2_row = w2_row / twice_area;

	float la12 = a12 / twice_area;
	float la20 = a20 / twice_area;
	float la01 = a01 / twice_area;

	float lb12 = b12 / twice_area;
	float lb20 = b20 / twice_area;
	float lb01 = b01 / twice_area;

	for (p.y = minY; p.y <= maxY; p.y += stepSize)
	{
		float l0 = lambda0_row;
		float l1 = lambda1_row;
		float l2 = lambda2_row;

		for (p.x = minX; p.x <= maxX; p.x += stepSize)
		{
			// NOTE (2014-12-30):  This is a (theoretically) slightly cheaper way to handle 
			// this branch, as it results in a single comparison rather than 3 separate 
			// comparisons. However, profiling does not show any difference at the moment. 
			// I am leaving this here, commented out, in case it proves useful at a later time.
//			u32 mask = *(u32 *)&l0 | *(u32 *)&l1 | *(u32 *)&l2;
//			if (~mask & 0x80000000)
			if (l0 >= 0.0f && l1 >= 0.0f && l2 >= 0.0f) 
			{
#ifdef RASTERIZER_SLOW_PATH
				float depth = v0.w + l1*(v1.w - v0.w) + l2*(v2.w - v0.w);
				float r = c0.r + l1*(c1.r - c0.r) + l2*(c2.r - c0.r);
				float g = c0.g + l1*(c1.g - c0.g) + l2*(c2.g - c0.g);
				float b = c0.b + l1*(c1.b - c0.b) + l2*(c2.b - c0.b);
				float a = c0.a + l1*(c1.a - c0.a) + l2*(c2.a - c0.a);
				Color color(r, g, b, a);

				SetPixel(backbuffer, (int)(p.x + 0.5f), (int)(p.y + 0.5f), color.rgba(), depth);
#else
				SetPixel(backbuffer, p, c0, c1, c2, v0.w, v1.w, v2.w, l0, l1, l2);
#endif
			}

			l0 += la12;
			l1 += la20;
			l2 += la01;
		}

		lambda0_row += lb12;
		lambda1_row += lb20;
		lambda2_row += lb01;
	}
}

static Color faceColors[6] = {
	Color(1.0f, 1.0f, 1.0f, 1.0f),
	Color(1.0f, 0.0f, 0.0f, 1.0f),
	Color(0.0f, 1.0f, 0.0f, 1.0f),
	Color(0.0f, 0.0f, 1.0f, 1.0f),
	Color(1.0f, 1.0f, 0.0f, 1.0f),
	Color(1.0f, 0.0f, 1.0f, 1.0f),
};

struct Triangle
{
	vec4 v0;
	vec4 v1;
	vec4 v2;
};

static void RenderMesh(win32_backbuffer * backbuffer, const vec4 * vertices, int vertexCount, const mat4x4& transform)
{
	static const int MAX_VERTICES = 1024;
	assert(vertexCount < MAX_VERTICES);

	float r = (gScreenWidth / 2.0f);
	float l = -(gScreenWidth / 2.0f);
	float t = -(gScreenHeight / 2.0f);
	float b = (gScreenHeight / 2.0f);
	float n = -1000.0f;
	float f = 1000.0f;
	static const mat4x4 frustumMatrix = FrustumMatrix(r, l, t, b, n, f);

	Triangle tris[(MAX_VERTICES / 3) + 1];
	for (int i = 0; i < vertexCount; i += 3)
	{
		int triIndex = i / 3;
		tris[triIndex].v0 = frustumMatrix * transform * vertices[i];
		tris[triIndex].v1 = frustumMatrix * transform * vertices[i + 1];
		tris[triIndex].v2 = frustumMatrix * transform * vertices[i + 2];
	}

	std::sort(&tris[0], &tris[vertexCount / 3], [](const Triangle& a, const Triangle& b) {
		return min3(a.v0.w, a.v1.w, a.v2.w) < min3(b.v0.w, b.v1.w, b.v2.w);
	});

	for (int i = 0; i < vertexCount / 3; ++i)
	{ 
		Rasterize(backbuffer, tris[i].v0, tris[i].v1, tris[i].v2, 
			faceColors[i % 6], faceColors[i % 6], faceColors[i % 6]);
	}
}

static const int meshVertexCount = 36;
static const vec4 meshVerts[meshVertexCount] = {
	vec4( -1.0f, -1.0f, -1.0f, 1.0f ),
	vec4( -1.0f, -1.0f,  1.0f, 1.0f ),
	vec4( -1.0f,  1.0f,  1.0f, 1.0f ),
//
	vec4(  1.0f,  1.0f, -1.0f, 1.0f ),
	vec4( -1.0f, -1.0f, -1.0f, 1.0f ),
	vec4( -1.0f,  1.0f, -1.0f, 1.0f ),
//
	vec4(  1.0f, -1.0f,  1.0f, 1.0f ),
	vec4( -1.0f, -1.0f, -1.0f, 1.0f ),
	vec4(  1.0f, -1.0f, -1.0f, 1.0f ),
//
	vec4(  1.0f,  1.0f, -1.0f, 1.0f ),
	vec4(  1.0f, -1.0f, -1.0f, 1.0f ),
	vec4( -1.0f, -1.0f, -1.0f, 1.0f ),
//
	vec4( -1.0f, -1.0f, -1.0f, 1.0f ),
	vec4( -1.0f,  1.0f,  1.0f, 1.0f ),
	vec4( -1.0f,  1.0f, -1.0f, 1.0f ),
//
	vec4(  1.0f, -1.0f,  1.0f, 1.0f ),
	vec4( -1.0f, -1.0f,  1.0f, 1.0f ),
	vec4( -1.0f, -1.0f, -1.0f, 1.0f ),
//
	vec4( -1.0f,  1.0f,  1.0f, 1.0f ),
	vec4( -1.0f, -1.0f,  1.0f, 1.0f ),
	vec4(  1.0f, -1.0f,  1.0f, 1.0f ),
//
	vec4(  1.0f,  1.0f,  1.0f, 1.0f ),
	vec4(  1.0f, -1.0f, -1.0f, 1.0f ),
	vec4(  1.0f,  1.0f, -1.0f, 1.0f ),
//		
	vec4(  1.0f, -1.0f, -1.0f, 1.0f ),
	vec4(  1.0f,  1.0f,  1.0f, 1.0f ),
	vec4(  1.0f, -1.0f,  1.0f, 1.0f ),
//		
	vec4(  1.0f,  1.0f,  1.0f, 1.0f ),
	vec4(  1.0f,  1.0f, -1.0f, 1.0f ),
	vec4( -1.0f,  1.0f, -1.0f, 1.0f ),
//
	vec4(  1.0f,  1.0f,  1.0f, 1.0f ),
	vec4( -1.0f,  1.0f, -1.0f, 1.0f ),
	vec4( -1.0f,  1.0f,  1.0f, 1.0f ),
//
	vec4(  1.0f,  1.0f,  1.0f, 1.0f ),
	vec4( -1.0f,  1.0f,  1.0f, 1.0f ),
	vec4(  1.0f, -1.0f,  1.0f, 1.0f ),
};

static void MakeCubeMesh(Mesh * mesh)
{
	mesh->vertexCount = 8;
	mesh->vertices = (vec4 *)calloc(sizeof(vec4), mesh->vertexCount);
	mesh->vertices[0] = vec4( 1.0f,  1.0f,  1.0f, 1.0f);
	mesh->vertices[1] = vec4(-1.0f,  1.0f,  1.0f, 1.0f);
	mesh->vertices[2] = vec4(-1.0f, -1.0f,  1.0f, 1.0f);
	mesh->vertices[3] = vec4(-1.0f, -1.0f, -1.0f, 1.0f);
	mesh->vertices[4] = vec4( 1.0f, -1.0f, -1.0f, 1.0f);
	mesh->vertices[5] = vec4( 1.0f,  1.0f, -1.0f, 1.0f);
	mesh->vertices[6] = vec4(-1.0f,  1.0f, -1.0f, 1.0f);
	mesh->vertices[7] = vec4( 1.0f, -1.0f,  1.0f, 1.0f);

	mesh->indexCount = 6*2*3;  // Six sides, two triangles per side, 3 vertices each.
	mesh->indices = (uint *)calloc(sizeof(uint), mesh->indexCount);
	mesh->indices[0] = 3;
	mesh->indices[1] = 2;
	mesh->indices[2] = 1;

	mesh->indices[3] = 5;
	mesh->indices[4] = 3;
	mesh->indices[5] = 6;

	mesh->indices[6] = 7;
	mesh->indices[7] = 3;
	mesh->indices[8] = 4;

	mesh->indices[9]  = 5;
	mesh->indices[10] = 4;
	mesh->indices[11] = 3;

	mesh->indices[12] = 3;
	mesh->indices[13] = 1;
	mesh->indices[14] = 6;

	mesh->indices[15] = 7;
	mesh->indices[16] = 2;
	mesh->indices[17] = 3;

	mesh->indices[18] = 1;
	mesh->indices[19] = 2;
	mesh->indices[20] = 7;

	mesh->indices[21] = 0;
	mesh->indices[22] = 4;
	mesh->indices[23] = 5;

	mesh->indices[24] = 4;
	mesh->indices[25] = 0;
	mesh->indices[26] = 7;

	mesh->indices[27] = 0;
	mesh->indices[28] = 5;
	mesh->indices[29] = 6;

	mesh->indices[30] = 0;
	mesh->indices[31] = 6;
	mesh->indices[32] = 1;

	mesh->indices[33] = 0;
	mesh->indices[34] = 1;
	mesh->indices[35] = 7;
}

static void 
RenderMesh(win32_backbuffer * backbuffer, Mesh * mesh, mat4x4 transform)
{
	float r = (gScreenWidth / 2.0f);
	float l = -(gScreenWidth / 2.0f);
	float t = -(gScreenHeight / 2.0f);
	float b = (gScreenHeight / 2.0f);
	float n = -1000.0f;
	float f = 1000.0f;
	static const mat4x4 frustumMatrix = FrustumMatrix(r, l, t, b, n, f);

	assert(mesh->vertexCount > 0);
	vec4 * xformedVerts = (vec4 *)malloc(sizeof(vec4) * mesh->vertexCount);

	for (u32 i = 0; i < mesh->vertexCount; ++i)
	{
		xformedVerts[i] = frustumMatrix * transform * mesh->vertices[i];
	}

	assert(mesh->indexCount % 3 == 0);
	for (u32 i = 0; i < mesh->indexCount; i += 3)
	{
		uint idxA = mesh->indices[i];
		uint idxB = mesh->indices[i + 1];
		uint idxC = mesh->indices[i + 2];
		Rasterize(backbuffer, xformedVerts[idxA], xformedVerts[idxB], xformedVerts[idxC], faceColors[i % 6], faceColors[i % 6], faceColors[i % 6]);
	}
}

static void 
RenderTestMesh(win32_backbuffer * backbuffer, Mesh * mesh)
{
	static const vec3 scale = { 1.0f, 1.0f, 1.0f };
	static const vec3 position = { 0.0f, 0.0f, 0.0f };
	static const vec4 axis(0.0f, 1.0f, 1.0f, 0.0f);

	static float angle = 0.0f;
	quaternion rotation = RotationAroundAxis(angle, axis);
	angle += 0.1f;
	mat4x4 transform = MakeTransformMatrix(rotation, scale, position);

	static const vec3 cameraPosition = { 0.0f, 0.0f, 15.0f }; 
	quaternion cameraRotation = { 0.0f, 0.0f, 0.0f, 1.0f };

	transform = MakeTransformMatrix(cameraRotation, scale, cameraPosition) * transform;

	RenderMesh(backbuffer, mesh, transform);
}

static void RenderTest(win32_backbuffer * backbuffer)
{
	static const vec3 scale = { 1.0f, 1.0f, 1.0f };
	static const vec3 position = { 0.0f, 0.0f, 0.0f };

	static float angle = 0.0f;
	vec4 axis = vec4( 1.0f, 0.5f, 0.0f, 0.0f );
	quaternion rotation = RotationAroundAxis(angle, axis);
	mat4x4 transform = MakeTransformMatrix(rotation, scale, position);
	angle += 0.1f;

	static const vec3 scale2 = { 0.5f, 1.5f, 0.5f };
	static const vec3 position2 = { 1.0f, 1.0f, 5.0f };
	static const quaternion rot2 = RotationAroundAxis((float)M_PI_4, axis);
	static const mat4x4 transform2 = MakeTransformMatrix(rot2, scale2, position2);

	static float cameraAngle = 0.0f;
	vec4 upAxis = vec4( 0.0f, -1.0f, 0.0f, 0.0f );
	quaternion cameraRotation = RotationAroundAxis(cameraAngle, upAxis);
	vec3 cameraPosition = { 0.0f, 0.0f, 15.0f };
	mat4x4 cameraTransform = MakeTransformMatrix(cameraRotation, scale, cameraPosition);
	cameraAngle -= 0.05f;

	RenderMesh(backbuffer, &meshVerts[0], meshVertexCount, cameraTransform * transform);
	RenderMesh(backbuffer, &meshVerts[0], meshVertexCount, cameraTransform * transform2);

	for (int i = 0; i < 4; ++i)
	{
		quaternion q = {0, 0, 0, 1}; //RotationAroundAxis((float)M_PI_4 * i, upAxis);
		vec3 p = { -1.75, -1.75, 2.5f * (i + 1) };
		mat4x4 t = MakeTransformMatrix(q, scale, p);
		RenderMesh(backbuffer, &meshVerts[0], meshVertexCount, cameraTransform * t);
	}
}

static u64
Win32TimerFrequency()
{
	LARGE_INTEGER out;
	QueryPerformanceFrequency(&out);
	return out.QuadPart;
}

static u64
Win32GetPerformanceTimer()
{
	LARGE_INTEGER out;
	QueryPerformanceCounter(&out);
	return out.QuadPart;
}

int CALLBACK 
WinMain(HINSTANCE hInstance,
		HINSTANCE hPrevInstance,
		LPSTR lpCmdLine,
		int nCmdShow)
{
	WNDCLASSEX windowClass = { 0 };
	windowClass.cbSize = sizeof(WNDCLASSEX);
	windowClass.style = CS_HREDRAW | CS_VREDRAW | CS_OWNDC;
	windowClass.lpfnWndProc = Win32MainWindowProc;
	windowClass.hInstance = hInstance;
	windowClass.lpszClassName = TEXT("SoftwareRendererMainWindowClass");

	RegisterClassEx(&windowClass);

	RECT clientSize;
	clientSize.top = 0;
	clientSize.left = 0;
	clientSize.right = gScreenWidth;
	clientSize.bottom = gScreenHeight;

	AdjustWindowRect(&clientSize, WS_OVERLAPPEDWINDOW, FALSE);

	HWND window = CreateWindowEx(0, windowClass.lpszClassName, 
								 TEXT("Software Renderer"), 
								 WS_OVERLAPPEDWINDOW, 
								 CW_USEDEFAULT, CW_USEDEFAULT, 
								 clientSize.right - clientSize.left, 
								 clientSize.bottom - clientSize.top, 
								 NULL, NULL, hInstance, NULL);
	assert(window);

	ShowWindow(window, nCmdShow);

	win32_backbuffer backbuffer;
	Win32AllocateDIBSection(&backbuffer, gScreenWidth, gScreenHeight);

	gRunning = true;

	TestMatrixMultiply();

	u64 frameTargetMS = 1000 / 30;
	u64 perfTicksPerMS = Win32TimerFrequency() / 1000;
	u64 lastTick = Win32GetPerformanceTimer();

	Mesh cube;
	MakeCubeMesh(&cube);

	MSG msg;
	while (gRunning)
	{
		while (PeekMessage(&msg, window, 0, 0, 1))
		{
			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}

		RECT clientRect;
		GetClientRect(window, &clientRect);
		int x = clientRect.left;
		int y = clientRect.top;
		int width = clientRect.right - clientRect.left;
		int height = clientRect.bottom - clientRect.top;
		
//		RenderTest(&backbuffer);
		RenderTestMesh(&backbuffer, &cube);

		Win32RedrawWindow(window, x, y, width, height, &backbuffer);
		Win32ClearBackbuffer(&backbuffer, RGBA32(0, 0, 0, 0));

		u64 elapsedSinceFrameStart = Win32GetPerformanceTimer() - lastTick;
		float frameTime = (float)elapsedSinceFrameStart / (float)perfTicksPerMS;

		{
			wchar_t buffer[1024] = {0};
			swprintf(buffer, _countof(buffer), L"frame time: %fms\n", frameTime);
			OutputDebugStringW(buffer);
		}

		u64 msToSleep = frameTargetMS - (elapsedSinceFrameStart / perfTicksPerMS);

		if (msToSleep < frameTargetMS) Sleep((DWORD)msToSleep);

		lastTick = Win32GetPerformanceTimer();
	}

	return 0;
}