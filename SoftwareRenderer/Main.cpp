#include <Windows.h>
#include <cassert>

#define _USE_MATH_DEFINES 1
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cfloat>

#define STBI_ONLY_PNG
#include "stb/stb_image.h"

#include "types.h"

// For std::sort
#include <algorithm>
// For loader
#include <vector>

// Disable warnings about C runtime functions
// In a production setting this is a bad idea, but this is just proof-of-concept.
#pragma warning(disable: 4996)

static const int gScreenWidth = 1280;
static const int gScreenHeight = 720;
static bool gRunning;

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

int ReadColor(Color * color, const char * c)
{
	color->a = 1.0f; // Colors are given as RGB, no alpha channel.

	char * cur = (char *)c;

	color->r = (float)strtod(cur, &cur);
	color->g = (float)strtod(cur, &cur);
	color->b = (float)strtod(cur, &cur);

	return cur - c;
}

Texture * LoadTexture(char * filename)
{
	FILE * file = fopen(filename, "rb");
	if (file == NULL) 
	{
		return NULL;
	}
	int x;
	int y;
	int comp; 
	u8 * data = stbi_load_from_file(file, &x, &y, &comp, 4);

	if (data == NULL) {
		assert(0);
		return NULL;
	}

	Texture * texture = (Texture *)malloc(sizeof(Texture));
	texture->width = x;
	texture->height = y;
	u32 * texData = (u32 *)malloc(4 * x * y);

	u32 * pixels = (u32 *)data;
	for (int i = 0; i < x * y; ++i)
	{
		texData[i] = pixels[i];
	}

	texture->data = (char *)texData;

	stbi_image_free(data);

	return texture;
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

	char * bytes = (char *)calloc(1, fileLength + 1);
	fileLength = fread(bytes, 1, fileLength, file);
	bytes[fileLength] = '\0';

	// Clear material struct
	memset(material, 0, sizeof(Material));

	const char * cur = bytes;
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
			char * end;
			material->specularIntensity = (float)strtod(lineStart + 2, &end);
			assert(end != lineStart + 2);
		}
		else if (!strncmp(lineStart, "illum", 5))
		{
			char * end;
			material->illumType = (u8)strtol(lineStart + 5, &end, 10);
			assert(end != lineStart + 5);
		}
		else if (*lineStart == 'K')
		{
			lineStart++;
			if      (*lineStart == 'a') ReadColor(&material->ambientColor, lineStart + 1);
			else if (*lineStart == 'd') ReadColor(&material->diffuseColor, lineStart + 1);
			else if (*lineStart == 's') ReadColor(&material->specularColor, lineStart + 1);
		}
		else if (!strncmp(lineStart, "map_K", 5))
		{
			// TODO:  Texture map.
			lineStart += 5;
			char * fnStart = (char *)lineStart + 2;
			int fnLength = lineEnd - fnStart;
			char filename[256] = { 0 };
			memcpy(filename, fnStart, fnLength);

			if (*lineStart == 'a')
			{
				material->ambientTexture = LoadTexture(filename);
			}
			else if (*lineStart == 'd')
			{
				material->diffuseTexture = LoadTexture(filename);
			}
		}
	}

	free(bytes);
}

void FindNext(char c, const char * str, char ** out)
{
	char * cur = (char *)str;
	while (*cur && *cur != c) ++cur;
	*out = cur;
}

void LoadMesh(char * filename, Mesh * mesh)
{
	FILE * file = fopen(filename, "r");
	fseek(file, 0, SEEK_END);
	u32 fileLength = ftell(file);  // This is a minimum, may overestimate length of text files.
	fseek(file, 0, SEEK_SET);

	char * bytes = (char *)calloc(1, fileLength + 1);
	fileLength = fread(bytes, 1, fileLength, file);
	bytes[fileLength] = '\0';

	fclose(file);

	memset(mesh, 0, sizeof(Mesh));

	std::vector<vec4> vertices;
	std::vector<vec3> uvws;
	std::vector<uint> indices;

	const char * cur = bytes;
	const char * lineStart;
	const char * lineEnd;

	while (ReadLine(&lineStart, &lineEnd, &cur))
	{
		if (*lineStart == 'v' && *(lineStart + 1) == ' ') 
		{
			// Vertex position
			float x;
			float y;
			float z;
			char * cur = (char *)lineStart + 1;
			x = (float)strtod(cur, &cur);
			y = (float)strtod(cur, &cur);
			z = (float)strtod(cur, &cur);

			vertices.push_back(vec4(x, y, z, 1.0f));
		}
		else if (*lineStart == 'v' && *(lineStart + 1) == 't')
		{
			// Vertex uvw
			float u;
			float v;
			float w;

			char * cur = (char *) lineStart + 2;
			u = (float)strtod(cur, &cur);
			v = (float)strtod(cur, &cur);
			w = (float)strtod(cur, &cur);

			assert(cur <= lineEnd);
			uvws.push_back(vec3(u, v, w));
		}
		else if (*lineStart == 'f')
		{
			// Face 
			++lineStart;
			char * lineCur = (char *)lineStart;

			uint indexA = -1;
			uint uvIndexA = -1;

			uint indexB = -1;
			uint uvIndexB = -1;

			uint indexC = -1;
			uint uvIndexC = -1;
			
			uint indexD = -1;
			uint uvIndexD = -1;

			indexA = (uint)strtol(lineCur, &lineCur, 10);
			++lineCur;
			uvIndexA = (uint)strtol(lineCur, &lineCur, 10);
			FindNext(' ', lineCur, &lineCur);
			assert(lineCur < lineEnd - 1);
			indexB = (uint)strtol(lineCur, &lineCur, 10);
			++lineCur;
			uvIndexB = (uint)strtol(lineCur, &lineCur, 10);
			FindNext(' ', lineCur, &lineCur);
			assert(lineCur < lineEnd - 1);
			indexC = (uint)strtol(lineCur, &lineCur, 10);
			++lineCur;
			uvIndexC = (uint)strtol(lineCur, &lineCur, 10);
			FindNext(' ', lineCur, &lineCur);

			char * last;
			indexD = (uint)strtol(lineCur, &last, 10);
			++last;
			uvIndexD = (uint)strtol(last, &last, 10);
			if (last != lineCur && last <= lineEnd)
			{
				// Quad
				indices.push_back(indexA - 1);
				indices.push_back(uvIndexA - 1);
				indices.push_back(indexB - 1);
				indices.push_back(uvIndexB - 1);
				indices.push_back(indexC - 1);
				indices.push_back(uvIndexC - 1);
				
				indices.push_back(indexA - 1);
				indices.push_back(uvIndexA - 1);
				indices.push_back(indexC - 1);
				indices.push_back(uvIndexC - 1);
				indices.push_back(indexD - 1);
				indices.push_back(uvIndexD - 1);
			}
			else
			{
				// Triangle
				indices.push_back(indexA - 1);
				indices.push_back(uvIndexA - 1);
				indices.push_back(indexB - 1);
				indices.push_back(uvIndexB - 1);
				indices.push_back(indexC - 1);
				indices.push_back(uvIndexC - 1);
			}
		}
	}

	free(bytes);

	mesh->vertexCount = vertices.size();
	mesh->uvwCount = uvws.size();
	mesh->indexCount = indices.size();

	// Allocate mesh data in a single block
	u32 verticesSize = mesh->vertexCount * sizeof(vec4);
	u32 uvwsSize = mesh->uvwCount * sizeof(vec3);
	u32 indicesSize = mesh->indexCount * sizeof(uint);
	u32 totalSize = verticesSize + uvwsSize + indicesSize;
	char * meshBytes = (char *)malloc(totalSize);

	assert(meshBytes);

	mesh->vertices = (vec4 *)meshBytes;
	mesh->uvws = (vec3 *)(meshBytes + verticesSize); 
	mesh->indices = (uint *)(meshBytes + verticesSize + uvwsSize);

	memcpy(mesh->vertices, &vertices[0], verticesSize);
	memcpy(mesh->uvws, &uvws[0], uvwsSize);
	memcpy(mesh->indices, &indices[0], indicesSize);
}

static Color
SampleTexture2D(Texture * texture, float u, float v)
{
	int width = texture->width;
	int height = texture->height;
	int x = (int)(u * width);
	int y = (int)(v * height);

	u32 * pixels = (u32 *)texture->data;

	return Color(pixels[y * width + x]);
}

// Rasterizer.
//
// Based on the triangle rasterizer described by Fabian Giesen here:
// https://fgiesen.wordpress.com/2013/02/08/triangle-rasterization-in-practice/
//

static vec3
HalfVector(vec3 a, vec3 b)
{
	vec3 h(
		a.x + b.x,
		a.y + b.y,
		a.z + b.z
	);
	return normalized(h);
}

struct VertexAttributes {
	vec3 uvw;
	float depth;
};

// MSVC will gladly generate conditional moves (or branches) rather than use
// SSE2 min/max for floats.
//
// So, we use intrinsics to force it.
#include <mmintrin.h>
inline float
minf(const float& a, const float& b)
{
	float rv;
	_mm_store_ss(&rv, _mm_min_ss(_mm_set_ss(a), _mm_set_ss(b)));
	return rv;
}

inline float
maxf(const float& a, const float& b)
{
	float rv;
	_mm_store_ss(&rv, _mm_max_ss(_mm_set_ss(a), _mm_set_ss(b)));
	return rv;
}

inline float
sse_floorf(const float& a)
{
	float rv;
	_mm_store_ss(&rv, _mm_floor_ps(_mm_set_ss(a)));
	return rv;
}

inline float
sse_ceilf(const float& a)
{
	float rv;
	_mm_store_ss(&rv, _mm_ceil_ps(_mm_set_ss(a)));
	return rv;
}

#define floorf(a) sse_floorf(a)
#define ceilf(a) sse_ceilf(a)

static void
SetPixel(win32_backbuffer * backbuffer, Material * material, vec4 p, 
		 VertexAttributes va0, VertexAttributes va1, VertexAttributes va2, 
		 float l0, float l1, float l2, 
		 vec3 surfaceNormal, vec3 viewNormal)
{
	static const vec3 sunDir(0.0f, -1.0f, 1.0f);
	int width = backbuffer->bmpInfo.bmiHeader.biWidth;
	int height = backbuffer->bmpInfo.bmiHeader.biHeight;

	int x = (int)(p.x + 0.5f);
	int y = (int)(p.y + 0.5f);

	int idx = y * width + x;

	float depth = va0.depth + l1*(va1.depth-va0.depth) + l2*(va2.depth-va0.depth);
	if (depth > backbuffer->depthBuffer[idx])
	{
		u32 * pixels = (u32 *)backbuffer->bmpMemory;

		float u = va0.uvw.x + l1*(va1.uvw.x - va0.uvw.x) + l2*(va2.uvw.x - va0.uvw.x);
		float v = va0.uvw.y + l1*(va1.uvw.y - va0.uvw.y) + l2*(va2.uvw.y - va0.uvw.y);
#if 0
		vec3 sunNorm = normalized(sunDir);
		float surfaceNormDotSun = maxf(dot(surfaceNormal, sunNorm), 0.0f);
		Color diffuse = material->diffuseColor * SampleTexture2D(material->diffuseTexture, u, v) * surfaceNormDotSun;

		vec3 h = HalfVector(sunNorm, viewNormal);
		float surfaceNormDotH = maxf(dot(surfaceNormal, h), 0.0f);
		float intensity = pow(surfaceNormDotH, material->specularIntensity);
		Color specular = material->specularColor * intensity * (dot(sunNorm, surfaceNormal) > 0.0f);
		Color color = 
			material->ambientColor * 0.025f + 
			diffuse + specular;
#else
		Color color = material->diffuseColor * SampleTexture2D(material->diffuseTexture, u, v);
#endif
		pixels[idx] = color.rgba();
		backbuffer->depthBuffer[idx] = depth;
	}
}

#define min3f(a, b, c) minf(a, minf(b, c))
#define max3f(a, b, c) maxf(a, maxf(b, c))

// Compute (twice) the area of the triangle abc.
static float
Orient2D(const vec4& a, const vec4& b, const vec4& c)
{
	return (b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x);
}

static bool IsTopLeft(const vec4& a, const vec4& b)
{
	return (abs(a.y - b.y) < 0.5f) // Top edge
		|| (b.x > a.x);            // Left edge
}

static void
Rasterize(win32_backbuffer * backbuffer, Material * material, vec4 v0, vec4 v1, vec4 v2, vec3 uv0, vec3 uv1, vec3 uv2, vec3 normal)
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
	float minX = min3f(v0.x, v1.x, v2.x);
	float minY = min3f(v0.y, v1.y, v2.y);
	float maxX = max3f(v0.x, v1.x, v2.x);
	float maxY = max3f(v0.y, v1.y, v2.y);

	// Clip bounding box to screen.
	minX = floorf(maxf(minX, 0));
	minY = floorf(maxf(minY, 0));
	maxX = ceilf(minf(maxX, gScreenWidth - 1.0f));
	maxY = ceilf(minf(maxY, gScreenHeight - 1.0f));

	float stepSize = 1.0f; // Pixels to step in each direction.

	// Per-step deltas for barycentric weights.
	float a01 = (v0.y - v1.y) * stepSize;
	float b01 = (v1.x - v0.x) * stepSize;
	float a12 = (v1.y - v2.y) * stepSize;
	float b12 = (v2.x - v1.x) * stepSize;
	float a20 = (v2.y - v0.y) * stepSize;
	float b20 = (v0.x - v2.x) * stepSize;

	vec4 p( minX, minY, 0.0f, 1.0f );

	// Fill-rule bias
	float bias0 = IsTopLeft(v1, v2) ? 0.0f : -1.0f * stepSize;
	float bias1 = IsTopLeft(v2, v0) ? 0.0f : -1.0f * stepSize;
	float bias2 = IsTopLeft(v0, v1) ? 0.0f : -1.0f * stepSize;
	
	// Calculate barycentric coordinates for first pixel.
	float w0_row = Orient2D(v1, v2, p) - bias0;
	float w1_row = Orient2D(v2, v0, p) - bias1;
	float w2_row = Orient2D(v0, v1, p) - bias2;

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

	VertexAttributes va0 = { uv0, v0.w };
	VertexAttributes va1 = { uv1, v1.w };
	VertexAttributes va2 = { uv2, v2.w };

	for (p.y = minY; p.y < maxY; p.y += stepSize)
	{
		float l0 = lambda0_row;
		float l1 = lambda1_row;
		float l2 = lambda2_row;

		for (p.x = minX; p.x < maxX; p.x += stepSize)
		{
			// Test whether any of the barycentric weights are negative.
			// Using a bitmask here is a definite win on larger meshes.
			u32 mask = *(u32 *)&l0 | *(u32 *)&l1 | *(u32 *)&l2;
			if (~mask & 0x80000000)
			{
				vec3 viewNormal(p.x, p.y, -1.0f);
				SetPixel(backbuffer, material, p, va0, va1, va2, l0, l1, l2, normal, normalized(viewNormal));
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

vec3 SurfaceNormal(const vec4& a, const vec4& b, const vec4& c)
{
	vec4 v1 = b - a;
	vec4 v2 = c - a;

	vec3 norm(
		v1.y*v2.z - v1.z*v2.y,
		v1.z*v2.x - v1.x*v2.z,
		v1.x*v2.y - v1.y*v2.x
	);

	return normalized(norm);
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

	assert(mesh->indexCount % 6 == 0);
	for (u32 i = 0; i < mesh->indexCount; i += 6)
	{
		uint idxA = mesh->indices[i];
		uint uvwA = mesh->indices[i + 1];
		uint idxB = mesh->indices[i + 2];
		uint uvwB = mesh->indices[i + 3];
		uint idxC = mesh->indices[i + 4];
		uint uvwC = mesh->indices[i + 5];

		vec4 v0 = xformedVerts[idxA];
		vec4 v1 = xformedVerts[idxB];
		vec4 v2 = xformedVerts[idxC];

		vec3 normal = SurfaceNormal(v0, v1, v2);

		Rasterize(backbuffer, mesh->material,
			v0, v1, v2,
			mesh->uvws[uvwA], mesh->uvws[uvwB], mesh->uvws[uvwC],
			normal);
	}
	free(xformedVerts);
}

static void 
RenderTestMesh(win32_backbuffer * backbuffer, Mesh * mesh, vec3 cameraPosition, quaternion cameraRotation)
{
	static const vec3 scale(1.0f, 1.0f, 1.0f);
	static const vec3 position(0.0f, -15.0f, 0.0f);
	static const vec4 axis(0.0f, 1.0f, 0.0f, 0.0f);
	static const vec4 xAxis(-1.0f, 0.0f, 0.0f, 0.0f);

	static float angle = 0.0f;
	quaternion rotation = RotationAroundAxis(angle, axis);
	angle += 0.01f;
	mat4x4 transform = MakeTransformMatrix(rotation, scale, position);

//	static const vec3 cameraPosition(0.0f, 0.0f, 250.0f); 
//	quaternion cameraRotation = RotationAroundAxis((float)M_PI * -1.25f, xAxis);

	transform = MakeTransformMatrix(cameraRotation, scale, cameraPosition) * transform;

	RenderMesh(backbuffer, mesh, transform);
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

	SetCurrentDirectory(TEXT("..\\Resources\\teapot\\"));

	ShowWindow(window, nCmdShow);

	win32_backbuffer backbuffer;
	Win32AllocateDIBSection(&backbuffer, gScreenWidth, gScreenHeight);

	gRunning = true;

	u64 frameTargetMS = 1000 / 30;
	u64 perfTicksPerMS = Win32TimerFrequency() / 1000;
	u64 lastTick = Win32GetPerformanceTimer();

	Mesh teapot;
	LoadMesh("teapot.obj", &teapot);
	Material mat = { 0 };
	LoadMaterial("default.mtl", &mat);
	teapot.material = &mat;

	vec3 cameraPosition(0.0f, 0.0f, 250.0f); 
	quaternion cameraRotation = { 0.0f, 0.0f, 0.0f, 1.0f };

	MSG msg;
	while (gRunning)
	{
		while (PeekMessage(&msg, window, 0, 0, 1))
		{
			switch(msg.message)
			{
			case WM_KEYDOWN:
				switch (msg.wParam)
				{
				case 'A':
					cameraPosition.x += 1.0f;
					break;
				case 'D':
					cameraPosition.x -= 1.0f;
					break;
				case 'W':
					cameraPosition.z -= 1.0f;
					break;
				case 'S':
					cameraPosition.z += 1.0f;
					break;

				case 'Q':
					cameraRotation = RotationAroundAxis((float)0.25, vec4(0.0f, 1.0f, 0.0f, 0.0f)) * cameraRotation;
					break;
				case 'E':
					cameraRotation = RotationAroundAxis((float)-0.25, vec4(0.0f, 1.0f, 0.0f, 0.0f)) * cameraRotation;
					break;
				}
				break;
			}
			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}

		RECT clientRect;
		GetClientRect(window, &clientRect);
		int x = clientRect.left;
		int y = clientRect.top;
		int width = clientRect.right - clientRect.left;
		int height = clientRect.bottom - clientRect.top;
		
		RenderTestMesh(&backbuffer, &teapot, cameraPosition, cameraRotation);

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