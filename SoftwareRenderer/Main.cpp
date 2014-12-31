#include <Windows.h>
#include <cassert>

#define _USE_MATH_DEFINES 1
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cfloat>

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

	char * cur = (char *)c;

	color->r = (float)strtod(cur, &cur);
	color->g = (float)strtod(cur, &cur);
	color->b = (float)strtod(cur, &cur);

	return cur - c;
}

Texture * LoadTexture(char * filename)
{
	return NULL;
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
			if      (*lineStart == 'a') ReadColor(&material->ambientColor, lineStart);
			else if (*lineStart == 'd') ReadColor(&material->diffuseColor, lineStart);
			else if (*lineStart == 's') ReadColor(&material->specularColor, lineStart);
		}
		else if (!strncmp(lineStart, "map_K", 5))
		{
			// TODO:  Texture map.
			lineStart += 5;
			char * fnStart = (char *)lineStart + 2;
			int fnLength = lineEnd - fnStart;
			char filename[256];
			memcpy(filename, fnStart, fnLength + 1);
			filename[fnLength] = '\0';

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
		}
		else if (*lineStart == 'f')
		{
			// Face 
			++lineStart;
			char * lineCur = (char *)lineStart;

			uint indexA = -1;
			uint indexB = -1;
			uint indexC = -1;
			uint indexD = -1;

			indexA = (uint)strtol(lineCur, &lineCur, 10);
			FindNext(' ', lineCur, &lineCur);
			assert(lineCur < lineEnd - 1);
			indexB = (uint)strtol(lineCur, &lineCur, 10);
			FindNext(' ', lineCur, &lineCur);
			assert(lineCur < lineEnd - 1);
			indexC = (uint)strtol(lineCur, &lineCur, 10);
			FindNext(' ', lineCur, &lineCur);

			char * last;
			indexD = (uint)strtol(lineCur, &last, 10);
			if (last != lineCur && last <= lineEnd)
			{
				// Quad
				indices.push_back(indexA - 1);
				indices.push_back(indexB - 1);
				indices.push_back(indexC - 1);
				
				indices.push_back(indexA - 1);
				indices.push_back(indexC - 1);
				indices.push_back(indexD - 1);
			}
			else
			{
				// Triangle
				indices.push_back(indexA - 1);
				indices.push_back(indexB - 1);
				indices.push_back(indexC - 1);
			}
		}
	}

	free(bytes);

	mesh->vertexCount = vertices.size();
	mesh->indexCount = indices.size();

	// Allocate mesh data in a single block
	u32 verticesSize = mesh->vertexCount * sizeof(vec4);
	u32 indicesSize = mesh->indexCount * sizeof(uint);
	u32 totalSize = verticesSize + indicesSize;
	char * meshBytes = (char *)malloc(totalSize);

	assert(meshBytes);

	mesh->vertices = (vec4 *)meshBytes;
	mesh->indices = (uint *)(meshBytes + verticesSize);

	memcpy(mesh->vertices, &vertices[0], verticesSize);
	memcpy(mesh->indices, &indices[0], indicesSize);
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
	minX = floorf(max(minX, 0));
	minY = floorf(max(minY, 0));
	maxX = ceilf(min(maxX, gScreenWidth - 1.0f));
	maxY = ceilf(min(maxY, gScreenHeight - 1.0f));

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
			// Test whether any of the barycentric weights are negative.
			// Using a bitmask here is a definite win on larger meshes.
			u32 mask = *(u32 *)&l0 | *(u32 *)&l1 | *(u32 *)&l2;
			if (~mask & 0x80000000)
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
		Rasterize(backbuffer, xformedVerts[idxA], xformedVerts[idxB], xformedVerts[idxC], //faceColors[0], faceColors[0], faceColors[0]);
			faceColors[(i) % 6], faceColors[(i + 1) % 6], faceColors[(i + 2) % 6]);
	}
}

static void 
RenderTestMesh(win32_backbuffer * backbuffer, Mesh * mesh)
{
	static const vec3 scale = { 1.0f, 1.0f, 1.0f };
	static const vec3 position = { 0.0f, -15.0f, 0.0f };
	static const vec4 axis(0.0f, 1.0f, 0.0f, 0.0f);
	static const vec4 xAxis(-1.0f, 0.0f, 0.0f, 0.0f);

	static float angle = 0.0f;
	quaternion rotation = RotationAroundAxis(angle, axis);
	angle += 0.01f;
	mat4x4 transform = MakeTransformMatrix(rotation, scale, position);

	static const vec3 cameraPosition = { 0.0f, 0.0f, 250.0f }; 
	quaternion cameraRotation = RotationAroundAxis((float)M_PI * -1.25f, xAxis);

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

	Mesh teapot;
	LoadMesh("C:\\Users\\Bryan Taylor\\Desktop\\teapot\\teapot.obj", &teapot);
//	LoadMesh("C:\\Users\\Bryan Taylor\\Desktop\\test.obj", &teapot);
	

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
		RenderTestMesh(&backbuffer, &teapot);

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