#include <Windows.h>
#include <cassert>

#define _USE_MATH_DEFINES 1
#include <cmath>
#include <cstdint>
#include <cstdio>

typedef unsigned int uint;
typedef uint8_t u8;
typedef uint32_t u32;
typedef uint64_t u64;

static const int gScreenWidth = 1280;
static const int gScreenHeight = 720;
static bool gRunning;

struct vec4
{
	float x;
	float y;
	float z;
	float w;
};

struct mat4x4
{
	float a0; float a1; float a2; float a3;
	float b0; float b1; float b2; float b3;
	float c0; float c1; float c2; float c3;
	float d0; float d1; float d2; float d3;
};

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
	float f = 1000.0f;
	mat4x4 test = FrustumMatrix(r, l, t, b, n, f);
	mat4x4 ident = { 0.0f };

	ident.a0 = 1.0f;
	ident.b1 = 1.0f;
	ident.c2 = 1.0f;
	ident.d3 = 1.0f;

	vec4 testVec = { -0.5f, -0.5f, 0.25f, 1.0f };

	OutputDebugStringW(L"Original:\n");

	DEBUGPrintMat4x4(test);

	OutputDebugStringW(L"\nAfter multiply by identity:\n");

	DEBUGPrintMat4x4(test * ident);

	OutputDebugStringW(L"\nTest vector:\n");

	DEBUGPrintVec4(testVec);

	OutputDebugStringW(L"\nTransformed vector:\n");

	DEBUGPrintVec4(test * testVec);

	OutputDebugStringW(L"\n90degree rotation matrix around Y:\n");

	mat4x4 rotMat = RotationMatrix(M_PI_2, RA_Y);
	DEBUGPrintMat4x4(rotMat);
	DEBUGPrintVec4(rotMat * testVec);

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
	u8 r;
	u8 g;
	u8 b;
	u8 a;

	Color(u8 _r, u8 _g, u8 _b, u8 _a) :
		r(_r),
		g(_g),
		b(_b),
		a(_a)
	{
	}

	Color(u32 rgba)
	{
		a = (rgba & 0xff000000) >> 24;
		r = (rgba & 0x00ff0000) >> 16;
		g = (rgba & 0x0000ff00) >> 8;
		b = (rgba & 0x000000ff);
	}
};

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

struct vec2i
{
	int x;
	int y;
};

struct vertex
{
	vec2i pos;
	Color color;
};

static void
SetPixel(win32_backbuffer * backbuffer, int x, int y, u32 color)
{
	int width = backbuffer->bmpInfo.bmiHeader.biWidth;
	int height = backbuffer->bmpInfo.bmiHeader.biHeight;

	assert(x >= 0 && x < width);
	assert(y >= 0 && y < height);

	u32 * pixels = (u32 *)backbuffer->bmpMemory;
	pixels[y * width + x] = color;
}

// Compute (twice) the area of the triangle abc.
static int
Orient2D(const vec2i& a, const vec2i& b, const vec2i& c)
{
	return (b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x);
}

static void
RasterizeTriangle(win32_backbuffer * backbuffer,
				  const vertex& v0, const vertex& v1, const vertex& v2)
{
	int minX = min(v0.pos.x, min(v1.pos.x, v2.pos.x));
	int minY = min(v0.pos.y, min(v1.pos.y, v2.pos.y));
	int maxX = max(v0.pos.x, max(v1.pos.x, v2.pos.x));
	int maxY = max(v0.pos.y, max(v1.pos.y, v2.pos.y));

	minX = max(0, minX);
	minY = max(0, minY);
	maxX = min(maxX, gScreenWidth - 1);
	maxY = min(maxY, gScreenHeight - 1);

	// Per-pixel step values for the barycentric weights.
	int a01 = v0.pos.y - v1.pos.y;
	int b01 = v1.pos.x - v0.pos.x;
	int a12 = v1.pos.y - v2.pos.y;
	int b12 = v2.pos.x - v1.pos.x;
	int a20 = v2.pos.y - v0.pos.y;
	int b20 = v0.pos.x - v2.pos.x;

	vec2i p = { minX, minY };

	int w0_row = Orient2D(v1.pos, v2.pos, p);
	int w1_row = Orient2D(v2.pos, v0.pos, p);
	int w2_row = Orient2D(v0.pos, v1.pos, p);

	for (p.y = minY; p.y <= maxY; ++p.y)
	{
		int w0 = w0_row;
		int w1 = w1_row;
		int w2 = w2_row;

		for (p.x = minX; p.x <= maxX; ++p.x)
		{
			// Calculate barycentric coordinates.
			// w0 + w1 + w2 = 2(Area of triangle v0v1v2)
//			int w0 = Orient2D(v1.pos, v2.pos, p);  // Among values in the triangle, this is largest when p = v0.pos
//			int w1 = Orient2D(v2.pos, v0.pos, p);  // Same for p = v1.pos
//			int w2 = Orient2D(v0.pos, v1.pos, p);  // Etc...

			// If any of w0, w1, w2 are negative, the point is outside the triangle.
			if ((w0 | w1 | w2) >= 0)
			{
				Color c0 = v0.color;
				Color c1 = v1.color;
				Color c2 = v2.color;

				// Interpolate colors.  Should work identically for UVs, depth, or any other linear interpolation.
				uint totalWeight = (w0 + w1 + w2);
				if (totalWeight > 0)
				{
					uint r = (c0.r*w0 + c1.r*w1 + c2.r*w2) / totalWeight;
					uint g = (c0.g*w0 + c1.g*w1 + c2.g*w2) / totalWeight;
					uint b = (c0.b*w0 + c1.b*w1 + c2.b*w2) / totalWeight;
					uint a = (c0.a*w0 + c1.a*w1 + c2.a*w2) / totalWeight;
				}

				SetPixel(backbuffer, p.x, p.y, RGBA32(0xff, 0xff, 0xff, 0xff));//RGBA32(r, g, b, a));
			}

			w0 += a12;
			w1 += a20;
			w2 += a01;
		}

		w0_row += b12;
		w1_row += b20;
		w2_row += b01;
	}
}

struct mat2x2
{
	float _00; float _01;
	float _10; float _11;
};

static void
RotateVec(vec2i * v, float angle)
{
	vec2i rv = *v;

	mat2x2 transform = {
		cos(angle), -sin(angle),
		sin(angle), cos(angle)
	};

	rv.x = (int)(v->x * transform._00 + v->y* transform._01);
	rv.y = (int)(v->x * transform._10 + v->y * transform._11);

	*v = rv;
}


static void
Rasterize(win32_backbuffer * backbuffer, const vec4& v0, const vec4& v1, const vec4& v2)
{
	vec2i a = { (int)((v0.x + 1.0f) * (gScreenWidth / 2.0f)), (int)((v0.y + 1.0f) * (gScreenHeight / 2.0f)) };
	vec2i b = { (int)((v1.x + 1.0f) * (gScreenWidth / 2.0f)), (int)((v1.y + 1.0f) * (gScreenHeight / 2.0f)) };
	vec2i c = { (int)((v2.x + 1.0f) * (gScreenWidth / 2.0f)), (int)((v2.y + 1.0f) * (gScreenHeight / 2.0f)) };

	vertex verts[3] = {
		{ a, Color(0xff, 0xff, 0xff, 0xff) },
		{ b, Color(0xff, 0xff, 0xff, 0xff) },
		{ c, Color(0xff, 0xff, 0xff, 0xff) },
	};

	RasterizeTriangle(backbuffer, verts[0], verts[1], verts[2]);
}

void RenderTest(win32_backbuffer * backbuffer)
{
	float r = (gScreenWidth / 2.0f);
	float l = -(gScreenWidth / 2.0f);
	float t = (gScreenHeight / 2.0f);
	float b = -(gScreenHeight / 2.0f);
	float n = 1.0f;
	float f = 100.0f;
	mat4x4 frustumMatrix = FrustumMatrix(r, l, t, b, n, f);

	static float angle = 0.0f;
	mat4x4 rotationMatrix = RotationMatrix(angle, RA_Y) * RotationMatrix(angle, RA_X);

	angle += 0.1f;

	vec4 verts[3] = {
		{ -100.0f, -100.0f, -100.0f, 1.0f },
		{ -100.0f, -100.0f,  100.0f, 1.0f },
		{ -100.0f,  100.0f,  100.0f, 1.0f },
/*		{  100.0f, -105.75f, -0.75f, 1.0f },
		{ -100.0f, -100.75f, -0.75f, 1.0f },
		{  100.0f, -100.75f, 10000.0f, 1.0f },*/
	};

	verts[0] = rotationMatrix * verts[0];
	verts[1] = rotationMatrix * verts[1];
	verts[2] = rotationMatrix * verts[2];

	vec4 transformedVerts[3];
	transformedVerts[0] = frustumMatrix * verts[0];
	transformedVerts[1] = frustumMatrix * verts[1];
	transformedVerts[2] = frustumMatrix * verts[2];

	Rasterize(backbuffer, transformedVerts[0], transformedVerts[1], transformedVerts[2]);
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

	HWND window = CreateWindowEx(0, windowClass.lpszClassName, 
								 TEXT("Software Renderer"), 
								 WS_OVERLAPPEDWINDOW, 
								 CW_USEDEFAULT, CW_USEDEFAULT, 
								 gScreenWidth, gScreenHeight, 
								 NULL, NULL, hInstance, NULL);
	assert(window);

	ShowWindow(window, nCmdShow);

	win32_backbuffer backbuffer;
	Win32AllocateDIBSection(&backbuffer, gScreenWidth, gScreenHeight);

	gRunning = true;

	vertex v0 = { { 0, gScreenHeight }, Color(0xff, 0, 0, 0xff) };
	vertex v1 = { { 0, 0 }, Color(0, 0, 0, 0xff) };
	vertex v2 = { { gScreenWidth, 0 }, Color(0xff, 0, 0, 0xff) };
	vertex v3 = { { gScreenWidth, gScreenHeight }, Color(0xff, 0xff, 0xff, 0xff) };

	
	TestMatrixMultiply();


	u64 frameTargetMS = 1000 / 15;
	u64 perfTicksPerMS = Win32TimerFrequency() / 1000;
	u64 lastTick = Win32GetPerformanceTimer();

	MSG msg;
	while (gRunning)
	{
		while (PeekMessage(&msg, window, 0, 0, 1))
		{
			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}

//		RasterizeTriangle(&backbuffer, v0, v1, v2);
//		RasterizeTriangle(&backbuffer, v0, v2, v3);

		RECT clientRect;
		GetClientRect(window, &clientRect);
		int x = clientRect.left;
		int y = clientRect.top;
		int width = clientRect.right - clientRect.left;
		int height = clientRect.bottom - clientRect.top;
		
		RenderTest(&backbuffer);

		Win32RedrawWindow(window, x, y, width, height, &backbuffer);
		Win32ClearBackbuffer(&backbuffer, RGBA32(0, 0, 0, 0));

		u64 elapsedSinceFrameStart = Win32GetPerformanceTimer() - lastTick;
		float frameTime = (float)elapsedSinceFrameStart / (float)perfTicksPerMS;

		{
			wchar_t buffer[1024] = {0};
			swprintf(buffer, _countof(buffer), L"frame time: %fms\n", frameTime);
//			OutputDebugStringW(buffer);
		}

//		u64 msToSleep = frameTargetMS - (elapsedSinceFrameStart / perfTicksPerMS);

//		assert(msToSleep < frameTargetMS);  // Check for underflow.

		Sleep(33);

		lastTick = Win32GetPerformanceTimer();
	}

	return 0;
}