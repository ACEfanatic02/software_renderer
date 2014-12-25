#include <Windows.h>
#include <cassert>

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
				uint r = (c0.r*w0 + c1.r*w1 + c2.r*w2) / totalWeight;
				uint g = (c0.g*w0 + c1.g*w1 + c2.g*w2) / totalWeight;
				uint b = (c0.b*w0 + c1.b*w1 + c2.b*w2) / totalWeight;
				uint a = (c0.a*w0 + c1.a*w1 + c2.a*w2) / totalWeight;

				SetPixel(backbuffer, p.x, p.y, RGBA32(r, g, b, a));
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

		RasterizeTriangle(&backbuffer, v0, v1, v2);
		RasterizeTriangle(&backbuffer, v0, v2, v3);

		RECT clientRect;
		GetClientRect(window, &clientRect);
		int x = clientRect.left;
		int y = clientRect.top;
		int width = clientRect.right - clientRect.left;
		int height = clientRect.bottom - clientRect.top;
		
		Win32RedrawWindow(window, x, y, width, height, &backbuffer);
		Win32ClearBackbuffer(&backbuffer, RGBA32(0, 0, 0, 0));

		u64 elapsedSinceFrameStart = Win32GetPerformanceTimer() - lastTick;
		float frameTime = (float)elapsedSinceFrameStart / (float)perfTicksPerMS;

		{
			wchar_t buffer[1024] = {0};
			swprintf(buffer, _countof(buffer), L"frame time: %fms\n", frameTime);
			OutputDebugStringW(buffer);
		}

//		u64 msToSleep = frameTargetMS - (elapsedSinceFrameStart / perfTicksPerMS);

//		assert(msToSleep < frameTargetMS);  // Check for underflow.

		Sleep(33);

		lastTick = Win32GetPerformanceTimer();
	}

	return 0;
}