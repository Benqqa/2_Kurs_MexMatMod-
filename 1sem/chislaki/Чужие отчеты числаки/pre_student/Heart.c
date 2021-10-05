#pragma warning(disable:4996)
#include "labengine.h"
#include <stdio.h>
#include <math.h>

typedef struct
{
	double x;
	double y;
} point_t;

typedef struct
{
	point_t xy1;
	point_t xy2;
} rect_t;

point_t Transform(point_t p, rect_t const* from, rect_t const* to) 
{
	point_t p2;
	p2.x = ((p.x - from->xy1.x)*(to->xy2.x - to->xy1.x)) / (from->xy2.x - from->xy1.x) + (to->xy1.x);
	p2.y = ((p.y - from->xy1.y)*(to->xy2.y - to->xy1.y)) / (from->xy2.y - from->xy1.y) + (to->xy1.y);
	return p2;
}

void DrawAxes(rect_t const* math, rect_t const* screen)
{
	int height, width;
	point_t point, p;
	height = LabGetHeight();
	width = LabGetWidth();
	p.x = 0.0;
	p.y = 0.0;
	LabSetColor(LABCOLOR_GREEN);
	point = Transform(p, math, screen);
	LabDrawLine(0, (int)point.y, width, (int)point.y);
	LabDrawLine((int)point.x, 0, (int)point.x, height);
	LabDrawFlush();
}

void DrawHeart(rect_t const* math, rect_t const* screen, double a)
{
	point_t point, p,p1;

	p.x = -1.81;
	while (p.x < 1.81)
	{	
		
		LabSetColorRGB((int)(fabs(p.x/3.6)*255), abs((int)(a/5 * 255)%510 - 255), (int)(fabs(1 - p.x / 3.6) * 255));			  //Default: LabSetColorRGB((int)(fabs(p.x/3.6)*255), abs((int)(a/5 * 255)%510 - 255), (int)(fabs(1 - p.x / 3.6) * 255));
		//LabSetColorRGB(r,g,b);
		//RGB == red, green, blue
		//r,g,b - integer,     0<=r,g,b<=255
		//purple - (128,0,128)
		//light-purple - (255,0,255)
		//yellow - (255,255,0)

		p.y = (pow(p.x, 2.0 / 3.0) + 0.9 * sqrt((3.3 - p.x * p.x)) *sin(a*3.1415* p.x) ); //производится вычисление следующей точки
		point = Transform(p, math, screen);
		LabDrawPoint((int)point.x, (int)point.y);
		p1.x = -p.x;
		p1.y = p.y;
		point = Transform(p1, math, screen);
		LabDrawPoint((int)point.x, (int)point.y);
		p.x += 0.007;
	}
	LabDrawFlush();
}


int main(void)
{
	int height, width, flag = 0;
	int i = 0, j = 0;
	double a=0;

	rect_t screen;
	rect_t math;

	if (LabInit())
	{
		{
			height = LabGetHeight();     //Здесь считывается размер окна и задаются координаты в системе экрана
			width = LabGetWidth();		 //В данном случае пересечение осей координат в левом верхнем углу окна,
			screen.xy1.x = 0;			 //а ось У направлена вниз. Ось Х осталась без изменений
			screen.xy1.y = 0;
			screen.xy2.x = width;
			screen.xy2.y = height;
		}
		
		{
			math.xy1.x = -5;   //Это математическая система координат. Точка пересечений осей лежит в центре экрана
			math.xy2.x = 5;    //и направлены они так, как вы уже привыкли: вправо для Х и вверх для У.
			math.xy1.y = 5;
			math.xy2.y = -5;
		}
    
		for (a = 0.001; a < 1000; a += 0.05)	 //Переменная а отвечает за плотность точек графика. Чем она меньше, тем ближе точки.
		{										 //Попробуйте ументьшить наращивание а с 0.05 до 0.01 и запустить программу.						
			LabDelay(20);
			DrawHeart(&math, &screen,a);
			LabClear();
		}

		LabInputKey();
		LabTerm();
	}
	return 0;
}

