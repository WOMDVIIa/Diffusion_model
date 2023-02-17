#include "C_rysowanie.h"
#include "bitmap_image.hpp"
//************************************************
C_rysowanie::C_rysowanie(const int x, const int y, double *tab, bool R, bool G, bool B) : rozmiar_tablicy_x(x), rozmiar_tablicy_y(y), red(0), green(0), blue(0)
{
	int licz = 0;
	tablica = new double *[rozmiar_tablicy_x];
	for (int i = 0; i < rozmiar_tablicy_x; i++)
	{
		tablica[i] = new double[rozmiar_tablicy_y];
		for (int j = 0; j < rozmiar_tablicy_y; j++)
		{
			tablica[i][j] = *(tab + licz);
			licz++;
		}
	}
	if (R) red = 255;
	if (G) green = 255;
	if (B) blue = 255;
}
//************************************************
C_rysowanie::C_rysowanie()	: rozmiar_tablicy_x(0), rozmiar_tablicy_y(0)
{
	cout << "Konstruktor domniemany klasy C_rysowanie" << endl;

	tablica = new double *[rozmiar_tablicy_x];
	for (int i = 0; i < rozmiar_tablicy_x; i++)
	{
		tablica[i] = new double[rozmiar_tablicy_y];
		for (int j = 0; j < rozmiar_tablicy_y; j++)
		{
			tablica[i][j] = NULL;
		}
	}
}
//************************************************
void C_rysowanie::wyrysuj(string filename)
{
	double boundary_min = 1.5;
	double boundary_max = 3.5;

	double colour;

	bitmap_image image(rozmiar_tablicy_x, rozmiar_tablicy_y);
	image_drawer draw(image);

	for (int i = 0; i < rozmiar_tablicy_x; i++)
	{
		for (int j = 0; j < rozmiar_tablicy_y; j++)
		{
			if (tablica[i][j] > boundary_min)
			{
				if (tablica[i][j] < boundary_max)
				{
					colour = 255 * (tablica[i][j] - boundary_min) / (boundary_max - boundary_min);
					draw.pen_color(255, colour, colour);
				}
				else
				{
					draw.pen_color(255, 255, 255); //bia³y
				}
			}
			else
			{
				draw.pen_color(red, green, blue);
			}

			draw.plot_pixel(i, j);
		}
	}

	image.save_image(filename);
}
//************************************************
C_rysowanie::~C_rysowanie()
{
	for (int i = 0; i < rozmiar_tablicy_x; i++)
		delete [] tablica[i];
	delete [] tablica;
}