#pragma once

#include <iostream>
#include <string>

using namespace std;
//************************************************
class C_rysowanie
{
public:
	C_rysowanie(const int wymiar_x, const int wymiar_y, double *tab, bool R, bool G, bool B);
	C_rysowanie();
	~C_rysowanie();
	//************************************************
	void wyrysuj(string filename);

private:
	double **tablica;
	const int rozmiar_tablicy_x, rozmiar_tablicy_y;
	int red, green, blue;
};