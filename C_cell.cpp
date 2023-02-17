#include "C_cell.h"
#include <cmath>

const double	C_cell::dx						= 0.5e-5 /*Fe-Al-Ni*/ /*1.0e-3*/ /*Pb_bi*/ /*1.0e-4*/ /*Fe-Mn JPN*/ /*1.16e-3*/ /*Moja próbka*/ /*1.0e-6*/ /*Si or Ni*/  /*7.5e-4*/ /*Al-Cu*/ /*1e-9*/,				//	[m] // 2021-05-13:		1e-10
				C_cell::volume					= dx * dx * dx,			//	[m3]
				//C_cell::density					= 7.2e3,				//	[kg / m3]
				C_cell::density = /*7.515270e3,*/ /*Fe-C-Si*/ 7.469350e3, /*Fe-Mn*/  /*7.76136e3,*/ /*Fe-Ni*/ /*10.0e3,*/ /*Pb-Bi*/ /*2.9912e3,*/ /*Al-Cu*/				//	[kg / m3]
				C_cell::nominal_mass			= C_cell::density * C_cell::volume,
				C_cell::fraction				= 0.8,
				C_cell::tan_alpha_coefficient	= 5 * 3.141592654 / 180,
				C_cell::liq_fcc_solid_fraction_on_edge[2] = {0.0, 1.0};
//const int		C_cell::tab_calculations_gaps[fixed_calculations_distance] = {5, 10, 20, 40, 80};
//const int		C_cell::tab_calculations_gaps[fixed_calculations_distance] = {5, 10, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200};		// dzia³a³o dobrze
const int		C_cell::tab_calculations_gaps[fixed_calculations_distance] = {5, 10, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500};
bool			C_cell::flag = true;
double			C_cell::tab_NEWS_int_length[angle_count][fill_count][NEWS_plus_interface_length];
//*********************************
C_cell::C_cell()
{
	element_mass = new double[no_of_phases][no_of_elements];
	element_chemical_potential = new double[no_of_phases][no_of_elements];
	element_concentration = new double[no_of_phases][no_of_elements];
	if (flag)
	{
		fill_NEWS_int_length_table();
		flag = false;
	}
	//cout << "Konstruktor domniemany klasy C_cell\n";
}
//*********************************
C_cell::C_cell(const double *tab_conc, list_of_states *state, double *temper)
{
	element_mass = new double[no_of_phases][no_of_elements];
	element_chemical_potential = new double[no_of_phases][no_of_elements];
	element_concentration = new double[no_of_phases][no_of_elements];

	set_initial_conditions(tab_conc, state, temper);
}
//*********************************
C_cell::~C_cell()
{
	delete [] element_concentration;
	delete [] element_chemical_potential;
	delete [] element_mass;
}
//*********************************
inline void C_cell::set_state(list_of_states *state)
{
	cell_state = *state;
}
//*********************************
inline double C_cell::hypotenuse(double *a, double *b)
{
	return sqrt(*a * *a + *b * *b);
}
//*********************************
inline double C_cell::hypotenuse(double *a)
{
	return sqrt(*a * *a + 1);
}
//*********************************
double C_cell::mass_percent_to_kilogram(const double* mass_percent)
{
	return density * volume * ((*mass_percent) / 100);
}
//*********************************
double C_cell::mass_percent_to_kilogram(const double* mass_percent, double *cell_mass)
{
	return *cell_mass * ((*mass_percent) / 100);
}
//*********************************
void C_cell::set_initial_conditions(const double *tab_conc, list_of_states *state, double *temper)
{
	temperature = *temper;				// [K]
	set_state(state);

	if (cell_state == FCC)
	{
		element_mass[cell_state][no_of_elements - 1] = density * volume;			// masa pierwiastka bazowego wyliczona jako masa ca³ej komórki - nastêpnie odejmowane bêd¹ masy pierwiastków stopowych
		for (int i = 0; i < (no_of_elements - 1); i++)
		{
			element_mass[cell_state][i] = mass_percent_to_kilogram(tab_conc + i);
			element_mass[cell_state][no_of_elements - 1] -= element_mass[cell_state][i];		// tu odejmowane s¹ masy pierwiastków stopowych od masy pierwiastka bazowego
			element_mass[LIQUID][i] = 0;
		}
		element_mass[LIQUID][no_of_elements - 1] = 0;
	}
	else if (cell_state == LIQUID)
	{
		element_mass[cell_state][no_of_elements - 1] = density * volume;
		for (int i = 0; i < (no_of_elements - 1); i++)
		{
			element_mass[cell_state][i] = mass_percent_to_kilogram(tab_conc + i + no_of_elements - 1);
			element_mass[cell_state][no_of_elements - 1] -= element_mass[cell_state][i];
			element_mass[FCC][i] = 0;
		}
		element_mass[FCC][no_of_elements - 1] = 0;
	}
	else
	{
		element_mass[FCC][no_of_elements - 1] = density * volume * fraction;
		element_mass[LIQUID][no_of_elements - 1] = density * volume * (1 - fraction);

		for (int i = 0; i < (no_of_elements - 1); i++)
		{
			element_mass[FCC][i] = mass_percent_to_kilogram(tab_conc + i) * fraction;
			element_mass[FCC][no_of_elements - 1] -= element_mass[FCC][i];

			element_mass[LIQUID][i] = mass_percent_to_kilogram(tab_conc + i + no_of_elements - 1) * (1 - fraction);
			element_mass[LIQUID][no_of_elements - 1] -= element_mass[LIQUID][i];
		}
	}
}
//*********************************
void C_cell::set_initial_conditions(const double *tab_conc, list_of_states *state, double *solid_fraction, double *temper)
{
	temperature = *temper;				// [K]
	set_state(state);

	if (cell_state == FCC)
	{
		element_mass[cell_state][no_of_elements - 1] = density * volume;			// masa pierwiastka bazowego wyliczona jako masa ca³ej komórki - nastêpnie odejmowane bêd¹ masy pierwiastków stopowych
		for (int i = 0; i < (no_of_elements - 1); i++)
		{
			element_mass[cell_state][i] = mass_percent_to_kilogram(tab_conc + i);
			element_mass[cell_state][no_of_elements - 1] -= element_mass[cell_state][i];		// tu odejmowane s¹ masy pierwiastków stopowych od masy pierwiastka bazowego
			element_mass[LIQUID][i] = 0;
		}
		element_mass[LIQUID][no_of_elements - 1] = 0;
	}
	else if (cell_state == LIQUID)
	{
		element_mass[cell_state][no_of_elements - 1] = density * volume;
		for (int i = 0; i < (no_of_elements - 1); i++)
		{
			element_mass[cell_state][i] = mass_percent_to_kilogram(tab_conc + i + no_of_elements - 1);
			element_mass[cell_state][no_of_elements - 1] -= element_mass[cell_state][i];
			element_mass[FCC][i] = 0;
		}
		element_mass[FCC][no_of_elements - 1] = 0;
	}
	else
	{
		element_mass[FCC][no_of_elements - 1] = density * volume * *solid_fraction;
		element_mass[LIQUID][no_of_elements - 1] = density * volume * (1 - *solid_fraction);

		for (int i = 0; i < (no_of_elements - 1); i++)
		{
			element_mass[FCC][i] = mass_percent_to_kilogram(tab_conc + i) * *solid_fraction;
			element_mass[FCC][no_of_elements - 1] -= element_mass[FCC][i];

			element_mass[LIQUID][i] = mass_percent_to_kilogram(tab_conc + i + no_of_elements - 1) * (1 - *solid_fraction);
			element_mass[LIQUID][no_of_elements - 1] -= element_mass[LIQUID][i];
		}
	}

	//for (int i = 0; i < 2; i++)
	//{
	//	for (int j = 0; j < 3; j++)
	//	{
	//		cout << "\nfaza_" << i << "_ele_" << j << "_masa = " << element_mass[i][j] << endl;
	//	}
	//}
	//cout << endl;
}
//*********************************
void C_cell::set_initial_conditions_mass_from_file(const double *tab_conc, list_of_states *state, double *solid_fraction, double *temper, double *cell_mass)
{
	temperature = *temper;				// [K]
	set_state(state);

	if (cell_state == FCC)
	{
		element_mass[cell_state][no_of_elements - 1] = *cell_mass;			// masa pierwiastka bazowego wyliczona jako masa ca³ej komórki - nastêpnie odejmowane bêd¹ masy pierwiastków stopowych
		for (int i = 0; i < (no_of_elements - 1); i++)
		{
			element_mass[cell_state][i] = mass_percent_to_kilogram(tab_conc + i, cell_mass);
			element_mass[cell_state][no_of_elements - 1] -= element_mass[cell_state][i];		// tu odejmowane s¹ masy pierwiastków stopowych od masy pierwiastka bazowego
			element_mass[LIQUID][i] = 0;
		}
		element_mass[LIQUID][no_of_elements - 1] = 0;
	}
	else if (cell_state == LIQUID)
	{
		element_mass[cell_state][no_of_elements - 1] = *cell_mass;
		for (int i = 0; i < (no_of_elements - 1); i++)
		{
			element_mass[cell_state][i] = mass_percent_to_kilogram(tab_conc + i + no_of_elements - 1, cell_mass);
			element_mass[cell_state][no_of_elements - 1] -= element_mass[cell_state][i];
			element_mass[FCC][i] = 0;
		}
		element_mass[FCC][no_of_elements - 1] = 0;
	}
	else
	{
		element_mass[FCC][no_of_elements - 1] = *cell_mass * *solid_fraction;
		//cout << "\nfcc_fe_masa = " << element_mass[FCC][no_of_elements - 1] << endl;
		element_mass[LIQUID][no_of_elements - 1] = *cell_mass * (1 - *solid_fraction);
		//cout << "\nliq_fe_masa = " << element_mass[LIQUID][no_of_elements - 1] << endl;

		for (int i = 0; i < (no_of_elements - 1); i++)
		{
			element_mass[FCC][i] = mass_percent_to_kilogram(tab_conc + i, cell_mass) * *solid_fraction;
			//cout << "\nfcc_element_masa = " << element_mass[FCC][i] << endl;
			element_mass[FCC][no_of_elements - 1] -= element_mass[FCC][i];

			element_mass[LIQUID][i] = mass_percent_to_kilogram(tab_conc + i + no_of_elements - 1, cell_mass) * (1 - *solid_fraction);
			//cout << "\nliq_fe_masa = " << element_mass[LIQUID][i] << endl;
			element_mass[LIQUID][no_of_elements - 1] -= element_mass[LIQUID][i];
		}
	}

	//for (int i = 0; i < 2; i++)
	//{
	//	for (int j = 0; j < 3; j++)
	//	{
	//		cout << "\nfaza_" << i << "_ele_" << j << "_masa = " << element_mass[i][j] << endl;
	//	}
	//}
	//cout << endl;
}
//*********************************
void C_cell::decrease_countdowns()
{
	if (immunity_countdown > 0)
	{
		immunity_countdown--;
	}
	//*********************************
	//if (immunity_countdown > 175)
	//{
	//	calculations_countdown -= 5;
	//}
	//else if (immunity_countdown > 0)
	//{
	//	calculations_countdown -= 50;
	//}
	//else
	//{
	calculations_countdown--;
	//}
}
//*********************************
void C_cell::calculate_phase_fraction_and_masses()
{
	if (cell_state == FCC)
	{
		phase_fraction	= 1.0;
		mass[FCC]		= mass[no_of_phases] = calculate_phase_mass(&element_mass[FCC][0]);
		mass[LIQUID]	= 0.0;
	}
	else if (cell_state == LIQUID)
	{
		phase_fraction	= 0.0;
		mass[LIQUID]	= mass[no_of_phases] = calculate_phase_mass(&element_mass[LIQUID][0]);
		mass[FCC]		= 0.0;
	}
	else
	{
		mass[FCC]			= calculate_phase_mass(&element_mass[FCC][0]);
		mass[LIQUID]		= calculate_phase_mass(&element_mass[LIQUID][0]);
		mass[no_of_phases]	= mass[FCC] + mass[LIQUID];
		
		phase_fraction		= mass[FCC] / mass[no_of_phases];
	}

	decrease_countdowns();

	//cout << "w klasie, masa = " << mass[no_of_phases] << endl;
}
//*********************************
double C_cell::calculate_phase_mass(double *tab_ele_mass)
{
	double temp_phase_mass = 0.0;

	for (int i = 0; i < no_of_elements; i++)
	{
		temp_phase_mass += tab_ele_mass[i];
	}

	return temp_phase_mass;
}
//*********************************
void C_cell::calculate_concentration_and_chemical_potential_coefficient()
{
	//mass[no_of_phases] = 0.0;

	if (cell_state != INTERFACE)
	{
		for (int j = 0; j < no_of_phases; j++)
		{
			for (int i = 0; i < no_of_elements; i++)
			{
				element_concentration[j][i] = element_mass[j][i] / volume;
				//mass[no_of_phases] += element_mass[j][i];
				//cout << "concentration, phase " << j << ", element " << i << " = " << element_concentration[j][i] << endl;
			}
		}
		//cout << endl;
	}
	else
	{
		for (int i = 0; i < no_of_elements; i++)
		{
			element_concentration[FCC][i] = element_mass[FCC][i] / volume / phase_fraction;
			//mass[no_of_phases] += element_mass[FCC][i];
			//cout << "concentration, phase FCC" << ", element " << i << " = " << element_concentration[FCC][i] << endl;
		}

		for (int i = 0; i < no_of_elements; i++)
		{
			element_concentration[LIQUID][i] = element_mass[LIQUID][i] / volume / (1 - phase_fraction);
			//mass[no_of_phases] += element_mass[LIQUID][i];
			//cout << "concentration, phase LIQ" << ", element " << i << " = " << element_concentration[LIQUID][i] << endl;
		}
		//cout << endl;
	}

	calculate_chemical_potential_coefficient();
}
//*********************************
void C_cell::calculate_chemical_potential_coefficient()
{
	//chemical_potential_coefficient = 100000 * (mass[no_of_phases] / nominal_mass - 1) * (mass[no_of_phases] / nominal_mass - 1) * (mass[no_of_phases] / nominal_mass - 1);
	//chemical_potential_coefficient = 10e5 * pow((mass[no_of_phases] / nominal_mass - 1) * 40, 3);
	chemical_potential_coefficient = old_chemical_potential_coefficient = mass[no_of_phases] / nominal_mass - 1;
	
	
	
	if (mass[no_of_phases] > nominal_mass)
	{
		// 2021-03-25 po pierwszych obliczeniach: mo¿liwe, ¿e za bardzo ta zmiana wypacza równowagê 
		//chemical_potential_coefficient = 1 / ( ( (mass[no_of_phases] / nominal_mass) - 1 ) * 10 + 1); 
		// 2021-04-29 by³o: chemical_potential_coefficient = 1 / (((mass[no_of_phases] / nominal_mass) - 1) * 800 + 1); ale jest fala na masie komórki - niestabilnie
		// 2021-09-10 Ÿle liczy w poruwnaniu do artyku³u
		//chemical_potential_coefficient = 1 / (((mass[no_of_phases] / nominal_mass) - 1) * 3000 + 1);		// 2021-05-13  3000 
		//chemical_potential_coefficient = (mass[no_of_phases] / nominal_mass) - 1;
		//chemical_potential_coefficient = 1;
 	   //chemical_potential_coefficient = mass[no_of_phases] - nominal_mass;
  		//chemical_potential_coefficient = 10e5 * pow((mass[no_of_phases] / nominal_mass - 1) * 2, 3);
	}
	else
	{
		//chemical_potential_coefficient = (1 - (mass[no_of_phases] / nominal_mass) ) * 10 + 1;
		//chemical_potential_coefficient = (1 - (mass[no_of_phases] / nominal_mass)) * 3000 + 1;				// 2021-05-13  3000 
		//chemical_potential_coefficient = (mass[no_of_phases] / nominal_mass) - 1;				// 2021-09-10
		//chemical_potential_coefficient = - (mass[no_of_phases] - nominal_mass) * (mass[no_of_phases] - nominal_mass);
		//chemical_potential_coefficient = 0;
	}
}
//*********************************
void C_cell::set_NEWS_int_lenght_pointer()
{
	int temp_normal = interface_normal / 5 - 1,
		temp_fraction = round(100 * phase_fraction) - 1;

	if (temp_normal != 17)
	{
		switch (origin_corner)
		{
		case NE:																						//				 __________
			pointer_solid_fraction_on_edge[0] = &tab_NEWS_int_length[temp_normal][temp_fraction][3];	// N = SW(S)	|	   \   |
			pointer_solid_fraction_on_edge[1] = &tab_NEWS_int_length[temp_normal][temp_fraction][2];	// E = SW(W)	|		\  |
			pointer_solid_fraction_on_edge[2] = &tab_NEWS_int_length[temp_normal][temp_fraction][1];	// W = SW(E)	|		 \ |
			pointer_solid_fraction_on_edge[3] = &tab_NEWS_int_length[temp_normal][temp_fraction][0];	// S = SW(N)	|_________\|
			break;

		case NW:																						//				 __________
			pointer_solid_fraction_on_edge[0] = &tab_NEWS_int_length[temp_normal][temp_fraction][3];	// N = SW(S)	|	/	   |
			pointer_solid_fraction_on_edge[1] = &tab_NEWS_int_length[temp_normal][temp_fraction][1];	// E = SW(E)	|  /	   |
			pointer_solid_fraction_on_edge[2] = &tab_NEWS_int_length[temp_normal][temp_fraction][2];	// W = SW(W)	| /		   |
			pointer_solid_fraction_on_edge[3] = &tab_NEWS_int_length[temp_normal][temp_fraction][0];	// S = SW(N)	|/_________|
			break;

		case SE:																						//				 __________
			pointer_solid_fraction_on_edge[0] = &tab_NEWS_int_length[temp_normal][temp_fraction][0];	// N = SW(N)	|	      /|
			pointer_solid_fraction_on_edge[1] = &tab_NEWS_int_length[temp_normal][temp_fraction][2];	// E = SW(W)	|		 / |
			pointer_solid_fraction_on_edge[2] = &tab_NEWS_int_length[temp_normal][temp_fraction][1];	// W = SW(E)	|	    /  |
			pointer_solid_fraction_on_edge[3] = &tab_NEWS_int_length[temp_normal][temp_fraction][3];	// S = SW(S)	|______/___|
			break;

		default: //case SW																								 __________
			pointer_solid_fraction_on_edge[0] = &tab_NEWS_int_length[temp_normal][temp_fraction][0];	// N = SW(N)	|\		   |
			pointer_solid_fraction_on_edge[1] = &tab_NEWS_int_length[temp_normal][temp_fraction][1];	// E = SW(E)	| \		   |
			pointer_solid_fraction_on_edge[2] = &tab_NEWS_int_length[temp_normal][temp_fraction][2];	// W = SW(W)	|  \	   |
			pointer_solid_fraction_on_edge[3] = &tab_NEWS_int_length[temp_normal][temp_fraction][3];	// S = SW(S)	|___\______|
		}
	}
	else // temp_normal = 17 => interface_normal = 90
	{
		switch (origin_corner)
		{
		case NE:																						//				 _______o
			pointer_solid_fraction_on_edge[0] = &tab_NEWS_int_length[temp_normal][temp_fraction][3];	// N = SW(S)	|_______|
			pointer_solid_fraction_on_edge[1] = &tab_NEWS_int_length[temp_normal][temp_fraction][2];	// E = SW(W)	|	  	|
			pointer_solid_fraction_on_edge[2] = &tab_NEWS_int_length[temp_normal][temp_fraction][1];	// W = SW(E)	|_______|
			pointer_solid_fraction_on_edge[3] = &tab_NEWS_int_length[temp_normal][temp_fraction][0];	// S = SW(N)
			break;

		case NW:																						//				o_______
			pointer_solid_fraction_on_edge[0] = &tab_NEWS_int_length[temp_normal][temp_fraction][2];	// N = SW(W)	| |		|
			pointer_solid_fraction_on_edge[1] = &tab_NEWS_int_length[temp_normal][temp_fraction][0];	// E = SW(N)	| |		|
			pointer_solid_fraction_on_edge[2] = &tab_NEWS_int_length[temp_normal][temp_fraction][3];	// W = SW(S)	|_|_____|
			pointer_solid_fraction_on_edge[3] = &tab_NEWS_int_length[temp_normal][temp_fraction][1];	// S = SW(E)
			break;

		case SE:																						//				 _______
			pointer_solid_fraction_on_edge[0] = &tab_NEWS_int_length[temp_normal][temp_fraction][1];	// N = SW(E)	|	  |	|
			pointer_solid_fraction_on_edge[1] = &tab_NEWS_int_length[temp_normal][temp_fraction][3];	// E = SW(S)	|	  |	|
			pointer_solid_fraction_on_edge[2] = &tab_NEWS_int_length[temp_normal][temp_fraction][0];	// W = SW(N)	|_____|_|o
			pointer_solid_fraction_on_edge[3] = &tab_NEWS_int_length[temp_normal][temp_fraction][2];	// S = SW(W)
			break;

		default: //case SW																								  _______
			pointer_solid_fraction_on_edge[0] = &tab_NEWS_int_length[temp_normal][temp_fraction][0];	// N = SW(N)	 |		 |
			pointer_solid_fraction_on_edge[1] = &tab_NEWS_int_length[temp_normal][temp_fraction][1];	// E = SW(E)	 |_______|
			pointer_solid_fraction_on_edge[2] = &tab_NEWS_int_length[temp_normal][temp_fraction][2];	// W = SW(W)	o|_______|
			pointer_solid_fraction_on_edge[3] = &tab_NEWS_int_length[temp_normal][temp_fraction][3];	// S = SW(S)	
		}
	}

	pointer_interface_length = &tab_NEWS_int_length[temp_normal][temp_fraction][4];
}
//*********************************
void C_cell::fill_NEWS_int_length_table()
{
	double tan_alpha = 0;
	double a, b = 0;

	for (int i = 0; i < angle_count_minus_one; i++)
	{
				tan_alpha = tan( (1 + i) * tan_alpha_coefficient);
				tan_alpha = round(tan_alpha * 1e7) * 1e-7;
				//cout << "tangens(" << (i + 1) * 5 << ") = " << tan_alpha << "\n\n";

		for (int j = 0; j < half_fill_count; j++)
		{
			b = sqrt(2 * (j + 1) / tan_alpha / 100);

			if (b > 1)
			{
				a = tan_alpha;
				b = (j + 1) * 0.01 - a / 2;

				tab_NEWS_int_length[i][j][0] = b;						// N
				tab_NEWS_int_length[i][j][1] = 0;						// E
				tab_NEWS_int_length[i][j][2] = 1;						// W
				tab_NEWS_int_length[i][j][3] = a + b;					// S
				tab_NEWS_int_length[i][j][4] = hypotenuse(&a);			// interface length

				tab_NEWS_int_length[i][98 - j][0] = 1 - (a + b);		// N
				tab_NEWS_int_length[i][98 - j][1] = 0;					// E (1 - 1)
				tab_NEWS_int_length[i][98 - j][2] = 1;					// W (1 - 0)
				tab_NEWS_int_length[i][98 - j][3] = 1 - b;				// S
				tab_NEWS_int_length[i][98 - j][4] = hypotenuse(&a);	// interface length
			}
			else
			{
				a = b * tan_alpha;

				if (a > 1)
				{
					a = 1 / tan_alpha;
					b = (j + 1) * 0.01 - a / 2;

					tab_NEWS_int_length[i][j][0] = 0;						// N
					tab_NEWS_int_length[i][j][1] = b;						// E
					tab_NEWS_int_length[i][j][2] = a + b;					// W
					tab_NEWS_int_length[i][j][3] = 1;						// S
					tab_NEWS_int_length[i][j][4] = hypotenuse(&a);			// interface length

					tab_NEWS_int_length[i][98 - j][0] = 0;					// N (1 - 1)
					tab_NEWS_int_length[i][98 - j][1] = 1 - (a + b);		// E
					tab_NEWS_int_length[i][98 - j][2] = 1 - b;				// W
					tab_NEWS_int_length[i][98 - j][3] = 1;					// S (1 - 0)
					tab_NEWS_int_length[i][98 - j][4] = hypotenuse(&a);	// interface length
				}
				else
				{
					// a, b < 1
					tab_NEWS_int_length[i][j][0] = 0;							// N
					tab_NEWS_int_length[i][j][1] = 0;							// E
					tab_NEWS_int_length[i][j][2] = b;							// W
					tab_NEWS_int_length[i][j][3] = a;							// S
					tab_NEWS_int_length[i][j][4] = hypotenuse(&a, &b);			// interface length

					tab_NEWS_int_length[i][98 - j][0] = 1 - a;					// N
					tab_NEWS_int_length[i][98 - j][1] = 1 - b;					// E
					tab_NEWS_int_length[i][98 - j][2] = 1;						// W (1 - 0)
					tab_NEWS_int_length[i][98 - j][3] = 1;						// S (1 - 0)
					tab_NEWS_int_length[i][98 - j][4] = hypotenuse(&a, &b);	// interface length
				}
			}
		}
	}

	for (int j = 0; j < fill_count; j++)
	{
		tab_NEWS_int_length[angle_count_minus_one][j][0] = 0;					// N
		tab_NEWS_int_length[angle_count_minus_one][j][1] = (j + 1) * 0.01;		// E
		tab_NEWS_int_length[angle_count_minus_one][j][2] = (j + 1) * 0.01;		// W
		tab_NEWS_int_length[angle_count_minus_one][j][3] = 1;					// S
		tab_NEWS_int_length[angle_count_minus_one][j][4] = 1;					// interface length
	}
}
//*********************************
void C_cell::present_yourself()
{
	cout << "----------------------------------\n";
	cout << "faza FCC:\n";
	
	double fcc_mass = 0.0;
	for (int i = 0; i < no_of_elements; i++)
	{
		cout << "mass_ele_al" << i << " = " << element_mass[0][i] << ", ";
		fcc_mass += element_mass[0][i];
	}
	cout << "\nmasa_fcc = " << mass[FCC] << endl;
	cout << "masa_fcc (liczone w present) = " << fcc_mass << endl;

	cout << "\nfaza LIQ:\n";
	double liq_mass = 0.0;
	for (int i = 0; i < no_of_elements; i++)
	{
		cout << "mass_ele_al" << i << " = " << element_mass[1][i] << ", ";
		liq_mass += element_mass[1][i];
	}
	cout << "\nmasa_liq = " << mass[LIQUID] << endl;
	cout << "masa_liq (liczone w present) = " << liq_mass << endl;

	cout << endl;
	cout << "temperatura = " << temperature << endl;
	cout << "masa komorki = " << mass[no_of_phases] << endl;
	cout << "masa komorki (liczone w present) = " << fcc_mass + liq_mass << endl;
	cout << "phase fraction = " << phase_fraction << endl;

	cout << endl;
	cout << "normalna do interfejsu = " << interface_normal << endl;
	cout << "naroznik poczatkowy = " << origin_corner << endl;
	for (int i = 0; i < 4; i++)
	{
		cout << "faza stala na stronie " << i << " = " << *pointer_solid_fraction_on_edge[i] << endl;		
	}
	
	cout << endl;
	for (int j = 0; j < no_of_phases; j++)
	{
		for (int i = 0; i < no_of_elements; i++)
		{
			cout << "potencjal: faza [" << j << "], pierwiastek [" << i << "] = " << element_chemical_potential[j][i] << endl;
		}
	}

	cout << endl;
	for (int j = 0; j < no_of_phases; j++)
	{
		for (int i = 0; i < no_of_elements; i++)
		{
			cout << "stezenie: faza [" << j << "], pierwiastek [" << i << "] = " << element_concentration[j][i] << endl;
		}
	}
	cout << "----------------------------------\n";
}



//void C_cell::get_potentials()
//{
//	//cout << "mass_C = " << mass_C << "\t\t mass_Si = " << mass_Si << endl;
//	//chemical_potential_pointer = system->get_chemical_potentials(temperature, mass_C, mass_Si);
//	//return system->get_chemical_potentials(temperature, mass_C, mass_Si, mass_Fe);
//	chemical_potential_pointer = system->get_chemical_potentials(temperature, mass_C, mass_Si, mass_Fe);
//}
//*********************************
//void C_cell::get_temp_volume()
//{
//	system->get_volume(&temp_volume);
//
//	system->print_equilibrium();
//	cout << "temp vol = " << temp_volume << endl;
//}