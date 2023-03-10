//#include "C_TQ_single_property_old.h"

#include "C_TQ_XX.h"
#include "C_TQ_Fe.h"
#include "C_TQ_single_phase.h"
#include "C_TQ_chem_pot.h"
#include "C_TQ_chem_pot_two_ph.h"
#include "C_TQ_single_property.h"
#include "C_cell.h"
#include <cmath>
#include <fstream>
#include <ctime>

using namespace std;

//--------------------------------------------------------------------
//------------------   PODSTAWOWE DANE O UK?ADZIE   ------------------
const int		no_of_dimensions					= 2,
				no_of_elements						= C_cell::no_of_elements,			// ilo?? pierwiastk?w 1/4
				no_of_phases						= C_cell::no_of_phases;				// konkretne fazy wyst?puj?ce w uk?adzie s? definiowane w C_TQ_Fe.cpp
#define			alloying_elements					"Al", "Ni" /*"Cu"*/ /*"C", "Si"*/ /*"Bi"*/ /*, "Mo",*/ /*"Ni"*/				// ilo?? pierwiastk?w 2/4
char			*element[no_of_elements]			= {alloying_elements, main_element};		// pierwiastek bazowy na ko?cu listy

//----------   DANE ZWI?ZANE Z OBLICZENIAMI: SIATKA, CZAS   ----------

const int		no_of_cells_x						= 100,											// [#] // 2021-03-23 -> 100
				no_of_cells_y						= 1 /*no_of_cells_x*/,							// [#] // 2021-03-23 -> 100
				cell_border_thickness				= 2,											// [#] // 2021-03-23 -> 2
				no_of_cells_x_total					= no_of_cells_x + 2 * cell_border_thickness,	// [#]
				no_of_cells_y_total					= no_of_cells_y + 2 * cell_border_thickness,	// [#]
				cell_loop_x_start					= cell_border_thickness,
				cell_loop_x_stop					= no_of_cells_x + cell_border_thickness,
				cell_loop_y_start					= cell_border_thickness,
				cell_loop_y_stop					= no_of_cells_y + cell_border_thickness,
				immunity_duration					= 100,		// wcze?niej dzia?a?o 200
				calculations_gap					= 100;				
				const double	dx = C_cell::dx,									// [m]
				dx_squared = C_cell::dx * C_cell::dx,						// [m2]
				//dt_multiplier = 1.00000001,										// [dimensionless]
				dt_multiplier = 1.10095,										// [dimensionless]
				initial_dt = 1e-4  /*5.0e-1*/;				//5.0e-6;	/*0.00016285044;*/
				//initial_dt							= 0.00063063085;
double			dt									= initial_dt,	//39.3132600195989e-5,		//		2021-05-18	1e-5 /*1e-7*/,								// [s]	// 2021-05-06 by?o 1e-10; 2021-05-13 by?o 1e-11         2021-05-14 5e-7
				dt_dx_squared						= dt * dx_squared;								// [s * m2]
//const double	dt									= 6e-7 /*1e-7*/,								// [s]	// 2021-05-06 by?o 1e-10; 2021-05-13 by?o 1e-11
const double	gas_constant						= 8.314462618,
				one_over_sqrt2						= 1 / 1.4121,
				half_circle_over_pi					= 180 / 3.14159,
				interface_to_FCC_threshold			= 0.99,
				interface_to_LIQ_threshold			= 1 - interface_to_FCC_threshold,
				interface_trigger_epsilon			= 0.011;

//---------------------   DANE DO OBLICZE? TC   ----------------------

const double	density											= C_cell::density,
				//mobility[no_of_elements]						= {6.01362E-18, 6.01362E-18, 1.87e-19};		// Fe-Al-Ni pierwsza pr?ba
				//mobility[no_of_elements]						= { 6.01362E-18, 1.6E-19, 1.87e-19 };		// Fe-Al-Ni druga pr?ba
				mobility[no_of_elements]						= { 6.01362E-18, 1.6E-19, 1.87e-19 };		// Fe-Al-Ni trzecia pr?ba
				//mobility[no_of_elements]						= { 3.585E-13, 3.585E-13 };						//e-12,e-15,e-12// [kg * mol / (J * s * m)] jako? to trzeba rozbi? -------- to dzia?a 3e-17
				//mobility[no_of_elements]						= { 1.85415E-13, 1.85415E-13};						//e-12,e-15,e-12// [kg * mol / (J * s * m)] jako? to trzeba rozbi? -------- to dzia?a 3e-17
				//mobility[no_of_elements]						= { 1e-10, 1e-12, 1e-16 };						//e-12,e-15,e-12// [kg * mol / (J * s * m)] jako? to trzeba rozbi?		// 2021-05-13:	= {1e-18, 1e-18, 1e-18};		tyle powinno by? = {1e-10, 1e-12, 1e-16};
//tu bylo				//mobility[no_of_elements]						= { 6.01362E-18, 3.00681E-14};											// [s * mol / kg]
				//mobility[no_of_elements]						= { 4.64e-15, 1.94e-17, 6.83e-20 };			// do pierwszych dobrych potr?jnych!								// [s * mol / kg]
				//mobility[no_of_elements]						= {1e-16, 1.5e-16, 2e-16, 2.5e-16, 3e-19};	// ilo?? pierwiastk?w 3/4	//e-12,e-15,e-12// [kg * mol / (J * s * m)] jako? to trzeba rozbi?
double			mobility_coefficient[no_of_phases + 1][no_of_elements],											// phases -> 
				mobility_coefficients_table[no_of_cells_y + 1][no_of_cells_x + 1][no_of_phases + 1][no_of_elements],
				temperature = 1000.0 + 273.15,										// [K]
				
				initial_conc[no_of_phases][no_of_elements - 1] = { 12.44, 18.08, 12.44, 18.08 },		// Beta 17		
				//initial_conc[no_of_phases][no_of_elements - 1] = { 17.35, 17.02, 17.35, 17.02 },		// Beta 17 kawalek dalej		
				initial_conc2[no_of_phases][no_of_elements - 1] = { 31.9, 0.2, 31.9, 0.2 };			// Beta 5
				//initial_conc2[no_of_phases][no_of_elements - 1] = { 14.34, 68.859, 14.34, 68.859 };	// Na podstawie Gamma 2
				//initial_conc[no_of_phases][no_of_elements - 1] = { 18.5, 0.2, 18.5, 0.2 },			// Beta 3
				//initial_conc2[no_of_phases][no_of_elements - 1] = { 30.8, 69.0, 30.8, 69.0 };			// Beta 15
	

				// by?o poni?ej
				//initial_conc[no_of_phases][no_of_elements - 1]	= { 0.8, /*, 0.5, 0.4*/ /*||*/ 0.8/*, 0.7, 0.6*/ };	// ilo?? pierwiastk?w 4/4	// [wt. %, ph0_el0, ph0_el1, ph1_el0, ph1_el1]			
				//initial_conc3[no_of_phases][no_of_elements - 1] = { 29.359, /*69.0,*/ /*0.05, 0.5, 0.4*/ /*||*/ 29.359/*, 69.0*//*, 0.05*/ /*, 0.05, 0.7, 0.6*/ };	// ilo?? pierwiastk?w 4/4	// [wt. %, ph0_el0, ph0_el1, ph1_el0, ph1_el1]			


				//initial_conc[no_of_phases][no_of_elements - 1]	= {0.615, 3.0/*, 0.5, 0.4*/,  /*||*/ 0.615, 3.0/*, 0.7, 0.6*/};	// ilo?? pierwiastk?w 4/4	// [wt. %, ph0_el0, ph0_el1, ph1_el0, ph1_el1]
				//initial_conc[no_of_phases][no_of_elements - 1]	= {3.7, 2.0/*, 0.5, 0.4*/,  /*||*/ 1.0, 0.5/*, 0.7, 0.6*/};	// ilo?? pierwiastk?w 4/4	// gonienie interfejsu [wt. %, ph0_el0, ph0_el1, ph1_el0, ph1_el1]
				//initial_conc[no_of_phases][no_of_elements - 1]	= {1.0, 1.5, /*0.5, */0.05, 0.055/*, 0.7*/};		// [wt. %, ph0_el0, ph0_el1, ph1_el0, ph1_el1]
				//initial_conc[no_of_phases][no_of_elements - 1]	= {0.1, 0.12, /*0.5, */2.0, 2.5/*, 0.7*/};		// [wt. %, ph0_el0, ph0_el1, ph1_el0, ph1_el1]

//-------------   ZMIENNE / DANE GLOBALNE (POMOCNICZE)   -------------

double			glob_var_mass;
const double	mass_accuracy_magnitude = 1e5;
int				no_of_mass_errors = 0;
const int		color_depth = 1024,
				no_of_screen_loop_stepps = 60,				// ilo?? razy ile program daje znac, ?e nadal liczy; ca?kowita ilo?? krok?w (iloczyn trzech loop_stepps)
				no_of_file_loop_stepps = 5,					// (ilo?? plik?w = file * screen)
				no_of_dt_loop_stepps = 5,					// tyle razy zostanie wykonana p?tla calculations zamin zostanie zapisany plik (przy zmiennym dt tyle razy zmieni si? zanim plik zostanie zapisany)
				no_of_calculations_loop_stepps = 50000;	// interwa? pomi?dzy wypisaniami do pliku 100000
				//no_of_screen_loop_stepps = 10,				// ilo?? razy ile program daje znac, ?e nadal liczy; ca?kowita ilo?? krok?w (iloczyn trzech loop_stepps)
				//no_of_file_loop_stepps = 10,					// (ilo?? plik?w = file * screen)
				//no_of_dt_loop_stepps = 1,					// tyle razy zostanie wykonana p?tla calculations zamin zostanie zapisany plik (przy zmiennym dt tyle razy zmieni si? zanim plik zostanie zapisany)
				//no_of_calculations_loop_stepps = 100;	// interwa? pomi?dzy wypisaniami do pliku 100000
string			txt_filename_mass,
				txt_filename_fraction;
fstream			conc_field_txt,
				summary_file;
time_t			start_time,
				end_time;
				
//----------------------------   TEMPY   -----------------------------

int				current_phase;
double			temp_mass_for_txt_output,
				temp_mass_change_state;

//--------------------------   WSKA?NIKI   ---------------------------

double			*chemical_potential_pointer;					// [J/mol]
C_cell			*pointer_Ctab[no_of_cells_y_total][no_of_cells_x_total];

//---------------------------   TABLICE   ----------------------------

double			tab_temperature[no_of_cells_y][no_of_cells_x],																// [K]
				tab_initial_concentration[no_of_cells_y][no_of_cells_x][no_of_phases][no_of_elements - 1],					// [wt. %]
				tab_initial_fraction_solid[no_of_cells_y][no_of_cells_x],													//	# <0; 1>
				tab_initial_mass_from_file[no_of_cells_y][no_of_cells_x],
				tab_delta_mass[no_of_cells_y + 1][no_of_cells_x + 1][no_of_phases][no_of_elements][no_of_dimensions],		// [kg]
				tab_delta_mass_through_interface[no_of_cells_y][no_of_cells_x][no_of_elements],								// [kg]
				tab_min_max_for_txt_file[no_of_elements][2];																// tab[C, Si, mass][min, max]

C_cell::list_of_states tab_initial_state[no_of_cells_y][no_of_cells_x];												// [FCC / LIQ / INTERFACE]

//-------------------------   OBIEKTY KLAS   -------------------------

C_cell (*Ctab)[no_of_cells_x];

//--------------------   OBIEKTY DO OBLICZE? TC   --------------------

C_TQ_chem_pot LIQ_calculations(phase_one, alloying_elements);
C_TQ_chem_pot FCC_calculations(phase_zero, alloying_elements);

//-------------------------   TESTOWANIE   ---------------------------

//--------------------   FUNKCJE DO TESTOWANIA   ---------------------

//---------------------------   FUNKCJE   ----------------------------

double calculate_mobility(C_cell *cell, int *phase, int *element)
{

	////dla Fe-Ni - pierwotne, s?abe wyznaczenie ale dzia?aj?ce 
	//if (*element == 0) // Ni
	//{
	//	double temp_wtpctNi = (cell->element_concentration[*phase][*element] / C_cell::density * 100);
	//	
	//	//double aa = 1.13E-23 * pow(temp_wtpctNi, 3) - 2.2985E-22 * pow(temp_wtpctNi, 2) + 7.7983E-21 * temp_wtpctNi + 2.2323E-20;
	//	//cout << "Ni = " << aa << endl;

	//	//return 2.96973E-21 * temp_wtpctNi + 4.96483E-20;
	//	// ver_2 return 1.13E-23 * pow(temp_wtpctNi, 3) - 2.2985E-22 * pow(temp_wtpctNi, 2) + 7.7983E-21 * temp_wtpctNi + 2.2323E-20;
	//	// ver_1.1 return 1.10761E-23 * pow(temp_wtpctNi, 3) - 2.24323E-22 * pow(temp_wtpctNi, 2) + 5.71988E-21 * temp_wtpctNi + 4.8819E-20;
	//	//return 2.36281E-23 * pow(temp_wtpctNi, 3) - 1.19255E-21 * pow(temp_wtpctNi, 2) + 2.47366E-20 * temp_wtpctNi - 6.46946E-22; // versja dopasowana do wynik?w dictry przy DFe const
	//	//return 4.0274E-24 * pow(temp_wtpctNi, 3) -7.2766E-23 * pow(temp_wtpctNi, 2) -7.4618E-23 * temp_wtpctNi + 8.3370E-20; // versja z danych Smithels przy DFe var v1
	//	//return 1.6269E-24 * pow(temp_wtpctNi, 3) + 1.0427E-22 * pow(temp_wtpctNi, 2) - 1.9248E-21 * temp_wtpctNi + 9.5058E-20; // versja z danych Smithels przy DFe var v2
	//	return -3.76158E-37 * pow(temp_wtpctNi, 3) + 1.05685E-22 * pow(temp_wtpctNi, 2) + 2.50758E-21 * temp_wtpctNi + 8.14994E-20; // versja z danych Smithels przy DFe var v3 - mno?nik 1.2
	//}
	//else // Fe
	//{
	//	double temp_wt_perc = 100 - (cell->element_concentration[*phase][*element] / C_cell::density * 100);
	//	//return -6.06424E-21 * (cell->element_concentration[*phase][*element] / C_cell::density * 100) + 6.93155E-19;
	//	//double temp_wtpctNi = 100 - (cell->element_concentration[*phase][*element] / C_cell::density * 100);

	//	//double bb = 1.4756E-20 * temp_wtpctNi + 7.7342E-20;
	//	//cout << "Fe = " << bb << endl;

	//	// ver_2 return 1.4756E-20 * temp_wtpctNi + 7.7342E-20;
	//	//return 8.65360743099875E-20;
	//	//return 1.6008E-20 * temp_wt_perc + 8.6731E-20; // DFe var v1
	//	//return 5.4770E-22 * pow(temp_wt_perc, 2) - 3.5374E-21 * temp_wt_perc + 8.7723E-20; // DFe var v2
	//	return -3.76158E-37 * pow(temp_wt_perc, 3) + 1.05685E-22 * pow(temp_wt_perc, 2) + 2.50758E-21 * temp_wt_perc + 8.14994E-20; // DFe var v3 - mno?nik 1.2
	//}

	// --------------------------------------------------------

	//// dla Fe-Si - Arrhenius
	//if (*element == 0) // Si
	//{
	//	return (0.000007 * exp(-243000 / (gas_constant * cell->temperature) ) ) / (gas_constant * cell->temperature);
	//}
	//else // Fe
	//{
	//	return (0.0004085 * exp(-311100 / (gas_constant * cell->temperature))) / (gas_constant * cell->temperature);
	//}

	// --------------------------------------------------------

	// dla Fe-Mn - JPN
	if (*element == 0) // Mn
	{
		double temp_wt_perc = (cell->element_concentration[*phase][*element] / C_cell::density * 100);
		if (temp_wt_perc < 6.933)
		{
			return (1.4892e-13 * pow(temp_wt_perc, 4) - 1.95095e-12 * pow(temp_wt_perc, 3) + 5.68074e-12 * pow(temp_wt_perc, 2) + 6.07599e-12 * temp_wt_perc + 3.79908e-11) / 10000 / (gas_constant * cell->temperature);
		}
		else
		{
			return (7.06717e-16 * pow(temp_wt_perc, 4) - 4.73619e-14 * pow(temp_wt_perc, 3) + 1.28557e-12 * pow(temp_wt_perc, 2) - 1.39268e-11 * temp_wt_perc + 9.59941e-11) / 10000 / (gas_constant * cell->temperature);
		}
		//return (1.3066E-12 * temp_wt_perc + 3.9023E-11) / 10000 / (gas_constant * cell->temperature); // v1
		//return (-6.1619E-13	* pow(temp_wt_perc, 2) + 5.5787E-12	* temp_wt_perc + 3.8023E-11) / 10000 / (gas_constant * cell->temperature); // v2
		//return (1.6822E-13 * pow(temp_wt_perc, 3) - 2.7776E-12 * pow(temp_wt_perc, 2) + 1.2478E-11 * temp_wt_perc + 3.8023E-11) / 10000 / (gas_constant * cell->temperature); // v3
		//return (7.5421E-14 * pow(temp_wt_perc, 4) - 9.0036E-13 * pow(temp_wt_perc, 3) + 1.5771E-12 * pow(temp_wt_perc, 2) + 8.5571E-12 * temp_wt_perc + 3.7931E-11) / 10000 / (gas_constant * cell->temperature); // v4
		//return (1.7576E-13 * pow(temp_wt_perc, 4) - 2.2164E-12 * pow(temp_wt_perc, 3) + 6.1651E-12 * pow(temp_wt_perc, 2) + 6.7667E-12 * temp_wt_perc + 3.7496E-11) / 10000 / (gas_constant * cell->temperature); // v5
		//return (1.41683e-13 * pow(temp_wt_perc, 4) - 1.83976e-12 * pow(temp_wt_perc, 3) + 5.21609e-12 * pow(temp_wt_perc, 2) + 6.30582e-12 * temp_wt_perc + 3.79458E-11) / 10000 / (gas_constant * cell->temperature); // v6
	}
	else // Fe
	{
		double temp_wt_perc = 100 - (cell->element_concentration[*phase][*element] / C_cell::density * 100);
		//cout << "\nwt% = " << temp_wt_perc << endl;
		//cout << "mob = " << (-5.70817e-14 * pow(temp_wt_perc, 3) + 6.37286e-13 * pow(temp_wt_perc, 2) + 5.10276e-13 * temp_wt_perc + 2.5213e-11) << endl;
		//cout << "mob ext2 = " << (-5.70817e-14 * pow(temp_wt_perc, 3) + 6.37286e-13 * pow(temp_wt_perc, 2) + 5.10276e-13 * temp_wt_perc + 2.5213e-11) / 10000 << endl;
		//cout << "mob ext = " << (-5.70817e-14 * pow(temp_wt_perc, 3) + 6.37286e-13 * pow(temp_wt_perc, 2) + 5.10276e-13 * temp_wt_perc + 2.5213e-11) / 10000 / (gas_constant * cell->temperature) << endl;
		//cout << "R * T = " << (gas_constant * cell->temperature) << endl;
		return (-5.70817e-14 * pow(temp_wt_perc, 3) + 6.37286e-13 * pow(temp_wt_perc, 2) + 5.10276e-13 * temp_wt_perc + 2.5213e-11) / 10000 / (gas_constant * cell->temperature);
		//return (2.2321E-12 * temp_wt_perc + 2.5065E-11) / 10000 / (gas_constant * cell->temperature);
		//return (1.41683e-13 * pow(temp_wt_perc, 4) - 1.83976e-12 * pow(temp_wt_perc, 3) + 5.21609e-12 * pow(temp_wt_perc, 2) + 6.30582e-12 * temp_wt_perc + 3.79458E-11) / 10000 / (gas_constant * cell->temperature); // v6
	}

	// --------------------------------------------------------

	////dla Fe-C
	//if (*element == 0) // C
	//{
	//	double temp_wt_perc = (cell->element_concentration[*phase][*element] / C_cell::density * 100);
	//	{
	//		double A_of_T = -2.69593e-5 + 7.41913e-8 * cell->temperature - 6.85119e-11 * pow(cell->temperature, 2) + 2.12558e-14 * pow(cell->temperature, 3);
	//		double B_of_T = 1.2217e-5 - 3.20908e-8 * cell->temperature + 2.78898e-11 * pow(cell->temperature, 2) - 7.97634e-15 * pow(cell->temperature, 3);
	//		//cout << "mobility C = " << (A_of_T + B_of_T * temp_wt_perc) / 10000 / (gas_constant * cell->temperature) << endl;
	//		return (A_of_T + B_of_T * temp_wt_perc) / 10000 / (gas_constant * cell->temperature); // dla zale?no?ci od cC
	//		//return (A_of_T + B_of_T * 0.75) / 10000 / (gas_constant * cell->temperature);	// niezale?ne od cC
	//	}
	//}
	//else // Fe
	//{
	//	//cout << "mobility Fe = " << (4.085 * exp(-311100 / (gas_constant * cell->temperature))) / 10000 / (gas_constant * cell->temperature) << endl;
	//	return (4.085 * exp(-311100 / (gas_constant * cell->temperature) ) ) / 10000 / (gas_constant * cell->temperature);
	//}

	//// dla Fe-C-Si
	//if (*element == 0) // C
	//{
	//	double temp_wt_perc = (cell->element_concentration[*phase][*element] / C_cell::density * 100);
	//	{
	//		//double A_of_T = -2.69593e-5 + 7.41913e-8 * cell->temperature - 6.85119e-11 * pow(cell->temperature, 2) + 2.12558e-14 * pow(cell->temperature, 3);
	//		//double B_of_T = 1.2217e-5 - 3.20908e-8 * cell->temperature + 2.78898e-11 * pow(cell->temperature, 2) - 7.97634e-15 * pow(cell->temperature, 3);
	//		////cout << "mobility C = " << (A_of_T + B_of_T * temp_wt_perc) / 10000 / (gas_constant * cell->temperature) << endl;
	//		//return (A_of_T + B_of_T * temp_wt_perc) / 10000 / (gas_constant * cell->temperature); // dla zale?no?ci od cC
	//		////return (A_of_T + B_of_T * 0.75) / 10000 / (gas_constant * cell->temperature);	// niezale?ne od cC

	//		//return 4.99293E-15; // v1
	//		//return 3.10532E-15; // v2
	//		//return 3.72639E-15; // +20% do mniejszego v3 (wcze?niej +25%)
	//		return 3.88165E-15; // +25%


	//	}
	//}
	//else if (*element == 1) // Si
	//{
	//	//return (0.000007 * exp(-243000 / (gas_constant * cell->temperature) ) ) / (gas_constant * cell->temperature);
	//	return 1.62097E-19;
	//}
	//else // Fe	
	//{
	//	//cout << "mobility Fe = " << (4.085 * exp(-311100 / (gas_constant * cell->temperature))) / 10000 / (gas_constant * cell->temperature) << endl;
	//	//return (4.085 * exp(-311100 / (gas_constant * cell->temperature) ) ) / 10000 / (gas_constant * cell->temperature);
	//	return 1.93731E-20;
	//}
}
//*********************************
void calculate_mobility_coefficients_table()
{
	for (int m = 0; m < no_of_cells_y; m++)
	{
		for (int n = 0; n < no_of_cells_x; n++)
		{
			for (int j = 0; j < no_of_phases + 1; j++)
			{
				for (int i = 0; i < no_of_elements; i++)
				{
					mobility_coefficients_table[m][n][j][i] = - calculate_mobility(&Ctab[m][n], &j, &i) * dx * dt /** pow(2, j) *// 2;
				}
			}
		}
	}

	for (int m = 0; m < no_of_cells_y; m++)
	{
		for (int j = 0; j < no_of_phases + 1; j++)
		{
			for (int i = 0; i < no_of_elements; i++)
			{
				mobility_coefficients_table[m][no_of_cells_x][j][i] = mobility_coefficients_table[m][no_of_cells_x - 1][j][i];
			}
		}
	}

	for (int n = 0; n < no_of_cells_x; n++)
	{
		for (int j = 0; j < no_of_phases + 1; j++)
		{
			for (int i = 0; i < no_of_elements; i++)
			{
				mobility_coefficients_table[no_of_cells_y][n][j][i] = mobility_coefficients_table[no_of_cells_y - 1][n][j][i];
			}
		}
	}
}
//*********************************
void calculate_mobility_coefficient()
{
	for (int j = 0; j < no_of_phases + 1; j++)
	{
		for (int i = 0; i < no_of_elements; i++)
		{
			//mobility_coefficient[i] = -mobility[i] * dx * dt / 2;
			mobility_coefficient[j][i] = -mobility[i] * dx * dt /** pow(2, j)*/ / 2;
		}
	}
}
//*********************************
void initial_temperature()
{
	for (int j = 0; j < no_of_cells_y; j++)
	{
		for (int i = 0; i < (no_of_cells_x); i++)
		{
			//tab_temperature[j][i] = temperature + double(i) * 45.0;
			tab_temperature[j][i] = temperature/* + 2 * i*/;
		}
	}
}
//*********************************
void initial_fraction_solid()
{
	string	local_txt_filename = "fraction_solid_input.txt";
	fstream txt_file_fraction_solid_input;
	txt_file_fraction_solid_input.open(local_txt_filename, ios::in);

	for (int j = 0; j < no_of_cells_y; j++)
	{
		for (int i = 0; i < no_of_cells_x; i++)
		{
			txt_file_fraction_solid_input >> tab_initial_fraction_solid[j][i];
		}
	}
	txt_file_fraction_solid_input.close();
}
//*********************************
void initial_state()
{
	for (int j = 0; j < no_of_cells_y; j++)
	{
		for (int i = 0; i < (no_of_cells_x); i++)
		{
			if (tab_initial_fraction_solid[j][i] == 0)
			{
				tab_initial_state[j][i] = C_cell::list_of_states::LIQUID;
			}
			else if (tab_initial_fraction_solid[j][i] == 1)
			{
				tab_initial_state[j][i] = C_cell::list_of_states::FCC;
			}
			else tab_initial_state[j][i] = C_cell::list_of_states::INTERFACE;
		}
	}
}
//*********************************
void initial_concentration()
{
	for (int k = 0; k < no_of_cells_y; k++)
	{
		//for (int j = 0; j < (no_of_cells_x); j++)	// dla 1 materia?u
		for (int j = 0; j < 50; j++)	// dla 2 materia??w
		//for (int j = 0; j < 58; j++)		// dla 3 materia??w
		{
			for (int i = 0; i < no_of_elements - 1; i++)
			{
				if (tab_initial_state[k][j] == C_cell::list_of_states::FCC)
				{
					tab_initial_concentration[k][j][C_cell::list_of_states::FCC][i]		= initial_conc[C_cell::list_of_states::FCC][i]/* + j * 0.40*/;
					tab_initial_concentration[k][j][C_cell::list_of_states::LIQUID][i]	= 0;
				}
				else if (tab_initial_state[k][j] == C_cell::list_of_states::LIQUID)
				{
					tab_initial_concentration[k][j][C_cell::list_of_states::FCC][i]		= 0;
					tab_initial_concentration[k][j][C_cell::list_of_states::LIQUID][i]	= initial_conc[C_cell::list_of_states::LIQUID][i];
				}
				else // INTERFACE
				{
					tab_initial_concentration[k][j][C_cell::list_of_states::FCC][i]		= initial_conc[C_cell::list_of_states::FCC][i];
					tab_initial_concentration[k][j][C_cell::list_of_states::LIQUID][i]	= initial_conc[C_cell::list_of_states::LIQUID][i];
				}
			}
		}

		// 2 r??ne materia?y

		for (int j = 50; j < (no_of_cells_x); j++)	// dla 2 materia??w
		//for (int j = 58; j < 62; j++)	// dla 3 materia??w
		{
			for (int i = 0; i < no_of_elements - 1; i++)
			{
				if (tab_initial_state[k][j] == C_cell::list_of_states::FCC)
				{
					tab_initial_concentration[k][j][C_cell::list_of_states::FCC][i] = initial_conc2[C_cell::list_of_states::FCC][i];
					tab_initial_concentration[k][j][C_cell::list_of_states::LIQUID][i] = 0;
				}
				else if (tab_initial_state[k][j] == C_cell::list_of_states::LIQUID)
				{
					tab_initial_concentration[k][j][C_cell::list_of_states::FCC][i] = 0;
					tab_initial_concentration[k][j][C_cell::list_of_states::LIQUID][i] = initial_conc2[C_cell::list_of_states::LIQUID][i];
				}
				else // INTERFACE
				{
					tab_initial_concentration[k][j][C_cell::list_of_states::FCC][i] = initial_conc2[C_cell::list_of_states::FCC][i];
					tab_initial_concentration[k][j][C_cell::list_of_states::LIQUID][i] = initial_conc2[C_cell::list_of_states::LIQUID][i];
				}
			}
		}

		//// 3 r??ne materia?y

		//for (int j = 62; j < (no_of_cells_x); j++)
		//{
		//	for (int i = 0; i < no_of_elements - 1; i++)
		//	{
		//		if (tab_initial_state[k][j] == C_cell::list_of_states::FCC)
		//		{
		//			tab_initial_concentration[k][j][C_cell::list_of_states::FCC][i] = initial_conc3[C_cell::list_of_states::FCC][i];
		//			tab_initial_concentration[k][j][C_cell::list_of_states::LIQUID][i] = 0;
		//		}
		//		else if (tab_initial_state[k][j] == C_cell::list_of_states::LIQUID)
		//		{
		//			tab_initial_concentration[k][j][C_cell::list_of_states::FCC][i] = 0;
		//			tab_initial_concentration[k][j][C_cell::list_of_states::LIQUID][i] = initial_conc3[C_cell::list_of_states::LIQUID][i];
		//		}
		//		else // INTERFACE
		//		{
		//			tab_initial_concentration[k][j][C_cell::list_of_states::FCC][i] = initial_conc3[C_cell::list_of_states::FCC][i];
		//			tab_initial_concentration[k][j][C_cell::list_of_states::LIQUID][i] = initial_conc3[C_cell::list_of_states::LIQUID][i];
		//		}
		//	}
		//}
		
		
		// ?rednia przy nieparzystej ilo?ci kom?rek
		//for (int i = 0; i < no_of_elements - 1; i++)
		//{
		//	tab_initial_concentration[k][50][C_cell::list_of_states::FCC][i] = (initial_conc[C_cell::list_of_states::FCC][i] + initial_conc2[C_cell::list_of_states::FCC][i]) / 2;
		//	tab_initial_concentration[k][50][C_cell::list_of_states::LIQUID][i] = 0;
		//}

	}
}
//*********************************
void initial_concentration_from_file()
{
	string	local_txt_filename = "concentration_input.ppm";
	fstream txt_file_concentration_input;
	txt_file_concentration_input.open(local_txt_filename, ios::in);

	for (int i = 0; i < no_of_elements - 1; i++)
	{
		for (int k = 0; k < no_of_cells_y; k++)
		{
			for (int j = 0; j < (no_of_cells_x); j++)
			{
				if (tab_initial_state[k][j] == C_cell::list_of_states::FCC)
				{
					txt_file_concentration_input >>	tab_initial_concentration[k][j][C_cell::list_of_states::FCC][i];
						
					//tab_initial_concentration[k][j][C_cell::list_of_states::FCC][i]		= tab_initial_concentration[k][j][C_cell::list_of_states::FCC][i] / (C_cell::volume * 100);
					tab_initial_concentration[k][j][C_cell::list_of_states::LIQUID][i]	= 0;
				}
				else if (tab_initial_state[k][j] == C_cell::list_of_states::LIQUID)
				{
					txt_file_concentration_input >>	tab_initial_concentration[k][j][C_cell::list_of_states::LIQUID][i];

					tab_initial_concentration[k][j][C_cell::list_of_states::FCC][i]		= 0;
					//tab_initial_concentration[k][j][C_cell::list_of_states::LIQUID][i]	= tab_initial_concentration[k][j][C_cell::list_of_states::LIQUID][i] / (C_cell::volume * 100);
				}
				else // INTERFACE
				{
					txt_file_concentration_input >>	tab_initial_concentration[k][j][C_cell::list_of_states::FCC][i];
					txt_file_concentration_input >>	tab_initial_concentration[k][j][C_cell::list_of_states::LIQUID][i];

					tab_initial_concentration[k][j][C_cell::list_of_states::FCC][i] = tab_initial_concentration[k][j][C_cell::list_of_states::FCC][i] / tab_initial_fraction_solid[k][j];
					tab_initial_concentration[k][j][C_cell::list_of_states::LIQUID][i] = tab_initial_concentration[k][j][C_cell::list_of_states::LIQUID][i] / (1 - tab_initial_fraction_solid[k][j]);

					//tab_initial_concentration[k][j][C_cell::list_of_states::LIQUID][i]	= tab_initial_concentration[k][j][C_cell::list_of_states::FCC][i] * (1 - tab_initial_fraction_solid[k][j]);
					//tab_initial_concentration[k][j][C_cell::list_of_states::FCC][i]		= tab_initial_concentration[k][j][C_cell::list_of_states::FCC][i] * tab_initial_fraction_solid[k][j];
				}
			}
		}
	}

	//*********************************

	double local_temp = 0;

	for (int k = 0; k < no_of_cells_y; k++)
	{
		for (int j = 0; j < (no_of_cells_x); j++)
		{
			txt_file_concentration_input >> tab_initial_mass_from_file[k][j];
			local_temp = tab_initial_mass_from_file[k][j] / C_cell::volume / 100;

			for (int l = 0; l < no_of_phases; l++)
			{
				for (int i = 0; i < no_of_elements - 1; i++)
				{
					tab_initial_concentration[k][j][l][i] = tab_initial_concentration[k][j][l][i] / local_temp;
				}
			}
		}
	}

	txt_file_concentration_input.close();
}
//*********************************
void C_cell_initialization()
{
	//utworzenie tablicy obiekt?w c_cell
	Ctab = new C_cell [no_of_cells_y][no_of_cells_x];

	for (int j = 0; j < no_of_cells_y; j++)
	{
		for (int i = 0; i < no_of_cells_x; i++)
		{
			Ctab[j][i].set_initial_conditions(&tab_initial_concentration[j][i][0][0], &tab_initial_state[j][i], &tab_initial_fraction_solid[j][i], &tab_temperature[j][i]);
		}
	}
}
//*********************************
void C_cell_initialization_from_file()
{
	//utworzenie tablicy obiekt?w c_cell
	Ctab = new C_cell [no_of_cells_y][no_of_cells_x];

	for (int j = 0; j < no_of_cells_y; j++)
	{
		for (int i = 0; i < no_of_cells_x; i++)
		{
			Ctab[j][i].set_initial_conditions_mass_from_file(&tab_initial_concentration[j][i][0][0], &tab_initial_state[j][i], &tab_initial_fraction_solid[j][i], &tab_temperature[j][i], &tab_initial_mass_from_file[j][i]);
		}
	}
}
//*********************************
void set_C_cell_pointers()
{
	for (int j = 0; j < no_of_cells_y; j++)
	{
		for (int i = 0; i < no_of_cells_x; i++)
		{
			pointer_Ctab[j + cell_border_thickness][i + cell_border_thickness] = &Ctab[j][i];
		}
	}

	//-----------------------------

	for (int j = 0; j < cell_border_thickness; j++)
	{
		for (int i = 0; i < no_of_cells_x; i++)
		{
			// warunki lustrzane
			pointer_Ctab[1 - j][i + cell_border_thickness]										= &Ctab[j][i];							// S
			pointer_Ctab[j + cell_border_thickness + no_of_cells_y][i + cell_border_thickness]	= &Ctab[no_of_cells_y - 1 - j][i];		// N
		}
	}

	//-----------------------------

	for (int j = 0; j < no_of_cells_y; j++)
	{
		for (int i = 0; i < cell_border_thickness; i++)
		{
			// warunki lustrzane
			pointer_Ctab[j + cell_border_thickness][1 - i]										= &Ctab[j][i];							// W
			pointer_Ctab[j + cell_border_thickness][i + cell_border_thickness + no_of_cells_x]	= &Ctab[j][no_of_cells_x - 1 - i];		// E
		}
	}

	//-----------------------------

	for (int j = 0; j < cell_border_thickness; j++)
	{
		for (int i = 0; i < cell_border_thickness; i++)
		{
			pointer_Ctab[1 - j][1 - i] = &Ctab[j][i];																													// SW
			pointer_Ctab[1 - j][cell_border_thickness + no_of_cells_x + i] = &Ctab[j][no_of_cells_x - 1 - i];															// SE
			pointer_Ctab[cell_border_thickness + no_of_cells_y + j][1 - i] = &Ctab[no_of_cells_y - 1 - j][i];															// NW
			pointer_Ctab[cell_border_thickness + no_of_cells_y + j][cell_border_thickness + no_of_cells_x + i] = &Ctab[no_of_cells_y - 1 - j][no_of_cells_x - 1 - i];	// NE

		}
	}
}
//*********************************
void set_calculatoins_gap();
void calculate_phase_fraction_and_masses();
void calculate_chemical_potentials_and_concentrations_in_entire_cell_table();
//*********************************
void initialization()
{
	//time(&start_time);
	initial_temperature();
	initial_fraction_solid();
	initial_state();
	initial_concentration();
	C_cell_initialization();
	//initial_concentration_from_file();
	//C_cell_initialization_from_file();

	set_calculatoins_gap();															// dodane ?eby liczy? D,B (c_Fe, c_Ni)
	calculate_phase_fraction_and_masses();											// dodane ?eby liczy? D,B (c_Fe, c_Ni)
	calculate_chemical_potentials_and_concentrations_in_entire_cell_table();		// dodane ?eby liczy? D,B (c_Fe, c_Ni)
	
	//dla var D
	//calculate_mobility_coefficients_table();										// dodane ?eby liczy? D,B (c_Fe, c_Ni)

	// dla D const
	calculate_mobility_coefficient();

	set_C_cell_pointers();
}
//*********************************
void reinitiate_C_TQ_LIQ_object()
{
	LIQ_calculations.reinitiate();
	LIQ_calculations.isolate_stable_phase("LIQ");
}
//*********************************
void reinitiate_C_TQ_FCC_object()
{
	FCC_calculations.reinitiate();
	FCC_calculations.isolate_stable_phase(phase_zero);
}
//*********************************
void reinitiate_C_TQ_objects()
{
	reinitiate_C_TQ_LIQ_object();
	reinitiate_C_TQ_FCC_object();
}
//*********************************
void calculate_phase_fraction_and_masses()
{
	for (int j = 0; j < no_of_cells_y; j++)
	{
		for (int i = 0; i < no_of_cells_x; i++)
		{
			Ctab[j][i].calculate_phase_fraction_and_masses();
		}
	}
}
//*********************************
void calculate_interface_normal()
{
	for (int j = cell_loop_y_start; j < cell_loop_y_stop; j++)
	{
		for (int i = cell_loop_x_start; i < cell_loop_x_stop; i++)
		{
			if (pointer_Ctab[j][i]->cell_state == C_cell::list_of_states::FCC)
			{
				pointer_Ctab[j][i]->interface_numerator = 1;
				pointer_Ctab[j][i]->interface_denominator = 1;

				for (int k = 0; k < 4; k++)
				{
					pointer_Ctab[j][i]->pointer_solid_fraction_on_edge[k] = &C_cell::liq_fcc_solid_fraction_on_edge[1];
				}
			}
			else if (pointer_Ctab[j][i]->cell_state == C_cell::list_of_states::LIQUID)
			{
				pointer_Ctab[j][i]->interface_numerator = 1;
				pointer_Ctab[j][i]->interface_denominator = 1;

				for (int k = 0; k < 4; k++)
				{
					pointer_Ctab[j][i]->pointer_solid_fraction_on_edge[k] = &C_cell::liq_fcc_solid_fraction_on_edge[0];
				}			
			}
			else	// INTERFACE
			{
				pointer_Ctab[j][i]->interface_numerator = ( pointer_Ctab[j + 1][i]->phase_fraction
							+ (pointer_Ctab[j + 1][i - 1]->phase_fraction + pointer_Ctab[j + 1][i + 1]->phase_fraction) * one_over_sqrt2
							+ pointer_Ctab[j + 2][i]->phase_fraction * 0.5)
							- ( pointer_Ctab[j - 1][i]->phase_fraction
							+ (pointer_Ctab[j - 1][i - 1]->phase_fraction + pointer_Ctab[j - 1][i + 1]->phase_fraction) * one_over_sqrt2
							+ pointer_Ctab[j - 2][i]->phase_fraction * 0.5);

				if (pointer_Ctab[j][i]->interface_numerator == 0) pointer_Ctab[j][i]->interface_numerator = 1E-10;

				pointer_Ctab[j][i]->interface_denominator = ( pointer_Ctab[j][i + 1]->phase_fraction
							+ (pointer_Ctab[j + 1][i + 1]->phase_fraction + pointer_Ctab[j - 1][i + 1]->phase_fraction) * one_over_sqrt2
							+ pointer_Ctab[j][i + 2]->phase_fraction * 0.5)
							- ( pointer_Ctab[j][i - 1]->phase_fraction
							+ (pointer_Ctab[j + 1][i - 1]->phase_fraction + pointer_Ctab[j - 1][i - 1]->phase_fraction) * one_over_sqrt2
							+ pointer_Ctab[j][i - 2]->phase_fraction * 0.5);

				if (pointer_Ctab[j][i]->interface_denominator == 0) pointer_Ctab[j][i]->interface_denominator = 1E-10;

				pointer_Ctab[j][i]->interface_normal = atan(pointer_Ctab[j][i]->interface_numerator / pointer_Ctab[j][i]->interface_denominator) * half_circle_over_pi;
				pointer_Ctab[j][i]->interface_normal = round(pointer_Ctab[j][i]->interface_normal / 5) * 5;

				//--------------------------------------------------------------------

				if (pointer_Ctab[j][i]->interface_normal < 0)
				{
					pointer_Ctab[j][i]->interface_normal = - pointer_Ctab[j][i]->interface_normal;
				}

				//-----------------------------

				if (pointer_Ctab[j][i]->interface_normal == 0)
				{
					if (pointer_Ctab[j][i]->interface_denominator > 0)							//	 _______
					{																			//	|	  |	|
						pointer_Ctab[j][i]->origin_corner = C_cell::list_of_corners::SE;		//	|	  |	|
						pointer_Ctab[j][i]->interface_normal = 90;								//	|_____|_|o
					}
					else	// denominator < 0														o_______
					{																			//	| |		|
						pointer_Ctab[j][i]->origin_corner = C_cell::list_of_corners::NW;		//	| |		|
						pointer_Ctab[j][i]->interface_normal = 90;								//	|_|_____|
					}
				}
				else if (pointer_Ctab[j][i]->interface_normal == 90)
				{
					if (pointer_Ctab[j][i]->interface_numerator > 0)							//	 _______o
					{																			//	|_______|
						pointer_Ctab[j][i]->origin_corner = C_cell::list_of_corners::NE;		//	|	  	|
						pointer_Ctab[j][i]->interface_normal = 90;								//	|_______|
					}
					else	// numeratorr < 0														 _______
					{																			//	| 		|
						pointer_Ctab[j][i]->origin_corner = C_cell::list_of_corners::SW;		//	|_______|
						pointer_Ctab[j][i]->interface_normal = 90;								// o|_______|
					}
				}
				else // <> 0 oraz <> 90 czyli <5;85>
				{
					if ((pointer_Ctab[j][i])->interface_numerator > 0)
					{
						if ((pointer_Ctab[j][i])->interface_denominator > 0)
						{
							// + / +
							pointer_Ctab[j][i]->origin_corner = C_cell::list_of_corners::NE;
						}
						else
						{
							// + / -
							pointer_Ctab[j][i]->origin_corner = C_cell::list_of_corners::NW;
						}
					}
					else if ((pointer_Ctab[j][i])->interface_denominator > 0)
						{
							// - / +
							pointer_Ctab[j][i]->origin_corner = C_cell::list_of_corners::SE;
						}
						else
						{
							// - / -
							pointer_Ctab[j][i]->origin_corner = C_cell::list_of_corners::SW;
						}
				}

				pointer_Ctab[j][i]->set_NEWS_int_lenght_pointer();
			}
		}
	}
}
//*********************************
void LIQ_to_INT_mass_change_in_rising_interface(C_cell *triggered_cell, const double *local_epsilon)
{
	for (int element = 0; element < no_of_elements; element++)
	{
		temp_mass_change_state = triggered_cell->element_mass[C_cell::list_of_states::LIQUID][element] * *local_epsilon;
		triggered_cell->element_mass[C_cell::list_of_states::FCC][element]	+= temp_mass_change_state;
		triggered_cell->element_mass[C_cell::list_of_states::LIQUID][element]	-= temp_mass_change_state;
		triggered_cell->immunity_countdown = immunity_duration;
		triggered_cell->calculations_countdown = 0;
	}
}
//*********************************
void FCC_to_INT_mass_change_in_rising_interface(C_cell *triggered_cell, const double *local_epsilon)
{
	for (int element = 0; element < no_of_elements; element++)
	{
		temp_mass_change_state = triggered_cell->element_mass[C_cell::list_of_states::FCC][element] * *local_epsilon;
		triggered_cell->element_mass[C_cell::list_of_states::LIQUID][element]	+= temp_mass_change_state;
		triggered_cell->element_mass[C_cell::list_of_states::FCC][element]	-= temp_mass_change_state;
		triggered_cell->immunity_countdown = immunity_duration;
		triggered_cell->calculations_countdown = 0;
	}
}
//*********************************
void iterface_rise_in_triggered_cell(C_cell *triggered_cell, C_cell *triggering_cell)
{
	if (triggered_cell->cell_state == C_cell::list_of_states::LIQUID)
	{
		if (triggering_cell->phase_fraction > interface_trigger_epsilon)
		{
			LIQ_to_INT_mass_change_in_rising_interface(triggered_cell, &interface_trigger_epsilon);
		}
		else
		{
			LIQ_to_INT_mass_change_in_rising_interface(triggered_cell, &triggering_cell->phase_fraction);
		}
	}
	else // used to be FCC
	{
		double LIQ_fraction = 1 - triggering_cell->phase_fraction;
		if (LIQ_fraction > interface_trigger_epsilon)
		{
			FCC_to_INT_mass_change_in_rising_interface(triggered_cell, &interface_trigger_epsilon);
		}
		else
		{
			FCC_to_INT_mass_change_in_rising_interface(triggered_cell, &LIQ_fraction);
		}
	}

	triggered_cell->cell_state = C_cell::list_of_states::INTERFACE;
	triggered_cell->calculate_phase_fraction_and_masses();
	triggered_cell->origin_corner		= triggering_cell->origin_corner;
	triggered_cell->interface_normal	= triggering_cell->interface_normal;
	triggered_cell->set_NEWS_int_lenght_pointer();
	triggered_cell->interface_distance = 0;
	triggered_cell->calculations_gap = &C_cell::tab_calculations_gaps[triggered_cell->interface_distance];
	triggered_cell->calculations_countdown = 0;
}
//*********************************
void check_triggered_cell_trigger_if_needed(C_cell *triggered_cell, C_cell *triggering_cell/*, int *edge_in_triggering_cell*/)
{
	if (triggered_cell->cell_state != C_cell::list_of_states::INTERFACE)
	{
		if (triggered_cell->immunity_countdown == 0)
		{
			iterface_rise_in_triggered_cell(triggered_cell, triggering_cell);		
		}
	}
}
//*********************************
void trigger_interface_movement()
{
	for (int y = cell_loop_y_start; y < cell_loop_y_stop; y++)
	{
		for (int x = cell_loop_x_start; x < cell_loop_x_stop; x++)
		{
			if (pointer_Ctab[y][x]->cell_state == C_cell::list_of_states::INTERFACE)
			{
				if (pointer_Ctab[y][x]->immunity_countdown == 0)
				{
					for (int i = 0; i < 4; i++)
					{
						if ((*pointer_Ctab[y][x]->pointer_solid_fraction_on_edge[i] > 0) && (*pointer_Ctab[y][x]->pointer_solid_fraction_on_edge[i] < 1))
						{
							switch (i)
							{
							case 0: // N triggers into S
								check_triggered_cell_trigger_if_needed(pointer_Ctab[y + 1][x], pointer_Ctab[y][x]);
								break;
							case 1: // E triggers into W
								check_triggered_cell_trigger_if_needed(pointer_Ctab[y][x + 1], pointer_Ctab[y][x]);
								break;
							case 2: // W triggers into E
								check_triggered_cell_trigger_if_needed(pointer_Ctab[y][x - 1], pointer_Ctab[y][x]);
								break;
							default: // (case 3:) S triggers into N
								check_triggered_cell_trigger_if_needed(pointer_Ctab[y - 1][x], pointer_Ctab[y][x]);
							}						
						}
					}
				}
			}
		}
	}
}
//*********************************
void FCC_potentials_in_C_cell_obj(C_cell *cell)
{
	reinitiate_C_TQ_FCC_object();

	chemical_potential_pointer = FCC_calculations.get_chemical_potentials(&cell->temperature, &cell->element_mass[C_cell::list_of_states::FCC][0]);
	FCC_calculations.get_mass(&glob_var_mass);
	glob_var_mass = round(mass_accuracy_magnitude * glob_var_mass / cell->mass[C_cell::list_of_states::FCC]);

	int local_counter = 0;

	while (mass_accuracy_magnitude - glob_var_mass)
	{
		local_counter++;
		if (local_counter > 2)
		{
			cout << local_counter << endl;
			FCC_calculations.calculations_with_print(&cell->temperature, &cell->element_mass[C_cell::list_of_states::FCC][0]);
		}

		cell->temperature += 1e-10;
		chemical_potential_pointer = FCC_calculations.get_chemical_potentials(&cell->temperature, &cell->element_mass[C_cell::list_of_states::FCC][0]);
		FCC_calculations.get_mass(&glob_var_mass);
		glob_var_mass = round(mass_accuracy_magnitude * glob_var_mass / cell->mass[C_cell::list_of_states::FCC]);

		no_of_mass_errors++;
	}

	for (int i = 0; i < no_of_elements; i++)
	{
		cell->element_chemical_potential[C_cell::list_of_states::FCC][i] = chemical_potential_pointer[i];
	}

	cell->calculations_countdown = *cell->calculations_gap;
}
//*********************************
void LIQ_potentials_in_C_cell_obj(C_cell *cell)
{
	reinitiate_C_TQ_LIQ_object();

	chemical_potential_pointer = LIQ_calculations.get_chemical_potentials(&cell->temperature, &cell->element_mass[C_cell::list_of_states::LIQUID][0]);
	LIQ_calculations.get_mass(&glob_var_mass);
	glob_var_mass = round(mass_accuracy_magnitude * abs(glob_var_mass / cell->mass[C_cell::list_of_states::LIQUID]) );
	
	int local_counter = 0;

	while (mass_accuracy_magnitude - glob_var_mass)
	{
		local_counter++;
		if (local_counter > 2)
		{
			cout << local_counter << endl;
			LIQ_calculations.calculations_with_print(&cell->temperature, &cell->element_mass[C_cell::list_of_states::LIQUID][0]);
		}

		cell->temperature += 1e-10;
		chemical_potential_pointer = LIQ_calculations.get_chemical_potentials(&cell->temperature, &cell->element_mass[C_cell::list_of_states::LIQUID][0]);
		LIQ_calculations.get_mass(&glob_var_mass);
		glob_var_mass = round(mass_accuracy_magnitude * glob_var_mass / cell->mass[C_cell::list_of_states::LIQUID]);
		
		no_of_mass_errors++;
	}

	for (int i = 0; i < no_of_elements; i++)
	{
		cell->element_chemical_potential[C_cell::list_of_states::LIQUID][i] = chemical_potential_pointer[i];
	}

	cell->calculations_countdown = *cell->calculations_gap;
}
//*********************************
void calculate_chemical_potentials_in_single_cell(C_cell* cell)
{
	if (cell->cell_state == C_cell::list_of_states::FCC)
	{
		FCC_potentials_in_C_cell_obj(cell);
	}
	else if (cell->cell_state == C_cell::list_of_states::LIQUID)
	{
		LIQ_potentials_in_C_cell_obj(cell);
	}
	else // INTERFACE
	{
		FCC_potentials_in_C_cell_obj(cell);
		LIQ_potentials_in_C_cell_obj(cell);
	}
}
//*********************************
void calculate_chemical_potentials_and_concentrations_in_entire_cell_table()
{
	for (int j = 0; j < no_of_cells_y; j++)
	{
		for (int i = 0; i < no_of_cells_x; i++)
		{
			if (Ctab[j][i].calculations_countdown <= 0)
			{
				calculate_chemical_potentials_in_single_cell(&Ctab[j][i]);
			}
			Ctab[j][i].calculate_concentration_and_chemical_potential_coefficient();
		}
	}
}
//*********************************
void calculate_delta_mass_with_interface_cell(C_cell *current_cell, C_cell *next_cell, int *y, int *x, int *phase, int *dimension)
{
	if (*phase == C_cell::list_of_states::FCC)
	{
		for (int element = 0; element < no_of_elements; element++)
		{
			if (*dimension == 0) // current E (1), next W (2)
			{
				// delta m_________________________________________________________________________________________________ = (- B * dx * dt / 2)__________ * (current_cell solid fraction on E_________________  + next_cell solid fraction on W___________________) / 2
				tab_delta_mass[*y - cell_border_thickness + 1][*x - cell_border_thickness + 1][*phase][element][*dimension] = mobility_coefficient[*phase][element] * (*(current_cell->pointer_solid_fraction_on_edge[1]) + *(next_cell->pointer_solid_fraction_on_edge[2]) ) / 2

				//	* ( C_0________________________________________________ + C_X+_____________________________________________)
					* (current_cell->element_concentration[*phase][element] + next_cell->element_concentration[*phase][element])

				//  * ( mu_X+________________________________________________ * mu_X+_coefficient________________________ - mu_0_____________________________________________________ * mu_0_coefficient____________________________)
					* (next_cell->element_chemical_potential[*phase][element] * next_cell->chemical_potential_coefficient - current_cell->element_chemical_potential[*phase][element] * current_cell->chemical_potential_coefficient);
			}
			else // dimension == 1; current N (0), next S (3)
			{
				// delta m________________________________________________________________________________________________ = (- B * dx * dt / 2)__________ * (current_cell solid fraction on N_________________  + next_cell solid fraction on S___________________) / 2
				tab_delta_mass[*y - cell_border_thickness + 1][*x - cell_border_thickness + 1][*phase][element][*dimension] = mobility_coefficient[*phase][element] * (*(current_cell->pointer_solid_fraction_on_edge[0]) + *(next_cell->pointer_solid_fraction_on_edge[3]) ) / 2

				//	* ( C_0________________________________________________ + C_X+_____________________________________________)
					* (current_cell->element_concentration[*phase][element] + next_cell->element_concentration[*phase][element])

				//  * ( mu_X+________________________________________________ * mu_X+_coefficient_______________________  - mu_0_____________________________________________________ * mu_0_coefficient____________________________)
					* (next_cell->element_chemical_potential[*phase][element] * next_cell->chemical_potential_coefficient - current_cell->element_chemical_potential[*phase][element] * current_cell->chemical_potential_coefficient);
			}
		}
	}
	else // (*phase == C_cell::list_of_states::LIQUID) <- teoretycznie LIQUID lub INTERFACE ale ta funkcja nie powinna by? nigdy wywolana z faza jako INTERFACE
	{
		for (int element = 0; element < no_of_elements; element++)
		{
			if (*dimension == 0) // current E (1), next W (2)
			{
				// delta m_________________________________________________________________________________________________ = (- B * dx * dt / 2)__________ * ( (1 - current_cell solid fraction on E___________________) + (1 - next_cell solid fraction on W___________________)  ) / 2
				tab_delta_mass[*y - cell_border_thickness + 1][*x - cell_border_thickness + 1][*phase][element][*dimension] = mobility_coefficient[*phase][element] * ( (1 - *(current_cell->pointer_solid_fraction_on_edge[1]) ) + (1 - *(next_cell->pointer_solid_fraction_on_edge[2]) )  ) / 2

				//	* ( C_0________________________________________________ + C_X+_____________________________________________)
					* (current_cell->element_concentration[*phase][element] + next_cell->element_concentration[*phase][element])

				//  * ( mu_X+________________________________________________ * mu_X+_coefficient___________________________  - mu_0_____________________________________________________ * mu_0_coefficient_____________________________)
					* (next_cell->element_chemical_potential[*phase][element] * next_cell-> chemical_potential_coefficient - current_cell->element_chemical_potential[*phase][element] * current_cell-> chemical_potential_coefficient);
			}
			else // dimension == 1; current N (0), next S (3)
			{
				// delta m_________________________________________________________________________________________________ = (- B * dx * dt / 2)__________ * ( (1 - current_cell solid fraction on N___________________) + (1 - next_cell solid fraction on S___________________) ) / 2
				tab_delta_mass[*y - cell_border_thickness + 1][*x - cell_border_thickness + 1][*phase][element][*dimension] = mobility_coefficient[*phase][element] * ( (1 - *(current_cell->pointer_solid_fraction_on_edge[0]) ) + (1 - *(next_cell->pointer_solid_fraction_on_edge[3]) ) ) / 2

				//	* ( C_0________________________________________________ + C_X+_____________________________________________)
					* (current_cell->element_concentration[*phase][element] + next_cell->element_concentration[*phase][element])

				//  * ( mu_X+________________________________________________ * mu_X+_coefficient________________________ - mu_0_____________________________________________________ * mu_0_coefficient____________________________)
					* (next_cell->element_chemical_potential[*phase][element] * next_cell->chemical_potential_coefficient - current_cell->element_chemical_potential[*phase][element] * current_cell->chemical_potential_coefficient);
			}
		}		
	}
}
//*********************************
void calculate_delta_mass(C_cell *current_cell, C_cell *next_cell, int *y, int *x, int *phase, int *dimension)
{
	for (int element = 0; element < no_of_elements; element++)
	{
		double *delta_mass_here = &tab_delta_mass[*y - cell_border_thickness + 1][*x - cell_border_thickness + 1][*phase][element][*dimension];

		//// 2021-04-29 dla zale?no?ci D od T i c
		//// delta m_________________________________________________________________________________________________ = (- B * dx * dt / 2)__________ 
		//tab_delta_mass[*y - cell_border_thickness + 1][*x - cell_border_thickness + 1][*phase][element][*dimension] = mobility_coefficients_table[*y - 1][*x - 1][*phase][element]

		//	//	* ( C_0________________________________________________ + C_X+_____________________________________________)
		//	* (current_cell->element_concentration[*phase][element] + next_cell->element_concentration[*phase][element])

		//	//  * ( mu_X+________________________________________________ * (1 - mu_X+_coefficient________________________) - mu_0_____________________________________________________ * (1 - mu_0_coefficient____________________________) )
		//	* (next_cell->element_chemical_potential[*phase][element] * (1 - next_cell->chemical_potential_coefficient) - current_cell->element_chemical_potential[*phase][element] * (1 - current_cell->chemical_potential_coefficient));


		// 2021-10-20 dla sta?ych D
		// delta m_________________________________________________________________________________________________ = (- B * dx * dt / 2)__________ 
		tab_delta_mass[*y - cell_border_thickness + 1][*x - cell_border_thickness + 1][*phase][element][*dimension] = mobility_coefficient[*phase][element]

			//	* ( C_0________________________________________________ + C_X+_____________________________________________)
			* (current_cell->element_concentration[*phase][element] + next_cell->element_concentration[*phase][element])

		//  * ( mu_X+________________________________________________ * (1 - mu_X+_coefficient________________________) - mu_0_____________________________________________________ * (1 - mu_0_coefficient____________________________) )
			* (next_cell->element_chemical_potential[*phase][element] * (1 - next_cell->chemical_potential_coefficient) - current_cell->element_chemical_potential[*phase][element] * (1 - current_cell->chemical_potential_coefficient) );



		//// 2021-09-13 zmiana strumienia w zale?no?ci czy jest wychodz?cy czy przychodz?cy z kom?rki o masie wi?kszej lub mniejszej ni? nominalna
		//// delta m_________________________________________________________________________________________________ = (- B * dx * dt / 2)__________ 
		//*delta_mass_here = mobility_coefficient[*phase][element]

		////	* ( C_0________________________________________________ + C_X+_____________________________________________)
		//	* (current_cell->element_concentration[*phase][element] + next_cell->element_concentration[*phase][element])

		////  * ( mu_X+________________________________________________ - mu_0_____________________________________________________
		//	* (next_cell->element_chemical_potential[*phase][element] - current_cell->element_chemical_potential[*phase][element]);

		//if (current_cell->chemical_potential_coefficient > next_cell->chemical_potential_coefficient)
		//{
		//	if (*delta_mass_here < 0)
		//	{
		//		*delta_mass_here = 0;
		//	}
		//}








		//// delta m_________________________________________________________________________________________________ = (- B * dx * dt / 2)__________ 
		//tab_delta_mass[*y - cell_border_thickness + 1][*x - cell_border_thickness + 1][*phase][element][*dimension] = mobility_coefficient[*phase][element]

		//	//	* ( C_0________________________________________________ + C_X+_____________________________________________)
		//	* (current_cell->element_concentration[*phase][element] + next_cell->element_concentration[*phase][element])

		//	//  * ( mu_X+________________________________________________ - mu_0_____________________________________________________) * mu_0_coefficient____________________________)
		//	* (next_cell->element_chemical_potential[*phase][element] - current_cell->element_chemical_potential[*phase][element]);// * (next_cell->chemical_potential_coefficient - current_cell->chemical_potential_coefficient);
	
		//if (tab_delta_mass[*y - cell_border_thickness + 1][*x - cell_border_thickness + 1][*phase][element][*dimension] > 0)
		//{
		//	tab_delta_mass[*y - cell_border_thickness + 1][*x - cell_border_thickness + 1][*phase][element][*dimension] *= next_cell->chemical_potential_coefficient / current_cell->chemical_potential_coefficient;
		//}
		//else
		//{
		//	tab_delta_mass[*y - cell_border_thickness + 1][*x - cell_border_thickness + 1][*phase][element][*dimension] *= current_cell->chemical_potential_coefficient / next_cell->chemical_potential_coefficient;
		//}
		//// do dupy to, zrobi? dwie sk?adowe, jedna jak by?a z napr??eniami maj?cymi wp?yw na potencja? chemiczny, druga relaksacja napr??e? - samoistne przechodzenie masy z rejon?w o wi?kszej masie do mniejszej
	}
}
//*********************************
void zero_delta_mass_values(int *y, int *x, int *phase, int *dimension)
{
	for (int element = 0; element < no_of_elements; element++)
	{
		tab_delta_mass[*y - cell_border_thickness + 1][*x - cell_border_thickness + 1][*phase][element][*dimension] = 0;
	}	
}
//*********************************
inline void zero_delta_mass_through_interface(int *y, int *x)
{
	for (int element = 0; element < no_of_elements; element++)
	{
		tab_delta_mass_through_interface[*y - cell_border_thickness][*x - cell_border_thickness][element] = 0;
	}
}
//*********************************
void calculate_delta_mass_through_interface(C_cell *current_cell, int *y, int *x)
{
	for (int element = 0; element < no_of_elements; element++)
	{	// we wzorze s? wsp??rz?dne w tablicy z g?ry przyj?te jako 0 lub 1, to s? poszczeg?lne fazy, od tego jak wygl?da ta funkcja zale?y znak przy przep?ywie pomi?dzy LIQ a FCC
		// teraz: (-) oznacza, ?e TRACI faza LIQ, (+) oznacza, ?e ZYSKUJE faza LIQ
		// delta m_______________________________________________________________________________________ = (- B * dx * dt / 2)____________________________________________________
		tab_delta_mass_through_interface[*y - cell_border_thickness][*x - cell_border_thickness][element] = mobility_coefficient[C_cell::list_of_states::INTERFACE][element] * *current_cell->pointer_interface_length

		//	* ( C_(phase_0 = FCC)_____________________________ + C_(phase_1 = LIQ_______________________________)
			* (current_cell->element_concentration[0][element] + current_cell->element_concentration[1][element])

		//  * ( mu_(phase_1 = LIQ)_________________________________ * mu_coefficient______________________________ - mu_(phase_0 = FCC)__________________________________ * mu_coefficient_______________________________)
			* (current_cell->element_chemical_potential[1][element] * current_cell->chemical_potential_coefficient - current_cell->element_chemical_potential[0][element] * current_cell-> chemical_potential_coefficient);
	}
}
//*********************************
void check_current_cell_and_neighbour_then_calculate_delta_mass(C_cell *current_cell, C_cell *next_cell, int *y, int *x, int dimension)
{
	if (current_cell->cell_state == C_cell::list_of_states::FCC)
	{
		if ( (*y > cell_border_thickness - 1) && (*x > cell_border_thickness - 1) )
		{
			zero_delta_mass_through_interface(y, x);
		}

		if (next_cell->cell_state == C_cell::list_of_states::FCC) // FCC-FCC
		{
			current_phase = C_cell::list_of_states::FCC;
			calculate_delta_mass(current_cell, next_cell, y, x, &current_phase, &dimension);

			current_phase = C_cell::list_of_states::LIQUID;
			zero_delta_mass_values(y, x, &current_phase, &dimension);
		}
		else if (next_cell->cell_state == C_cell::list_of_states::INTERFACE) // FCC-INTERFACE
		{
			current_phase = C_cell::list_of_states::FCC;
			calculate_delta_mass_with_interface_cell(current_cell, next_cell, y, x, &current_phase, &dimension);

			current_phase = C_cell::list_of_states::LIQUID;
			zero_delta_mass_values(y, x, &current_phase, &dimension);				
		}
		else // FCC-LIQ
		{
			current_phase = C_cell::list_of_states::FCC;
			zero_delta_mass_values(y, x, &current_phase, &dimension);

			current_phase = C_cell::list_of_states::LIQUID;
			zero_delta_mass_values(y, x, &current_phase, &dimension);
		}
	}
	else if (current_cell->cell_state == C_cell::list_of_states::LIQUID)
	{
		if ( (*y > cell_border_thickness - 1) && (*x > cell_border_thickness - 1) )
		{
			zero_delta_mass_through_interface(y, x);
		}

		if (next_cell->cell_state == C_cell::list_of_states::LIQUID) // LIQ-LIQ
		{
			current_phase = C_cell::list_of_states::LIQUID;
			calculate_delta_mass(current_cell, next_cell, y, x, &current_phase, &dimension);

			current_phase = C_cell::list_of_states::FCC;
			zero_delta_mass_values(y, x, &current_phase, &dimension);				
		}
		else if (next_cell->cell_state == C_cell::list_of_states::INTERFACE) // LIQ-INTERFACE
		{
			current_phase = C_cell::list_of_states::LIQUID;
			calculate_delta_mass_with_interface_cell(current_cell, next_cell, y, x, &current_phase, &dimension);

			current_phase = C_cell::list_of_states::FCC;
			zero_delta_mass_values(y, x, &current_phase, &dimension);
		}
		else // LIQ-FCC
		{
			current_phase = C_cell::list_of_states::FCC;
			zero_delta_mass_values(y, x, &current_phase, &dimension);

			current_phase = C_cell::list_of_states::LIQUID;
			zero_delta_mass_values(y, x, &current_phase, &dimension);					
		}
	}
	else // INTERFACE
	{
		if ( (*y > cell_border_thickness - 1) && (*x > cell_border_thickness - 1) )
		{
			calculate_delta_mass_through_interface(current_cell, y, x);
		}

		if (next_cell->cell_state == C_cell::list_of_states::FCC) // INTERFACE-FCC
		{
			current_phase = C_cell::list_of_states::FCC;
			calculate_delta_mass_with_interface_cell(current_cell, next_cell, y, x, &current_phase, &dimension);

			current_phase = C_cell::list_of_states::LIQUID;
			zero_delta_mass_values(y, x, &current_phase, &dimension);								
		}
		else if (next_cell->cell_state == C_cell::list_of_states::LIQUID) // INTERFACE-LIQ
		{
			current_phase = C_cell::list_of_states::LIQUID;
			calculate_delta_mass_with_interface_cell(current_cell, next_cell, y, x, &current_phase, &dimension);

			current_phase = C_cell::list_of_states::FCC;
			zero_delta_mass_values(y, x, &current_phase, &dimension);				
		}
		else // INTERFACE-INTERFACE
		{
			current_phase = C_cell::list_of_states::FCC;
			calculate_delta_mass_with_interface_cell(current_cell, next_cell, y, x, &current_phase, &dimension);

			current_phase = C_cell::list_of_states::LIQUID;
			calculate_delta_mass_with_interface_cell(current_cell, next_cell, y, x, &current_phase, &dimension);
		}
	}
}
//*********************************
void check_current_cell_and_neighbour_then_calculate_delta_mass_no_interface_calculations(C_cell *current_cell, C_cell *next_cell, int *y, int *x, int dimension)
{
	// ta funkcja wywo?ywana jako 2. w kolejno?ci (wymiana masy z kom?rk? powy?ej), przep?yw przez interfejs zosta? ju? przeliczony przy pierwszym wywo?aniu z kom?rk? na prawo
	if (current_cell->cell_state == C_cell::list_of_states::FCC)
	{
		if (next_cell->cell_state == C_cell::list_of_states::FCC) // FCC-FCC
		{
			current_phase = C_cell::list_of_states::FCC;
			calculate_delta_mass(current_cell, next_cell, y, x, &current_phase, &dimension);

			current_phase = C_cell::list_of_states::LIQUID;
			zero_delta_mass_values(y, x, &current_phase, &dimension);
		}
		else if (next_cell->cell_state == C_cell::list_of_states::INTERFACE) // FCC-INTERFACE
		{
			current_phase = C_cell::list_of_states::FCC;
			calculate_delta_mass_with_interface_cell(current_cell, next_cell, y, x, &current_phase, &dimension);

			current_phase = C_cell::list_of_states::LIQUID;
			zero_delta_mass_values(y, x, &current_phase, &dimension);				
		}
		else // FCC-LIQ
		{
			current_phase = C_cell::list_of_states::FCC;
			zero_delta_mass_values(y, x, &current_phase, &dimension);

			current_phase = C_cell::list_of_states::LIQUID;
			zero_delta_mass_values(y, x, &current_phase, &dimension);
		}
	}
	else if (current_cell->cell_state == C_cell::list_of_states::LIQUID)
	{
		if (next_cell->cell_state == C_cell::list_of_states::LIQUID) // LIQ-LIQ
		{
			current_phase = C_cell::list_of_states::LIQUID;
			calculate_delta_mass(current_cell, next_cell, y, x, &current_phase, &dimension);

			current_phase = C_cell::list_of_states::FCC;
			zero_delta_mass_values(y, x, &current_phase, &dimension);				
		}
		else if (next_cell->cell_state == C_cell::list_of_states::INTERFACE) // LIQ-INTERFACE
		{
			current_phase = C_cell::list_of_states::LIQUID;
			calculate_delta_mass_with_interface_cell(current_cell, next_cell, y, x, &current_phase, &dimension);

			current_phase = C_cell::list_of_states::FCC;
			zero_delta_mass_values(y, x, &current_phase, &dimension);
		}
		else // LIQ-FCC
		{
			current_phase = C_cell::list_of_states::FCC;
			zero_delta_mass_values(y, x, &current_phase, &dimension);

			current_phase = C_cell::list_of_states::LIQUID;
			zero_delta_mass_values(y, x, &current_phase, &dimension);					
		}
	}
	else // INTERFACE
	{
		if (next_cell->cell_state == C_cell::list_of_states::FCC) // INTERFACE-FCC
		{
			current_phase = C_cell::list_of_states::FCC;
			calculate_delta_mass_with_interface_cell(current_cell, next_cell, y, x, &current_phase, &dimension);

			current_phase = C_cell::list_of_states::LIQUID;
			zero_delta_mass_values(y, x, &current_phase, &dimension);								
		}
		else if (next_cell->cell_state == C_cell::list_of_states::LIQUID) // INTERFACE-LIQ
		{
			current_phase = C_cell::list_of_states::LIQUID;
			calculate_delta_mass_with_interface_cell(current_cell, next_cell, y, x, &current_phase, &dimension);

			current_phase = C_cell::list_of_states::FCC;
			zero_delta_mass_values(y, x, &current_phase, &dimension);				
		}
		else // INTERFACE-INTERFACE
		{
			current_phase = C_cell::list_of_states::FCC;
			calculate_delta_mass_with_interface_cell(current_cell, next_cell, y, x, &current_phase, &dimension);

			current_phase = C_cell::list_of_states::LIQUID;
			calculate_delta_mass_with_interface_cell(current_cell, next_cell, y, x, &current_phase, &dimension);
		}
	}
}
//*********************************
void calculate_delta_mass_in_entire_cell_table()
{
	for (int y = 0 + cell_border_thickness - 1; y < no_of_cells_y + cell_border_thickness; y++)
	{
		for (int x = 0 + cell_border_thickness - 1; x < no_of_cells_x + cell_border_thickness; x++)
		{
			check_current_cell_and_neighbour_then_calculate_delta_mass(pointer_Ctab[y][x], pointer_Ctab[y][x + 1], &y, &x, 0);
			check_current_cell_and_neighbour_then_calculate_delta_mass_no_interface_calculations(pointer_Ctab[y][x], pointer_Ctab[y + 1][x], &y, &x, 1);
		}
	}

	//for (int y = 0; y < no_of_cells_y; y++)
	//{
	//	for (int x = 0; x < no_of_cells_x; x++)
	//	{
	//		if (pointer_Ctab[y + cell_border_thickness][x + cell_border_thickness]->old_chemical_potential_coefficient > 0.0005)
	//		{
	//			for (int element = 0; element < no_of_elements; element++)
	//			{
	//				for (int dimension = 0; dimension < no_of_dimensions; dimension++)
	//				{
	//					if (tab_delta_mass[y + 1][x + 1][0][element][dimension] < 0)
	//						tab_delta_mass[y + 1][x + 1][0][element][dimension] = 0;
	//				}
	//				
	//				if (tab_delta_mass[y][x + 1][0][element][1] > 0)
	//					tab_delta_mass[y][x + 1][0][element][1] = 0;

	//				if (tab_delta_mass[y + 1][x][0][element][0] > 0)
	//					tab_delta_mass[y + 1][x][0][element][0] = 0;
	//			}
	//		}
	//	}
	//}
}
//*********************************
void single_timestep_calculations_no_triggering()
{
	calculate_phase_fraction_and_masses();
	calculate_interface_normal();
	//trigger_interface_movement();
	calculate_chemical_potentials_and_concentrations_in_entire_cell_table();
	//calculate_delta_mass_in_entire_cell_table();
}
//*********************************
void single_timestep_calculations()
{
	calculate_phase_fraction_and_masses();
	calculate_interface_normal();
	trigger_interface_movement();
	calculate_chemical_potentials_and_concentrations_in_entire_cell_table();
	calculate_delta_mass_in_entire_cell_table();
}
//*********************************
void update_mass_after_calculations()
{
	for (int y = cell_loop_y_start; y < cell_loop_y_stop; y++)
	{
		for (int x = cell_loop_x_start; x < cell_loop_x_stop; x++)
		{
			for (int phase = 0; phase < no_of_phases; phase++)
			{
				for (int element = 0; element < no_of_elements; element++)
				{
					for (int direction = 0; direction < no_of_dimensions; direction++)
					{
						pointer_Ctab[y][x]->element_mass[phase][element] -= tab_delta_mass[y - 1][x - 1][phase][element][direction];
					}

					pointer_Ctab[y][x]->element_mass[phase][element] += tab_delta_mass[y - 1][x - 2][phase][element][0];
					pointer_Ctab[y][x]->element_mass[phase][element] += tab_delta_mass[y - 2][x - 1][phase][element][1];
				}
			}

			for (int element = 0; element < no_of_elements; element++)
			{
				pointer_Ctab[y][x]->element_mass[C_cell::list_of_states::LIQUID][element]	+= tab_delta_mass_through_interface[y - 2][x - 2][element];
				pointer_Ctab[y][x]->element_mass[C_cell::list_of_states::FCC][element]		-= tab_delta_mass_through_interface[y - 2][x - 2][element];
			}
		}
	}	
}
//*********************************
void relax_stress()
{
	for (int y = cell_loop_y_start; y < cell_loop_y_stop; y++)
	{
		for (int x = cell_loop_x_start; x < cell_loop_x_stop; x++)
		{
			double temp_mass_difference = (pointer_Ctab[y][x]->mass[no_of_phases] - pointer_Ctab[y][x + 1]->mass[no_of_phases]) / C_cell::nominal_mass / 6;

			//if (x == 51)
			//{
			//	cout << "xx";
			//}

			for (int element = 0; element < no_of_elements; element++)
			{
				if (temp_mass_difference < 0)
				{
					double temp_delta_mass = pointer_Ctab[y][x + 1]->element_mass[0][element] * temp_mass_difference;

					pointer_Ctab[y][x + 1]->element_mass[0][element]	+= temp_delta_mass;
					pointer_Ctab[y][x]->element_mass[0][element]		-= temp_delta_mass;
				}
				else
				{
					double temp_delta_mass = pointer_Ctab[y][x]->element_mass[0][element] * temp_mass_difference;

					pointer_Ctab[y][x]->element_mass[0][element]		-= temp_delta_mass;
					pointer_Ctab[y][x + 1]->element_mass[0][element]	+= temp_delta_mass;
				}
			}

		}
	}
}
//*********************************
void check_interface_vanishing()
{
	for (int y = cell_loop_y_start; y < cell_loop_y_stop; y++)
	{
		for (int x = cell_loop_x_start; x < cell_loop_x_stop; x++)
		{
			if (pointer_Ctab[y][x]->cell_state == C_cell::list_of_states::INTERFACE)
			{
				if (pointer_Ctab[y][x]->immunity_countdown == 0)
				{
					if (pointer_Ctab[y][x]->phase_fraction > interface_to_FCC_threshold)	// INT TO FCC
					{
						//*********************************
						for (int i = 0; i < 4; i++)
						{
							if (*pointer_Ctab[y][x]->pointer_solid_fraction_on_edge[i] == 0)
							{
								switch (i)
								{
								case 0: // N triggers into S
									check_triggered_cell_trigger_if_needed(pointer_Ctab[y + 1][x], pointer_Ctab[y][x]);
									break;
								case 1: // E triggers into W
									check_triggered_cell_trigger_if_needed(pointer_Ctab[y][x + 1], pointer_Ctab[y][x]);
									break;
								case 2: // W triggers into E
									check_triggered_cell_trigger_if_needed(pointer_Ctab[y][x - 1], pointer_Ctab[y][x]);
									break;
								default: // (case 3:) S triggers into N
									check_triggered_cell_trigger_if_needed(pointer_Ctab[y - 1][x], pointer_Ctab[y][x]);
								}						
							}
						}
						//*********************************
						for (int element = 0; element < no_of_elements; element++)
						{
							pointer_Ctab[y][x]->element_mass[C_cell::list_of_states::FCC][element] += pointer_Ctab[y][x]->element_mass[C_cell::list_of_states::LIQUID][element];
							pointer_Ctab[y][x]->element_mass[C_cell::list_of_states::LIQUID][element] = 0;
						}

						pointer_Ctab[y][x]->cell_state = C_cell::list_of_states::FCC;
						pointer_Ctab[y][x]->phase_fraction = 1;
						pointer_Ctab[y][x]->immunity_countdown = immunity_duration;
						pointer_Ctab[y][x]->interface_distance = 1;
						pointer_Ctab[y][x]->calculations_gap = &C_cell::tab_calculations_gaps[pointer_Ctab[y][x]->interface_distance];
						pointer_Ctab[y][x]->calculations_countdown = 0;
					}
					else if (pointer_Ctab[y][x]->phase_fraction < interface_to_LIQ_threshold)	// INT TO LIQ
					{
						//*********************************
						for (int i = 0; i < 4; i++)
						{
							if (*pointer_Ctab[y][x]->pointer_solid_fraction_on_edge[i] == 1)
							{
								switch (i)
								{
								case 0: // N triggers into S
									check_triggered_cell_trigger_if_needed(pointer_Ctab[y + 1][x], pointer_Ctab[y][x]);
									break;
								case 1: // E triggers into W
									check_triggered_cell_trigger_if_needed(pointer_Ctab[y][x + 1], pointer_Ctab[y][x]);
									break;
								case 2: // W triggers into E
									check_triggered_cell_trigger_if_needed(pointer_Ctab[y][x - 1], pointer_Ctab[y][x]);
									break;
								default: // (case 3:) S triggers into N
									check_triggered_cell_trigger_if_needed(pointer_Ctab[y - 1][x], pointer_Ctab[y][x]);
								}						
							}
						}
						//*********************************
						for (int element = 0; element < no_of_elements; element++)
						{
							pointer_Ctab[y][x]->element_mass[C_cell::list_of_states::LIQUID][element] += pointer_Ctab[y][x]->element_mass[C_cell::list_of_states::FCC][element];
							pointer_Ctab[y][x]->element_mass[C_cell::list_of_states::FCC][element] = 0;
						}
						
						pointer_Ctab[y][x]->cell_state = C_cell::list_of_states::LIQUID;
						pointer_Ctab[y][x]->phase_fraction = 0;
						pointer_Ctab[y][x]->immunity_countdown = immunity_duration;
						pointer_Ctab[y][x]->interface_distance = 1;
						pointer_Ctab[y][x]->calculations_gap = &C_cell::tab_calculations_gaps[pointer_Ctab[y][x]->interface_distance];
						pointer_Ctab[y][x]->calculations_countdown = 0;
					}
				}
			}
		}
	}	
}
//*********************************
void set_file_name(int *file_no)
{
	txt_filename_mass = "c:/Users/rz_m/Nextcloud/wydobywanie/wydobywanie/!output/";
	txt_filename_mass += "mass_";
	if (*file_no < 100)
	{
		txt_filename_mass += "0";
		if (*file_no < 10)
		{
			txt_filename_mass += "0";
		}
	}
	txt_filename_mass += to_string(*file_no);
	txt_filename_mass += ".ppm";
	//cout << txt_filename_mass << endl;

	//*********************************

	txt_filename_fraction = "c:/Users/rz_m/Nextcloud/wydobywanie/wydobywanie/!output/";
	txt_filename_fraction += "fraction_";
	if (*file_no < 100)
	{
		txt_filename_fraction += "0";
		if (*file_no < 10)
		{
			txt_filename_fraction += "0";
		}
	}
	txt_filename_fraction += to_string(*file_no);
	txt_filename_fraction += ".ppm";
	//cout << txt_filename_fraction << endl;
}
//*********************************
void export_mass_to_file(double *calculated_total_time)
{
	conc_field_txt.open(txt_filename_mass, ios::out);
	conc_field_txt.precision(8);		// 2021-05-06 by?o 5	//2021-05-18 by?o 6

	for (int k = 0; k < no_of_elements - 1; k++)
	{
		for (int j = 0; j < no_of_cells_y; j++)
		{
			for (int i = 0; i < no_of_cells_x; i++)
			{
				if (Ctab[j][i].cell_state == C_cell::list_of_states::FCC)
				{
					conc_field_txt << Ctab[j][i].element_concentration[C_cell::list_of_states::FCC][k] * Ctab[j][i].phase_fraction << "\t";
				}
				else if (Ctab[j][i].cell_state == C_cell::list_of_states::LIQUID)
				{
					conc_field_txt << Ctab[j][i].element_concentration[C_cell::list_of_states::LIQUID][k] * (1 - Ctab[j][i].phase_fraction) << "\t";
				}
				else //INTERFACE
				{
					conc_field_txt << Ctab[j][i].element_concentration[C_cell::list_of_states::FCC][k] * Ctab[j][i].phase_fraction << "\t";
					conc_field_txt << Ctab[j][i].element_concentration[C_cell::list_of_states::LIQUID][k] * (1 - Ctab[j][i].phase_fraction) << "\t";
				}
				//conc_field_txt << Ctab[j][i].element_concentration[0][k] * Ctab[j][i].phase_fraction + Ctab[j][i].element_concentration[1][k] * (1 - Ctab[j][i].phase_fraction) << "\t";
			}
			conc_field_txt << "\n";
		}
	}

	for (int j = 0; j < no_of_cells_y; j++)
	{
		for (int i = 0; i < no_of_cells_x; i++)
		{
			conc_field_txt << Ctab[j][i].mass[no_of_phases] << "\t";
		}
		conc_field_txt << "\n";
	}

	conc_field_txt	<< "\nSimulated time\t[s]  = " << *calculated_total_time << "\n"
					<< "\Simulated time\t[d]  = " << *calculated_total_time / 86400 << "\n"
					<< "current dt\t[s]  = " << dt / dt_multiplier << "\n"
					<< "next dt\t\t[s]  = " << dt << "\n"
					<< "Total length\t[mm] = " << dx * no_of_cells_x * 1000;

	conc_field_txt.close();
}
//*********************************
void export_fraction_to_file()
{
	conc_field_txt.open(txt_filename_fraction, ios::out);
	conc_field_txt.precision(5);

	for (int j = 0; j < no_of_cells_y; j++)
	{
		for (int i = 0; i < no_of_cells_x; i++)
		{
			conc_field_txt << Ctab[j][i].phase_fraction << "\t";
			//conc_field_txt << Ctab[j][i].cell_state << "\t";
			//conc_field_txt << Ctab[j][i].interface_distance << "\t";
		}
		conc_field_txt << "\n";
	}

	conc_field_txt.close();
}
//*********************************
void interface_distance_setting(int *distance, int *j, int *i)
{
	if (Ctab[*j][*i].cell_state == C_cell::list_of_states::INTERFACE)
	{
		*distance = 0;
	}

	if (*distance < Ctab[*j][*i].interface_distance)
	{
		Ctab[*j][*i].interface_distance = *distance;
	}
	else
	{
		*distance = Ctab[*j][*i].interface_distance;
	}

	if (*distance < (C_cell::fixed_calculations_distance - 1) )
	{
		(*distance)++;
	}

	//cout << Ctab[*j][*i].interface_distance << "\t";	
}
//*********************************
void temp()
{
	for (int j = 0; j < no_of_cells_y; j++)
	{
		for (int i = 0; i < no_of_cells_x; i++)
		{
			cout << Ctab[j][i].interface_distance << "\t";	
		}
		cout << endl;
	}
	cout << endl;
}
//*********************************
void set_calculatoins_gap()
{
	int distance = C_cell::fixed_calculations_distance - 1;

	for (int j = 0; j < no_of_cells_y; j++)
	{
		distance = C_cell::fixed_calculations_distance - 1;

		for (int i = 0; i < no_of_cells_x; i++)
		{
			Ctab[j][i].interface_distance = C_cell::fixed_calculations_distance - 1;

			interface_distance_setting(&distance, &j, &i);
		}
	}
	//temp();

	//*********************************

	for (int j = 0; j < no_of_cells_y; j++)
	{
		distance = C_cell::fixed_calculations_distance - 1;

		for (int i = no_of_cells_x - 1; i >= 0; i--)
		{
			interface_distance_setting(&distance, &j, &i);
		}
	}
	//temp();

	//*********************************

	for (int i = 0; i < no_of_cells_x; i++)
	{
		distance = C_cell::fixed_calculations_distance - 1;

		for (int j = 0; j < no_of_cells_y; j++)
		{
			interface_distance_setting(&distance, &j, &i);
		}
	}
	//temp();
	
	//*********************************

	for (int i = 0; i < no_of_cells_x; i++)
	{
		distance = C_cell::fixed_calculations_distance - 1;

		for (int j = no_of_cells_y - 1; j >= 0; j--)
		{
			interface_distance_setting(&distance, &j, &i);
			Ctab[j][i].calculations_gap = &C_cell::tab_calculations_gaps[Ctab[j][i].interface_distance];
			if (*Ctab[j][i].calculations_gap < Ctab[j][i].calculations_countdown)
			{
				Ctab[j][i].calculations_countdown = *Ctab[j][i].calculations_gap;
			}
		}
	}
	//temp();
}
//*********************************
void export_header_file()
{
	set_calculatoins_gap();
	single_timestep_calculations_no_triggering();
	txt_filename_mass = "c:/Users/rz_m/Nextcloud/wydobywanie/wydobywanie/!output/";
	txt_filename_mass += "plik_!header";
	txt_filename_mass += ".ppm";
	cout << txt_filename_mass << endl;

	conc_field_txt.open(txt_filename_mass, ios::out);
	conc_field_txt.precision(10);
	conc_field_txt << no_of_cells_y << " = no of cells y" << endl;
	conc_field_txt << no_of_cells_x << " = no of cells x" << endl;
	conc_field_txt << dx << " = dx [m]" << endl;
	conc_field_txt << dx * no_of_cells_x << " = total length [m]" << endl;
	conc_field_txt << no_of_elements << " = no of elements" << endl;
	for (int i = 0; i < C_cell::no_of_elements; i++)
	{
		conc_field_txt << element[i] << ", ";
	}
	conc_field_txt  << " = elements" << endl;
	//conc_field_txt << no_of_file_loop_stepps * no_of_screen_loop_stepps << endl;

	for (int i = 0; i < no_of_elements - 1; i++)
	{
		conc_field_txt << Ctab[0][0].element_concentration[0][i] * Ctab[0][0].phase_fraction + Ctab[0][0].element_concentration[1][i] * (1 - Ctab[0][0].phase_fraction) << endl;
	}
	conc_field_txt << Ctab[0][0].mass[no_of_phases] << endl;

	for (int k = 0; k < no_of_elements - 1; k++)
	{
		for (int j = 0; j < no_of_cells_y; j++)
		{
			for (int i = 0; i < no_of_cells_x; i++)
			{
				if (Ctab[j][i].cell_state == C_cell::list_of_states::FCC)
				{
					conc_field_txt << Ctab[j][i].element_concentration[C_cell::list_of_states::FCC][k] * Ctab[j][i].phase_fraction << "\t";
				}
				else if (Ctab[j][i].cell_state == C_cell::list_of_states::LIQUID)
				{
					conc_field_txt << Ctab[j][i].element_concentration[C_cell::list_of_states::LIQUID][k] * (1 - Ctab[j][i].phase_fraction) << "\t";
				}
				else //INTERFACE
				{
					conc_field_txt << Ctab[j][i].element_concentration[C_cell::list_of_states::FCC][k] * Ctab[j][i].phase_fraction << "\t";
					conc_field_txt << Ctab[j][i].element_concentration[C_cell::list_of_states::LIQUID][k] * (1 - Ctab[j][i].phase_fraction) << "\t";
				}

				//conc_field_txt << Ctab[j][i].element_concentration[0][k] * Ctab[j][i].phase_fraction + Ctab[j][i].element_concentration[1][k] * (1 - Ctab[j][i].phase_fraction) << "\t";
			}
			conc_field_txt << "\n";
		}
	}

	for (int j = 0; j < no_of_cells_y; j++)
	{
		for (int i = 0; i < no_of_cells_x; i++)
		{
			conc_field_txt << Ctab[j][i].mass[no_of_phases] << "\t";
		}
		conc_field_txt << "\n";
	}

	conc_field_txt <<	no_of_screen_loop_stepps << " = screen stepps" << endl;
	conc_field_txt <<	no_of_file_loop_stepps << " = file stepps" << endl;
	conc_field_txt <<	no_of_dt_loop_stepps << " = dt stepps" << endl;
	conc_field_txt <<	no_of_calculations_loop_stepps << " = calculations stepps" << endl;
	conc_field_txt <<	temperature << " = temperature" << endl;

	for (int i = 0; i < no_of_elements; i++)
	{
		conc_field_txt << mobility[i] << " = mobility, element " << i << endl;
	}

	conc_field_txt.close();
}
//*********************************
void calculations_and_output()
{
	time(&start_time);

	export_header_file();
	// dodatek do test?w
	int a = -1;
	//int calculated_timestepps = 0;
	//const double continuation_time = 44537778;
	const double continuation_time = 0.0;
	double calculated_total_time = 0.0;
	set_file_name(&a);
	export_fraction_to_file();

	for (int timestep_screen = 0; timestep_screen < no_of_screen_loop_stepps; timestep_screen++)
	{
		cout << "timestep screen = " << timestep_screen << endl;
		
		for (int timestep_file = 0; timestep_file < no_of_file_loop_stepps; timestep_file++)
		{
			set_calculatoins_gap();

			for (int timestep_dt = 0; timestep_dt < no_of_dt_loop_stepps; timestep_dt++)
			{
				for (int timestep_calculations = 0; timestep_calculations < no_of_calculations_loop_stepps; timestep_calculations++)
				{
					single_timestep_calculations();
					update_mass_after_calculations();
					relax_stress();
					check_interface_vanishing();
				}
				//Ctab[0][100].present_yourself();

				dt *= dt_multiplier;
				if (dt > 0.05)
				{
					dt = 0.05;
				}

				dt_dx_squared = dt * dx_squared;
				//calculate_mobility_coefficients_table();
				calculate_mobility_coefficient();
			}

			//cout << Ctab[0][49].chemical_potential_coefficient << endl;
			//cout << Ctab[0][50].chemical_potential_coefficient << endl;
			//cout << "Fe 49: " << tab_delta_mass[1][49][0][2][0] << endl;
			//cout << "Fe 50: " << tab_delta_mass[1][50][0][2][0] << endl;
			//cout << "Fe 51: " << tab_delta_mass[1][51][0][2][0] << endl;
			//calculated_timestepps += no_of_calculations_loop_stepps;
			calculated_total_time = continuation_time + initial_dt * (1 - pow(dt_multiplier, 0 + no_of_dt_loop_stepps * (timestep_screen * no_of_file_loop_stepps + timestep_file + 1) ) ) / (1 - dt_multiplier) * no_of_calculations_loop_stepps;
			//calculated_total_time = continuation_time + (timestep_screen + 1) * (timestep_file + 1) * (no_of_dt_loop_stepps + 1) - 1;

			int timestep = timestep_screen * no_of_file_loop_stepps + timestep_file;
			set_file_name(&timestep);
			export_mass_to_file(&calculated_total_time);
			export_fraction_to_file();
		}
	}



	time(&end_time);
}
//*********************************
void export_summary_file()
{
	string local_txt_filename;

	local_txt_filename = "c:/Users/rz_m/Nextcloud/wydobywanie/wydobywanie/!output/";
	local_txt_filename += "plik_!summary";
	local_txt_filename += ".txt";
	cout << local_txt_filename << endl;

	summary_file.open(local_txt_filename, ios::out);
	summary_file.precision(1);
	summary_file << scientific;

	summary_file << (int) difftime(end_time, start_time) << "\t [sec]\t= calculations duration" << "\n";
	summary_file << (double) (no_of_calculations_loop_stepps * no_of_file_loop_stepps * no_of_screen_loop_stepps * no_of_cells_x * no_of_cells_y) << " [#]\t= stepps * X * Y" << "\n\n";

	summary_file << no_of_calculations_loop_stepps << "\t [#]\t= calculations stepps" << endl;
	summary_file << no_of_file_loop_stepps << "\t [#]\t= file stepps" << endl;
	summary_file << no_of_screen_loop_stepps << "\t [#]\t= screen stepps" << endl;
	summary_file << (double) (no_of_calculations_loop_stepps * no_of_file_loop_stepps * no_of_screen_loop_stepps) << " [#]\t= total stepps" << endl << endl;

	summary_file << no_of_cells_x << "\t [#]\t= number of X axis cells" << endl;
	summary_file << no_of_cells_y << "\t [#]\t= number of Y axis cells" << endl << endl;

	summary_file << no_of_mass_errors << "\t [#]\t= mass miscalculations" << endl;
	summary_file.close();
}

//--------------------------------------------------------------------
//----------------------------   MAIN   ------------------------------

int main()
{
	cout << "Czesc Marek\n" << endl;

	//C_TQ_chem_pot obiekt(phase_zero,  "Bi");
	//double tab[] = { 10.0, 10.0 };
	//obiekt.enter_conditions_ele_mass(tab);
	//obiekt.enter_condition_temperature(300.0);
	//obiekt.compute_equilibrium();
	//obiekt.print_equilibrium();

	//obiekt.enter_condition_temperature(320.0);
	//obiekt.compute_equilibrium();
	//obiekt.print_equilibrium();

	//chemical_potential_pointer = obiekt.get_chemical_potentials(300.0, 10.0, 10.0);
	//cout << *chemical_potential_pointer;

	

	initialization();
	calculations_and_output();
	export_summary_file();


	// Z PROMOTOREM
	// 21-04-2020 sprawdzi? jednostki zmiany masy, po pierwszym kroku masa staje si? ujemna, pomys?: gramy a kilogramy w TQ, mobilno??
	// 21-05-2020 update: zmiana masy dzia?a ale poprzez sztuczn?, niesprawdzon? zmian? mobilno?ci, zweryfikowa? warto?ci mobilno?ci w oparciu o wzory (promotor) oraz odnie?? si? do uwagi "gramy na kg w TQ"

	// 14-07-2020 jak ma si? rzeczywista pr?dko?? interfejsu do rzeczywistej pr?dko?ci dyfuzji atomu w gradiencie potencja??w chemicznych? - dzi?ki temu odpowied? na pytanie o krotno?? mobilno?ci na froncie 
	
	// SUGESTIE
	// 30-06-2020 total interface length?
	// 13-07-2020 mo?e check interface vanishing robi? ?adziej, wszystkie sprawdzania ruchu interfejsu?
	// 13-07-2020 WA?NE - jest triggerowanie przy przechodzeniu r?wnoleg?ym interfejsu, problem si? pojawi, je?li triggerowana kom?rka b?dzie immune, mo?e doda? sprawdzenie co jaki? czas czy nie trzeba "dotriggerowa?"

	delete [] Ctab;
	//system("c:/temp/123.txt");
	//system("c:/Users/rz_m/Nextcloud/wydobywanie/wydobywanie/!output/plik_!summarytxt");
	system("pause");
}

// log:
// 2021-05-06: problem z nieustalaniem r?wnej masy, rozwi?zanie: zwi?kszy? czu?o?? na odchylenie od r?wnowagi przy mu_coef: zwi?kszy? mno?nik ALE to powoduje powstawanie fal na masie kom?rki,
// albo zmniejszy? krok albo inaczej rozwi?za? problem
// 2021-05-06: podj?to rozwi?zanie ze zmniejszeniem ktoku czaswego - sukces, mo?na wyliczy? obszar i czas oblicze?!
// 2021-05-13: pierwsze obliczenia czasu: 5 min oblicze? 2.4e-5 sec czasu rzeczywistego 
// dx-7 dt-6
// 2021-05-13: zrobiono rzeczywiste mobilno?ci, powalczy? z czasem oblicze? (zmienne dt)
// 2021-05-14: dt = 5e-6 wygl?da, ?e dzia?a
// 2021-05-19 wyniki ze sta?ym dt=1e-5 i ze zmiennym nie pokrywaj? si?, przetestuj? wyniki dla sta?ego dt 1e-5 i 5e-6 czy s? identyczne, trzeba znale?? maksymalne dt, kt?re daje dobre wyniki 
// po sprawdzeniu dt = 1e-5, 5e-6, 2.5e-6 : wyniki bardzo do siebie zbli?one, przy sta?ym 1e-5 problemu nie ma ale mo?e on by? szybko zauwa?alny dla niewiele wzrastaj?cego dt ju? w pierwszych krokach
// dlatego zacz?? nie od 1e-5 a od 5e-6
// q=1,0046 jest wystarczaj?ce po por?wnaniu z q=1,0025
// 2021-05-19 wznawianie oblicze? dla var dt dzia? bez zarzutu
// 2021-05-24 dodano kolejn? pentl? redukuj?c? ilo?? plik?w a nie sprawiaj?c?, ?e rzadziej dochodzi do zwi?kszenia dt