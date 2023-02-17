#pragma once

#include "C_TQ_chem_pot.h"	// <------------------------------------ sprawdziæ po co ten include 09-04-2019 na koñcu pliku jest wskaŸnik typu chem_pot wiêc pewnie dlatego
//********************************
class C_cell
{
public:
	enum				  {	no_of_elements = 3};		// iloœæ pierwiastków w programie, bardzo wa¿ne miejsce
	enum				  {	no_of_phases = 2};
	enum				  { fixed_calculations_distance	= 27};
	enum list_of_states	  {	FCC = 0,
							LIQUID = 1,
							INTERFACE = 2	};
	enum list_of_corners  {	NE = 0,
							SE = 1,
							SW = 2,
							NW = 3	};
	enum NEWS_tab		  {	NEWS_count = 4,
							NEWS_plus_interface_length = NEWS_count + 1,
							angle_count = 90 / 5,
							angle_count_minus_one = angle_count - 1,
							half_fill_count = 50,
							fill_count = 99	};
	//********************************
	C_cell();
	C_cell(const double *tab_conc, list_of_states *state, double *temper);
	~C_cell();
	//********************************
	inline void set_state(list_of_states *state);
	inline double hypotenuse(double *a, double *b);
	inline double hypotenuse(double *a);
	double mass_percent_to_kilogram(const double *mass_percent);
	double mass_percent_to_kilogram(const double *mass_percent, double *cell_mass);
	void set_initial_conditions(const double *tab_conc, list_of_states *state, double *temper);
	void set_initial_conditions(const double *tab_conc, list_of_states *state, double *solid_fraction, double *temper);
	void set_initial_conditions_mass_from_file(const double *tab_conc, list_of_states *state, double *solid_fraction, double *temper, double *cell_mass);
	void decrease_countdowns();
	void calculate_phase_fraction_and_masses();
	double calculate_phase_mass(double *tab_ele_mass);
	void calculate_concentration_and_chemical_potential_coefficient();
	void calculate_chemical_potential_coefficient();
	void set_NEWS_int_lenght_pointer();
	void fill_NEWS_int_length_table();
	void present_yourself();
	//void get_potentials();
	//void get_temp_volume();

	//********************************	
public:
	static const double	dx,
						volume,
						density,
						nominal_mass,
						fraction,
						tan_alpha_coefficient,
						liq_fcc_solid_fraction_on_edge[2];
	static double		tab_NEWS_int_length[angle_count][fill_count][NEWS_plus_interface_length];
	static const int	tab_calculations_gaps[fixed_calculations_distance];

	static bool	flag;

	bool	whether_to_check_edge[4] = {false, false, false, false}, //NEWS
			trigger[5] = {false, false, false, false, false};

	int		immunity_countdown		= 0,
			calculations_countdown	= 0,
			interface_distance		= 0;

	const int	*calculations_gap;
	
	double	temperature,
			phase_fraction,
			mass[no_of_phases + 1],									// indeks jak w enum list_of_states: FCC = 0, LIQ = 1, ostatnia wartosc: mass[no_of_phases] = masa calej komorki
			chemical_potential_coefficient,
			interface_numerator,
			interface_denominator,
			interface_normal;									// <5; 90> co 5

	const double	*pointer_solid_fraction_on_edge[4],
					*pointer_interface_length;

	double	(*element_mass)[no_of_elements],					// [kg]
			(*element_chemical_potential)[no_of_elements],		// []
			(*element_concentration)[no_of_elements];			// [kg / m3]

	list_of_states	cell_state;
	list_of_corners origin_corner;


	//C_TQ_chem_pot *system;
};