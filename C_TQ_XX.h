#define		WIN64
#define		main_element		"Fe"
#define		phase_zero			"BCC"
#define		phase_one			"LIQ"
#pragma once

#include <iostream>
#include <string>
#include "tqroot.h"

using namespace std;
//********************************
class C_TQ_XX
{
public:
	static string separator;

	tc_elements_strings* list_ele;
	TC_INT* number_ele;

	C_TQ_XX(char* ele0, char* ele1 = 0, char* ele2 = 0, char* ele3 = 0, char* ele4 = 0,
		char* ele5 = 0, char* ele6 = 0, char* ele7 = 0, char* ele8 = 0);
	C_TQ_XX() { cout << "Konstruktor domniemany klasy C_TQ_XX" << endl; };
	~C_TQ_XX();
	//********************************
	void enter_condition_temperature(TC_FLOAT temp);
	void enter_condition_temperature(TC_FLOAT* temp);
	void enter_conditions_ele_mass(TC_FLOAT ele0_mass, TC_FLOAT ele1_mass,
		TC_FLOAT ele2_mass = 0, TC_FLOAT ele3_mass = 0, TC_FLOAT ele4_mass = 0, TC_FLOAT ele5_mass = 0,
		TC_FLOAT ele6_mass = 0, TC_FLOAT ele7_mass = 0, TC_FLOAT ele8_mass = 0);
	void enter_conditions_ele_mass(TC_FLOAT* ele0_mass, TC_FLOAT* ele1_mass,
		TC_FLOAT* ele2_mass = 0, TC_FLOAT* ele3_mass = 0, TC_FLOAT* ele4_mass = 0, TC_FLOAT* ele5_mass = 0,
		TC_FLOAT* ele6_mass = 0, TC_FLOAT* ele7_mass = 0, TC_FLOAT* ele8_mass = 0);
	void enter_conditions_ele_mass(TC_FLOAT* tab_ele_mass);
	//********************************
	void compute_equilibrium();
	TC_INT compute_equilibrium(TC_INT a);
	void print_equilibrium();
	void reinitiate();

	//double& get_value(TC_FLOAT temp, TC_INT no_of_phase, TC_FLOAT ele0_mass = 3.5, TC_FLOAT ele1_mass = 0,		// ------------------------------------------------- dodano do prog Jacka
	//	TC_FLOAT ele2_mass = 0, TC_FLOAT ele3_mass = 0, TC_FLOAT ele4_mass = 0, TC_FLOAT ele5_mass = 0,
	//	TC_FLOAT ele6_mass = 0, TC_FLOAT ele7_mass = 0, TC_FLOAT ele8_mass = 0);
//********************************
protected:
	TC_INT* iwsg, * iwse, ierr, number_of_elements, number_of_phases;

	


	static TC_STRING used_database;
	static char* tc_installation_directory,
		* log_file_directory;

	char** table_of_elements;
	//char  (*table_of_phases)[TC_STRLEN_PHASES];

	TC_FLOAT* element_mass,
		temperature;

	static TC_FLOAT	n, p;

	TC_INT* element_index,
		//*phase_index,			//	09-04-2019 nie wiem czy potrzebne
		* num_condition_mass,
		num_condition_t,
		//num_condition_n,
		num_condition_p,
		phase_index;			//	sprawdziæ: 09-04-2019 u¿ywane do pozostawienia jednej fazy "ENTERED" w C_TQ_single_phase
//********************************
	void database_list();
	void initialization();
	void db_open();
	void common_constructors_part();
	void input_elements_to_the_system();
	void set_system_phases();
	void get_data();
	void get_elements_index();
	void get_phases_number();
	void set_cond_p();															// zmiana kodu przy przekszta³ceniu warunków z N = 1, W% na masy wszystkich (³¹cznie z XX) pierwiastków

	void set_conditions_ele_mass();
	void get_phase_index(char* ph_name);

	//void remove_condition_temperature();

	//TC_FLOAT returned_value; // ------------------------------------------------- dodano do prog Jacka
};