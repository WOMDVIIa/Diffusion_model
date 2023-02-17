#define WIN64
#pragma once

#include <iostream>
#include <string>
#include "tqroot.h"

using namespace std;
//********************************
class C_TQ_Fe_old
{
public:
	C_TQ_Fe_old(char *ele0, char *ele1 = 0, char *ele2 = 0, char *ele3 = 0, char *ele4 = 0,
					 char *ele5 = 0, char *ele6 = 0, char *ele7 = 0, char *ele8 = 0);
	C_TQ_Fe_old() {  cout << "Konstruktor domniemany klasy C_TQ_Fe" << endl; };
	~C_TQ_Fe_old();
	//********************************
	void enter_condition_temperature(TC_FLOAT temp);
	void enter_condition_temperature(TC_FLOAT *temp);
	void enter_conditions_ele_wtpc(TC_FLOAT ele0_wtpc, TC_FLOAT ele1_wtpc = 0,
		TC_FLOAT ele2_wtpc = 0, TC_FLOAT ele3_wtpc = 0, TC_FLOAT ele4_wtpc = 0, TC_FLOAT ele5_wtpc = 0,
		TC_FLOAT ele6_wtpc = 0, TC_FLOAT ele7_wtpc = 0, TC_FLOAT ele8_wtpc = 0);
	void enter_conditions_ele_wtpc(TC_FLOAT *ele0_wtpc, TC_FLOAT *ele1_wtpc = 0,
		TC_FLOAT *ele2_wtpc = 0, TC_FLOAT *ele3_wtpc = 0, TC_FLOAT *ele4_wtpc = 0, TC_FLOAT *ele5_wtpc = 0,
		TC_FLOAT *ele6_wtpc = 0, TC_FLOAT *ele7_wtpc = 0, TC_FLOAT *ele8_wtpc = 0);
	//********************************
	void compute_equilibrium();	
	void print_equilibrium();
	void reinitiate();

	//double& get_value(TC_FLOAT temp, TC_INT no_of_phase, TC_FLOAT ele0_wtpc = 3.5, TC_FLOAT ele1_wtpc = 0,		// ------------------------------------------------- dodano do prog Jacka
	//	TC_FLOAT ele2_wtpc = 0, TC_FLOAT ele3_wtpc = 0, TC_FLOAT ele4_wtpc = 0, TC_FLOAT ele5_wtpc = 0,
	//	TC_FLOAT ele6_wtpc = 0, TC_FLOAT ele7_wtpc = 0, TC_FLOAT ele8_wtpc = 0);
//********************************
protected:
	static string separator;
	TC_INT *iwsg, *iwse, ierr, number_of_elements, number_of_phases;
	
	static TC_STRING used_database;
	static char	*tc_installation_directory,
				*log_file_directory;

	char* *table_of_elements;
	//char  (*table_of_phases)[TC_STRLEN_PHASES];

	TC_FLOAT	*element_wtpc,
				temperature;

	static TC_FLOAT	n, p;

	TC_INT		*element_index,
				//*phase_index,			//	sprawdziæ
				*num_condition_wtp,
				num_condition_t,
				num_condition_n,
				num_condition_p,
				phase_index;			//	sprawdziæ
	//********************************
	void initialization();
	void db_open();
	void common_constructors_part();
		void input_elements_to_the_system();
		void set_system_phases();
		void get_data();
		void get_elements_index();
		void get_phases_number();
		void set_cond_p_n();
	
	void set_conditions_ele_wtpc();
	void get_phase_index(char *ph_name);

	//void remove_condition_temperature();

	//TC_FLOAT returned_value; // ------------------------------------------------- dodano do prog Jacka
};