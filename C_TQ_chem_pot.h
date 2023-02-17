#pragma once

#include "C_TQ_single_phase.h"
//********************************
class C_TQ_chem_pot : public C_TQ_single_phase
{
public:
	C_TQ_chem_pot(char *phase, char *ele0, char *ele1 = 0, char *ele2 = 0, char *ele3 = 0, char *ele4 = 0,
					 char *ele5 = 0, char *ele6 = 0, char *ele7 = 0, char *ele8 = 0);
	C_TQ_chem_pot() {  cout << "Konstruktor domniemany klasy C_TQ_chem_pot" << endl; };
	~C_TQ_chem_pot();
	//********************************
	double* get_chemical_potentials(double temper, TC_FLOAT ele0_wtpc,
		TC_FLOAT ele1_wtpc = 0, TC_FLOAT ele2_wtpc = 0, TC_FLOAT ele3_wtpc = 0,
		TC_FLOAT ele4_wtpc = 0, TC_FLOAT ele5_wtpc = 0,	TC_FLOAT ele6_wtpc = 0,
		TC_FLOAT ele7_wtpc = 0, TC_FLOAT ele8_wtpc = 0);
	double* C_TQ_chem_pot::get_chemical_potentials(double *temper, double *tab_ele_mass);
	void C_TQ_chem_pot::calculations_with_print(double *temper, double *tab_ele_mass);
	void get_mass(double *vol);

	TC_INT tempest;
	//********************************
protected:
	TC_FLOAT *chemical_potentials_table;
	//********************************
	void retrieve_chemical_potentials();
};