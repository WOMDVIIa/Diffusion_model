#pragma once

#include "C_TQ_Fe.h"
//********************************
class C_TQ_chem_pot_two_ph : protected C_TQ_Fe
{
public:
	C_TQ_chem_pot_two_ph(char *ele0, char *ele1 = 0, char *ele2 = 0, char *ele3 = 0, char *ele4 = 0,
					 char *ele5 = 0, char *ele6 = 0, char *ele7 = 0, char *ele8 = 0);
	C_TQ_chem_pot_two_ph() { cout << "Konstruktor domniemany klasy C_TQ_chem_pot" << endl; };
	~C_TQ_chem_pot_two_ph();
	//********************************
	double* get_chemical_potentials(double temp, TC_FLOAT ele0_wtpc,
		TC_FLOAT ele1_wtpc = 0, TC_FLOAT ele2_wtpc = 0, TC_FLOAT ele3_wtpc = 0,
		TC_FLOAT ele4_wtpc = 0, TC_FLOAT ele5_wtpc = 0,	TC_FLOAT ele6_wtpc = 0,
		TC_FLOAT ele7_wtpc = 0, TC_FLOAT ele8_wtpc = 0);
//********************************
protected:
	TC_FLOAT *chemical_potentials_table;
	//********************************
	void retrieve_chemical_potentials();
};