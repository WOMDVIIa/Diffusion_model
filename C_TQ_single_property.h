#pragma once

#include "C_TQ_single_phase.h"
//********************************
class C_TQ_single_property : public C_TQ_single_phase
{
public:
	C_TQ_single_property(char *mnem_prop, char *phase, char *ele0, char *ele1 = 0, char *ele2 = 0, char *ele3 = 0, char *ele4 = 0,
					 char *ele5 = 0, char *ele6 = 0, char *ele7 = 0, char *ele8 = 0);
	C_TQ_single_property() {  cout << "Konstruktor domniemany klasy C_TQ_single_property" << endl; };
	~C_TQ_single_property();
	//********************************
	double& get_value(TC_FLOAT temp, TC_FLOAT ele0_wtpc, TC_FLOAT ele1_wtpc = 0,
		TC_FLOAT ele2_wtpc = 0, TC_FLOAT ele3_wtpc = 0, TC_FLOAT ele4_wtpc = 0, TC_FLOAT ele5_wtpc = 0,
		TC_FLOAT ele6_wtpc = 0, TC_FLOAT ele7_wtpc = 0, TC_FLOAT ele8_wtpc = 0);
//********************************
protected:
	char *mnemotic_of_property;
	TC_FLOAT returned_value;

	double value;
};