#pragma once

#include "C_TQ_Fe_old.h"
//********************************
class C_TQ_liquidus : protected C_TQ_Fe_old
{
public:
	C_TQ_liquidus(char *ele0, char *ele1 = 0, char *ele2 = 0, char *ele3 = 0, char *ele4 = 0,
					 char *ele5 = 0, char *ele6 = 0, char *ele7 = 0, char *ele8 = 0);
	C_TQ_liquidus() {  cout << "Konstruktor domniemany klasy C_TQ_liquidus" << endl; };
	~C_TQ_liquidus();
	//********************************
	double& get_liquidus(TC_FLOAT ele0_wtpc, TC_FLOAT ele1_wtpc = 0,
		TC_FLOAT ele2_wtpc = 0, TC_FLOAT ele3_wtpc = 0, TC_FLOAT ele4_wtpc = 0, TC_FLOAT ele5_wtpc = 0,
		TC_FLOAT ele6_wtpc = 0, TC_FLOAT ele7_wtpc = 0, TC_FLOAT ele8_wtpc = 0);
//********************************
protected:
	TC_FLOAT returned_liquidus;
	TC_FLOAT delta_temp;
	bool sulphur_present;
	//********************************
	bool check_sulphur_presence();
	void set_liquid_fix_1();
	void find_liquidus_manually();
};