#include "C_TQ_single_property_old.h"
//********************************
C_TQ_single_property_old::C_TQ_single_property_old(char *mnem_prop, char *phase, char *ele0, char *ele1, char *ele2, char *ele3,
				 char *ele4, char *ele5, char *ele6, char *ele7, char *ele8)
				 : C_TQ_single_phase_old(phase, ele0, ele1, ele2, ele3, ele4, ele5, ele6, ele7, ele8)
{
	mnemotic_of_property = mnem_prop;
}
//********************************
C_TQ_single_property_old::~C_TQ_single_property_old()
{
}
//********************************
double& C_TQ_single_property_old::get_value(TC_FLOAT temp, TC_FLOAT ele0_wtpc, TC_FLOAT ele1_wtpc,
	TC_FLOAT ele2_wtpc, TC_FLOAT ele3_wtpc, TC_FLOAT ele4_wtpc, TC_FLOAT ele5_wtpc, TC_FLOAT ele6_wtpc,
	TC_FLOAT ele7_wtpc, TC_FLOAT ele8_wtpc)
{
	enter_condition_temperature(temp);
	enter_conditions_ele_wtpc(&ele0_wtpc, &ele1_wtpc, &ele2_wtpc, &ele3_wtpc,
		&ele4_wtpc, &ele5_wtpc, &ele6_wtpc, &ele7_wtpc, &ele8_wtpc);
	
	compute_equilibrium();
	//print_equilibrium();

	tq_get1("V", -1, -1, &volume, iwsg, iwse);
	tq_get1("M", -1, -1, &mass, iwsg, iwse);

	//cout << "mass [kg] = " << mass << endl;
	//cout << "vol [m3] = " << volume << endl;

	returned_value = mass / volume;

	return returned_value;
}