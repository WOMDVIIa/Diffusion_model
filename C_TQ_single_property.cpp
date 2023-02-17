#include "C_TQ_single_property.h"
//********************************
C_TQ_single_property::C_TQ_single_property(char *mnem_prop, char *phase, char *ele0, char *ele1, char *ele2, char *ele3,
				 char *ele4, char *ele5, char *ele6, char *ele7, char *ele8)
				 : 	C_TQ_single_phase(phase, ele0, ele1, ele2, ele3, ele4, ele5, ele6, ele7, ele8)
{
	mnemotic_of_property = mnem_prop;
}
//********************************
C_TQ_single_property::~C_TQ_single_property()
{
}
//********************************
double& C_TQ_single_property::get_value(TC_FLOAT temp, TC_FLOAT ele0_mass, TC_FLOAT ele1_mass,
	TC_FLOAT ele2_mass, TC_FLOAT ele3_mass, TC_FLOAT ele4_mass, TC_FLOAT ele5_mass, TC_FLOAT ele6_mass,
	TC_FLOAT ele7_mass, TC_FLOAT ele8_mass)
{
	enter_condition_temperature(temp);
	enter_conditions_ele_mass(&ele0_mass, &ele1_mass, &ele2_mass, &ele3_mass,
		&ele4_mass, &ele5_mass, &ele6_mass, &ele7_mass, &ele8_mass);
	
	compute_equilibrium();
	print_equilibrium();

	tq_get1(mnemotic_of_property, -1, -1, &returned_value, iwsg, iwse);
	//tq_get1("V", -1, -1, &returned_value, iwsg, iwse);
	//returned_value = returned_value / value;
	return returned_value;
}