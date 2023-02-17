#include "C_TQ_liquidus.h"
//********************************
C_TQ_liquidus::C_TQ_liquidus(char *ele0, char *ele1, char *ele2, char *ele3,
				 char *ele4, char *ele5, char *ele6, char *ele7, char *ele8)
				 : 	C_TQ_Fe_old(ele0, ele1, ele2, ele3, ele4, ele5, ele6, ele7, ele8)
{
	sulphur_present = check_sulphur_presence();
}
//********************************
C_TQ_liquidus::~C_TQ_liquidus()
{
}
//*********************************
bool C_TQ_liquidus::check_sulphur_presence()
{
	for(int i = 0; i < number_of_elements; i++)
	{
		if(table_of_elements[i][0] == 'S')
			if(!table_of_elements[i][1])
				return true;
	}

	return false;
}
//*********************************
double& C_TQ_liquidus::get_liquidus(TC_FLOAT ele0_wtpc, TC_FLOAT ele1_wtpc,
	TC_FLOAT ele2_wtpc, TC_FLOAT ele3_wtpc, TC_FLOAT ele4_wtpc, TC_FLOAT ele5_wtpc, TC_FLOAT ele6_wtpc,
	TC_FLOAT ele7_wtpc, TC_FLOAT ele8_wtpc)
{
	enter_conditions_ele_wtpc(&ele0_wtpc, &ele1_wtpc, &ele2_wtpc, &ele3_wtpc,
		&ele4_wtpc, &ele5_wtpc, &ele6_wtpc, &ele7_wtpc, &ele8_wtpc);

	if(!sulphur_present)
	{
		set_liquid_fix_1();
		compute_equilibrium();
	}
	else
		find_liquidus_manually();

	tq_get1("T", -1, -1, &returned_liquidus, iwsg, iwse);
	print_equilibrium();

	return returned_liquidus;
}
//*********************************
void C_TQ_liquidus::set_liquid_fix_1()
{
	get_phase_index("LIQ");
	tq_csp(phase_index, "FIXED", 1.0, iwsg, iwse);
}
//*********************************
void C_TQ_liquidus::find_liquidus_manually()
{
	get_phase_index("LIQ");
	int no_stable_phases;
	temperature = 2000;
	delta_temp = 100;
	double driving_force, eps = 5e-3;
	bool liq_present, above_liquidus = true;

	while( (delta_temp > eps) || ( -delta_temp > eps) )
	{
		enter_condition_temperature(temperature);
		compute_equilibrium();
	
		if(delta_temp < 0)
			delta_temp *= -1;
		no_stable_phases = 0;
		liq_present = false;
		for(int i = 0; i < number_of_phases; i++)
		{
			tq_get1("DG", i + 1, -1, &driving_force, iwsg, iwse);
			if(!driving_force)
			{
				no_stable_phases++;
				if(phase_index == i + 1)
					liq_present = true;
			}
		}
		if( (liq_present) && (no_stable_phases == 1) )
		{
			if(above_liquidus)
				delta_temp /= -1;
			else
				delta_temp /= -2;
		}
		else
		{
			delta_temp /= 2;
			above_liquidus = false;
		}

		temperature += delta_temp;
	}
}