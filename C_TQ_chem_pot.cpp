#include "C_TQ_chem_pot.h"
//********************************
C_TQ_chem_pot::C_TQ_chem_pot(char *phase, char *ele0, char *ele1, char *ele2, char *ele3,
				 char *ele4, char *ele5, char *ele6, char *ele7, char *ele8)
				 : 	C_TQ_single_phase(phase, ele0, ele1, ele2, ele3, ele4, ele5, ele6, ele7, ele8)
{
	chemical_potentials_table = new TC_FLOAT[number_of_elements];
}
//********************************
C_TQ_chem_pot::~C_TQ_chem_pot()
{
	delete [] chemical_potentials_table;
}
//********************************
double* C_TQ_chem_pot::get_chemical_potentials(double temper, TC_FLOAT ele0_mass,
		TC_FLOAT ele1_mass, TC_FLOAT ele2_mass, TC_FLOAT ele3_mass,	TC_FLOAT ele4_mass,
		TC_FLOAT ele5_mass,	TC_FLOAT ele6_mass,	TC_FLOAT ele7_mass, TC_FLOAT ele8_mass)
{
	enter_condition_temperature(&temper);

	enter_conditions_ele_mass(&ele0_mass, &ele1_mass, &ele2_mass, &ele3_mass,		// masa wprowadzona w kg przeliczona przez tq na g
		&ele4_mass, &ele5_mass, &ele6_mass, &ele7_mass, &ele8_mass);
	
	compute_equilibrium();
	//print_equilibrium();
	retrieve_chemical_potentials();

	return chemical_potentials_table;
}
//*********************************
double* C_TQ_chem_pot::get_chemical_potentials(double *temper, double *tab_ele_mass)
{
	enter_condition_temperature(temper);
	enter_conditions_ele_mass(tab_ele_mass);

	////tu dopisane do testów
	//tempest = compute_equilibrium(tempest);
	//if (tempest != 0)
	//{
	//	cout << "test" << endl;
	//	return 0;
	//}

	compute_equilibrium();
	//print_equilibrium();	
	retrieve_chemical_potentials();

	return chemical_potentials_table;
}
//*********************************
void C_TQ_chem_pot::calculations_with_print(double *temper, double *tab_ele_mass)
{
	enter_condition_temperature(temper);
	enter_conditions_ele_mass(tab_ele_mass);
	compute_equilibrium();
	print_equilibrium();
	//retrieve_chemical_potentials();
}	
//********************************
void C_TQ_chem_pot::retrieve_chemical_potentials()
{
	for(int i = 0; i < number_of_elements; i++)
	{
		tq_get1("MU", -1, element_index[i], &chemical_potentials_table[i], iwsg, iwse);
	}
}
//*********************************
void C_TQ_chem_pot::get_mass(double *vol)
{
	tq_get1("M", -1, -1, vol, iwsg, iwse);
}