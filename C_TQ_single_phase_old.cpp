#include "C_TQ_single_phase_old.h"
//********************************
C_TQ_single_phase_old::C_TQ_single_phase_old(char *phase, char *ele0, char *ele1, char *ele2, char *ele3,
				 char *ele4, char *ele5, char *ele6, char *ele7, char *ele8)
				 : C_TQ_Fe_old(ele0, ele1, ele2, ele3, ele4, ele5, ele6, ele7, ele8)
{
	isolate_stable_phase(phase);
}
//********************************
C_TQ_single_phase_old::~C_TQ_single_phase_old()
{
}
//*********************************
void C_TQ_single_phase_old::isolate_stable_phase(char *stable_phase)
{
	suspend_all_phases();
	get_phase_index(stable_phase);
	tq_csp(phase_index, "entered", 0.0, iwsg, iwse);
}
//*********************************
void C_TQ_single_phase_old::suspend_all_phases()
{
	for(int i = 0; i < number_of_phases; i++)
	{
		tq_csp(i + 1, "sus", 0.0, iwsg, iwse);
	}
}