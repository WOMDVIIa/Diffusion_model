#pragma once

#include "C_TQ_Fe_old.h"
//********************************
class C_TQ_single_phase_old : public C_TQ_Fe_old
{
public:
	C_TQ_single_phase_old(char *phase, char *ele0, char *ele1 = 0, char *ele2 = 0, char *ele3 = 0, char *ele4 = 0,
					 char *ele5 = 0, char *ele6 = 0, char *ele7 = 0, char *ele8 = 0);
	C_TQ_single_phase_old() {  cout << "Konstruktor domniemany klasy C_TQ_single_phase" << endl; };
	~C_TQ_single_phase_old();
//********************************
	void isolate_stable_phase(char *stable_phase);
	void suspend_all_phases();
};