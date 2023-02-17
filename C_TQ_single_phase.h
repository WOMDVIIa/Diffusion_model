#pragma once

#include "C_TQ_XX.h"
//********************************
class C_TQ_single_phase : public C_TQ_XX
{
public:
	C_TQ_single_phase(char *phase, char *ele0, char *ele1 = 0, char *ele2 = 0, char *ele3 = 0, char *ele4 = 0,
					 char *ele5 = 0, char *ele6 = 0, char *ele7 = 0, char *ele8 = 0);
	C_TQ_single_phase() {  cout << "Konstruktor domniemany klasy C_TQ_single_phase" << endl; };
	~C_TQ_single_phase();
//********************************
	void isolate_stable_phase(char *stable_phase);
	void suspend_all_phases();
};