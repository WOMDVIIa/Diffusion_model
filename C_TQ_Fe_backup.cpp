#include "C_TQ_Fe.h"
//*********************************
string C_TQ_Fe::separator(50, '-');
char* C_TQ_Fe::tc_installation_directory = "C:\\Program Files (x86)\\Thermo-Calc\\Thermo-Calc\\tc\\vers";
char* C_TQ_Fe::log_file_directory = "C:\\Users\\Marek\\AppData\\Local\\Temp\\!tq";
TC_STRING C_TQ_Fe::used_database = "TCFE5";
TC_FLOAT C_TQ_Fe::n = 1;
TC_FLOAT C_TQ_Fe::p = 101325;
//*********************************
C_TQ_Fe::C_TQ_Fe(char *ele0, char *ele1, char *ele2, char *ele3, char *ele4,
				 char *ele5, char *ele6, char *ele7, char *ele8)
{
	initialization();
	db_open();

	if(!ele1)
	{
		number_of_elements = 2;
		table_of_elements = new char*[number_of_elements];
		table_of_elements[0] = ele0;
		table_of_elements[1] = "Fe";
	}
	else if(!ele2)
	{
		number_of_elements = 3;
		table_of_elements = new char*[number_of_elements];
		table_of_elements[0] = ele0;
		table_of_elements[1] = ele1;
		table_of_elements[2] = "Fe";		
	}
	else if(!ele3)
	{
		number_of_elements = 4;
		table_of_elements = new char*[number_of_elements];
		table_of_elements[0] = ele0;
		table_of_elements[1] = ele1;
		table_of_elements[2] = ele2;
		table_of_elements[3] = "Fe";
	}
	else if(!ele4)
	{
		number_of_elements = 5;
		table_of_elements = new char*[number_of_elements];
		table_of_elements[0] = ele0;
		table_of_elements[1] = ele1;
		table_of_elements[2] = ele2;
		table_of_elements[3] = ele3;
		table_of_elements[4] = "Fe";		
	}
	else if(!ele5)
	{
		number_of_elements = 6;
		table_of_elements = new char*[number_of_elements];
		table_of_elements[0] = ele0;
		table_of_elements[1] = ele1;
		table_of_elements[2] = ele2;
		table_of_elements[3] = ele3;
		table_of_elements[4] = ele4;
		table_of_elements[5] = "Fe";		
	}
	else if(!ele6)
	{
		number_of_elements = 7;
		table_of_elements = new char*[number_of_elements];
		table_of_elements[0] = ele0;
		table_of_elements[1] = ele1;
		table_of_elements[2] = ele2;
		table_of_elements[3] = ele3;
		table_of_elements[4] = ele4;
		table_of_elements[5] = ele5;
		table_of_elements[6] = "Fe";		
	}
	else if(!ele7)
	{
		number_of_elements = 8;
		table_of_elements = new char*[number_of_elements];
		table_of_elements[0] = ele0;
		table_of_elements[1] = ele1;
		table_of_elements[2] = ele2;
		table_of_elements[3] = ele3;
		table_of_elements[4] = ele4;
		table_of_elements[5] = ele5;
		table_of_elements[6] = ele6;
		table_of_elements[7] = "Fe";		
	}
	else if(!ele8)
	{
		number_of_elements = 9;
		table_of_elements = new char*[number_of_elements];
		table_of_elements[0] = ele0;
		table_of_elements[1] = ele1;
		table_of_elements[2] = ele2;
		table_of_elements[3] = ele3;
		table_of_elements[4] = ele4;
		table_of_elements[5] = ele5;
		table_of_elements[6] = ele6;
		table_of_elements[7] = ele7;
		table_of_elements[8] = "Fe";		
	}
	else
	{
		number_of_elements = 10;
		table_of_elements = new char*[number_of_elements];
		table_of_elements[0] = ele0;
		table_of_elements[1] = ele1;
		table_of_elements[2] = ele2;
		table_of_elements[3] = ele3;
		table_of_elements[4] = ele4;
		table_of_elements[5] = ele5;
		table_of_elements[6] = ele6;
		table_of_elements[7] = ele7;
		table_of_elements[8] = ele8;
		table_of_elements[9] = "Fe";		
	}

	common_constructors_part();
}
//*********************************
//C_TQ_Fe::C_TQ_Fe(int n_ele, char *elements[]) : number_of_elements(n_ele + 1)
//{
//	initialization();
//	db_open();
//	//*********************************
//	//Getting elements list
//	table_of_elements = new char*[number_of_elements]; //Pointer to the table of pointers
//	for(int i = 0; i < number_of_elements - 1; i++)
//	{
//		table_of_elements[i] = elements[i];
//	}
//	table_of_elements[number_of_elements - 1] = "Fe";
//	//*********************************
//	common_constructors_part();
//}
//*********************************
C_TQ_Fe::~C_TQ_Fe()
{
	delete [] num_condition_wtp;
	delete [] element_index;
	delete [] element_wtpc;
	delete [] table_of_elements;
	delete [] iwsg;
	delete [] iwse;
}
//*********************************
void C_TQ_Fe::initialization()
{
	iwsg = new TC_INT[TC_NWSG];		//TQ workspace
	iwse = new TC_INT[TC_NWSE];		//TQ workspace
	//*********************************
	tq_ini3(tc_installation_directory, log_file_directory, TC_NWSG, TC_NWSE, iwsg, iwse);
	if(tq_sg2err(&ierr) )
	{
		cout << "error initializing \n";
		system("pause");
	}
}
//*********************************
void C_TQ_Fe::db_open()
{
	tq_opdb(used_database, iwsg, iwse);
}
//*********************************
void C_TQ_Fe::common_constructors_part()
{
	input_elements_to_the_system();
	set_system_phases();
	get_data();
	//*********************************
	//Input of weight percent and other conditions
	element_wtpc = new TC_FLOAT[number_of_elements - 1];
	element_index = new TC_INT[number_of_elements];
	num_condition_wtp = new TC_INT[number_of_elements - 1];

	get_elements_index();
	get_phases_number();
	set_cond_p_n();
}
//*********************************
void C_TQ_Fe::input_elements_to_the_system()
{
	for(int i = 0; i < number_of_elements; i++)
	{
		tq_defel(table_of_elements[i], iwsg, iwse);
		if(tq_sg2err(&ierr) )
		{
			cout << "There is no such element! Try again. \n";
			tq_reserr();
		}
	}
}
//*********************************
void C_TQ_Fe::set_system_phases()
{
	tq_rejph("*", iwsg, iwse);
	tq_resph("FCC LIQ BCC GRA CEM", iwsg, iwse);
	//tq_resph("LIQ GRA", iwsg, iwse);
}
//*********************************
void C_TQ_Fe::get_data()
{
	tq_gdat(iwsg, iwse);
	if(tq_sg2err(&ierr) )
	{
		cout << "cannot get data from database \n";
	}
}
//*********************************
void C_TQ_Fe::get_elements_index()
{
	for(int i = 0; i < number_of_elements; i++)
	{
		tq_gsci(&element_index[i], table_of_elements[i], iwsg, iwse);
	}
}
//*********************************
void C_TQ_Fe::get_phases_number()
{
	tq_gnp(&number_of_phases, iwsg, iwse);
}
//*********************************
void C_TQ_Fe::set_cond_p_n()
{
	tq_setc("P", -1, -1, p, &num_condition_p, iwsg, iwse);
	tq_setc("N", -1, -1, n, &num_condition_n, iwsg, iwse);
}
//*********************************
//*** END OF CONSTRUCTORS PART ****
//*********************************
void C_TQ_Fe::enter_conditions_ele_wtpc(TC_FLOAT ele0_wtpc, TC_FLOAT ele1_wtpc,
	TC_FLOAT ele2_wtpc, TC_FLOAT ele3_wtpc, TC_FLOAT ele4_wtpc, TC_FLOAT ele5_wtpc, TC_FLOAT ele6_wtpc,
	TC_FLOAT ele7_wtpc, TC_FLOAT ele8_wtpc)
{
	switch(number_of_elements)
	{
	case 10:
		element_wtpc[8] = ele8_wtpc;
	case 9:
		element_wtpc[7] = ele7_wtpc;
	case 8:
		element_wtpc[6] = ele6_wtpc;
	case 7:
		element_wtpc[5] = ele5_wtpc;
	case 6:
		element_wtpc[4] = ele4_wtpc;
	case 5:
		element_wtpc[3] = ele3_wtpc;
	case 4:
		element_wtpc[2] = ele2_wtpc;
	case 3:
		element_wtpc[1] = ele1_wtpc;
	case 2:
		element_wtpc[0] = ele0_wtpc;
	default:
		break;
	}
	set_conditions_ele_wtpc();
}
//*********************************
void C_TQ_Fe::enter_conditions_ele_wtpc(TC_FLOAT *ele0_wtpc, TC_FLOAT *ele1_wtpc,
	TC_FLOAT *ele2_wtpc, TC_FLOAT *ele3_wtpc, TC_FLOAT *ele4_wtpc, TC_FLOAT *ele5_wtpc, TC_FLOAT *ele6_wtpc,
	TC_FLOAT *ele7_wtpc, TC_FLOAT *ele8_wtpc)
{
	switch(number_of_elements)
	{
	case 10:
		element_wtpc[8] = *ele8_wtpc;
	case 9:
		element_wtpc[7] = *ele7_wtpc;
	case 8:
		element_wtpc[6] = *ele6_wtpc;
	case 7:
		element_wtpc[5] = *ele5_wtpc;
	case 6:
		element_wtpc[4] = *ele4_wtpc;
	case 5:
		element_wtpc[3] = *ele3_wtpc;
	case 4:
		element_wtpc[2] = *ele2_wtpc;
	case 3:
		element_wtpc[1] = *ele1_wtpc;
	case 2:
		element_wtpc[0] = *ele0_wtpc;
	default:
		break;
	}
	set_conditions_ele_wtpc();
}
//*********************************
void C_TQ_Fe::set_conditions_ele_wtpc()
{
	for(int i = 0; i < number_of_elements - 1; i++)
	{
		tq_setc("W%", -1, element_index[i], element_wtpc[i], &num_condition_wtp[i], iwsg, iwse);
	}
}
//*********************************
void C_TQ_Fe::enter_condition_temperature(TC_FLOAT temp)
{
	temperature = temp;
	tq_setc("T", -1, -1, temperature, &num_condition_t, iwsg, iwse);
}
//*********************************
void C_TQ_Fe::enter_condition_temperature(TC_FLOAT *temp)
{
	temperature = *temp;
	tq_setc("T", -1, -1, temperature, &num_condition_t, iwsg, iwse);
}
//*********************************
void C_TQ_Fe::compute_equilibrium()
{
	tq_ce("", 0, 0, 0, iwsg, iwse);
}
//*********************************
//**** SET CONDITIONS AND C-E *****
//*********************************
void C_TQ_Fe::print_equilibrium()
{
	tq_sio("OUTPUT", 6);

	tq_le(iwsg, iwse);
	cout << "\n" << separator << "\n\n";
}
//*********************************
void C_TQ_Fe::get_phase_index(char *ph_name)
{
	tq_gpi(&phase_index, ph_name, iwsg, iwse);

	//table_of_phases = new char[number_of_phases][TC_STRLEN_PHASES];
	//phase_index = new TC_INT[number_of_phases];

	//for(int i = 0; i < number_of_phases; i++)
	//{
	//	tq_gpn(i + 1, table_of_phases[i], TC_STRLEN_PHASES, iwsg, iwse);
	//	
	//	tq_gpi(&phase_index[i], table_of_phases[i], iwsg, iwse);
	//	cout << table_of_phases[i] << endl;
	//}

	//for(int i = 0; i < number_of_phases; i++)
	//{
	//	cout << "Phase: " << table_of_phases[i] << ",\t" "has index: " << phase_index[i] << endl;
	//}
}
//*********************************
//void C_TQ_Fe::remove_condition_temperature()
//{
//	tq_remc(num_condition_t, iwsg, iwse);
//}