#include "C_TQ_XX.h"
//*********************************
string C_TQ_XX::separator(50, '-');
//char* C_TQ_XX::tc_installation_directory =  "c:\\Program Files\\Thermo-Calc\\2018a\\";
char* C_TQ_XX::tc_installation_directory = "c:\\Program Files\\Thermo-Calc\\2019a\\";
//char* C_TQ_XX::tc_installation_directory = "c:\\Program Files\\Thermo - Calc\\2019a\\";
// "c:\\Program Files\\Thermo-Calc\\4.1\\";
/*C:\\Program Files(x86)\\Thermo - Calc\\Thermo - Calc\\tc\\vers*/
char* C_TQ_XX::log_file_directory = "C:\\temp\\!tq";
// "C:\\Users\\WOMDVIIa\\AppData\\Local\\Temp\\!tq";
TC_STRING C_TQ_XX::used_database = "TCFE7";
//TC_STRING C_TQ_XX::used_database = "TCAL4";
//TC_STRING C_TQ_XX::used_database = "TCBIN";
TC_FLOAT C_TQ_XX::n = 1;
TC_FLOAT C_TQ_XX::p = 101325;
//*********************************
C_TQ_XX::C_TQ_XX(char* ele0, char* ele1, char* ele2, char* ele3, char* ele4,
	char* ele5, char* ele6, char* ele7, char* ele8)
{
	initialization();
	db_open();

	if (!ele1)
	{
		number_of_elements = 2;
		table_of_elements = new char* [number_of_elements];
		table_of_elements[0] = ele0;
		table_of_elements[1] = main_element;
	}
	else if (!ele2)
	{
		number_of_elements = 3;
		table_of_elements = new char* [number_of_elements];
		table_of_elements[0] = ele0;
		table_of_elements[1] = ele1;
		table_of_elements[2] = main_element;
	}
	else if (!ele3)
	{
		number_of_elements = 4;
		table_of_elements = new char* [number_of_elements];
		table_of_elements[0] = ele0;
		table_of_elements[1] = ele1;
		table_of_elements[2] = ele2;
		table_of_elements[3] = main_element;
	}
	else if (!ele4)
	{
		number_of_elements = 5;
		table_of_elements = new char* [number_of_elements];
		table_of_elements[0] = ele0;
		table_of_elements[1] = ele1;
		table_of_elements[2] = ele2;
		table_of_elements[3] = ele3;
		table_of_elements[4] = main_element;
	}
	else if (!ele5)
	{
		number_of_elements = 6;
		table_of_elements = new char* [number_of_elements];
		table_of_elements[0] = ele0;
		table_of_elements[1] = ele1;
		table_of_elements[2] = ele2;
		table_of_elements[3] = ele3;
		table_of_elements[4] = ele4;
		table_of_elements[5] = main_element;
	}
	else if (!ele6)
	{
		number_of_elements = 7;
		table_of_elements = new char* [number_of_elements];
		table_of_elements[0] = ele0;
		table_of_elements[1] = ele1;
		table_of_elements[2] = ele2;
		table_of_elements[3] = ele3;
		table_of_elements[4] = ele4;
		table_of_elements[5] = ele5;
		table_of_elements[6] = main_element;
	}
	else if (!ele7)
	{
		number_of_elements = 8;
		table_of_elements = new char* [number_of_elements];
		table_of_elements[0] = ele0;
		table_of_elements[1] = ele1;
		table_of_elements[2] = ele2;
		table_of_elements[3] = ele3;
		table_of_elements[4] = ele4;
		table_of_elements[5] = ele5;
		table_of_elements[6] = ele6;
		table_of_elements[7] = main_element;
	}
	else if (!ele8)
	{
		number_of_elements = 9;
		table_of_elements = new char* [number_of_elements];
		table_of_elements[0] = ele0;
		table_of_elements[1] = ele1;
		table_of_elements[2] = ele2;
		table_of_elements[3] = ele3;
		table_of_elements[4] = ele4;
		table_of_elements[5] = ele5;
		table_of_elements[6] = ele6;
		table_of_elements[7] = ele7;
		table_of_elements[8] = main_element;
	}
	else
	{
		number_of_elements = 10;
		table_of_elements = new char* [number_of_elements];
		table_of_elements[0] = ele0;
		table_of_elements[1] = ele1;
		table_of_elements[2] = ele2;
		table_of_elements[3] = ele3;
		table_of_elements[4] = ele4;
		table_of_elements[5] = ele5;
		table_of_elements[6] = ele6;
		table_of_elements[7] = ele7;
		table_of_elements[8] = ele8;
		table_of_elements[9] = main_element;
	}

	common_constructors_part();
}
//*********************************
//C_TQ_XX::C_TQ_XX(int n_ele, char *elements[]) : number_of_elements(n_ele + 1)
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
//	table_of_elements[number_of_elements - 1] = "XX";
//	//*********************************
//	common_constructors_part();
//}
//*********************************
C_TQ_XX::~C_TQ_XX()
{
	delete[] num_condition_mass;
	delete[] element_index;
	delete[] element_mass;
	delete[] table_of_elements;
	delete[] iwse;
	delete[] iwsg;
}
//*********************************
void C_TQ_XX::database_list()
{
	//tc_databases_strings *qwe;
	//TC_INT *temp;
	list_ele = new tc_elements_strings[200];
	number_ele = new TC_INT[100];

	tq_lide(list_ele, number_ele, iwsg, iwse);

	for (int i = 0; i < 200; i++)
	{
		cout << list_ele->element[i];
	}
	cout << endl << *number_ele << endl;

	delete[] list_ele;
	delete[] number_ele;
}
//*********************************
void C_TQ_XX::initialization()
{
	iwsg = new TC_INT[TC_NWSG];		//TQ workspace
	iwse = new TC_INT[TC_NWSE];		//TQ workspace
	//*********************************
	tq_ini3(tc_installation_directory, log_file_directory, TC_NWSG, TC_NWSE, iwsg, iwse);
	// tu by³a próba snl
	if (tq_sg2err(&ierr))
	{
		cout << "error initializing \n";
		system("pause");
	}
}
//*********************************
void C_TQ_XX::db_open()
{
	tq_opdb(used_database, iwsg, iwse);

	//C_TQ_XX::database_list();
}
//*********************************
void C_TQ_XX::common_constructors_part()
{
	input_elements_to_the_system();
	set_system_phases();
	get_data();
	//*********************************
	//Input of weight percent and other conditions
	element_mass = new TC_FLOAT[number_of_elements];								// zmiana kodu przy przekszta³ceniu warunków z N = 1, W% na masy wszystkich (³¹cznie z XX) pierwiastków
	element_index = new TC_INT[number_of_elements];
	num_condition_mass = new TC_INT[number_of_elements];							// zmiana kodu przy przekszta³ceniu warunków z N = 1, W% na masy wszystkich (³¹cznie z XX) pierwiastków

	get_elements_index();
	get_phases_number();
	//tq_snl(500, 1e-9, 1e-30, 'N', iwsg, iwse);
	set_cond_p();
}
//*********************************
void C_TQ_XX::input_elements_to_the_system()
{
	for (int i = 0; i < number_of_elements; i++)
	{
		tq_defel(table_of_elements[i], iwsg, iwse);
		if (tq_sg2err(&ierr))
		{
			cout << "There is no such element! Try again. \n";
			tq_reserr();
		}
	}
}
//*********************************
void C_TQ_XX::set_system_phases()
{
	tq_rejph("*", iwsg, iwse);
	//tq_resph("FCC LIQ", iwsg, iwse);	// tu jest miejsce, gdzie definiuje siê jakie fazy s¹ uwzglêdniane w obliczeniach w tej klasie
	//tq_resph("LIQ FCC", iwsg, iwse);	// tu jest miejsce, gdzie definiuje siê jakie fazy s¹ uwzglêdniane w obliczeniach w tej klasie
	//tq_resph("FCC_L12 LIQ", iwsg, iwse);	// tu jest miejsce, gdzie definiuje siê jakie fazy s¹ uwzglêdniane w obliczeniach w tej klasie
	tq_resph("BCC LIQ", iwsg, iwse);
}
//*********************************
void C_TQ_XX::get_data()
{
	tq_gdat(iwsg, iwse);
	if (tq_sg2err(&ierr))
	{
		cout << "cannot get data from database \n";
	}
}
//*********************************
void C_TQ_XX::get_elements_index()
{
	for (int i = 0; i < number_of_elements; i++)
	{
		tq_gsci(&element_index[i], table_of_elements[i], iwsg, iwse);
	}
}
//*********************************
void C_TQ_XX::get_phases_number()
{
	tq_gnp(&number_of_phases, iwsg, iwse);
}
//*********************************
void C_TQ_XX::set_cond_p()
{
	tq_setc("P", -1, -1, p, &num_condition_p, iwsg, iwse);
	//tq_setc("N", -1, -1, n, &num_condition_n, iwsg, iwse);						// zmiana kodu przy przekszta³ceniu warunków z N = 1, W% na masy wszystkich (³¹cznie z XX) pierwiastków
}
//*********************************
//*** END OF CONSTRUCTORS PART ****
//*********************************
void C_TQ_XX::enter_conditions_ele_mass(TC_FLOAT ele0_mass, TC_FLOAT ele1_mass,
	TC_FLOAT ele2_mass, TC_FLOAT ele3_mass, TC_FLOAT ele4_mass, TC_FLOAT ele5_mass, TC_FLOAT ele6_mass,
	TC_FLOAT ele7_mass, TC_FLOAT ele8_mass)
{
	switch (number_of_elements)													// zmiana kodu przy przekszta³ceniu warunków z N = 1, W% na masy wszystkich (³¹cznie z XX) pierwiastków
	{
	case 9:
		element_mass[8] = ele8_mass;
	case 8:
		element_mass[7] = ele7_mass;
	case 7:
		element_mass[6] = ele6_mass;
	case 6:
		element_mass[5] = ele5_mass;
	case 5:
		element_mass[4] = ele4_mass;
	case 4:
		element_mass[3] = ele3_mass;
	case 3:
		element_mass[2] = ele2_mass;
	case 2:
		element_mass[1] = ele1_mass;
	case 1:
		element_mass[0] = ele0_mass;
	default:
		break;
	}
	set_conditions_ele_mass();
}
//*********************************
void C_TQ_XX::enter_conditions_ele_mass(TC_FLOAT* ele0_mass, TC_FLOAT* ele1_mass,
	TC_FLOAT* ele2_mass, TC_FLOAT* ele3_mass, TC_FLOAT* ele4_mass, TC_FLOAT* ele5_mass, TC_FLOAT* ele6_mass,
	TC_FLOAT* ele7_mass, TC_FLOAT* ele8_mass)
{
	switch (number_of_elements)													// zmiana kodu przy przekszta³ceniu warunków z N = 1, W% na masy wszystkich (³¹cznie z XX) pierwiastków
	{
	case 9:
		element_mass[8] = *ele8_mass;
	case 8:
		element_mass[7] = *ele7_mass;
	case 7:
		element_mass[6] = *ele6_mass;
	case 6:
		element_mass[5] = *ele5_mass;
	case 5:
		element_mass[4] = *ele4_mass;
	case 4:
		element_mass[3] = *ele3_mass;
	case 3:
		element_mass[2] = *ele2_mass;
	case 2:
		element_mass[1] = *ele1_mass;
	case 1:
		element_mass[0] = *ele0_mass;
	default:
		break;
	}
	set_conditions_ele_mass();
}
//*********************************
void C_TQ_XX::enter_conditions_ele_mass(TC_FLOAT* tab_ele_mass)
{
	for (int i = 0; i < number_of_elements; i++)
	{
		element_mass[i] = /*1e20 * */tab_ele_mass[i];
	}
	set_conditions_ele_mass();
}
//*********************************
void C_TQ_XX::set_conditions_ele_mass()
{
	for (int i = 0; i < number_of_elements; i++)									// zmiana kodu przy przekszta³ceniu warunków z N = 1, W% na masy wszystkich (³¹cznie z XX) pierwiastków
	{
		tq_setc("M", -1, element_index[i], element_mass[i], &num_condition_mass[i], iwsg, iwse);	// zmiana kodu przy przekszta³ceniu warunków z N = 1, W% na masy wszystkich (³¹cznie z XX) pierwiastków
	}
}
//*********************************
void C_TQ_XX::enter_condition_temperature(TC_FLOAT temp)
{
	temperature = temp;
	tq_setc("T", -1, -1, temperature, &num_condition_t, iwsg, iwse);
}
//*********************************
void C_TQ_XX::enter_condition_temperature(TC_FLOAT* temp)
{
	temperature = *temp;
	tq_setc("T", -1, -1, temperature, &num_condition_t, iwsg, iwse);
}
//*********************************
void C_TQ_XX::compute_equilibrium()
{
	//for (int i = 0; i < number_of_elements; i++)
	//{
	//	cout << "element [" << i << "].mass = " << element_mass[i] << endl;
	//}
	tq_reserr();
	tq_ce("", 0, 0, 0, iwsg, iwse);
	if (tq_sg2err(&ierr))
	{
		tq_reserr();
		tq_sio("OUTPUT", 6);
		tq_ls(iwsg, iwse);
		cout << "c-e failed \n";
		for (int i = 0; i < number_of_elements; i++)
		{
			cout << "element [" << i << "].mass = " << element_mass[i] << endl;
		}
	}
}
//*********************************
TC_INT C_TQ_XX::compute_equilibrium(TC_INT a)
{
	a = 0;

	tq_reserr();
	tq_ce("", 0, 0, 0, iwsg, iwse);
	if (tq_sg2err(&ierr))
	{
		a = ierr;
		//tq_reserr();
		tq_sio("OUTPUT", 6);
		tq_ls(iwsg, iwse);
		cout << "c-e failed \n";
		for (int i = 0; i < number_of_elements; i++)
		{
			cout << "element [" << i << "].mass = " << element_mass[i] << endl;
		}
	}
	return a;
}
//*********************************
//**** SET CONDITIONS AND C-E *****
//*********************************
void C_TQ_XX::print_equilibrium()
{
	tq_sio("OUTPUT", 6);

	tq_lc(iwsg, iwse);
	cout << "\n" << separator << "\n\n";
	tq_le(iwsg, iwse);
	cout << "\n" << separator << "\n\n";
}
//*********************************
void C_TQ_XX::reinitiate()
{
	tq_pini(iwsg, iwse);
	set_cond_p();
}
//*********************************
void C_TQ_XX::get_phase_index(char* ph_name)
{
	tq_gpi(&phase_index, ph_name, iwsg, iwse);
}
//*********************************
//double& C_TQ_XX::get_value(TC_FLOAT temp, TC_INT no_of_phase, TC_FLOAT ele0_mass, TC_FLOAT ele1_mass,		// ------------------------------------------------- dodano do prog Jacka
//	TC_FLOAT ele2_mass, TC_FLOAT ele3_mass, TC_FLOAT ele4_mass, TC_FLOAT ele5_mass, TC_FLOAT ele6_mass,
//	TC_FLOAT ele7_mass, TC_FLOAT ele8_mass)
//{
//	enter_condition_temperature(temp);
//	enter_conditions_ele_mass(&ele0_mass, &ele1_mass, &ele2_mass, &ele3_mass,
//		&ele4_mass, &ele5_mass, &ele6_mass, &ele7_mass, &ele8_mass);
//	
//	compute_equilibrium();
//	print_equilibrium();
//
//	tq_get1("W", no_of_phase, 1, &returned_value, iwsg, iwse);
//	return returned_value;
//}