#include "C_TQ_Fe.h"
#include "C_TQ_single_phase.h"
#include "C_TQ_chem_pot.h"
#include "C_TQ_single_property_old.h"
#include "C_cell.h"
#include <fstream>
//*********************************
//const int		system_size_x		= 5,
//				no_of_fcc_cells		= 2;
//
//double			initial_concentration_C_fcc = 0.5,				// [wt %]
//				initial_concentration_C_fcc_max = 0.7,
//				initial_concentration_Si_fcc = 0.5,
//				initial_concentration_Si_fcc_max = 0.1,
//				initial_concentration_C_liq = 4.0,
//				initial_concentration_Si_liq = 1.5;
//
//double			initial_C_tab[system_size_x],
//				initial_Si_tab[system_size_x];
//
//double			flux_C_right[system_size_x - 1],
//				flux_Si_right[system_size_x - 1],
//				flux_Fe_right[system_size_x - 1];
//				//flux_C_right_recalculated[system_size_x - 1],
//				//flux_Si_right_recalculated[system_size_x - 1],
//				//flux_Fe_right_recalculated[system_size_x - 1];
////double			plus_, minus_, abs_sum_;
//
//const double	/*dx = 1e-2,*/						//	[m]  duplikat, istnieje w C_cell
//				dx_squared = C_cell::dx * C_cell::dx,
//				dt = 1e-4,							//	[s]
//				mobility_C_fcc = 5e-10,				//	[]
//				mobility_Si_fcc = 1e-10,			//	[]
//				mobility_Fe_fcc = 1e-14,			//	[]
//				mobility_C_liq = 1,				//	[]
//				mobility_Si_liq = 1,			//	[]
//				mobility_Fe_liq = 1,			//	[]
//				mobility_C_interface = 1,		//	[]
//				mobility_Si_interface = 1,		//	[]
//				mobility_Fe_interface = 1;		//	[]
//*********************************
//void boundary_conditions(C_cell *tablica)
//{
//	tablica[0].mass_C = tablica[1].mass_C;
//	tablica[0].mass_Si = tablica[1].mass_Si;
//	tablica[0].mass_Fe = tablica[1].mass_Fe;
//
//	tablica[system_size_x - 1].mass_C =	 tablica[system_size_x - 2].mass_C;
//	tablica[system_size_x - 1].mass_Si = tablica[system_size_x - 2].mass_Si;
//	tablica[system_size_x - 1].mass_Fe = tablica[system_size_x - 2].mass_Fe;
//}
//*********************************
int main()
{
	string	txt_filename;
	fstream conc_field_txt;

	for (int i = 1; i < system_size_x - 1; i++)
	{
		initial_C_tab[i] = initial_concentration_C_fcc + (i - 1) * (initial_concentration_C_fcc_max - initial_concentration_C_fcc) / (system_size_x - 3);
		initial_Si_tab[i] = initial_concentration_Si_fcc + (i - 1) * (initial_concentration_Si_fcc_max - initial_concentration_Si_fcc) / (system_size_x - 3);
	}

	C_cell test;

	test.set_initial_conditions(&initial_C_tab[1], &initial_Si_tab[1], "FCC");

	test.get_potentials();

	C_TQ_single_property_old test_dens("V", "FCC", "C", "Si");

	double temp_value = test_dens.get_value(1450, initial_C_tab[1], initial_Si_tab[1]);

	cout << "density = " << temp_value;


	C_cell wektor[system_size_x];
	cout << endl;

	//for (int i = 0; i < system_size_x; i++)
	//{
	//	wektor[i].present_yourself();
	//}
	//cout << "\n*****************************\n\n";

	//for (int i = 1; i < system_size_x - 1; i++)
	//{
	//	if (i < no_of_fcc_cells)
	//	{
	//		wektor[i].set_initial_conditions(&initial_concentration_C_fcc, &initial_concentration_Si_fcc, "FCC");
	//	}
	//	else
	//		wektor[i].set_initial_conditions(&initial_concentration_C_liq, &initial_concentration_Si_liq, "LIQ");
	//}
	
	for (int i = 0; i < system_size_x; i++)
	{
		wektor[i].set_initial_conditions(&initial_C_tab[i], &initial_Si_tab[i], "FCC");
	}
	boundary_conditions(wektor);

	//*************************************************

	for (int k = 0; k < 400; k++)
	{
		cout << "k = " << k << endl;

		for (int i = 0; i < system_size_x; i++)
		{
			wektor[i].system->reinitiate();
		}
		//*************************************************
		for (int j = 0; j < 5000; j++) // max 5000 - co tyle trzeba reinicjowaæ bo pamiêæ w TQ siê zapycha
		{
			for (int i = 0; i < system_size_x; i++)
			{
				wektor[i].get_potentials();
			}

			for (int i = 0; i < system_size_x - 1; i++)
			{
				flux_C_right[i]		= -mobility_C_fcc * (wektor[i + 1].chemical_potential_pointer[0] - wektor[i].chemical_potential_pointer[0]) / C_cell::dx;
				flux_Si_right[i]	= -mobility_Si_fcc * (wektor[i + 1].chemical_potential_pointer[1] - wektor[i].chemical_potential_pointer[1]) / C_cell::dx;
				flux_Fe_right[i]	= -mobility_Fe_fcc * (wektor[i + 1].chemical_potential_pointer[2] - wektor[i].chemical_potential_pointer[2]) / C_cell::dx;

				//plus_ = minus_ = abs_sum_ = 0.0;
				//if (flux_C_right[i] > 0)
				//{
				//	plus_ += flux_C_right[i];
				//}
				//else
				//	minus_ += abs(flux_C_right[i]);
				//if (flux_Si_right[i] > 0)
				//{
				//	plus_ += flux_Si_right[i];
				//}
				//else
				//	minus_ += abs(flux_Si_right[i]);
				//if (flux_Fe_right[i] > 0)
				//{
				//	plus_ += flux_Fe_right[i];
				//}
				//else
				//	minus_ += abs(flux_Fe_right[i]);
				//abs_sum_ = abs(flux_C_right[i]) + abs(flux_Si_right[i]) + abs(flux_Fe_right[i]);
				//if ((flux_C_right[i] + flux_Si_right[i] + flux_Fe_right[i]) > 0)
				//{
				//	if (flux_C_right[i] > 0)
				//	{
				//		flux_C_right_recalculated[i] = flux_C_right[i] * minus_ / (abs_sum_ - minus_);
				//	}
				//	else
				//		flux_C_right_recalculated[i] = flux_C_right[i];
				//}
				//else
				//{
				//	if (flux_C_right[i] < 0)
				//	{
				//		flux_C_right_recalculated[i] = flux_C_right[i] * plus_ / (abs_sum_ - plus_);
				//	}
				//	else
				//		flux_C_right_recalculated[i] = flux_C_right[i];
				//}
				////*************************************************
				//if ((flux_C_right[i] + flux_Si_right[i] + flux_Fe_right[i]) > 0)
				//{
				//	if (flux_Si_right[i] > 0)
				//	{
				//		flux_Si_right_recalculated[i] = flux_Si_right[i] * minus_ / (abs_sum_ - minus_);
				//	}
				//	else
				//		flux_Si_right_recalculated[i] = flux_Si_right[i];
				//}
				//else
				//{
				//	if (flux_Si_right[i] < 0)
				//	{
				//		flux_Si_right_recalculated[i] = flux_Si_right[i] * plus_ / (abs_sum_ - plus_);
				//	}
				//	else
				//		flux_Si_right_recalculated[i] = flux_Si_right[i];
				//}
				////*************************************************
				//if ((flux_C_right[i] + flux_Si_right[i] + flux_Fe_right[i]) > 0)
				//{
				//	if (flux_Fe_right[i] > 0)
				//	{
				//		flux_Fe_right_recalculated[i] = flux_Fe_right[i] * minus_ / (abs_sum_ - minus_);
				//	}
				//	else
				//		flux_Fe_right_recalculated[i] = flux_Fe_right[i];
				//}
				//else
				//{
				//	if (flux_Fe_right[i] < 0)
				//	{
				//		flux_Fe_right_recalculated[i] = flux_Fe_right[i] * plus_ / (abs_sum_ - plus_);
				//	}
				//	else
				//		flux_Fe_right_recalculated[i] = flux_Fe_right[i];
				//}
				////*************************************************
				//double sprawdzenie = flux_C_right_recalculated[i] + flux_Si_right_recalculated[i] + flux_Fe_right_recalculated[i];
				//if (abs(sprawdzenie) > 1e-10)
				//{
				//	cout << "sprawdzenie wykazalo sume przeliczonych strumieni rozna od zera\n"
				//		<< "sprawdzenie = " << sprawdzenie << endl;
				//}
				//cout << "flux_C_right = " << flux_C_right[i] << "\t\tflux_C_right_recalculated = " << flux_C_right_recalculated[i] << endl
				//	<< "flux_Si_right = " << flux_Si_right[i] << "\t\tflux_Si_right_recalculated = " << flux_Si_right_recalculated[i] << endl
				//	<< "flux_Fe_right = " << flux_Fe_right[i] << "\t\tflux_Fe_right_recalculated = " << flux_Fe_right_recalculated[i] << endl
				//	<< "sprawdzenie = " << sprawdzenie << endl;
				//cout << "\n*****************************\n\n";
			}

			for (int i = 1; i < system_size_x - 1; i++)
			{
				wektor[i].mass_C	= wektor[i].mass_C + (flux_C_right[i - 1] - flux_C_right[i]) * dt * dx_squared;
				wektor[i].mass_Si	= wektor[i].mass_Si + (flux_Si_right[i - 1] - flux_Si_right[i]) * dt * dx_squared;
				wektor[i].mass_Fe	= wektor[i].mass_Fe + (flux_Fe_right[i - 1] - flux_Fe_right[i]) * dt * dx_squared;
				
				//wektor[i].mass_C = wektor[i].mass_C + (flux_C_right_recalculated[i - 1] - flux_C_right_recalculated[i]) * dt * dx_squared;
				//wektor[i].mass_Si = wektor[i].mass_Si + (flux_Si_right_recalculated[i - 1] - flux_Si_right_recalculated[i]) * dt * dx_squared;
				//wektor[i].mass_Fe = wektor[i].mass_Fe + (flux_Fe_right_recalculated[i - 1] - flux_Fe_right_recalculated[i]) * dt * dx_squared;
			}

			boundary_conditions(wektor);
		}	// KONIEC PÊTLI WEWNÊTRZNEJ

		//*************************************************

		txt_filename = "conc_field_";
		txt_filename += to_string(k);
		txt_filename += ".txt";

		conc_field_txt.open(txt_filename, ios::out);

		conc_field_txt << "timestep\t" << k << endl;
		conc_field_txt << "i\t" << "conc_C\t" << "conc_Si\t" << "conc_Fe\n";
		
		for (int i = 0; i < system_size_x; i++)
		{
			conc_field_txt << i << "\t";
			for (int j = 0; j < 3; j++)
			{
				conc_field_txt << wektor[i].chemical_potential_pointer[j] << "\t";
			}
			conc_field_txt << endl;
		}

		conc_field_txt.close();

	}	// KONIEC PÊTLI ZEWNÊTRZNEJ

	//*************************************************

	system("pause");
}