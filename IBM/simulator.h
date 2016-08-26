#ifndef __SIMULATOR_H_
#define __SIMULATOR_H_

void simulation(Parameter *Initial_Parameter, vector<Species*> com, double destruct_param){
	string output_final_occupancy; //output_final_occupancy stores the lambda, tau , gammaH, final result as a fraction of occupancy and the duration

	ostringstream save_final_occupancy;

	string output_time_occupancy; // stores the changes in occupancy through time

	ostringstream save_time_occupancy;

	//INITIALIZE SPECIES PARAMETER
	int M = Initial_Parameter->get_M();
	double muR = Initial_Parameter->get_muR();
	double tau=Initial_Parameter->get_tau();
	double lambda=Initial_Parameter->get_lambda();
	double gammaH = Initial_Parameter->get_gammaH();
	double size =Initial_Parameter->get_size();
	double T = Initial_Parameter->get_T();

	//find species with longest dispersal range
	double new_max_ls, max_ls = 0.0;
	for (int i = 0; i < M; i++) {
		new_max_ls = std::max(max_ls, com[i]->get_l());
		max_ls = new_max_ls;
	}

	// SIMLUATE ALL DIFFERENT LEVEL OF FRAGMENTATION AND DENSITY OF RESOURCE PRODUCTION

		for(int j = 0; j < 1 ; j++){
			save_final_occupancy << tau << "," << lambda <<","<<gammaH;

			cout<<"\n \nlambda is " <<lambda<< ", tau is "<< tau<< ", gammaH is "<< gammaH<<endl;
            cout.flush();

			save_time_occupancy << tau <<"," << lambda << "," << gammaH;


			//INITIALIZE HABITAT GRAIN: PATCHES LOCATION AND TYPE
			vector<Grain*> *habitat = initialize_habitat(Initial_Parameter,gammaH, lambda);

			//INITIALIZE SPECIES INCIDENCE
			double resource_density = PI*lambda*lambda*tau*gammaH/muR;
			//double resource_density = PI*lambda*lambda*tau/ muR;

			vector<double> prevalence = initialize_incidence(Initial_Parameter, resource_density, com);

			//INITIALIZE GRID AND PATCH PER GRID
			Grid *grid = initialize(Initial_Parameter, resource_density, prevalence, habitat, com, max_ls);
			PatchPerGrid *ppg = initializePatchGrid(grid);

			//INITIALIZE COLONIZATION RATE AND COLONIZATION PER GRID
			double total_colonization = Total_colonization(grid, com, Initial_Parameter, max_ls);
			ColonizationPerGrid *cpg= initializeColonizationGrid(grid);

			//outputcpg(cpg);

			//INITIALIZE VARIABLE CONTAINING BIRTH DEATH AND COLONIZATION RATE
			double total_birth = PI*lambda*lambda*tau*gammaH*size*size;
			uInt npatches = grid->count_patches();
			double total_death = npatches*muR;

			Variable *var = new Variable(total_birth, total_death, total_colonization, M, muR);

			initialize_temp_result(var,grid, M);//this contains the number of patches and species occupancy, it is a private of var.

			cout<<"number of resource units "<<npatches<<" initialization done\n";

			double sum_occ = 0.0; //this is used to test if the species persist and is used very late (declare here to avoid cross initialization)
			double t = 0.0;
			double discrete_time = 0.0;
			double step_time = 50.0; // for transient time, occupancy is recorded for each  time unit.

      int half_T = (int) floor(Initial_Parameter->get_T()/2);
			double choose_event=0.0;
			double total_event=0.0;
			double dt=0;
            //double step = 0.0;
			vector<double> final_result(M+1);
			for (int m = 0 ; m < M+1; m++) final_result[m]= 0.0;

			uInt total_abundance = var->get_total_abundance();

			if(npatches == 0){
				save_time_occupancy << "{" << t;
				for( int m = 0; m < M ; m++){
					final_result[m+1]= 0.0;
					cout<<final_result[m+1]<<"\t";
					save_final_occupancy<<","<<final_result[m+1];
					save_time_occupancy << ";" << final_result[m+1];
				}
				save_final_occupancy<<"\n";
				save_time_occupancy <<"} \n";
				goto end;
			}

			// save_time_occupancy << ",{" << discrete_time;
			// for (int m = 1; m < M+1; m++) {
			// 	assert(var->fraction_occupancy(m)<=1. && var->fraction_occupancy(m)>=0.);
			// 	save_time_occupancy << ";" << var->fraction_occupancy(m);
			// }
			// save_time_occupancy <<"}";


			while(t<T && total_abundance != 0){
				total_event = var->get_birth()+ var->get_death() + var->get_col();//calculates the rate that an event happens
				dt = get_exp(total_event);
				t+=dt;
				choose_event= total_event*get_random();

				if(choose_event <= var->get_birth()) {
					// cout<<"birth \n";
					Patch *new_p = new_patch(habitat,Initial_Parameter, lambda, com);
					update_birth(new_p, grid, ppg, cpg,Initial_Parameter, com, var, max_ls);

					// double totcol_b = total_colonization_check(grid, com,  Initial_Parameter,  max_ls);
					// cout<<" birth total colonization var " << var->get_col() << " and from calculation " << totcol_b<<endl;
				}
				else{
					if(choose_event <= var->get_birth() + var->get_death()) {
					// cout<<"death \n ";
					vector<uInt> chosen_death = select_patch_death(grid,ppg, var);
					update_death(chosen_death, grid, ppg, cpg, Initial_Parameter, com,var);
					npatches =grid->count_patches();
					if(npatches == 0){
						save_time_occupancy << "{" << discrete_time;
						for( int m = 0; m < M ; m++){
							final_result[m+1]= 0.0;
							cout<<final_result[m+1]<<"\t";
							save_final_occupancy<<","<<final_result[m+1];
							save_time_occupancy << ";" << final_result[m+1];
						}
						save_final_occupancy<<"\n";
						save_time_occupancy<<"} \n";
						goto end;
					}

					// double totcol_d = total_colonization_check(grid, com,  Initial_Parameter, max_ls);
					// cout<<" death total colonization var " << var->get_col() << " and from calculation " << totcol_d<<endl;
					}
					else {
						// cout<<"colonization \n " ;
					vector<uInt> chosen_col = select_patch_col(grid, cpg,var,M);
					update_colonization(chosen_col, grid, cpg, Initial_Parameter, var);
					update_colonization_rest(chosen_col, grid, cpg, Initial_Parameter, com, var);

					// double totcol_c = total_colonization_check(grid, com , Initial_Parameter, max_ls);
					// cout<<" colonization total colonization var " << var->get_col() << " and from calculation " << totcol_c<<endl;
					}
				}

				// if (discrete_time < t) {
				// 	save_time_occupancy <<",{ " << discrete_time;
				// 	for (int m = 1; m < M+1;  m++) {
				// 		assert(var->fraction_occupancy(m)<=1. && var->fraction_occupancy(m)>=0.);
				// 		save_time_occupancy << ";" << var->fraction_occupancy(m);
				// 	}
				// 	save_time_occupancy << "}";
				// 	discrete_time+= step_time;
				// }
				//
				// total_abundance = var->get_total_abundance();
				//
				// if(t > half_T + final_result[0] || total_abundance == 0 ){
				// 	//if (t < half_T+step+1){
				// 		//step+= 1.0;
				// 		for (int m = 0 ; m<M+1 ; m++){
				// 			assert(var->fraction_occupancy(m)<=1.);
				// 			final_result[m]+= var->fraction_occupancy(m);//(double) var->get_temp_result(m+1)/var->get_temp_result(0);
				// 		}
				// 	//}
				// }

				//system("read");
				//for(int m = 0; m < M+1; m++) cout<< final_result[m]<<endl;

			}// first while closes

			// save_time_occupancy << "\n";

			// INITIATE HABITAT DESTRUCTION
			int case123;

			std::cout << "enter 1 if reducing tau, 2 if reducing lambda, 3 is reducing gammaH: ";
			std::cin >> case123;
			grid->export_all_resource_units(1);
			switch (case123) {
				case 1:            // reduce tau
					tau *= destruct_param;
					destroy_patches_tau(grid, destruct_param);
					grid->export_all_resource_units(2);
					break;
				case 2:            // reduce lambda
					lambda *= sqrt(destruct_param);
					destroy_patches_lambda(grid, habitat, destruct_param, lambda);
					grid->export_all_resource_units(2);
					break;
				case 3:            // perturb gammaH
					gammaH *= destruct_param;
					destroy_patches_gammaH(grid, habitat, destruct_param);
					grid->export_all_resource_units(2);
					export_all_patches(habitat, 2);
					break;
			}
			update_general_variable(grid, var, cpg, ppg, com, Initial_Parameter, destruct_param, tau, lambda, gammaH, size, muR, max_ls);

			discrete_time = ceil(T);
			while(t<2*T && total_abundance != 0){
				total_event = var->get_birth()+ var->get_death() + var->get_col();//calculates the rate that an event happens
				dt = get_exp(total_event);
				t+=dt;
				choose_event= total_event*get_random();

				if(choose_event <= var->get_birth()) {
					// cout<<"birth \n";
					Patch *new_p = new_patch(habitat,Initial_Parameter, lambda, com);
					update_birth(new_p, grid, ppg, cpg,Initial_Parameter, com, var, max_ls);

					// double totcol_b = total_colonization_check(grid, com,  Initial_Parameter,  max_ls);
					// cout<<" birth total colonization var " << var->get_col() << " and from calculation " << totcol_b<<endl;
				}
				else{
					if(choose_event <= var->get_birth() + var->get_death()) {
					// cout<<"death \n ";
					vector<uInt> chosen_death = select_patch_death(grid,ppg, var);
					update_death(chosen_death, grid, ppg, cpg, Initial_Parameter, com,var);
					npatches =grid->count_patches();
					if(npatches == 0){
						save_time_occupancy << "{" << discrete_time;
						for( int m = 0; m < M ; m++){
							final_result[m+1]= 0.0;
							cout<<final_result[m+1]<<"\t";
							save_final_occupancy<<","<<final_result[m+1];
							save_time_occupancy << ";" << final_result[m+1];
						}
						save_final_occupancy<<"\n";
						save_time_occupancy<<"} \n";
						goto end;
					}

					// double totcol_d = total_colonization_check(grid, com,  Initial_Parameter, max_ls);
					// cout<<" death total colonization var " << var->get_col() << " and from calculation " << totcol_d<<endl;
					}
					else {
						// cout<<"colonization \n " ;
					vector<uInt> chosen_col = select_patch_col(grid, cpg,var,M);
					update_colonization(chosen_col, grid, cpg, Initial_Parameter, var);
					update_colonization_rest(chosen_col, grid, cpg, Initial_Parameter, com, var);

					// double totcol_c = total_colonization_check(grid, com , Initial_Parameter, max_ls);
					// cout<<" colonization total colonization var " << var->get_col() << " and from calculation " << totcol_c<<endl;
					}
				}

				if (discrete_time < t) {
					save_time_occupancy <<",{ " << discrete_time;
					for (int m = 1; m < M+1;  m++) {
						assert(var->fraction_occupancy(m)<=1. && var->fraction_occupancy(m)>=0.);
						save_time_occupancy << ";" << var->fraction_occupancy(m);
					}
					save_time_occupancy << "}";
					discrete_time+= step_time;
				}

				total_abundance = var->get_total_abundance();

				if(t > T+ half_T + final_result[0] || total_abundance == 0 ){
					//if (t < half_T+step+1){
						//step+= 1.0;
						for (int m = 0 ; m<M+1 ; m++){
							assert(var->fraction_occupancy(m)<=1.);
							final_result[m]+= var->fraction_occupancy(m);//(double) var->get_temp_result(m+1)/var->get_temp_result(0);
						}
					//}
				}

				//system("read");
				//for(int m = 0; m < M+1; m++) cout<< final_result[m]<<endl;

			}// second while closes


      cout << "Cleaning up\n";
      cout.flush();

			output(grid);
			for( int m = 0; m < M ; m++){
				sum_occ += final_result[m+1];
				final_result[m+1]= final_result[m+1]/final_result[0];
				cout<<final_result[m+1]<<"\t";
				save_final_occupancy<<","<<final_result[m+1];
			}
			save_final_occupancy<<"\n";

			end:


			// grid->export_all_resource_units();

			delete grid;
            /* Delete everything in habitat */
            for(uInt i = 0; i < (*habitat).size(); i++) {
                delete (*habitat)[i];
            }
			delete habitat;
			delete cpg;
			delete ppg;
			delete var;


			/* CREATE A CONDITION SUCH THAT IF THE SPECIES PERSIST, THE SIMULATION ENDS*/
			// if(sum_occ > 0.) break;
		} //for loop for destruction parameter (or landscape quality) closes
	//}


	double time_passed;
	time_passed= clock();

	cout<<"\n sim time is : " << time_passed/CLOCKS_PER_SEC << "seconds "<<endl;

	save_final_occupancy<<"\n "<<time_passed/CLOCKS_PER_SEC<<endl;

	output_final_occupancy.append(save_final_occupancy.str());
	output_time_occupancy.append(save_time_occupancy.str());

	// ostringstream ltau;
	// string addname;

	//ltau<<"_l_"<<lambda<<"_tau_"<<tau<<"_w_"<<w<<".csv";
	//  ltau<<".csv";
	// addname = ltau.str();
	string filename = Initial_Parameter->get_name().append(".csv");
	string filename2 = Initial_Parameter->get_name().append(".csv").insert(0,"time_occ_");
	// string filename = Initial_Parameter->get_name().append(addname);

	ofstream final_occ;
	final_occ.open (filename.c_str());
	final_occ << output_final_occupancy;
	final_occ.close();


	ofstream time_occ;
	time_occ.open(filename2.c_str());
	time_occ << output_time_occupancy;
	time_occ.close();
}//simulator closes

/*start testing if the colonization update works

double tt_col, death_rate;

cout<< " death rate "<< var->get_death()<< "\n" << "computed death rate "<< var->get_temp_result(0)<<"\n";

tt_col= total_colonization_check(grid,Initial_Parameter, SpChar);

cout<<"total colonization "<< var->get_col()<<"\n"<<"long total colonization "<< tt_col<<"\n";


system("pause");

end testing*/


#endif // __SIMULATOR_H_
