#ifndef __SIMULATOR_H_
#define __SIMULATOR_H_

void simulation(Parameter *Initial_Parameter, vector<Species*> com){
	string valiny; //valiny stores the lambda, tau , gammaH, final result as a fraction of occupancy and the duration

	ostringstream convert;

	//INITIALIZE SPECIES PARAMETER
	int M = Initial_Parameter->get_M();
	double muR = Initial_Parameter->get_muR();
	double lambda=Initial_Parameter->get_lambda();
	double gammaH = Initial_Parameter->get_gammaH();
	double size =Initial_Parameter->get_size();
	double T = Initial_Parameter->get_T();

	//find species with highest dispersal range
	double new_max_ls, max_ls = 0.0;
	for (int i = 0; i < M; i++) {
		new_max_ls = std::max(max_ls, com[i]->get_l());
		max_ls = new_max_ls;
	}

	double MIN_TAU = 0.1;// 0.75; //this needs to be adjusted with pp and N_STEP_L
	double MAX_TAU = 10;
	//double pp = 1.15; //this is the power used to increment tau;

	int N_STEP_T = 100; // number of discretization of tau

	double tau;
	// SIMLUATE ALL DIFFERENT LEVEL OF FRAGMENTATION AND DENSITY OF RESOURCE PRODUCTION

		for(int j = 0; j < N_STEP_T ; j++){
			tau =(double) MAX_TAU*(j+1)/N_STEP_T;
			//tau = (double)  (m_MIN_TAU/(PI*lambda*lambda*gammaH))*pow(pp,j);
			// tau = (double) MIN_TAU*pow(pp,j);
			//lambda = (double) MAX_LAMBDA*(i+1)/N_STEP_L;
			convert<<tau<<"," <<lambda<<","<<gammaH;

			cout<<"\n \n lambda is " <<lambda<< ", tau is "<< tau<< ", gammaH is "<< gammaH<<endl;
            cout.flush();

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

			double t = 0.0;
      int half_T = (int) floor(Initial_Parameter->get_T()/2);
			double choose_event=0.0;
			double total_event=0.0;
			double dt=0;
            //double step = 0.0;
			vector<double> final_result(M+1);
			for (int m = 0 ; m < M+1; m++) final_result[m]= 0.0;

			if(npatches == 0){
				for( int m = 0; m < M ; m++){
					final_result[m+1]= 0.0;
					cout<<final_result[m+1]<<"\t";
					convert<<","<<final_result[m+1];
				}
				convert<<"\n";
				goto end;
			}

			while(t<T){
				total_event = var->get_birth()+ var->get_death() + var->get_col();//calculates the rate that an event happens
				dt = get_exp(total_event);
				t+=dt;
				choose_event= total_event*get_random();

				if(choose_event <= var->get_birth()) {
					//cout<<"birth \n";
					Patch *new_p = new_patch(habitat,Initial_Parameter, lambda, com);
					update_birth(new_p, grid, ppg, cpg,Initial_Parameter, com, var, max_ls);

					//double totcol = Total_colonization(grid, Initial_Parameter, SpChar);
					//cout<<" total colonization var " << var->get_col() << " and from calculation " << totcol<<endl;
				}else{
					if(choose_event <= var->get_birth()+var->get_death()) {
					//cout<<"death \n ";
					vector<uInt> chosen_death = select_patch_death(grid,ppg, var);
					update_death(chosen_death, grid, ppg, cpg, Initial_Parameter, com,var);
					npatches =grid->count_patches();
					if(npatches == 0){
						for( int m = 0; m < M ; m++){
							final_result[m+1]= 0.0;
							cout<<final_result[m+1]<<"\t";
							convert<<","<<final_result[m+1];
						}
						convert<<"\n";
						goto end;
					}

					//double totcol = Total_colonization(grid, Initial_Parameter, SpChar);
					//cout<<" total colonization var " << var->get_col() << " and from calculation " << totcol<<endl;
					}
					else {
					//cout<<"colonization \n " ;
					vector<uInt> chosen_col = select_patch_col(grid, cpg,var,M);
					update_colonization(chosen_col, grid, cpg, Initial_Parameter, var);
					update_colonization_rest(chosen_col, grid, cpg, Initial_Parameter, com, var);

					//double totcol = Total_colonization(grid, Initial_Parameter, SpChar);
					//cout<<" total colonization var " << var->get_col() << " and from calculation " << totcol<<endl;
					}
				}

				if(t > half_T + final_result[0] ){
					//if (t < half_T+step+1){
						//step+= 1.0;
						for (int m = 0 ; m<M+1 ; m++){
							assert(var->fraction_occupancy(m)<=1);
							final_result[m]+= var->fraction_occupancy(m);//(double) var->get_temp_result(m+1)/var->get_temp_result(0);
						}
					//}
				}

				//for(int m = 0; m < M+1; m++) cout<< final_result[m]<<endl;

			}//while closes

			/*start testing if the colonization update works

			double tt_col, death_rate;

			cout<< " death rate "<< var->get_death()<< "\n" << "computed death rate "<< var->get_temp_result(0)<<"\n";

			tt_col= total_colonization_check(grid,Initial_Parameter, SpChar);

			cout<<"total colonization "<< var->get_col()<<"\n"<<"long total colonization "<< tt_col<<"\n";


			system("pause");

			end testing*/

            cout << "Cleaning up\n";
            cout.flush();

			//cout<< final_result[0]<<"\n";
			output(grid);
			for( int m = 0; m < M ; m++){
				final_result[m+1]= final_result[m+1]/final_result[0];
				cout<<final_result[m+1]<<"\t";
				convert<<","<<final_result[m+1];
			}
			convert<<"\n";

			end:

			//grid->export_all_resource_units();

			delete grid;
            /* Delete everything in habitat */
            for(uInt i = 0; i < (*habitat).size(); i++) {
                delete (*habitat)[i];
            }
			delete habitat;
			delete cpg;
			delete ppg;
			delete var;
		} //for loop for tau closes
	//}

    /* Free the species */

    //for(int m=0 ; m < M ; m++){
    //    delete com[m];
    //}

	double time_passed;
	time_passed= clock();

	cout<<"\n sim time is : " << time_passed/CLOCKS_PER_SEC << "seconds "<<endl;

	convert<<"\n "<<time_passed/CLOCKS_PER_SEC<<endl;

	valiny.append(convert.str());

	ostringstream ltau;
	string addname;

	//ltau<<"_l_"<<lambda<<"_tau_"<<tau<<"_w_"<<w<<".csv";
	//  ltau<<".csv";
	// addname = ltau.str();
	string filename = Initial_Parameter->get_name().append(".csv");
	// string filename = Initial_Parameter->get_name().append(addname);

	ofstream myfile;
	myfile.open (filename.c_str());
	myfile << valiny;
	myfile.close();
}//simulator closes


#endif // __SIMULATOR_H_
