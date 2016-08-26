#ifndef __SECTION_DEATH_H_
#define __SECTION_DEATH_H_

vector<uInt> select_patch_death(Grid *grid, PatchPerGrid *ppg, Variable *var){

	vector<uInt> result(3);

	vector<int> chosen_grid = ppg->rand_choose_gridD(var->get_temp_result(0));// gives the total number of patches

	int ncell = grid->get_number_of_cells();
	int ni= chosen_grid[0];//(int) chosen_grid % ncell;
	int nj= chosen_grid[1];//(int)(chosen_grid -ni)/ncell;

	assert(ni<ncell || nj <ncell);

	vector<Patch*> *cell = grid->get_cell(ni,nj);

	assert(ppg->get_grid(ni,nj) > 0);

	uInt chosen_patch = get_rand_integer(ppg->get_grid(ni,nj));//select the patch in grid ni nj
	if(chosen_patch>= cell->size()){
		cout<<"size_"<<cell->size()<<"_chosen_"<<chosen_patch<<"\n";
	}
		assert(chosen_patch<cell->size());

	result[0]=ni;
	result[1]=nj;
	result[2]=chosen_patch;

	return result;
}

void update_death(vector<uInt> select_patch_death, Grid *grid,
	PatchPerGrid *ppg, ColonizationPerGrid *cpg, Parameter *Initial_Parameter,
	vector<Species*> com, Variable *var){

	assert(select_patch_death.size() == 3);
	int ni = select_patch_death[0];
	int nj = select_patch_death[1];
	uInt chosen_patch = select_patch_death[2];

	vector<Patch*> *cell = grid->get_cell(ni, nj);

	Patch *old_p = (*cell)[chosen_patch];

	int M = Initial_Parameter->get_M();

	double temp_rate;

	for(int m = 0;  m< M;  m++){
		double old_col = old_p->get_colonization(m);

		var->modify_col(-old_col);
		cpg->modify_col_rate(ni, nj ,-old_col);
		var->modify_temp_result(m+1, -(old_p->occupancy[m]) );
		//var->modify_col(-old_p->colonization[m]);
		//cpg->modify_col_rate(ni, nj ,-old_p->colonization[m]);
	}

	var->modify_death(); //this updates the death rate by minus muR and # patches by minus one
	ppg->substract_one_nsp(ni,nj);

	int n_present = old_p->count_present_species();

	Patch *temp_p = new Patch(old_p->get_x(),old_p->get_y(),  old_p->get_quality(), M, com,  old_p->get_cx(),old_p->get_cy());

	for( int m = 0; m< M ; m++){
		temp_p->occupancy[m] = old_p->occupancy[m];
	}

	grid->delete_patch(ni, nj, chosen_patch);

	int ncell = grid->get_number_of_cells();
	int celli, cellj;

	if ( n_present!=0 ){
		for( int i = -1 ; i < 2 ; i++) {
			for (int j = -1 ; j < 2 ; j++ ){
				celli = (ni + i + ncell)%ncell;
				cellj =(nj + j + ncell)%ncell;
				vector<Patch*> *temp_grid = grid->get_cell(celli, cellj);

				for (uInt k = 0; k < temp_grid->size(); k ++){
					Patch *target = (*temp_grid)[k];
					for ( int m = 0 ; m < M; m++ ){
						if(temp_p->occupancy[m]==1 && target->occupancy[m]==0){
							temp_rate = colonization_rate(temp_p, target, com[m], Initial_Parameter);
							target->modify_col_patch(m, -temp_rate);
							//target->colonization[m]-= temp_rate ;
							cpg->modify_col_rate(celli, cellj, -temp_rate);
							var->modify_col(-temp_rate);
						}
					}
				}
			}
		}
	}

	delete temp_p;

	if(var->get_col()<0.0 || var->get_col()< chop_0){
		var->set_col_to_0();
	}
}


/*Patch *old_patch(Grid *grid, PatchPerGrid *ppg, Variable *var){

	int chosen_grid = ppg->rand_choose_gridD(var->get_temp_result(0));// gives the total number of patches

	int ncell = grid->get_number_of_cells();
	int ni= (int) chosen_grid%ncell;
	int nj= (int)(chosen_grid - ni)/ncell;

	vector<Patch*> *cell = grid->get_cell(ni,nj);

	uInt chosen_patch = get_rand_integer(ppg->get_grid(ni,nj));

	Patch * old_p = (*cell)[chosen_patch];

	return old_p;
}


double update_death(Patch *old_p, Grid *grid, ColonizationPerGrid *cpg, Parameter *Initial_Parameter, vector<Species*> SpChar, int ni , int nj){

	double temp_rate, total_colonization = 0.0;

	int M = Initial_Parameter->get_M();
	int celli, cellj;
	int ncell = cpg->n_cell();

	for( int i = -1 ; i < 2 ; i++) {
		for (int j = -1 ; j < 2 ; j++ ){
			celli = (ni + i + ncell)%ncell;
			cellj =(nj + j + ncell)%ncell;
			vector<Patch*> * temp_grid = grid->get_cell(celli, cellj);
			for (uInt k = 0; k < temp_grid->size(); k ++){
				Patch *target = (*temp_grid)[k];
				for ( int m = 0 ; m < M; m++ ){
					if(old_p->occupancy[m]==1 && target->occupancy[m]==0){
						temp_rate = colonization_rate(old_p, target , SpChar[m], Initial_Parameter);
						target->colonization[m]-= temp_rate ;
						cpg->modify_col_rate(celli, cellj, -temp_rate);
						total_colonization-=temp_rate;
					}
				}
			}
		}
	}

	return total_colonization;
}*/

#endif //__SECTION_DEATH_H_
