#ifndef __SECTION_COLONIZATION_H_
#define __SECTION_COLONIZATION_H_

vector<uInt> rand_choose_patch(Grid * grid, ColonizationPerGrid *cpg, int ni, int nj, int M){

	assert(0.0<=cpg->get_grid(ni,nj));

	double cpg_in_ni_nj = cpg->get_grid(ni,nj);

	double random = cpg_in_ni_nj*get_random();
	vector<uInt> result(2);
	double sum = 0.0;
	vector<Patch*> *patches= grid->get_cell(ni,nj);

	for( uInt k = 0 ; k < patches->size() ; k++){
		Patch *p = (*patches)[k];
		for( int m = 0 ; m < M ; m++){
			sum+= p->get_colonization(m);
			if (random <= sum){
				result[0]= m;
				result[1]= k ;goto end;
			}
		}
	}
end:
	return result;
}

vector<uInt> select_patch_col(Grid *grid, ColonizationPerGrid *cpg, Variable *var, int M){

	vector<uInt> result(4);

	vector<int> chosen_cell = cpg->rand_choose_gridC(var->get_col());// get_col gives the total colonization rate

	int ncell = grid->get_number_of_cells();
	int ni= chosen_cell[0];
	int nj= chosen_cell[1];

	assert(ni<ncell && nj <ncell);

	vector<uInt> chosen_patch_sp = rand_choose_patch(grid, cpg, ni, nj, M);

	int sp_id = chosen_patch_sp[0];//species identity
	assert(sp_id<M);
	uInt chosen_patch = chosen_patch_sp[1];// patch number in the grid
	vector<Patch*> *cell =grid->get_cell(ni,nj);

	/*if(cell->size()== 0){
		cout<<"grid size: "<< cell->size()<< "chosen patch: " << chosen_patch<<"\n";
		cout<<"ni: "<<ni<<", nj: "<<nj<<"\n";
		outputcpg(cpg);
	}
	*/
	assert(chosen_patch < cell->size());


	result[0]=(uInt) ni;
	result[1]=(uInt) nj;
	result[2]=chosen_patch;
	result[3]= sp_id;

	return result;
}

void update_colonization(vector<uInt> select_patch_col, Grid *grid, ColonizationPerGrid *cpg, Parameter *Initial_Parameter, Variable *var){
	assert(select_patch_col.size() == 4);
  int ni = select_patch_col[0];
	int nj = select_patch_col[1];
	uInt chosen_patch = select_patch_col[2];
	int sp_id = select_patch_col[3];

	vector<Patch*> *cell = grid->get_cell(ni,nj);

  assert((*cell).size() > chosen_patch);

	Patch *col_p = (*cell)[chosen_patch];
	col_p->occupancy[sp_id]=1;

	double last_col = col_p->get_colonization(sp_id);
	cpg->modify_col_rate(ni, nj, -last_col);
	var->modify_col(-last_col);

	col_p->set_col_to_0(sp_id);

	//updates to colonization rate for the colonized patch
	int M = Initial_Parameter->get_M();
	for( int m = 0 ; m < M ; m++){
		if(col_p->occupancy[m]==0){
			double old_local_col = col_p->get_colonization(m);
			double new_local_col = (Initial_Parameter->get_z()==0.0 )? 0.0 : Initial_Parameter->get_z()*old_local_col;

			//col_p->colonization[m]=new_local_col;

			double temp_local = new_local_col- old_local_col;
			col_p->modify_col_patch(m,temp_local);
			cpg-> modify_col_rate(ni, nj, temp_local );
			var->modify_col(temp_local); //- because modify_col substracts
		}
	}

	var->modify_temp_result(sp_id+1, 1);

	if (var->get_col()<0.0 || var->get_col()<chop_0){
		var->set_col_to_0();
	}
}

void update_colonization_rest(vector<uInt> select_patch_col, Grid *grid,
	ColonizationPerGrid *cpg, Parameter * Initial_Parameter,
	vector<Species*> com,Variable *var){
	int ni = select_patch_col[0];
	int nj = select_patch_col[1];
	uInt chosen_patch = select_patch_col[2];
	int sp_id = select_patch_col[3];

	vector<Patch*> *cell = grid->get_cell(ni,nj);

	Patch *col_p = (*cell)[chosen_patch];

	int ncell = grid->get_number_of_cells();
	int celli, cellj;
	double temp_rate;

	for( int i = -1 ; i < 2 ; i++) {
		for (int j = -1 ; j < 2 ; j++ ){
			celli = (ni + i + ncell)%ncell;
			cellj =(nj + j + ncell)%ncell;
			vector<Patch*> *temp_grid = grid->get_cell(celli, cellj);
			for (uInt k = 0; k < temp_grid->size(); k ++){
				Patch *target = (*temp_grid)[k];
				if(target->occupancy[sp_id]==0){
					temp_rate = colonization_rate(col_p, target , com[sp_id], Initial_Parameter);
					target->modify_col_patch(sp_id, temp_rate);
					//target->colonization[sp_id]+= temp_rate ;
					cpg->modify_col_rate(celli, cellj, temp_rate);
					var->modify_col(temp_rate);
				}
			}
		}
	}

	if (var->get_col()<0.0 || var->get_col()<chop_0){
		var->set_col_to_0();
	}

}


#endif //__SECTION_COLONIZATION_H_
