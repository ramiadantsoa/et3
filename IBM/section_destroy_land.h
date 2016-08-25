#ifndef __SECTION_DESTROY_LAND_H_
#define __SECTION_DESTROY_LAND_H_


void destroy_patches_tau(Grid *grid, double destruct_param){
	int ncell = grid->get_number_of_cells();

	assert(destruct_param > 0. && destruct_param <1);
	double rand;

	for(int i = 0; i<ncell; i++) {
		for(int j = 0; j<ncell; j++) {
			uInt np = grid->count_per_grid(i,j);
			for(uInt k = 0; k < np; k++) {
				rand = get_random();
				if(rand > destruct_param){
					grid->delete_patch(i,j, (np-1)-k);
				}
			}
		}
	}
	grid->export_all_resource_units(2);
}

void destroy_patches_lambda(Grid *grid, vector<Grain*> *habitat, double destruct_param, double lambda){
	int ncell = grid->get_number_of_cells();
	double size = grid->get_size();
	assert(destruct_param > 0. && destruct_param <1);

	double x, y, cx, cy, dx, dy;

	for(int i = 0; i<ncell; i++) {
		for(int j = 0; j<ncell; j++) {
			uInt np = grid->count_per_grid(i,j);
			for(uInt k = 0; k < np; k++) {
				Patch *p = grid->get_ru(i,j, np-1-k);
				// Patch *p = (*patches)[(np-1)-k];
				x = p->get_x();
				y = p->get_y();
				cx = p->get_cx();
				cy = p->get_cy();

				dx = min(abs(cx-x), abs(cx - (size + x)));
				dy = min(abs(cy-y), abs(cy - (size + y)));
				if(sqrt(dx*dx + dy*dy) > lambda ){
					grid->delete_patch(i,j, np-1-k);
				}
			}
		}
	}
}

void destroy_patches_gammaH(Grid *grid, vector<Grain*> *habitat, double destruct_param){

	int ngrains = habitat->size();
	int ngrains_deleted = ceil(ngrains*(1-destruct_param));

	assert(ngrains != ngrains_deleted);

	int ncell = grid->get_number_of_cells();
	double size = grid->get_size();
	assert(destruct_param > 0. && destruct_param <1);

	double x_pos, y_pos, cx, cy;

	for(int ii = 0; ii < ngrains_deleted; ii++){
		Grain *g = (*habitat)[ii];
		x_pos = g->get_x_pos();
		y_pos = g->get_y_pos();


		for(int i = 0; i<ncell; i++) {
			for(int j = 0; j<ncell; j++) {
				uInt np = grid->count_per_grid(i,j);
				// vector<Patch*> *patches = grid -> get_cell(i,j);
				// Go through all patches in this cell
				for(unsigned int k = 0; k<np; k++) {
					Patch *p = grid->get_ru(i,j, np-1-k);

					cx = p->get_cx();
					cy = p->get_cy();

					if(cx == x_pos && cy == y_pos ){
						//  delete that patch
						grid->delete_patch(i,j, np-1-k);
					}
				}
			}
		}
	}
	//NEED TO DELETE THE GRAIN g
	habitat->erase(habitat->begin(), habitat->begin()+ngrains_deleted);
}
//
//
void update_general_variable(Grid *grid, Variable *var, ColonizationPerGrid *cpg, PatchPerGrid *ppg, vector<Species*> com, Parameter *Initial_Parameter, double destruct_param, double tau, double lambda, double gammaH, double size, double muR, double max_ls){
	uInt new_np = grid->count_patches();
	std::cout<<new_np;
	assert(new_np >= 1);

	double new_total_birth = PI*lambda*lambda*tau*gammaH*size*size;
	double new_total_death = new_np*muR;
	double new_total_colonization = Total_colonization(grid, com, Initial_Parameter, max_ls);


	var->modify_total_birth(new_total_birth);
	var->modify_total_death(new_total_death);
	var->modify_total_colonization(new_total_colonization);

	// reinitializing colonization per grid
	int ncell = grid->get_number_of_cells();
	double new_col_rate;
	for(int i = 0; i<ncell; i++) {
		for(int j = 0; j<ncell; j++){
			new_col_rate = grid->colrate_per_grid(i,j);
			cpg->assign_col_rate(i,j, new_col_rate);
			new_np = grid->count_per_grid(i,j);
			ppg->change_nsp(i,j, new_np);
		}
	}
}


#endif //__SECTION_DESTROY_LAND_H_
