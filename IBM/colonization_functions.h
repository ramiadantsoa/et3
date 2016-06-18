#ifndef __COLONIZATION_FUNCTION_H_
#define __COLONIZATION_FUNCTION_H_

double colonization_rate(Patch *source, Patch *target, Species *sp, Parameter *param){
	int present_sp = target->count_present_species();
	double col_a = sp->get_col_a();

	double s_x = source->get_x();
	double s_y = source->get_y();

	double t_x =target->get_x();
	double t_y =target->get_y();

	double size= param->get_size();
	double est=param->get_est();
	double z = param->get_z();

	double fitn = (est==1) ? target->get_fitness(sp->get_sp_id()) : source->get_fitness(sp->get_sp_id());

	double result = col_a*competition(present_sp,z)*fitn*dispersal((s_x-t_x),(s_y-t_y),sp->get_l(),size);

	return result;
}


double total_colonization_rate(Patch *target, Species *sp, Grid *g, Parameter *param, double max_l){

	int loc_i = (int)floor ((target->get_x())/ max_l);
	int loc_j = (int)floor ((target->get_y())/ max_l);

	int celli, cellj;
	int ncells = g->get_number_of_cells() ;
	int sp_id =sp->get_sp_id();

	double size = param->get_size();

	double result=0.0;

	for ( int i = -1 ; i < 2 ; i++ ){
		for ( int j = -1 ; j < 2 ; j++){
			celli = (loc_i + i+ ncells)%ncells;
			cellj =(loc_j + j + ncells)%ncells;

			vector<Patch*> *grid = g->get_cell(celli, cellj);

			for (uInt k = 0; k < grid->size(); k ++){
				Patch *source = (*grid)[k];
				if(source->occupancy[sp_id]==1){
				result+= colonization_rate(source, target, sp, param);
				}
			}
		}
	}
	return result;
}

double total_colonization_check(Grid* grid, Parameter *Initial_Parameter, vector<Species*> com, double max_l){
	int ncell = grid->get_number_of_cells();
	int M = Initial_Parameter->get_M();

	double temp0;
	double total_colonization=0.0;

	for (int i= 0 ; i<ncell ; i++){
		for ( int j = 0 ; j<ncell ; j++ ){
			vector<Patch*> *cell = grid->get_cell(i,j);
			for(unsigned int k = 0 ; k < cell->size() ; k++ ){
				Patch *target = (*cell)[k];
				for( int m = 0 ; m < M ; m++){
					if(target->occupancy[m] == 0){
						temp0 = total_colonization_rate(target, com[m], grid, Initial_Parameter, max_l);
						//target->colonization[m] = temp0;
						total_colonization+= temp0;
					}
				}
			}
		}
	}

	return total_colonization;
}

double Total_colonization(Grid *grid, vector<Species*> com, Parameter *param, double max_l){

	double M = param->get_M();
	int ncell = grid->get_number_of_cells();
	double temp0, total_colonization = 0.0;

	for (int i= 0 ; i<ncell ; i++){
		for ( int j = 0 ; j<ncell ; j++ ){
			vector<Patch*> *cell = grid->get_cell(i,j);
			for(unsigned int k = 0 ; k < cell->size() ; k++ ){
				Patch *target = (*cell)[k];
				for( int m = 0 ; m < M ; m++){
					if(target->occupancy[m] == 0){
						temp0 = total_colonization_rate(target, com[m], grid, param, max_l);
						target->modify_col_patch(m, temp0);
						//target->colonization[m] = temp0;
						total_colonization+= temp0;
					}
				}
			}
		}
	}
	return total_colonization;
}

#endif //__COLONIZATION_FUNCTION_H_
