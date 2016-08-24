#ifndef __SECTION_BIRTH_H_
#define __SECTION_BIRTH_H_

Patch *new_patch(vector<Grain*> *habitat, Parameter *Initial_Parameter, double lambda, vector<Species*> com){

	uInt ngrain = habitat->size();
	double size = Initial_Parameter->get_size();
	int M = Initial_Parameter->get_M();

	uInt chosen_grain = get_rand_integer(ngrain);

	Grain *grain = (*habitat)[chosen_grain];

	double _lambda = lambda*sqrt(get_random());
	double _angle = 2*PI*get_random();
	double q = get_random_patch(grain->get_q(), Initial_Parameter->get_aggreg());

	double x_pos =  grain->get_x_pos();
	double x = x_pos+_lambda*cos(_angle);
		if(x<0.0)
			x+=size;
		if(x> size)
			x-=size;

	double y_pos =  grain->get_y_pos();
	double y = y_pos +_lambda*sin(_angle);
		if(y<0.0)
			y+=size;
		if(y> size)
			y-=size;

		Patch *new_p = new Patch(x, y, q, M, com, x_pos, y_pos);

		return new_p;
}

void update_birth(Patch *new_p, Grid *grid, PatchPerGrid *ppg,
	ColonizationPerGrid *cpg, Parameter *Initial_Parameter, vector<Species*> com,
	 Variable *var, double max_l){

	double x = new_p->get_x();
	double y = new_p->get_y();

	int ni = (int)floor (x/max_l);
	int nj = (int) floor (y/max_l);
	double temp_rate, sum = 0.0;

	for ( int m = 0; m < Initial_Parameter->get_M() ; m++){
		new_p->occupancy[m] = 0;
		temp_rate = total_colonization_rate(new_p, com[m], grid, Initial_Parameter, max_l);
		new_p->modify_col_patch(m, temp_rate);
		//new_p->colonization[m] = temp_rate;
		sum+=temp_rate;
	}

	grid->add_patch(new_p);
	var->modify_birth();// add one to temp_result 0 and muR to total death
	ppg->add_one_sp(ni, nj); //update the number of patch per grid
	cpg->modify_col_rate( ni, nj, sum); //update the colonization rate per grid

	var->modify_col(sum);

}

#endif //__SECTION_BIRTH_H_
