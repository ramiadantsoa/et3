#ifndef __INITIALIZE_INCIDENCE_H_
#define __INITIALIZE_INCIDENCE_H_

vector<double> initialize_incidence(Parameter *Initial_Parameter , double resource_density, vector<Species*> com){
	double initial_density=0.0;
	int M = Initial_Parameter->get_M();
	double muR = Initial_Parameter->get_muR();
	double col_a;

	vector<double> prevalence(M);

	for(int i = 0; i<M; i++){
		col_a = com[i]->get_col_a();
		if(Initial_Parameter->get_est()==0){
			if(Initial_Parameter->get_z()==1){
				initial_density= 1-(muR*muR/(resource_density*col_a));
			}
			else{
				initial_density =(1-(muR*muR/(resource_density*col_a)))/M;
			}
		}
		else{
			if(Initial_Parameter->get_z()==1){
				initial_density = 0.5;
			}
			else{
				//initial_density= 0.2;
				initial_density =(1-(muR*muR/(resource_density*col_a)))/M;
			}
		}
		prevalence[i]= initial_density;
	}
	return prevalence;
}


#endif //INITIALIZE_INCIDENCE_H_
