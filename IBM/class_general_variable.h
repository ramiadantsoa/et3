#ifndef __GENERAL_VARIABLE_H_
#define __GENERAL_VARIABLE_H_

class Variable{
private:
	double total_birth;
	double total_death;
	double total_colonization;
	uInt *temp_result;
	double muR;
	int M;

public:
	Variable(double _birth, double _death, double _col, int _M, double _muR){
		this->total_birth= _birth;
		this->total_death= _death;
		this->total_colonization= _col;
		this->muR= _muR;
		this->M = _M;
		this->temp_result = new uInt[_M+1];// this should contain number of patches, and the occupancy of each species??
		for( int m = 0 ; m < _M+1 ; m++){
			this->temp_result[m]= 0;
		}
	}

	~Variable(){
		delete[] this->temp_result;
	}

	double get_birth(){
		return this->total_birth;
	}

	double get_death(){
		return this->total_death;
	}

	double get_col(){
		return this->total_colonization;
	}

	uInt get_temp_result(int i){
		return this->temp_result[i];
	}

	void modify_birth(){
		this->total_death+=muR;
		this->temp_result[0]+=1;
	}

	void modify_death(){
		this->total_death-=muR;
		this->temp_result[0]-=1;
	}

	void modify_col(double temp_col){
		this->total_colonization+=temp_col;
	}

	void set_col_to_0(){
		this->total_colonization = 0.0;
	}

	void modify_temp_result(int i, uInt change){
		this->temp_result[i]+= change;
	}

	double fraction_occupancy(int i){
		return (double) this->temp_result[i]/this->temp_result[0];
	}

};


void initialize_temp_result(Variable *var, Grid *grid, int M){
	uInt temp_value;

	var->modify_temp_result(0,grid->count_patches());// the first entry is for the total number of resource units (called patches here)

	for ( int m = 0 ; m< M ; m++){
		temp_value = grid->count_species()[m];
		var->modify_temp_result(m+1,temp_value);
	}
}

void modify_temp_result_death(Patch * p, Variable *var, int M){
	short int temp_occup;
	for( int m = 0 ; m <M; m++){
		temp_occup = p->occupancy[m];
		var->modify_temp_result(m+1,-temp_occup);
	}
}

void output_temp_result(Variable *var, int M){
	for ( int m = 0 ; m<M+1 ;  m++){
		cout<<var->get_temp_result(m)<<" \t";
	}
	cout<<endl;
}


#endif //__GENERAL_VARIABLE_H_