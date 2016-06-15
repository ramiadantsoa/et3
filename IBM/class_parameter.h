#ifndef __PARAMETER_H_
#define __PARAMETER_H_
class Parameter{
private:
	double size;
	int M;
	double l_s;
	double z;
	double T;
	short int establishment;
	short int resource_aggreg;
	double muR;
	double col_a;
	double diff_a;
	string name;
	int replicates;
	double lambda;
	double gammaH;
	
public:
	Parameter(double _size, short int _M, double _l_s, double _z, double _T, short int _establishment, short int _resource_aggregation, double _muR, int _repl,double _col_a,double _diff_a,double _lambda, double _gH){
		this->size =_size;
		this->M=_M;
		this->l_s=_l_s;
		this->z= _z;
		this->T = _T;
		this->establishment=_establishment;
		this->lambda=_lambda;
		this->gammaH = _gH;
		this->resource_aggreg= _resource_aggregation;
		this->muR = _muR;
		this->replicates = _repl;
		this->col_a=_col_a;
		this->diff_a = _diff_a;
	}

	double get_lambda(){
		return this->lambda;
	}

	double get_gammaH() {
		return this->gammaH;
	}

	int get_replicates(){
		return this->replicates;
	}

	double get_size(){
		return this->size;
	}

	short int get_M(){
		return this->M;
	}

	double get_l_s(){
		return this->l_s;
	}

	double get_z(){
		return this->z;
	}

	double get_T(){
		return this->T;
	}

	short int get_est(){
		return this->establishment;
	}

	short int get_aggeg(){
		return this->resource_aggreg;
	}

	double get_muR(){
		return this->muR;
	}

	double get_col_a(){
		return this->col_a;
	}

	double get_diff_a(){
		return this->diff_a;
	}

	string get_name(){
		ostringstream tostring;
		string est =establishment==1? "Est" : "Fec";
		int model = 3+ 4*(1- (int) this->z) + this->resource_aggreg; 
		tostring<<"task_1_"<<est<<"_model_"<<model<<"_M_"<<this->M<<"_size_"<<this->size<<"_T_"<<this->T<<"_ls_"<< this->l_s<<"_diff_a_"<<this->diff_a<<"_muR_"<<this->muR<<"_lambda_"<<this->lambda << "_gH_" << this->gammaH <<"_replicates_"<<this->replicates;
		name = tostring.str();
		return this->name;
	}

};

#endif // __PARAMETER_H_