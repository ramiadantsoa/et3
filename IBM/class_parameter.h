#ifndef __PARAMETER_H_
#define __PARAMETER_H_
class Parameter{
private:
	int com_id;
	int M;
	double z;
	short int establishment;
	double resource_aggreg;
	double muR;
	double lambda;
	double gammaH;
	double size;
	double T;
	int replicates;
	string name;


public:
	Parameter(int _com_id, int _M, double _z, short int _establishment,
		double _resource_aggregation, double _muR, double _lambda, double _gammaH,
		double _size, double _T, int _repl){
		this->com_id =_com_id;
		this->M =_M;
		this->z= _z;
		this->establishment=_establishment;
		this->resource_aggreg= _resource_aggregation;
		this->muR = _muR;
		this->lambda=_lambda;
		this->gammaH = _gammaH;
		this->size =_size;
		this->T = _T;
		this->replicates = _repl;
		this->com_id = _com_id;
}

	short int get_M(){
		return this->M;
	}

	double get_z(){
		return this->z;
	}

	short int get_est(){
		return this->establishment;
	}

	double get_aggreg(){
		return this->resource_aggreg;
	}

	double get_muR(){
		return this->muR;
	}

	double get_lambda(){
		return this->lambda;
	}

	double get_gammaH() {
		return this->gammaH;
	}

	double get_size(){
		return this->size;
	}

	double get_T(){
		return this->T;
	}

	int get_replicates(){
		return this->replicates;
	}

	string get_name(){
		ostringstream tostring;
		string est =establishment==1? "Est" : "Fec";
		tostring<<"sim_"<<est<<"_z_"<<this->z<<"_aggreg_"<<this->resource_aggreg
		<<"_M_"<<this->M<<"_com_id_"<<com_id<<"_size_"<<this->size<<"_T_"<<this->T<<"_muR_"<<this->muR
		<<"_lambda_"<<this->lambda << "_gH_" << this->gammaH <<"_replicates_"
		<<this->replicates;
		name = tostring.str();
		return this->name;
	}

};

#endif // __PARAMETER_H_
