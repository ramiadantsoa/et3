#ifndef __SPECIES_H_
#define __SPECIES_H_
class Species{
private:
	double col_a; //colonisation rate parameter
	double opt_q;
	double nu;
	double l_s; //dispersal length
	short int sp_id;

public:
	Species(double _col_a, double _opt_q, double _nu, double _l_s, short int _sp_id){
		this-> col_a = _col_a;
		this-> l_s = _l_s;
		this-> opt_q = _opt_q;
		this-> nu = _nu;
		this -> sp_id = _sp_id;
	}

	double get_sp_log_bessel(){
		return get_log_bessel(this->nu);
	}

	double get_col_a(){
		return this->col_a;
	}

	double get_l(){
		return this->l_s;
	}

	double get_opt_q(){
		return this->opt_q;
	}

	double get_nu(){
		return this->nu;
	}

	short int get_sp_id(){
		return this->sp_id;
	}
};

#endif // __SPECIES_H_
