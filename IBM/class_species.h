#ifndef __SPECIES_H_
#define __SPECIES_H_
class Species{
private:
	double col_a; //colonisation rate parameter
	double opt_q;
	double nu;
	double l_s; //dispersal length

public:
	Species(double _col_a, double _opt_q, double _nu, double _l_s){
		this-> col_a = _col_a;
		this-> l_s = _l_s;
		this-> opt_q = _opt_q;
		this-> nu = _nu;
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
};

#endif // __SPECIES_H_
