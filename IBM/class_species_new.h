#ifndef __SPECIES_H_
#define __SPECIES_H_
class Species{
private:
	double col_a; //colonisation rate parameter
	double diff_a;
	double l_s; //dispersal length
	double opt_q;
	double nu_inv; //inverse of nu
	int sp_id;
	double log_bessel;

public:
	Species(double _col_a, double _diff_a, double _l_s, int _k, int _M, double _nu){
		this-> col_a = (_k==0)? _col_a*_diff_a : _col_a;
		this-> l_s = _l_s;
		this-> opt_q= (double) _k/(_M-1);
		this-> nu_inv= (_k==0)? 0 :(_M-1)*(_M-1)/_w; //it is the inverse of nu
		this-> sp_id = _k;
		this-> log_bessel = (_nu == -1)? -1: boost::math::cyl_bessel_i(order,1/nu);

	}

	double get_log_bessel(){
		return this->log_bessel;
	}
	int get_sp_id(){
		return this->sp_id;
	}

	double get_col_a(){
		return this->col_a;
	}

	/*THIS IS REALLY OLD WHEN THE DISPERSAL KERNEL HAS A TAIL
	double get_eps(){
		double pre_eps = 0.0;//gg(this->nu);
		return this->eps = pre_eps;
	}*/

	double get_l(){
		return this->l_s;
	}

	double get_opt_q(){
		return this->opt_q;
	}
};

#endif // __SPECIES_H_
