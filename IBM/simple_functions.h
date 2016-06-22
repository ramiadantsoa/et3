#ifndef __SIMPLE_FUNCTIONS_H_
#define __SIMPLE_FUNCTIONS_H_

double get_log_bessel(double nu){
  /* Compute the logarithm of modifed Bessel of first order*/
  assert(nu == -1 || nu > 0.0015);
  long double result;
  if (nu == - 1) {
    result = -1; // this is symbolical, in fact if nu = -1 (meaning infinity = generalist), the fitness is directly 1 in class patch
  }
  else{
    double temp = boost::math::cyl_bessel_i(0,1/nu);
    assert(temp >= 0);
    // std::cout<< " the bessel is " << temp << endl;
    result = log(temp);
    // boost::multiprecision::cpp_dec_float_100 temp = boost::math::cyl_bessel_i(0,1/nu);
    // assert(result >= 0);
    // std::cout<< " the new value is " << temp << endl;
    // double result = 0;//log(temp);
  }
  return result;
}

//double get_random() {
    /* Give a random number from the unit range (0,1) */
//    return dsfmt_genrand_open_open(&rng);
//}


double get_random(){
  boost::uniform_real<double> dist(0,1);   // /uniform Distribution
  boost::variate_generator< boost::mt19937&, boost::uniform_real<>
 > rand(rng,dist);    // Variate generator
 return rand();
}

// Generate a random number between -aggreg and aggreg
double get_random_patch(double grain_q, double aggreg){
  double patch_q;
  double y = grain_q + aggreg*(2*get_random() - 1); // map (0,1) to (-aggreg + grain_q, aggreg + grain_q))
  patch_q = fmod(y + 1, 1);
 return patch_q;
}

/* Generate a random number from Exponential distribution with parameter 'lambda' */
double get_exp(double lambda) {
    double U = get_random();
	double result = - log(U)/lambda;
	return result;
}

/* Generate a random integer between 0 and value */
uInt get_rand_integer( uInt value){
	double temp0 = value*get_random();
	double result = 0.0;

	while(temp0 > result){ // this is slow and needs a better implementation (using floor)
		result+=1;
	}

	return (unsigned int) result-1;
}

/* Generate a random integer from a Poisson distribution with parameter 'mean' */
uInt get_poisson(double mean){
	double p =1.0;
	uInt k =0;

	if(mean<50){
		double a= exp(-mean);
		while( p>a){
			k+=1;
			p*=get_random();
		}
		k-=1;
	}
	else{
		double u = get_random();
		double v = get_random();
		double temp = sqrt(-2*log(u))*cos(2*PI*v);
		k = (uInt)floor(mean+sqrt(mean)*temp);
	}
	return (uInt) k;
}

/* Calculate z^n_occ */
double competition (int present_sp, double z){
	double result ;

	if( present_sp == 0){
		result =1;
	}
	else {
		if( z == 0.0){
			result = 0.0;
		}
		else{
			result = pow(z,present_sp);
		}
	}
	return result;
}

/* Calculate distance between (0,0) and (x,y) and the amount of propagules reaching (x,y) */
double dispersal (double x, double y , double L, double size){
	double result ;

	double x_coord = (x*x <= (size-x)*(size -x)) ? x*x : (size-x)*(size -x);
	double y_coord = (y*y <= (size-y)*(size -y)) ? y*y : (size-y)*(size -y);

	if(sqrt(x_coord +y_coord) <= L ){
		result = 1/(PI*L*L); // + eps/(size*size); because eps=0
	}
	else {
		result =0.0;// eps /(size*size);for now since eps =0
	}
	return result;
}

/*double d_sp(double L, int k, double gen_spe_ratio){
	double result;

	if(k==0){
		result= gen_spe_ratio*L;
	}
	else{
		result = L;
	}
	return result;
}*/


/*void avoid_repl(double times){
	double random=0.0;
	for (int r = 0; r<(int)times; r++){
			random=get_random();
	}
}*/

#endif // __SIMPLE_FUNCTIONS_H_
