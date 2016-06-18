#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <time.h>
#include <sstream>
#include <string>
#include <assert.h>

//#include "dSFMT/dSFMT.h"
#include <boost/math/special_functions/bessel.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#define PI 3.1415926535897932385

#define chop_0 0.000000000000001

using namespace std;

typedef unsigned long int uInt;


//const int rng_seed = (const int)time(NULL);
//dsfmt_t rng;

//void init_rng() {
//    dsfmt_init_gen_rand(&rng, rng_seed);
//}

boost::mt19937 rng(std::time(0));

#include "simple_functions.h"

/*CLASS DECLARATION AND DEFINITION*/
#include "class_parameter.h"
#include "class_species.h"
#include "class_grain.h"
#include "initialize_incidence.h"
#include "class_patch.h"
#include "class_grid.h"
#include "class_patch_per_grid.h"
#include "class_colonization_per_grid.h"
#include "class_general_variable.h"

/*COLONIZATION AND TOTAL COLONIZATION FUNCTION AND RANDOM CHOOOSE PATCH*/

#include "colonization_functions.h"
#include "section_birth.h"
#include "section_death.h"
#include "section_colonization.h"

/*  CONTAINS THE MAIN SIMULATOR */

#include "simulator.h"

int main(){
	double z =0.0, aggreg = 0.5, size = 4, simtime = 1, muR = 0.1;
	short int est = 1,  M = 5, replicates = 0, com_id = 2;
	double lambda = 0.2, gammaH = 1;

  std::time_t seed_time = std::time(0);
  rng.seed(seed_time);

  //std::cout << "random number is " << get_random()<< std::endl;

	// double log_bessel, nu;
  // std::cout << "enter nu ";
  // std::cin >> nu;
  // log_bessel = get_log_bessel(nu);
  // std::cout << "the log bessel of " << nu <<" is "<< log_bessel << std::endl;

	/*cout<<"Compet z: ";
	cin>>compet_z;
	cout<<"Aggregation: ";
	cin>>aggreg;
	cout<<"lambda: ";
	cin>>lambda;
	cout<<"gammaH: ";
	cin>>gammaH;
	cout<<"l_s: ";
	cin>>l_s;
	cout<<"diff_a:";
	cin>>diff_a;
	cout << "enter muR: ";
	cin>>muR;
  */
	//cout<<"\n";
	//cout<<"Community size: ";
	//cin>> M;
	M = 3;

/*	cout<<"size: ";
	cin>>size;
	cout<<"time: ";
	cin>>time;
*/
	short int rep = 0;

	// std::cout << "community size ";
	// std::cin >> M;
	// std::cout << "community ID ";
	// std::cin >> com_id;

	std::vector<double> pre_com(M*5);
	for (int i = 0; i < M; i++) {
		// std::cout << " \n for species" << i << std::endl;
		// std::cout << "enter col_a "<< std::endl;
		// std::cin >> com[4*i];

		pre_com[5*i] = 0.5; //colonization rate
		pre_com[5*i+1] = get_random(); //optimal q
		pre_com[5*i+2] = (i+1)*get_random(); // niche width nu, by convention, nu = -1 for generalist
		pre_com[5*i+3] = 0.2; // dispersal range
		pre_com[5*i+4] = i; // species id

		// std::cout << "enter optimal q " << std::endl;
		// std::cin >> com[4*i+1];
		//
		// std::cout << "enter nu " << std::endl;
		// std::cin >> com[4*i+2];
		//
		// std::cout << "enter dispersal l " << std::endl;
		// std::cin >> com[4*i+3];
		//
	}

	Parameter *param = new Parameter(com_id, M, z, est, aggreg, muR, lambda, gammaH,
		size, simtime, replicates);

	// std::vector<double> v(M*4);
	// v = param->get_com();

	//for (int i = 0; i < M*4; i++) {
		// std::cout << "value" << v[i] << "_and_log_bessel "<< get_log_bessel(v[i])<< std::endl;
	// }


	//INITIALIZE SPECIES CHARACTERISTICS
	vector<Species*> com(M);
	/* order of variable col, optimal q, nu, and dispersal */
	for(int i=0 ; i < M ; i++){
		Species *sp = new Species(pre_com[5*i], pre_com[5*i+1], pre_com[5*i+2], pre_com[5*i+3], pre_com[5*i+4]);
		com[i] = sp;
	}


	// for(int i=0 ; i < M ; i++){
	// 	std::cout << "the opt_q " << com[i]->get_opt_q()<< std::endl;
	// 	std::cout << "the log_bessel " << com[i]->get_sp_log_bessel()<< std::endl;
	// }

	//
	// Patch *p1 = new Patch(0.2,0.2,0.2,M,com);
	// Patch *p2 = new Patch(0.2,0.2,0.4,M,com);
	// Patch *p3 = new Patch(0.2,0.2,0.6,M,com);
	// Patch *p4 = new Patch(0.2,0.2,0.8,M,com);
	//
	// for (int i = 0; i < M; i++) {
	// 	std::cout << "the niche with of  species " << i << " is " << com[i]->get_nu() << std::endl;
	// 	std::cout << "the optimal q of  species " << i << " is " << com[i]->get_opt_q() << std::endl;
	// 	std::cout << " in patch 1 the fitness "<< i <<" is " << p1->get_fitness(i) << std::endl;
	// 	std::cout << " in patch 2 the fitness "<< i <<" is " << p2->get_fitness(i) << std::endl;
	// 	std::cout << " in patch 3 the fitness "<< i <<" is " << p3->get_fitness(i) << std::endl;
	// 	std::cout << " in patch 4 the fitness "<< i <<" is " << p4->get_fitness(i) << std::endl << std::endl;
	//
	// }
	//
	// double new_max_ls, max_ls = 0.0;
	// for (int i = 0; i < M; i++) {
	// 	new_max_ls = std::max(max_ls, com[i]->get_l());
	// 	max_ls = new_max_ls;
	// }
	// std::cout << "the max ls is " << max_ls << std::endl;

	double agg, midq, newq;
	std::cout << " enter aggreg ";
	cin>> agg;
	std::cout << "enter mid q ";
	cin>> midq;
	std::cout << "the new q is " << get_random_patch(midq, agg) <<  std::endl;

	// delete p1;
	// delete p2;



	//int comsize = (com.size())/4;
	//Parameter *param = new Parameter(size, M ,l_s, compet_z, simtime, est, aggreg, muR, rep,lambda, gammaH, com);
	simulation(param,com);

	delete param;
	for(int i=0 ; i < M ; i++){
		delete com[i]; // = sp;
	}

	return 0;
}
