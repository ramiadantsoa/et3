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
// #include <boost/multiprecision/cpp_dec_float.hpp>


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

/*  DESTROYING THE LANDSCAPE */

#include "section_destroy_land.h"

/*  CONTAINS THE MAIN SIMULATOR */

#include "simulator.h"


// typedef number<cpp_dec_float<100> > cpp_dec_float_100;

int main(){
	double z = 1.0, aggreg = 0.0, size = 6, simtime = 1000, muR = 0.1;
	short int est = 1,  M = 5, replicates = 0, com_id = 1;
	double tau = 40, lambda = 0.2, gammaH = 2.0;

  std::time_t seed_time = std::time(0);
  rng.seed(123456);
  // rng.seed(seed_time);

	// Parameter input
	// std::cout << "enter size ";
	// std::cin >> size;
	// std::cout << " \nenter patch size ";
	// std::cin >> lambda;
	// std::cout << " \nenter patch density ";
	// std::cin >> gammaH;
	// std::cout << "\nenter time ";
	// std::cin >> simtime;
	//
	// // Model types
	// std::cout << "\nenter competition ";
	// std::cin >> z;
	// std::cout << "\nenter aggregation "; //note that aggregation now takes values from 0 to 0.5, when 0 means aggregated, 0.5 means uniform.
	// std::cin >> aggreg;
	// //
	// // Enter community size and id
	// std::cout << "\nM ";
	// std::cin >> M;

	std::vector<double> pre_com(M*5);

	// For the moment assume a fixed community
// This is for a single generalist
	pre_com[0] = 0.5;
	pre_com[1] = 0.5;
	pre_com[2] = 0.5;
	pre_com[3] = 0.5;
	pre_com[4] = 0;

// This is for all (4) the specialists
	// colonization rate
	pre_com[0+5] = pre_com[0+2*5] = pre_com[0+3*5] = pre_com[0+4*5] = 0.5;
	// optimal q
	pre_com[1+5] = 0.;
	pre_com[1+2*5] = 0.25;
	pre_com[1+3*5] = 0.5;
	pre_com[1+4*5] = 0.75;
	// niche width
	pre_com[2+5] = pre_com[2+2*5] = pre_com[2+3*5] = pre_com[2+4*5] = 0.125;
	// dispersal range
	pre_com[3+5] = pre_com[3+2*5] = pre_com[3+3*5] = pre_com[3+4*5] = 0.125;
	// species id
	pre_com[4+5] = 1;
	pre_com[4+2*5] = 2;
	pre_com[4+3*5] = 3;
	pre_com[4+4*5] = 4;

	//INITIALIZE SPECIES CHARACTERISTICS
	vector<Species*> com(M);
	/* order of variable col, optimal q, nu, and dispersal */
	for(int i=0 ; i < M ; i++){
		Species *sp = new Species(pre_com[5*i], pre_com[5*i+1], pre_com[5*i+2], pre_com[5*i+3], pre_com[5*i+4]);
		com[i] = sp;
	}

	double destruct_param = 0.25;

	for (int rep = 0; rep < 1; rep++) {
		Parameter *param = new Parameter(com_id, M, z, est, aggreg, muR, tau, lambda, gammaH,
			size, simtime, rep);
		simulation(param,com, destruct_param);
		delete param;
	}

	for(int i = 0 ; i < M ; i++){
		delete com[i]; // = sp;
	}

	return 0;
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

// double agg, midq, newq;
// std::cout << " enter aggreg ";
// cin>> agg;
// std::cout << "enter mid q ";
// cin>> midq;
// std::cout << "the new q is " << get_random_patch(midq, agg) <<  std::endl;

// delete p1;
// delete p2;



//int comsize = (com.size())/4;
//Parameter *param = new Parameter(size, M ,l_s, compet_z, simtime, est, aggreg, muR, rep,lambda, gammaH, com);



// Enter species parameters
// std::vector<double> pre_com(M*5);
// for (int i = 0; i < M; i++) {
//
// 	pre_com[5*i] = 0.5; //colonization rate
// 	pre_com[5*i+1] = get_random(); //optimal q
// 	std::cout << "enter niche width nu " ;
// 	std::cin >> pre_com[5*i+2];//(i+1)*get_random(); // niche width nu, by convention, nu = -1 for generalist
// 	pre_com[5*i+3] = 0.2; // dispersal range
// 	pre_com[5*i+4] = i; // species id
//
// }
