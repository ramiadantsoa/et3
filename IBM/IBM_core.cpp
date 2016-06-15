#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <time.h>
#include <sstream>
#include <string>
#include <assert.h>

#include "dSFMT/dSFMT.h"

#define PI 3.1415926535897932385

#define chop_0 0.000000000000001

using namespace std;

typedef unsigned long int uInt;

const int rng_seed = (const int)time(NULL);
dsfmt_t rng;

void init_rng() {
    dsfmt_init_gen_rand(&rng, rng_seed);
}

#include "simple_functions.h"

/*CLASS DECLARATION AND DEFINITION*/
#include "class_parameter.h"
#include "class_species.h"
#include "class_grain.h"
#include "initialize_incidence.h"
#include "class_patch_new.h"
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
	double compet_z =0.0, size, time, muR = 0.1, col_a = 0.5, diff_a = 1, w;
	short int est = 1, aggreg, M;
	double l_s = 0.5, lambda, gammaH;


	init_rng();

	cout<<"Compet z: ";
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

	//cout<<"\n";
	//cout<<"Community size: ";
	//cin>> M;
	M = 5;
	w = 2.0;

	cout<<"size: ";
	cin>>size;
	cout<<"time: ";
	cin>>time;

	short int rep = 0;

	Parameter *param = new Parameter(size, M ,l_s, compet_z, time, est, aggreg, muR, rep,col_a,diff_a,lambda, gammaH);
	simulation(param,w);
	delete param;


	return 0;
}
