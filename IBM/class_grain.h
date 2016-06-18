#ifndef __GRAIN_H_
#define __GRAIN_H_

class Grain{
private:
	double x_pos;
	double y_pos;
	double lambda;
	double q;

public:
	Grain(double _x, double _y, double _lambda, double _q){
		this->x_pos =_x;
		this->y_pos = _y;
		this->lambda = _lambda;
		this->q = _q;
	}

	double get_x_pos(){
		return this->x_pos;
	}

	double get_y_pos(){
		return this->y_pos;
	}

	double get_lambda(){
		return this->lambda;
	}

	double get_q(){
		return this->q;
	}
};

vector<Grain*> *initialize_habitat (Parameter *Initial_Parameter , double gammaH, double lambda){
	double SIZE = Initial_Parameter->get_size();
	int cSIZE= (int) ceil(SIZE);
	//double tau = Initial_Parameter->get_tau();
	int n_grains = ceil(SIZE*SIZE*gammaH);
	assert(n_grains >= 1);
	double x_loc, y_loc, q;

	// string all_patches; // stores the position and the type of all patches
	// ostringstream convert_all_patches;
	// string patch_file_name= "patch data.csv";

	vector<Grain*> *Habitat = new vector<Grain*>;

	for(int i = 0; i < n_grains ; i++){
		x_loc= get_random()*SIZE;
		y_loc =get_random()*SIZE;
		q =get_random();

		Grain *grain =new Grain(x_loc, y_loc,lambda,q);

		Habitat->push_back(grain);

		// convert_all_patches<<x_loc<<","<<y_loc<<","<<q<<"\n";
	}

	// all_patches.append(convert_all_patches.str());
	// ofstream patch_file;
	// patch_file.open(patch_file_name.c_str());
	// patch_file<<all_patches;
	// patch_file.close();

	return Habitat;
}


#endif // __GRAIN_H_
