#ifndef __GRID_H_
#define __GRID_H_

class Grid {
    public:
        Grid(double _size, double _max_l, int _M) {
            this->size = _size;
            this->max_l = _max_l;
            this->ncells = (int) ceil(size/max_l);
            this->cells = new vector<Patch*>[ncells*ncells]; // create a new array
            this->M = _M;
        }

        ~Grid() {
            for ( int i = 0; i < this->ncells; i++) {
                for (int j = 0; j < this->ncells; j++) {
                    vector<Patch*> *v = get_cell(i,j);
                    for(unsigned int k = 0; k < v->size(); k++) {
                        Patch *p = (*v)[k];
                        delete p;
                    }
                }
            }

            delete[] this->cells;
        }

        int get_number_of_cells() {
            return this->ncells;
        }

        void add_patch(Patch *p) {
            int i = (int) floor(p->get_x() / this->max_l);
            int j = (int) floor(p->get_y() / this->max_l);

            vector<Patch*> *vec = get_cell(i, j);
            vec->push_back(p);
        }

		void delete_patch( int i , int j , uInt k ){
			 vector<Patch*> *vec = get_cell(i, j);
             Patch *p = (*vec)[k];
			 vec->erase(vec->begin()+k); //CAN BE OPTIMIZED
             delete p;
		}


        vector<Patch*> *get_cell(int i, int j) {
            return &this->cells[i + j*ncells];
        }

		Patch *get_patch (int i , uInt k){
			vector<Patch*> *patches = &this->cells[i];
			return (*patches)[k];
		}


        uInt count_patches() {
            uInt number_of_patches = 0;
            for(int i = 0; i<ncells; i++) {
                for(int j = 0; j<ncells; j++) {
                    number_of_patches += get_cell(i,j)->size();
                }
            }
            return number_of_patches;
        }

        vector<uInt> count_species() {
            vector<uInt> count_vec(this->M, 0);
            for (int i = 0; i<ncells; i++) {
                for(int j = 0; j<ncells; j++) {
                    vector<Patch*> *patches = get_cell(i,j);
                    // Go through all patches in this cell
                    for(unsigned int k = 0; k<patches->size(); k++) {
                        Patch *p = (*patches)[k];
                        for(int m = 0; m < this->M; m++) {
                            count_vec[m] += p->occupancy[m];
                        }
                    }

                }
            }
            return count_vec;
        }

		uInt count_per_grid( int i, int j){
			uInt n_in_grid = 0;
			vector<Patch*> *patches = get_cell(i,j);
			n_in_grid =patches->size();
            return n_in_grid;
		}

		double colrate_per_grid (int i, int j){
			double col_grid =0.0;
			vector<Patch*> *patches = get_cell(i,j);

			for(unsigned int k = 0; k<patches->size(); k++) {
				Patch *p = (*patches)[k];
                for(int m = 0; m < this->M; m++) {
					col_grid += p->get_colonization(m);
				}
			}
			return col_grid;
		}

		/*void export_all_resource_units(){
			string all_resource_units; // stores the position and the type of all patches
			ostringstream convert_all_resource_units;
			string resource_units_file_name= "resource unit data.csv";

			 for (int i = 0; i<ncells; i++) {
                for(int j = 0; j<ncells; j++) {
                    vector<Patch*> *patches = get_cell(i,j);
                    // Go through all patches in this cell
                    for(unsigned int k = 0; k<patches->size(); k++) {
                        Patch *p = (*patches)[k];
						convert_all_resource_units<<p->get_x()<<","<<p->get_y()<<","<<p->get_quality()<<",";
                        for(int m = 0; m < this->M; m++) {
                            convert_all_resource_units<<p->occupancy[m]<<",";
                        }
						convert_all_resource_units<<"\n";
                    }

                }
            }

			all_resource_units.append(convert_all_resource_units.str());
			ofstream resource_units_file;
			resource_units_file.open(resource_units_file_name.c_str());
			resource_units_file<<all_resource_units;
			resource_units_file.close();

		}*/


    private:
        double size; /* The size of the side */
        double max_l; /* maximum dispersal range */
        int ncells; /* The grid has size ncells*ncells */
        int M;
        vector<Patch*> *cells;
};

Grid *initialize(Parameter *param, double resourceDensity,
  vector<double> occupancyFractions, vector<Grain*> *habitat, vector<Species*> com,
  double max_ls) {
	int ngrain = habitat->size(); // represent the number of forest patch that can produces resource unit

	int M = param->get_M();
	double size =param->get_size();

	unsigned int npatches = get_poisson(resourceDensity*size*size);//initialize the number of resource units (called patches)

	Grid *g = new Grid(size, max_ls, M);

    for (unsigned int i = 0; i<npatches; i++){

		// Draw random integer from size Habitat to choose in which habitat grain the new patch will be
		uInt chosen_grain = get_rand_integer(ngrain);

		//cout<<" the chosen grain is " << chosen_grain<<endl;

		Grain *grain = (*habitat)[chosen_grain];

		double lambda =grain->get_lambda();
		double x_pos =grain->get_x_pos();
		double y_pos =grain->get_y_pos();

    double _lambda = lambda*sqrt(get_random());
    double _angle = 2*PI*get_random();

		double x = x_pos +_lambda*cos(_angle);
		if(x<0.0){
			x+= size;
		}
		if( x > size ){
			x-= size;
		}

		double y = y_pos +_lambda*sin(_angle);
		if(y<0.0){
			y+= size ;
		}
		if( y > size ){
			y-=size;
		}

		double q =  get_random_patch( grain->get_q(), param->get_aggreg());

		//cout<<" patch quality "<< q<<endl;

    Patch *p = new Patch(x, y, q, M, com);

		//cout<< "x " << p->get_x() << " y " << p->get_y() << " q " << p->get_quality()<<" and fitness "<<p->get_fitness(1)<<"\n";

        //INITIALIZE SPECIES OCCUPANCY
		for (int m = 0; m < M; m++) {
            int is_occupied = 0;
            double epsilon = get_random();
            if (epsilon < occupancyFractions[m]) {
                is_occupied = 1;
            }
            p->occupancy[m] = is_occupied;
            //p->colonization[m] = 0.0;
        }

        g->add_patch(p);
    }

	//g->export_all_resource_units();

    return g;
}

void output(Grid *g) {
    uInt total_patches = g->count_patches();
    vector<uInt> species_counts = g->count_species();
    cout << total_patches;
    for (unsigned int i = 0; i<species_counts.size(); i++) {
        cout << ", \t" << species_counts[i];
    }
    cout << "\n";
}



#endif // __GRID_H_
